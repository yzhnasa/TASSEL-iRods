/*
 *  AlleleDepthBuilder
 */
package net.maizegenetics.dna.snp.depth;

import java.util.ArrayList;
import java.util.List;

import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.util.SuperByteMatrix;
import net.maizegenetics.util.SuperByteMatrixBuilder;

/**
 * Builder to store information on DNA read depths.
 *
 * @author Terry Casstevens
 */
public class AlleleDepthBuilder {

    private List<SuperByteMatrix> myDepths = null;
    private final boolean myIsHDF5;
    private final int myMaxNumAlleles;
    private final int myNumSites;
    private int myNumTaxa = 0;

    private AlleleDepthBuilder(boolean hdf5, int numSites, int maxNumAlleles) {
        myIsHDF5 = hdf5;
        myMaxNumAlleles = maxNumAlleles;
        myNumSites = numSites;

        if (myIsHDF5) {

        } else {
            myDepths = new ArrayList<>();
        }
    }

    private AlleleDepthBuilder(boolean hdf5, int numTaxa, int numSites, int maxNumAlleles) {
        myIsHDF5 = hdf5;
        myMaxNumAlleles = maxNumAlleles;
        myNumSites = numSites;
        myNumTaxa = numTaxa;

        if (myIsHDF5) {

        } else {
            myDepths = new ArrayList<>();
            for (int i = 0; i < myNumTaxa; i++) {
                myDepths.add(SuperByteMatrixBuilder.getInstance(myNumSites, myMaxNumAlleles));
            }
        }
    }

    public static AlleleDepthBuilder getInstance(int numTaxa, int numSites, int maxNumAlleles) {
        return new AlleleDepthBuilder(numTaxa, numSites, maxNumAlleles);
    }

    public static AlleleDepthBuilder getNucleotideInstance(int numTaxa, int numSites) {
        return new AlleleDepthBuilder(numTaxa, numSites, NucleotideAlignmentConstants.NUMBER_NUCLEOTIDE_ALLELES);
    }

    public AlleleDepthBuilder setDepth(int taxon, int site, byte allele, byte value) {
        myDepths[allele].set(taxon, site, value);
        return this;
    }

    /**
     * Set allele for the all sites and alleles for a taxon simultaneously.
     *
     * @param taxon Index of taxon
     * @param value array[sites][allele] of all values
     * @return
     */
    public AlleleDepthBuilder setDepth(int taxon, byte[][] value) {
        return this;
    }

    public AlleleDepth build() {
        if (myIsHDF5) {
            return null;
        } else {
            SuperByteMatrix[] temp = new SuperByteMatrix[myDepths.size()];
            temp = myDepths.toArray(temp);
            myDepths = null;
            return new MemoryAlleleDepth(temp);
        }
    }
    }
}
