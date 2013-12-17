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
        return new AlleleDepthBuilder(false, numTaxa, numSites, maxNumAlleles);
    }

    public static AlleleDepthBuilder getNucleotideInstance(int numTaxa, int numSites) {
        return new AlleleDepthBuilder(false, numTaxa, numSites, NucleotideAlignmentConstants.NUMBER_NUCLEOTIDE_ALLELES);
    }

    public static AlleleDepthBuilder getHDF5NucleotideInstance(int numSites) {
        return new AlleleDepthBuilder(true, numSites, NucleotideAlignmentConstants.NUMBER_NUCLEOTIDE_ALLELES);
    }

    /**
     * Set allele value for taxon, site, and allele. Value will be translated
     * using AlleleDepthUtil.
     *
     * @param taxon taxon
     * @param site site
     * @param allele allele
     * @param value value
     *
     * @return builder
     */
    public AlleleDepthBuilder setDepth(int taxon, int site, byte allele, int value) {
        if (myIsHDF5) {
            throw new IllegalStateException("AlleleDepthBuilder: setDepth: use addTaxon for HDF5 files.");
        }
        myDepths.get(taxon).set(site, allele, AlleleDepthUtil.depthIntToByte(value));
        return this;
    }

    /**
     * Set allele for the all sites and alleles for a taxon simultaneously.
     * Value will be translated using AlleleDepthUtil.
     *
     * @param taxon Index of taxon
     * @param value array[sites][allele] of all values
     *
     * @return builder
     */
    public AlleleDepthBuilder setDepth(int taxon, int[][] value) {
        if (myIsHDF5) {
            throw new IllegalStateException("AlleleDepthBuilder: setDepth: use addTaxon for HDF5 files.");
        }
        int numSites = value.length;
        if (numSites != myNumSites) {
            throw new IllegalArgumentException("AlleleDepthBuilder: setDepth: value number of sites: " + numSites + " should have: " + myNumSites);
        }
        int numAlleles = value[0].length;
        if (numAlleles != myMaxNumAlleles) {
            throw new IllegalArgumentException("AlleleDepthBuilder: setDepth: value number of alleles: " + numAlleles + " should have: " + myMaxNumAlleles);
        }
        for (int s = 0; s < myNumSites; s++) {
            for (int a = 0; a < myMaxNumAlleles; a++) {
                setDepth(taxon, s, (byte) a, value[s][a]);
            }
        }
        return this;
    }

    /**
     * Set allele value for taxon, site, and allele. Value should have already
     * been translated using AlleleDepthUtil.
     *
     * @param taxon taxon
     * @param site site
     * @param allele allele
     * @param value value
     *
     * @return builder
     */
    public AlleleDepthBuilder setDepth(int taxon, int site, byte allele, byte value) {
        if (myIsHDF5) {
            throw new IllegalStateException("AlleleDepthBuilder: setDepth: use addTaxon for HDF5 files.");
        }
        myDepths.get(taxon).set(site, allele, value);
        return this;
    }

    /**
     * Set allele for the all sites and alleles for a taxon simultaneously.
     * Values should have already been translated using AlleleDepthUtil.
     *
     * @param taxon Index of taxon
     * @param value array[sites][allele] of all values
     *
     * @return builder
     */
    public AlleleDepthBuilder setDepth(int taxon, byte[][] value) {
        if (myIsHDF5) {
            throw new IllegalStateException("AlleleDepthBuilder: setDepth: use addTaxon for HDF5 files.");
        }
        int numSites = value.length;
        if (numSites != myNumSites) {
            throw new IllegalArgumentException("AlleleDepthBuilder: setDepth: value number of sites: " + numSites + " should have: " + myNumSites);
        }
        int numAlleles = value[0].length;
        if (numAlleles != myMaxNumAlleles) {
            throw new IllegalArgumentException("AlleleDepthBuilder: setDepth: value number of alleles: " + numAlleles + " should have: " + myMaxNumAlleles);
        }
        for (int s = 0; s < myNumSites; s++) {
            for (int a = 0; a < myMaxNumAlleles; a++) {
                setDepth(taxon, s, (byte) a, value[s][a]);
            }
        }
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
