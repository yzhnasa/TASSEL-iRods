/*
 *  AlleleDepthBuilder
 */
package net.maizegenetics.dna.snp.depth;

import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.util.SuperByteMatrix;
import net.maizegenetics.util.SuperByteMatrixBuilder;

/**
 * Builder to store information on DNA read depths. TODO this needs methods to
 * support incremental depth building. What about HDF5 depth building?
 *
 * @author Terry Casstevens
 */
public class AlleleDepthBuilder {

    private SuperByteMatrix[] myDepths = null;

    private AlleleDepthBuilder(int numTaxa, int numSites, int maxNumAlleles) {
        myDepths = new SuperByteMatrix[maxNumAlleles];
        for (int i = 0; i < maxNumAlleles; i++) {
            myDepths[i] = SuperByteMatrixBuilder.getInstance(numTaxa, numSites);
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
        SuperByteMatrix[] temp = myDepths;
        myDepths = null;
        return new MemoryAlleleDepth(temp);
    }
}
