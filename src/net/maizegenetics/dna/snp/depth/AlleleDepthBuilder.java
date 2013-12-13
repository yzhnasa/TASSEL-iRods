/*
 *  AlleleDepthBuilder
 */
package net.maizegenetics.dna.snp.depth;

import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;

/**
 * Builder to store information on DNA read depths.
 * TODO this needs methods to support incremental depth building.  What about HDF5 depth building?
 * @author Terry Casstevens
 */
public class AlleleDepthBuilder {

    private byte[][][] myDepths = null;

    private AlleleDepthBuilder(int numTaxa, int numSites, int maxNumAlleles) {
        myDepths = new byte[numTaxa][numSites][maxNumAlleles];
    }

    public static AlleleDepthBuilder getInstance(int numTaxa, int numSites, int maxNumAlleles) {
        return new AlleleDepthBuilder(numTaxa, numSites, maxNumAlleles);
    }

    public static AlleleDepthBuilder getNucleotideInstance(int numTaxa, int numSites) {
        return new AlleleDepthBuilder(numTaxa, numSites, NucleotideAlignmentConstants.NUMBER_NUCLEOTIDE_ALLELES);
    }

    public AlleleDepthBuilder setDepth(int taxon, int site, byte allele, byte value) {
        myDepths[taxon][site][allele] = value;
        return this;
    }

    /**
     * Set allele for the all sites and alleles for a taxon simultaneously.
     * @param taxon Index of taxon
     * @param value array[sites][allele] of all values
     * @return
     */
    public AlleleDepthBuilder setDepth(int taxon, byte[][] value) {
        myDepths[taxon] = value;
        return this;
    }

    public AlleleDepth build() {
        byte[][][] temp = myDepths;
        myDepths = null;
        return new MemoryAlleleDepth(temp);
    }
}
