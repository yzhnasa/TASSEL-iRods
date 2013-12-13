/*
 *  MemoryAlleleDepth
 */
package net.maizegenetics.dna.snp.depth;

/**
 * In memory allele depth class.  In memory allele depth occupies 6-fold more memory than the allele calls.
 *
 * @author Terry Casstevens
 */
class MemoryAlleleDepth implements AlleleDepth {

    private final byte[][][] myDepths;

    MemoryAlleleDepth(byte[][][] depths) {
        myDepths = depths;
    }

    @Override
    public byte[] getDepthForAlleles(int taxon, int site) {
        return myDepths[taxon][site];
    }
}
