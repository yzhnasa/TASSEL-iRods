/*
 *  MemoryAlleleDepth
 */
package net.maizegenetics.pal.alignment.depth;

/**
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
