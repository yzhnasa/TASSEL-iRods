/*
 *  MemoryAlleleDepth
 */
package net.maizegenetics.dna.snp.depth;

import net.maizegenetics.util.SuperByteMatrix;

/**
 * In memory allele depth class. In memory allele depth occupies 6-fold more
 * memory than the allele calls.
 *
 * @author Terry Casstevens
 */
class MemoryAlleleDepth implements AlleleDepth {

    private final SuperByteMatrix[] myDepths;
    private final int myNumAlleles;

    MemoryAlleleDepth(SuperByteMatrix[] depths) {
        myDepths = depths;
        myNumAlleles = depths.length;
    }

    @Override
    public byte[] getDepthForAlleles(int taxon, int site) {
        byte[] result = new byte[myNumAlleles];
        for (int a = 0; a < myNumAlleles; a++) {
            result[a] = myDepths[a].get(taxon, site);
        }
        return result;
    }

    @Override
    public byte getDepthForAllele(int taxon, int site, int allele) {
        return myDepths[allele].get(taxon, site);
    }
}
