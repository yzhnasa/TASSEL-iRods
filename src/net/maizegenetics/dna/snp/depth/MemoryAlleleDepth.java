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
    public int[] getDepthForAlleles(int taxon, int site) {
        int[] result = new int[myNumAlleles];
        for (int a = 0; a < myNumAlleles; a++) {
            result[a] = getDepthForAllele(taxon, site, a);
        }
        return result;
    }

    @Override
    public int getDepthForAllele(int taxon, int site, int allele) {
        return AlleleDepthUtil.depthByteToInt(myDepths[taxon].get(site, allele));
    }

    @Override
    public int getDepth(int taxon, int site) {
        int result = 0;
        for (int a = 0; a < myNumAlleles; a++) {
            result += getDepthForAllele(taxon, site, a);
        }
        return result;
    }
}
