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
class MemoryAlleleDepth extends AbstractAlleleDepth {

    private final SuperByteMatrix[] myDepths;

    MemoryAlleleDepth(SuperByteMatrix[] depths) {
        super(6);
        myDepths = depths;
    }

    @Override
    public byte depthForAlleleByte(int taxon, int site, int allele) {
        return myDepths[taxon].get(site, allele);
    }
}
