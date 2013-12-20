/*
 *  AbstractAlleleDepth
 */
package net.maizegenetics.dna.snp.depth;

/**
 *
 * @author Terry Casstevens
 */
public abstract class AbstractAlleleDepth implements AlleleDepth {
    
    private final int myNumAlleles;
    
    AbstractAlleleDepth(int numAlleles) {
        myNumAlleles = numAlleles;
    }
    
    @Override
    public int[] depthForAlleles(int taxon, int site) {
        int[] result = new int[myNumAlleles];
        for (int a = 0; a < myNumAlleles; a++) {
            result[a] = depthForAllele(taxon, site, a);
        }
        return result;
    }
    
    @Override
    public int depthForAllele(int taxon, int site, int allele) {
        return AlleleDepthUtil.depthByteToInt(depthForAlleleByte(taxon, site, allele));
    }
    
    @Override
    public int depth(int taxon, int site) {
        int result = 0;
        for (int a = 0; a < myNumAlleles; a++) {
            result += depthForAllele(taxon, site, a);
        }
        return result;
    }
    
    @Override
    public int depthForTaxon(int taxon) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    @Override
    public int depthForSite(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
}
