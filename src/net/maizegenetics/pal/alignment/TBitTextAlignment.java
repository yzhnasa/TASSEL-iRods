/*
 * TBitTextAlignment
 */
package net.maizegenetics.pal.alignment;

import net.maizegenetics.pal.ids.IdGroup;

/**
 *
 * @author terry
 */
public class TBitTextAlignment extends TBitAlignment {

    protected TBitTextAlignment(Alignment a, int maxNumAlleles, boolean retainRareAlleles, boolean isFinalized) {
        super(a, maxNumAlleles, retainRareAlleles, isFinalized);
    }

    protected TBitTextAlignment(IdGroup idGroup, byte[][] data, GeneticMap map, byte[] reference, String[][] alleleStates, int[] variableSites, int maxNumAlleles, Locus[] loci, int[] lociOffsets, String[] snpIDs, boolean retainRareAlleles, boolean isFinalized) {
        super(idGroup, data, map, reference, alleleStates, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles, isFinalized);
    }

    @Override
    public String getBaseAsString(int taxon, int site) {
        byte[] temp = getBaseArray(taxon, site);
        return myAlleleStates[site][temp[0]] + ":" + myAlleleStates[site][temp[1]];
    }

    @Override
    public String[] getBaseAsStringArray(int taxon, int site) {
        byte[] temp = getBaseArray(taxon, site);
        return new String[]{myAlleleStates[site][temp[0]], myAlleleStates[site][temp[1]]};
    }

    @Override
    public String getBaseAsString(int site, byte value) {
        return myAlleleStates[site][value];
    }
}
