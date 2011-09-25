/*
 * SBitTextAlignment
 */
package net.maizegenetics.pal.alignment;

import net.maizegenetics.pal.ids.IdGroup;

/**
 *
 * @author terry
 */
public class SBitTextAlignment extends SBitAlignment {

    protected SBitTextAlignment(Alignment a, int maxNumAlleles, boolean retainRareAlleles) {
        super(a, maxNumAlleles, retainRareAlleles);
    }

    protected SBitTextAlignment(IdGroup idGroup, byte[][] data, GeneticMap map, byte[] reference, String[][] alleleStates, int[] variableSites, int maxNumAlleles, Locus[] loci, int[] lociOffsets, String[] snpIDs, boolean retainRareAlleles) {
        super(idGroup, data, map, reference, alleleStates, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles);
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
