/*
 * TBitNucleotideAlignment
 */
package net.maizegenetics.pal.alignment;

import net.maizegenetics.pal.ids.IdGroup;

/**
 *
 * @author terry
 */
public class TBitNucleotideAlignment extends TBitAlignment {

    protected TBitNucleotideAlignment(Alignment a, int maxNumAlleles, boolean retainRareAlleles) {
        super(a, maxNumAlleles, retainRareAlleles);
    }

    protected TBitNucleotideAlignment(IdGroup idGroup, byte[][] data, GeneticMap map, byte[] reference, String[][] alleleStates, int[] variableSites, int maxNumAlleles, Locus[] loci, int[] lociOffsets, String[] snpIDs, boolean retainRareAlleles) {
        super(idGroup, data, map, reference, alleleStates, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles);
    }

    @Override
    public String getBaseAsString(int taxon, int site) {
        return NucleotideAlignmentConstants.getNucleotideIUPAC(getBase(taxon, site));
    }
    
    public String getDiploidAsString(int site, byte value) {
        return NucleotideAlignmentConstants.getNucleotideIUPAC(value);
    }
    
    @Override
    public boolean isIndel(int site) {
        int[][] alleles = getAllelesSortedByFrequency(site);
        int numAlleles = Math.min(alleles[0].length, 2);
        for (int i = 0; i < numAlleles; i++) {
            if ((alleles[0][i] == 4) || (alleles[0][i] == 5)) {
                return true;
            }
        }
        return false;
    }

}
