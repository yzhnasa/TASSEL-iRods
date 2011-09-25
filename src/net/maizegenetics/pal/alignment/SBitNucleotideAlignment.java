/*
 * SBitNucleotideAlignment
 */
package net.maizegenetics.pal.alignment;

import net.maizegenetics.pal.ids.IdGroup;

/**
 *
 * @author terry
 */
public class SBitNucleotideAlignment extends SBitAlignment {

    protected SBitNucleotideAlignment(Alignment a, int maxNumAlleles, boolean retainRareAlleles) {
        super(a, maxNumAlleles, retainRareAlleles);
    }

    protected SBitNucleotideAlignment(IdGroup idGroup, byte[][] data, GeneticMap map, byte[] reference, String[][] alleleStates, int[] variableSites, int maxNumAlleles, Locus[] loci, int[] lociOffsets, String[] snpIDs, boolean retainRareAlleles) {
        super(idGroup, data, map, reference, alleleStates, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles);
    }

    @Override
    public String getBaseAsString(int taxon, int site) {
        return NucleotideAlignmentConstants.getNucleotideIUPAC(getBase(taxon, site));
    }

    @Override
    public String[] getBaseAsStringArray(int taxon, int site) {
        return new String[]{getBaseAsString(taxon, site)};
    }
}
