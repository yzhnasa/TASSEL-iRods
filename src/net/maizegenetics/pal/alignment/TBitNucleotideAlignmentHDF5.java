/*
 * TBitNucleotideAlignmentHDF5
 */
package net.maizegenetics.pal.alignment;

import ch.systemsx.cisd.hdf5.IHDF5Reader;
import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.util.OpenBitSet;

/**
 *
 * @author terry
 */
public class TBitNucleotideAlignmentHDF5 extends BitAlignmentHDF5 {

    private static final long serialVersionUID = -5197800047652332969L;

    protected TBitNucleotideAlignmentHDF5(IHDF5Reader hdf5, IdGroup idGroup, byte[][] alleles, GeneticMap map, byte[] reference, int[] variableSites, int maxNumAlleles, Locus[] loci, int[] lociOffsets, String[] snpIDs, boolean retainRareAlleles) {
        super(hdf5, idGroup, alleles, map, reference, NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles);
    }

    @Override
    public byte[] getBaseArray(int taxon, int site) {
        OpenBitSet[] tBitData = getCachedTaxon(taxon);
        byte[] result = new byte[2];
        result[0] = Alignment.UNKNOWN_ALLELE;
        result[1] = Alignment.UNKNOWN_ALLELE;
        try {
            int count = 0;
            for (int i = 0; i < myMaxNumAlleles; i++) {
                if (tBitData[i].fastGet(site)) {
                    if (count == 0) {
                        result[1] = myAlleles[site][i];
                    }
                    result[count++] = myAlleles[site][i];
                }
            }

            // Check For Rare Allele
            if (retainsRareAlleles() && tBitData[myMaxNumAlleles].fastGet(site)) {
                if (count == 0) {
                    result[1] = Alignment.RARE_ALLELE;
                }
                result[count] = Alignment.RARE_ALLELE;
            }

        } catch (IndexOutOfBoundsException e) {
            throw new IllegalStateException("TBitNucleotideAlignmentHDF5: getBaseArray: bit sets indicate more than two alleles for taxon: " + taxon + "   site: " + site);
        }
        return result;
    }

    @Override
    public String getBaseAsString(int taxon, int site) {
        return NucleotideAlignmentConstants.getNucleotideIUPAC(getBase(taxon, site));
    }

    @Override
    public String getDiploidAsString(int site, byte value) {
        return NucleotideAlignmentConstants.getNucleotideIUPAC(value);
    }

    @Override
    public boolean isIndel(int site) {
        int[][] alleles = getAllelesSortedByFrequency(site);
        int numAlleles = Math.min(alleles[0].length, 2);
        for (int i = 0; i < numAlleles; i++) {
            if ((alleles[0][i] == NucleotideAlignmentConstants.INSERT_ALLELE) || (alleles[0][i] == NucleotideAlignmentConstants.GAP_ALLELE)) {
                return true;
            }
        }
        return false;
    }
}
