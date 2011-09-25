/*
 * SBitAlignment
 */
package net.maizegenetics.pal.alignment;

import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;
import net.maizegenetics.util.UnmodifiableBitSet;

/**
 * This data alignment is optimized for operations
 * involving lots of SNPs from a taxon - imputation, genetic distance, kinship, diversity, etc.
 * It is not optimized for LD or association mapping.
 *
 * @author terry
 */
public class SBitAlignment extends AbstractAlignment {

    private OpenBitSet[][] myData;
    private int myNumDataRows;
    private boolean myIsDirty = false;

    protected SBitAlignment(Alignment a, int maxNumAlleles, boolean retainRareAlleles) {
        super(a, maxNumAlleles, retainRareAlleles);
        loadAlleles(a);
    }

    protected SBitAlignment(IdGroup idGroup, byte[][] data, GeneticMap map, byte[] reference, String[][] alleleStates, int[] variableSites, int maxNumAlleles, Locus[] loci, int[] lociOffsets, String[] snpIDs, boolean retainRareAlleles) {
        super(idGroup, data, map, reference, alleleStates, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles);
        loadAlleles(data);
    }

    public static SBitAlignment getInstance(Alignment a) {
        return SBitAlignment.getInstance(a, a.getMaxNumAlleles(), a.retainsRareAlleles());
    }

    public static SBitAlignment getInstance(Alignment a, int maxNumAlleles, boolean retainRareAlleles) {
        String[][] alleleStates = a.getAlleleEncodings();
        if ((alleleStates == null) || (alleleStates.length == 0)) {
            throw new IllegalArgumentException("SBitAlignment: init: allele states can not be empty.");
        }
        if (alleleStates.length == 1) {
            return new SBitAlignment(a, maxNumAlleles, retainRareAlleles);
        } else {
            return new SBitTextAlignment(a, maxNumAlleles, retainRareAlleles);
        }
    }

    public static SBitAlignment getInstance(IdGroup idGroup, byte[][] data, GeneticMap map, byte[] reference, String[][] alleleStates, int[] variableSites, int maxNumAlleles, Locus[] loci, int[] lociOffsets, String[] snpIDs, boolean retainRareAlleles) {
        if ((alleleStates == null) || (alleleStates.length == 0)) {
            throw new IllegalArgumentException("SBitAlignment: init: allele states can not be empty.");
        }
        if (alleleStates.length == 1) {
            return new SBitAlignment(idGroup, data, map, reference, alleleStates, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles);
        } else {
            return new SBitTextAlignment(idGroup, data, map, reference, alleleStates, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles);
        }
    }

    public static SBitAlignment getNucleotideInstance(IdGroup idGroup, byte[][] data, GeneticMap map, byte[] reference, int[] variableSites, int maxNumAlleles, Locus[] loci, int[] lociOffsets, String[] snpIDs, boolean retainRareAlleles) {
        return new SBitNucleotideAlignment(idGroup, data, map, reference, NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles);
    }

    public static SBitAlignment getNucleotideInstance(IdGroup idGroup, String[][] data, GeneticMap map, byte[] reference, int[] variableSites, int maxNumAlleles, Locus[] loci, int[] lociOffsets, String[] snpIDs, boolean retainRareAlleles) {

        if ((maxNumAlleles < 1) || (maxNumAlleles > NucleotideAlignmentConstants.NUMBER_NUCLEOTIDE_ALLELES)) {
            throw new IllegalArgumentException("SBitAlignment: getNucleotideInstance: max number of alleles must be between 1 and 14 inclusive: " + maxNumAlleles);
        }

        if ((data == null) || (data.length == 0)) {
            throw new IllegalArgumentException("SBitAlignment: getNucleotideInstance: data can not be empty.");
        }

        if (data.length != idGroup.getIdCount()) {
            throw new IllegalArgumentException("SBitAlignment: getNucleotideInstance: data rows not equal to number of identifers.");
        }

        byte[][] dataBytes = AlignmentUtils.getDataBytes(data, NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES, NucleotideAlignmentConstants.NUMBER_NUCLEOTIDE_ALLELES);

        return SBitAlignment.getNucleotideInstance(idGroup, dataBytes, map, reference, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles);

    }

    public static SBitAlignment getNucleotideInstance(IdGroup idGroup, String[] data, GeneticMap map, byte[] reference, int[] variableSites, int maxNumAlleles, Locus[] loci, int[] lociOffsets, String[] snpIDs, boolean retainRareAlleles) {

        if ((maxNumAlleles < 1) || (maxNumAlleles > NucleotideAlignmentConstants.NUMBER_NUCLEOTIDE_ALLELES)) {
            throw new IllegalArgumentException("SBitAlignment: getNucleotideInstance: max number of alleles must be between 1 and 14 inclusive: " + maxNumAlleles);
        }

        if ((data == null) || (data.length == 0)) {
            throw new IllegalArgumentException("SBitAlignment: getNucleotideInstance: data can not be empty.");
        }

        if (data.length != idGroup.getIdCount()) {
            throw new IllegalArgumentException("SBitAlignment: getNucleotideInstance: data rows not equal to number of identifers.");
        }

        byte[][] dataBytes = AlignmentUtils.getDataBytes(data);

        return SBitAlignment.getNucleotideInstance(idGroup, dataBytes, map, reference, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles);

    }

    public static SBitAlignment getInstance(IdGroup idGroup, String[][] data, GeneticMap map, byte[] reference, int[] variableSites, int maxNumAlleles, Locus[] loci, int[] lociOffsets, String[] snpIDs, boolean retainRareAlleles) {

        if ((maxNumAlleles < 1) || (maxNumAlleles > 14)) {
            throw new IllegalArgumentException("SBitAlignment: getInstance: max number of alleles must be between 1 and 14 inclusive: " + maxNumAlleles);
        }

        if ((data == null) || (data.length == 0)) {
            throw new IllegalArgumentException("SBitAlignment: getInstance: data can not be empty.");
        }

        if (data.length != idGroup.getIdCount()) {
            throw new IllegalArgumentException("SBitAlignment: getInstance: data rows not equal to number of identifers.");
        }

        String[][] alleleStates = AlignmentUtils.getAlleleStates(data, maxNumAlleles);

        byte[][] dataBytes = AlignmentUtils.getDataBytes(data, alleleStates, maxNumAlleles);

        return SBitAlignment.getInstance(idGroup, dataBytes, map, reference, alleleStates, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles);

    }

    private void loadAlleles(byte[][] data) {

        myNumDataRows = myMaxNumAlleles;
        if (retainsRareAlleles()) {
            myNumDataRows++;
        }
        myData = new OpenBitSet[myNumDataRows][myNumSites];
        for (int al = 0; al < myNumDataRows; al++) {
            for (int s = 0; s < myNumSites; s++) {
                myData[al][s] = new OpenBitSet(getSequenceCount());
            }
        }
        byte[] cb = new byte[2];
        for (int s = 0; s < myNumSites; s++) {
            for (int t = 0, n = getSequenceCount(); t < n; t++) {
                cb[0] = (byte) ((data[t][s] >>> 4) & 0xf);
                cb[1] = (byte) (data[t][s] & 0xf);
                for (int i = 0; i < 2; i++) {
                    if (cb[i] != Alignment.UNKNOWN_ALLELE) {
                        boolean isRare = true;
                        for (int j = 0; j < myMaxNumAlleles; j++) {
                            if (cb[i] == myAlleles[s][j]) {
                                myData[j][s].fastSet(t);
                                isRare = false;
                                break;
                            }
                        }
                        if (isRare && retainsRareAlleles()) {
                            myData[myMaxNumAlleles][s].fastSet(t);
                        }
                    }
                }
            }
        }

    }

    private void loadAlleles(Alignment a) {

        myNumDataRows = myMaxNumAlleles;
        if (retainsRareAlleles()) {
            myNumDataRows++;
        }
        myData = new OpenBitSet[myNumDataRows][myNumSites];
        for (int al = 0; al < myNumDataRows; al++) {
            for (int s = 0; s < myNumSites; s++) {
                myData[al][s] = new OpenBitSet(getSequenceCount());
            }
        }
        for (int s = 0; s < myNumSites; s++) {
            for (int t = 0, n = getSequenceCount(); t < n; t++) {
                byte[] cb = a.getBaseArray(t, s);
                for (int i = 0; i < 2; i++) {
                    if (cb[i] != Alignment.UNKNOWN_ALLELE) {
                        boolean isRare = true;
                        for (int j = 0; j < myMaxNumAlleles; j++) {
                            if (cb[i] == myAlleles[s][j]) {
                                myData[j][s].fastSet(t);
                                isRare = false;
                                break;
                            }
                        }
                        if (isRare && retainsRareAlleles()) {
                            myData[myMaxNumAlleles][s].fastSet(t);
                        }
                    }
                }
            }
        }

    }

    @Override
    public byte getBase(int taxon, int site) {
        byte[] temp = getBaseArray(taxon, site);
        return (byte) ((temp[0] << 4) | temp[1]);
    }

    @Override
    public byte[] getBaseArray(int taxon, int site) {
        byte[] result = new byte[2];
        result[0] = Alignment.UNKNOWN_ALLELE;
        result[1] = Alignment.UNKNOWN_ALLELE;
        try {
            int count = 0;
            for (int i = 0; i < myMaxNumAlleles; i++) {
                if (myData[i][site].fastGet(taxon)) {
                    if (count == 0) {
                        result[1] = myAlleles[site][i];
                    }
                    result[count++] = myAlleles[site][i];
                }
            }

            // Check For Rare Allele
            if (retainsRareAlleles() && myData[myMaxNumAlleles][site].fastGet(taxon)) {
                if (count == 0) {
                    result[1] = Alignment.RARE_ALLELE;
                }
                result[count] = Alignment.RARE_ALLELE;
            }

        } catch (IndexOutOfBoundsException e) {
            throw new IllegalStateException("SBitAlignment: getBaseArray: bit sets indicate more than two alleles for taxon: " + taxon + "   site: " + site);
        }
        return result;
    }

    @Override
    public BitSet getAllelePresenceForAllTaxa(int site, int alleleNumber) {
        return UnmodifiableBitSet.getInstance(myData[alleleNumber][site]);
    }
}
