/*
 * TBitAlignment
 */
package net.maizegenetics.pal.alignment;

import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;
import net.maizegenetics.util.UnmodifiableBitSet;

/**
 *
 * @author terry
 */
public class TBitAlignment extends AbstractAlignment {

    private OpenBitSet[][] myData;
    private int myNumDataRows;

    protected TBitAlignment(Alignment a, int maxNumAlleles, boolean retainRareAlleles) {
        super(a, maxNumAlleles, retainRareAlleles);
        loadAlleles(a);
    }

    protected TBitAlignment(IdGroup idGroup, byte[][] data, GeneticMap map, byte[] reference, String[][] alleleStates, int[] variableSites, int maxNumAlleles, Locus[] loci, int[] lociOffsets, String[] snpIDs, boolean retainRareAlleles) {
        super(idGroup, data, map, reference, alleleStates, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles);
        loadAlleles(data);
    }

    public static TBitAlignment getInstance(Alignment a) {
        return TBitAlignment.getInstance(a, a.getMaxNumAlleles(), a.retainsRareAlleles());
    }

    public static TBitAlignment getInstance(Alignment a, int maxNumAlleles, boolean retainRareAlleles) {

        if ((a instanceof TBitAlignment) && (a.getMaxNumAlleles() == maxNumAlleles) && (a.retainsRareAlleles() == retainRareAlleles)) {
            return (TBitAlignment) a;
        }

        String[][] alleleStates = a.getAlleleEncodings();
        if ((alleleStates == null) || (alleleStates.length == 0)) {
            throw new IllegalStateException("TBitAlignment: init: allele states should not be empty.");
        }

        boolean isNucleotide = false;
        if (alleleStates.length == 1) {
            isNucleotide = true;
            if (alleleStates[0].length == NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES[0].length) {
                for (int i = 0; i < alleleStates.length; i++) {
                    if (!alleleStates[0][i].equals(NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES[0][i])) {
                        isNucleotide = false;
                    }
                }
            }

        }

        if (isNucleotide) {
            return new TBitNucleotideAlignment(a, maxNumAlleles, retainRareAlleles);
        } else if (alleleStates.length == 1) {
            return new TBitAlignment(a, maxNumAlleles, retainRareAlleles);
        } else {
            return new TBitTextAlignment(a, maxNumAlleles, retainRareAlleles);
        }

    }

    public static TBitAlignment getInstance(IdGroup idGroup, byte[][] data, GeneticMap map, byte[] reference, String[][] alleleStates, int[] variableSites, int maxNumAlleles, Locus[] loci, int[] lociOffsets, String[] snpIDs, boolean retainRareAlleles) {
        if ((alleleStates == null) || (alleleStates.length == 0)) {
            throw new IllegalArgumentException("TBitAlignment: init: allele states can not be empty.");
        }
        if (alleleStates.length == 1) {
            return new TBitAlignment(idGroup, data, map, reference, alleleStates, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles);
        } else {
            return new TBitTextAlignment(idGroup, data, map, reference, alleleStates, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles);
        }
    }

    public static TBitAlignment getNucleotideInstance(IdGroup idGroup, byte[][] data, GeneticMap map, byte[] reference, int[] variableSites, int maxNumAlleles, Locus[] loci, int[] lociOffsets, String[] snpIDs, boolean retainRareAlleles) {
        return new TBitNucleotideAlignment(idGroup, data, map, reference, NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles);
    }

    public static TBitAlignment getNucleotideInstance(IdGroup idGroup, String[][] data, GeneticMap map, byte[] reference, int[] variableSites, int maxNumAlleles, Locus[] loci, int[] lociOffsets, String[] snpIDs, boolean retainRareAlleles) {

        if ((maxNumAlleles < 1) || (maxNumAlleles > NucleotideAlignmentConstants.NUMBER_NUCLEOTIDE_ALLELES)) {
            throw new IllegalArgumentException("TBitAlignment: getNucleotideInstance: max number of alleles must be between 1 and 14 inclusive: " + maxNumAlleles);
        }

        if ((data == null) || (data.length == 0)) {
            throw new IllegalArgumentException("TBitAlignment: getNucleotideInstance: data can not be empty.");
        }

        if (data.length != idGroup.getIdCount()) {
            throw new IllegalArgumentException("TBitAlignment: getNucleotideInstance: data rows not equal to number of identifers.");
        }

        byte[][] dataBytes = AlignmentUtils.getDataBytes(data, NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES, NucleotideAlignmentConstants.NUMBER_NUCLEOTIDE_ALLELES);

        return TBitAlignment.getNucleotideInstance(idGroup, dataBytes, map, reference, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles);

    }

    public static TBitAlignment getNucleotideInstance(IdGroup idGroup, String[] data, GeneticMap map, byte[] reference, int[] variableSites, int maxNumAlleles, Locus[] loci, int[] lociOffsets, String[] snpIDs, boolean retainRareAlleles) {

        if ((maxNumAlleles < 1) || (maxNumAlleles > NucleotideAlignmentConstants.NUMBER_NUCLEOTIDE_ALLELES)) {
            throw new IllegalArgumentException("TBitAlignment: getNucleotideInstance: max number of alleles must be between 1 and 14 inclusive: " + maxNumAlleles);
        }

        if ((data == null) || (data.length == 0)) {
            throw new IllegalArgumentException("TBitAlignment: getNucleotideInstance: data can not be empty.");
        }

        if (data.length != idGroup.getIdCount()) {
            throw new IllegalArgumentException("TBitAlignment: getNucleotideInstance: data rows not equal to number of identifers.");
        }

        byte[][] dataBytes = AlignmentUtils.getDataBytes(data);

        return TBitAlignment.getNucleotideInstance(idGroup, dataBytes, map, reference, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles);

    }

    public static TBitAlignment getInstance(IdGroup idGroup, String[][] data, GeneticMap map, byte[] reference, int[] variableSites, int maxNumAlleles, Locus[] loci, int[] lociOffsets, String[] snpIDs, boolean retainRareAlleles) {

        if ((maxNumAlleles < 1) || (maxNumAlleles > 14)) {
            throw new IllegalArgumentException("TBitAlignment: getInstance: max number of alleles must be between 1 and 14 inclusive: " + maxNumAlleles);
        }

        if ((data == null) || (data.length == 0)) {
            throw new IllegalArgumentException("TBitAlignment: getInstance: data can not be empty.");
        }

        if (data.length != idGroup.getIdCount()) {
            throw new IllegalArgumentException("TBitAlignment: getInstance: data rows not equal to number of identifers.");
        }

        String[][] alleleStates = AlignmentUtils.getAlleleStates(data, maxNumAlleles);

        byte[][] dataBytes = AlignmentUtils.getDataBytes(data, alleleStates, maxNumAlleles);

        return TBitAlignment.getInstance(idGroup, dataBytes, map, reference, alleleStates, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles);

    }

    private void loadAlleles(byte[][] data) {

        myNumDataRows = myMaxNumAlleles;
        if (retainsRareAlleles()) {
            myNumDataRows++;
        }
        myData = new OpenBitSet[myNumDataRows][getSequenceCount()];
        for (int al = 0; al < myNumDataRows; al++) {
            for (int t = 0; t < getSequenceCount(); t++) {
                myData[al][t] = new OpenBitSet(myNumSites);
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
                                myData[j][t].fastSet(s);
                                isRare = false;
                                break;
                            }
                        }
                        if (isRare && retainsRareAlleles()) {
                            myData[myMaxNumAlleles][t].fastSet(s);
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
        myData = new OpenBitSet[myNumDataRows][getSequenceCount()];
        for (int al = 0; al < myNumDataRows; al++) {
            for (int t = 0; t < getSequenceCount(); t++) {
                myData[al][t] = new OpenBitSet(myNumSites);
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
                                myData[j][t].fastSet(s);
                                isRare = false;
                                break;
                            }
                        }
                        if (isRare && retainsRareAlleles()) {
                            myData[myMaxNumAlleles][t].fastSet(s);
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
                if (myData[i][taxon].fastGet(site)) {
                    if (count == 0) {
                        result[1] = myAlleles[site][i];
                    }
                    result[count++] = myAlleles[site][i];
                }
            }

            // Check For Rare Allele
            if (retainsRareAlleles() && myData[myMaxNumAlleles][taxon].fastGet(site)) {
                if (count == 0) {
                    result[1] = Alignment.RARE_ALLELE;
                }
                result[count] = Alignment.RARE_ALLELE;
            }

        } catch (IndexOutOfBoundsException e) {
            e.printStackTrace();
            throw new IllegalStateException("TBitAlignment: getBaseArray: bit sets indicate more than two alleles for taxon: " + taxon + "   site: " + site);
        }
        return result;
    }

    @Override
    public BitSet getAllelePresenceForAllSites(int taxon, int alleleNumber) {
        return UnmodifiableBitSet.getInstance(myData[alleleNumber][taxon]);
    }

    @Override
    public long[] getAllelePresenceForSitesBlock(int taxon, int alleleNumber, int startBlock, int endBlock) {
        long[] result = new long[endBlock - startBlock];
        System.arraycopy(myData[alleleNumber][taxon].getBits(), startBlock, result, 0, endBlock - startBlock);
        return result;
    }

    @Override
    public boolean isHeterozygous(int taxon, int site) {
        int count = 0;
        for (int i = 0; i < myNumDataRows; i++) {
            if (myData[i][taxon].fastGet(site)) {
                count++;
                if (count == 2) {
                    return true;
                }
            }
        }
        return false;
    }

    @Override
    public int getTotalGametesNotMissingForTaxon(int taxon) {

        OpenBitSet temp = new OpenBitSet(getSequenceCount());
        for (int i = 0; i < myNumDataRows; i++) {
            temp.or(myData[i][taxon]);
        }
        return ((int) temp.cardinality()) * 2;

    }

    @Override
    public int getHeterozygousCountForTaxon(int taxon) {

        int result = 0;
        for (int i = 0; i < myNumDataRows; i++) {
            for (int j = i + 1; j < myNumDataRows; j++) {
                result += (int) OpenBitSet.intersectionCount(myData[i][taxon], myData[j][taxon]);
            }
        }
        return result;

    }

    @Override
    public boolean isSBitFriendly() {
        return false;
    }

    @Override
    public boolean isTBitFriendly() {
        return true;
    }
}
