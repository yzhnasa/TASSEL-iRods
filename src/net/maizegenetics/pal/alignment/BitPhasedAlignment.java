/*
 * BitPhasedAlignment
 */
package net.maizegenetics.pal.alignment;

import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;
import net.maizegenetics.util.UnmodifiableBitSet;

/**
 * This data alignment is optimized for operations involving lots of SNPs from a
 * taxon - imputation, genetic distance, kinship, diversity, etc. It is not
 * optimized for LD or association mapping.
 *
 * @author terry
 */
public class BitPhasedAlignment extends AbstractAlignment {

    private OpenBitSet[][] mySBitData0;
    private OpenBitSet[][] mySBitData1;
    private OpenBitSet[][] myTBitData0;
    private OpenBitSet[][] myTBitData1;
    private int myNumDataRows;

    protected BitPhasedAlignment(Alignment a, int maxNumAlleles, boolean retainRareAlleles) {
        super(a, maxNumAlleles, retainRareAlleles);
        loadAlleles(a);
    }

    protected BitPhasedAlignment(IdGroup idGroup, byte[][] data, GeneticMap map, byte[] reference, String[][] alleleStates, int[] variableSites, int maxNumAlleles, Locus[] loci, int[] lociOffsets, String[] snpIDs, boolean retainRareAlleles) {
        super(idGroup, data, map, reference, alleleStates, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles);
        loadAlleles(data);
    }

    public static BitPhasedAlignment getInstance(Alignment a) {
        return BitPhasedAlignment.getInstance(a, a.getMaxNumAlleles(), a.retainsRareAlleles());
    }

    public static BitPhasedAlignment getInstance(Alignment a, int maxNumAlleles, boolean retainRareAlleles) {
        String[][] alleleStates = a.getAlleleEncodings();
        if ((alleleStates == null) || (alleleStates.length == 0)) {
            throw new IllegalArgumentException("SBitPhasedAlignment: init: allele states can not be empty.");
        }
        if (alleleStates.length == 1) {
            return new BitPhasedAlignment(a, maxNumAlleles, retainRareAlleles);
        } else {
            return new BitPhasedTextAlignment(a, maxNumAlleles, retainRareAlleles);
        }
    }

    public static BitPhasedAlignment getInstance(IdGroup idGroup, byte[][] data, GeneticMap map, byte[] reference, String[][] alleleStates, int[] variableSites, int maxNumAlleles, Locus[] loci, int[] lociOffsets, String[] snpIDs, boolean retainRareAlleles) {
        if ((alleleStates == null) || (alleleStates.length == 0)) {
            throw new IllegalArgumentException("SBitPhasedAlignment: init: allele states can not be empty.");
        }
        if (alleleStates.length == 1) {
            return new BitPhasedAlignment(idGroup, data, map, reference, alleleStates, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles);
        } else {
            return new BitPhasedTextAlignment(idGroup, data, map, reference, alleleStates, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles);
        }
    }

    public static BitPhasedAlignment getNucleotideInstance(IdGroup idGroup, byte[][] data, GeneticMap map, byte[] reference, int[] variableSites, int maxNumAlleles, Locus[] loci, int[] lociOffsets, String[] snpIDs, boolean retainRareAlleles) {
        return new BitPhasedNucleotideAlignment(idGroup, data, map, reference, NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles);
    }

    public static BitPhasedAlignment getNucleotideInstance(IdGroup idGroup, String[][] data, GeneticMap map, byte[] reference, int[] variableSites, int maxNumAlleles, Locus[] loci, int[] lociOffsets, String[] snpIDs, boolean retainRareAlleles) {

        if ((maxNumAlleles < 1) || (maxNumAlleles > NucleotideAlignmentConstants.NUMBER_NUCLEOTIDE_ALLELES)) {
            throw new IllegalArgumentException("SBitPhasedAlignment: getNucleotideInstance: max number of alleles must be between 1 and 14 inclusive: " + maxNumAlleles);
        }

        if ((data == null) || (data.length == 0)) {
            throw new IllegalArgumentException("SBitPhasedAlignment: getNucleotideInstance: data can not be empty.");
        }

        if (data.length != idGroup.getIdCount()) {
            throw new IllegalArgumentException("SBitPhasedAlignment: getNucleotideInstance: data rows not equal to number of identifers.");
        }

        byte[][] dataBytes = AlignmentUtils.getDataBytes(data, NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES, NucleotideAlignmentConstants.NUMBER_NUCLEOTIDE_ALLELES);

        return BitPhasedAlignment.getNucleotideInstance(idGroup, dataBytes, map, reference, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles);

    }

    public static BitPhasedAlignment getNucleotideInstance(IdGroup idGroup, String[] data, GeneticMap map, byte[] reference, int[] variableSites, int maxNumAlleles, Locus[] loci, int[] lociOffsets, String[] snpIDs, boolean retainRareAlleles) {

        if ((maxNumAlleles < 1) || (maxNumAlleles > NucleotideAlignmentConstants.NUMBER_NUCLEOTIDE_ALLELES)) {
            throw new IllegalArgumentException("SBitPhasedAlignment: getNucleotideInstance: max number of alleles must be between 1 and 14 inclusive: " + maxNumAlleles);
        }

        if ((data == null) || (data.length == 0)) {
            throw new IllegalArgumentException("SBitPhasedAlignment: getNucleotideInstance: data can not be empty.");
        }

        if (data.length != idGroup.getIdCount()) {
            throw new IllegalArgumentException("SBitPhasedAlignment: getNucleotideInstance: data rows not equal to number of identifers.");
        }

        byte[][] dataBytes = AlignmentUtils.getDataBytes(data);

        return BitPhasedAlignment.getNucleotideInstance(idGroup, dataBytes, map, reference, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles);

    }

    public static BitPhasedAlignment getInstance(IdGroup idGroup, String[][] data, GeneticMap map, byte[] reference, int[] variableSites, int maxNumAlleles, Locus[] loci, int[] lociOffsets, String[] snpIDs, boolean retainRareAlleles) {

        if ((maxNumAlleles < 1) || (maxNumAlleles > 14)) {
            throw new IllegalArgumentException("SBitPhasedAlignment: getInstance: max number of alleles must be between 1 and 14 inclusive: " + maxNumAlleles);
        }

        if ((data == null) || (data.length == 0)) {
            throw new IllegalArgumentException("SBitPhasedAlignment: getInstance: data can not be empty.");
        }

        if (data.length != idGroup.getIdCount()) {
            throw new IllegalArgumentException("SBitPhasedAlignment: getInstance: data rows not equal to number of identifers.");
        }

        String[][] alleleStates = AlignmentUtils.getAlleleStates(data, maxNumAlleles);

        byte[][] dataBytes = AlignmentUtils.getDataBytes(data, alleleStates, maxNumAlleles);

        return BitPhasedAlignment.getInstance(idGroup, dataBytes, map, reference, alleleStates, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles);

    }

    private void loadAlleles(byte[][] data) {

        myNumDataRows = myMaxNumAlleles;
        if (retainsRareAlleles()) {
            myNumDataRows++;
        }
        mySBitData0 = new OpenBitSet[myNumDataRows][myNumSites];
        mySBitData1 = new OpenBitSet[myNumDataRows][myNumSites];
        for (int al = 0; al < myNumDataRows; al++) {
            for (int s = 0; s < myNumSites; s++) {
                mySBitData0[al][s] = new OpenBitSet(getSequenceCount());
                mySBitData1[al][s] = new OpenBitSet(getSequenceCount());
            }
        }
        byte[] cb = new byte[2];
        for (int s = 0; s < myNumSites; s++) {
            for (int t = 0, n = getSequenceCount(); t < n; t++) {
                cb[0] = (byte) ((data[t][s] >>> 4) & 0xf);
                cb[1] = (byte) (data[t][s] & 0xf);

                if (cb[0] != Alignment.UNKNOWN_ALLELE) {
                    boolean isRare = true;
                    for (int j = 0; j < myMaxNumAlleles; j++) {
                        if (cb[0] == myAlleles[s][j]) {
                            mySBitData0[j][s].fastSet(t);
                            isRare = false;
                            break;
                        }
                    }
                    if (isRare && retainsRareAlleles()) {
                        mySBitData0[myMaxNumAlleles][s].fastSet(t);
                    }
                }

                if (cb[1] != Alignment.UNKNOWN_ALLELE) {
                    boolean isRare = true;
                    for (int j = 0; j < myMaxNumAlleles; j++) {
                        if (cb[1] == myAlleles[s][j]) {
                            mySBitData1[j][s].fastSet(t);
                            isRare = false;
                            break;
                        }
                    }
                    if (isRare && retainsRareAlleles()) {
                        mySBitData1[myMaxNumAlleles][s].fastSet(t);
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
        mySBitData0 = new OpenBitSet[myNumDataRows][myNumSites];
        mySBitData1 = new OpenBitSet[myNumDataRows][myNumSites];
        for (int al = 0; al < myNumDataRows; al++) {
            for (int s = 0; s < myNumSites; s++) {
                mySBitData0[al][s] = new OpenBitSet(getSequenceCount());
                mySBitData1[al][s] = new OpenBitSet(getSequenceCount());
            }
        }
        for (int s = 0; s < myNumSites; s++) {
            for (int t = 0, n = getSequenceCount(); t < n; t++) {
                byte[] cb = a.getBaseArray(t, s);

                if (cb[0] != Alignment.UNKNOWN_ALLELE) {
                    boolean isRare = true;
                    for (int j = 0; j < myMaxNumAlleles; j++) {
                        if (cb[0] == myAlleles[s][j]) {
                            mySBitData0[j][s].fastSet(t);
                            isRare = false;
                            break;
                        }
                    }
                    if (isRare && retainsRareAlleles()) {
                        mySBitData0[myMaxNumAlleles][s].fastSet(t);
                    }
                }

                if (cb[1] != Alignment.UNKNOWN_ALLELE) {
                    boolean isRare = true;
                    for (int j = 0; j < myMaxNumAlleles; j++) {
                        if (cb[1] == myAlleles[s][j]) {
                            mySBitData1[j][s].fastSet(t);
                            isRare = false;
                            break;
                        }
                    }
                    if (isRare && retainsRareAlleles()) {
                        mySBitData1[myMaxNumAlleles][s].fastSet(t);
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

            for (int i = 0; i < myMaxNumAlleles; i++) {
                if (mySBitData0[i][site].fastGet(taxon)) {
                    result[0] = myAlleles[site][i];
                    break;
                }
            }

            for (int i = 0; i < myMaxNumAlleles; i++) {
                if (mySBitData1[i][site].fastGet(taxon)) {
                    result[1] = myAlleles[site][i];
                    break;
                }
            }

            // Check For Rare Allele
            if (retainsRareAlleles() && mySBitData0[myMaxNumAlleles][site].fastGet(taxon)) {
                result[0] = Alignment.RARE_ALLELE;
            }
            if (retainsRareAlleles() && mySBitData1[myMaxNumAlleles][site].fastGet(taxon)) {
                result[1] = Alignment.RARE_ALLELE;
            }

        } catch (IndexOutOfBoundsException e) {
            throw new IllegalStateException("SBitPhasedAlignment: getBaseArray: bit sets indicate more than two alleles for taxon: " + taxon + "   site: " + site);
        }
        return result;
    }

    @Override
    public boolean isPhased() {
        return true;
    }

    @Override
    public boolean isSBitFriendly() {
        return true;
    }

    @Override
    public boolean isTBitFriendly() {
        return false;
    }

    @Override
    public int getTotalNumAlleles() {
        return myNumDataRows;
    }

    @Override
    public BitSet getPhasedAllelePresenceForAllSites(int taxon, boolean firstParent, int alleleNumber) {
        if (myTBitData0 != null) {
            if (firstParent) {
                return UnmodifiableBitSet.getInstance(myTBitData0[alleleNumber][taxon]);
            } else {
                return UnmodifiableBitSet.getInstance(myTBitData1[alleleNumber][taxon]);
            }
        } else {
            throw new IllegalStateException("BitPhasedAlignment: getPhasedAllelePresenceForAllSites: This alignment hasn't been optimized for Taxa Operations.");
        }
    }

    @Override
    public BitSet getPhasedAllelePresenceForAllTaxa(int site, boolean firstParent, int alleleNumber) {
        if (mySBitData0 != null) {
            if (firstParent) {
                return UnmodifiableBitSet.getInstance(mySBitData0[alleleNumber][site]);
            } else {
                return UnmodifiableBitSet.getInstance(mySBitData1[alleleNumber][site]);
            }
        } else {
            throw new IllegalStateException("BitPhasedAlignment: getPhasedAllelePresenceForAllTaxa: This alignment hasn't been optimized for Site Operations.");
        }
    }

    @Override
    public long[] getPhasedAllelePresenceForSitesBlock(int taxon, boolean firstParent, int alleleNumber, int startBlock, int endBlock) {
        if (myTBitData0 != null) {
            BitSet temp = getAllelePresenceForAllSites(taxon, alleleNumber);
            long[] result = new long[endBlock - startBlock];
            System.arraycopy(temp.getBits(), startBlock, result, 0, endBlock - startBlock);
            return result;
        } else {
            throw new IllegalStateException("BitPhasedAlignment: getPhasedAllelePresenceForSitesBlock: This alignment hasn't been optimized for Taxa Operations.");
        }
    }

    @Override
    public BitSet getAllelePresenceForAllTaxa(int site, int alleleNumber) {
        if (mySBitData0 != null) {
            OpenBitSet temp = new OpenBitSet(getSequenceCount());
            temp.or(mySBitData0[alleleNumber][site]);
            temp.or(mySBitData1[alleleNumber][site]);
            return temp;
        } else {
            throw new IllegalStateException("BitPhasedAlignment: getAllelePresenceForAllTaxa: This alignment hasn't been optimized for Site Operations.");
        }
    }

    @Override
    public BitSet getAllelePresenceForAllSites(int taxon, int alleleNumber) {
        if (myTBitData0 != null) {
            OpenBitSet temp = new OpenBitSet(getSequenceCount());
            temp.or(myTBitData0[alleleNumber][taxon]);
            temp.or(myTBitData1[alleleNumber][taxon]);
            return temp;
        } else {
            throw new IllegalStateException("BitPhasedAlignment: getAllelePresenceForAllSites: This alignment hasn't been optimized for Taxa Operations.");
        }
    }

    @Override
    public long[] getAllelePresenceForSitesBlock(int taxon, int alleleNumber, int startBlock, int endBlock) {
        if (myTBitData0 != null) {
            BitSet temp = getAllelePresenceForAllSites(taxon, alleleNumber);
            long[] result = new long[endBlock - startBlock];
            System.arraycopy(temp.getBits(), startBlock, result, 0, endBlock - startBlock);
            return result;
        } else {
            throw new IllegalStateException("BitPhasedAlignment: getAllelePresenceForSitesBlock: This alignment hasn't been optimized for Taxa Operations.");
        }
    }
}
