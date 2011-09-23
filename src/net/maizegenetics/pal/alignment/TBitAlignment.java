/*
 * TBitAlignment
 */
package net.maizegenetics.pal.alignment;

import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.Identifier;
import net.maizegenetics.pal.ids.SimpleIdGroup;
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
    private boolean myIsDirty = false;

    protected TBitAlignment(Alignment a, int maxNumAlleles, boolean retainRareAlleles, boolean isFinalized) {
        super(a, maxNumAlleles, retainRareAlleles, isFinalized);
        loadAlleles(a);
    }

    protected TBitAlignment(IdGroup idGroup, byte[][] data, GeneticMap map, byte[] reference, String[][] alleleStates, int[] variableSites, int maxNumAlleles, Locus[] loci, int[] lociOffsets, String[] snpIDs, boolean retainRareAlleles, boolean isFinalized) {
        super(idGroup, data, map, reference, alleleStates, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles, isFinalized);
        loadAlleles(data);
    }

    public static TBitAlignment getInstance(Alignment a) {
        return TBitAlignment.getInstance(a, a.getMaxNumAlleles(), a.retainsRareAlleles(), true);
    }

    public static TBitAlignment getInstance(Alignment a, int maxNumAlleles, boolean retainRareAlleles, boolean isFinalized) {
        String[][] alleleStates = a.getAlleleEncodings();
        if ((alleleStates == null) || (alleleStates.length == 0)) {
            throw new IllegalArgumentException("TBitAlignment: init: allele states can not be empty.");
        }
        if (alleleStates.length == 1) {
            return new TBitAlignment(a, maxNumAlleles, retainRareAlleles, isFinalized);
        } else {
            return new TBitTextAlignment(a, maxNumAlleles, retainRareAlleles, isFinalized);
        }
    }

    public static TBitAlignment getInstance(IdGroup idGroup, byte[][] data, GeneticMap map, byte[] reference, String[][] alleleStates, int[] variableSites, int maxNumAlleles, Locus[] loci, int[] lociOffsets, String[] snpIDs, boolean retainRareAlleles, boolean isFinalized) {
        if ((alleleStates == null) || (alleleStates.length == 0)) {
            throw new IllegalArgumentException("TBitAlignment: init: allele states can not be empty.");
        }
        if (alleleStates.length == 1) {
            return new TBitAlignment(idGroup, data, map, reference, alleleStates, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles, isFinalized);
        } else {
            return new TBitTextAlignment(idGroup, data, map, reference, alleleStates, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles, isFinalized);
        }
    }

    public static TBitAlignment getNucleotideInstance(IdGroup idGroup, String[][] data, GeneticMap map, byte[] reference, int[] variableSites, int maxNumAlleles, Locus[] loci, int[] lociOffsets, String[] snpIDs, boolean retainRareAlleles, boolean isFinalized) {

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

        return TBitAlignment.getNucleotideInstance(idGroup, dataBytes, map, reference, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles, isFinalized);

    }

    public static TBitAlignment getNucleotideInstance(IdGroup idGroup, String[] data, GeneticMap map, byte[] reference, int[] variableSites, int maxNumAlleles, Locus[] loci, int[] lociOffsets, String[] snpIDs, boolean retainRareAlleles, boolean isFinalized) {

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

        return TBitAlignment.getNucleotideInstance(idGroup, dataBytes, map, reference, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles, isFinalized);

    }

    public static TBitAlignment getNucleotideInstance(IdGroup idGroup, byte[][] data, GeneticMap map, byte[] reference, int[] variableSites, int maxNumAlleles, Locus[] loci, int[] lociOffsets, String[] snpIDs, boolean retainRareAlleles, boolean isFinalized) {
        return new TBitNucleotideAlignment(idGroup, data, map, reference, NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles, isFinalized);
    }

    public static TBitAlignment getInstance(IdGroup idGroup, String[][] data, GeneticMap map, byte[] reference, int[] variableSites, int maxNumAlleles, Locus[] loci, int[] lociOffsets, String[] snpIDs, boolean retainRareAlleles, boolean isFinalized) {

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

        return TBitAlignment.getInstance(idGroup, dataBytes, map, reference, alleleStates, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles, isFinalized);

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
                    if (cb[i] == Alignment.UNKNOWN_ALLELE) {
                        myData[myMaxNumAlleles][t].fastSet(s);
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
    public BitSet getAllelePresensceForAllSites(int taxon, int alleleNumber) {
        if (isFinalized()) {
            return UnmodifiableBitSet.getInstance(myData[alleleNumber][taxon]);
        } else {
            return myData[alleleNumber][taxon];
        }
    }

    @Override
    public long[] getAllelePresensceForSitesBlock(int taxon, int alleleNumber, int startBlock, int endBlock) {
        long[] result = new long[endBlock - startBlock];
        System.arraycopy(myData[alleleNumber][taxon].getBits(), startBlock, result, 0, endBlock - startBlock);
        return result;
    }

    //
    // MutableAlignment Methods...
    //
    @Override
    public void setBase(int taxon, int site, byte newBase) {

        if (isFinalized()) {
            throw new IllegalStateException("TBitAlignment: This Alignment has been finalized and can not be changed.");
        }

        for (int x = 0; x < myNumDataRows; x++) {
            myData[x][taxon].fastClear(site);
        }

        byte[] cb = new byte[2];
        cb[0] = (byte) ((newBase >>> 4) & 0xf);
        cb[1] = (byte) (newBase & 0xf);
        for (int i = 0; i < 2; i++) {
            if (cb[i] == Alignment.UNKNOWN_ALLELE) {
                myData[myMaxNumAlleles][taxon].fastSet(site);
            } else {
                boolean isRare = true;
                for (int j = 0; j < myMaxNumAlleles; j++) {
                    if (cb[i] == myAlleles[site][j]) {
                        myData[j][taxon].fastSet(site);
                        isRare = false;
                        break;
                    } else if (myAlleles[site][j] == Alignment.UNKNOWN_ALLELE) {
                        myAlleles[site][j] = cb[i];
                        myData[j][taxon].fastSet(site);
                        isRare = false;
                        break;
                    }
                }
                if (isRare) {
                    if (retainsRareAlleles()) {
                        myData[myMaxNumAlleles + 1][taxon].fastSet(site);
                    } else {
                        myData[myMaxNumAlleles][taxon].fastSet(site);
                    }
                }
            }
        }

        myIsDirty = true;

    }

    @Override
    public void setBaseRange(int taxon, int startSite, byte[] newBases) {

        if (isFinalized()) {
            throw new IllegalStateException("TBitAlignment: This Alignment has been finalized and can not be changed.");
        }

        for (int i = 0, n = newBases.length; i < n; i++) {
            setBase(taxon, i + startSite, newBases[i]);
        }

    }

    @Override
    public void addTaxon(Identifier id) {

        if (isFinalized()) {
            throw new IllegalStateException("TBitAlignment: This Alignment has been finalized and can not be changed.");
        }

        Identifier[] ids = new Identifier[myIdGroup.getIdCount() + 1];
        for (int i = 0; i < myIdGroup.getIdCount(); i++) {
            ids[i] = myIdGroup.getIdentifier(i);
        }
        ids[myIdGroup.getIdCount()] = id;
        myIdGroup = new SimpleIdGroup(ids);

        OpenBitSet[][] temp = myData;
        myData = new OpenBitSet[myNumDataRows][getSequenceCount()];
        for (int a = 0; a < myNumDataRows; a++) {
            for (int t = 0; t < getSequenceCount() - 1; t++) {
                myData[a][t] = temp[a][t];
            }
            myData[a][getSequenceCount()] = new OpenBitSet(myNumSites);
        }

        myIsDirty = true;

    }

    @Override
    public void removeTaxon(int taxon) {

        if (isFinalized()) {
            throw new IllegalStateException("TBitAlignment: This Alignment has been finalized and can not be changed.");
        }

        int numIds = myIdGroup.getIdCount() - 1;
        Identifier[] ids = new Identifier[numIds];
        int offset = 0;
        for (int i = 0; i < numIds; i++) {
            if (taxon == i) {
                offset = 1;
                ids[i] = myIdGroup.getIdentifier(i + offset);
            } else {
                ids[i] = myIdGroup.getIdentifier(i + offset);
            }
        }
        myIdGroup = new SimpleIdGroup(ids);

        OpenBitSet[][] temp = myData;
        myData = new OpenBitSet[myNumDataRows][getSequenceCount()];
        offset = 0;
        for (int a = 0; a < myNumDataRows; a++) {
            for (int t = 0; t < getSequenceCount(); t++) {
                if (taxon == t) {
                    offset = 1;
                    myData[a][t] = temp[a][t + offset];
                } else {
                    myData[a][t] = temp[a][t + offset];
                }

            }
        }

        myIsDirty = true;

    }

    @Override
    public void clean() {
        if (isFinalized()) {
            throw new IllegalStateException("TBitAlignment: This Alignment has been finalized and can not be changed.");
        }
        // Nothing needs to be done?
        myIsDirty = false;
    }

    @Override
    public boolean isDirty() {
        return myIsDirty;
    }
}
