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

    protected SBitAlignment(Alignment a, int maxNumAlleles, boolean retainRareAlleles, boolean isFinalized) {
        super(a, maxNumAlleles, retainRareAlleles, isFinalized);
        loadAlleles(a);
    }

    protected SBitAlignment(IdGroup idGroup, byte[][] data, GeneticMap map, byte[] reference, String[][] alleleStates, int[] variableSites, int maxNumAlleles, Locus[] loci, int[] lociOffsets, String[] snpIDs, boolean retainRareAlleles, boolean isFinalized) {
        super(idGroup, data, map, reference, alleleStates, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles, isFinalized);
        loadAlleles(data);
    }

    public static SBitAlignment getInstance(Alignment a) {
        return SBitAlignment.getInstance(a, a.getMaxNumAlleles(), a.retainsRareAlleles(), true);
    }

    public static SBitAlignment getInstance(Alignment a, int maxNumAlleles, boolean retainRareAlleles, boolean isFinalized) {
        String[][] alleleStates = a.getAlleleEncodings();
        if ((alleleStates == null) || (alleleStates.length == 0)) {
            throw new IllegalArgumentException("SBitAlignment: init: allele states can not be empty.");
        }
        if (alleleStates.length == 1) {
            return new SBitAlignment(a, maxNumAlleles, retainRareAlleles, isFinalized);
        } else {
            return new SBitTextAlignment(a, maxNumAlleles, retainRareAlleles, isFinalized);
        }
    }

    public static SBitAlignment getInstance(IdGroup idGroup, byte[][] data, GeneticMap map, byte[] reference, String[][] alleleStates, int[] variableSites, int maxNumAlleles, Locus[] loci, int[] lociOffsets, String[] snpIDs, boolean retainRareAlleles, boolean isFinalized) {
        if ((alleleStates == null) || (alleleStates.length == 0)) {
            throw new IllegalArgumentException("SBitAlignment: init: allele states can not be empty.");
        }
        if (alleleStates.length == 1) {
            return new SBitAlignment(idGroup, data, map, reference, alleleStates, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles, isFinalized);
        } else {
            return new SBitTextAlignment(idGroup, data, map, reference, alleleStates, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles, isFinalized);
        }
    }

    public static SBitAlignment getNucleotideInstance(IdGroup idGroup, byte[][] data, GeneticMap map, byte[] reference, int[] variableSites, int maxNumAlleles, Locus[] loci, int[] lociOffsets, String[] snpIDs, boolean retainRareAlleles, boolean isFinalized) {
        return new SBitNucleotideAlignment(idGroup, data, map, reference, NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles, isFinalized);
    }

    public static SBitAlignment getNucleotideInstance(IdGroup idGroup, String[][] data, GeneticMap map, byte[] reference, int[] variableSites, int maxNumAlleles, Locus[] loci, int[] lociOffsets, String[] snpIDs, boolean retainRareAlleles, boolean isFinalized) {

        if ((maxNumAlleles < 1) || (maxNumAlleles > NucleotideAlignmentConstants.NUMBER_NUCLEOTIDE_ALLELES)) {
            throw new IllegalArgumentException("SBitAlignment: getInstance: max number of alleles must be between 1 and 14 inclusive: " + maxNumAlleles);
        }

        if ((data == null) || (data.length == 0)) {
            throw new IllegalArgumentException("SBitAlignment: getInstance: data can not be empty.");
        }

        if (data.length != idGroup.getIdCount()) {
            throw new IllegalArgumentException("SBitAlignment: getInstance: data rows not equal to number of identifers.");
        }

        byte[][] dataBytes = AlignmentUtils.getDataBytes(data, NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES, NucleotideAlignmentConstants.NUMBER_NUCLEOTIDE_ALLELES);

        return SBitAlignment.getNucleotideInstance(idGroup, dataBytes, map, reference, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles, isFinalized);

    }

    public static SBitAlignment getInstance(IdGroup idGroup, String[][] data, GeneticMap map, byte[] reference, int[] variableSites, int maxNumAlleles, Locus[] loci, int[] lociOffsets, String[] snpIDs, boolean retainRareAlleles, boolean isFinalized) {

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

        return SBitAlignment.getInstance(idGroup, dataBytes, map, reference, alleleStates, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles, isFinalized);

    }

    private void loadAlleles(byte[][] data) {

        myNumDataRows = myMaxNumAlleles + 1;
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
                    if (cb[i] == Alignment.UNKNOWN_ALLELE) {
                        myData[myMaxNumAlleles][s].fastSet(t);
                    } else {
                        boolean isRare = true;
                        for (int j = 0; j < myMaxNumAlleles; j++) {
                            if (cb[i] == myAlleles[s][j]) {
                                myData[j][s].fastSet(t);
                                isRare = false;
                                break;
                            }
                        }
                        if (isRare) {
                            if (retainsRareAlleles()) {
                                // Record as Rare Allele
                                myData[myMaxNumAlleles + 1][s].fastSet(t);
                            } else {
                                // Change to Unknown
                                myData[myMaxNumAlleles][s].fastSet(t);
                            }
                        }
                    }
                }
            }
        }

    }

    private void loadAlleles(Alignment a) {

        myNumDataRows = myMaxNumAlleles + 1;
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
                    if (cb[i] == Alignment.UNKNOWN_ALLELE) {
                        myData[myMaxNumAlleles][s].fastSet(t);
                    } else {
                        boolean isRare = true;
                        for (int j = 0; j < myMaxNumAlleles; j++) {
                            if (cb[i] == myAlleles[s][j]) {
                                myData[j][s].fastSet(t);
                                isRare = false;
                                break;
                            }
                        }
                        if (isRare) {
                            if (retainsRareAlleles()) {
                                // Record as Rare Allele
                                myData[myMaxNumAlleles + 1][s].fastSet(t);
                            } else {
                                // Change to Unknown
                                myData[myMaxNumAlleles][s].fastSet(t);
                            }
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
            if (retainsRareAlleles() && myData[myMaxNumAlleles + 1][site].fastGet(taxon)) {
                if (count == 0) {
                    result[1] = Alignment.RARE_ALLELE;
                }
                result[count] = Alignment.RARE_ALLELE;
            }

            // Check For Unknown
            if (myData[myMaxNumAlleles][site].fastGet(taxon)) {
                result[1] = Alignment.UNKNOWN_ALLELE;
            }

        } catch (IndexOutOfBoundsException e) {
            throw new IllegalStateException("SBitAlignment: getBaseArray: bit sets indicate more than two alleles for taxon: " + taxon + "   site: " + site);
        }
        return result;
    }

    @Override
    public BitSet getAllelePresenceForAllTaxa(int site, int alleleNumber) {
        if (isFinalized()) {
            return UnmodifiableBitSet.getInstance(myData[alleleNumber][site]);
        } else {
            return myData[alleleNumber][site];
        }
    }

    //
    // MutableAlignment Methods...
    //
    @Override
    public void setBase(int taxon, int site, byte newBase) {

        if (isFinalized()) {
            throw new IllegalStateException("SBitAlignment: This Alignment has been finalized and can not be changed.");
        }

        for (int x = 0; x < myNumDataRows; x++) {
            myData[x][site].fastClear(taxon);
        }

        byte[] cb = new byte[2];
        cb[0] = (byte) ((newBase >>> 4) & 0xf);
        cb[1] = (byte) (newBase & 0xf);
        for (int i = 0; i < 2; i++) {
            if (cb[i] == Alignment.UNKNOWN_ALLELE) {
                myData[myMaxNumAlleles][site].fastSet(taxon);
            } else {
                boolean isRare = true;
                for (int j = 0; j < myMaxNumAlleles; j++) {
                    if (cb[i] == myAlleles[site][j]) {
                        myData[j][site].fastSet(taxon);
                        isRare = false;
                        break;
                    } else if (myAlleles[site][j] == Alignment.UNKNOWN_ALLELE) {
                        myAlleles[site][j] = cb[i];
                        myData[j][site].fastSet(taxon);
                        isRare = false;
                        break;
                    }
                }
                if (isRare) {
                    if (retainsRareAlleles()) {
                        myData[myMaxNumAlleles + 1][site].fastSet(taxon);
                    } else {
                        myData[myMaxNumAlleles][site].fastSet(taxon);
                    }
                }
            }
        }

        myIsDirty = true;
    }

    @Override
    public void setBaseRange(int taxon, int startSite, byte[] newBases) {

        if (isFinalized()) {
            throw new IllegalStateException("SBitAlignment: This Alignment has been finalized and can not be changed.");
        }

        for (int i = 0, n = newBases.length; i < n; i++) {
            setBase(taxon, i + startSite, newBases[i]);
        }

    }

    @Override
    public void addSite(int site) {

        if (isFinalized()) {
            throw new IllegalStateException("SBitAlignment: This Alignment has been finalized and can not be changed.");
        }

        if (myAlleleStates.length != 1) {
            throw new IllegalStateException("SBitAlignment: addSite: Currently no way to add site to Alignment that has multiple allele state mappings.");
        }

        OpenBitSet[][] temp = myData;
        myNumSites++;
        myData = new OpenBitSet[myNumDataRows][myNumSites];
        int offset = 0;
        for (int a = 0; a < myNumDataRows; a++) {
            for (int s = 0; s < myNumSites; s++) {
                if (site == s) {
                    myData[a][s] = new OpenBitSet(getSequenceCount());
                    offset = 1;
                } else {
                    myData[a][s] = temp[a][s - offset];
                }

            }
        }

        byte[][] tempAlleles = myAlleles;
        myAlleles = new byte[myNumSites][];
        offset = 0;
        for (int i = 0; i < myNumSites; i++) {
            if (site == i) {
                for (int j = 0; j < myMaxNumAlleles; j++) {
                    myAlleles[i][j] = Alignment.UNKNOWN_ALLELE;
                }
                offset = 1;
            } else {
                myAlleles[i] = tempAlleles[i - offset];
            }
        }

        myIsDirty = true;

    }

    @Override
    public void removeSite(int site) {

        if (isFinalized()) {
            throw new IllegalStateException("SBitAlignment: This Alignment has been finalized and can not be changed.");
        }

        OpenBitSet[][] temp = myData;
        myNumSites--;
        myData = new OpenBitSet[myNumDataRows][myNumSites];
        int offset = 0;
        for (int a = 0; a < myNumDataRows; a++) {
            for (int s = 0; s < myNumSites; s++) {
                if (site == s) {
                    offset = 1;
                    myData[a][s] = temp[a][s + offset];
                } else {
                    myData[a][s] = temp[a][s + offset];
                }

            }
        }

        byte[][] tempAlleles = myAlleles;
        myAlleles = new byte[myNumSites][];
        offset = 0;
        for (int i = 0; i < myNumSites; i++) {
            if (site == i) {
                offset = 1;
                myAlleles[i] = tempAlleles[i + offset];
            } else {
                myAlleles[i] = tempAlleles[i + offset];
            }
        }

        if (myAlleleStates.length != 1) {
            String[][] tempAlleleStates = myAlleleStates;
            myAlleleStates = new String[myNumSites][];
            offset = 0;
            for (int i = 0; i < myNumSites; i++) {
                if (site == i) {
                    offset = 1;
                    myAlleleStates[i] = tempAlleleStates[i + offset];
                } else {
                    myAlleleStates[i] = tempAlleleStates[i + offset];
                }
            }
        }

        myIsDirty = true;

    }

    @Override
    public void clean() {
        if (isFinalized()) {
            throw new IllegalStateException("SBitAlignment: This Alignment has been finalized and can not be changed.");
        }
        // Nothing needs to be done?
        myIsDirty = false;
    }

    @Override
    public boolean isDirty() {
        return myIsDirty;
    }
}
