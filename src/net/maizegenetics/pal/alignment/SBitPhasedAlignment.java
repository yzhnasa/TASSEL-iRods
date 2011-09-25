/*
 * SBitPhasedAlignment
 */
package net.maizegenetics.pal.alignment;

import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.util.OpenBitSet;

/**
 * This data alignment is optimized for operations
 * involving lots of SNPs from a taxon - imputation, genetic distance, kinship, diversity, etc.
 * It is not optimized for LD or association mapping.
 *
 * @author terry
 */
public class SBitPhasedAlignment extends AbstractAlignment {

    private OpenBitSet[][] myData0;
    private OpenBitSet[][] myData1;
    private int myNumDataRows;
    private boolean myIsDirty = false;

    protected SBitPhasedAlignment(Alignment a, int maxNumAlleles, boolean retainRareAlleles) {
        super(a, maxNumAlleles, retainRareAlleles);
        loadAlleles(a);
    }

    protected SBitPhasedAlignment(IdGroup idGroup, byte[][] data, GeneticMap map, byte[] reference, String[][] alleleStates, int[] variableSites, int maxNumAlleles, Locus[] loci, int[] lociOffsets, String[] snpIDs, boolean retainRareAlleles) {
        super(idGroup, data, map, reference, alleleStates, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles);
        loadAlleles(data);
    }

    public static SBitPhasedAlignment getInstance(Alignment a) {
        return SBitPhasedAlignment.getInstance(a, a.getMaxNumAlleles(), a.retainsRareAlleles());
    }

    public static SBitPhasedAlignment getInstance(Alignment a, int maxNumAlleles, boolean retainRareAlleles) {
        String[][] alleleStates = a.getAlleleEncodings();
        if ((alleleStates == null) || (alleleStates.length == 0)) {
            throw new IllegalArgumentException("SBitPhasedAlignment: init: allele states can not be empty.");
        }
        if (alleleStates.length == 1) {
            return new SBitPhasedAlignment(a, maxNumAlleles, retainRareAlleles);
        } else {
            return new SBitPhasedTextAlignment(a, maxNumAlleles, retainRareAlleles);
        }
    }

    public static SBitPhasedAlignment getInstance(IdGroup idGroup, byte[][] data, GeneticMap map, byte[] reference, String[][] alleleStates, int[] variableSites, int maxNumAlleles, Locus[] loci, int[] lociOffsets, String[] snpIDs, boolean retainRareAlleles) {
        if ((alleleStates == null) || (alleleStates.length == 0)) {
            throw new IllegalArgumentException("SBitPhasedAlignment: init: allele states can not be empty.");
        }
        if (alleleStates.length == 1) {
            return new SBitPhasedAlignment(idGroup, data, map, reference, alleleStates, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles);
        } else {
            return new SBitPhasedTextAlignment(idGroup, data, map, reference, alleleStates, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles);
        }
    }

    public static SBitPhasedAlignment getNucleotideInstance(IdGroup idGroup, byte[][] data, GeneticMap map, byte[] reference, int[] variableSites, int maxNumAlleles, Locus[] loci, int[] lociOffsets, String[] snpIDs, boolean retainRareAlleles) {
        return new SBitPhasedNucleotideAlignment(idGroup, data, map, reference, NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles);
    }

    public static SBitPhasedAlignment getNucleotideInstance(IdGroup idGroup, String[][] data, GeneticMap map, byte[] reference, int[] variableSites, int maxNumAlleles, Locus[] loci, int[] lociOffsets, String[] snpIDs, boolean retainRareAlleles) {

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

        return SBitPhasedAlignment.getNucleotideInstance(idGroup, dataBytes, map, reference, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles);

    }

    public static SBitPhasedAlignment getNucleotideInstance(IdGroup idGroup, String[] data, GeneticMap map, byte[] reference, int[] variableSites, int maxNumAlleles, Locus[] loci, int[] lociOffsets, String[] snpIDs, boolean retainRareAlleles) {

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

        return SBitPhasedAlignment.getNucleotideInstance(idGroup, dataBytes, map, reference, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles);

    }

    public static SBitPhasedAlignment getInstance(IdGroup idGroup, String[][] data, GeneticMap map, byte[] reference, int[] variableSites, int maxNumAlleles, Locus[] loci, int[] lociOffsets, String[] snpIDs, boolean retainRareAlleles) {

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

        return SBitPhasedAlignment.getInstance(idGroup, dataBytes, map, reference, alleleStates, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles);

    }

    private void loadAlleles(byte[][] data) {

        myNumDataRows = myMaxNumAlleles;
        if (retainsRareAlleles()) {
            myNumDataRows++;
        }
        myData0 = new OpenBitSet[myNumDataRows][myNumSites];
        myData1 = new OpenBitSet[myNumDataRows][myNumSites];
        for (int al = 0; al < myNumDataRows; al++) {
            for (int s = 0; s < myNumSites; s++) {
                myData0[al][s] = new OpenBitSet(getSequenceCount());
                myData1[al][s] = new OpenBitSet(getSequenceCount());
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
                            myData0[j][s].fastSet(t);
                            isRare = false;
                            break;
                        }
                    }
                    if (isRare && retainsRareAlleles()) {
                        myData0[myMaxNumAlleles][s].fastSet(t);
                    }
                }

                if (cb[1] != Alignment.UNKNOWN_ALLELE) {
                    boolean isRare = true;
                    for (int j = 0; j < myMaxNumAlleles; j++) {
                        if (cb[1] == myAlleles[s][j]) {
                            myData1[j][s].fastSet(t);
                            isRare = false;
                            break;
                        }
                    }
                    if (isRare && retainsRareAlleles()) {
                        myData1[myMaxNumAlleles][s].fastSet(t);
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
        myData0 = new OpenBitSet[myNumDataRows][myNumSites];
        myData1 = new OpenBitSet[myNumDataRows][myNumSites];
        for (int al = 0; al < myNumDataRows; al++) {
            for (int s = 0; s < myNumSites; s++) {
                myData0[al][s] = new OpenBitSet(getSequenceCount());
                myData1[al][s] = new OpenBitSet(getSequenceCount());
            }
        }
        for (int s = 0; s < myNumSites; s++) {
            for (int t = 0, n = getSequenceCount(); t < n; t++) {
                byte[] cb = a.getBaseArray(t, s);

                if (cb[0] != Alignment.UNKNOWN_ALLELE) {
                    boolean isRare = true;
                    for (int j = 0; j < myMaxNumAlleles; j++) {
                        if (cb[0] == myAlleles[s][j]) {
                            myData0[j][s].fastSet(t);
                            isRare = false;
                            break;
                        }
                    }
                    if (isRare && retainsRareAlleles()) {
                        myData0[myMaxNumAlleles][s].fastSet(t);
                    }
                }

                if (cb[1] != Alignment.UNKNOWN_ALLELE) {
                    boolean isRare = true;
                    for (int j = 0; j < myMaxNumAlleles; j++) {
                        if (cb[1] == myAlleles[s][j]) {
                            myData1[j][s].fastSet(t);
                            isRare = false;
                            break;
                        }
                    }
                    if (isRare && retainsRareAlleles()) {
                        myData1[myMaxNumAlleles][s].fastSet(t);
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
                if (myData0[i][site].fastGet(taxon)) {
                    result[0] = myAlleles[site][i];
                    break;
                }
            }

            for (int i = 0; i < myMaxNumAlleles; i++) {
                if (myData1[i][site].fastGet(taxon)) {
                    result[1] = myAlleles[site][i];
                    break;
                }
            }

            // Check For Rare Allele
            if (retainsRareAlleles() && myData0[myMaxNumAlleles][site].fastGet(taxon)) {
                result[0] = Alignment.RARE_ALLELE;
            }
            if (retainsRareAlleles() && myData1[myMaxNumAlleles][site].fastGet(taxon)) {
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
}
