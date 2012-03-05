/*
 * SBitAlignment
 */
package net.maizegenetics.pal.alignment;

import java.util.Random;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.SimpleIdGroup;
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

    protected SBitAlignment(Alignment a, int maxNumAlleles, boolean retainRareAlleles) {
        super(a, maxNumAlleles, retainRareAlleles);
        long currentTime = System.currentTimeMillis();
        loadAlleles(a);
        long prevTime = currentTime;
        currentTime = System.currentTimeMillis();
        System.out.println("Time to load alleles: " + ((currentTime - prevTime) / 1000));
    }

    protected SBitAlignment(IdGroup idGroup, byte[][] data, GeneticMap map, byte[] reference, String[][] alleleStates, int[] variableSites, int maxNumAlleles, Locus[] loci, int[] lociOffsets, String[] snpIDs, boolean retainRareAlleles) {
        super(idGroup, data, map, reference, alleleStates, variableSites, maxNumAlleles, loci, lociOffsets, snpIDs, retainRareAlleles);
        long currentTime = System.currentTimeMillis();
        loadAlleles(data);
        long prevTime = currentTime;
        currentTime = System.currentTimeMillis();
        System.out.println("Time to load alleles: " + ((currentTime - prevTime) / 1000));
    }

    // TESTING ONLY
    public static SBitAlignment getInstance() {

        int numTaxa = 11321;
        int numSites = 72425;

        String[] ids = new String[numTaxa];
        for (int i = 0; i < numTaxa; i++) {
            ids[i] = "Taxa" + i;
        }
        IdGroup idGroup = new SimpleIdGroup(ids);

        byte[][] data = new byte[numTaxa][numSites];
        Random random = new Random();
        for (int t = 0; t < numTaxa; t++) {
            for (int s = 0; s < numSites; s++) {
                byte temp = (byte) random.nextInt(6);
                byte value = (byte) ((temp << 4) | temp);
                data[t][s] = value;
            }
        }

        int[] variableSites = new int[numSites];
        String[] snpIDs = new String[numSites];
        for (int i = 0; i < numSites; i++) {
            variableSites[i] = i;
            snpIDs[i] = "SNPID_" + i;
        }

        return SBitAlignment.getNucleotideInstance(idGroup, data, null, null, variableSites, Alignment.DEFAULT_MAX_NUM_ALLELES, new Locus[]{new Locus("10", "10", 0, 0, null, null)}, new int[]{0}, snpIDs, true);

    }

    public static SBitAlignment getInstance(Alignment a) {
        return SBitAlignment.getInstance(a, a.getMaxNumAlleles(), a.retainsRareAlleles());
    }

    public static SBitAlignment getInstance(Alignment a, int maxNumAlleles, boolean retainRareAlleles) {

        if ((a instanceof SBitAlignment) && (a.getMaxNumAlleles() == maxNumAlleles) && (a.retainsRareAlleles() == retainRareAlleles)) {
            return (SBitAlignment) a;
        }

        String[][] alleleStates = a.getAlleleEncodings();
        if ((alleleStates == null) || (alleleStates.length == 0)) {
            throw new IllegalStateException("SBitAlignment: init: allele states should not be empty.");
        }

        if ((a instanceof SBitNucleotideAlignment) || (a instanceof TBitNucleotideAlignment)) {
            return new SBitNucleotideAlignment(a, maxNumAlleles, retainRareAlleles);
        } else if (alleleStates.length == 1) {
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
        int numSeqs = getSequenceCount();
        myData = new OpenBitSet[myNumDataRows][myNumSites];
        for (int al = 0; al < myNumDataRows; al++) {
            for (int s = 0; s < myNumSites; s++) {
                myData[al][s] = new OpenBitSet(numSeqs);
            }
        }
        ExecutorService pool = Executors.newFixedThreadPool(10);
        for (int s = 0; s < myNumSites; s++) {
            pool.execute(new ProcessSite(data, myData, s));
        }

        try {
            pool.shutdown();
            if (!pool.awaitTermination(120, TimeUnit.SECONDS)) {
                throw new IllegalStateException("ImportUtils: readFromHapmap: processing threads timed out.");
            }
        } catch (Exception e) {
            throw new IllegalStateException("ImportUtils: readFromHapmap: processing threads problem.");
        }

    }

    private class ProcessSite implements Runnable {

        private OpenBitSet[][] myData;
        private byte[][] myOrigData;
        private int mySite;

        public ProcessSite(byte[][] origData, OpenBitSet[][] data, int site) {
            myData = data;
            myOrigData = origData;
            mySite = site;
        }

        public void run() {
            int numSeqs = getSequenceCount();
            byte[] cb = new byte[2];
            for (int t = 0; t < numSeqs; t++) {
                cb[0] = (byte) ((myOrigData[t][mySite] >>> 4) & 0xf);
                cb[1] = (byte) (myOrigData[t][mySite] & 0xf);
                for (int i = 0; i < 2; i++) {
                    if (cb[i] != Alignment.UNKNOWN_ALLELE) {
                        boolean isRare = true;
                        for (int j = 0; j < myMaxNumAlleles; j++) {
                            if (cb[i] == myAlleles[mySite][j]) {
                                myData[j][mySite].fastSet(t);
                                isRare = false;
                                break;
                            }
                        }
                        if (isRare && retainsRareAlleles()) {
                            myData[myMaxNumAlleles][mySite].fastSet(t);
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
        int numSeqs = getSequenceCount();
        myData = new OpenBitSet[myNumDataRows][myNumSites];
        for (int al = 0; al < myNumDataRows; al++) {
            for (int s = 0; s < myNumSites; s++) {
                myData[al][s] = new OpenBitSet(numSeqs);
            }
        }
        for (int s = 0; s < myNumSites; s++) {
            for (int t = 0; t < numSeqs; t++) {
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

    @Override
    public int getTotalGametesNotMissing(int site) {

        OpenBitSet temp = new OpenBitSet(getSequenceCount());
        for (int i = 0; i < myNumDataRows; i++) {
            temp.or(myData[i][site]);
        }
        return ((int) temp.cardinality()) * 2;

    }

    @Override
    public int getMinorAlleleCount(int site) {

        if ((myMaxNumAlleles < 2) || (myAlleles[site][1] == Alignment.UNKNOWN_ALLELE)) {
            return 0;
        }

        OpenBitSet temp = new OpenBitSet(getSequenceCount());
        for (int i = 0; i < myNumDataRows; i++) {
            if (i != 1) {
                temp.or(myData[i][site]);
            }
        }
        temp.flip(0, temp.size());
        temp.and(myData[1][site]);

        return (int) temp.cardinality() + (int) myData[1][site].cardinality();

    }

    @Override
    public double getMinorAlleleFrequency(int site) {
        int minorAlleleCount = getMinorAlleleCount(site);
        if (minorAlleleCount == 0) {
            return 0.0;
        }
        return (double) minorAlleleCount / (double) getTotalGametesNotMissing(site);
    }

    @Override
    public double getMajorAlleleFrequency(int site) {
        int majorAlleleCount = getMajorAlleleCount(site);
        if (majorAlleleCount == 0) {
            return 0.0;
        }
        return (double) majorAlleleCount / (double) getTotalGametesNotMissing(site);
    }

    @Override
    public int getMajorAlleleCount(int site) {

        if (myAlleles[site][0] == Alignment.UNKNOWN_ALLELE) {
            return 0;
        }

        OpenBitSet temp = new OpenBitSet(getSequenceCount());
        for (int i = 1; i < myNumDataRows; i++) {
            temp.or(myData[i][site]);
        }
        temp.flip(0, temp.size());
        temp.and(myData[0][site]);

        return (int) temp.cardinality() + (int) myData[0][site].cardinality();

    }

    @Override
    public boolean isHeterozygous(int taxon, int site) {
        int count = 0;
        for (int i = 0; i < myNumDataRows; i++) {
            if (myData[i][site].fastGet(taxon)) {
                count++;
                if (count == 2) {
                    return true;
                }
            }
        }
        return false;
    }

    @Override
    public int getHeterozygousCount(int site) {

        int result = 0;
        for (int i = 0; i < myNumDataRows; i++) {
            for (int j = i + 1; j < myNumDataRows; j++) {
                result += (int) OpenBitSet.intersectionCount(myData[i][site], myData[j][site]);
            }
        }
        return result;

    }

    @Override
    public boolean isPolymorphic(int site) {
        boolean nonZero = false;
        for (int i = 0; i < myNumDataRows; i++) {
            int numTaxa = (int) myData[i][site].cardinality();
            if (numTaxa != 0) {
                if (nonZero) {
                    return true;
                }
                nonZero = true;
            }
        }
        return false;
    }

    @Override
    public Object[][] getDiploidCounts() {

        if (myAlleleStates.length != 1) {
            return super.getDiploidCounts();
        }

        long[][] counts = new long[16][16];
        for (int site = 0; site < myNumSites; site++) {
            for (int i = 0; i < myMaxNumAlleles; i++) {
                byte indexI = myAlleles[site][i];
                counts[indexI][indexI] += myData[i][site].cardinality();
                for (int j = i + 1; j < myMaxNumAlleles; j++) {
                    byte indexJ = myAlleles[site][j];
                    long ijHet = OpenBitSet.intersectionCount(myData[i][site], myData[j][site]);
                    if (indexI < indexJ) {
                        counts[indexI][indexJ] += ijHet;
                    } else {
                        counts[indexJ][indexI] += ijHet;
                    }
                    counts[indexI][indexI] -= ijHet;
                    counts[indexJ][indexJ] -= ijHet;
                }
            }
        }

        int numAlleles = 0;
        long unknownCount = (long) getSequenceCount() * (long) myNumSites;
        for (byte x = 0; x < 16; x++) {
            for (byte y = x; y < 16; y++) {
                if (counts[x][y] != 0) {
                    numAlleles++;
                    unknownCount -= counts[x][y];
                }
            }
        }

        if (unknownCount > 0) {
            numAlleles++;
        }

        Object[][] result = new Object[2][numAlleles];
        int nextResult = 0;
        for (byte x = 0; x < 16; x++) {
            for (byte y = x; y < 16; y++) {
                if (counts[x][y] != 0) {
                    byte value = (byte) ((x << 4) | y);
                    result[0][nextResult] = getDiploidAsString(0, value);
                    result[1][nextResult++] = counts[x][y];
                }
            }
        }

        if (unknownCount > 0) {
            result[0][nextResult] = getDiploidAsString(0, UNKNOWN_DIPLOID_ALLELE);
            result[1][nextResult] = unknownCount;
        }

        boolean change = true;
        while (change) {

            change = false;

            for (int k = 0; k < numAlleles - 1; k++) {

                if ((Long) result[1][k] < (Long) result[1][k + 1]) {

                    Object temp = result[0][k];
                    result[0][k] = result[0][k + 1];
                    result[0][k + 1] = temp;

                    Object tempCount = result[1][k];
                    result[1][k] = result[1][k + 1];
                    result[1][k + 1] = tempCount;

                    change = true;
                }
            }

        }

        return result;
    }

    @Override
    public Object[][] getDiploidssSortedByFrequency(int site) {

        if (myAlleleStates.length != 1) {
            return super.getDiploidssSortedByFrequency(site);
        }

        int[][] counts = new int[16][16];
        for (int i = 0; i < myMaxNumAlleles; i++) {
            byte indexI = myAlleles[site][i];
            counts[indexI][indexI] += (int) myData[i][site].cardinality();
            for (int j = i + 1; j < myMaxNumAlleles; j++) {
                byte indexJ = myAlleles[site][j];
                int ijHet = (int) OpenBitSet.intersectionCount(myData[i][site], myData[j][site]);
                if (indexI < indexJ) {
                    counts[indexI][indexJ] += ijHet;
                } else {
                    counts[indexJ][indexI] += ijHet;
                }
                counts[indexI][indexI] -= ijHet;
                counts[indexJ][indexJ] -= ijHet;
            }
        }

        int numAlleles = 0;
        int unknownCount = getSequenceCount();
        for (byte x = 0; x < 16; x++) {
            for (byte y = x; y < 16; y++) {
                if (counts[x][y] != 0) {
                    numAlleles++;
                    unknownCount -= counts[x][y];
                }
            }
        }

        if (unknownCount > 0) {
            numAlleles++;
        }

        Object[][] result = new Object[2][numAlleles];
        int nextResult = 0;
        for (byte x = 0; x < 16; x++) {
            for (byte y = x; y < 16; y++) {
                if (counts[x][y] != 0) {
                    byte value = (byte) ((x << 4) | y);
                    result[0][nextResult] = getDiploidAsString(0, value);
                    result[1][nextResult++] = counts[x][y];
                }
            }
        }

        if (unknownCount > 0) {
            result[0][nextResult] = getDiploidAsString(0, UNKNOWN_DIPLOID_ALLELE);
            result[1][nextResult] = unknownCount;
        }

        boolean change = true;
        while (change) {

            change = false;

            for (int k = 0; k < numAlleles - 1; k++) {

                if ((Integer) result[1][k] < (Integer) result[1][k + 1]) {

                    Object temp = result[0][k];
                    result[0][k] = result[0][k + 1];
                    result[0][k + 1] = temp;

                    Object tempCount = result[1][k];
                    result[1][k] = result[1][k + 1];
                    result[1][k + 1] = tempCount;

                    change = true;
                }
            }

        }

        return result;

    }

    @Override
    public int[][] getAllelesSortedByFrequency(int site) {

        int[] counts = new int[16];
        for (int i = 0; i < myNumDataRows; i++) {
            byte indexI;
            if ((retainsRareAlleles()) && (i == myMaxNumAlleles)) {
                indexI = Alignment.RARE_ALLELE;
            } else {
                indexI = myAlleles[site][i];
            }
            counts[indexI] += (int) myData[i][site].cardinality() * 2;
            for (int j = i + 1; j < myNumDataRows; j++) {
                byte indexJ;
                if ((retainsRareAlleles()) && (j == myMaxNumAlleles)) {
                    indexJ = Alignment.RARE_ALLELE;
                } else {
                    indexJ = myAlleles[site][j];
                }
                int ijHet = (int) OpenBitSet.intersectionCount(myData[i][site], myData[j][site]);
                counts[indexI] -= ijHet;
                counts[indexJ] -= ijHet;
            }
        }

        int numAlleles = 0;
        for (byte x = 0; x < Alignment.UNKNOWN_ALLELE; x++) {
            if (counts[x] != 0) {
                numAlleles++;
            }
        }

        int current = 0;
        int[][] result = new int[2][numAlleles];
        for (byte x = 0; x < Alignment.UNKNOWN_ALLELE; x++) {
            if (counts[x] != 0) {
                result[0][current] = x;
                result[1][current++] = counts[x];
            }
        }

        boolean change = true;
        while (change) {

            change = false;

            for (int k = 0; k < numAlleles - 1; k++) {

                if (result[1][k] < result[1][k + 1]) {

                    int temp = result[0][k];
                    result[0][k] = result[0][k + 1];
                    result[0][k + 1] = temp;

                    int tempCount = result[1][k];
                    result[1][k] = result[1][k + 1];
                    result[1][k + 1] = tempCount;

                    change = true;
                }
            }

        }

        return result;

    }
}
