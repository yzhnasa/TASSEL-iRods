/*
 * SBitAlignment
 */
package net.maizegenetics.pal.alignment;

import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Reader;
import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.SimpleIdGroup;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;
import net.maizegenetics.util.UnmodifiableBitSet;

/**
 * This data alignment reads from HDF5 files a native SBit alignment
 *
 * @author Ed & Terry
 */
public class SBitAlignmentNucleotideHDF5 extends AbstractAlignment {

    private OpenBitSet[] myData;
    private int myNumDataRows;
    private IHDF5Reader h5 = null;
    private int cachedRow = 0;
    private int[] myVariableSites;
    private int myNumWords = 0;
    private String myLocusPath;
    private Locus[] myLoci;

    public SBitAlignmentNucleotideHDF5(String theHDF5file, Locus[] loci, IdGroup idGroup,
            String[][] alleleStates, int[] positions, byte[][] alleles) {
        super(idGroup, alleleStates);
        h5 = HDF5Factory.openForReading(theHDF5file);
        myLoci = loci;
        myNumWords = h5.getIntAttribute("/", "numWords");
        myNumSites = positions.length;
        myVariableSites = positions;
        myAlleles = alleles;
        myMaxNumAlleles = myAlleles[0].length;
        myLocusPath = "L:" + loci[0].getName();
        cachSiteRow(0);
    }

    public static SBitAlignmentNucleotideHDF5 getInstance(String theHDF5file, String locusName) {
        IHDF5Reader h5 = HDF5Factory.openForReading(theHDF5file);
        IdGroup idg = new SimpleIdGroup(h5.readStringArray("taxaNames"));
        String[][] alleleStates = {h5.readStringArray("alleleStates"),};
        String locusPath = "L:" + locusName;
        int[] positions = h5.readIntArray(locusPath + "/positions");
        byte[][] alleles = h5.readByteMatrix(locusPath + "/alleles");
        Locus[] theLocus = new Locus[]{new Locus(locusName, locusName, 0, positions[positions.length - 1], null, null)};
        h5.close();
        return new SBitAlignmentNucleotideHDF5(theHDF5file, theLocus, idg,
                alleleStates, positions, alleles);
    }

    private void cachSiteRow(int site) {
        myData = new OpenBitSet[myMaxNumAlleles];
        for (int aNum = 0; aNum < myMaxNumAlleles; aNum++) {
            myData[aNum] = new OpenBitSet(h5.readLongMatrixBlockWithOffset(myLocusPath + "/" + aNum, 1, myNumWords, site, 0)[0], myNumWords);
        }
        cachedRow = site;
    }

    @Override
    public byte getBase(int taxon, int site) {
        byte[] temp = getBaseArray(taxon, site);
        return (byte) ((temp[0] << 4) | temp[1]);
    }

    @Override
    public byte[] getBaseArray(int taxon, int site) {
        if (site != cachedRow) {
            cachSiteRow(site);
        }
        byte[] result = new byte[2];
        result[0] = Alignment.UNKNOWN_ALLELE;
        result[1] = Alignment.UNKNOWN_ALLELE;
        try {
            int count = 0;
            for (int i = 0; i < myMaxNumAlleles; i++) {
                if (myData[i].fastGet(taxon)) {
                    if (count == 0) {
                        result[1] = myAlleles[site][i];
                    }
                    result[count++] = myAlleles[site][i];
                }
            }

            // Check For Rare Allele
            if (retainsRareAlleles() && myData[myMaxNumAlleles].fastGet(taxon)) {
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
        if (site != cachedRow) {
            cachSiteRow(site);
        }
        return UnmodifiableBitSet.getInstance(myData[alleleNumber]);
    }

    @Override
    public int getTotalGametesNotMissing(int site) {
        if (site != cachedRow) {
            cachSiteRow(site);
        }
        OpenBitSet temp = new OpenBitSet(getSequenceCount());
        for (int i = 0; i < myNumDataRows; i++) {
            temp.or(myData[i]);
        }
        return ((int) temp.cardinality()) * 2;

    }

    @Override
    public int getMinorAlleleCount(int site) {
        if (site != cachedRow) {
            cachSiteRow(site);
        }
        if ((myMaxNumAlleles < 2) || (myAlleles[site][1] == Alignment.UNKNOWN_ALLELE)) {
            return 0;
        }

        OpenBitSet temp = new OpenBitSet(getSequenceCount());
        for (int i = 0; i < myNumDataRows; i++) {
            if (i != 1) {
                temp.or(myData[i]);
            }
        }
        temp.flip(0, temp.size());
        temp.and(myData[1]);

        return (int) temp.cardinality() + (int) myData[1].cardinality();

    }

    @Override
    public double getMinorAlleleFrequency(int site) {
        if (site != cachedRow) {
            cachSiteRow(site);
        }
        int minorAlleleCount = getMinorAlleleCount(site);
        if (minorAlleleCount == 0) {
            return 0.0;
        }
        return (double) minorAlleleCount / (double) getTotalGametesNotMissing(site);
    }

    @Override
    public double getMajorAlleleFrequency(int site) {
        if (site != cachedRow) {
            cachSiteRow(site);
        }
        int majorAlleleCount = getMajorAlleleCount(site);
        if (majorAlleleCount == 0) {
            return 0.0;
        }
        return (double) majorAlleleCount / (double) getTotalGametesNotMissing(site);
    }

    @Override
    public int getMajorAlleleCount(int site) {
        if (site != cachedRow) {
            cachSiteRow(site);
        }
        if (myAlleles[site][0] == Alignment.UNKNOWN_ALLELE) {
            return 0;
        }

        OpenBitSet temp = new OpenBitSet(getSequenceCount());
        for (int i = 1; i < myNumDataRows; i++) {
            temp.or(myData[i]);
        }
        temp.flip(0, temp.size());
        temp.and(myData[0]);

        return (int) temp.cardinality() + (int) myData[0].cardinality();

    }

    @Override
    public boolean isHeterozygous(int taxon, int site) {
        if (site != cachedRow) {
            cachSiteRow(site);
        }
        int count = 0;
        for (int i = 0; i < myNumDataRows; i++) {
            if (myData[i].fastGet(taxon)) {
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
        if (site != cachedRow) {
            cachSiteRow(site);
        }
        int result = 0;
        for (int i = 0; i < myNumDataRows; i++) {
            for (int j = i + 1; j < myNumDataRows; j++) {
                result += (int) OpenBitSet.intersectionCount(myData[i], myData[j]);
            }
        }
        return result;

    }

    @Override
    public boolean isPolymorphic(int site) {
        if (site != cachedRow) {
            cachSiteRow(site);
        }
        boolean nonZero = false;
        for (int i = 0; i < myNumDataRows; i++) {
            int numTaxa = (int) myData[i].cardinality();
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
            if (site != cachedRow) {
                cachSiteRow(site);
            }
            for (int i = 0; i < myMaxNumAlleles; i++) {
                byte indexI = myAlleles[site][i];
                counts[indexI][indexI] += myData[i].cardinality();
                for (int j = i + 1; j < myMaxNumAlleles; j++) {
                    byte indexJ = myAlleles[site][j];
                    long ijHet = OpenBitSet.intersectionCount(myData[i], myData[j]);
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
        if (site != cachedRow) {
            cachSiteRow(site);
        }
        if (myAlleleStates.length != 1) {
            return super.getDiploidssSortedByFrequency(site);
        }

        int[][] counts = new int[16][16];
        for (int i = 0; i < myMaxNumAlleles; i++) {
            byte indexI = myAlleles[site][i];
            counts[indexI][indexI] += (int) myData[i].cardinality();
            for (int j = i + 1; j < myMaxNumAlleles; j++) {
                byte indexJ = myAlleles[site][j];
                int ijHet = (int) OpenBitSet.intersectionCount(myData[i], myData[j]);
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
        if (site != cachedRow) {
            cachSiteRow(site);
        }
        int[] counts = new int[16];
        for (int i = 0; i < myNumDataRows; i++) {
            byte indexI;
            if ((retainsRareAlleles()) && (i == myMaxNumAlleles)) {
                indexI = Alignment.RARE_ALLELE;
            } else {
                indexI = myAlleles[site][i];
            }
            counts[indexI] += (int) myData[i].cardinality() * 2;
            for (int j = i + 1; j < myNumDataRows; j++) {
                byte indexJ;
                if ((retainsRareAlleles()) && (j == myMaxNumAlleles)) {
                    indexJ = Alignment.RARE_ALLELE;
                } else {
                    indexJ = myAlleles[site][j];
                }
                int ijHet = (int) OpenBitSet.intersectionCount(myData[i], myData[j]);
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
    public String getBaseAsString(int taxon, int site) {
        return NucleotideAlignmentConstants.getNucleotideIUPAC(getBase(taxon, site));
    }

    @Override
    public int getPositionInLocus(int site) {
        try {
            return myVariableSites[site];
        } catch (Exception e) {
            return site;
        }
    }

    @Override
    public Locus getLocus(int site) {
        return myLoci[0];
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
            if ((alleles[0][i] == 4) || (alleles[0][i] == 5)) {
                return true;
            }
        }
        return false;
    }

    public static void main(String[] args) {
        String infile = "/Users/terry/TASSELTutorialData/data/mdp_genotype.hmp.txt";
        String outfile = "/Users/terry/TASSELTutorialData/data/mdp_genotype.hmp.h5";
        BitAlignment a = (BitAlignment) ImportUtils.readFromHapmap(infile, null);
        ExportUtils.writeToHDF5(a, outfile);
        SBitAlignmentNucleotideHDF5 sbah2 = SBitAlignmentNucleotideHDF5.getInstance(outfile, "10");
        String outHap = "/Users/terry/TASSELTutorialData/data/testing_h5.hmp.txt";
        ExportUtils.writeToHapmap(sbah2, false, outHap, '\t', null);
        long time = System.currentTimeMillis();
        for (int i = 0; i < sbah2.getSiteCount(); i++) {
            if (a.getBase(i % 100, i) != sbah2.getBase(i % 100, i)) {
                System.out.println("Error:" + i);
            }
        }
        System.out.println("Accessing time:" + (System.currentTimeMillis() - time));

    }
}
