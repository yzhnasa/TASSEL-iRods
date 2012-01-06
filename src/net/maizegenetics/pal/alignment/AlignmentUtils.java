/*
 * AlignmentUtils
 */
package net.maizegenetics.pal.alignment;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

/**
 *
 * @author terry
 */
public class AlignmentUtils {

    private static final Integer ONE = new Integer(1);

    private AlignmentUtils() {
        // utility class
    }

    public static int[][] getAllelesSortedByFrequency(byte[][] data, int site) {

        int[] stateCnt = new int[16];
        for (int i = 0; i < data.length; i++) {
            byte first = (byte) ((data[i][site] >>> 4) & 0xf);
            byte second = (byte) (data[i][site] & 0xf);
            if (first != Alignment.UNKNOWN_ALLELE) {
                stateCnt[first]++;
            }
            if (second != Alignment.UNKNOWN_ALLELE) {
                stateCnt[second]++;
            }
        }

        int count = 0;
        for (int j = 0; j < 16; j++) {
            if (stateCnt[j] != 0) {
                count++;
            }
        }

        int result[][] = new int[2][count];
        int index = 0;
        for (int k = 0; k < 16; k++) {
            if (stateCnt[k] != 0) {
                result[0][index] = k;
                result[1][index] = stateCnt[k];
                index++;
            }
        }

        boolean change = true;
        while (change) {

            change = false;

            for (int k = 0; k < count - 1; k++) {

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

    public static byte[] getAlleles(byte[][] data, int site) {
        int[][] alleles = getAllelesSortedByFrequency(data, site);
        int resultSize = alleles[0].length;
        byte[] result = new byte[resultSize];
        for (int i = 0; i < resultSize; i++) {
            result[i] = (byte) alleles[0][i];
        }
        return result;
    }

    public static Object[][] getAllelesSortedByFrequency(String[][] data, int site) {

        Map<String, Integer> stateCnt = new HashMap();
        for (int i = 0; i < data.length; i++) {
            String[] temp = data[i][site].split(":");
            String first;
            String second;
            if ((temp == null) || (temp.length == 0)) {
                first = second = Alignment.UNKNOWN_ALLELE_STR;
            } else if (temp.length == 1) {
                first = second = temp[0].trim();
            } else {
                first = temp[0].trim();
                second = temp[1].trim();
            }
            if (!first.equalsIgnoreCase(Alignment.UNKNOWN_ALLELE_STR)) {
                Integer count = (Integer) stateCnt.get(first);
                if (count == null) {
                    stateCnt.put(first, ONE);
                } else {
                    stateCnt.put(first, count + 1);
                }
            }
            if (!second.equalsIgnoreCase(Alignment.UNKNOWN_ALLELE_STR)) {
                Integer count = (Integer) stateCnt.get(second);
                if (count == null) {
                    stateCnt.put(second, ONE);
                } else {
                    stateCnt.put(second, count + 1);
                }
            }
        }

        int count = stateCnt.size();

        Object[][] result = new Object[2][count];
        Iterator itr = stateCnt.keySet().iterator();
        int index = 0;
        while (itr.hasNext()) {
            String key = (String) itr.next();
            result[0][index] = key;
            result[1][index] = (Integer) stateCnt.get(key);
            index++;
        }

        boolean change = true;
        while (change) {

            change = false;

            for (int k = 0; k < count - 1; k++) {

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

    public static List<String> getAlleles(String[][] data, int site) {
        Object[][] alleles = getAllelesSortedByFrequency(data, site);
        if ((alleles == null) || (alleles.length == 0)) {
            return null;
        }
        int resultSize = alleles[0].length;
        List result = new ArrayList();
        for (int i = 0; i < resultSize; i++) {
            result.add((String) alleles[0][i]);
        }
        return result;
    }

    public static String[][] getAlleleStates(String[][] data, int maxNumAlleles) {

        int numSites = data[0].length;

        String[][] alleleStates = new String[numSites][16];
        for (int i = 0; i < numSites; i++) {
            for (int j = 0; j < 16; j++) {
                if (j == Alignment.RARE_ALLELE) {
                    alleleStates[i][j] = Alignment.RARE_ALLELE_STR;
                } else {
                    alleleStates[i][j] = Alignment.UNKNOWN_ALLELE_STR;
                }
            }
        }

        for (int site = 0; site < numSites; site++) {
            List alleles = AlignmentUtils.getAlleles(data, site);
            if (alleles != null) {
                int numAlleles = Math.min(alleles.size(), maxNumAlleles);
                for (int k = 0; k < numAlleles; k++) {
                    alleleStates[site][k] = (String) alleles.get(k);
                }
            }
        }

        return alleleStates;

    }

    public static byte[][] getDataBytes(String[] data) {

        int numTaxa = data.length;

        int numSites = data[0].length();

        byte[][] dataBytes = new byte[numTaxa][numSites];

        for (int site = 0; site < numSites; site++) {
            for (int taxon = 0; taxon < numTaxa; taxon++) {
                dataBytes[taxon][site] = NucleotideAlignmentConstants.getNucleotideDiploidByte(data[taxon].charAt(site));
            }
        }

        return dataBytes;

    }

    public static byte[][] getDataBytes(String[][] data, String[][] alleleStates, int maxNumAlleles) {

        int numTaxa = data.length;

        int numSites = data[0].length;

        byte[][] dataBytes = new byte[numTaxa][numSites];

        if (alleleStates.length == 1) {
            for (int site = 0; site < numSites; site++) {
                setDataBytes(data, alleleStates[0], maxNumAlleles, numTaxa, site, dataBytes);
            }
        } else {
            for (int site = 0; site < numSites; site++) {
                setDataBytes(data, alleleStates[site], maxNumAlleles, numTaxa, site, dataBytes);
            }
        }

        return dataBytes;

    }

    private static void setDataBytes(String[][] data, String[] alleleStates, int maxNumAlleles, int numTaxa, int site, byte[][] dataBytes) {

        for (int taxon = 0; taxon < numTaxa; taxon++) {
            if (data[taxon][site].equalsIgnoreCase(Alignment.UNKNOWN_ALLELE_STR)) {
                dataBytes[taxon][site] = Alignment.UNKNOWN_ALLELE;
            } else {
                dataBytes[taxon][site] = Alignment.RARE_ALLELE;
                for (int k = 0; k < maxNumAlleles; k++) {
                    if (alleleStates[k].equals(data[taxon][site])) {
                        dataBytes[taxon][site] = (byte) k;
                        break;
                    }
                }
            }
        }

    }

    /**
     * Counts the number of pairwise allele matches. For identical
     * homozygotes, this equals n * n, where n is the ploidy level.
     * The probability that alleles drawn from the two loci are
     * identical in state, P(IBS), equals this count divided by n*n.
     * If either genotype is unknown returns -1.
     *
     * @param b1        the genotype for the first site
     * @param b2        the genotype for the second site
     *
     * @return  a count of the number of times pairs of alleles match
     * examining all possible pairwise combinations
     */
    public static int getDiploidIdentity(byte b1, byte b2) {
        throw new UnsupportedOperationException();
    }

    /**
     * remove sites based on minimum frequency (the count of good bases, INCLUDING GAPS)
     * and based on the proportion of good alleles (including gaps) different from consensus
     *
     * @param aa the AnnotatedAlignment to filter
     * @param minimumProportion minimum proportion of sites different from the consensus
     * @param minimumCount      minimum number of sequences with a good bases (not N or ?), where GAP IS CONSIDERED A GOOD BASE
     */
    public static Alignment removeSitesBasedOnFreqIgnoreMissing(Alignment aa, double minimumProportion, int minimumCount) {
        int[] includeSites = getIncludedSitesBasedOnFreqIgnoreMissing(aa, minimumProportion, minimumCount);
        Alignment mlaa = FilterAlignment.getInstance(aa, includeSites);
        return mlaa;
    }

    /**
     * get sites to be included based on minimum frequency (the count of good
     * bases, INCLUDING GAPS) and based on the proportion of good sites (INCLUDING
     * GAPS) different from consensus
     *
     * @param aa the AnnotatedAlignment to filter
     * @param minimumProportion minimum proportion of sites different from the consensus
     * @param minimumCount      minimum number of sequences with a good base or a gap (but not N or ?)
     */
    public static int[] getIncludedSitesBasedOnFreqIgnoreMissing(Alignment aa, double minimumProportion, int minimumCount) {
        ArrayList<Integer> includeAL = new ArrayList<Integer>();
        for (int i = 0, n = aa.getSiteCount(); i < n; i++) {

            int totalNonMissing = aa.getTotalCountNotMissing(i);

            double obsMinProp = aa.getMinorAlleleFrequency(i);

            if ((totalNonMissing > 0) && (totalNonMissing >= minimumCount) && (obsMinProp >= minimumProportion)) {
                includeAL.add(i);
            }
        }
        int[] includeSites = new int[includeAL.size()];
        for (int i = 0; i < includeAL.size(); i++) {
            includeSites[i] = includeAL.get(i);
        }
        return includeSites;
    }

    /**
     * Remove sites based on site position (excluded sites are <firstSite and >lastSite)
     * This not effect any prior exclusions.
     * 
     * @param aa the AnnotatedAlignment to filter
     * @param firstSite first site to keep in the range
     * @param lastSite  last site to keep in the range
     */
    public static Alignment removeSitesOutsideRange(Alignment aa, int firstSite, int lastSite) {
        if ((firstSite < 0) || (firstSite > lastSite)) {
            return null;
        }
        if (lastSite > aa.getSiteCount() - 1) {
            return null;
        }
        return FilterAlignment.getInstance(aa, firstSite, lastSite);
    }
}
