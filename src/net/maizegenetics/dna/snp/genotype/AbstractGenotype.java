/*
 *  AbstractGenotype
 */
package net.maizegenetics.dna.snp.genotype;

import net.maizegenetics.dna.snp.Alignment;
import net.maizegenetics.dna.snp.AlignmentUtils;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;

import org.apache.log4j.Logger;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

/**
 *
 * @author Terry Casstevens
 */
abstract class AbstractGenotype implements Genotype {

    private static final Logger myLogger = Logger.getLogger(AbstractGenotype.class);
    private static final int DEFAULT_MAX_NUM_ALLELES = NucleotideAlignmentConstants.NUMBER_NUCLEOTIDE_ALLELES;
    protected final int myTaxaCount;
    protected final int mySiteCount;
    private final String[][] myAlleleEncodings;
    private final boolean myIsPhased;
    private final AlleleFreqCache myAlleleFreqCache;

    AbstractGenotype(int numTaxa, int numSites, boolean phased, String[][] alleleEncodings, int maxNumAlleles) {
        myTaxaCount = numTaxa;
        mySiteCount = numSites;
        myIsPhased = phased;
        myAlleleEncodings = alleleEncodings;
        myAlleleFreqCache = new AlleleFreqCache(this, maxNumAlleles);
    }

    AbstractGenotype(int numTaxa, int numSites, boolean phased, String[][] alleleEncodings) {
        this(numTaxa, numSites, phased, alleleEncodings, DEFAULT_MAX_NUM_ALLELES);
    }

    @Override
    public byte[] genotypeArray(int taxon, int site) {
        return AlignmentUtils.getDiploidValues(genotype(taxon, site));
    }

    @Override
    public byte[] genotypeRange(int taxon, int startSite, int endSite) {
        byte[] result = new byte[endSite - startSite];
        for (int i = startSite; i < endSite; i++) {
            result[i - startSite] = genotype(taxon, i);
        }
        return result;
    }

    @Override
    public byte[] genotypeAllSites(int taxon) {
        byte[] result = new byte[mySiteCount];
        for (int i = 0; i < mySiteCount; i++) {
            result[i] = genotype(taxon, i);
        }
        return result;
    }

    @Override
    public String genotypeAsString(int taxon, int site) {
        String[][] alleleStates = alleleDefinitions();
        byte[] temp = genotypeArray(taxon, site);
        return alleleStates[0][temp[0]] + ":" + alleleStates[0][temp[1]];
    }

    @Override
    public String genotypeAsStringRange(int taxon, int startSite, int endSite) {
        StringBuilder builder = new StringBuilder();
        for (int i = startSite; i < endSite; i++) {
            if (i != startSite) {
                builder.append(";");
            }
            builder.append(genotypeAsString(taxon, i));
        }
        return builder.toString();
    }

    @Override
    public String genotypeAsStringRow(int taxon) {
        return genotypeAsStringRange(taxon, 0, mySiteCount);
    }

    @Override
    public String[] genotypeAsStringArray(int taxon, int site) {
        String[][] alleleStates = alleleDefinitions();
        byte[] temp = genotypeArray(taxon, site);
        return new String[]{alleleStates[0][temp[0]], alleleStates[0][temp[1]]};
    }

    @Override
    public int[][] getAllelesSortedByFrequency(int site) {
        return myAlleleFreqCache.getAllelesSortedByFrequency(site);
    }

    @Override
    public boolean isHeterozygous(int taxon, int site) {
        byte[] values = genotypeArray(taxon, site);
        if (values[0] == values[1]) {
            return false;
        } else {
            return true;
        }
    }

    @Override
    public int getHeterozygousCount(int site) {
        int result = 0;
        for (int i = 0, n = myTaxaCount; i < n; i++) {
            if (isHeterozygous(i, site)) {
                result++;
            }
        }
        return result;
    }

    @Override
    public boolean isPolymorphic(int site) {

        byte first = Alignment.UNKNOWN_ALLELE;
        for (int i = 0, n = myTaxaCount; i < n; i++) {
            byte[] current = genotypeArray(i, site);
            if (current[0] != Alignment.UNKNOWN_ALLELE) {
                if (first == Alignment.UNKNOWN_ALLELE) {
                    first = current[0];
                } else if (first != current[0]) {
                    return true;
                }
            }
            if (current[1] != Alignment.UNKNOWN_ALLELE) {
                if (first == Alignment.UNKNOWN_ALLELE) {
                    first = current[1];
                } else if (first != current[1]) {
                    return true;
                }
            }
        }

        return false;

    }

    @Override
    public boolean isAllPolymorphic() {

        for (int i = 0, n = mySiteCount; i < n; i++) {
            if (!isPolymorphic(i)) {
                return false;
            }
        }

        return true;

    }

    @Override
    public boolean isPhased() {
        return myIsPhased;
    }

    @Override
    public boolean retainsRareAlleles() {
        return true;
    }

    @Override
    public String[][] alleleDefinitions() {
        return myAlleleEncodings;
    }

    @Override
    public String[] alleleDefinitions(int site) {
        if (myAlleleEncodings.length == 1) {
            return myAlleleEncodings[0];
        } else {
            return myAlleleEncodings[site];
        }
    }

    @Override
    public String genotypeAsString(int site, byte value) {
        return alleleDefinitions(site)[value];
    }

    @Override
    public String getDiploidAsString(int site, byte value) {
        String[] alleleStates = alleleDefinitions(site);
        return alleleStates[(value >>> 4) & 0xf] + ":" + alleleStates[value & 0xf];
    }

    @Override
    public int getMaxNumAlleles() {
        return DEFAULT_MAX_NUM_ALLELES;
    }

    @Override
    public int getTotalGametesNotMissing(int site) {

        int result = 0;
        for (int i = 0, n = myTaxaCount; i < n; i++) {
            byte[] current = genotypeArray(i, site);
            if (current[0] != Alignment.UNKNOWN_ALLELE) {
                result++;
            }
            if (current[1] != Alignment.UNKNOWN_ALLELE) {
                result++;
            }
        }
        return result;

    }

    @Override
    public int getTotalNotMissing(int site) {

        int result = 0;
        for (int i = 0, n = myTaxaCount; i < n; i++) {
            byte current = genotype(i, site);
            if (current != Alignment.UNKNOWN_DIPLOID_ALLELE) {
                result++;
            }
        }
        return result;

    }

    @Override
    public byte[] getMajorAlleleForAllSites() {
        byte[] result = new byte[mySiteCount];
        for (int i = 0; i < mySiteCount; i++) {
            result[i] = getMajorAllele(i);
        }
        return result;
    }

    @Override
    public byte[] getMinorAlleleForAllSites() {
        byte[] result = new byte[mySiteCount];
        for (int i = 0; i < mySiteCount; i++) {
            result[i] = getMinorAllele(i);
        }
        return result;
    }

    @Override
    public int getMinorAlleleCount(int site) {

        int[][] alleles = getAllelesSortedByFrequency(site);

        if (alleles[0].length >= 2) {
            return alleles[1][1];
        } else {
            return 0;
        }

    }

    @Override
    public byte getMinorAllele(int site) {
        int[][] alleles = getAllelesSortedByFrequency(site);

        if (alleles[0].length >= 2) {
            return (byte) alleles[0][1];
        } else {
            return Alignment.UNKNOWN_ALLELE;
        }
    }

    @Override
    public String getMinorAlleleAsString(int site) {
        return genotypeAsString(site, getMinorAllele(site));
    }

    @Override
    public byte[] getMinorAlleles(int site) {
        int[][] alleles = getAllelesSortedByFrequency(site);
        int resultSize = alleles[0].length - 1;
        byte[] result = new byte[resultSize];
        for (int i = 0; i < resultSize; i++) {
            result[i] = (byte) alleles[0][i + 1];
        }
        return result;
    }

    @Override
    public int getMajorAlleleCount(int site) {

        int[][] alleles = getAllelesSortedByFrequency(site);

        if (alleles[0].length >= 1) {
            return alleles[1][0];
        } else {
            return 0;
        }

    }

    @Override
    public byte getMajorAllele(int site) {
        int[][] alleles = getAllelesSortedByFrequency(site);

        if (alleles[0].length >= 1) {
            return (byte) alleles[0][0];
        } else {
            return Alignment.UNKNOWN_ALLELE;
        }
    }

    @Override
    public String getMajorAlleleAsString(int site) {
        return genotypeAsString(site, getMajorAllele(site));
    }

    @Override
    public double getMajorAlleleFrequency(int site) {

        int[][] alleles = getAllelesSortedByFrequency(site);

        int numAlleles = alleles[0].length;
        if (numAlleles >= 1) {
            int totalNonMissing = 0;
            for (int i = 0; i < numAlleles; i++) {
                totalNonMissing = totalNonMissing + alleles[1][i];
            }
            return (double) alleles[1][0] / (double) totalNonMissing;
        } else {
            return 0.0;
        }

    }

    @Override
    public double getMinorAlleleFrequency(int site) {

        int[][] alleles = getAllelesSortedByFrequency(site);

        int numAlleles = alleles[0].length;
        if (numAlleles >= 2) {
            int totalNonMissing = 0;
            for (int i = 0; i < numAlleles; i++) {
                totalNonMissing = totalNonMissing + alleles[1][i];
            }
            return (double) alleles[1][1] / (double) totalNonMissing;
        } else {
            return 0.0;
        }

    }

    @Override
    public Object[][] getDiploidsSortedByFrequency(int site) {

        Integer ONE_INTEGER = 1;

        Map<String, Integer> diploidValueCounts = new HashMap<String, Integer>();
        for (int r = 0; r < myTaxaCount; r++) {
            String current = genotypeAsString(r, site);
            Integer num = diploidValueCounts.get(current);
            if (num == null) {
                diploidValueCounts.put(current, ONE_INTEGER);
            } else {
                diploidValueCounts.put(current, ++num);
            }
        }

        Object[][] result = new Object[2][diploidValueCounts.size()];

        int i = 0;
        Iterator<String> itr = diploidValueCounts.keySet().iterator();
        while (itr.hasNext()) {
            String key = itr.next();
            Integer count = diploidValueCounts.get(key);
            result[0][i] = key;
            result[1][i++] = count;
        }

        boolean change = true;
        while (change) {

            change = false;

            for (int k = 0, n = diploidValueCounts.size() - 1; k < n; k++) {

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
    public Object[][] getDiploidCounts() {

        Map<String, Long> diploidValueCounts = new HashMap<String, Long>();
        for (int c = 0; c < mySiteCount; c++) {
            Object[][] diploids = getDiploidsSortedByFrequency(c);
            for (int i = 0; i < diploids[0].length; i++) {
                String current = (String) diploids[0][i];
                Long count = (long) ((Integer) diploids[1][i]).intValue();
                Long num = diploidValueCounts.get(current);
                if (num == null) {
                    diploidValueCounts.put(current, count);
                } else {
                    diploidValueCounts.put(current, (num + count));
                }
            }
        }

        Object[][] result = new Object[2][diploidValueCounts.size()];

        int i = 0;
        Iterator itr = diploidValueCounts.keySet().iterator();
        while (itr.hasNext()) {
            String key = (String) itr.next();
            Long count = diploidValueCounts.get(key);
            result[0][i] = key;
            result[1][i++] = count;
        }

        boolean change = true;
        while (change) {

            change = false;

            for (int k = 0, n = diploidValueCounts.size() - 1; k < n; k++) {

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
    public Object[][] getMajorMinorCounts() {

        String[][] alleleStates = alleleDefinitions();

        if (alleleStates.length != 1) {
            return new Object[0][0];
        }

        long[][] counts = new long[16][16];

        for (int site = 0; site < mySiteCount; site++) {
            byte[] alleles = getAlleles(site);
            if ((alleles == null) || alleles.length == 0) {
                // do nothing
            } else if (alleles.length == 1) {
                counts[alleles[0]][alleles[0]]++;
            } else {
                counts[alleles[0]][alleles[1]]++;
            }
        }

        int numAlleles = 0;
        for (byte x = 0; x < 16; x++) {
            for (byte y = 0; y < 16; y++) {
                if (counts[x][y] != 0) {
                    numAlleles++;
                }
            }
        }

        Object[][] result = new Object[2][numAlleles];
        int nextResult = 0;
        for (byte x = 0; x < 16; x++) {
            for (byte y = 0; y < 16; y++) {
                if (counts[x][y] != 0) {
                    result[0][nextResult] = genotypeAsString(0, x) + ":" + genotypeAsString(0, y);
                    result[1][nextResult++] = counts[x][y];
                }
            }
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
    public int getTotalGametesNotMissingForTaxon(int taxon) {

        int result = 0;
        for (int i = 0, n = mySiteCount; i < n; i++) {
            byte[] current = genotypeArray(taxon, i);
            if (current[0] != Alignment.UNKNOWN_ALLELE) {
                result++;
            }
            if (current[1] != Alignment.UNKNOWN_ALLELE) {
                result++;
            }
        }
        return result;

    }

    @Override
    public int getHeterozygousCountForTaxon(int taxon) {
        int result = 0;
        for (int i = 0, n = mySiteCount; i < n; i++) {
            if (isHeterozygous(taxon, i)) {
                result++;
            }
        }
        return result;
    }

    @Override
    public int getTotalNotMissingForTaxon(int taxon) {

        int result = 0;
        for (int i = 0, n = mySiteCount; i < n; i++) {
            byte current = genotype(taxon, i);
            if (current != Alignment.UNKNOWN_DIPLOID_ALLELE) {
                result++;
            }
        }
        return result;

    }

    @Override
    public byte[] getAlleles(int site) {
        int[][] alleles = getAllelesSortedByFrequency(site);
        int resultSize = alleles[0].length;
        byte[] result = new byte[resultSize];
        for (int i = 0; i < resultSize; i++) {
            result[i] = (byte) alleles[0][i];
        }
        return result;
    }

    @Override
    public int getSiteCount() {
        return mySiteCount;
    }

    @Override
    public int getTaxaCount() {
        return myTaxaCount;
    }

    @Override
    public byte[] getGenotypeForAllSites(int taxon) {
        int numSites = getSiteCount();
        byte[] result = new byte[numSites];
        for (int i = 0; i < numSites; i++) {
            result[i] = genotype(taxon, i);
        }
        return result;
    }

    @Override
    public byte[] getGenotypeForSiteRange(int taxon, int start, int end) {
        int numSites = end - start;
        byte[] result = new byte[numSites];
        for (int i = start; i < end; i++) {
            result[i] = genotype(taxon, i);
        }
        return result;
    }

    @Override
    public byte[] getGenotypeForAllTaxa(int site) {
        int numTaxa = getTaxaCount();
        byte[] result = new byte[numTaxa];
        for (int i = 0; i < numTaxa; i++) {
            result[i] = genotype(i, site);
        }
        return result;
    }
}
