/*
 *  NucleotideGenotype
 */
package net.maizegenetics.pal.alignment.genotype;

import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Map;

import net.maizegenetics.pal.alignment.NucleotideAlignmentConstants;
import net.maizegenetics.util.SuperByteMatrix;

/**
 *
 * @author Terry Casstevens
 */
public class NucleotideGenotype extends ByteGenotype {

    private static final int SHIFT_AMOUNT = 10;
    /**
     * Byte representations of DNA sequences are stored in blocks of 65536 sites
     */
    private static final int NUM_SITES_TO_CACHE = 1 << SHIFT_AMOUNT;
    public static final int SITE_BLOCK_MASK = ~(NUM_SITES_TO_CACHE - 1);
    private static final int MAX_CACHE_SIZE = 50;
    private final Map<Integer, SiteStats> myCachedSites = new LinkedHashMap<Integer, SiteStats>((3 * MAX_CACHE_SIZE) / 2) {
        @Override
        protected boolean removeEldestEntry(Map.Entry eldest) {
            return size() > MAX_CACHE_SIZE;
        }
    };
    private final int myMaxNumThreads = Runtime.getRuntime().availableProcessors();
    private int myNumRunningThreads = 0;

    private static int getStartSite(int site) {
        return site & SITE_BLOCK_MASK;
    }

    NucleotideGenotype(SuperByteMatrix genotype, boolean phased) {
        super(genotype, phased, NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES);
    }

    private SiteStats getCachedSiteStats(int site) {
        int startSite = getStartSite(site);
        SiteStats result = myCachedSites.get(startSite);
        if (result == null) {
            if (myNumRunningThreads < myMaxNumThreads) {
                myNumRunningThreads++;
                new Thread(new LookAheadSiteStats(startSite + NUM_SITES_TO_CACHE * 10)).start();
            }
            result = calculateAlleleFreq(startSite);
            myCachedSites.put(startSite, result);
        }

        if ((site == startSite) && (myNumRunningThreads < myMaxNumThreads)) {
            myNumRunningThreads++;
            new Thread(new LookAheadSiteStats(startSite + NUM_SITES_TO_CACHE * 20)).start();
        }
        return result;
    }

    @Override
    public int[][] getAllelesSortedByFrequency(int site) {
        return getCachedSiteStats(site).getAllelesSortedByFrequency(site);
    }

    private SiteStats calculateAlleleFreq(int startSite) {
        int numSites = Math.min(NUM_SITES_TO_CACHE, getSiteCount() - startSite);
        int numTaxa = getTaxaCount();
        int[][] alleleFreq = new int[numSites][6];
        for (int taxon = 0; taxon < numTaxa; taxon++) {
            for (int s = 0; s < numSites; s++) {
                byte[] b = getBaseArray(taxon, s + startSite);
                if (b[0] < 6) {
                    alleleFreq[s][b[0]]++;
                }
                if (b[1] < 6) {
                    alleleFreq[s][b[1]]++;
                }
            }
        }

        int[][][] alleleCounts = new int[numSites][2][];
        for (int s = 0; s < numSites; s++) {
            int[] cntAndAllele = new int[6];
            for (byte i = 0; i < 6; i++) {
                cntAndAllele[i] = (alleleFreq[s][i] << 4) | (5 - i);  //size | allele (the 5-i is to get the sort right, so if case of ties A is first)
            }
            Arrays.sort(cntAndAllele);  //ascending quick sort
            int numAlleles = 0;
            for (byte i = 5; i >= 0; i--) {
                if (cntAndAllele[i] <= 0xF) {
                    numAlleles = 5 - i;
                    break;
                }
            }
            alleleCounts[s][0] = new int[numAlleles];
            alleleCounts[s][1] = new int[numAlleles];
            for (int i = 5, n = 5 - numAlleles; i > n; i--) {
                alleleCounts[s][0][5 - i] = (byte) (5 - (0xF & cntAndAllele[i]));
                alleleCounts[s][1][5 - i] = alleleFreq[s][alleleCounts[s][0][5 - i]];
            }
        }

        return new SiteStats(startSite, alleleCounts);
    }

    private class LookAheadSiteStats implements Runnable {

        private final int myStartSite;

        public LookAheadSiteStats(int site) {
            myStartSite = getStartSite(site);
        }

        @Override
        public void run() {
            try {
                if (myStartSite >= mySiteCount) {
                    return;
                }
                SiteStats result = myCachedSites.get(myStartSite);
                if (result == null) {
                    result = calculateAlleleFreq(myStartSite);
                    myCachedSites.put(myStartSite, result);
                }
            } finally {
                myNumRunningThreads--;
            }
        }
    }

    private class SiteStats {

        private final int myStartSite;
        private final int[][][] myAlleleCounts;

        SiteStats(int startSite, int[][][] alleleCounts) {
            myStartSite = startSite;
            myAlleleCounts = alleleCounts;
        }

        public int[][] getAllelesSortedByFrequency(int site) {
            return myAlleleCounts[site - myStartSite];
        }
    }
}
