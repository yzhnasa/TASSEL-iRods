// LinkageDisequilibrium.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)
package net.maizegenetics.pal.popgen;

import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.report.TableReport;
import net.maizegenetics.pal.statistics.FisherExact;

import net.maizegenetics.util.ProgressListener;

import java.io.Serializable;
import java.io.StringWriter;

import net.maizegenetics.pal.alignment.SBitAlignment;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;

/**
 * This class calculates D' and r^2 estimates of linkage disequilibrium.  It also
 * calculates the significance of the LD by either Fisher Exact or the multinomial
 * permutation test.  This class can work with either normal alignments of annotated alignments.
 * The alignments should be stripped of invariable numSites.
 *
 * 2 state estimates of D' and r^2 can be found reviewed and discussed in Weir 1996
 *
 * multi-state loci (>=3) require an averaging approach.  In TASSEL 3 in 2010, Buckler
 * removed these approach as the relative magnitudes and meaningfulness of these approaches
 * has never been clear.  Additionally with the moving away from SSR to SNPs these methods are less
 * relevant.  Researchers should convert to biallelic - either by ignoring rarer classes or collapsing
 * rarer states.
 *
 * @version $Id: LinkageDisequilibrium.java,v 1
 *
 * @author Ed Buckler
 */
public class LinkageDisequilibrium extends Thread implements Serializable, TableReport {

    public enum testDesign {

        All, SlidingWindow, SiteByAll
    };
    private Alignment myAlignment;
    private SBitAlignment mySBitAlignment;
    private boolean myRapidPermute = true;
    private int myMinTaxaForEstimate = 2;
    private int myWindowSize = 50;
    private int myTestSite = -1;  // this is only set when one versus all numSites is calculated.
    private int myTotalTests = 0;
    private testDesign myCurrDesign = testDesign.SlidingWindow;
    private float[] myRSqr, myDPrime, myPVal;
    private int[] myiRow, myiCol, mySampleSize;
    private ProgressListener myListener = null;
    private FisherExact myFisherExact;
    private boolean myIsAccumulativeReport = false;
    private int myNumAccumulativeBins = 100;
    private float myAccumulativeInterval;
    private int[] myAccumulativeRValueBins;
    private static String NotImplemented = "NotImplemented";
    private static String NA = "N/A";
    private static Integer IntegerTwo = Integer.valueOf(2);

    
    //
    // I've commented out all the unused constructors.  I will remove
    // once I determine there is no need for them.  -Terry
    //
    
    /**
     * compute LD based on an alignment
     *
     *  @param alignment  Alignment or AnnotationAlignment (this should only contain
     *                    polymorphic numSites)
     *  @param myRapidPermute Use a rapid approach to P-value estimation (see Contigency Table)
     *  @param numberOfPermutations The number of permutations to determine P values
     *  @param windowSize The size of the LD window, determined by user.
     */
    //public LinkageDisequilibrium(Alignment alignment, int numberOfPermutations, int windowSize, boolean rapidPermute, testDesign LDType) {
    //    myAlignment = alignment;
    //    if (myAlignment instanceof SBitAlignment) {
    //        mySBitAlignment = (SBitAlignment) myAlignment;
    //    } else {
    //        mySBitAlignment = SBitAlignment.getInstance(myAlignment, 2, false);
    //    }
    //    myFisherExact = new FisherExact(myAlignment.getSequenceCount() + 10);
    //    myRapidPermute = rapidPermute;
        // this.numberOfPermutations = numberOfPermutations;
    //    myWindowSize = windowSize;
    //    myCurrDesign = LDType;
    //}

    /**
     * Compute LD based on alignment.
     * 
     * @param alignment  Alignment or AnnotationAlignment (this should only contain
     *                    polymorphic numSites)
     * @param myRapidPermute Use a rapid approach to P-value estimation (see Contigency Table)
     * @param numberOfPermutations The number of permutations to determine P values
     * @param windowSize The size of the LD window, determined by user.
     * @param myTestSite
     */
    public LinkageDisequilibrium(Alignment alignment, int numberOfPermutations, int windowSize, boolean rapidPermute, testDesign LDType, int testSite, ProgressListener listener, boolean isAccumulativeReport, int numAccumulateIntervals) {
        myAlignment = alignment;
        if (myAlignment instanceof SBitAlignment) {
            mySBitAlignment = (SBitAlignment) myAlignment;
        } else {
            mySBitAlignment = SBitAlignment.getInstance(myAlignment, 2, false);
        }
        myFisherExact = new FisherExact(myAlignment.getSequenceCount() + 10);
        myRapidPermute = rapidPermute;
        myWindowSize = windowSize;
        myCurrDesign = LDType;
        myTestSite = testSite;
        myListener = listener;
        myIsAccumulativeReport = isAccumulativeReport;
        if (myIsAccumulativeReport) {
            myNumAccumulativeBins = numAccumulateIntervals;
        }
    }

    /**
     * compute LD based on an alignment
     *
     *  @param alignment  Alignment or AnnotationAlignment (this should only contain
     *                    polymorphic numSites)
     *  @param myRapidPermute Use a rapid approach to P-value estimation (see Contigency Table)
     *  @param numberOfPermutations The number of permutations to determine P values
     *   @param myMinTaxaForEstimate After excluding missing data, the minimum number of taxa needed for LD estimation
     */
    //public LinkageDisequilibrium(Alignment alignment, boolean rapidPermute, int numberOfPermutations, int minTaxaForEstimate) {
    //    myAlignment = alignment;
    //    if (myAlignment instanceof SBitAlignment) {
    //        mySBitAlignment = (SBitAlignment) myAlignment;
    //    } else {
    //        mySBitAlignment = SBitAlignment.getInstance(myAlignment, 2, false);
    //    }
    //    myFisherExact = new FisherExact(myAlignment.getSequenceCount() + 10);
    //    myRapidPermute = rapidPermute;
    //    // this.numberOfPermutations = numberOfPermutations;
    //    myMinTaxaForEstimate = minTaxaForEstimate;
    //}

    /**
     * compute LD based on an alignment
     *
     *  @param alignment  Alignment or AnnotationAlignment (this should only contain
     *                    polymorphic numSites)
     *  @param myRapidPermute Use a rapid approach to P-value estimation (see Contigency Table)
     *  @param numberOfPermutations The number of permutations to determine P values
     *   @param myMinTaxaForEstimate After excluding missing data, the minimum number of taxa needed for LD estimation
     */
    //public LinkageDisequilibrium(Alignment alignment, boolean rapidPermute, int numberOfPermutations,
    //        int minTaxaForEstimate, int windowOfSites, testDesign currDesign) {
    //    this(alignment, rapidPermute, numberOfPermutations, minTaxaForEstimate,
    //            windowOfSites, currDesign, -1);
    //}

    /**
     * compute LD based on an alignment
     *
     *  @param alignment  Alignment or AnnotationAlignment (this should only contain
     *                    polymorphic numSites)
     *  @param myRapidPermute Use a rapid approach to P-value estimation (see Contigency Table)
     *  @param numberOfPermutations The number of permutations to determine P values
     *  @param myMinTaxaForEstimate After excluding missing data, the minimum number of taxa needed for LD estimation
     */
    //public LinkageDisequilibrium(Alignment alignment, boolean rapidPermute, int numberOfPermutations,
    //        int minTaxaForEstimate, int windowOfSites, testDesign currDesign, int testSite) {
    //    myAlignment = alignment;
    //    if (myAlignment instanceof SBitAlignment) {
    //        mySBitAlignment = (SBitAlignment) myAlignment;
    //    } else {
    //        mySBitAlignment = SBitAlignment.getInstance(myAlignment, 2, false);
    //    }
    //    myFisherExact = new FisherExact(myAlignment.getSequenceCount() + 10);
    //    myRapidPermute = rapidPermute;
        // this.numberOfPermutations = numberOfPermutations;
    //    myMinTaxaForEstimate = minTaxaForEstimate;
    //    myWindowSize = windowOfSites;
    //    myCurrDesign = currDesign;
    //    myTestSite = testSite;
    //}

    /**
     * starts the thread to calculate LD
     */
    public void run() {
        initMatrices();
        designLDTests();
        calculateBitLDForInbred(true);
    }

    /**
     * Sets of the tests to be done.  Whether the test is done or not is determined by
     * what getIndex returns, which is based on the testDesign
     */
    private void designLDTests() {

        int sites = myAlignment.getSiteCount();

        if (myCurrDesign == testDesign.SlidingWindow) {

            for (int r = 0; r < sites; r++) {
                int start = Math.max(0, r - myWindowSize);
                int end = Math.max(0, r - 1);

                for (int c = start; c <= end; c++) {
                    if (c != r) {
                        int index = getIndex(r, c);
                        if (index > -1) {
                            myiRow[index] = r;
                            myiCol[index] = c;
                        }
                    }
                }
            }

        } else {

            for (int r = 1; r < sites; r++) {
                for (int c = 0; c < r; c++) {
                    int index = getIndex(r, c);
                    if (index > -1) {
                        myiRow[index] = r;
                        myiCol[index] = c;
                    }
                }
            }

        }
    }

    private void initMatrices() {

        int numSites = myAlignment.getSiteCount();
        if (myCurrDesign == testDesign.All) {
            myTotalTests = numSites * (numSites - 1) / 2;
        } else if (myCurrDesign == testDesign.SlidingWindow) {
            int n = Math.min(numSites - 1, myWindowSize);
            myTotalTests = ((n * (n + 1)) / 2) + (numSites - n - 1) * n;
        } else if (myCurrDesign == testDesign.SiteByAll) {
            myTotalTests = numSites;
        }

        if (myIsAccumulativeReport) {
            myAccumulativeInterval = 1.0f / (float) myNumAccumulativeBins;
            myAccumulativeRValueBins = new int[myNumAccumulativeBins + 1];
            myiRow = new int[myTotalTests];
            myiCol = new int[myTotalTests];
        } else {
            myRSqr = new float[myTotalTests];
            myDPrime = new float[myTotalTests];
            myPVal = new float[myTotalTests];
            myiRow = new int[myTotalTests];
            myiCol = new int[myTotalTests];
            mySampleSize = new int[myTotalTests];
        }

    }

    private void calculateBitLDForInbred(boolean collapseMinor) {  //only calculates disequilibrium for inbreds

        int[][] contig;
        for (int currTest = 0; currTest < myTotalTests; currTest++) {
            int r = myiRow[currTest];
            int c = myiCol[currTest];
            int currentProgress = (int) ((double) 100.0 * ((double) currTest / (double) myTotalTests));
            fireProgress(currentProgress);
            contig = new int[2][2];
            BitSet rMj = mySBitAlignment.getAllelePresenceForAllTaxa(r, 0);
            BitSet rMn = mySBitAlignment.getAllelePresenceForAllTaxa(r, 1);
            BitSet cMj = mySBitAlignment.getAllelePresenceForAllTaxa(c, 0);
            BitSet cMn = mySBitAlignment.getAllelePresenceForAllTaxa(c, 1);
            int n = 0;
            n += contig[0][0] = (int) OpenBitSet.intersectionCount(rMj, cMj);
            n += contig[1][0] = (int) OpenBitSet.intersectionCount(rMn, cMj);
            n += contig[0][1] = (int) OpenBitSet.intersectionCount(rMj, cMn);
            n += contig[1][1] = (int) OpenBitSet.intersectionCount(rMn, cMn);

            float rVal = (float) calculateRSqr(contig[0][0], contig[1][0], contig[0][1], contig[1][1], myMinTaxaForEstimate);

            if (myIsAccumulativeReport) {
                if (rVal == Float.NaN) {
                    myAccumulativeRValueBins[myNumAccumulativeBins]++;
                } else {
                    int index = (int) Math.floor(rVal / myAccumulativeInterval);
                    myAccumulativeRValueBins[index]++;
                }
            } else {
                mySampleSize[currTest] = n;
                myRSqr[currTest] = myDPrime[currTest] = myPVal[currTest] = Float.NaN;
                myRSqr[currTest] = rVal;
                myDPrime[currTest] = (float) calculateDPrime(contig[0][0], contig[1][0], contig[0][1], contig[1][1], myMinTaxaForEstimate);
                if (Double.isNaN(myRSqr[currTest]) || Double.isNaN(myDPrime[currTest])) {
                    myPVal[currTest] = Float.NaN;
                } else {
                    myPVal[currTest] = (float) myFisherExact.getTwoTailedP(contig[0][0], contig[1][0], contig[0][1], contig[1][1]);
                }
            }

        } //end of currTest
    }

    private void calculateLDForInbred(boolean collapseMinor) {  //only calculates disequilibrium for inbreds
        int n;
        FisherExact fisherExact = new FisherExact(myAlignment.getSequenceCount() + 10);
        int[][] contig;
        for (int currTest = 0; currTest < myTotalTests; currTest++) {
            int r = myiRow[currTest];
            int c = myiCol[currTest];
            byte rowMajor = (byte) myAlignment.getMajorAllele(r);
            byte rowMinor = (byte) myAlignment.getMinorAllele(r);
            byte colMajor = (byte) myAlignment.getMajorAllele(c);
            byte colMinor = (byte) myAlignment.getMinorAllele(c);
            //double currentProgress = 100 * r * r / (myAlignment.getSiteCount() * myAlignment.getSiteCount());
            int currentProgress = 100 * currTest / myTotalTests;
            fireProgress((int) currentProgress);
            contig = new int[2][2];
            if ((rowMinor == Alignment.UNKNOWN_ALLELE) || (colMinor == Alignment.UNKNOWN_ALLELE)) {
                myRSqr[currTest] = myDPrime[currTest] = myPVal[currTest] = Float.NaN;
                mySampleSize[currTest] = 0;
            } else {
                n = 0;
                for (int sample = 0; sample < myAlignment.getSequenceCount(); sample++) {
                    byte x = myAlignment.getBase(sample, r);
                    byte y = myAlignment.getBase(sample, c);
                    if ((x == Alignment.UNKNOWN_ALLELE) || (y == Alignment.UNKNOWN_ALLELE)) {
                        continue;
                    }
                    int x1, y1;
                    if (x == rowMajor) {
                        x1 = 0;
                    } else if (collapseMinor || (x == rowMinor)) {
                        x1 = 1;
                    } else {
                        continue;
                    }
                    if (y == colMajor) {
                        y1 = 0;
                    } else if (collapseMinor || (y == colMinor)) {
                        y1 = 1;
                    } else {
                        continue;
                    }
                    contig[x1][y1]++;
                    n++;
                } //end of sample
                mySampleSize[currTest] = n;
                myRSqr[currTest] = (float) calculateRSqr(contig[0][0], contig[1][0], contig[0][1], contig[1][1], myMinTaxaForEstimate);
                myDPrime[currTest] = (float) calculateDPrime(contig[0][0], contig[1][0], contig[0][1], contig[1][1], myMinTaxaForEstimate);
                if (Double.isNaN(myRSqr[currTest]) || Double.isNaN(myDPrime[currTest])) {
                    myPVal[currTest] = Float.NaN;
                } else {
                    myPVal[currTest] = (float) fisherExact.getTwoTailedP(contig[0][0], contig[1][0], contig[0][1], contig[1][1]);
                }
            }
        } //end of currTest
    }

    static double calculateDPrime(int countAB, int countAb, int countaB, int countab, int minTaxaForEstimate) {
        //this is the normalized D' is Weir Genetic Data Analysis II 1986 p120
        double freqR, freqC, freq, countR, countC, nonmissingSampleSize;
        nonmissingSampleSize = countAB + countAb + countaB + countab;
        if (nonmissingSampleSize < minTaxaForEstimate) {
            return Double.NaN;
        }
        countR = countab + countAb;
        countC = countab + countaB;
        freqR = (nonmissingSampleSize - countR) / nonmissingSampleSize;
        freqC = (nonmissingSampleSize - countC) / nonmissingSampleSize;
        // if((freqR==0)||(freqC==0)||(freqR==1)||(freqC==1)) return -999;  //changed by ed 8-13-2004
        if ((freqR == 0) || (freqC == 0) || (freqR == 1) || (freqC == 1)) {
            return Double.NaN;
        }
        freq = ((double) countAB / nonmissingSampleSize) - (freqR * freqC);
        if (freq < 0) {
            return freq / Math.max(-freqR * freqC, -(1 - freqR) * (1 - freqC));
        } else {
            return freq / Math.min((1 - freqR) * freqC, (1 - freqC) * freqR);
        }  //check these equations
    }

    static double calculateRSqr(int countAB, int countAb, int countaB, int countab, int minTaxaForEstimate) {
        //this is the Hill & Robertson measure as used in Awadella Science 1999 286:2524
        double freqA, freqB, rsqr, nonmissingSampleSize;
        nonmissingSampleSize = countAB + countAb + countaB + countab;
        if (nonmissingSampleSize < minTaxaForEstimate) {
            return Double.NaN;
        }
        freqA = (double) (countAB + countAb) / nonmissingSampleSize;
        freqB = (double) (countAB + countaB) / nonmissingSampleSize;

        //Through missing data & incomplete datasets some alleles can be fixed this returns missing value
        if ((freqA == 0) || (freqB == 0) || (freqA == 1) || (freqB == 1)) {
            return Double.NaN;
        }

        rsqr = ((double) countAB / nonmissingSampleSize) * ((double) countab / nonmissingSampleSize);
        rsqr -= ((double) countaB / nonmissingSampleSize) * ((double) countAb / nonmissingSampleSize);
        rsqr *= rsqr;
        rsqr /= freqA * (1 - freqA) * freqB * (1 - freqB);
        return rsqr;
    }

    public static double[] getLDForSitePair(Alignment a1, int site1, Alignment a2, int site2,
            int minComp, int minMinor, FisherExact fisherExact) {
        double[] results = null;
        byte rowMajor = (byte) a1.getMajorAllele(site1);
        byte rowMinor = (byte) a1.getMinorAllele(site1);
        byte colMajor = (byte) a2.getMajorAllele(site2);
        byte colMinor = (byte) a2.getMinorAllele(site2);
        int[][] contig = new int[2][2];
        if ((rowMinor == Alignment.UNKNOWN_ALLELE) || (colMinor == Alignment.UNKNOWN_ALLELE)) {
            return null;
        }
        int n = 0;
        for (int sample = 0; sample < a1.getSequenceCount(); sample++) {
            byte x = a1.getBase(sample, site1);
            byte y = a2.getBase(sample, site2);
            if ((x == Alignment.UNKNOWN_ALLELE) || (y == Alignment.UNKNOWN_ALLELE)) {
                continue;
            }
            int x1, y1;
            if (x == rowMajor) {
                x1 = 0;
            } else if (x == rowMinor) {
                x1 = 1;
            } else {
                continue;
            }
            if (y == colMajor) {
                y1 = 0;
            } else if (y == colMinor) {
                y1 = 1;
            } else {
                continue;
            }
            contig[x1][y1]++;
            n++;
        } // end of sample
        // System.out.println(site1+" "+site2+" contig"+Arrays.deepToString(contig));
        if ((n < minComp) || (contig[0][1] + contig[1][1] < minMinor) || (contig[1][0] + contig[1][1] < minMinor)) {
            //          System.out.println("null contig"+Arrays.deepToString(contig));
            return null;
        }
        results = new double[4];
        // System.out.println("contig"+Arrays.deepToString(contig));
        results[0] = n;
        results[1] = (float) calculateRSqr(contig[0][0], contig[1][0], contig[0][1], contig[1][1], minComp);
        results[2] = (float) calculateDPrime(contig[0][0], contig[1][0], contig[0][1], contig[1][1], minComp);
        if (Double.isNaN(results[1]) || Double.isNaN(results[2])) {
            results[3] = Float.NaN;
            // System.out.println("NaN contig"+Arrays.deepToString(contig)+" "+results[0]+" "+results[1]+" "+results[2]);
        } else {
            results[3] = (float) fisherExact.getTwoTailedP(contig[0][0], contig[1][0], contig[0][1], contig[1][1]);
        }
        // if(results[3]<0.00001) System.out.println("contig"+Arrays.deepToString(contig)+" "+results[0]+" "+results[1]+" "+results[2]+" "+results[3]);
        // if(results[3]<0.1) System.out.printf("%d %d %d %d %f \n",contig[0][0], contig[1][0], contig[0][1], contig[1][1], results[3]);
        return results;
    }

    private int getIndex(int r, int c) {
        int index = -1; // -1 is not found
        if (r < c) {
            int t = r;
            r = c;
            c = t;
        }
        if ((c < 0) || (c == r) || (r >= myAlignment.getSiteCount())) {
            return -1;
        }
        if (myCurrDesign == testDesign.All) {
            index = (r * (r - 1) / 2) + c;
        } else if (myCurrDesign == testDesign.SlidingWindow) {
            int sites = myAlignment.getSiteCount();
            int n = Math.min(sites - 1, myWindowSize);
            if ((r - c) <= n) {
                if (r < n) {
                    index = (r * (r - 1) / 2) + c;
                } else {
                    index = (n * (n - 1) / 2) + ((r - n) * n) + (c - r + n);
                }
            }
        } else if (myCurrDesign == testDesign.SiteByAll) {
            if (r == myTestSite) {
                index = c;
            } else if (c == myTestSite) {
                index = r;
            }
        }
        if (index >= myTotalTests) {
            return -1;
        }

        return index;
    }

    /**
     * Returns P-value estimate for a given pair of numSites.  If there were only 2 alleles
     * at each locus, then the Fisher Exact P-value (one-tail) is returned.  If more states
     * then the permutaed Monte Carlo test is used.
     * @param r is site 1
     * @param c is site 2
     * @return P-value
     */
    public double getP(int r, int c) {
        int i = getIndex(r, c);
        if (i > -1) {
            return (double) myPVal[i];
        } else {
            return Double.NaN;
        }
    }

    /**
     * Get number of gametes included in LD calculations (after missing data was excluded)
     * @param r is site 1
     * @param c is site 2
     * @return number of gametes
     *
     */
    private int getN(int r, int c) {
        int i = getIndex(r, c);
        if (i > -1) {
            return mySampleSize[i];
        } else {
            return 0;
        }
    }

    /**
     * Returns D' estimate for a given pair of numSites
     * @param r is site 1
     * @param c is site 2
     * @return D'
     */
    public double getDPrime(int r, int c) {
        int i = getIndex(r, c);
        if (i > -1) {
            return (double) myDPrime[i];
        } else {
            return Double.NaN;
        }
    }

    /**
     * Returns r^2 estimate for a given pair of numSites
     * @param r is site 1
     * @param c is site 2
     * @return D'
     */
    public double getRSqr(int r, int c) {
        int i = getIndex(r, c);
        if (i > -1) {
            return (double) myRSqr[i];
        } else {
            return Double.NaN;
        }
    }

    /**
     * Returns the counts of the numSites in the alignment
     */
    public int getSiteCount() {
        return myAlignment.getSiteCount();
    }

    /**
     * Returns an annotated aligment if one was used for this LD
     * this could be used to access information of locus position
     */
    public Alignment getAlignment() {
        return myAlignment;
    }

    /**
     * Returns representation of the LD results as a string
     */
    public String toString() {
        String delimit = "\t";
        StringWriter sw = new StringWriter();
        Object[] colNames = getTableColumnNames();
        for (int j = 0; j < colNames.length; j++) {
            sw.write(colNames[j].toString());
            sw.write(delimit);
        }
        sw.write("\n");

        for (int r = 0; r < myTotalTests; r++) {
            Object[] theRow = getRow(r);
            for (int i = 0; i < theRow.length; i++) {
                sw.write(theRow[i].toString());
                sw.write(delimit);
            }
        }
        return sw.toString();
    }

    //Implementation of TableReport Interface
    /**
     * @return column names for the table
     */
    public Object[] getTableColumnNames() {
        String[] annotatedLabels = null;
        if (myIsAccumulativeReport) {
            annotatedLabels = new String[]{"R2BinMin", "R2BinMax", "Count"};
        } else {
            annotatedLabels = new String[]{"Locus1", "Position1", "Site1",
                "NumberOfStates1", "States1", "Frequency1", "Locus2", "Position2",
                "Site2", "NumberOfStates2", "States2", "Frequency2", "Dist_bp", "R^2", "DPrime", "pDiseq", "N"};
        }
        return annotatedLabels;
    }

    /**
     * get the data elements
     *
     * @return the data elements
     */
    public Object[][] getTableData() {
        return getTableData(0, getRowCount() - 1);
    }

    /**
     * Get Table Data in specified range inclusive.
     * Indices start at 0.
     *
     * @param start start position
     * @param end end position
     * @return
     */
    public Object[][] getTableData(int start, int end) {
        if ((start < 0) || (end >= getRowCount())) {
            throw new IndexOutOfBoundsException("getTableData: start: " + start + "  end: " + end);
        }
        if (end < start) {
            return null;
        }
        Object[][] temp = new Object[end - start + 1][];
        for (int i = start; i <= end; i++) {
            temp[i] = getRow(i);
        }

        return temp;

    }

    /**
     * Returns specified row.
     *
     * @param row row number
     *
     * @return row
     */
    public Object[] getRow(int row) {

        if (myIsAccumulativeReport) {
            Object[] data = new Object[3];
            if (row == myNumAccumulativeBins) {
                data[0] = Float.NaN;
                data[1] = Float.NaN;
                data[2] = Integer.valueOf(myAccumulativeRValueBins[row]);
            } else {
                float start = myAccumulativeInterval * (float) row;
                data[0] = Float.valueOf(start);
                data[1] = Float.valueOf(start + myAccumulativeInterval);
                data[2] = Integer.valueOf(myAccumulativeRValueBins[row]);
            }
            return data;
        } else {
            int labelOffset = 0;
            Object[] data = new Object[17];

            int r = myiRow[row];
            int c = myiCol[row];

            String rState = (char) myAlignment.getMajorAllele(r) + ":" + (char) myAlignment.getMinorAllele(r);
            Integer rStr = Integer.valueOf(r);

            String cState = (char) myAlignment.getMajorAllele(c) + ":" + (char) myAlignment.getMinorAllele(c);
            Integer cStr = Integer.valueOf(c);

            data[labelOffset++] = myAlignment.getLocusName(r);
            data[labelOffset++] = Integer.valueOf(myAlignment.getPositionInLocus(r));
            data[labelOffset++] = rStr;

            data[labelOffset++] = IntegerTwo;
            data[labelOffset++] = rState;
            data[labelOffset++] = NotImplemented;
            data[labelOffset++] = myAlignment.getLocusName(c);
            data[labelOffset++] = Integer.valueOf(myAlignment.getPositionInLocus(c));
            data[labelOffset++] = cStr;

            data[labelOffset++] = IntegerTwo;
            data[labelOffset++] = cState;
            data[labelOffset++] = NotImplemented;
            if (myAlignment.getLocusName(r).equals(myAlignment.getLocusName(c))) {
                data[labelOffset++] = Integer.valueOf(Math.abs(myAlignment.getPositionInLocus(r) - myAlignment.getPositionInLocus(c)));
            } else {
                data[labelOffset++] = NA;
            }
            data[labelOffset++] = Double.valueOf(getRSqr(r, c));
            data[labelOffset++] = Double.valueOf(getDPrime(r, c));
            data[labelOffset++] = Double.valueOf(getP(r, c));
            data[labelOffset++] = Integer.valueOf(getN(r, c));

            return data;
        }

    }

    /**
     * @return the title of the table
     */
    public String getTableTitle() {
        return "Linkage Disequilibrium";
    }

    // interface ExtendedTableReport
    public int getColumnCount() {
        return getTableColumnNames().length;
    }

    public int getRowCount() {
        if (myIsAccumulativeReport) {
            return myNumAccumulativeBins + 1;
        } else {
            return myTotalTests;
        }
    }

    public int getElementCount() {
        return getRowCount() * getColumnCount();
    }

    public Object getValueAt(int row, int col) {
        return getRow(row)[col];
    }

    protected void fireProgress(int percent) {

        if (percent < 0) {
            percent = -percent;
        }

        if (myListener != null) {
            myListener.progress(percent, null);
        }

    }
}
