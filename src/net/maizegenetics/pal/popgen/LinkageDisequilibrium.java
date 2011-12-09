// LinkageDisequilibrium.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)
package net.maizegenetics.pal.popgen;

import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.datatype.DataType;
import net.maizegenetics.pal.report.TableReport;
import net.maizegenetics.pal.statistics.FisherExact;

import net.maizegenetics.util.ProgressListener;

import java.io.Serializable;
import java.io.StringWriter;

/**
 * This class calculates D' and r^2 estimates of linkage disequilibrium.  It also
 * calculates the significance of the LD by either Fisher Exact or the multinomial
 * permutation test.  This class can work with either normal alignments of annotated alignments.
 * The alignments should be stripped of invariable sites.
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
    protected Alignment theAlignment;
    boolean rapidPermute = true;
//    int numberOfPermutations = 1000;
    int minTaxaForEstimate = 2;
    int windowOfSites = 50;
    int testSite = -1;  //this is only set when one versus all sites is calculated.
    int totalTests = 0;
    testDesign currDesign = testDesign.SlidingWindow;
    float[] rsqr, dprime, pval; //4*3 = 12 bytes per entry
    int[] irow, icol, sampleSize;  //4*3 = 12 bytes per entry
    boolean isDense = true;
    private ProgressListener myListener = null;

    /**
     * compute LD based on an alignment
     *
     *  @param alignment  Alignment or AnnotationAlignment (this should only contain
     *                    polymorphic sites)
     *  @param rapidPermute Use a rapid approach to P-value estimation (see Contigency Table)
     *  @param numberOfPermutations The number of permutations to determine P values
     *  @param windowSize The size of the LD window, determined by user.
     */
    public LinkageDisequilibrium(Alignment alignment, int numberOfPermutations, int windowSize, boolean rapidPermute, testDesign LDType) {
        this.theAlignment = alignment;
        this.rapidPermute = rapidPermute;
        // this.numberOfPermutations = numberOfPermutations;
        this.windowOfSites = windowSize;
        this.currDesign = LDType;
    }

    /**
     * Compute LD based on alignment.
     * 
     * @param alignment  Alignment or AnnotationAlignment (this should only contain
     *                    polymorphic sites)
     * @param rapidPermute Use a rapid approach to P-value estimation (see Contigency Table)
     * @param numberOfPermutations The number of permutations to determine P values
     * @param windowSize The size of the LD window, determined by user.
     * @param testSite
     */
    public LinkageDisequilibrium(Alignment alignment, int numberOfPermutations, int windowSize, boolean rapidPermute, testDesign LDType, int testSite, ProgressListener listener) {
        this.theAlignment = alignment;
        this.rapidPermute = rapidPermute;
        this.windowOfSites = windowSize;
        this.currDesign = LDType;
        this.testSite = testSite;
        myListener = listener;
    }

    /**
     * compute LD based on an alignment
     *
     *  @param alignment  Alignment or AnnotationAlignment (this should only contain
     *                    polymorphic sites)
     *  @param rapidPermute Use a rapid approach to P-value estimation (see Contigency Table)
     *  @param numberOfPermutations The number of permutations to determine P values
     *   @param minTaxaForEstimate After excluding missing data, the minimum number of taxa needed for LD estimation
     */
    public LinkageDisequilibrium(Alignment alignment, boolean rapidPermute, int numberOfPermutations, int minTaxaForEstimate) {
        this.theAlignment = alignment;
        this.rapidPermute = rapidPermute;
        // this.numberOfPermutations = numberOfPermutations;
        this.minTaxaForEstimate = minTaxaForEstimate;
    }

    /**
     * compute LD based on an alignment
     *
     *  @param alignment  Alignment or AnnotationAlignment (this should only contain
     *                    polymorphic sites)
     *  @param rapidPermute Use a rapid approach to P-value estimation (see Contigency Table)
     *  @param numberOfPermutations The number of permutations to determine P values
     *   @param minTaxaForEstimate After excluding missing data, the minimum number of taxa needed for LD estimation
     */
    public LinkageDisequilibrium(Alignment alignment, boolean rapidPermute, int numberOfPermutations,
            int minTaxaForEstimate, int windowOfSites, testDesign currDesign) {
        this(alignment, rapidPermute, numberOfPermutations, minTaxaForEstimate,
                windowOfSites, currDesign, -1);
    }

    /**
     * compute LD based on an alignment
     *
     *  @param alignment  Alignment or AnnotationAlignment (this should only contain
     *                    polymorphic sites)
     *  @param rapidPermute Use a rapid approach to P-value estimation (see Contigency Table)
     *  @param numberOfPermutations The number of permutations to determine P values
     *  @param minTaxaForEstimate After excluding missing data, the minimum number of taxa needed for LD estimation
     */
    public LinkageDisequilibrium(Alignment alignment, boolean rapidPermute, int numberOfPermutations,
            int minTaxaForEstimate, int windowOfSites, testDesign currDesign, int testSite) {
        this.theAlignment = alignment;
        this.rapidPermute = rapidPermute;
        // this.numberOfPermutations = numberOfPermutations;
        this.minTaxaForEstimate = minTaxaForEstimate;
        this.windowOfSites = windowOfSites;
        this.currDesign = currDesign;
        this.testSite = testSite;
    }

    /**
     * starts the thread to calculate LD
     */
    public void run() {
        initMatrices();
        designLDTests();
        calculateLDForInbred(true);
    }

    /**
     * Sets of the tests to be done.  Whether the test is done or not is determined by
     * what getIndex returns, which is based on the testDesign
     */
    private void designLDTests() {

        int sites = theAlignment.getSiteCount();

        if (currDesign == testDesign.SlidingWindow) {

            for (int r = 0; r < sites; r++) {
                //int half = Math.round((float)windowOfSites / 2.0f);
                int start = Math.max(0, r - windowOfSites);
                int end = Math.max(0, r - 1);

                for (int c = start; c <= end; c++) {
                    if (c != r) {
                        int index = getIndex(r, c);
                        if (index > -1) {
                            irow[index] = r;
                            icol[index] = c;
                        }
                    }
                }
            }

        } else {

            for (int r = 1; r < sites; r++) {
                for (int c = 0; c < r; c++) {
                    int index = getIndex(r, c);
                    if (index > -1) {
                        irow[index] = r;
                        icol[index] = c;
                    }
                }
            }

        }
    }

    private void initMatrices() {
        int sites = theAlignment.getSiteCount();
        if (currDesign == testDesign.All) {
            totalTests = sites * (sites - 1) / 2;
        } else if (currDesign == testDesign.SlidingWindow) {
            //int half = Math.round((float)windowOfSites / 2.0f);
            int n = Math.min(sites - 1, windowOfSites);
            totalTests = ((n * (n + 1)) / 2) + (sites - n - 1) * n;
        } else if (currDesign == testDesign.SiteByAll) {
            totalTests = theAlignment.getSiteCount();
        }
        rsqr = new float[totalTests];
        dprime = new float[totalTests];
        pval = new float[totalTests];
        irow = new int[totalTests];
        icol = new int[totalTests];
        sampleSize = new int[totalTests];
    }

    private void calculateLDForInbred(boolean collapseMinor) {  //only calculates disequilibrium for inbreds
        int n;
        FisherExact fisherExact = new FisherExact(theAlignment.getSequenceCount() + 10);
        int[][] contig;
        for (int currTest = 0; currTest < totalTests; currTest++) {
            int r = irow[currTest];
            int c = icol[currTest];
            byte rowMajor = (byte) theAlignment.getMajorAllele(r);
            byte rowMinor = (byte) theAlignment.getMinorAllele(r);
            byte colMajor = (byte) theAlignment.getMajorAllele(c);
            byte colMinor = (byte) theAlignment.getMinorAllele(c);
            //double currentProgress = 100 * r * r / (theAlignment.getSiteCount() * theAlignment.getSiteCount());
            int currentProgress = 100 * currTest / totalTests;
            fireProgress((int) currentProgress);
            contig = new int[2][2];
            if ((rowMinor == DataType.UNKNOWN_BYTE) || (colMinor == DataType.UNKNOWN_BYTE)) {
                rsqr[currTest] = dprime[currTest] = pval[currTest] = Float.NaN;
                sampleSize[currTest] = 0;
            } else {
                n = 0;
                for (int sample = 0; sample < theAlignment.getSequenceCount(); sample++) {
                    byte x = theAlignment.getBase(sample, r);
                    byte y = theAlignment.getBase(sample, c);
                    if ((x == DataType.UNKNOWN_BYTE) || (y == DataType.UNKNOWN_BYTE)) {
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
                }//end of sample
                sampleSize[currTest] = n;
                rsqr[currTest] = (float) calculateRSqr(contig[0][0], contig[1][0], contig[0][1], contig[1][1], minTaxaForEstimate);
                dprime[currTest] = (float) calculateDPrime(contig[0][0], contig[1][0], contig[0][1], contig[1][1], minTaxaForEstimate);
                if (Double.isNaN(rsqr[currTest]) || Double.isNaN(dprime[currTest])) {
                    pval[currTest] = Float.NaN;
                } else {
                    pval[currTest] = (float) fisherExact.getTwoTailedP(contig[0][0], contig[1][0], contig[0][1], contig[1][1]);
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
        if ((rowMinor == DataType.UNKNOWN_BYTE) || (colMinor == DataType.UNKNOWN_BYTE)) {
            return null;
        }
        int n = 0;
        for (int sample = 0; sample < a1.getSequenceCount(); sample++) {
            byte x = a1.getBase(sample, site1);
            byte y = a2.getBase(sample, site2);
            if ((x == DataType.UNKNOWN_BYTE) || (y == DataType.UNKNOWN_BYTE)) {
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
        }//end of sample
//        System.out.println(site1+" "+site2+" contig"+Arrays.deepToString(contig));
        if ((n < minComp) || (contig[0][1] + contig[1][1] < minMinor) || (contig[1][0] + contig[1][1] < minMinor)) {
            //          System.out.println("null contig"+Arrays.deepToString(contig));
            return null;
        }
        results = new double[4];
//        System.out.println("contig"+Arrays.deepToString(contig));
        results[0] = n;
        results[1] = (float) calculateRSqr(contig[0][0], contig[1][0], contig[0][1], contig[1][1], minComp);
        results[2] = (float) calculateDPrime(contig[0][0], contig[1][0], contig[0][1], contig[1][1], minComp);
        if (Double.isNaN(results[1]) || Double.isNaN(results[2])) {
            results[3] = Float.NaN;
//            System.out.println("NaN contig"+Arrays.deepToString(contig)+" "+results[0]+" "+results[1]+" "+results[2]);
        } else {
            results[3] = (float) fisherExact.getTwoTailedP(contig[0][0], contig[1][0], contig[0][1], contig[1][1]);
        }
//        if(results[3]<0.00001) System.out.println("contig"+Arrays.deepToString(contig)+" "+results[0]+" "+results[1]+" "+results[2]+" "+results[3]);
//        if(results[3]<0.1) System.out.printf("%d %d %d %d %f \n",contig[0][0], contig[1][0], contig[0][1], contig[1][1], results[3]);
        return results;
    }

    private int getIndex(int r, int c) {
        int index = -1; //-1 is not found
//        if((r==c)||(r<0)||(c<0)||(r>=theAlignment.getSiteCount())||(c>=theAlignment.getSiteCount())) return -1;
        if (r < c) {
            int t = r;
            r = c;
            c = t;
        }
        if ((c < 0) || (c == r) || (r >= theAlignment.getSiteCount())) {
            return -1;
        }
        if (currDesign == testDesign.All) {
            index = (r * (r - 1) / 2) + c;
        } else if (currDesign == testDesign.SlidingWindow) {
            int sites = theAlignment.getSiteCount();
            //int half = Math.round((float)windowOfSites / 2.0f);
            int n = Math.min(sites - 1, windowOfSites);
            if ((r - c) <= n) {
                if (r < n) {
                    index = (r * (r - 1) / 2) + c;
                } else {
                    index = (n * (n - 1) / 2) + ((r - n) * n) + (c - r + n);
                }
            }
        } else if (currDesign == testDesign.SiteByAll) {
            if (r == testSite) {
                index = c;
            } else if (c == testSite) {
                index = r;
            }
        }
        if (index >= totalTests) {
            return -1;
        }

        return index;
    }

    /** Returns P-value estimate for a given pair of sites.  If there were only 2 alleles
     *  at each locus, then the Fisher Exact P-value (one-tail) is returned.  If more states
     *  then the permutaed Monte Carlo test is used.
     *  @param r is site 1
     *  @param c is site 2
     *  @return P-value
     */
    public double getP(int r, int c) {
        int i = getIndex(r, c);
        if (i > -1) {
            return (double) pval[i];
        } else {
            return Double.NaN;
        }
    }

    /** Get number of gametes included in LD calculations (after missing data was excluded)
     *  @param r is site 1
     *  @param c is site 2
     *  @return number of gametes
     *
     */
    private int getN(int r, int c) {
        int i = getIndex(r, c);
        if (i > -1) {
            return sampleSize[i];
        } else {
            return 0;
        }
    }

    /** Returns D' estimate for a given pair of sites
     *  @param r is site 1
     *  @param c is site 2
     *  @return D'
     */
    public double getDPrime(int r, int c) {
        int i = getIndex(r, c);
        if (i > -1) {
            return (double) dprime[i];
        } else {
            return Double.NaN;
        }
    }

    /** Returns r^2 estimate for a given pair of sites
     *  @param r is site 1
     *  @param c is site 2
     *  @return D'
     */
    public double getRSqr(int r, int c) {
        int i = getIndex(r, c);
        if (i > -1) {
            return (double) rsqr[i];
        } else {
            return Double.NaN;
        }
    }

    /**
     * Returns the counts of the sites in the alignment
     */
    public int getSiteCount() {
        return theAlignment.getSiteCount();
    }

    /**
     * Returns an annotated aligment if one was used for this LD
     * this could be used to access information of locus position
     */
    public Alignment getAnnotatedAlignment() {
        return theAlignment;
    }

    /** returns representation of the LD results as a string */
    public String toString() {
        String delimit = "\t";
        StringWriter sw = new StringWriter();
        Object[] colNames = getTableColumnNames();
        for (int j = 0; j < colNames.length; j++) {
            sw.write(colNames[j].toString());
            sw.write(delimit);
        }
        sw.write("\n");

        for (int r = 0; r < totalTests; r++) {
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
        String[] annotatedLabels = {"Locus1", "Position1", "Site1",
            "NumberOfStates1", "States1", "Frequency1", "Locus2", "Position2",
            "Site2", "NumberOfStates2", "States2", "Frequency2", "Dist_bp", "R^2", "DPrime", "pDiseq", "N"};
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

        String NotImplemented = "NotImplemented";
        String NA = "N/A";
        Object[] data;
        java.text.NumberFormat nf = new java.text.DecimalFormat();
        nf.setMaximumFractionDigits(8);
        int labelOffset = 0;
        data = new String[17];

        int r = irow[row];
        int c = icol[row];

        String rState = (char) theAlignment.getMajorAllele(r) + ":" + theAlignment.getMinorAllele(r);
        //String rState = theAlignment.getAlleles(c)
        String rStr = String.valueOf(r);

        String cState = (char) theAlignment.getMajorAllele(c) + ":" + theAlignment.getMinorAllele(c);
        String cStr = String.valueOf(c);

        data[labelOffset++] = theAlignment.getLocusName(r);
        data[labelOffset++] = String.valueOf(theAlignment.getPositionInLocus(r));
        data[labelOffset++] = rStr;

        data[labelOffset++] = String.valueOf(2);
        data[labelOffset++] = rState;
        data[labelOffset++] = NotImplemented;
        data[labelOffset++] = theAlignment.getLocusName(c);
        data[labelOffset++] = String.valueOf(theAlignment.getPositionInLocus(c));
        data[labelOffset++] = cStr;

        data[labelOffset++] = String.valueOf(2);
        data[labelOffset++] = cState;
        data[labelOffset++] = NotImplemented;
        if (theAlignment.getLocusName(r).equals(theAlignment.getLocusName(c))) {
            data[labelOffset++] = String.valueOf(Math.abs(theAlignment.getPositionInLocus(r) - theAlignment.getPositionInLocus(c)));
        } else {
            data[labelOffset++] = NA;
        }
        data[labelOffset++] = nf.format(getRSqr(r, c));
        data[labelOffset++] = nf.format(getDPrime(r, c));
        data[labelOffset++] = nf.format(getP(r, c));
        data[labelOffset++] = String.valueOf(getN(r, c));

        return data;

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
        return this.totalTests;
    }

    public int getElementCount() {
        return getRowCount() * getColumnCount();
    }

    /**
     * Remove the rows passed
     * This method is not implemented as the primary reason for implementing ExtendedTableReport
     * is to 
     *
     * @param rows
     */
    public void deleteRows(int[] rows) {
        System.out.println("LinkageDisequilibrium.deleteRows(int[] rows is not implemented)");
    }

    /**
     * Remove the columns passed
     *
     * @param cols
     */
    public void deleteColumns(int[] cols) {
        System.out.println("LinkageDisequilibrium.deleteColumns(int[] cols) is not implemented");
    }

    public Object getValueAt(int row, int col) {
        return getRow(row)[col];
    }

    protected void fireProgress(int percent) {

        if (myListener != null) {
            myListener.progress(percent, null);
        }

    }
}
