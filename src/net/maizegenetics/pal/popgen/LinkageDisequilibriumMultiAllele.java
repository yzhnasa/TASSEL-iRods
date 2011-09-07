// LinkageDisequilibriumMultiAllele.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)
package net.maizegenetics.pal.popgen;

import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.datatype.DataType;
import net.maizegenetics.pal.report.TableReport;
import net.maizegenetics.pal.statistics.ContigencyTable;
import net.maizegenetics.pal.statistics.FisherExact;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.io.Serializable;
import java.io.StringWriter;

import java.util.Vector;

/**
 * This class calculates D' and r^2 estimates of linkage disequilibrium.  It also
 * calculates the significance of the LD by either Fisher Exact or the multinomial
 * permutation test.  This class can work with either normal alignments of annotated alignments.
 * The alignments should be stripped of invariable sites.
 *
 * 2 state estimates of D' and r^2 can be found reviewed and discussed in Weir 1996
 *
 * multi-state loci (>=3) require an averaging approach.  These should not be used for
 * popgen parameter estimates, unless you know specifically that it works for multistate loci.
 * The estimate of D' is the approach used by Farnir 2000 Genome Research 10:220-227
 * that follows Hedrick 1987.  r^2 was estimated in a similar way.
 *
 * @version $Id: LinkageDisequilibriumMultiAllele.java,v 1
 *
 * @author Ed Buckler
 */
public class LinkageDisequilibriumMultiAllele extends Thread implements Serializable, TableReport {
    /*this class converts the alignment into a matrix for internal use, as this is 2.5 times faster than
    making lots of calls to the original matrix
     */

    protected Alignment theAlignment;
    Vector[] stateVector;
    boolean rapidPermute = true;
    int numberOfPermutations = 1000;
    int minTaxaForEstimate = 2;
    double[][] diseq, pDiseq;
    private double currentProgress;

    /**
     * compute LD based on an alignment.  The default is to used used rapid permutations
     * that provides slightly biased P-values, and 1000 permutations to evaluate P-values.
     *
     *  @param alignment  Alignment or AnnotationAlignment (this should only contain
     *                    polymorphic sites)
     */
    public LinkageDisequilibriumMultiAllele(Alignment alignment) {
        this.theAlignment = alignment;
    }

    /**
     * compute LD based on an alignment
     *
     *  @param alignment  Alignment or AnnotationAlignment (this should only contain
     *                    polymorphic sites)
     *  @param rapidPermute Use a rapid approach to P-value estimation (see Contigency Table)
     *  @param numberOfPermutations The number of permutations to determine P values
     */
    public LinkageDisequilibriumMultiAllele(Alignment alignment, boolean rapidPermute, int numberOfPermutations) {
        this.theAlignment = alignment;
        this.rapidPermute = rapidPermute;
        this.numberOfPermutations = numberOfPermutations;
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
    public LinkageDisequilibriumMultiAllele(Alignment alignment, boolean rapidPermute, int numberOfPermutations, int minTaxaForEstimate) {
        this.theAlignment = alignment;
        this.rapidPermute = rapidPermute;
        this.numberOfPermutations = numberOfPermutations;
        this.minTaxaForEstimate = minTaxaForEstimate;
    }

    /**
     * compute LD based on an site pattern, needs to be implemented
     *
     * @param sp site pattern
     */
    //  public LinkageDisequilibriumMultiAllele(SitePattern sp) {
    //  }
    /**
     * starts the thread to calculate LD
     */
    public void run() {
        // Date theDate = new Date();
        // System.out.println("Start LD run now!");
        calculateMultiProbDiseq();
        // Date stopDate = new Date();
        // System.out.println("Stop LD run Time="+(theDate.getTime()-stopDate.getTime()));
    }

    /**
     *Determines the number of states for each site, and creates a byte matrix will all
     *the states renumbered from 0 to states-1.  It also sets GAP and UNKNOWN_CHARACTERS to
     *-99.
     *
     * @return matrix of renumbered states
     */
    byte[][] determineNumberOfStates() {
        // System.out.println("States starting to be loaded into S");
        stateVector = new Vector[theAlignment.getSiteCount()];
        byte[][] S = new byte[theAlignment.getSequenceCount()][theAlignment.getSiteCount()];
        char c;
        for (int i = 0; i < theAlignment.getSiteCount(); i++) {
            stateVector[i] = new Vector();
            for (int j = 0; j < theAlignment.getSequenceCount(); j++) {
                c = theAlignment.getBaseChar(j, i);
                Character cObj = Character.valueOf(c);
                DataType dt;
                dt = theAlignment.getDataType();
                //if ((c == DataType.UNKNOWN_CHARACTER || dt.isUnknownChar(c))) {
                if (c == DataType.UNKNOWN_CHARACTER) {
                    S[j][i] = -99;
                } else if (!stateVector[i].contains(cObj)) {
                    stateVector[i].add(cObj);
                    S[j][i] = (byte) stateVector[i].indexOf(cObj);
                } else {
                    S[j][i] = (byte) stateVector[i].indexOf(cObj);
                }
            }
        }
        // System.out.println("States loaded into S");
        return S;
    }

    private void calculateMultiProbDiseq() {  //only calculates disequilibrium for inbreds
        byte S[][];
        S = this.determineNumberOfStates();
        int rows, cols, n;
        ContigencyTable contigencyTable = new ContigencyTable(theAlignment.getSequenceCount() + 10);
        FisherExact fisherExact = new FisherExact(theAlignment.getSequenceCount() + 10);

        int[][] contig;

        diseq = new double[theAlignment.getSiteCount()][theAlignment.getSiteCount()];
        pDiseq = new double[theAlignment.getSiteCount()][theAlignment.getSiteCount()];

        for (int r = 0; r < theAlignment.getSiteCount(); r++) {
            currentProgress = 100 * r * r / (theAlignment.getSiteCount() * theAlignment.getSiteCount());
            rows = stateVector[r].size();
            // System.out.println("r="+rows);
            pDiseq[r][r] = 0.0;
            diseq[r][r] = 1.0;
            // System.out.println(r+"'s states="+rows);
            for (int c = 0; c < r; c++) {
                cols = stateVector[c].size();
                contig = new int[rows][cols];
                n = 0;
                for (int sample = 0; sample < theAlignment.getSequenceCount(); sample++) {//rChar=theAlignment.getData(sample,r);
                    //cChar=theAlignment.getData(sample,c);
                    if ((S[sample][r] != -99) && (S[sample][c] != -99)) {//System.out.println("sample="+sample+"  S[sample][r]="+S[sample][r]+"  S[sample][c]="+S[sample][c]);
                        contig[S[sample][r]][S[sample][c]]++;
                        n++;
                    }
                }//end of sample
                //System.out.println("r="+r+" c="+c+" AA="+contig[0][0]+" Aa="+contig[0][1]+" aA="+contig[1][0]+" aa="+contig[1][1]);
                if ((rows == 2) && (cols == 2)) {
                    diseq[r][c] = calculateRSqrDisequilibrium(contig[0][0], contig[1][0], contig[0][1], contig[1][1]);
                    diseq[c][r] = calculateDisequilibrium(contig[0][0], contig[1][0], contig[0][1], contig[1][1]);
                    pDiseq[r][c] = n;
                    if (Double.isNaN(diseq[r][c]) || Double.isNaN(diseq[r][c])) {
                        pDiseq[c][r] = Double.NaN;
                    } else {
                        pDiseq[c][r] = fisherExact.getTwoTailedP(contig[0][0], contig[1][0], contig[0][1], contig[1][1]);
                    }
                } else if ((rows >= 2) && (cols >= 2)) //multiallelic, then use a permuted contigency table test
                {
                    diseq[r][c] = calculateRSqrDisequilibrium(contig, rows, cols);
                    diseq[c][r] = calculateDisequilibrium(contig, rows, cols);
                    contigencyTable.setMatrix(contig);
                    pDiseq[r][c] = n;
                    if (Double.isNaN(diseq[r][c]) || Double.isNaN(diseq[r][c])) {
                        pDiseq[c][r] = Double.NaN;
                    } else {
                        if (this.rapidPermute) {
                            pDiseq[c][r] = contigencyTable.calcRapidMonteCarloExactTest(numberOfPermutations);
                        } else {
                            pDiseq[c][r] = contigencyTable.calcMonteCarloExactTest(numberOfPermutations);
                        }
                    }
                } else //One locus must be monomorphic, do not calculate LD
                {
                    diseq[c][r] = Double.NaN;
                    diseq[r][c] = Double.NaN;
                    pDiseq[r][c] = n;
                    pDiseq[c][r] = Double.NaN;
                }
            // System.out.println("pDiseq["+r+"]["+c+"]="+pDiseq[r][c]);
            } //end of c
        } //end of r
    }

    double calculateDisequilibrium(int countAA, int countAa, int countaA, int countaa) {
        //this is the normalized D' is Weir Genetic Data Analysis II 1986 p120
        double freqR, freqC, freq, countR, countC, nonmissingSampleSize;
        nonmissingSampleSize = countAA + countAa + countaA + countaa;
        if (nonmissingSampleSize <= minTaxaForEstimate) {
            return Double.NaN;
        }
        countR = countaa + countAa;
        countC = countaa + countaA;
        freqR = (nonmissingSampleSize - countR) / nonmissingSampleSize;
        freqC = (nonmissingSampleSize - countC) / nonmissingSampleSize;
        // if((freqR==0)||(freqC==0)||(freqR==1)||(freqC==1)) return -999;  //changed by ed 8-13-2004
        if ((freqR == 0) || (freqC == 0) || (freqR == 1) || (freqC == 1)) {
            return Double.NaN;
        }
        freq = ((double) countAA / nonmissingSampleSize) - (freqR * freqC);
        if (freq < 0) {
            return freq / Math.max(-freqR * freqC, -(1 - freqR) * (1 - freqC));
        } else {
            return freq / Math.min((1 - freqR) * freqC, (1 - freqC) * freqR);
        }  //check these equations
    }

    double calculateRSqrDisequilibrium(int countAB, int countAb, int countaB, int countab) {
        //this is the Hill & Robertson measure as used in Awadella Science 1999 286:2524
        double freqA, freqB, rsqr, nonmissingSampleSize;
        nonmissingSampleSize = countAB + countAb + countaB + countab;
        if (nonmissingSampleSize <= minTaxaForEstimate) {
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

    double calculateDisequilibrium(int[][] contig, int rows, int cols) {
        //this is the D' approach used by Farnir 2000 Genome Research 10:220-227
        //this follow Hedrick 1987
        double Dij, Dmax, D_prime_ij, D_prime = 0, pi, qj, nonmissingSampleSize = 0;
        double[] p_margin = new double[rows];
        double[] q_margin = new double[cols];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                nonmissingSampleSize += contig[i][j];
                p_margin[i] += contig[i][j];
                q_margin[j] += contig[i][j];
            }
        }
        if (nonmissingSampleSize <= minTaxaForEstimate) {
            return Double.NaN;
        }
        for (int i = 0; i < rows; i++) {//Through missing data & incomplete datasets some alleles can be fixed this return missing value
            //       if(p_margin[i]==0) continue;   //changed by ed 8-13-2004
            if (p_margin[i] == 0) {
                return Double.NaN;
            }
            pi = p_margin[i] / nonmissingSampleSize;
            for (int j = 0; j < cols; j++) {//Through missing data & incomplete datasets some alleles can be fixed this return missing value
                // if(q_margin[j]==0) continue;   //changed by ed 8-13-2004
                if (q_margin[j] == 0) {
                    return Double.NaN;
                }
                qj = q_margin[j] / nonmissingSampleSize;
                Dij = (contig[i][j] / nonmissingSampleSize) - (pi * qj);
                if (Dij < 0) {
                    Dmax = Math.min(pi * qj, (1 - pi) * (1 - qj));
                } else {
                    Dmax = Math.min(pi * (1 - qj), (1 - pi) * qj);
                }  //check these equations
                D_prime_ij = Dij / Dmax;
                D_prime += (pi * qj * Math.abs(D_prime_ij));
            }
        }
        return D_prime;
    }

    double calculateRSqrDisequilibrium(int[][] contig, int rows, int cols) {
        //this is the D' approach used by Farnir 2000 Genome Research 10:220-227
        //this follows Hedrick 1987
        //but I have really made this up myself fusing Garnir with Hill
        double r_sqrsum = 0, pi, qj;
        int countAB, countAb, countaB, countab, nonmissingSampleSize = 0;
        int[] p_margin = new int[rows];
        int[] q_margin = new int[cols];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                nonmissingSampleSize += contig[i][j];
                p_margin[i] += contig[i][j];
                q_margin[j] += contig[i][j];
            }
        }
        if (nonmissingSampleSize <= minTaxaForEstimate) {
            return Double.NaN;
        }
        for (int i = 0; i < rows; i++) {//Through missing data & incomplete datasets some alleles can be fixed this return missing value
            //       if(p_margin[i]==0) continue;   //changed by ed 8-13-2004
            if (p_margin[i] == 0) {
                return Double.NaN;
            }
            pi = (double) p_margin[i] / (double) nonmissingSampleSize;
            for (int j = 0; j < cols; j++) {//Through missing data & incomplete datasets some alleles can be fixed this return missing value
                // if(q_margin[j]==0) continue;   //changed by ed 8-13-2004
                if (q_margin[j] == 0) {
                    return Double.NaN;
                }
                qj = (double) q_margin[j] / (double) nonmissingSampleSize;
                countAB = contig[i][j];
                countAb = p_margin[i] - countAB;
                countaB = q_margin[j] - countAB;
                countab = nonmissingSampleSize - countAB - countAb - countaB;
                r_sqrsum += pi * qj * calculateRSqrDisequilibrium(countAB, countAb, countaB, countab);
            }
        }
        return r_sqrsum;
    }

    /** Returns P-value estimate for a given pair of sites.  If there were only 2 alleles
     *  at each locus, then the Fisher Exact P-value (one-tail) is returned.  If more states
     *  then the permutaed Monte Carlo test is used.
     *  @param r is site 1
     *  @param c is site 2
     *  @return P-value
     */
    public double getP(int r, int c) {
        if (r <= c) {
            return pDiseq[r][c];
        } else {
            return pDiseq[c][r];
        }
    }

    /** Get number of gametes included in LD calculations (after missing data was excluded)
     *  @param r is site 1
     *  @param c is site 2
     *  @return number of gametes
     *
     */
    private int getN(int r, int c) {
        if (c <= r) {
            return (int) pDiseq[r][c];
        } else {
            return (int) pDiseq[c][r];
        }
    }

    /** Returns D' estimate for a given pair of sites
     *  @param r is site 1
     *  @param c is site 2
     *  @return D'
     */
    public double getDPrime(int r, int c) {
        if (r <= c) {
            return diseq[r][c];
        } else {
            return diseq[c][r];
        }
    }

    /** Returns r^2 estimate for a given pair of sites
     *  @param r is site 1
     *  @param c is site 2
     *  @return D'
     */
    public double getRSqr(int r, int c) {
        if (c <= r) {
            return diseq[r][c];
        } else {
            return diseq[c][r];
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
        StringWriter sw = new StringWriter();
        print(this, new PrintWriter(sw));
        return sw.toString();
    }

    /** print the LD to the PrintWrite */
    public void print(LinkageDisequilibriumMultiAllele ld, PrintWriter out) {

        out.print("LocusName1\t Chromosome1\t ChromoPosition1\t ");
        out.print("Site1\t NumberOfStates1\t States1\t Frequency1\t ");
        out.print("LocusName2\t Chromosome2\t ChromoPosition2\t ");
        out.print("Site2\t NumberOfStates2\t States2\t Frequency2\t ");
        out.print("R^2\t DPrime\t pDiseq\t N");
        out.println();
        for (int r = 0; r < theAlignment.getSiteCount(); r++) {
            String rState = getStatesForPrint(r);
            for (int c = 0; c <= r; c++) {
                String cState = getStatesForPrint(c);
                out.print(theAlignment.getLocusName(r) + "\t" + "\t" + "\t");
                out.print(theAlignment.getPositionInLocus(r) + "\t" + stateVector[r].size() + "\t" + rState + "\t" + "NotDone" + "\t");
                out.print(theAlignment.getLocusName(c) + "\t" + "\t" + "\t");
                out.print(theAlignment.getPositionInLocus(c) + "\t" + stateVector[c].size() + "\t" + cState + "\t" + "NotDone" + "\t");
                out.print(ld.getRSqr(r, c) + "\t" + ld.getDPrime(r, c) + "\t" + ld.getP(r, c) + "\t" + ld.getN(r, c));
                out.println();
            }
        }
        out.println();
        out.println("Taxa Included:");
        for (int i = 0; i < theAlignment.getSequenceCount(); i++) {
            out.println(theAlignment.getIdGroup().getIdentifier(i).getName());
        }

    }


//This could be removed if the Datatype class had a toString method that chooses between outputting states and characters
    private String getStatesForPrint(int site) {
        String s = "[";
        Character C;
        DataType dt;
        dt = theAlignment.getDataType();
        for (int i = 0; i < stateVector[site].size(); i++) {
            C = (Character) stateVector[site].get(i);
            s += dt.getFormattedString(C.charValue()) + " ";
        }
        s += "]";
        return s;
    }

    //Implementation of TableReport Interface
    /**
     * @return column names for the table
     */
    public Object[] getTableColumnNames() {
        String[] basicLabels = {"Site1", "NumberOfStates1", "States1", "Frequency1",
            "Site2", "NumberOfStates2", "States2", "Frequency2", "R^2", "DPrime", "pDiseq", "N"};
        String[] annotatedLabels = {"LocusName1", "Chromosome1", "ChromoPosition1", "Site1",
            "NumberOfStates1", "States1", "Frequency1", "LocusName2", "Chromosome2", "ChromoPosition2",
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
        data = new String[19];

        int x = (int) ((-1.0 + Math.sqrt(8.0 * row + 1)) / 2.0);
        int r = x + 1;
        int c = (int) (row - 0.5 * x * x - 0.5 * x);

        String rState = getStatesForPrint(r);
        String rStr = String.valueOf(r);

        String cState = getStatesForPrint(c);
        String cStr = String.valueOf(c);

        data[labelOffset++] = theAlignment.getLocusName(r);
        //Need to fix this...
        //data[labelOffset++] = String.valueOf(theAlignment.getChromosome(r));
        data[labelOffset++] = "chromosome";
        data[labelOffset++] = String.valueOf(theAlignment.getPositionInLocus(r));
        data[labelOffset++] = rStr;

        data[labelOffset++] = String.valueOf(stateVector[r].size());
        data[labelOffset++] = rState;
        data[labelOffset++] = NotImplemented;
        data[labelOffset++] = theAlignment.getLocusName(c);
        //Need to fix this...
        //data[labelOffset++] = String.valueOf(theAlignment.getChromosome(c));
        data[labelOffset++] = "chromosome";
        data[labelOffset++] = String.valueOf(theAlignment.getPositionInLocus(c));
        data[labelOffset++] = cStr;

        data[labelOffset++] = String.valueOf(stateVector[c].size());
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

    public void writeReport(String delimit, File saveFile) {

        if (saveFile == null) {
            return;
        }

        FileWriter fw = null;
        BufferedWriter bw = null;

        try {

            fw = new FileWriter(saveFile);
            bw = new BufferedWriter(fw);
            String NotImplemented = "NotImplemented";
            String NA = "N/A";
            java.text.NumberFormat nf = new java.text.DecimalFormat();
            nf.setMaximumFractionDigits(8);

            Object[] colNames = getTableColumnNames();
            int cols = colNames.length;
            for (int j = 0; j < cols; j++) {
                bw.write(colNames[j].toString());
                if (j < (cols - 1)) {
                    bw.write(delimit);
                }
            }
            bw.write("\n");

            for (int r = 0; r < theAlignment.getSiteCount(); r++) {
                String rState = getStatesForPrint(r);
                String rStr = String.valueOf(r);
                for (int c = 0; c < r; c++) {
                    String cState = getStatesForPrint(c);
                    String cStr = String.valueOf(c);
                    bw.write(theAlignment.getLocusName(r));
                    bw.write(delimit);
                    //bw.write(String.valueOf(theAlignment.getChromosome(r)));
                    bw.write("chromosome");
                    bw.write(delimit);
                    bw.write(String.valueOf(theAlignment.getPositionInLocus(r)));
                    bw.write(delimit);
                    bw.write(rStr);
                    bw.write(delimit);
                    bw.write(String.valueOf(stateVector[r].size()));
                    bw.write(delimit);
                    bw.write(rState);
                    bw.write(delimit);
                    bw.write(NotImplemented);
                    bw.write(delimit);
                    bw.write(theAlignment.getLocusName(c));
                    bw.write(delimit);
                    //bw.write(String.valueOf(theAlignment.getChromosome(c)));
                    bw.write("chromosome");
                    bw.write(delimit);
                    bw.write(String.valueOf(theAlignment.getPositionInLocus(c)));
                    bw.write(delimit);
                    bw.write(cStr);
                    bw.write(delimit);
                    bw.write(String.valueOf(stateVector[c].size()));
                    bw.write(delimit);
                    bw.write(cState);
                    bw.write(delimit);
                    bw.write(NotImplemented);
                    if (theAlignment.getLocusName(r).equals(theAlignment.getLocusName(c))) {
                        bw.write(delimit);
                        bw.write(String.valueOf(Math.abs(theAlignment.getPositionInLocus(r) - theAlignment.getPositionInLocus(c))));
                    } else {
                        bw.write(delimit);
                        bw.write(NA);
                    }
                    bw.write(delimit);
                    bw.write(nf.format(getRSqr(r, c)));
                    bw.write(delimit);
                    bw.write(nf.format(getDPrime(r, c)));
                    bw.write(delimit);
                    bw.write(nf.format(getP(r, c)));
                    bw.write(delimit);
                    bw.write(String.valueOf(getN(r, c)));
                    bw.write("\n");
                }

            }

        } catch (Exception e) {
            System.out.println("TableReportUtils: writeReport: problem writing file: " + e.getMessage());
        } finally {
            try {
                bw.close();
                fw.close();
            } catch (Exception e) {
                // do nothing
            }
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
        int total = theAlignment.getSiteCount() * (theAlignment.getSiteCount() - 1) / 2;
        return total;
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

}

