/*
 * TASSEL - Trait Analysis by a aSSociation Evolution & Linkage
 * Copyright (C) 2003 Ed Buckler
 *
 * This software evaluates linkage disequilibrium nucletide diversity and 
 * associations. For more information visit http://www.maizegenetics.net
 *
 * This software is distributed under GNU general public license and without
 * any warranty ot technical support.
 *
 * You can redistribute and/or modify it under the terms of GNU General 
 * public license. 
 *
 */
package net.maizegenetics.pal.distance;

import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.datatype.DataType;
import net.maizegenetics.util.ProgressListener;

/**
 * This class calculates an identity by state matrix.  It is scaled so only non-missing data is used.
 *
 * @author Ed Buckler
 * @version 1.0
 */
public class IBSDistanceMatrix extends DistanceMatrix {

    private ProgressListener myListener = null;
    private int numSeqs;
    private Alignment theAlignment;
    /**
     * Holds the average numbers of sites in the comparisons
     */
    private double avgTotalSites;
    /**
     * Caches one row at a time to improve performance
     */
    private int lastCachedRowNum = 0;
    private byte[] lastCachedRow = null;

    /**
     * compute observed distances.  Missing sites are ignored.  This is no weighting for ambigous bases.
     *
     * @param theAlignment Alignment used to computed proportion that
     */
    public IBSDistanceMatrix(Alignment theAlignment) {
        this(theAlignment, null);
    }

    public IBSDistanceMatrix(Alignment theAlignment, ProgressListener listener) {
        super();
        myListener = listener;
        numSeqs = theAlignment.getSequenceCount();
        this.theAlignment = theAlignment;
        setIdGroup(theAlignment.getIdGroup());
        computeDistances();
    }

    private void computeDistances() {
        avgTotalSites = 0;
        int count = 0;
        double[] params;
        double[][] distance = new double[numSeqs][numSeqs];
        for (int i = 0; i < numSeqs; i++) {
            distance[i][i] = 0;
            for (int j = i + 1; j < numSeqs; j++) {
                params = calculateDistance(i, j);
                distance[i][j] = params[0];
                distance[j][i] = params[0];
                avgTotalSites += params[1];
                // System.err.println("i="+i+" j="+j+"  params[0]="+params[0]+"  params[1]="+params[1]);
                count++;
            }
            fireProgress((int) (((double) (i + 1) / (double) numSeqs) * 100.0));
        }
        lastCachedRow = null;
        setDistances(distance);
        avgTotalSites /= (double) count;
        // System.err.println("avgTotalSites="+avgTotalSites);
    }

    private double[] calculateDistance(int s1, int s2) {

        int siteCount = theAlignment.getSiteCount();

        if ((lastCachedRow == null) || (s1 != lastCachedRowNum)) {
            lastCachedRowNum = s1;
            lastCachedRow = theAlignment.getBaseRange(s1, 0, siteCount - 1);
        }

        byte[] s2Row = theAlignment.getBaseRange(s2, 0, siteCount - 1);

        double[] params = new double[2];
        int numIdentical = 0, numDifferent = 0;
        for (int i = 0; i < siteCount; i++) {
            byte s1b = lastCachedRow[i];
            byte s2b = s2Row[i];
            if ((s1b != DataType.UNKNOWN_CHARACTER) && (s2b != DataType.UNKNOWN_CHARACTER)) {
                if (s1b == s2b) {
                    numIdentical++;
                } else {
                    numDifferent++;
                }
            }
            // if(i<12) System.out.println(i+"  s1="+s1+" b:"+s1b+"  s2="+s2+" b:"+s2b+"  id="+numIdentical+" nonid="+numDifferent);
        }
        params[1] = numIdentical + numDifferent;
        if (params[1] == 0) {
            params[0] = Double.NaN;
        } else {
            params[0] = (double) numDifferent / params[1];
        }
        // System.out.println("  s1="+s1+"  s2="+s2+"  id="+numIdentical+" nonid="+numDifferent);
        return params;
    }

    public double getAverageTotalSites() {
        return avgTotalSites;
    }

    public double[][] getDistance() {
        return super.getDistances();
    }

    public String toString(int d) {
        double[][] distance = this.getDistance();
        /*Return a string representation of this matrix with 'd'
        displayed digits*/
        String newln = System.getProperty("line.separator");
        String outPut = new String();
        String num = new String();
        int i, j;
        java.text.NumberFormat nf = new java.text.DecimalFormat();
        nf.setMaximumFractionDigits(5);
        for (i = 0; i < distance.length; i++) {
            for (j = 0; j < distance[i].length; j++) {

                //Numeric x = new Numeric(distance[i][j]);
                num = nf.format(d);
                //num = x.toString(d);
                //ed change that screws up formatting
                //num=""+this.element[i][j];
                outPut = outPut + num + (char) 9;
            }
            outPut = outPut + newln;
        }
        return outPut;
    }

    public String toString() {
        return this.toString(6);
    }

    protected void fireProgress(int percent) {

        if (myListener != null) {
            myListener.progress(percent, null);
        }

    }
}
