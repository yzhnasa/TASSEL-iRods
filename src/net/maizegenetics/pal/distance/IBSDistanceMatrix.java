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
import net.maizegenetics.pal.alignment.TBitAlignment;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.BitUtil;
import net.maizegenetics.util.OpenBitSet;
import net.maizegenetics.util.ProgressListener;

/**
 * This class calculates an identity by state matrix.  It is scaled so only non-missing data is used.
 * Class needs to be updated to use TBit alignment, and then bit calculations
 *
 * @author Ed Buckler
 * @version 1.0
 */
public class IBSDistanceMatrix extends DistanceMatrix {

    private ProgressListener myListener = null;
    private int numSeqs;
    private Alignment theAlignment;
    private TBitAlignment theTBA=null;
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
        if(theAlignment instanceof TBitAlignment) {
            theTBA=(TBitAlignment)theAlignment;
        } else {
            theTBA=TBitAlignment.getInstance(theAlignment,2,false);
        }
    //  this should have an option to only use the 2 or 3 most common alleles
        setIdGroup(theAlignment.getIdGroup());
        long time=System.currentTimeMillis();
//        computeDistances();       
//        System.out.println("Old Distance Time:"+(System.currentTimeMillis()-time));
//        time=System.currentTimeMillis();
//        computeBitDistances();
//        System.out.println("NewBit Distance Time:"+(System.currentTimeMillis()-time));
        time=System.currentTimeMillis();
        computeHetBitDistances();
        System.out.println("NewBitHet Distance Time:"+(System.currentTimeMillis()-time));
    }
    
    /**
     * Only work for inbreds
     */
    private void computeBitDistances() {
        avgTotalSites = 0;
        int count = 0;
        double[] params;
        double[][] distance = new double[numSeqs][numSeqs];
        for (int i = 0; i < numSeqs; i++) {
            distance[i][i] = 0;
            BitSet iMj=theTBA.getAllelePresenceForAllSites(i, 0);
            BitSet iMn=theTBA.getAllelePresenceForAllSites(i, 1);
            //there are lots of approaches to deal with the hets
            for (int j = i + 1; j < numSeqs; j++) {
                BitSet jMj=theTBA.getAllelePresenceForAllSites(j, 0);
                BitSet jMn=theTBA.getAllelePresenceForAllSites(j, 1);
                long same=OpenBitSet.intersectionCount(iMj, jMj)+OpenBitSet.intersectionCount(iMn, jMn);
                long diff=OpenBitSet.intersectionCount(iMj, jMn)+OpenBitSet.intersectionCount(iMn, jMj);
                double identity=(double)same/(double)(same+diff);
                double dist=1-identity;
                distance[i][j] = distance[j][i] = dist;
                avgTotalSites += (same+diff);  //this assumes not hets
                count++;
            }
            fireProgress((int) (((double) (i + 1) / (double) numSeqs) * 100.0));
        }
        lastCachedRow = null;
        setDistances(distance);
        avgTotalSites /= (double) count;
    }
    
    
    /**
     * This is a cleanest, fastest and most accurate way to calculate distance.
     * It includes a simple approach for hets.
     * This method is actually 10% faster than the inbred approach.
     * Key reason for the speed increase is only 3 bit counts versus the 4 above.
     */
    private void computeHetBitDistances() {
        avgTotalSites = 0;
        int count = 0;
        double[] params;
        double[][] distance = new double[numSeqs][numSeqs];
        for (int i = 0; i < numSeqs; i++) {
            distance[i][i] = 0;
            long[] iMj=theTBA.getAllelePresenceForAllSites(i, 0).getBits();
            long[] iMn=theTBA.getAllelePresenceForAllSites(i, 1).getBits();
            for (int j = i + 1; j < numSeqs; j++) {
                long[] jMj=theTBA.getAllelePresenceForAllSites(j, 0).getBits();
                long[] jMn=theTBA.getAllelePresenceForAllSites(j, 1).getBits();
                int sameCnt=0, diffCnt=0, hetCnt=0;
                for(int x=0; x<iMj.length; x++) {
                    long same=(iMj[x]&jMj[x])|(iMn[x]&jMn[x]);
                    long diff=(iMj[x]&jMn[x])|(iMn[x]&jMj[x]);
                    long hets=same&diff;
                    sameCnt+=BitUtil.pop(same);
                    diffCnt+=BitUtil.pop(diff);
                    hetCnt+=BitUtil.pop(hets);
                }
                double identity=(double)(sameCnt+(hetCnt/2))/(double)(sameCnt+diffCnt+hetCnt);
                double dist=1-identity;
                distance[i][j] = distance[j][i] = dist;
                avgTotalSites += sameCnt+diffCnt+hetCnt;  //this assumes not hets
                count++;
            }
            fireProgress((int) (((double) (i + 1) / (double) numSeqs) * 100.0));
        }
        lastCachedRow = null;
        setDistances(distance);
        avgTotalSites /= (double) count;
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
                count++;
            }
            fireProgress((int) (((double) (i + 1) / (double) numSeqs) * 100.0));
        }
        lastCachedRow = null;
        setDistances(distance);
        avgTotalSites /= (double) count;
    }

    private double[] calculateDistance(int s1, int s2) {

        int siteCount = theAlignment.getSiteCount();

        if ((lastCachedRow == null) || (s1 != lastCachedRowNum)) {
            lastCachedRowNum = s1;
            lastCachedRow = theAlignment.getBaseRange(s1, 0, siteCount);
        }

        byte[] s2Row = theAlignment.getBaseRange(s2, 0, siteCount);

        double[] params = new double[2];
        int numIdentical = 0, numDifferent = 0;
        for (int i = 0; i < siteCount; i++) {
            byte s1b = lastCachedRow[i];
            byte s2b = s2Row[i];
            if ((s1b != Alignment.UNKNOWN_DIPLOID_ALLELE) && (s2b != Alignment.UNKNOWN_DIPLOID_ALLELE)) {
                if (s1b == s2b) {
                    numIdentical++;
                } else {
                    numDifferent++;
                }
            }
        }
        params[1] = numIdentical + numDifferent;
        if (params[1] == 0) {
            params[0] = Double.NaN;
        } else {
            params[0] = (double) numDifferent / params[1];
        }
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
