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

import net.maizegenetics.pal.alignment.SitePattern;

/**
 * Title:        TASSEL
 * Description:  A java program to deal with diversity
 * Copyright:    Copyright (c) 2000
 * Company:      USDA-ARS/NCSU
 *
 * @author Ed Buckler
 * @version 1.0
 */
public class GeneralizedAlignmentDistanceMatrix extends DistanceMatrix {

    private int numSeqs;
    private IgnoreMissingPairwiseDistance pwd;
    private SitePattern sitePattern;
    /**
     * PairwiseDistance mode:  missing is counted and goes in total (net.maizegenetics.pal DEFAULT)
     */
    public static final int INCLUDE_MISSING = 1;
    /**
     * PairwiseDistance mode:  missing is ignored
     */
    public static final int IGNORE_MISSING = 2;
    /**
     * PairwiseDistance mode:  synonmyous site of codon alignment are calculated, missing is ignored
     */
    public static final int SYNONYMOUS_IGNORE_MISSING = 3;
    /**
     * PairwiseDistance mode:  nonsynonmyous site of codon alignment are calculated, missing is ignored
     */
    public static final int NONSYNONYMOUS_IGNORE_MISSING = 4;
    int pairwiseDistMode;
    /**
     * Holds the average numbers of sites in the comparisons
     */
    double avgTotalSites;

    /**
     * compute observed distances
     *
     * @param sp               site pattern
     * @param pairwiseDistMode PairwiseDistance calculation mode (for example IGNORE_MISSING);
     */
    public GeneralizedAlignmentDistanceMatrix(SitePattern sp, int pairwiseDistMode) {
        super();
        numSeqs = sp.getSequenceCount();
        sitePattern = sp;
        setIdGroup(sp.getIdGroup());

        this.pairwiseDistMode = pairwiseDistMode;
        switch (pairwiseDistMode) {
            case INCLUDE_MISSING: {
                pwd = new IncludeMissingPairwiseDistance(sp);
                break;
            }
            case IGNORE_MISSING: {
                pwd = new IgnoreMissingPairwiseDistance(sp);
                break;
            }
        }

//		pwd = new PairwiseDistance(sp);

        computeDistances();
    }

    /**
     * recompute maximum-likelihood distances under new site pattern
     *
     * @param sp site pattern
     */
    public void recompute(SitePattern sp) {
        pwd.updateSitePattern(sp);

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
                params = pwd.getDistance(i, j);
                distance[i][j] = params[0];
                distance[j][i] = params[0];
                avgTotalSites += params[1];
//                    System.err.println("i="+i+" j="+j+"  params[0]="+params[0]+"  params[1]="+params[1]);
                count++;
            }
        }
        setDistances(distance);


        avgTotalSites /= (double) count;
//      System.err.println("avgTotalSites="+avgTotalSites);
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
}
