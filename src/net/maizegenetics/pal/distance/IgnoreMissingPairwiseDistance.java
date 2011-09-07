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
import net.maizegenetics.pal.math.UnivariateMinimum;
import net.maizegenetics.pal.substmodel.SubstitutionModel;
import net.maizegenetics.pal.util.BranchLimits;

import java.io.Serializable;

/**
 * Title:        TASSEL
 * Description:  A java program to deal with diversity
 * Copyright:    Copyright (c) 2000
 * Company:      USDA-ARS/NCSU
 * @author Ed Buckler
 * @version 1.0
 */
public class IgnoreMissingPairwiseDistance implements Serializable {

    /** last estimated distance */
    public double distance;
    /** last estimate standard error of a distance */
    public double distanceSE;

    public IgnoreMissingPairwiseDistance(SitePattern sp) {
        updateSitePattern(sp);
    }

    /**
     * update model of substitution
     *
     * @param m model of substitution
     */
    public void updateModel(SubstitutionModel m) {
        of.updateModel(m);
    }

    /**
     * update site pattern
     *
     * @param sp site pattern
     */
    public void updateSitePattern(SitePattern sp) {
        sitePattern = sp;
        numSites = sp.getSiteCount();
        numPatterns = sp.getNumPatterns();
        numStates = sp.getDataType().getNumStates();
        weight = sp.getWeight();

        jcratio = ((double) numStates - 1.0) / (double) numStates;

        if (modelBased) {
            of.updateSitePattern(sp);
        }
    }

    /**
     * compute distance between two sequences in the given alignment
     *
     * @param s1 number of first sequence
     * @param s2 number of second sequence
     *
     * @return estimated distance (observed or ML, depending on constructor used)
     */
    public double[] getDistance(int s1, int s2) {
        return getDistance(sitePattern.getPattern()[s1], sitePattern.getPattern()[s2]);
    }

    /**
     * compute distance between two sequences (not necessarly
     * in the given alignment but with the same weights in the site pattern)
     *
     * @param s1 site pattern of first sequence
     * @param s2 site pattern of second sequence
     *
     * @return estimated distance (observed or ML, depending on constructor used)
     */
    public double[] getDistance(byte[] s1, byte[] s2) {
        double dist;
        double[] params;
        params = getObservedDistance(s1, s2);
        dist = params[0];

        if (modelBased && dist != 0.0) {
            // Apply generalized JC correction if possible
            double start = 1.0 - dist / jcratio;
            if (start > 0.0) {
                start = -jcratio * Math.log(start);
            } else {
                start = dist;
            }

            // Determine ML distance
            of.setSequences(s1, s2);
            if (start > BranchLimits.MAXARC || start < BranchLimits.MINARC) {
                // Don't use start value
                dist = um.findMinimum(of, BranchLimits.FRACDIGITS);
            } else {
                // Use start value
                dist = um.findMinimum(start, of, BranchLimits.FRACDIGITS);
            }
        }

        if (modelBased) {
            double f2x = um.f2minx;

            if (1.0 / (BranchLimits.MAXARC * BranchLimits.MAXARC) < f2x) {
                distanceSE = Math.sqrt(1.0 / f2x);
            } else {
                distanceSE = BranchLimits.MAXARC;
            }
        } else {
            distanceSE = 0.0;
        }

        params[0] = distance = dist;

        return params;
    }
    //
    // Private stuff
    //
    protected int numSites;
    protected int numPatterns;
    protected int numStates;
    protected int[] weight;
    protected double jcratio;
    protected boolean modelBased;
    protected SitePattern sitePattern;
    protected UnivariateMinimum um;
    protected SequencePairLikelihood of;

    protected boolean isDifferent(int s1, int s2) {

        // Check for identity
        if (s1 == s2) {
            return false;
        }

        // The remaining pairs are all different
        return true;
    }

    protected double[] getObservedDistance(byte[] seqPat1, byte[] seqPat2) {
        //this method also handles when weights are set to zero
        double[] params = new double[2];
        int diff = 0;
        int totalSitesPresent = 0;
        for (int i = 0; i < numPatterns; i++) {
            if (weight[i] > 0) {
                // ? is considered identical to anything
                if (seqPat1[i] != numStates && seqPat2[i] != numStates) {
                    totalSitesPresent += weight[i];
                    if (isDifferent(seqPat1[i], seqPat2[i])) {
                        diff += weight[i];
                    }
                }
            }
        }
        if (totalSitesPresent < 1) {
            totalSitesPresent = 1;
        }  //this could be removed if the downstream stuff can handle NaN
        params[0] = (double) diff / (double) totalSitesPresent;
        params[1] = (double) totalSitesPresent;
        return params;
    }
}