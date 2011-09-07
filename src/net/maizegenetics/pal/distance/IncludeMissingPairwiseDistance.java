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
 * @author Ed Buckler
 * @version 1.0
 */
public class IncludeMissingPairwiseDistance extends IgnoreMissingPairwiseDistance {

    public IncludeMissingPairwiseDistance(SitePattern sp) {
        super(sp);
    }

    protected boolean isDifferent(int s1, int s2) {

        // Check for identity
        if (s1 == s2) {
            return false;
        }

        // The remaining pairs are all different
        return true;
    }

    protected double[] getObservedDistance(byte[] seqPat1, byte[] seqPat2) {//this method also handles when weights are set to zero
        double[] params = new double[2];
        int diff = 0;
        int totalSitesPresent = 0;
        for (int i = 0; i < numPatterns; i++) {
            if (weight[i] > 0) {
                totalSitesPresent += weight[i];
                // ? is considered identical to anything
                if (seqPat1[i] != numStates && seqPat2[i] != numStates) {
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