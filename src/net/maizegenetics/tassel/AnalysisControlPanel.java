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
package net.maizegenetics.tassel;

import net.maizegenetics.baseplugins.*;
import net.maizegenetics.baseplugins.genomicselection.RidgeRegressionEmmaPlugin;

import java.awt.Dimension;
import java.awt.FlowLayout;

/**
 * Title:        TASSEL
 * Description:  A java program to deal with diversity
 * Copyright:    Copyright (c) 2000
 * Company:      USDA-ARS/NCSU
 *
 * @author Ed Buckler
 * @version 1.0
 */
public class AnalysisControlPanel extends AbstractControlPanel {

    public AnalysisControlPanel(TASSELMainFrame theQAF, DataTreePanel theDTP) {
        super(theQAF, theDTP);
        try {
            addPlugin(new SequenceDiversityPlugin(theTASSELMainFrame, true));
            addPlugin(new LinkageDisequilibriumPlugin(theTASSELMainFrame, true));
            addPlugin(new CreateTreePlugin(theTASSELMainFrame, true));
            // addPlugin(new ExtractSNPAssaysPlugin(theTASSELMainFrame, true));
            // addPlugin(new LogisticRegressionAssocPlugin(theTASSELMainFrame,true));
            addPlugin(new KinshipPlugin(theTASSELMainFrame, true));
            addPlugin(new FixedEffectLMPlugin(theTASSELMainFrame, true));
            addPlugin(new MLMPlugin(theTASSELMainFrame, true));
            addPlugin(new RidgeRegressionEmmaPlugin(theTASSELMainFrame,true));
            jbInit();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    private void jbInit() throws Exception {

        this.setLayout(new FlowLayout());
        this.setFont(new java.awt.Font("Dialog", 1, 12));
        this.setMaximumSize(new Dimension(32767, 200));
        this.setMinimumSize(new Dimension(745, 34));
        this.setPreferredSize(new Dimension(745, 34));

    }
}