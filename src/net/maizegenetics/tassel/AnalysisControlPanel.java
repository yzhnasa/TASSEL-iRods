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
            ConvertSBitTBitPlugin convertAlignment = new ConvertSBitTBitPlugin(theQAF, theDTP);
            addPlugin(new SequenceDiversityPlugin(theTASSELMainFrame, true));
            addPlugin(new LinkageDisequilibriumPlugin(theTASSELMainFrame, true));
            addPlugin(new CreateTreePlugin(theTASSELMainFrame, true, convertAlignment));
            // addPlugin(new ExtractSNPAssaysPlugin(theTASSELMainFrame, true));
            // addPlugin(new LogisticRegressionAssocPlugin(theTASSELMainFrame,true));
            addPlugin(new KinshipPlugin(theTASSELMainFrame, true));
            addPlugin(new FixedEffectLMPlugin(theTASSELMainFrame, true));
            addPlugin(new MLMPlugin(theTASSELMainFrame, true));
            addPlugin(new RidgeRegressionEmmaPlugin(theTASSELMainFrame, true));
            addPlugin(new GenotypeSummaryPlugin(theTASSELMainFrame, true));
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }
}