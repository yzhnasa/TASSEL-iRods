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

import net.maizegenetics.baseplugins.ArchaeopteryxPlugin;
import net.maizegenetics.baseplugins.TableDisplayPlugin;
import net.maizegenetics.baseplugins.TreeDisplayPlugin;
import net.maizegenetics.baseplugins.Grid2dDisplayPlugin;
import net.maizegenetics.baseplugins.LinkageDiseqDisplayPlugin;
import net.maizegenetics.baseplugins.chart.ChartDisplayPlugin;
import net.maizegenetics.baseplugins.QQDisplayPlugin;
import net.maizegenetics.baseplugins.ManhattanDisplayPlugin;

/**
 * Title:        TASSEL
 * Description:  A java program to deal with diversity
 * Copyright:    Copyright (c) 2000
 * Company:      USDA-ARS/NCSU
 *
 * @author Ed Buckler
 * @version 1.0
 */
public class ResultControlPanel extends AbstractControlPanel {

    public ResultControlPanel(TASSELMainFrame theQAF, DataTreePanel theDTP) {
        super(theQAF, theDTP);

        try {

            addPlugin(new TableDisplayPlugin(theTASSELMainFrame, true));
            addPlugin(new TreeDisplayPlugin(theTASSELMainFrame, true));
            addPlugin(new Grid2dDisplayPlugin(theTASSELMainFrame, true));
            addPlugin(new LinkageDiseqDisplayPlugin(theTASSELMainFrame, true));
            addPlugin(new ChartDisplayPlugin(theTASSELMainFrame, true));
            addPlugin(new QQDisplayPlugin(theTASSELMainFrame, true));
            addPlugin(new ManhattanDisplayPlugin(theTASSELMainFrame, true));

        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }
}
