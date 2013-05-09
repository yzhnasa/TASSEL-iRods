package net.maizegenetics.tassel;

import net.maizegenetics.baseplugins.*;
import net.maizegenetics.baseplugins.numericaltransform.NumericalTransformPlugin;

/**
 * Title:        TASSEL
 * Description:  A java program to deal with diversity
 * Copyright:    Copyright (c) 2000
 * Company:      USDA-ARS/NCSU
 *
 * @author Ed Buckler
 * @version 1.0
 */
public class DataControlPanel extends AbstractControlPanel {

    private PlinkLoadPlugin myPlinkLoadPlugin;
    private FlapjackLoadPlugin myFlapjackLoadPlugin;

    public DataControlPanel(TASSELMainFrame theQAF, DataTreePanel theDTP) {
        super(theQAF, theDTP);

        try {

            myPlinkLoadPlugin = new PlinkLoadPlugin(theTASSELMainFrame, true);
            myPlinkLoadPlugin.addListener(theDataTreePanel);

            myFlapjackLoadPlugin = new FlapjackLoadPlugin(theTASSELMainFrame, true);
            myFlapjackLoadPlugin.addListener(theDataTreePanel);

            addPlugin(new FileLoadPlugin(theTASSELMainFrame, true, myPlinkLoadPlugin, myFlapjackLoadPlugin));
            addPlugin(new ExportPlugin(theTASSELMainFrame, true));
            addPlugin(new ConvertSBitTBitPlugin(theTASSELMainFrame, true));
            addPlugin(new FilterAlignmentPlugin(theTASSELMainFrame, true));
            addPlugin(new FilterSiteNamePlugin(theTASSELMainFrame, true));
            addPlugin(new FilterTaxaAlignmentPlugin(theTASSELMainFrame, true));
            addPlugin(new FilterTaxaPropertiesPlugin(theTASSELMainFrame, true));
            addPlugin(new FilterTraitsPlugin(theTASSELMainFrame, true));
            //todo the kinship filtering needs to be reworked
            //kinship need to changed to matrices, and MLM should just work without filtering
            //addPlugin(new GenotypeTransformPlugin(theTASSELMainFrame, true));
            addPlugin(new GenotypeImputationPlugin(theTASSELMainFrame, true));
            addPlugin(new NumericalTransformPlugin(theTASSELMainFrame, true));
            addPlugin(new SynonymizerPlugin(theTASSELMainFrame, true));
            addPlugin(new UnionAlignmentPlugin(theTASSELMainFrame, true));
            addPlugin(new IntersectionAlignmentPlugin(theTASSELMainFrame, true));
            addPlugin(new MergeAlignmentsPlugin(theTASSELMainFrame, true));
            addPlugin(new SeparatePlugin(theTASSELMainFrame, true));

        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }
}
