package net.maizegenetics.tassel;


import net.maizegenetics.baseplugins.*;
import net.maizegenetics.baseplugins.gdpc.ConvertGDPCToPALPlugin;
import net.maizegenetics.baseplugins.gdpc.GDPCImportPlugin;
import net.maizegenetics.baseplugins.gdpc.GDPCPlugin;
import net.maizegenetics.baseplugins.numericaltransform.NumericalTransformPlugin;

import java.awt.*;



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
    
    private GDPCPlugin myGDPCPlugin;
    
    private GDPCImportPlugin myGDPCImportPlugin;

    private PlinkLoadPlugin myPlinkLoadPlugin;

    private FlapjackLoadPlugin myFlapjackLoadPlugin;
    
    
    public DataControlPanel(TASSELMainFrame theQAF, DataTreePanel theDTP) {
        super(theQAF, theDTP);
        
        try {
            jbInit();
            
            myGDPCPlugin = new GDPCPlugin(theTASSELMainFrame, true);
            
            myGDPCImportPlugin = new GDPCImportPlugin();
            
            ConvertGDPCToPALPlugin convertPlugin = new ConvertGDPCToPALPlugin();
            myGDPCPlugin.addGDPCListener(convertPlugin);
            myGDPCImportPlugin.addListener(convertPlugin);
            convertPlugin.addListener(theDataTreePanel);
            addMainDisplayPlugin(myGDPCPlugin);

            myPlinkLoadPlugin = new PlinkLoadPlugin(theTASSELMainFrame, true);
            myPlinkLoadPlugin.addListener(theDataTreePanel);

            myFlapjackLoadPlugin = new FlapjackLoadPlugin(theTASSELMainFrame, true);
            myFlapjackLoadPlugin.addListener(theDataTreePanel);

            addPlugin(new FileLoadPlugin(theTASSELMainFrame, true, myGDPCImportPlugin, myPlinkLoadPlugin, myFlapjackLoadPlugin));
            addPlugin(new ExportPlugin(theTASSELMainFrame, true));
            addPlugin(new FilterAlignmentPlugin(theTASSELMainFrame, true));
            addPlugin(new FilterTaxaAlignmentPlugin(theTASSELMainFrame, true));
            addPlugin(new FilterTraitsPlugin(theTASSELMainFrame, true));
            //todo the kinship filtering needs to be reworked
            //kinship need to changed to matrices, and MLM should just work without filtering
            //addPlugin(new GenotypeTransformPlugin(theTASSELMainFrame, true));
            addPlugin(new GenotypeImputationPlugin(theTASSELMainFrame, true));
            addPlugin(new NumericalTransformPlugin(theTASSELMainFrame, true));
            addPlugin(new SynonymizerPlugin(theTASSELMainFrame, true));
            addPlugin(new UnionAlignmentPlugin(theTASSELMainFrame, true));
            addPlugin(new IntersectionAlignmentPlugin(theTASSELMainFrame, true));
            addPlugin(new SeparatePlugin(theTASSELMainFrame, true));
            
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }
    
    
    void jbInit() throws Exception {
        
        this.setLayout(new FlowLayout());
        this.setFont(new java.awt.Font("Dialog", 1, 12));
        this.setMaximumSize(new Dimension(32767, 200));
        this.setMinimumSize(new Dimension(700, 34));
        
    }
    
    public void setDataTreePanel(DataTreePanel theDataTreePanel) {
        this.theDataTreePanel = theDataTreePanel;
    }
    
    
    public boolean saveGDPCSettings() {
        myGDPCPlugin.saveSettings();
        return true;
    }
}