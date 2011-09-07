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


import net.maizegenetics.baseplugins.TableDisplayPlugin;
import net.maizegenetics.baseplugins.TreeDisplayPlugin;
import net.maizegenetics.baseplugins.Grid2dDisplayPlugin;
import net.maizegenetics.baseplugins.LinkageDiseqDisplayPlugin;
import net.maizegenetics.baseplugins.chart.ChartDisplayPlugin;
import net.maizegenetics.baseplugins.QQDisplayPlugin;

import java.awt.event.ActionEvent;
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
            addPlugin(new TableDisplayPlugin(theTASSELMainFrame,true));
            addPlugin(new TreeDisplayPlugin(theTASSELMainFrame,true));
            addPlugin(new Grid2dDisplayPlugin(theTASSELMainFrame, true));
            addPlugin(new LinkageDiseqDisplayPlugin(theTASSELMainFrame, true));
            addPlugin(new ChartDisplayPlugin(theTASSELMainFrame, true));
            addPlugin(new QQDisplayPlugin(theTASSELMainFrame, true));
            addPlugin(new ManhattanDisplayPlugin(theTASSELMainFrame, true));
            
            // If jbInit is needed later, we should implement the functions
            // as Plugins.  -Terry
            //jbInit();
            
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    /*
    void jbInit() {
        this.setLayout(new FlowLayout());
        this.setFont(new java.awt.Font("Dialog", 1, 12));
        this.setMaximumSize(new Dimension(32767, 200));
        this.setMinimumSize(new Dimension(645, 34));
        this.setPreferredSize(new Dimension(645, 34));

        URL imageURL;



        imageURL = ResultControlPanel.class.getResource("images/GenePlot.gif");
        ImageIcon genePlotIcon = null;
        if(imageURL != null)genePlotIcon = new ImageIcon(imageURL);
        if(genePlotIcon != null)geneButton.setIcon(genePlotIcon);
        geneButton.setText("Gene Plot");
        geneButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(ActionEvent e) {
                geneButton_actionPerformed(e);
            }
        });
        geneButton.setFont(new java.awt.Font("Dialog", 1, 12));
        geneButton.setToolTipText("Plot Gene Structure");
        geneButton.setMargin(new Insets(2, 2, 2, 2));
        geneButton.setBackground(Color.white);

        ImageIcon barChartIcon = null;

        VisualGenoType.setBackground(Color.white);
        VisualGenoType.setToolTipText("Visual Geno Type");
        VisualGenoType.setText("VGtype");
        VisualGenoType.setFont(new java.awt.Font("Dialog", 1, 12));
        imageURL = ResultControlPanel.class.getResource("images/VizGenotype_25x17.gif");
        ImageIcon vizGenIcon = null;
        if(imageURL != null)vizGenIcon = new ImageIcon(imageURL);
        if(barChartIcon != null)VisualGenoType.setIcon(vizGenIcon);
        VisualGenoType.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(ActionEvent e) {
                VisualGenoType_actionPerformed(e);
            }
        });

        hapViewButton.setBackground(Color.white);
        hapViewButton.setToolTipText("Haplotype Viewer");
        hapViewButton.setText("HapView");
        hapViewButton.setFont(new java.awt.Font("Dialog", 1, 12));
        imageURL = ResultControlPanel.class.getResource("images/HapViewer_23x19.gif");
        ImageIcon hapViewerIcon = null;
        if(imageURL != null)hapViewerIcon  = new ImageIcon(imageURL);
        if(barChartIcon != null)hapViewButton.setIcon(hapViewerIcon);
        hapViewButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(ActionEvent e) {
                hapViewButton_actionPerformed(e);
            }
        });

        this.add(geneButton, null);
        this.add(VisualGenoType, null);
        this.add(hapViewButton, null);
    }
     */



    void geneButton_actionPerformed(ActionEvent e) {
   /*     Object[] data = theDataTreePanel.getSelectedData();
        String[] names = theDataTreePanel.getSelectedNames();
        Vector trVector = new Vector();
        net.maizegenetics.pal.report.TableReport[] tr = null;
        net.maizegenetics.pal.alignment.AnnotationAlignment aa = null;
        if (data[0] instanceof net.maizegenetics.pal.alignment.AnnotationAlignment) {
            JOptionPane.showMessageDialog(theTASSELMainFrame, TASSELMainFrame.GENOTYPE_DATA_NEEDED);
        }else{
            for (int i = 0; i < data.length; i++) {
                if (data[i] instanceof net.maizegenetics.pal.report.TableReport) {
                    trVector.add((net.maizegenetics.pal.report.TableReport) data[i]);
                }
                if (data[i] instanceof net.maizegenetics.pal.alignment.AnnotationAlignment) {
                    aa = (net.maizegenetics.pal.alignment.AnnotationAlignment) data[i];
                }
                // todo currently there is no way to import annotation data; GDPC should eventually provide annotation data
                //
                //todo this module needs to be cleaned up/fixed
            }
            if (aa != null) {
                boolean annotatedSites = false;
                // check to see if this alignment has annotation information
                for (int i = 0; i < aa.getSiteCount(); i++) {
                    if (aa.getPositionType(i) != 0) {
                        annotatedSites = true;
                        continue;
                    }
                }
                if (annotatedSites == false) {
                    JOptionPane.showMessageDialog(theTASSELMainFrame, "Selected data does not have annotation information");

                } else {
                    if (trVector.size() > 0) {
                        tr = new net.maizegenetics.pal.report.TableReport[trVector.size()];
                        for (int i = 0; i < trVector.size(); i++) {
                            tr[i] = (net.maizegenetics.pal.report.TableReport) trVector.get(i);
                        }
                    }
                    GeneDialog theGeneDialog = new GeneDialog(theTASSELMainFrame, aa, tr);
                    theGeneDialog.setVisible(true);
                }
            }
        }
        */
    }

    void VisualGenoType_actionPerformed(ActionEvent e) {
/*        Object[] data = theDataTreePanel.getSelectedData();
        Vector trVector = new Vector();
        net.maizegenetics.pal.alignment.AnnotationAlignment aa = null;
        if(data.length == 0){
            JOptionPane.showMessageDialog(theTASSELMainFrame, TASSELMainFrame.GENOTYPE_DATA_NEEDED);
        }
        for (int i = 0; i < data.length; i++) {
            //if(data[i] instanceof net.maizegenetics.pal.report.TableReport)
            //{trVector.add((net.maizegenetics.pal.report.TableReport)data[i]);}
            if (data[i] instanceof net.maizegenetics.pal.alignment.AnnotationAlignment) {
                aa = (net.maizegenetics.pal.alignment.AnnotationAlignment) data[i];
                VisualGenotypePlot theVGtype = new VisualGenotypePlot(theTASSELMainFrame, aa);
                theVGtype.setVisible(true);
            }else{
                JOptionPane.showMessageDialog(theTASSELMainFrame, TASSELMainFrame.GENOTYPE_DATA_NEEDED);
            }
        }  */
    }

    void hapViewButton_actionPerformed(ActionEvent e) {
  /*      Object[] data = theDataTreePanel.getSelectedData();
        if (data[0] instanceof net.maizegenetics.pal.alignment.AnnotationAlignment) {
            HaplotypeFrame hf = new HaplotypeFrame(theTASSELMainFrame, theDataTreePanel, (net.maizegenetics.pal.alignment.AnnotationAlignment) data[0]);
        }
        else {
            JOptionPane.showMessageDialog(theTASSELMainFrame, TASSELMainFrame.GENOTYPE_DATA_NEEDED);
        }
      */
    }
}