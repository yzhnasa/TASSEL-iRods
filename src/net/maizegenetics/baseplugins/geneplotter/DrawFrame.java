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
package net.maizegenetics.baseplugins.geneplotter;

import net.maizegenetics.pal.alignment.Alignment;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyVetoException;
import java.beans.VetoableChangeListener;

/**
 * Title:        TASSEL
 * Description:  A java program to deal with diversity
 * Copyright:    Copyright (c) 2000
 * Company:      USDA-ARS/NCSU
 * @author Ed Buckler
 * @version 1.0
 */
public class DrawFrame extends JInternalFrame implements VetoableChangeListener {

    net.maizegenetics.pal.report.TableReport theTableSource;
    BorderLayout borderLayout1 = new BorderLayout();
    GenePainter genePanel;
    GenePlotter[] graphPanel = new GenePlotter[4];
    JPanel jPanel1 = new JPanel();
    JButton SettingsButton = new JButton();
    BorderLayout borderLayout2 = new BorderLayout();
    Alignment aa;
    JFrame mainFrame;
    boolean flag;
    String[] numericColumnHeaders;
    double[][] numericData;
    int windows, option = 0;
    GeneDialog thegeneDialog;
    private JDesktopPane desktopPane = new JDesktopPane();
    int[] optionArray;
    int referenceDataIndex;
    GenePaintSettings thegenepaint = null;
    GenePlotSettings thegeneplot = null;

    public DrawFrame(JFrame frame, Alignment aa, net.maizegenetics.pal.report.TableReport tr, boolean flag, int windows, GeneDialog geneDialog) {
        this.aa = aa;
        this.theTableSource = tr;
        this.mainFrame = frame;
        this.flag = flag;
        this.windows = windows;
        this.thegeneDialog = geneDialog;
        addVetoableChangeListener(this);

        if (flag) {
            genePanel = new GenePainter(aa);
            thegenepaint = new GenePaintSettings(mainFrame, aa, genePanel);
        } else if (theTableSource != null) {
            findNumericData(theTableSource.getTableColumnNames(), theTableSource.getTableData());
            sortNumericData(referenceDataIndex);
            graphPanel[windows - 1] = new GenePlotter(aa, numericData, numericColumnHeaders, referenceDataIndex);
            thegeneplot = new GenePlotSettings(mainFrame, aa, graphPanel[windows - 1], numericColumnHeaders, numericData);
        }
        try {
            jbInit();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private void jbInit() throws Exception {
        jPanel1.setLayout(borderLayout2);
        this.setIconifiable(true);
        this.setMaximizable(true);
        this.getContentPane().setLayout(borderLayout1);
        SettingsButton.setText("Settings");
        SettingsButton.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                SettingsButton_actionPerformed(e);
            }
        });
        if (flag) {
            this.setTitle("Gene Paint");
            genePanel.setBorder(BorderFactory.createEtchedBorder());
            this.getContentPane().add(jPanel1, null);
            jPanel1.add(genePanel, BorderLayout.CENTER);
            jPanel1.add(SettingsButton, BorderLayout.EAST);
            this.show();
        } else if (theTableSource != null) {
            this.setClosable(true);
            this.setTitle("Gene Plot : " + getGenePlotTitle());
            graphPanel[windows - 1].setBorder(BorderFactory.createEtchedBorder());
            this.getContentPane().add(jPanel1, null);
            jPanel1.add(graphPanel[windows - 1], BorderLayout.CENTER);
            jPanel1.add(SettingsButton, BorderLayout.EAST);
            this.show();
        }
    }

    void SettingsButton_actionPerformed(ActionEvent e) {
        if (flag) {
//          GenePaintSettings thegenepaint = new GenePaintSettings(theTASSELMainFrameenePanel);
            thegenepaint.setVisible(true);
            repaint();
        } else if (theTableSource != null) {
//           GenePlotSettings thegeneplot = new GenePlotSettings(theTASSELMainFrameraphPanel[windows-1], numericColumnHeaders, numericData);
            thegeneplot.setVisible(true);
            repaint();
        }
    }

    void CloseButton_actionPerformed(ActionEvent e) {
        this.setVisible(false);
    }

    private boolean findNumericData(Object[] header, Object[][] data) {
        //This code is inherited from Jack Liu PowerSSR, and it puts the cols in the row.  Essentially it transposes the matrix

        if (header.length <= 1) {
            return false;
        }
        if (data.length < 1) {
            return false;
        }
        for (int i = 0; i < data.length; i++) {
            if (data[i].length != header.length) {
                return false;
            }
        }
        String[] tempColumnHeader = new String[header.length];
        double[][] tempColumnData = new double[header.length][data.length];
        int[] badDataCount = new int[header.length];
        try {
            for (int i = 0; i < data.length; i++) {
                for (int j = 0; j < header.length; j++) {
                    tempColumnHeader[j] = header[j].toString();
                    try {
                        double d = Double.valueOf(data[i][j].toString()).doubleValue();
                        tempColumnData[j][i] = d;
                    } catch (NumberFormatException ex) {
                        tempColumnData[j][i] = Double.NaN;
                        badDataCount[j]++;
                    }
                }
            }
        } catch (NumberFormatException ex) {
            return false;
        }
        int goodCol = 0;
        for (int j = 0; j < header.length; j++) {
            if (badDataCount[j] < (data.length * 0.9)) {
                goodCol++;
            }
        }
        numericColumnHeaders = new String[goodCol];
        numericData = new double[goodCol][data.length];
        goodCol = 0;
        for (int j = 0; j < header.length; j++) {
            if (badDataCount[j] < (data.length / 2)) {
                numericColumnHeaders[goodCol] = tempColumnHeader[j];
                if (numericColumnHeaders[goodCol].equals("Site")) {
                    referenceDataIndex = goodCol;
                }
                numericData[goodCol] = tempColumnData[j];
                goodCol++;
            }
        }
        return true;
    }

    void sortNumericData(int sortCol) {
        //this numeric data puts the
        boolean swap = true;
        double temp;
        while (swap == true) {
            swap = false;
            for (int i = 1; i < numericData[0].length; i++) {
                if (numericData[sortCol][i - 1] > numericData[sortCol][i]) {
                    swap = true;
                    for (int j = 0; j < numericData.length; j++) {
                        temp = numericData[j][i - 1];
                        numericData[j][i - 1] = numericData[j][i];
                        numericData[j][i] = temp;
                    }
                }
            }
        }
    }

    public Dimension getSize() {
        return new Dimension(200, 200);
    }

    public void vetoableChange(PropertyChangeEvent pce) throws PropertyVetoException {
        if (pce.getPropertyName().equals(IS_CLOSED_PROPERTY)) {
            boolean changed = ((Boolean) pce.getNewValue()).booleanValue();
            if (changed) {
                int option = JOptionPane.showOptionDialog(this, "Close " +
                        getTitle() + "?",
                        "Close Confirmation",
                        JOptionPane.YES_NO_OPTION,
                        JOptionPane.QUESTION_MESSAGE,
                        null, null, null);
                if (option != JOptionPane.YES_OPTION) {
                    throw new PropertyVetoException("Cancelled", null);
                } else {
                    thegeneDialog.WindowClosed();
                }
            }
        }
    }

    public String getGenePlotTitle() {
        return (JOptionPane.showInputDialog(" Set Title for new Window "));
    }
}