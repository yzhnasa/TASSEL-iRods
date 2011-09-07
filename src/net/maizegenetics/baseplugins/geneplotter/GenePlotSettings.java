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
import javax.swing.border.TitledBorder;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.FocusEvent;

/**
 * Title:        TASSEL
 * Description:  A java program to deal with diversity
 * Copyright:    Copyright (c) 2000
 * Company:      USDA-ARS/NCSU
 * @author Ed Buckler
 * @version 1.0
 */
public class GenePlotSettings extends JDialog {

    Alignment aa;
    Color linecolor;
    String scales[] = {"500", "1000", "1500", "2000", "2500"};
    String Icons[] = {"o", "*", "+", "."};
    BorderLayout borderLayout1 = new BorderLayout();
    GenePlotter geneplot;
    JButton SetLineColor = new JButton();
    JCheckBox isBaseNumbers = new JCheckBox();
    JLabel BaseIncLabel = new JLabel();
    JComboBox BaseInc = new JComboBox(scales);
    JButton OkButton = new JButton();
    JLabel selectNumDataLabel = new JLabel();
    JComboBox SelectNumData = new JComboBox();
    String[] numericColumnHeaders;
    String[] thresholdvalues;
    double[][] numericData;
    JCheckBox DrawLineCheckbox = new JCheckBox();
    JCheckBox ThresholdCheckbox = new JCheckBox();
    JCheckBox SelIconsCheckbox = new JCheckBox();
    JComboBox SelectIcon = new JComboBox(Icons);
    boolean isdrawline = true, isthreshold = false, isicon = false, isbasenumber = true;
    int numdata = 0, iconstyle = 1, baseinc = 0;
    float value;
    JTextField ThresholdTextBox = new JTextField();
    TitledBorder titledBorder1;
    GridBagLayout gridBagLayout1 = new GridBagLayout();

    public GenePlotSettings(JFrame frame, Alignment aa, GenePlotter geneplot, String[] numericColumnHeaders, double[][] numericData) {
        super(frame, "Gene Plot Settings", true);
        this.aa = aa;
        this.geneplot = geneplot;
        this.numericColumnHeaders = numericColumnHeaders;
        this.numericData = numericData;
        SelectNumData = new JComboBox(numericColumnHeaders);
        try {
            jbInit();
            pack();
        } catch (Exception ex) {
        }
    }

    void jbInit() throws Exception {
        titledBorder1 = new TitledBorder("");
        SetLineColor.setEnabled(true);
        JPanel settingspanel = new JPanel();
        settingspanel.setPreferredSize(new Dimension(500, 425));
        settingspanel.setToolTipText("");
        settingspanel.setLayout(gridBagLayout1);
        this.setResizable(true);
        geneplot.setBackground(SystemColor.menu);
        geneplot.setForeground(Color.lightGray);
        SetLineColor.setBorder(BorderFactory.createRaisedBevelBorder());
        SetLineColor.setText("Set LineColor");
        SetLineColor.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                SetLineColor_actionPerformed(e);
            }
        });
        isBaseNumbers.setSelected(true);
        isBaseNumbers.setText("      Want Base Numbers ?");
        isBaseNumbers.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                isBaseNumbers_actionPerformed(e);
            }
        });
        BaseIncLabel.setToolTipText("");
        BaseIncLabel.setText("   Base Increment");
        OkButton.setText("OK");
        OkButton.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                OkButton_actionPerformed(e);
            }
        });
        BaseInc.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                BaseInc_actionPerformed(e);
            }
        });
        selectNumDataLabel.setText("   Select Numeric Data");
        SelectNumData.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                SelectNumData_actionPerformed(e);
            }
        });
        DrawLineCheckbox.setSelected(true);
        DrawLineCheckbox.setText("    Draw Line");
        DrawLineCheckbox.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                DrawLineCheckbox_actionPerformed(e);
            }
        });
        ThresholdCheckbox.setText("   Select Threshold");
        ThresholdCheckbox.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                ThresholdCheckbox_actionPerformed(e);
            }
        });
        SelIconsCheckbox.setText("   Select Icons");
        SelIconsCheckbox.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                SelIconsCheckbox_actionPerformed(e);
            }
        });
        SelectIcon.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                SelectIcon_actionPerformed(e);
            }
        });

        ThresholdTextBox.setBorder(titledBorder1);
        ThresholdTextBox.addFocusListener(new java.awt.event.FocusAdapter() {

            public void focusLost(FocusEvent e) {
                ThresholdTextBox_focusLost(e);
            }
        });

        settingspanel.add(SelectNumData, new GridBagConstraints(2, 0, 2, 1, 1.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(12, 0, 15, 57), 123, 4));
        settingspanel.add(DrawLineCheckbox, new GridBagConstraints(0, 1, 2, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(14, 95, 8, 6), 75, 0));
        settingspanel.add(SetLineColor, new GridBagConstraints(2, 1, 2, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(8, 0, 0, 56), 85, 18));
        settingspanel.add(ThresholdCheckbox, new GridBagConstraints(0, 2, 2, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(21, 96, 7, 0), 52, 0));
        settingspanel.add(selectNumDataLabel, new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(2, 83, 0, 29), 25, 35));
        settingspanel.add(SelIconsCheckbox, new GridBagConstraints(0, 3, 1, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(25, 95, 0, 0), 55, 0));
        settingspanel.add(SelectIcon, new GridBagConstraints(2, 3, 2, 1, 1.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(27, 8, 0, 55), 119, 5));
        settingspanel.add(isBaseNumbers, new GridBagConstraints(0, 4, 1, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(35, 43, 0, 12), 37, 0));
        settingspanel.add(BaseIncLabel, new GridBagConstraints(1, 4, 2, 1, 0.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(23, 0, 0, 0), 26, 23));
        settingspanel.add(BaseInc, new GridBagConstraints(3, 4, 1, 1, 1.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(31, 30, 0, 21), 20, 3));
        settingspanel.add(ThresholdTextBox, new GridBagConstraints(2, 2, 2, 1, 1.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(24, 0, 0, 51), 145, 0));
        settingspanel.add(OkButton, new GridBagConstraints(0, 5, 1, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(7, 157, 9, 17), 20, 0));
        this.getContentPane().add(settingspanel, BorderLayout.CENTER);
    }

    void SetLineColor_actionPerformed(ActionEvent e) {
        linecolor = JColorChooser.showDialog(null, "Choose Exon Color", linecolor);
        if (linecolor == null) {
            linecolor = Color.blue;
        }
    }

    void isBaseNumbers_actionPerformed(ActionEvent e) {
        isbasenumber = isBaseNumbers.isSelected();
    }

    void BaseInc_actionPerformed(ActionEvent e) {
        baseinc = BaseInc.getSelectedIndex();
    }

    void SelectNumData_actionPerformed(ActionEvent e) {
        numdata = SelectNumData.getSelectedIndex();
    }

    void SelectIcon_actionPerformed(ActionEvent e) {
        iconstyle = SelectIcon.getSelectedIndex();
    }

    void ThresholdCheckbox_actionPerformed(ActionEvent e) {
        isthreshold = ThresholdCheckbox.isSelected();
    }

    void SelIconsCheckbox_actionPerformed(ActionEvent e) {
        isicon = SelIconsCheckbox.isSelected();
    }

    void DrawLineCheckbox_actionPerformed(ActionEvent e) {
        isdrawline = DrawLineCheckbox.isSelected();
    }

    void OkButton_actionPerformed(ActionEvent e) {
        geneplot.isBaseNumbers(isbasenumber);
        geneplot.setBaseNumberInc(baseinc);
        geneplot.SetNumData(numdata);
        geneplot.setIcon(isicon);
        geneplot.setIconStyle(iconstyle);
        geneplot.setThreshold(isthreshold);
        geneplot.isDrawLine(isdrawline);
        geneplot.setLineColor(linecolor);
        if (isthreshold) {
            geneplot.SetThresholdData(value);
        }
        setVisible(false);
    }

    void ThresholdTextBox_focusLost(FocusEvent e) {
        try {
            if (isthreshold) {
                value = Float.parseFloat(ThresholdTextBox.getText());
            }
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }
}