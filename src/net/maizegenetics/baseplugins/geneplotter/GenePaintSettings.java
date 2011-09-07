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

/**
 * Title:        TASSEL
 * Description:  A java program to deal with diversity
 * Copyright:    Copyright (c) 2000
 * Company:      USDA-ARS/NCSU
 * @author Ed Buckler
 * @version 1.0
 */
public class GenePaintSettings extends JDialog {

    Alignment aa;
    private Color exoncolor;
    private Color linecolor;
    String scales[] = {"500", "1000", "1500", "2000", "5000"};
    String widths[] = {"1", "2", "3", "4", "5"};
    JComboBox Scaleinc = new JComboBox(scales);
    BorderLayout borderLayout1 = new BorderLayout();
    JButton SetExonColor = new JButton();
    JLabel ScaleInc = new JLabel();
    JLabel ExonWidth = new JLabel();
    JCheckBox CheckBox = new JCheckBox();
    JComboBox exonwidth = new JComboBox(widths);
    JComboBox scaleinc = new JComboBox(scales);
    JButton OkButton = new JButton();
    GenePainter genepaint;
    JButton SetLineColor = new JButton();
    JLabel LineWidth = new JLabel();
    JComboBox linewidth = new JComboBox(widths);
    GridBagLayout gridBagLayout1 = new GridBagLayout();
    JFrame mainFrame = null;

    public GenePaintSettings(JFrame frame, Alignment aa, GenePainter genepaint) {
        super(frame, "Gene Paint Settings", true);
        this.mainFrame = frame;
        this.aa = aa;
        this.genepaint = genepaint;
        try {
            jbInit();
            pack();
        } catch (Exception ex) {
        }
    }

    void jbInit() throws Exception {
        JPanel settingspanel = new JPanel();
        settingspanel.setPreferredSize(new Dimension(500, 425));
        settingspanel.setToolTipText("");
        SetExonColor.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                SetExonColor_actionPerformed(e);
            }
        });
        OkButton.setText("OK");
        OkButton.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                OkButton_actionPerformed(e);
            }
        });
        SetExonColor.setPreferredSize(new Dimension(200, 20));
        SetExonColor.setText("Set Exon color");
        settingspanel.setLayout(gridBagLayout1);
        ScaleInc.setText("  Scale Increment");
        ExonWidth.setToolTipText("");
        ExonWidth.setText("    Exon Width");
        this.setResizable(true);
        CheckBox.setSelected(true);
        CheckBox.setText("    Want  Base  Numbers ?");
        CheckBox.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                CheckBox_actionPerformed(e);
            }
        });
        scaleinc.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                scaleinc_actionPerformed(e);
            }
        });
        exonwidth.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                exonwidth_actionPerformed(e);
            }
        });
        SetLineColor.setText("Set Line Color");
        SetLineColor.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                SetLineColor_actionPerformed(e);
            }
        });
        LineWidth.setText("    Line Width");
        LineWidth.setToolTipText("");
        linewidth.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                linewidth_actionPerformed(e);
            }
        });
        this.getContentPane().add(settingspanel, BorderLayout.SOUTH);
        settingspanel.add(SetExonColor, new GridBagConstraints(0, 0, 3, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(26, 130, 0, 79), -9, 24));
        settingspanel.add(OkButton, new GridBagConstraints(0, 6, 3, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(19, 162, 37, 105), 84, 0));
        settingspanel.add(ScaleInc, new GridBagConstraints(0, 5, 2, 1, 0.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 137, 0, 0), 35, 35));
        settingspanel.add(scaleinc, new GridBagConstraints(1, 5, 2, 1, 1.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(20, 0, 12, 78), 26, 0));
        settingspanel.add(ExonWidth, new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(10, 130, 0, 0), 52, 35));
        settingspanel.add(SetLineColor, new GridBagConstraints(0, 2, 3, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(6, 130, 0, 75), 89, 21));
        settingspanel.add(exonwidth, new GridBagConstraints(1, 1, 2, 1, 1.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(24, 0, 15, 82), 44, 0));
        settingspanel.add(CheckBox, new GridBagConstraints(0, 4, 3, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(20, 147, 0, 64), 38, 0));
        settingspanel.add(LineWidth, new GridBagConstraints(0, 3, 2, 1, 0.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(16, 138, 0, 0), 56, 35));
        settingspanel.add(linewidth, new GridBagConstraints(2, 3, 1, 1, 1.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(30, 0, 15, 72), 44, 0));
    }

    void SetExonColor_actionPerformed(ActionEvent e) {
        exoncolor = JColorChooser.showDialog(mainFrame, "Choose Exon Color", exoncolor);
        if (exoncolor == null) {
            exoncolor = Color.blue;
        }
        genepaint.setExonColor(exoncolor);
    }

    void exonwidth_actionPerformed(ActionEvent e) {
        genepaint.setExonWidth(exonwidth.getSelectedIndex() + 1);
    }

    void SetLineColor_actionPerformed(ActionEvent e) {
        linecolor = JColorChooser.showDialog(mainFrame, "Choose  Line Color", linecolor);
        if (linecolor == null) {
            linecolor = Color.red;
        }
        genepaint.setLineColor(linecolor);
    }

    void linewidth_actionPerformed(ActionEvent e) {
        genepaint.setLineWidth(linewidth.getSelectedIndex() + 1);
    }

    void CheckBox_actionPerformed(ActionEvent e) {
        genepaint.setBaseNumbers(CheckBox.isSelected());
    }

    void scaleinc_actionPerformed(ActionEvent e) {
        genepaint.setBaseNumberInc(scaleinc.getSelectedIndex());
    }

    void OkButton_actionPerformed(ActionEvent e) {
        setVisible(false);
        repaint();
    }
}