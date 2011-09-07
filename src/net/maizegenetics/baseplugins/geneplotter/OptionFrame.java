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

public class OptionFrame extends JInternalFrame {
  JPanel jPanel1 = new JPanel();
  JButton SaveSVGButton = new JButton();
  JButton SaveJPGButton = new JButton();
  JButton addWindowButton = new JButton();
  GeneDialog thegeneDialog;
  JComboBox TableSourceCombo = new JComboBox();
  net.maizegenetics.pal.report.TableReport[] theTableSource;
  String[] trOptions;
  GridBagLayout gridBagLayout1 = new GridBagLayout();

  public OptionFrame(GeneDialog thegeneDialog, net.maizegenetics.pal.report.TableReport[] tr) {
    this.thegeneDialog=thegeneDialog;
    this.theTableSource=tr;
    if (theTableSource == null)
        TableSourceCombo = new JComboBox();
    else {
    trOptions = new String[theTableSource.length];
    for (int i=0;i<theTableSource.length;i++)
          trOptions[i]="Diversity Data : "+(i+1);
    TableSourceCombo = new JComboBox(trOptions);
    }
    try {
      jbInit();
    }
    catch(Exception e) {
      e.printStackTrace();
    }
  }
  private void jbInit() throws Exception {
    jPanel1.setLayout(gridBagLayout1);
    SaveSVGButton.setText("Save SVG");
    SaveSVGButton.addActionListener(new java.awt.event.ActionListener() {
      public void actionPerformed(ActionEvent e) {
        SaveSVGButton_actionPerformed(e);
      }
    });
    SaveJPGButton.setText("Save JPG");
    SaveJPGButton.addActionListener(new java.awt.event.ActionListener() {
      public void actionPerformed(ActionEvent e) {
        SaveJPGButton_actionPerformed(e);
      }
    });
    addWindowButton.setText("Add GenePlot Window (Max 4)");
    addWindowButton.addActionListener(new java.awt.event.ActionListener() {
      public void actionPerformed(ActionEvent e) {
        addWindowButton_actionPerformed(e);
      }
    });
    this.setBorder(BorderFactory.createLineBorder(Color.black));
    this.getContentPane().setBackground(Color.pink);
    TableSourceCombo.addActionListener(new java.awt.event.ActionListener() {
      public void actionPerformed(ActionEvent e) {
        TableSourceCombo_actionPerformed(e);
      }
    });
    jPanel1.setMinimumSize(new Dimension(100, 50));
    jPanel1.setPreferredSize(new Dimension(787, 50));
    this.setTitle("Options Panel");
    TableSourceCombo.setPreferredSize(new Dimension(26, 25));
    jPanel1.add(SaveSVGButton,  new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 25, 24, 0), 0, 2));
    jPanel1.add(addWindowButton,  new GridBagConstraints(3, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 37, 17, 48), 31, 9));
    jPanel1.add(SaveJPGButton,  new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 63, 26, 0), 0, 0));
    jPanel1.add(TableSourceCombo,  new GridBagConstraints(2, 0, 1, 1, 1.0, 0.0
            ,GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(0, 59, 17, 0), 149, 5));
    this.getContentPane().add(jPanel1, BorderLayout.CENTER);
  }

  void addWindowButton_actionPerformed(ActionEvent e) {
       if(addWindowButton.isEnabled()) thegeneDialog.AddNewWindow();
  }

  void disableAddWindowButton () {
           addWindowButton.setEnabled(false);
  }

  void enableAddWindowButton () {
          addWindowButton.setEnabled(true);
  }

  void TableSourceCombo_actionPerformed(ActionEvent e) {
       thegeneDialog.SettrData(TableSourceCombo.getSelectedIndex());
  }

  void SaveJPGButton_actionPerformed(ActionEvent e) {
       thegeneDialog.SaveJPG();
  }

  void SaveSVGButton_actionPerformed(ActionEvent e) {
       thegeneDialog.SaveSVG();
  }
}