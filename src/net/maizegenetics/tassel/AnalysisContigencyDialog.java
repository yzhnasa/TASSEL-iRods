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
// This snippet creates a new dialog box
// with buttons on the bottom.

//Title:
//Version:
//Copyright:
//Author:
//Company:
//Description:

package net.maizegenetics.tassel;

import net.maizegenetics.pal.statistics.ContigencyTable;
import net.maizegenetics.pal.statistics.FisherExact;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.File;
import java.util.Vector;
import net.maizegenetics.prefs.TasselPrefs;

public class AnalysisContigencyDialog extends JDialog {
    int rows = 2, columns = 2, reps;
    JTextField[][] entryTextField;
    //JdbNavComboBox[] popJdbNavComboBox;
    // QueryDataSet[] popQueryDataSet;
    Vector[] stockIDInPop;
    JPanel mainPanel = new JPanel();
    JPanel controlPanel = new JPanel();
    JPanel resultsPanel = new JPanel();
    JButton okButton = new JButton();
    JButton cancelButton = new JButton();

    JTextField txtLeftTailed = new JTextField();
    JTextField txtRightTailed = new JTextField();
    JTextField txtTwoTailed = new JTextField();
    JPanel buttonPanel = new JPanel();
    GridLayout buttonPanelGridLayout = new GridLayout();
    JLabel lblContingencyVsFisher = new JLabel();
    JPanel contigPanel = new JPanel();
    GridLayout contigPanelGridLayout = new GridLayout();
//  JdbComboBox rowComboBox = new JdbComboBox();
    JLabel lblRows = new JLabel();
    // Database database1 = new Database();
    ButtonGroup theButtonGroup = new ButtonGroup();
    String[] comboPoss = {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10"};
    JComboBox columnComboBox = new JComboBox(comboPoss);
    JComboBox rowComboBox = new JComboBox(comboPoss);
    JLabel lblColumn = new JLabel();
    JTextField repTextField = new JTextField();
    JLabel lblReps = new JLabel();
    JLabel resultLabel = new JLabel();
    GridBagLayout gridBagLayout2 = new GridBagLayout();

    public AnalysisContigencyDialog(JFrame frame, String title, boolean modal) {
        super(frame, title, modal);
        try {
            jbInit();
        } catch (Exception e) {
            e.printStackTrace();
        }
        pack();
    }


    private void jbInit() throws Exception {
        lblColumn.setText("Columns");
        columnComboBox.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                columnComboBox_actionPerformed(e);
            }
        });

        String[] petStrings = {"Bird", "Cat", "Dog", "Rabbit", "Pig"};

//Create the combo box, select item at index 4.
//Indices start at 0, so 4 specifies the pig.
        JComboBox petList = new JComboBox(petStrings);


        //   columnComboBox.setItems(new String[] {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10"});
        columnComboBox.setSelectedIndex(1);
        buttonPanel.setLayout(buttonPanelGridLayout);
        controlPanel.setLayout(gridBagLayout2);
        okButton.setText("Calculate");
        okButton.addActionListener(new ContigencyDialog_okButton_actionAdapter(this));
        cancelButton.setText("Cancel");
        buttonPanelGridLayout.setHgap(4);
        cancelButton.addActionListener(new ContigencyDialog_cancelButton_actionAdapter(this));
        this.addWindowListener(new ContigencyDialog_this_windowAdapter(this));
        mainPanel.setLayout(new GridBagLayout());

        lblContingencyVsFisher.setFont(new java.awt.Font("Dialog", 1, 18));
        lblContingencyVsFisher.setText("Contigency Or Fisher Exact Test");

        contigPanel.setLayout(contigPanelGridLayout);
        contigPanelGridLayout.setColumns(2);
        contigPanelGridLayout.setRows(3);

        rowComboBox.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                rowComboBox_actionPerformed(e);
            }
        });


        buttonPanel.add(okButton, null);
        buttonPanel.add(cancelButton, null);

        // setup results panel
        resultsPanel.setLayout(new GridLayout(4, 2));

        JLabel lblLeftTailed = new JLabel("Left-Tailed P: ", SwingConstants.RIGHT);
        lblLeftTailed.setLabelFor(txtLeftTailed);

        JLabel lblRightTailed = new JLabel("Right-Tailed P: ", SwingConstants.RIGHT);
        lblRightTailed.setLabelFor(txtRightTailed);

        JLabel lblTwoTailed = new JLabel("Two-Tailed P: ", SwingConstants.RIGHT);
        lblTwoTailed.setLabelFor(lblTwoTailed);

        resultsPanel.add(lblLeftTailed);
        resultsPanel.add(txtLeftTailed);
        resultsPanel.add(lblRightTailed);
        resultsPanel.add(txtRightTailed);
        resultsPanel.add(lblTwoTailed);
        resultsPanel.add(txtTwoTailed);

//    rowComboBox.setItems(new String[] {"1", "2", "3", "4", "5", "6", "7", "8", "9", "10"});
        rowComboBox.setSelectedIndex(1);
        lblRows.setText("Rows");
        repTextField.setText("5000");
        lblReps.setText("Reps");
        resultLabel.setFont(new java.awt.Font("Dialog", 1, 12));
        resultLabel.setForeground(Color.red);
        resultLabel.setText("Enter Data & Hit Calculate");

        // setup the panel containing controls for setting number of rows and columns
        controlPanel.add(lblContingencyVsFisher, new GridBagConstraints(0, 0, 3, 1, 0.0, 0.0
                , GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 13, 15, 164), 6, 19));
        controlPanel.add(rowComboBox, new GridBagConstraints(0, 1, 1, 1, 1.0, 0.0
                , GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(0, 32, 0, 0), 63, 5));
        controlPanel.add(lblRows, new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
                , GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(31, 32, 0, 0), 72, 12));
        controlPanel.add(columnComboBox, new GridBagConstraints(1, 1, 1, 1, 1.0, 0.0
                , GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(0, 13, 0, 0), 63, 5));
        controlPanel.add(lblColumn, new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0
                , GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(28, 13, 0, 0), 60, 12));
        controlPanel.add(repTextField, new GridBagConstraints(2, 1, 1, 1, 1.0, 0.0
                , GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 28, 0, 107), 60, 6));
        controlPanel.add(lblReps, new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0
                , GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(28, 28, 0, 141), 31, 12));
        controlPanel.add(resultLabel, new GridBagConstraints(0, 2, 3, 1, 0.0, 0.0
                , GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(17, 27, 0, 90), 204, 8));

        // control panel holds the panel holding the data input fields
        controlPanel.add(contigPanel, new GridBagConstraints(0, 3, 3, 1, 1.0, 1.0
                , GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(4, 20, 4, 20), 395, 237));

        // layout the highest level panel
        mainPanel.add(controlPanel, new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
                , GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(0, 0, 0, 0), 0, 0));
        mainPanel.add(resultsPanel, new GridBagConstraints(0, 1, 2, 2, 1.0, 0.0
                , GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(4, 20, 4, 20), 0, 0));
        mainPanel.add(buttonPanel, new GridBagConstraints(0, 3, 1, 1, 1.0, 0.0
                , GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(4, 8, 4, 8), 0, 0));


//    init_Pops();
        getContentPane().add(mainPanel);
    }

    void init_Pops() {
        try {
            if(rows == 2 && columns == 2){
                repTextField.setEnabled(false);
            }else{
                repTextField.setEnabled(true);
            }
            contigPanel.removeAll();
            contigPanelGridLayout.setRows(rows);
            contigPanelGridLayout.setColumns(columns);
            entryTextField = new JTextField[rows][columns];
            for (int r = 0; r < rows; r++)
                for (int c = 0; c < columns; c++) {//popCheckbox[i]=new JCheckBox("Pop:"+i);
                    entryTextField[r][c] = new JTextField("0");
                    entryTextField[r][c].setHorizontalAlignment(JTextField.CENTER);
                    contigPanel.add(entryTextField[r][c], null);
                }
        } catch (Exception e) {
            System.out.println("Contigency:" + e);
        }
        pack();
    }

    // OK
    void okButton_actionPerformed(ActionEvent e) {
        Integer ii = new Integer(0);
        double p, leftTailedP, rightTailedP, twoTailedP;

        int[][] contig = new int[rows][columns];
        int sum = 0;
        try {
            reps = ii.parseInt(repTextField.getText());
            for (int r = 0; r < rows; r++)
                for (int c = 0; c < columns; c++) {
                    contig[r][c] = ii.parseInt(entryTextField[r][c].getText());
                    sum += contig[r][c];
                }
        }//end of try
        catch (Exception ee) {
            resultLabel.setText("One of the values was bad");
        }
        if ((rows == 2) && (columns == 2)) {
            FisherExact fisherExact = new FisherExact(sum + 10);
            p = fisherExact.getCumlativeP(contig[0][0], contig[1][0], contig[0][1], contig[1][1]);
            leftTailedP = fisherExact.getLeftTailedP(contig[0][0], contig[1][0], contig[0][1], contig[1][1]);
            rightTailedP = fisherExact.getRightTailedP(contig[0][0], contig[1][0], contig[0][1], contig[1][1]);
            twoTailedP = fisherExact.getTwoTailedP(contig[0][0], contig[1][0], contig[0][1], contig[1][1]);

            resultLabel.setText("Fisher Exact P=" + p);
            txtLeftTailed.setText(new Double(leftTailedP).toString());
            txtRightTailed.setText(new Double(rightTailedP).toString());
            txtTwoTailed.setText(new Double(twoTailedP).toString());
        } else {
            ContigencyTable contigencyTable = new ContigencyTable(sum + 10);
            contigencyTable.setMatrix(contig);
            p = contigencyTable.calcMonteCarloExactTest(reps);
            resultLabel.setText("Contigency Table Permutation P=" + p);
            txtLeftTailed.setText("");
            txtRightTailed.setText("");
            txtTwoTailed.setText("");
        }
    }

    public void showDialog() {
        getSettings();
        setVisible(true);
    }

    void getSettings() {
        //allStockCheckBox.setSelected(theSettings.allStockIncluded);
    }

    File getOutputFile() {
        JFileChooser filerSave = new JFileChooser(TasselPrefs.getSaveDir());
        filerSave.showSaveDialog(this);
        File saveFile = filerSave.getSelectedFile();
        return saveFile;
    }

    // Cancel
    void cancelButton_actionPerformed(ActionEvent e) {
        dispose();
    }

    void this_windowClosing(WindowEvent e) {
        dispose();
    }

    void rowComboBox_actionPerformed(ActionEvent e) {
        rows = rowComboBox.getSelectedIndex() + 1;
//  System.out.println("rows="+rows);
        init_Pops();
    }

    void columnComboBox_actionPerformed(ActionEvent e) {
        columns = columnComboBox.getSelectedIndex() + 1;
//  System.out.println("rows="+rows);
        init_Pops();
    }

class ContigencyDialog_okButton_actionAdapter implements ActionListener {
    AnalysisContigencyDialog adaptee;

    ContigencyDialog_okButton_actionAdapter(AnalysisContigencyDialog adaptee) {
        this.adaptee = adaptee;
    }

    public void actionPerformed(ActionEvent e) {
        adaptee.okButton_actionPerformed(e);
    }
}

class ContigencyDialog_cancelButton_actionAdapter implements ActionListener {
    AnalysisContigencyDialog adaptee;

    ContigencyDialog_cancelButton_actionAdapter(AnalysisContigencyDialog adaptee) {
        this.adaptee = adaptee;
    }

    public void actionPerformed(ActionEvent e) {
        adaptee.cancelButton_actionPerformed(e);
    }
}

class ContigencyDialog_this_windowAdapter extends WindowAdapter {
    AnalysisContigencyDialog adaptee;

    ContigencyDialog_this_windowAdapter(AnalysisContigencyDialog adaptee) {
        this.adaptee = adaptee;
    }

    public void windowClosing(WindowEvent e) {
        adaptee.this_windowClosing(e);
    }
}
}


