/*
 * LinkageDisequilibriumPlugin.java
 *
 * Created on December 22, 2006, 5:02 PM
 *
 */
package net.maizegenetics.baseplugins;

import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginEvent;
import net.maizegenetics.pal.popgen.LinkageDisequilibrium;

import javax.swing.*;

import java.awt.*;
import java.awt.event.ActionEvent;

import java.net.URL;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 *
 * @author Ed Buckler
 */
public class LinkageDisequilibriumPlugin extends AbstractPlugin {

    boolean isRapidAnalysis = true;
    int permutationNumber = 1000;
    int windowSize = 50;
    LinkageDisequilibrium.testDesign LDType = LinkageDisequilibrium.testDesign.SlidingWindow;
    int testSite = -1;

    /** Creates a new instance of LinkageDisequilibriumPlugin */
    public LinkageDisequilibriumPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    public DataSet performFunction(DataSet input) {

        try {
            List<Datum> alignInList = input.getDataOfType(Alignment.class);
            if (alignInList.size() < 1) {
                JOptionPane.showMessageDialog(getParentFrame(), "Invalid selection.  Please select sequence or marker alignment.");
                return null;
            }

            if (isInteractive()) {
                LinkageDiseqDialog myDialog = new LinkageDiseqDialog(isRapidAnalysis, permutationNumber);
                myDialog.setLocationRelativeTo(getParentFrame());
                myDialog.setVisible(true);
                if (myDialog.isCancel()) {
                    return null;
                }
                isRapidAnalysis = myDialog.getRapidLDAnalysis();
                windowSize = myDialog.getWindowSize();
                LDType = myDialog.getLDType();
            }

            List result = new ArrayList();
            Iterator<Datum> itr = alignInList.iterator();
            while (itr.hasNext()) {
                DataSet tds = null;
                Datum current = itr.next();
                tds = processDatum(current);
                if (tds != null) {
                    result.add(tds);
                    fireDataSetReturned(new PluginEvent(tds, LinkageDisequilibriumPlugin.class));
                }
            }

            return DataSet.getDataSet(result, this);
        } finally {
            fireProgress(100);
        }

    }

    private DataSet processDatum(Datum input) {
        Alignment aa = (Alignment) input.getData();
        net.maizegenetics.pal.popgen.LinkageDisequilibrium theLD = new net.maizegenetics.pal.popgen.LinkageDisequilibrium(aa, permutationNumber, windowSize, isRapidAnalysis, LDType, testSite, this);
        try {
            theLD.run();
            Datum td = new Datum("LD:" + input.getName(), theLD, "LD analysis");
            DataSet tds = new DataSet(td, this);
            return tds;
        } catch (Error e) {
            e.printStackTrace();
            StringBuilder builder = new StringBuilder();
            builder.append("Unable to run Linkage Disequilibrium analysis ");
            builder.append(e);
            String str = builder.toString();
            JOptionPane.showMessageDialog(getParentFrame(), str);
        }
        return null;
    }

    public boolean isRapidAnalysis() {
        return isRapidAnalysis;
    }

    public void setRapidAnalysis(boolean rapidAnalysis) {
        isRapidAnalysis = rapidAnalysis;
    }

    public int isPermutationNumber() {
        return permutationNumber;
    }

    public void setPermutationNumber(int permutationNumber) {
        this.permutationNumber = permutationNumber;
    }

    public void setLDType(LinkageDisequilibrium.testDesign type) {
        LDType = type;
    }

    public LinkageDisequilibrium.testDesign getLDType() {
        return LDType;
    }

    public void setWinSize(int winSize) {
        windowSize = winSize;
    }

    public int getWinSize() {
        return windowSize;
    }

    public void setTestSite(int site) {
        testSite = site;
    }

    public int getTestSite() {
        return testSite;
    }

    /**
     * Icon for this plugin to be used in buttons, etc.
     *
     * @return ImageIcon
     */
    public ImageIcon getIcon() {
        URL imageURL = LinkageDisequilibriumPlugin.class.getResource("images/LDPlot.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    /**
     * Button name for this plugin to be used in buttons, etc.
     *
     * @return String
     */
    public String getButtonName() {
        return "Link. Diseq.";
    }

    /**
     * Tool Tip Text for this plugin
     *
     * @return String
     */
    public String getToolTipText() {
        return "Linkage Disequilibrium";
    }
}

/**
 * Title:        TASSEL
 * Description:  A java program to deal with diversity
 * Copyright:    Copyright (c) 2000
 * Company:      USDA-ARS/NCSU
 * @author Ed Buckler
 * @version 1.0
 */
class LinkageDiseqDialog extends JDialog {

    boolean rapidPermutations = true;
//    int numberPermutations = 1000;
    boolean runAnalysis = false;
    JPanel panel1 = new JPanel();
    JCheckBox rapidCheckBox = new JCheckBox();
//    JTextField permNumberTextField = new JTextField();
    JLabel jLabel1 = new JLabel();
    JButton runButton = new JButton();
    JButton closeButton = new JButton();
    GridBagLayout gridBagLayout1 = new GridBagLayout();
    ButtonGroup matrixSelection = new ButtonGroup();
    JRadioButton denseMatrixButton = new JRadioButton("Full Matrix LD");
    JRadioButton slidingWindowButton = new JRadioButton("Sliding Window LD");
    JTextField windowSizeTextField = new JTextField();
    JLabel windowSizeLabel = new JLabel();
    int windowSize = 50;

    public LinkageDiseqDialog(boolean rapidPermutations, int numberPermutations) {
        super((Frame) null, "Linkage Disequilibrium", true);
        this.rapidPermutations = rapidPermutations;
//        this.numberPermutations = numberPermutations;
        try {
            jbInit();
            pack();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    void jbInit() throws Exception {
        panel1.setLayout(gridBagLayout1);
        rapidCheckBox.setSelected(rapidPermutations);
        rapidCheckBox.setText("Rapid Permutations (slightly biased p-values)");
//        permNumberTextField.setText(numberPermutations + "");
        jLabel1.setText("Permutations");

        windowSizeTextField.setText(windowSize + "");
        windowSizeLabel.setText("LD Window Size");

        denseMatrixButton.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                denseMatrixButton_actionPerformed(e);
            }
        });

        slidingWindowButton.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                slidingWindowButton_actionPerformed(e);
            }
        });

        matrixSelection.add(denseMatrixButton);
        matrixSelection.add(slidingWindowButton);
        slidingWindowButton.setSelected(true);

        runButton.setText("Run");
        runButton.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                runButton_actionPerformed(e);
            }
        });
        closeButton.setText("Close");
        closeButton.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                closeButton_actionPerformed(e);
            }
        });
        getContentPane().add(panel1);
        panel1.add(rapidCheckBox, new GridBagConstraints(0, 0, 2, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(102, 44, 0, 33), 82, 0));
//        panel1.add(permNumberTextField, new GridBagConstraints(0, 1, 1, 1, 1.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(16, 44, 0, 6), 60, -1));
//        panel1.add(jLabel1, new GridBagConstraints(1, 1, 1, 1, 0.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(13, 0, 0, 13), 60, 9));
        panel1.add(runButton, new GridBagConstraints(0, 4, 1, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(82, 55, 26, 0), 32, 5));
        panel1.add(closeButton, new GridBagConstraints(1, 4, 1, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(83, 20, 26, 156), 23, 5));
        panel1.add(windowSizeTextField, new GridBagConstraints(0, 3, 1, 1, 1.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(16, 44, 0, 6), 60, -1));
        panel1.add(windowSizeLabel, new GridBagConstraints(1, 3, 1, 1, 0.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(13, 0, 0, 13), 60, 9));
        panel1.add(denseMatrixButton, new GridBagConstraints(0, 1, 1, 1, 1.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(16, 44, 0, 0), 60, -1));
        panel1.add(slidingWindowButton, new GridBagConstraints(0, 2, 1, 1, 1.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 44, 0, 6), 60, -1));
    }

    private void denseMatrixButton_actionPerformed(ActionEvent e) {
        windowSizeTextField.setVisible(false);
        windowSizeLabel.setVisible(false);
    }

    private void slidingWindowButton_actionPerformed(ActionEvent e) {
        windowSizeTextField.setVisible(true);
        windowSizeLabel.setVisible(true);
    }

    public LinkageDisequilibrium.testDesign getLDType() {
        if (denseMatrixButton.isSelected()) {
            return LinkageDisequilibrium.testDesign.All;
        } else {
            return LinkageDisequilibrium.testDesign.SlidingWindow;
        }
    }

    /** Returns whether the run button was chosen*/
    public boolean isRunAnalysis() {
        return runAnalysis;
    }

    /** Returns whether the run button was chosen*/
    public boolean isCancel() {
        return !runAnalysis;
    }

    /** Returns whether the rapid permutation approach should be used */
    public boolean getRapidLDAnalysis() {
        return rapidPermutations;
    }

//    /** Return the number of permutations to use*/
//    public int getPermutationNumber() {
//        return numberPermutations;
//    }
    /** Retrun the window size*/
    public int getWindowSize() {
        return windowSize;
    }

    void runButton_actionPerformed(ActionEvent e) {
        rapidPermutations = rapidCheckBox.isSelected();
//        try {
//            numberPermutations = Integer.parseInt(permNumberTextField.getText());
//        } catch (Exception ee) {
//            permNumberTextField.setText("Set Integer");
//            return;
//        }
        try {
            windowSize = Integer.parseInt(windowSizeTextField.getText());
        } catch (Exception ee) {
            windowSizeTextField.setText("Set Integer");
        }
        runAnalysis = true;
        setVisible(false);
    }

    void closeButton_actionPerformed(ActionEvent e) {
        runAnalysis = false;
        setVisible(false);
    }
}
