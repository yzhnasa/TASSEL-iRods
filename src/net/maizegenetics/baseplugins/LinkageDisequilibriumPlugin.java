/*
 * LinkageDisequilibriumPlugin.java
 *
 * Created on December 22, 2006, 5:02 PM
 *
 */
package net.maizegenetics.baseplugins;

import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.FilterAlignment;

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
import java.util.List;

/**
 *
 * @author Ed Buckler
 */
public class LinkageDisequilibriumPlugin extends AbstractPlugin {

    private boolean myIsRapidAnalysis = true;
    private int myPermutationNumber = 1000;
    private int myWindowSize = 50;
    private LinkageDisequilibrium.testDesign myLDType = LinkageDisequilibrium.testDesign.SlidingWindow;
    private int myTestSite = -1;
    private int myNumAccumulateIntervals = 100;
    private boolean myIsAccumulateResults = false;
    private FilterAlignment myPossibleAlignmentForSiteList;
    private String myPossibleAlignmentName;
    private int[] myPossibleSiteList;

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

            Datum current = alignInList.get(0);

            if (alignInList.size() > 1) {
                try {
                    FilterAlignment temp = (FilterAlignment) alignInList.get(1).getData();
                    if (temp.getBaseAlignment() == current.getData()) {
                        myPossibleAlignmentForSiteList = temp;
                        myPossibleAlignmentName = alignInList.get(1).getName();
                    }
                } catch (Exception e) {
                    // do nothing
                }
            }

            if (isInteractive()) {
                LinkageDiseqDialog myDialog = new LinkageDiseqDialog(myIsRapidAnalysis, myPossibleAlignmentName);
                myDialog.setLocationRelativeTo(getParentFrame());
                myDialog.setVisible(true);
                if (myDialog.isCancel()) {
                    return null;
                }
                myLDType = myDialog.getLDType();

                myIsRapidAnalysis = myDialog.getRapidLDAnalysis();

                if (myLDType == LinkageDisequilibrium.testDesign.SlidingWindow) {
                    myWindowSize = myDialog.getWindowSize();
                }

                myIsAccumulateResults = myDialog.isAccumulateResults();
                if (myIsAccumulateResults) {
                    myNumAccumulateIntervals = myDialog.getNumAccumulateIntervals();
                }

                if (myLDType == LinkageDisequilibrium.testDesign.SiteByAll) {
                    myTestSite = myDialog.getTestSite();

                }

                if (myLDType == LinkageDisequilibrium.testDesign.SiteList) {
                    myPossibleSiteList = myPossibleAlignmentForSiteList.getBaseSitesShown();

                }
            }

            List result = new ArrayList();
            DataSet tds = null;

            tds = processDatum(current);
            if (tds != null) {
                result.add(tds);
                fireDataSetReturned(new PluginEvent(tds, LinkageDisequilibriumPlugin.class));
            }

            return DataSet.getDataSet(result, this);
        } finally {
            fireProgress(100);
        }

    }

    private DataSet processDatum(Datum input) {
        Alignment aa = (Alignment) input.getData();
        LinkageDisequilibrium theLD = new LinkageDisequilibrium(aa, myPermutationNumber, myWindowSize, myIsRapidAnalysis, myLDType, myTestSite, this, myIsAccumulateResults, myNumAccumulateIntervals, myPossibleSiteList);
        try {
            theLD.run();
            Datum td = new Datum("LD:" + input.getName(), theLD, "LD Analysis");
            DataSet tds = new DataSet(td, this);
            return tds;
        } catch (Exception e) {
            e.printStackTrace();
            StringBuilder builder = new StringBuilder();
            builder.append("Unable to run Linkage Disequilibrium analysis ");
            builder.append(e.getMessage());
            String str = builder.toString();
            JOptionPane.showMessageDialog(getParentFrame(), str);
        }
        return null;
    }

    public boolean isRapidAnalysis() {
        return myIsRapidAnalysis;
    }

    public void setRapidAnalysis(boolean rapidAnalysis) {
        myIsRapidAnalysis = rapidAnalysis;
    }

    public int isPermutationNumber() {
        return myPermutationNumber;
    }

    public void setPermutationNumber(int permutationNumber) {
        myPermutationNumber = permutationNumber;
    }

    public void setLDType(LinkageDisequilibrium.testDesign type) {
        myLDType = type;
    }

    public LinkageDisequilibrium.testDesign getLDType() {
        return myLDType;
    }

    public void setWinSize(int winSize) {
        myWindowSize = winSize;
    }

    public int getWinSize() {
        return myWindowSize;
    }

    public void setNumAccumulateIntervals(int numIntervals) {
        myNumAccumulateIntervals = numIntervals;
    }

    public int getNumAccumulateIntervals() {
        return myNumAccumulateIntervals;
    }

    public void setIsAccumulateResults(boolean accumulate) {
        myIsAccumulateResults = accumulate;
    }

    public boolean getIsAccumulateResults() {
        return myIsAccumulateResults;
    }

    public void setTestSite(int site) {
        myTestSite = site;
    }

    public int getTestSite() {
        return myTestSite;
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

    private boolean myIsRapidPermutations = true;
    // private int numberPermutations = 1000;
    private boolean myRunAnalysis = false;
    private JCheckBox myRapidCheckBox = new JCheckBox();
    // private JTextField permNumberTextField = new JTextField();
    // private JLabel jLabel1 = new JLabel();
    private JButton myRunButton = new JButton("Run");
    private JButton myCloseButton = new JButton("Close");
    private JRadioButton myDenseMatrixButton = new JRadioButton("Full Matrix LD");
    private JRadioButton mySlidingWindowButton = new JRadioButton("Sliding Window LD");
    private JTextField myWindowSizeTextField = new JTextField();
    private JLabel myWindowSizeLabel = new JLabel("LD Window Size");
    private int myWindowSize = 50;
    private JRadioButton mySiteByAllButton = new JRadioButton("Site By All LD");
    private JTextField mySiteByAllTextField = new JTextField();
    private JLabel mySiteByAllLabel = new JLabel("Site");
    private int myTestSite = 0;
    private JRadioButton mySiteListButton = new JRadioButton();
    private JRadioButton myAccumulativeResultsButton = new JRadioButton("Accumulate R2 Results");
    private JTextField myAccumulativeResultsTextField = new JTextField();
    private JLabel myAccumulativeResultsLabel = new JLabel("Number Of Intervals");
    private int myNumAccumulativeInterval = 100;
    private String myAlignmentForSiteList;

    public LinkageDiseqDialog(boolean rapidPermutations, String alignmentForSiteList) {
        super((Frame) null, "Linkage Disequilibrium", true);
        myIsRapidPermutations = rapidPermutations;
        myAlignmentForSiteList = alignmentForSiteList;
        // numberPermutations = numberPermutations;
        try {
            jbInit();
            pack();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    private void jbInit() throws Exception {
        JPanel panel1 = new JPanel();
        panel1.setLayout(new GridBagLayout());
        myRapidCheckBox.setSelected(myIsRapidPermutations);
        myRapidCheckBox.setText("Rapid Permutations (slightly biased p-values)");
        // permNumberTextField.setText(numberPermutations + "");
        // jLabel1.setText("Permutations");

        myDenseMatrixButton.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                hideOptions();
            }
        });

        myWindowSizeTextField.setText(String.valueOf(myWindowSize));

        mySlidingWindowButton.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                hideOptions();
                myWindowSizeTextField.setVisible(true);
                myWindowSizeLabel.setVisible(true);
            }
        });

        mySiteByAllTextField.setText(String.valueOf(myTestSite));

        mySiteByAllButton.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                hideOptions();
                mySiteByAllTextField.setVisible(true);
                mySiteByAllLabel.setVisible(true);
            }
        });

        mySiteListButton.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                hideOptions();
            }
        });

        if (myAlignmentForSiteList != null) {
            mySiteListButton.setText("Sites From: " + myAlignmentForSiteList);
            mySiteListButton.setVisible(true);
        } else {
            mySiteListButton.setVisible(false);
        }

        myAccumulativeResultsButton.setSelected(false);
        myAccumulativeResultsTextField.setText(String.valueOf(myNumAccumulativeInterval));

        ButtonGroup matrixSelection = new ButtonGroup();
        matrixSelection.add(myDenseMatrixButton);
        matrixSelection.add(mySlidingWindowButton);
        matrixSelection.add(mySiteByAllButton);
        matrixSelection.add(mySiteListButton);
        mySlidingWindowButton.setSelected(true);
        hideOptions();
        myWindowSizeTextField.setVisible(true);
        myWindowSizeLabel.setVisible(true);

        myRunButton.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                runButton_actionPerformed(e);
            }
        });

        myCloseButton.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                closeButton_actionPerformed(e);
            }
        });

        getContentPane().add(panel1);
        panel1.add(myRapidCheckBox, new GridBagConstraints(0, 0, 2, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(102, 44, 0, 33), 82, 0));
        // panel1.add(permNumberTextField, new GridBagConstraints(0, 1, 1, 1, 1.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(16, 44, 0, 6), 60, -1));
        // panel1.add(jLabel1, new GridBagConstraints(1, 1, 1, 1, 0.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(13, 0, 0, 13), 60, 9));
        panel1.add(myDenseMatrixButton, new GridBagConstraints(0, 1, 1, 1, 1.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(16, 44, 0, 0), 60, -1));

        panel1.add(mySlidingWindowButton, new GridBagConstraints(0, 2, 1, 1, 1.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(25, 44, 0, 6), 60, -1));
        panel1.add(myWindowSizeTextField, new GridBagConstraints(0, 3, 1, 1, 1.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(16, 44, 0, 6), 60, -1));
        panel1.add(myWindowSizeLabel, new GridBagConstraints(1, 3, 1, 1, 0.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(13, 0, 0, 13), 60, 9));

        panel1.add(mySiteByAllButton, new GridBagConstraints(0, 4, 1, 1, 1.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(25, 44, 0, 6), 60, -1));
        panel1.add(mySiteByAllTextField, new GridBagConstraints(0, 5, 1, 1, 1.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(16, 44, 0, 6), 60, -1));
        panel1.add(mySiteByAllLabel, new GridBagConstraints(1, 5, 1, 1, 0.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(13, 0, 0, 13), 60, 9));

        panel1.add(mySiteListButton, new GridBagConstraints(0, 6, 1, 1, 1.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(25, 44, 0, 6), 60, -1));

        panel1.add(myAccumulativeResultsButton, new GridBagConstraints(0, 7, 1, 1, 1.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(25, 44, 0, 6), 60, -1));
        panel1.add(myAccumulativeResultsTextField, new GridBagConstraints(0, 8, 1, 1, 1.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(16, 44, 0, 6), 60, -1));
        panel1.add(myAccumulativeResultsLabel, new GridBagConstraints(1, 8, 1, 1, 0.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(13, 0, 0, 13), 60, 9));

        panel1.add(myRunButton, new GridBagConstraints(0, 9, 1, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(82, 55, 26, 0), 32, 5));
        panel1.add(myCloseButton, new GridBagConstraints(1, 9, 1, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(83, 20, 26, 156), 23, 5));
    }

    private void hideOptions() {
        myWindowSizeTextField.setVisible(false);
        myWindowSizeLabel.setVisible(false);
        mySiteByAllTextField.setVisible(false);
        mySiteByAllLabel.setVisible(false);
    }

    public boolean isAccumulateResults() {
        return myAccumulativeResultsButton.isSelected();
    }

    public LinkageDisequilibrium.testDesign getLDType() {
        if (myDenseMatrixButton.isSelected()) {
            return LinkageDisequilibrium.testDesign.All;
        } else if (mySlidingWindowButton.isSelected()) {
            return LinkageDisequilibrium.testDesign.SlidingWindow;
        } else if (mySiteByAllButton.isSelected()) {
            return LinkageDisequilibrium.testDesign.SiteByAll;
        } else if (mySiteListButton.isSelected()) {
            return LinkageDisequilibrium.testDesign.SiteList;
        } else {
            throw new IllegalStateException("LinkageDisequilibriumPlugin: getLDType: No known LD Type selected.");
        }
    }

    /** Returns whether the run button was chosen*/
    public boolean isRunAnalysis() {
        return myRunAnalysis;
    }

    /** Returns whether the run button was chosen*/
    public boolean isCancel() {
        return !myRunAnalysis;
    }

    /** Returns whether the rapid permutation approach should be used */
    public boolean getRapidLDAnalysis() {
        return myIsRapidPermutations;
    }

    /** Return the window size */
    public int getWindowSize() {
        return myWindowSize;
    }

    public int getTestSite() {
        return myTestSite;
    }

    public int getNumAccumulateIntervals() {
        return myNumAccumulativeInterval;
    }

    void runButton_actionPerformed(ActionEvent e) {

        myIsRapidPermutations = myRapidCheckBox.isSelected();
        //        try {
        //            numberPermutations = Integer.parseInt(permNumberTextField.getText());
        //        } catch (Exception ee) {
        //            permNumberTextField.setText("Set Integer");
        //            return;
        //        }

        if (getLDType() == LinkageDisequilibrium.testDesign.SlidingWindow) {
            try {
                myWindowSize = Integer.parseInt(myWindowSizeTextField.getText());
            } catch (Exception ee) {
                myWindowSizeTextField.setText("Set Integer");
                return;
            }
        }

        if (getLDType() == LinkageDisequilibrium.testDesign.SiteByAll) {
            try {
                myTestSite = Integer.parseInt(mySiteByAllTextField.getText());
            } catch (Exception ee) {
                mySiteByAllTextField.setText("Set Integer");
                return;
            }
        }

        if (isAccumulateResults()) {
            try {
                myNumAccumulativeInterval = Integer.parseInt(myAccumulativeResultsTextField.getText());
            } catch (Exception ee) {
                myAccumulativeResultsTextField.setText("Set Integer");
                return;
            }
        }
        myRunAnalysis = true;
        setVisible(false);
    }

    void closeButton_actionPerformed(ActionEvent e) {
        myRunAnalysis = false;
        setVisible(false);
    }
}