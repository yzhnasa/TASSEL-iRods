/*
 * LinkageDiseqDisplayPlugin.java
 *
 * Created on December 22, 2006, 5:02 PM
 *
 */
package net.maizegenetics.baseplugins;

import net.maizegenetics.pal.gui.LinkageDisequilibriumComponent;
import net.maizegenetics.pal.popgen.LinkageDisequilibrium;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.net.URL;
import java.util.List;

/**
 *
 * @author Ed Buckler
 */
public class LinkageDiseqDisplayPlugin extends AbstractDisplayPlugin {

    private int upperCorner = LinkageDisequilibriumComponent.P_VALUE;
    private int lowerCorner = LinkageDisequilibriumComponent.RSQUARE;
    private boolean blockSchematic = true;
    private boolean chromosomalView = false;
    private boolean includeLabels = true;

    /** Creates a new instance of LinkageDiseqDisplayPlugin */
    public LinkageDiseqDisplayPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    public DataSet performFunction(DataSet input) {

        try {
            List<Datum> LDInList = input.getDataOfType(LinkageDisequilibrium.class);
            if (LDInList.size() != 1) {
                String message = "Invalid selection.  Please select one LD result.";
                if (isInteractive()) {
                    JOptionPane.showMessageDialog(getParentFrame(), message);
                } else {
                    System.out.println(message);
                }
                return null;
            }
            LinkageDisequilibrium theLD = (LinkageDisequilibrium) LDInList.get(0).getData();
            if (isInteractive()) {
                try {
                    LinkageDiseqDisplayDialog myDialog = new LinkageDiseqDisplayDialog(this, theLD);
                    myDialog.setLocationRelativeTo(getParentFrame());
                    myDialog.setVisible(true);
                } catch (Exception ex) {
                    ex.printStackTrace();
                    JOptionPane.showMessageDialog(this.getParentFrame(), "Unable to create LD plot " + ex);
                } catch (Error er) {
                    er.printStackTrace();
                    JOptionPane.showMessageDialog(this.getParentFrame(), "Unable to create LD plot " + er);
                }
            } else if (getSaveFile() != null) {
                LinkageDisequilibriumComponent ldc = new LinkageDisequilibriumComponent(theLD, blockSchematic, chromosomalView);
                ldc.setUpperCorner(upperCorner);
                ldc.setLowerCorner(lowerCorner);
                ldc.setSize(getImageWidth(), getImageHeight());
                ldc.setShowLabels(includeLabels);
                saveDataToFile(ldc, getSaveFile());
            }

            return null;
        } finally {
            fireProgress(100);
        }
    }

    public int getUpperCorner() {
        return upperCorner;
    }

    public void setUpperCorner(int upperCorner) {
        this.upperCorner = upperCorner;
    }

    public int getLowerCorner() {
        return lowerCorner;
    }

    public void setLowerCorner(int lowerCorner) {
        this.lowerCorner = lowerCorner;
    }

    public boolean isBlockSchematic() {
        return blockSchematic;
    }

    public void setBlockSchematic(boolean blockSchematic) {
        this.blockSchematic = blockSchematic;
    }

    public void setShowLabels(boolean includeLabels) {
        this.includeLabels = includeLabels;
    }

    public boolean isChromosomalView() {
        return chromosomalView;
    }

    public void setChromosomalView(boolean chromosomalView) {
        this.chromosomalView = chromosomalView;
    }

    /**
     * Icon for this plugin to be used in buttons, etc.
     *
     * @return ImageIcon
     */
    public ImageIcon getIcon() {
        URL imageURL = LinkageDiseqDisplayPlugin.class.getResource("images/LDPlot.gif");
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
        return "LD Plot";
    }

    /**
     * Tool Tip Text for this plugin
     *
     * @return String
     */
    public String getToolTipText() {
        return "Display Linkage Disequilibrium";
    }
}

class LinkageDiseqDisplayDialog extends JDialog {

    JPanel panel1 = new JPanel();
    JButton okButton = new JButton();
    JButton printButton = new JButton();
    JPanel linkPanel = new JPanel();
    LinkageDisequilibriumComponent ldFigurePanel;
    // net.maizegenetics.pal.popgen.LinkageDisequilibrium theLinkageDiseq;
    BorderLayout borderLayout1 = new BorderLayout();
    LinkageDisequilibrium theLinkageDisequilibrium;
    LinkageDiseqDisplayPlugin theLinkageDiseqDisplayPlugin;
    // TASSELMainFrame theTASSELMainFrame;
    BorderLayout borderLayout2 = new BorderLayout();
    JPanel jPanel1 = new JPanel();
    JButton saveButton = new JButton();
    JRadioButton upDPrimeRadioButton = new JRadioButton();
    JRadioButton upRSqrRadioButton = new JRadioButton();
    ButtonGroup theButtonGroup = new ButtonGroup();
    ButtonGroup lowerButtonGroup = new ButtonGroup();
    JRadioButton upPRadioButton = new JRadioButton();
    JRadioButton lowRSqrRadioButton = new JRadioButton();
    JRadioButton lowDPrimeRadioButton = new JRadioButton();
    JRadioButton lowPRadioButton = new JRadioButton();
    GridBagLayout gridBagLayout1 = new GridBagLayout();
    ButtonGroup geneChromoButtonGroup = new ButtonGroup();
    JRadioButton geneRadioButton = new JRadioButton();
    JRadioButton chromoRadioButton = new JRadioButton();
    JCheckBox schematicCheckBox = new JCheckBox();
    // JComboBox formatComboBox = new JComboBox();

    public LinkageDiseqDisplayDialog(LinkageDiseqDisplayPlugin theQAF, LinkageDisequilibrium theLinkageDisequilibrium) {
        super(theQAF.getParentFrame(), "Linkage Disequilibrium", false);

        this.theLinkageDiseqDisplayPlugin = theQAF;
        this.theLinkageDisequilibrium = theLinkageDisequilibrium;
        try {
            jbInit();
            ldFigurePanel = new LinkageDisequilibriumComponent(theLinkageDisequilibrium, true, false);
            linkPanel.add(ldFigurePanel, BorderLayout.CENTER);
            pack();
        } catch (Exception ex) {
            ex.printStackTrace();
            JOptionPane.showMessageDialog(this.getParent(), "Unable to create LD plot " + ex);
        }
        repaint();
    }

    private void jbInit() throws Exception {
        panel1.setLayout(borderLayout2);
        okButton.setText("Close");
        okButton.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                okButton_actionPerformed(e);
            }
        });
        linkPanel.setBorder(BorderFactory.createEtchedBorder());
        linkPanel.setLayout(borderLayout1);
        panel1.setPreferredSize(new Dimension(600, 600));
        saveButton.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                saveButton_actionPerformed(e);
            }
        });
        saveButton.setText("Save");
        upDPrimeRadioButton.setText("D\'");
        upDPrimeRadioButton.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                upDPrimeRadioButton_actionPerformed(e);
            }
        });
        upRSqrRadioButton.setSelected(true);
        upRSqrRadioButton.setText("r^2");
        upRSqrRadioButton.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                upRSqrRadioButton_actionPerformed(e);
            }
        });
        upPRadioButton.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                upPRadioButton_actionPerformed(e);
            }
        });
        upPRadioButton.setText("P");
        lowRSqrRadioButton.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                lowRSqrRadioButton_actionPerformed(e);
            }
        });
        lowRSqrRadioButton.setText("r^2");
        jPanel1.setLayout(gridBagLayout1);
        jPanel1.setToolTipText("");
        lowDPrimeRadioButton.setText("D\'");
        lowDPrimeRadioButton.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                lowDPrimeRadioButton_actionPerformed(e);
            }
        });
        lowPRadioButton.setSelected(true);
        lowPRadioButton.setText("P");
        lowPRadioButton.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                lowPRadioButton_actionPerformed(e);
            }
        });
        geneRadioButton.setSelected(true);
        geneRadioButton.setText("Gene View");
        geneRadioButton.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                geneRadioButton_actionPerformed(e);
            }
        });
        chromoRadioButton.setText("Chromo View");
        chromoRadioButton.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                chromoRadioButton_actionPerformed(e);
            }
        });
        schematicCheckBox.setSelected(true);
        schematicCheckBox.setText("Schematic");
        schematicCheckBox.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                schematicCheckBox_actionPerformed(e);
            }
        });
        // String[] s=AbstractDisplayPlugin.getPossibleGraphicOutFormats();
        // formatComboBox=new JComboBox(s);
        theButtonGroup.add(upRSqrRadioButton);
        theButtonGroup.add(upDPrimeRadioButton);
        theButtonGroup.add(upPRadioButton);
        lowerButtonGroup.add(lowRSqrRadioButton);
        lowerButtonGroup.add(lowDPrimeRadioButton);
        lowerButtonGroup.add(lowPRadioButton);
        getContentPane().add(panel1);
        panel1.add(linkPanel, BorderLayout.CENTER);
        this.getContentPane().add(jPanel1, BorderLayout.SOUTH);
        jPanel1.add(upPRadioButton, new GridBagConstraints(0, 0, 1, 2, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(3, 3, 3, 3), 0, 0));
        jPanel1.add(upRSqrRadioButton, new GridBagConstraints(1, 0, 1, 2, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(3, 3, 3, 3), 0, 0));
        jPanel1.add(upDPrimeRadioButton, new GridBagConstraints(2, 0, 1, 2, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(3, 3, 3, 3), 0, 0));
        jPanel1.add(lowRSqrRadioButton, new GridBagConstraints(1, 2, 1, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(3, 3, 3, 3), 0, 0));
        jPanel1.add(lowDPrimeRadioButton, new GridBagConstraints(2, 2, 1, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(3, 3, 3, 3), 3, -2));
        jPanel1.add(lowPRadioButton, new GridBagConstraints(0, 2, 1, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(3, 3, 3, 3), 0, 0));
        jPanel1.add(chromoRadioButton, new GridBagConstraints(3, 2, 1, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
        jPanel1.add(okButton, new GridBagConstraints(5, 1, 2, 2, 0.0, 0.0, GridBagConstraints.SOUTHEAST, GridBagConstraints.NONE, new Insets(0, 13, 12, 38), 0, 0));
        jPanel1.add(geneRadioButton, new GridBagConstraints(3, 0, 1, 2, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
        jPanel1.add(saveButton, new GridBagConstraints(4, 0, 1, 2, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
        // jPanel1.add(formatComboBox, new GridBagConstraints(5, 0, 1, 1, 0.0, 0.0
        // ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(2, 6, 0, 0), 0, 0));
        jPanel1.add(schematicCheckBox, new GridBagConstraints(4, 2, 1, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
        // jPanel1.add(svgButton, new GridBagConstraints(6, 0, 1, 1, 0.0, 0.0
        // ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
        geneChromoButtonGroup.add(chromoRadioButton);
        geneChromoButtonGroup.add(geneRadioButton);
    }

    void okButton_actionPerformed(ActionEvent e) {
        dispose();
    }/////////end calculateDisequilibrium

    void saveButton_actionPerformed(ActionEvent e) {
        //     String s=formatComboBox.getSelectedItem().toString();
        this.theLinkageDiseqDisplayPlugin.saveDataToFile(ldFigurePanel);
    }

    void upDPrimeRadioButton_actionPerformed(ActionEvent e) {
        ldFigurePanel.setUpperCorner(LinkageDisequilibriumComponent.DPRIME);
        repaint();
    }

    void upRSqrRadioButton_actionPerformed(ActionEvent e) {
        ldFigurePanel.setUpperCorner(LinkageDisequilibriumComponent.RSQUARE);
        repaint();
    }

    void upPRadioButton_actionPerformed(ActionEvent e) {
        ldFigurePanel.setUpperCorner(LinkageDisequilibriumComponent.P_VALUE);
        repaint();
    }

    void lowRSqrRadioButton_actionPerformed(ActionEvent e) {
        ldFigurePanel.setLowerCorner(LinkageDisequilibriumComponent.RSQUARE);
        repaint();
    }

    void lowPRadioButton_actionPerformed(ActionEvent e) {
        ldFigurePanel.setLowerCorner(LinkageDisequilibriumComponent.P_VALUE);
        repaint();
    }

    void lowDPrimeRadioButton_actionPerformed(ActionEvent e) {
        ldFigurePanel.setLowerCorner(LinkageDisequilibriumComponent.DPRIME);
        repaint();
    }

    void geneRadioButton_actionPerformed(ActionEvent e) {
        ldFigurePanel.setScaleOfView(false);
        repaint();
    }

    void chromoRadioButton_actionPerformed(ActionEvent e) {
        ldFigurePanel.setScaleOfView(true);
        repaint();
    }

    void schematicCheckBox_actionPerformed(ActionEvent e) {
        ldFigurePanel.setShowSchematic(schematicCheckBox.isSelected());
        repaint();
    }
}
