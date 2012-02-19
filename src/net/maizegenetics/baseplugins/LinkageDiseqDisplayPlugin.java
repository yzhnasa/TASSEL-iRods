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
import javax.swing.event.ChangeEvent;

import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.KeyEvent;

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
                LinkageDisequilibriumComponent ldc = new LinkageDisequilibriumComponent(theLD, blockSchematic, chromosomalView, theLD.getSiteCount(), Math.max(1, theLD.getSiteCount()), Math.max(1, theLD.getSiteCount()));
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

    final int myDefaultViewableSize = 35;
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
    JScrollBar verticalScrollBar = new JScrollBar(Adjustable.VERTICAL);
    JScrollBar horizontalScrollBar = new JScrollBar(Adjustable.HORIZONTAL);
//  JComboBox formatComboBox = new JComboBox();
//  Windowed Viewer
    JPanel viewerLocationOptionsPanel = new JPanel();
    //window size
    JLabel windowSizeLabel = new JLabel();
    JSlider windowSizeSlider;
    JTextField windowSizeText = new JTextField();
    //x axis
    JSlider windowXSlider;
//  Y axis panel
    JPanel yAxisPanel = new JPanel();
    //y axis
    JSlider windowYSlider;

    public LinkageDiseqDisplayDialog(LinkageDiseqDisplayPlugin theQAF, LinkageDisequilibrium theLinkageDisequilibrium) {
        super(theQAF.getParentFrame(), "Linkage Disequilibrium", false);

        this.theLinkageDiseqDisplayPlugin = theQAF;
        this.theLinkageDisequilibrium = theLinkageDisequilibrium;
        try {
            jbInit();
            if (theLinkageDisequilibrium.getSiteCount() > myDefaultViewableSize) {
                ldFigurePanel = new LinkageDisequilibriumComponent(theLinkageDisequilibrium, true, false, myDefaultViewableSize, theLinkageDisequilibrium.getSiteCount() / 2, theLinkageDisequilibrium.getSiteCount() / 2);
            } else {
                ldFigurePanel = new LinkageDisequilibriumComponent(theLinkageDisequilibrium, true, false, theLinkageDisequilibrium.getSiteCount(), theLinkageDisequilibrium.getSiteCount() / 2, theLinkageDisequilibrium.getSiteCount() / 2);
            }
            linkPanel.add(ldFigurePanel, BorderLayout.CENTER);
            pack();
        } catch (Exception ex) {
            ex.printStackTrace();
            JOptionPane.showMessageDialog(this.getParent(), "Unable to create LD plot " + ex);
        }
        repaint();
    }

    private void jbInit() throws Exception {

        jPanel1.setLayout(gridBagLayout1);
        jPanel1.setToolTipText("");

//  Window Viewer init

        if (theLinkageDisequilibrium.getSiteCount() > myDefaultViewableSize) {

            int numSites = theLinkageDisequilibrium.getSiteCount();

            viewerLocationOptionsPanel.setLayout(new GridBagLayout());
            yAxisPanel.setLayout(new GridBagLayout());

            //Window size
            windowSizeLabel.setText("Select viewable size");
            windowSizeSlider = new JSlider(1, Math.min(numSites, 300), myDefaultViewableSize);
            windowSizeText.setText(String.valueOf(windowSizeSlider.getValue()));

            windowSizeSlider.addChangeListener(new javax.swing.event.ChangeListener() {

                public void stateChanged(ChangeEvent ce) {
                    sizeSlider_actionPerformed(ce);
                }
            });

            windowSizeText.addKeyListener(new java.awt.event.KeyListener() {

                public void keyTyped(KeyEvent ke) {
                    //                countTextField_keyTyped(ke);
                }

                public void keyPressed(KeyEvent ke) {
                    //                throw new UnsupportedOperationException("Not supported yet.");
                }

                public void keyReleased(KeyEvent ke) {
                    sizeTextField_keyTyped(ke);
                }
            });

            //X coords
            windowXSlider = new JSlider(myDefaultViewableSize / 2, numSites - (myDefaultViewableSize / 2 + (myDefaultViewableSize % 2)), numSites / 2);

            windowXSlider.addChangeListener(new javax.swing.event.ChangeListener() {

                public void stateChanged(ChangeEvent ce) {
                    xSlider_actionPerformed(ce);
                }
            });

            //Y coords
            windowYSlider = new JSlider(myDefaultViewableSize / 2, numSites - (myDefaultViewableSize / 2 + (myDefaultViewableSize % 2)), numSites / 2);
            windowYSlider.setOrientation(JSlider.VERTICAL);
            windowYSlider.setInverted(true);

            windowYSlider.addChangeListener(new javax.swing.event.ChangeListener() {

                public void stateChanged(ChangeEvent ce) {
                    ySlider_actionPerformed(ce);
                }
            });

            //Add window size to panel
            viewerLocationOptionsPanel.add(windowSizeLabel, new GridBagConstraints(0, 0, 1, 1, 0.0, 1.0, GridBagConstraints.LINE_END, GridBagConstraints.NONE, new Insets(0, 5, 0, 5), 0, 0));
            viewerLocationOptionsPanel.add(windowSizeSlider, new GridBagConstraints(1, 0, 2, 1, 0.7, 1.0, GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 0, 0));
            viewerLocationOptionsPanel.add(windowSizeText, new GridBagConstraints(3, 0, 1, 1, 0.3, 1.0, GridBagConstraints.LINE_START, GridBagConstraints.HORIZONTAL, new Insets(0, 5, 0, 30), 0, 0));

            //Add X coords to panel
//            viewerLocationOptionsPanel.add(windowXLabel, new GridBagConstraints(0, 1, 1, 1, 0.0, 1.0, GridBagConstraints.LINE_END, GridBagConstraints.NONE, new Insets(0, 5, 0, 5), 0, 0));
//            viewerLocationOptionsPanel.add(windowXSlider, new GridBagConstraints(1, 1, 2, 1, 0.7, 1.0, GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 0, 0));
//            viewerLocationOptionsPanel.add(windowXText, new GridBagConstraints(3, 1, 1, 1, 0.3, 1.0, GridBagConstraints.LINE_START, GridBagConstraints.HORIZONTAL, new Insets(0, 5, 0, 5), 0, 0));
            jPanel1.add(windowXSlider, new GridBagConstraints(0, 0, 6, 1, 1.0, 1.0, GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(0, 60, 0, 85), 0, 0));

            //Add Y coords to panel
//            viewerLocationOptionsPanel.add(windowYLabel, new GridBagConstraints(0, 2, 1, 1, 0.0, 1.0, GridBagConstraints.LINE_END, GridBagConstraints.NONE, new Insets(0, 5, 0, 5), 0, 0));
//            viewerLocationOptionsPanel.add(windowYSlider, new GridBagConstraints(1, 2, 2, 1, 0.7, 1.0, GridBagConstraints.CENTER, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 0, 0));
//            viewerLocationOptionsPanel.add(windowYText, new GridBagConstraints(3, 2, 1, 1, 0.3, 1.0, GridBagConstraints.LINE_START, GridBagConstraints.HORIZONTAL, new Insets(0, 5, 0, 5), 0, 0));
            yAxisPanel.add(windowYSlider, new GridBagConstraints(0, 0, 1, 1, 0.0, 1.0, GridBagConstraints.CENTER, GridBagConstraints.VERTICAL, new Insets(10, 0, 30, 0), 0, 0));

            getContentPane().add(viewerLocationOptionsPanel, BorderLayout.NORTH);
            getContentPane().add(yAxisPanel, BorderLayout.EAST);
        }

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
        jPanel1.add(upPRadioButton, new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(3, 30, 3, 3), 0, 0));
        jPanel1.add(upRSqrRadioButton, new GridBagConstraints(1, 1, 1, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(3, 3, 3, 3), 0, 0));
        jPanel1.add(upDPrimeRadioButton, new GridBagConstraints(2, 1, 1, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(3, 3, 3, 3), 0, 0));
        jPanel1.add(lowRSqrRadioButton, new GridBagConstraints(1, 2, 1, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(3, 3, 3, 3), 0, 0));
        jPanel1.add(lowDPrimeRadioButton, new GridBagConstraints(2, 2, 1, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(3, 3, 3, 3), 0, 0));
        jPanel1.add(lowPRadioButton, new GridBagConstraints(0, 2, 1, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(3, 30, 3, 3), 0, 0));
        jPanel1.add(chromoRadioButton, new GridBagConstraints(3, 2, 1, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
        jPanel1.add(okButton, new GridBagConstraints(5, 2, 1, 1, 0.0, 0.0, GridBagConstraints.LINE_END, GridBagConstraints.NONE, new Insets(0, 0, 0, 30), 0, 0));
        jPanel1.add(geneRadioButton, new GridBagConstraints(3, 1, 1, 1, 0.0, 0.0, GridBagConstraints.LINE_START, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
        jPanel1.add(saveButton, new GridBagConstraints(5, 1, 1, 1, 0.0, 0.0, GridBagConstraints.LINE_END, GridBagConstraints.NONE, new Insets(0, 0, 0, 30), 0, 0));
//    jPanel1.add(formatComboBox, new GridBagConstraints(5, 0, 1, 1, 0.0, 0.0
//            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(2, 6, 0, 0), 0, 0));
        jPanel1.add(schematicCheckBox, new GridBagConstraints(4, 1, 1, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
//    jPanel1.add(svgButton, new GridBagConstraints(6, 0, 1, 1, 0.0, 0.0
//            ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
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

    void sizeSlider_actionPerformed(ChangeEvent ce) {
        windowSizeText.setText(String.valueOf(windowSizeSlider.getValue()));

        windowXSlider.setValue(Math.max(1, theLinkageDisequilibrium.getSiteCount() / 2));
        windowYSlider.setValue(Math.max(1, theLinkageDisequilibrium.getSiteCount() / 2));

        int ldMeasureLower = 0;
        int ldMeasureUpper = 0;

        if (lowRSqrRadioButton.isSelected()) {
            ldMeasureLower = LinkageDisequilibriumComponent.RSQUARE;
        } else if (lowDPrimeRadioButton.isSelected()) {
            ldMeasureLower = LinkageDisequilibriumComponent.DPRIME;
        } else if (upPRadioButton.isSelected()) {
            ldMeasureLower = LinkageDisequilibriumComponent.P_VALUE;
        }

        if (upRSqrRadioButton.isSelected()) {
            ldMeasureUpper = LinkageDisequilibriumComponent.RSQUARE;
        } else if (upDPrimeRadioButton.isSelected()) {
            ldMeasureUpper = LinkageDisequilibriumComponent.DPRIME;
        } else if (upPRadioButton.isSelected()) {
            ldMeasureUpper = LinkageDisequilibriumComponent.P_VALUE;
        }

        ldFigurePanel.setWindowSize(windowSizeSlider.getValue(), ldMeasureLower, ldMeasureUpper);

        windowXSlider.setMinimum(Math.max(1, windowSizeSlider.getValue() / 2));
        windowXSlider.setMaximum(theLinkageDisequilibrium.getSiteCount() - (windowSizeSlider.getValue() / 2 + (windowSizeSlider.getValue() % 2)));

        windowYSlider.setMinimum(Math.max(1, windowSizeSlider.getValue() / 2));
        windowYSlider.setMaximum(theLinkageDisequilibrium.getSiteCount() - (windowSizeSlider.getValue() / 2 + (windowSizeSlider.getValue() % 2)));

        repaint();
    }

    void sizeTextField_keyTyped(KeyEvent ke) {
        try {
            if (!windowSizeText.getText().equals("")) {
                int value = Integer.valueOf(windowSizeText.getText());
                if (value >= windowSizeSlider.getMinimum() && value <= windowSizeSlider.getMaximum()) {
                    windowSizeSlider.setValue(value);
                } else if (value <= windowSizeSlider.getMinimum()) {
                    windowSizeSlider.setValue(windowSizeSlider.getMinimum());
                    windowSizeText.setText(String.valueOf(windowSizeSlider.getMinimum()));
                } else if (value >= windowSizeSlider.getMaximum()) {
                    windowSizeSlider.setValue(windowSizeSlider.getMaximum());
                    windowSizeText.setText(String.valueOf(windowSizeSlider.getMaximum()));
                }
            }
        } catch (NumberFormatException nfe) {
            windowSizeText.setText(String.valueOf(windowSizeSlider.getValue()));
        }
    }

    void xSlider_actionPerformed(ChangeEvent ce) {
        int ldMeasureLower = 0;
        int ldMeasureUpper = 0;

        if (lowRSqrRadioButton.isSelected()) {
            ldMeasureLower = LinkageDisequilibriumComponent.RSQUARE;
        } else if (lowDPrimeRadioButton.isSelected()) {
            ldMeasureLower = LinkageDisequilibriumComponent.DPRIME;
        } else if (upPRadioButton.isSelected()) {
            ldMeasureLower = LinkageDisequilibriumComponent.P_VALUE;
        }

        if (upRSqrRadioButton.isSelected()) {
            ldMeasureUpper = LinkageDisequilibriumComponent.RSQUARE;
        } else if (upDPrimeRadioButton.isSelected()) {
            ldMeasureUpper = LinkageDisequilibriumComponent.DPRIME;
        } else if (upPRadioButton.isSelected()) {
            ldMeasureUpper = LinkageDisequilibriumComponent.P_VALUE;
        }

        ldFigurePanel.setWindowX(windowXSlider.getValue(), ldMeasureLower, ldMeasureUpper);
        repaint();
    }

    void ySlider_actionPerformed(ChangeEvent ce) {
        int ldMeasureLower = 0;
        int ldMeasureUpper = 0;

        if (lowRSqrRadioButton.isSelected()) {
            ldMeasureLower = LinkageDisequilibriumComponent.RSQUARE;
        } else if (lowDPrimeRadioButton.isSelected()) {
            ldMeasureLower = LinkageDisequilibriumComponent.DPRIME;
        } else if (upPRadioButton.isSelected()) {
            ldMeasureLower = LinkageDisequilibriumComponent.P_VALUE;
        }

        if (upRSqrRadioButton.isSelected()) {
            ldMeasureUpper = LinkageDisequilibriumComponent.RSQUARE;
        } else if (upDPrimeRadioButton.isSelected()) {
            ldMeasureUpper = LinkageDisequilibriumComponent.DPRIME;
        } else if (upPRadioButton.isSelected()) {
            ldMeasureUpper = LinkageDisequilibriumComponent.P_VALUE;
        }

        ldFigurePanel.setWindowY(windowYSlider.getValue(), ldMeasureLower, ldMeasureUpper);
        repaint();
    }

    public int getWindowSizeSelection() {
        return windowSizeSlider.getValue();
    }

    public int getWindowXSelection() {
        return windowXSlider.getValue();
    }

    public int getWindowYSelection() {
        return windowYSlider.getValue();
    }
}
