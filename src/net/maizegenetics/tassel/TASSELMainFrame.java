/*
 * TASSEL - Trait Analysis by a aSSociation Evolution & Linkage
 * Copyright (C) 2003 Ed Buckler
 *
 * This software evaluates linkage disequilibrium nucletide diversity and
 * associations. For more information visit http://www.maizegenetics.net
 *
 * This software is distributed under GNU general public license, version 2 and without
 * any warranty or technical support.
 *
 * You can redistribute and/or modify it under the terms of GNU General
 * public license.
 *
 */
//Title:      TASSELMainFrame
//Version:
//Copyright:  Copyright (c) 1997
//Author:     Ed Buckler
//Company:    NCSU
package net.maizegenetics.tassel;


import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;

import net.maizegenetics.prefs.TasselPrefs;

import net.maizegenetics.gui.PrintTextArea;

import javax.swing.*;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Event;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Image;
import java.awt.Insets;
import java.awt.Point;
import java.awt.Toolkit;
import java.awt.event.*;
import java.awt.image.ImageProducer;

import java.io.*;

import java.net.URL;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Map;

import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import net.maizegenetics.baseplugins.ExportPlugin;
import net.maizegenetics.gui.PrintHeapAction;
import net.maizegenetics.plugindef.PluginEvent;
import net.maizegenetics.plugindef.ThreadedPluginListener;
import net.maizegenetics.progress.ProgressPanel;
import net.maizegenetics.util.Utils;
import net.maizegenetics.wizard.Wizard;
import net.maizegenetics.wizard.WizardPanelDescriptor;
import net.maizegenetics.wizard.panels.TestPanel1Descriptor;
import net.maizegenetics.wizard.panels.TestPanel2Descriptor;
import net.maizegenetics.wizard.panels.TestPanel3Descriptor;
import net.maizegenetics.wizard.panels.TestPanel4Descriptor;

import org.apache.log4j.Logger;

/**
 * TASSELMainFrame
 *
 */
public class TASSELMainFrame extends JFrame {

    private static final Logger myLogger = Logger.getLogger(TASSELMainFrame.class);
    public static final String version = "4.1.9";
    public static final String versionDate = "October 18, 2012";
    DataTreePanel theDataTreePanel;
    DataControlPanel theDataControlPanel;
    AnalysisControlPanel theAnalysisControlPanel;
    ResultControlPanel theResultControlPanel;
    private String tasselDataFile = "TasselDataFile";
    //a variable to control when the progress bar was last updated
    private long lastProgressPaint = 0;
    private String dataTreeLoadFailed = "Unable to open the saved data tree.  The file format of this version is "
            + "incompatible with other versions.";
    static final String GENOTYPE_DATA_NEEDED = "Please select genotypic data from the data tree.";
    static final String RESULTS_DATA_NEEDED = "Please select results data from the data tree.";
    JFileChooser filerSave = new JFileChooser();
    JFileChooser filerOpen = new JFileChooser();
    JPanel mainPanel = new JPanel();
    JPanel dataTreePanelPanel = new JPanel();
    JPanel reportPanel = new JPanel();
    JPanel optionsPanel = new JPanel();
    JPanel optionsPanelPanel = new JPanel();
    JPanel modeSelectorsPanel = new JPanel();
    JPanel buttonPanel = new JPanel();
    JSplitPane dataTreeReportMainPanelsSplitPanel = new JSplitPane();
    JSplitPane dataTreeReportPanelsSplitPanel = new JSplitPane();
    JScrollPane reportPanelScrollPane = new JScrollPane();
    JTextArea reportPanelTextArea = new JTextArea();
    JScrollPane mainPanelScrollPane = new JScrollPane();
    JPanel mainDisplayPanel = new JPanel();
    //mainPanelTextArea corresponds to what is called Main Panel in the user documentation
    ThreadedJTextArea mainPanelTextArea = new ThreadedJTextArea();
    JTextField statusBar = new JTextField();
    JButton resultButton = new JButton();
    JButton dataButton = new JButton();
    JButton deleteButton = new JButton();
    JButton printButton = new JButton();
    JButton analysisButton = new JButton();
    JPopupMenu mainPopupMenu = new JPopupMenu();
    JMenuBar jMenuBar = new JMenuBar();
    JMenu fileMenu = new JMenu();
    JMenu toolsMenu = new JMenu();
    JMenuItem saveMainMenuItem = new JMenuItem();
    JCheckBoxMenuItem matchCheckBoxMenuItem = new JCheckBoxMenuItem();
    JMenuItem openCompleteDataTreeMenuItem = new JMenuItem();
    JMenuItem openDataMenuItem = new JMenuItem();
    JMenuItem saveAsDataTreeMenuItem = new JMenuItem();
    JMenuItem saveCompleteDataTreeMenuItem = new JMenuItem();
    JMenuItem saveDataTreeAsMenuItem = new JMenuItem();
    JMenuItem exitMenuItem = new JMenuItem();
    JMenu helpMenu = new JMenu();
    JMenuItem helpMenuItem = new JMenuItem();
    JMenuItem preferencesMenuItem = new JMenuItem();
    JMenuItem aboutMenuItem = new JMenuItem();
    PreferencesDialog thePreferencesDialog;
    String UserComments = "";
    private final ProgressPanel myProgressPanel = ProgressPanel.getInstance();
    JButton wizardButton = new JButton();
    ExportPlugin myExportPlugin = null;

    public TASSELMainFrame(boolean debug) {
        try {
            loadSettings();
            addMenuBar();
            theDataTreePanel = new DataTreePanel(this, true, debug);
            theDataTreePanel.setToolTipText("Data Tree Panel");
            theDataControlPanel = new DataControlPanel(this, theDataTreePanel);
            theAnalysisControlPanel = new AnalysisControlPanel(this, theDataTreePanel);
            theResultControlPanel = new ResultControlPanel(this, theDataTreePanel);
            theResultControlPanel.setToolTipText("Report Panel");
            initializeMyFrame();
            setIcon();
            initDataMode();

            this.setTitle("TASSEL (Trait Analysis by aSSociation, Evolution, and Linkage) " + this.version);
            
            myLogger.info("Tassel Version: " + version + "  Date: " + versionDate);
            myLogger.info("Max Available Memory Reported by JVM: " + Utils.getMaxHeapSizeMB() + " MB");
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private void setIcon() {
        URL url = this.getClass().getResource("Logo_small.png");
        if (url == null) {
            return;
        }
        Image img = null;
        try {
            img = createImage((ImageProducer) url.getContent());
        } catch (Exception e) {
        }
        if (img != null) {
            setIconImage(img);
        }
    }

    private void initWizard() {
        Wizard wizard = new Wizard(theDataTreePanel);
        wizard.getDialog().setTitle("TASSEL Wizard (In development)");

        WizardPanelDescriptor descriptor1 = new TestPanel1Descriptor();
        wizard.registerWizardPanel(TestPanel1Descriptor.IDENTIFIER, descriptor1);

        WizardPanelDescriptor descriptor2 = new TestPanel2Descriptor();
        wizard.registerWizardPanel(TestPanel2Descriptor.IDENTIFIER, descriptor2);

        WizardPanelDescriptor descriptor3 = new TestPanel3Descriptor();
        wizard.registerWizardPanel(TestPanel3Descriptor.IDENTIFIER, descriptor3);

        WizardPanelDescriptor descriptor4 = new TestPanel4Descriptor();
        wizard.registerWizardPanel(TestPanel4Descriptor.IDENTIFIER, descriptor4);

        wizard.setCurrentPanel(TestPanel1Descriptor.IDENTIFIER);

        wizard.getDialog().setLocationRelativeTo(this);

        int ret = wizard.showModalDialog();

        // System.out.println("Dialog return code is (0=Finish,1=Cancel,2=Error): " + ret);
        // System.out.println("Second panel selection is: " +
        //    (((TestPanel2)descriptor2.getPanelComponent()).getRadioButtonSelected()));

        // System.exit(0);
    }

    private JButton getHeapButton() {

        JButton heapButton = new JButton(PrintHeapAction.getInstance(this));
        heapButton.setText("Show Memory");
        heapButton.setToolTipText("Show Memory Usage");

        return heapButton;

    }

    private void initDataMode() {

        buttonPanel.removeAll();
        
        buttonPanel.add(theDataControlPanel, BorderLayout.NORTH);
        
        dataTreePanelPanel.removeAll();

        dataTreePanelPanel.add(theDataTreePanel, BorderLayout.CENTER);

        this.validate();

        repaint();

    }

    private void initAnalysisMode() {

        buttonPanel.removeAll();

        buttonPanel.add(theAnalysisControlPanel, null);

        dataTreePanelPanel.removeAll();

        dataTreePanelPanel.add(theDataTreePanel, BorderLayout.CENTER);

        this.validate();

        repaint();
    }

    private void initResultMode() {
        buttonPanel.removeAll();
        buttonPanel.add(theResultControlPanel, null);
        dataTreePanelPanel.removeAll();
        dataTreePanelPanel.add(theDataTreePanel, BorderLayout.CENTER);
        this.validate();
        repaint();
    }

    //Component initialization
    private void initializeMyFrame() throws Exception {
        this.getContentPane().setLayout(new BorderLayout());
        Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();

        // it is time for TASSEL to claim more (almost all) of the screen real estate for itself
        // this size was selected so as to encourage the user to resize to full screen,  thereby
        // insuring that all parts of the frame are visible.
        this.setSize(new Dimension(screenSize.width * 19 / 20, screenSize.height * 19 / 20));
        this.setTitle("TASSEL (Trait Analysis by aSSociation, Evolution, and Linkage)");
        this.addWindowListener(new java.awt.event.WindowAdapter() {

            public void windowClosing(WindowEvent e) {
                this_windowClosing(e);
            }
        });

        filerSave.setDialogType(JFileChooser.SAVE_DIALOG);
        mainPanel.setLayout(new BorderLayout());
        dataTreeReportPanelsSplitPanel.setOrientation(JSplitPane.VERTICAL_SPLIT);
        optionsPanelPanel.setLayout(new BorderLayout());
        dataTreePanelPanel.setLayout(new BorderLayout());
        dataTreePanelPanel.setToolTipText("Data Tree Panel");
        reportPanel.setLayout(new BorderLayout());
        reportPanelTextArea.setEditable(false);
        reportPanelTextArea.setToolTipText("Report Panel");
        mainPanelTextArea.setDoubleBuffered(true);
        mainPanelTextArea.setEditable(false);
        mainPanelTextArea.setFont(new java.awt.Font("Monospaced", 0, 12));
        mainPanelTextArea.setToolTipText("Main Panel");
        mainPanelTextArea.addMouseListener(new java.awt.event.MouseAdapter() {

            public void mouseClicked(MouseEvent e) {
                mainTextArea_mouseClicked(e);
            }
        });
        statusBar.setBackground(Color.lightGray);
        statusBar.setBorder(null);
        statusBar.setText("Program Status");
        modeSelectorsPanel.setLayout(new GridBagLayout());
        modeSelectorsPanel.setMinimumSize(new Dimension(380, 32));
        modeSelectorsPanel.setPreferredSize(new Dimension(700, 32));
        URL imageURL = TASSELMainFrame.class.getResource("images/help1.gif");
        ImageIcon helpIcon = null;
        if (imageURL != null) {
            helpIcon = new ImageIcon(imageURL);
        }
        JButton helpButton = new JButton();
        helpButton.setIcon(helpIcon);
        helpButton.setMargin(new Insets(0, 0, 0, 0));
        helpButton.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                helpButton_actionPerformed(e);
            }
        });
        helpButton.setBackground(Color.white);
        helpButton.setMinimumSize(new Dimension(20, 20));
        helpButton.setToolTipText("Help me!!");
        resultButton.setText("Results");
        resultButton.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                resultButton_actionPerformed(e);
            }
        });
        resultButton.setMargin(new Insets(2, 2, 2, 2));
        imageURL = TASSELMainFrame.class.getResource("images/Results.gif");
        ImageIcon resultsIcon = null;
        if (imageURL != null) {
            resultsIcon = new ImageIcon(imageURL);
        }
        if (resultsIcon != null) {
            resultButton.setIcon(resultsIcon);
        }
        resultButton.setPreferredSize(new Dimension(90, 25));
        resultButton.setMinimumSize(new Dimension(87, 25));
        resultButton.setMaximumSize(new Dimension(90, 25));
        resultButton.setBackground(Color.white);

        wizardButton.setText("Wizard");
        wizardButton.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                wizardButton_actionPerformed(e);
            }
        });
        wizardButton.setMargin(new Insets(2, 2, 2, 2));

        wizardButton.setPreferredSize(new Dimension(90, 25));
        wizardButton.setMinimumSize(new Dimension(87, 25));
        wizardButton.setMaximumSize(new Dimension(90, 25));
        wizardButton.setBackground(Color.white);

        dataButton.setBackground(Color.white);
        dataButton.setMaximumSize(new Dimension(90, 25));
        dataButton.setMinimumSize(new Dimension(87, 25));
        dataButton.setPreferredSize(new Dimension(90, 25));
        imageURL = TASSELMainFrame.class.getResource("images/DataSeq.gif");
        ImageIcon dataSeqIcon = null;
        if (imageURL != null) {
            dataSeqIcon = new ImageIcon(imageURL);
        }
        if (dataSeqIcon != null) {
            dataButton.setIcon(dataSeqIcon);
        }
        dataButton.setMargin(new Insets(2, 2, 2, 2));
        dataButton.setText("Data");
        dataButton.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {

                dataButton_actionPerformed(e);
            }
        });

        printButton.setBackground(Color.white);

        printButton.setToolTipText("Print selected datum");

        imageURL = TASSELMainFrame.class.getResource("images/print1.gif");

        ImageIcon printIcon = null;

        if (imageURL != null) {
            printIcon = new ImageIcon(imageURL);
        }

        if (printIcon != null) {
            printButton.setIcon(printIcon);
        }

        printButton.setMargin(new Insets(0, 0, 0, 0));

        printButton.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {

                printButton_actionPerformed(e);
            }
        });

        analysisButton.setText("Analysis");

        analysisButton.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {

                analysisButton_actionPerformed(e);
            }
        });

        analysisButton.setMargin(new Insets(2, 2, 2, 2));

        imageURL = TASSELMainFrame.class.getResource("images/Analysis.gif");

        ImageIcon analysisIcon = null;

        if (imageURL != null) {
            analysisIcon = new ImageIcon(imageURL);
        }

        if (analysisIcon != null) {
            analysisButton.setIcon(analysisIcon);
        }

        analysisButton.setPreferredSize(new Dimension(90, 25));

        analysisButton.setMinimumSize(new Dimension(87, 25));

        analysisButton.setMaximumSize(new Dimension(90, 25));

        analysisButton.setBackground(Color.white);

        // delete button added moved from data panel by yogesh.
        deleteButton.setOpaque(true);
        deleteButton.setForeground(Color.RED);
        deleteButton.setText("Delete");
        deleteButton.setFont(new java.awt.Font("Dialog", 1, 12));
        deleteButton.setToolTipText("Delete Dataset");
        deleteButton.setMargin(new Insets(2, 2, 2, 2));


        deleteButton.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                theDataTreePanel.deleteSelectedNodes();
            }
        });


        optionsPanel.setLayout(new BorderLayout(0, 0));
        
        buttonPanel.setLayout(new BorderLayout(0, 0));

        optionsPanel.setToolTipText("Options Panel");

        mainPopupMenu.setInvoker(this);
        saveMainMenuItem.setText("Save");
        matchCheckBoxMenuItem.setText("Match");
        matchCheckBoxMenuItem.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                matchCheckBoxMenuItem_actionPerformed(e);
            }
        });
        fileMenu.setText("File");
        toolsMenu.setText("Tools");
        saveCompleteDataTreeMenuItem.setText("Save Data Tree");
        saveCompleteDataTreeMenuItem.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                saveCompleteDataTreeMenuItem_actionPerformed(e);
            }
        });
        saveDataTreeAsMenuItem.setText("Save Data Tree As ...");
        saveDataTreeAsMenuItem.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                saveDataTreeMenuItem_actionPerformed(e);
            }
        });
        openCompleteDataTreeMenuItem.setText("Open Data Tree");
        openCompleteDataTreeMenuItem.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {

                openCompleteDataTreeMenuItem_actionPerformed(e);
            }
        });

        openDataMenuItem.setText("Open Data Tree...");
        openDataMenuItem.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {

                openDataMenuItem_actionPerformed(e);
            }
        });

        saveAsDataTreeMenuItem.setText("Save Selected As...");
        saveAsDataTreeMenuItem.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                ExportPlugin plugin = getExportPlugin();
                PluginEvent event = new PluginEvent(theDataTreePanel.getSelectedTasselDataSet());
                ProgressPanel progressPanel = getProgressPanel();
                ThreadedPluginListener thread = new ThreadedPluginListener(plugin, event);
                thread.start();
                progressPanel.addPlugin(plugin);
            }
        });

        exitMenuItem.setText("Exit");
        exitMenuItem.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                exitMenuItem_actionPerformed(e);
            }
        });

        helpMenu.setText("Help");
        helpMenuItem.setText("Help Manual");
        helpMenuItem.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                helpButton_actionPerformed(e);
            }
        });

        preferencesMenuItem.setText("Set Preferences");

        preferencesMenuItem.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {

                preferencesMenuItem_actionPerformed(e);

            }
        });

        aboutMenuItem.setText("About");

        aboutMenuItem.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {

                helpAbout_actionPerformed(e);

            }
        });

        this.getContentPane().add(dataTreeReportMainPanelsSplitPanel, BorderLayout.CENTER);

        dataTreeReportMainPanelsSplitPanel.add(optionsPanelPanel, JSplitPane.TOP);

        optionsPanelPanel.add(dataTreeReportPanelsSplitPanel, BorderLayout.CENTER);

        dataTreeReportPanelsSplitPanel.add(dataTreePanelPanel, JSplitPane.TOP);

        dataTreePanelPanel.add(theDataTreePanel, BorderLayout.CENTER);

        JSplitPane reportProgress = new JSplitPane(JSplitPane.VERTICAL_SPLIT);

        reportProgress.add(reportPanel, JSplitPane.TOP);

        reportPanel.add(reportPanelScrollPane, BorderLayout.CENTER);

        reportPanelScrollPane.getViewport().add(reportPanelTextArea, null);

        reportProgress.add(new JScrollPane(myProgressPanel), JSplitPane.BOTTOM);

        dataTreeReportPanelsSplitPanel.add(reportProgress, JSplitPane.BOTTOM);

        dataTreeReportMainPanelsSplitPanel.add(mainPanel, JSplitPane.BOTTOM);

//********************************************************************************************************
        // set up so that the TASSELMainFrame can display the new alignment viewer
        //mainPanel.add(mainPanelScrollPane, BorderLayout.CENTER);
        mainPanel.add(mainDisplayPanel, BorderLayout.CENTER);
        mainDisplayPanel.setLayout(new BorderLayout());
        mainPanelScrollPane.getViewport().add(mainPanelTextArea, null);
        mainDisplayPanel.add(mainPanelScrollPane, BorderLayout.CENTER);
//********************************************************************************************************

        mainPanelScrollPane.getViewport().add(mainPanelTextArea, null);

        this.getContentPane().add(statusBar, BorderLayout.SOUTH);

        this.getContentPane().add(optionsPanel, BorderLayout.NORTH);

        optionsPanel.add(modeSelectorsPanel, BorderLayout.NORTH);

        modeSelectorsPanel.add(resultButton, new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(2, 0, 1, 0), 0, 0));

        modeSelectorsPanel.add(dataButton, new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(2, 4, 1, 0), 0, 0));

        modeSelectorsPanel.add(analysisButton, new GridBagConstraints(1, 0, 1, 1, 0.0, 0.0, GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(2, 0, 1, 0), 0, 0));

        modeSelectorsPanel.add(helpButton, new GridBagConstraints(7, 0, 1, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(2, 0, 1, 2), 0, 0));

        modeSelectorsPanel.add(printButton, new GridBagConstraints(6, 0, 1, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(2, 0, 1, 0), 0, 0));

        modeSelectorsPanel.add(deleteButton, new GridBagConstraints(3, 0, 1, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(2, 20, 1, 2), 0, 0));

        modeSelectorsPanel.add(wizardButton, new GridBagConstraints(4, 0, 1, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(2, 10, 1, 0), 0, 0));

        modeSelectorsPanel.add(getHeapButton(), new GridBagConstraints(5, 0, 1, 1, 0.0, 0.0, GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(2, 30, 1, 0), 0, 0));

        optionsPanel.add(buttonPanel, BorderLayout.SOUTH);

        mainPopupMenu.add(matchCheckBoxMenuItem);

        mainPopupMenu.add(saveMainMenuItem);



        dataTreeReportMainPanelsSplitPanel.setDividerLocation(this.getSize().width / 4);

        dataTreeReportPanelsSplitPanel.setDividerLocation((int) (this.getSize().height / 3.5));

        reportProgress.setDividerLocation((int) (this.getSize().height / 3.5));
    }

    private ExportPlugin getExportPlugin() {
        if (myExportPlugin == null) {
            myExportPlugin = new ExportPlugin(this, true);
        }
        return myExportPlugin;
    }

    private void addMenuBar() {

        jMenuBar.add(fileMenu);

        jMenuBar.add(toolsMenu);

        jMenuBar.add(helpMenu);

        fileMenu.add(saveCompleteDataTreeMenuItem);

        fileMenu.add(openCompleteDataTreeMenuItem);

        fileMenu.add(saveDataTreeAsMenuItem);

        fileMenu.add(openDataMenuItem);

        fileMenu.add(saveAsDataTreeMenuItem);

        fileMenu.addSeparator();

        fileMenu.add(exitMenuItem);


        toolsMenu.add(preferencesMenuItem);

        helpMenu.add(helpMenuItem);

        helpMenu.add(aboutMenuItem);

        this.setJMenuBar(jMenuBar);

    }

    //Help | About action performed
    private void helpAbout_actionPerformed(ActionEvent e) {

        AboutBox dlg = new AboutBox(this);

        Dimension dlgSize = dlg.getPreferredSize();

        Dimension frmSize = getSize();

        Point loc = getLocation();

        dlg.setLocation((frmSize.width - dlgSize.width) / 2 + loc.x, (frmSize.height - dlgSize.height) / 2 + loc.y);

        dlg.setModal(true);

        dlg.setVisible(true);
    }

    public void sendMessage(String text) {
        statusBar.setForeground(Color.BLACK);

        statusBar.setText(text);

    }

    public void sendErrorMessage(String text) {
        statusBar.setForeground(Color.RED);

        statusBar.setText(text);
    }

    public void setMainText(String text) {

        mainPanelTextArea.start(text);

    }

    public void setMainText(StringBuffer text) {

        mainPanelTextArea.setDoubleBuffered(false);

        mainPanelTextArea.start(text.toString());

    }

    public void setNoteText(String text) {

        reportPanelTextArea.setText(text);

    }

    private void loadSettings() {
        filerOpen.setCurrentDirectory(new File(TasselPrefs.getOpenDir()));
        filerSave.setCurrentDirectory(new File(TasselPrefs.getSaveDir()));
    }

    /**
     * Provides a save filer that remembers the last location something was saved to
     */
    public File getSaveFile() {

        File saveFile = null;
        int returnVal = filerSave.showSaveDialog(this);

        if (returnVal == JFileChooser.APPROVE_OPTION) {
            saveFile = filerSave.getSelectedFile();

            TasselPrefs.putSaveDir(filerSave.getCurrentDirectory().getPath());
        }
        return saveFile;
    }

    /**
     * Provides a open filer that remember the last location something was opened from
     */
    public File getOpenFile() {

        File openFile = null;

        int returnVal = filerOpen.showOpenDialog(this);
        System.out.println("returnVal = " + returnVal);

        System.out.println("JFileChooser.OPEN_DIALOG " + JFileChooser.OPEN_DIALOG);

        if (returnVal == JFileChooser.OPEN_DIALOG || returnVal == JFileChooser.APPROVE_OPTION) {
            openFile = filerOpen.getSelectedFile();
            System.out.println("openFile = " + openFile);
            TasselPrefs.putOpenDir(filerOpen.getCurrentDirectory().getPath());
        }

        return openFile;
    }

    void this_windowClosing(WindowEvent e) {
        exitMenuItem_actionPerformed(null);
    }

    public void addDataSet(DataSet theDataSet, String defaultNode) {
        theDataTreePanel.addDataSet(theDataSet, defaultNode);
    }

    private void saveDataTree(String file) {

        Map dataToSerialize = new LinkedHashMap();
        Map dataFromTree = theDataTreePanel.getDataList();
        StringBuilder builder = new StringBuilder();

        Iterator itr = dataFromTree.keySet().iterator();
        while (itr.hasNext()) {

            Datum currentDatum = (Datum) itr.next();
            String currentNode = (String) dataFromTree.get(currentDatum);

            try {
                ByteArrayOutputStream out = new ByteArrayOutputStream();
                ObjectOutputStream oos = new ObjectOutputStream(out);
                oos.writeObject(currentDatum);
                oos.close();
                if (out.toByteArray().length > 0) {
                    dataToSerialize.put(currentDatum, currentNode);
                }
            } catch (Exception e) {
                myLogger.warn("saveDataTree: object: " + currentDatum.getName() + " type: " + currentDatum.getData().getClass().getName() + " does not serialize.");
                myLogger.warn("saveDataTree: message: " + e.getMessage());
                if (builder.length() == 0) {
                    builder.append("Due to error, these data sets could not\n");
                    builder.append("included in the saved file...");
                }
                builder.append("Data set: ");
                builder.append(currentDatum.getName());
                builder.append(" type: ");
                builder.append(currentDatum.getData().getClass().getName());
                builder.append("\n");
            }

        }

        try {

            File theFile = new File(Utils.addSuffixIfNeeded(file, ".zip"));
            FileOutputStream fos = new FileOutputStream(theFile);
            java.util.zip.ZipOutputStream zos = new ZipOutputStream(fos);

            Map data = theDataTreePanel.getDataList();
            ZipEntry thisEntry = new ZipEntry("DATA");
            zos.putNextEntry(thisEntry);
            ObjectOutputStream oos = new ObjectOutputStream(zos);
            oos.writeObject(dataToSerialize);
            oos.flush();
            zos.closeEntry();
            fos.close();
            sendMessage("Data saved to " + theFile.getAbsolutePath());

        } catch (Exception ee) {
            sendErrorMessage("Data could not be saved: " + ee);
            ee.printStackTrace();
        }

        if (builder.length() != 0) {
            JOptionPane.showMessageDialog(this, builder.toString(), "These data sets not saved...", JOptionPane.INFORMATION_MESSAGE);
        }

    }

    private boolean readDataTree(String file) {

        boolean loadedDataTreePanel = false;
        try {

            FileInputStream fis = null;
            ObjectInputStream ois = null;
            if (file.endsWith("zip")) {

                fis = new FileInputStream(file);
                java.util.zip.ZipInputStream zis = new java.util.zip.ZipInputStream(fis);
                zis.getNextEntry();
                ois = new ObjectInputStream(zis);

            } else {

                fis = new FileInputStream(file);
                ois = new ObjectInputStream(fis);

            }

            try {
                Map data = (Map) ois.readObject();
                Iterator itr = data.keySet().iterator();
                while (itr.hasNext()) {
                    Datum currentDatum = (Datum) itr.next();
                    String currentNode = (String) data.get(currentDatum);
                    theDataTreePanel.addDatum(currentNode, currentDatum);
                }
                loadedDataTreePanel = true;
            } catch (InvalidClassException ice) {
                JOptionPane.showMessageDialog(this, dataTreeLoadFailed, "Incompatible File Format", JOptionPane.INFORMATION_MESSAGE);
            } finally {
                fis.close();
            }

            if (loadedDataTreePanel) {
                sendMessage("Data loaded.");
            }

        } catch (FileNotFoundException fnfe) {
            JOptionPane.showMessageDialog(this, "File not found: " + file, "File not found", JOptionPane.INFORMATION_MESSAGE);
            sendErrorMessage("Data tree could not be loaded.");
            return false;
        } catch (Exception ee) {
            JOptionPane.showMessageDialog(this, dataTreeLoadFailed + ee, "Incompatible File Format", JOptionPane.INFORMATION_MESSAGE);
            sendErrorMessage("Data tree could not be loaded.");
            return false;
        }

        return loadedDataTreePanel;

    }

    private void wizardButton_actionPerformed(ActionEvent e) {
        initWizard();
    }

    private void dataButton_actionPerformed(ActionEvent e) {
        initDataMode();
    }

    private void analysisButton_actionPerformed(ActionEvent e) {
        initAnalysisMode();
    }

    private void resultButton_actionPerformed(ActionEvent e) {
        initResultMode();
    }

    private void printButton_actionPerformed(ActionEvent e) {
        DataSet ds = theDataTreePanel.getSelectedTasselDataSet();
        try {
            for (int i = 0; i < ds.getSize(); i++) {
                PrintTextArea pta = new PrintTextArea(this);
                pta.printThis(ds.getData(i).getData().toString());
            }
        } catch (Exception ee) {
            System.err.println("printButton_actionPerformed:" + ee);
        }
        this.statusBar.setText("Datasets were sent to the printer");
    }

    private void helpButton_actionPerformed(ActionEvent e) {
        HelpDialog theHelpDialog = new HelpDialog(this);
        theHelpDialog.setLocationRelativeTo(this);
        theHelpDialog.setVisible(true);
    }

    private void mainTextArea_mouseClicked(MouseEvent e) {

        // For this example event, we are checking for right-mouse click.

        if (e.getModifiers() == Event.META_MASK) {

            mainPanelTextArea.add(mainPopupMenu); // JPopupMenu must be added to the component whose event is chosen.

            // Make the jPopupMenu visible relative to the current mouse position in the container.

            mainPopupMenu.show(mainPanelTextArea, e.getX(), e.getY());
        }
    }

    private void matchCheckBoxMenuItem_actionPerformed(ActionEvent e) {
        //theSettings.matchChar = matchCheckBoxMenuItem.isSelected();
    }

    private void openCompleteDataTreeMenuItem_actionPerformed(ActionEvent e) {

        String dataFileName = this.tasselDataFile + ".zip";
        File dataFile = new File(dataFileName);
        if (dataFile.exists()) {
            readDataTree(dataFileName);
        } else if (new File("QPGADataFile").exists()) {
            // this exists to maintain backward compatibility with previous versions (pre-v0.99)
            readDataTree("QPGADataFile");
        } else {
            JOptionPane.showMessageDialog(this, "File: " + dataFile.getAbsolutePath() + " does not exist.\n"
                    + "Try using File/Open Data Tree...");
        }
    }

    private void openDataMenuItem_actionPerformed(ActionEvent e) {

        File f = getOpenFile();
        if (f != null) {
            readDataTree(f.getAbsolutePath());
        }
    }

    private void saveDataTreeMenuItem_actionPerformed(ActionEvent e) {

        File f = getSaveFile();
        if (f != null) {
            saveDataTree(f.getAbsolutePath());
        }
    }

    private void saveCompleteDataTreeMenuItem_actionPerformed(ActionEvent e) {
        saveDataTree(this.tasselDataFile + ".zip");
    }

    private void exitMenuItem_actionPerformed(ActionEvent e) {
        System.exit(0);
    }

    private void preferencesMenuItem_actionPerformed(ActionEvent e) {
        if (thePreferencesDialog == null) {
            thePreferencesDialog = new PreferencesDialog();
            thePreferencesDialog.pack();
        }

        thePreferencesDialog.setLocationRelativeTo(this);
        thePreferencesDialog.setVisible(true);
    }

    public void updateMainDisplayPanel(JPanel panel) {
        mainDisplayPanel.removeAll();
        mainDisplayPanel.add(panel, BorderLayout.CENTER);
        mainDisplayPanel.repaint();
        mainDisplayPanel.validate();
    }

    public void addMenu(JMenu menu) {
        jMenuBar.add(menu);
    }

    public DataTreePanel getDataTreePanel() {
        return theDataTreePanel;
    }

    public ProgressPanel getProgressPanel() {
        return myProgressPanel;
    }
}
