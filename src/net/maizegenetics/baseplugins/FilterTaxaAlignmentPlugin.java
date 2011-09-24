/*
 * FilterTaxaAlignmentPlugin.java
 *
 * Created on December 22, 2006, 5:02 PM
 *
 */
package net.maizegenetics.baseplugins;

import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.FilterAlignment;
import net.maizegenetics.pal.alignment.FilterPhenotype;
import net.maizegenetics.pal.alignment.Phenotype;

import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.Identifier;
import net.maizegenetics.pal.ids.SimpleIdGroup;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginEvent;

import javax.swing.*;

import java.awt.BorderLayout;
import java.awt.Frame;
import java.awt.event.ActionEvent;
import java.awt.event.KeyEvent;

import java.net.URL;

import java.util.ArrayList;
import java.util.List;

import org.apache.log4j.Logger;

/**
 *
 * @author Ed Buckler
 */
public class FilterTaxaAlignmentPlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(FilterTaxaAlignmentPlugin.class);
    private IdGroup myIdsToKeep = null;
    private IdGroup myIdsToRemove = null;

    /** Creates a new instance of FilterTaxaAlignmentPlugin */
    public FilterTaxaAlignmentPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    public DataSet performFunction(DataSet input) {

        try {

            List inputData = input.getDataSet();
            if (inputData.size() != 1) {
                if (isInteractive()) {
                    JOptionPane.showMessageDialog(getParentFrame(), "Invalid selection. Please select a single sequence or phenotype.");
                } else {
                    myLogger.error("performFunction: Please input a single sequence or phenotype.");
                }
                return null;
            }

            Datum inputDatum = (Datum) inputData.get(0);

            if (!(inputDatum.getData() instanceof Alignment) && !(inputDatum.getData() instanceof Phenotype)) {
                if (isInteractive()) {
                    JOptionPane.showMessageDialog(getParentFrame(), "Invalid selection. Please select a single sequence or phenotype.");
                } else {
                    myLogger.error("performFunction: Please input a single sequence or phenotype.");
                }
                return null;
            }

            Datum td = processDatum(inputDatum, isInteractive());
            if (td == null) {
                return null;
            }

            DataSet output = new DataSet(td, this);
            //I am setting the firing class as the metadata - so that the control panel know where the event is coming from
            fireDataSetReturned(new PluginEvent(output, FilterTaxaAlignmentPlugin.class));

            return output;

        } finally {
            fireProgress(100);
        }
    }

    private Datum processDatum(Datum inDatum, boolean isInteractive) {

        Object theData = inDatum.getData();

        if (isInteractive) {
            DataRowFilterDialog myDialog = null;
            if (theData instanceof Alignment) {
                myDialog = new DataRowFilterDialog((Alignment) theData);
            } else if (theData instanceof Phenotype) {
                myDialog = new DataRowFilterDialog((Phenotype) theData);
            } else {
                JOptionPane.showMessageDialog(getParentFrame(), "Invalid selection. Please select a single sequence or phenotype.");
                return null;
            }
            myDialog.setLocationRelativeTo(getParentFrame());
            myDialog.setVisible(true);
            if (!myDialog.isSelectionMade()) {
                return null;
            }
            myIdsToKeep = myDialog.getIdsToKeep();
            myDialog.dispose();
        }

        if (((myIdsToKeep == null) || (myIdsToKeep.getIdCount() == 0))
                && ((myIdsToRemove == null) || (myIdsToRemove.getIdCount() == 0))) {
            return null;
        }

        Object result = null;
        int count = 0;
        if (theData instanceof Alignment) {
            if (myIdsToKeep != null) {
                result = FilterAlignment.getInstance((Alignment) theData, myIdsToKeep);
            } else if (myIdsToRemove != null) {
                result = FilterAlignment.getInstanceRemoveIDs((Alignment) theData, myIdsToRemove);
            }
            count = ((Alignment) result).getIdGroup().getIdCount();
        } else if (theData instanceof Phenotype) {
            if (myIdsToKeep != null) {
                result = FilterPhenotype.getInstance((Phenotype) theData, myIdsToKeep, null);
            } else if (myIdsToRemove != null) {
                result = FilterPhenotype.getInstanceRemoveIDs((Phenotype) theData, myIdsToRemove);
            }
            count = ((FilterPhenotype) result).getRowCount();
        } else {
            myLogger.error("processDatum: Please input a single sequence or phenotype.  Unknown data type: " + theData.getClass().getName());
            return null;
        }

        String theName, theComment;
        theName = inDatum.getName() + "_" + count + " Rows";
        theComment = "Subset of " + count + " Taxa\n" + inDatum.getComment();
        return new Datum(theName, result, theComment);

    }

    public IdGroup getIdsToKeep() {
        return myIdsToKeep;
    }

    public void setIdsToKeep(IdGroup idsToKeep) {
        myIdsToKeep = idsToKeep;
    }

    public void setIdsToRemove(IdGroup idsToRemove) {
        myIdsToRemove = idsToRemove;
    }

    /**
     * Icon for this plugin to be used in buttons, etc.
     *
     * @return ImageIcon
     */
    public ImageIcon getIcon() {
        URL imageURL = FilterTaxaAlignmentPlugin.class.getResource("images/Filter_horizontal.gif");
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
        return "Taxa";
    }

    /**
     * Tool Tip Text for this plugin
     *
     * @return String
     */
    public String getToolTipText() {
        return "Select Taxa Within Dataset";
    }
}

/**
 * User: dallas & Ed 12/27/2006
 * Date: May 11, 2004
 * Time: 4:16:54 PM
 */
class DataRowFilterDialog extends JDialog {

    private IdGroup theIdsToKeep = null;
    private boolean isSelectionMade = false;
    private JPanel primaryPanel = new JPanel();
    private JScrollPane jScrollPane1;
    private JPanel buttonPanel = new JPanel();
    private JButton captureSelectionButton = new JButton();
    private JButton deleteSelectedButton = new JButton();
    private JButton closeButton = new JButton();
    private JTable table;

    /**
     * For filtering sequences, genotypes, and traits
     * @param etr an extended Table report with taxa in the first column (eg. Alignments)
     */
    public DataRowFilterDialog(Alignment etr) {
        super((Frame) null, true);
        Object[][] tableData = new Object[etr.getSequenceCount()][2];
        for (int i = 0; i < tableData.length; i++) {
            tableData[i][0] = etr.getIdGroup().getIdentifier(i);
            tableData[i][1] = etr.getBaseAsStringRow(i);
        }
        Object[] tableColumnNames = {"Taxa", "Sequence"};
        table = new JTable(tableData, tableColumnNames) {

            public boolean isCellEditable(int row, int column) {
                return false;
            }

            public boolean getRowSelectionAllowed() {
                return true;
            }

            public boolean getColumnSelectionAllowed() {
                return false;
            }
        };
        try {
            initDialog();
            pack();
        } catch (Exception ex) {
            ex.printStackTrace();
        }

    }

    public DataRowFilterDialog(Phenotype etr) {
        super((Frame) null, true);
        table = new JTable(etr.getTableData(), etr.getTableColumnNames()) {

            public boolean isCellEditable(int row, int column) {
                return false;
            }

            public boolean getRowSelectionAllowed() {
                return true;
            }

            public boolean getColumnSelectionAllowed() {
                return false;
            }
        };
        try {
            initDialog();
            pack();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    private void initDialog() throws Exception {
        primaryPanel.setLayout(new BorderLayout());
        jScrollPane1 = new JScrollPane(table);
        jScrollPane1.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
        jScrollPane1.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);

        captureSelectionButton.setText("Capture Selected");
        captureSelectionButton.setMnemonic(KeyEvent.VK_S);
        captureSelectionButton.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                captureSelectionButton_actionPerformed(e);
            }
        });
        deleteSelectedButton.setText("Capture Unselected");
        deleteSelectedButton.setMnemonic(KeyEvent.VK_D);
        deleteSelectedButton.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                deleteSelectedButton_actionPerformed(e);
            }
        });
        closeButton.setText("Close");
        closeButton.setMnemonic(KeyEvent.VK_C);
        closeButton.addActionListener(new java.awt.event.ActionListener() {

            public void actionPerformed(ActionEvent e) {
                closeButton_actionPerformed(e);
            }
        });

        primaryPanel.add(jScrollPane1, BorderLayout.CENTER);
        primaryPanel.add(buttonPanel, BorderLayout.SOUTH);
        buttonPanel.add(captureSelectionButton, null);
        buttonPanel.add(deleteSelectedButton, null);
        buttonPanel.add(closeButton, null);
        getContentPane().add(primaryPanel);
    }

    private void captureSelectionButton_actionPerformed(ActionEvent e) {
        makeIdsList(true);
        isSelectionMade = true;
        this.setVisible(false);
    }

    private void deleteSelectedButton_actionPerformed(ActionEvent e) {
        makeIdsList(false);
        isSelectionMade = true;
        this.setVisible(false);
    }

    private void makeIdsList(boolean isKeepSelected) {
        ArrayList<Identifier> ids = new ArrayList<Identifier>();
        for (int i = 0; i < table.getRowCount(); i++) {
            if (table.isRowSelected(i) == isKeepSelected) {
                ids.add((Identifier) table.getModel().getValueAt(i, 0));
            }
        }
        Identifier[] idsp = new Identifier[ids.size()];
        ids.toArray(idsp);
        theIdsToKeep = new SimpleIdGroup(idsp);
    }

    private void closeButton_actionPerformed(ActionEvent e) {
        this.setVisible(false);
        isSelectionMade = false;
    }

    public IdGroup getIdsToKeep() {
        return theIdsToKeep;
    }

    public boolean isSelectionMade() {
        return isSelectionMade;
    }
}
