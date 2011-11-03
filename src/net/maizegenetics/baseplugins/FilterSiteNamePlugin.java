/*
 * FilterSiteNamePlugin.java
 *
 * Created on November 2, 2011
 *
 */
package net.maizegenetics.baseplugins;

import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.FilterAlignment;
import net.maizegenetics.pal.ids.Identifier;

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

import javax.swing.table.AbstractTableModel;
import javax.swing.table.TableModel;

import org.apache.log4j.Logger;

/**
 *
 * @author terry
 */
public class FilterSiteNamePlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(FilterSiteNamePlugin.class);
    private int[] mySitesToKeep = null;

    /** Creates a new instance of FilterSiteNamePlugin */
    public FilterSiteNamePlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    public DataSet performFunction(DataSet input) {

        try {

            List inputData = input.getDataOfType(Alignment.class);
            if (inputData.size() != 1) {
                if (isInteractive()) {
                    JOptionPane.showMessageDialog(getParentFrame(), "Invalid selection. Please select a single alignment.");
                } else {
                    myLogger.error("performFunction: Please input a single alignment.");
                }
                return null;
            }

            Datum td = processDatum((Datum) inputData.get(0), isInteractive());
            if (td == null) {
                return null;
            }

            DataSet output = new DataSet(td, this);

            fireDataSetReturned(new PluginEvent(output, FilterSiteNamePlugin.class));

            return output;

        } finally {
            fireProgress(100);
        }
    }

    private Datum processDatum(Datum inDatum, boolean isInteractive) {

        Alignment alignment = (Alignment) inDatum.getData();

        if (isInteractive) {
            SiteNameFilterDialog myDialog = new SiteNameFilterDialog(alignment);
            myDialog.setLocationRelativeTo(getParentFrame());
            myDialog.setVisible(true);
            if (!myDialog.isSelectionMade()) {
                return null;
            }
            mySitesToKeep = myDialog.getSitesToKeep();
            myDialog.dispose();
        }

        if (((mySitesToKeep == null) || (mySitesToKeep.length == 0))) {
            return null;
        }

        Alignment result = FilterAlignment.getInstance(alignment, mySitesToKeep);

        String theName, theComment;
        theName = inDatum.getName() + "_" + result.getSiteCount() + "_Sites";
        theComment = "Subset of " + result.getSiteCount() + " from " + alignment.getSiteCount() + " Sites\n" + inDatum.getComment();
        return new Datum(theName, result, theComment);

    }

    public int[] getSitesToKeep() {
        return mySitesToKeep;
    }

    public void setIdsToKeep(int[] sitesToKeep) {
        mySitesToKeep = sitesToKeep;
    }

    /**
     * Icon for this plugin to be used in buttons, etc.
     *
     * @return ImageIcon
     */
    public ImageIcon getIcon() {
        URL imageURL = FilterSiteNamePlugin.class.getResource("images/Filter.gif");
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
        return "Site Names";
    }

    /**
     * Tool Tip Text for this plugin
     *
     * @return String
     */
    public String getToolTipText() {
        return "Select Site Names Within Dataset";
    }
}

/**
 */
class SiteNameFilterDialog extends JDialog {

    private int[] theSitesToKeep = null;
    private boolean isSelectionMade = false;
    private JPanel primaryPanel = new JPanel();
    private JScrollPane jScrollPane1;
    private JPanel buttonPanel = new JPanel();
    private JButton captureSelectionButton = new JButton();
    private JButton deleteSelectedButton = new JButton();
    private JButton closeButton = new JButton();
    private JTable table;

    /**
     * For filtering alignment by site names.
     *
     * @param alignment alignment
     */
    public SiteNameFilterDialog(final Alignment alignment) {
        super((Frame) null, true);
        TableModel tableModel = new AbstractTableModel() {

            @Override
            public int getRowCount() {
                return alignment.getSiteCount();
            }

            @Override
            public int getColumnCount() {
                return 1;
            }

            @Override
            public Object getValueAt(int rowIndex, int columnIndex) {
                return alignment.getSNPID(rowIndex);
            }

            @Override
            public String getColumnName(int column) {
                if (column == 0) {
                    return "Site Name";
                } else {
                    return super.getColumnName(column);
                }
            }
        };
        table = new JTable(tableModel) {

            @Override
            public boolean isCellEditable(int row, int column) {
                return false;
            }

            @Override
            public boolean getRowSelectionAllowed() {
                return true;
            }

            @Override
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
        if (isKeepSelected) {
            theSitesToKeep = new int[table.getSelectedRowCount()];
        } else {
            theSitesToKeep = new int[table.getRowCount() - table.getSelectedRowCount()];
        }
        int count = 0;
        for (int i = 0; i < table.getRowCount(); i++) {
            if (table.isRowSelected(i) == isKeepSelected) {
                theSitesToKeep[count++] = i;
            }
        }
    }

    private void closeButton_actionPerformed(ActionEvent e) {
        this.setVisible(false);
        isSelectionMade = false;
    }

    public int[] getSitesToKeep() {
        return theSitesToKeep;
    }

    public boolean isSelectionMade() {
        return isSelectionMade;
    }
}
