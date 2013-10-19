/*
 * FilterTaxaAlignmentPlugin.java
 *
 * Created on December 22, 2006, 5:02 PM
 *
 */
package net.maizegenetics.baseplugins;

import net.maizegenetics.gui.AbstractAvailableListModel;
import net.maizegenetics.gui.SelectFromAvailableDialog;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.FilterAlignment;
import net.maizegenetics.pal.alignment.FilterPhenotype;
import net.maizegenetics.pal.alignment.Phenotype;
import net.maizegenetics.pal.taxa.TaxaList;
import net.maizegenetics.pal.taxa.TaxaListBuilder;
import net.maizegenetics.pal.taxa.Taxon;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginEvent;
import org.apache.log4j.Logger;

import javax.swing.*;
import java.awt.*;
import java.net.URL;
import java.util.List;

/**
 *
 * @author Ed Buckler
 */
public class FilterTaxaAlignmentPlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(FilterTaxaAlignmentPlugin.class);
    private TaxaList myIdsToKeep = null;
    private TaxaList myIdsToRemove = null;

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
            TaxaList origIdGroup = null;
            SelectFromAvailableDialog dialog = null;
            if (theData instanceof Alignment) {
                final Alignment alignment = (Alignment) theData;
                origIdGroup = alignment.getTaxaList();
                AbstractAvailableListModel listModel = new AbstractAvailableListModel() {

                    @Override
                    public int getRealSize() {
                        return alignment.getSequenceCount();
                    }

                    @Override
                    public String getRealElementAt(int index) {
                        return alignment.getTaxaName(index);
                    }
                };
                dialog = new SelectFromAvailableDialog(getParentFrame(), "Taxa Filter", listModel);
            } else if (theData instanceof Phenotype) {
                final Phenotype phenotype = (Phenotype) theData;
                origIdGroup = phenotype.getTaxa();
                AbstractAvailableListModel listModel = new AbstractAvailableListModel() {

                    @Override
                    public int getRealSize() {
                        return phenotype.getNumberOfTaxa();
                    }

                    @Override
                    public String getRealElementAt(int index) {
                        return phenotype.getTaxon(index).getName();
                    }
                };
                dialog = new SelectFromAvailableDialog(getParentFrame(), "Taxa Filter", listModel);
            } else {
                JOptionPane.showMessageDialog(getParentFrame(), "Invalid selection. Please select a single sequence or phenotype.");
                return null;
            }
            dialog.setLocationRelativeTo(getParentFrame());
            dialog.setVisible(true);
            if (dialog.isCanceled()) {
                return null;
            }
            int[] indicesToKeep = dialog.getDesiredIndices();
            Taxon[] ids = new Taxon[indicesToKeep.length];
            for (int i = 0; i < indicesToKeep.length; i++) {
                ids[i] = origIdGroup.get(indicesToKeep[i]);
            }
            myIdsToKeep=new TaxaListBuilder().addAll(ids).build();
            dialog.dispose();
        }

        if (((myIdsToKeep == null) || (myIdsToKeep.getTaxaCount() == 0))
                && ((myIdsToRemove == null) || (myIdsToRemove.getTaxaCount() == 0))) {
            return null;
        }

        Object result = null;
        int count = 0;
        if (theData instanceof Alignment) {
            if (myIdsToKeep != null) {
                result = FilterAlignment.getInstance((Alignment) theData, myIdsToKeep, false);
            } else if (myIdsToRemove != null) {
                result = FilterAlignment.getInstanceRemoveIDs((Alignment) theData, myIdsToRemove);
            }
            count = ((Alignment) result).getTaxaCount();
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

    public TaxaList getIdsToKeep() {
        return myIdsToKeep;
    }

    public void setIdsToKeep(TaxaList idsToKeep) {
        myIdsToKeep = idsToKeep;
    }

    public void setIdsToRemove(TaxaList idsToRemove) {
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
        return "Taxa Names";
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
