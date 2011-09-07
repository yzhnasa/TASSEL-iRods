/*
 * ConvertGDPCToPALPlugin.java
 *
 * Created on January 30, 2007, 2:38 PM
 *
 */
package net.maizegenetics.baseplugins.gdpc;

import gov.usda.gdpc.DefaultGenotypeTable;
import gov.usda.gdpc.GenotypeExperiment;
import gov.usda.gdpc.GenotypeExperimentProperty;
import gov.usda.gdpc.GenotypeGroup;
import gov.usda.gdpc.GenotypeTable;
import gov.usda.gdpc.PhenotypeTable;
import gov.usda.gdpc.TaxonParentGroup;
import gov.usda.gdpc.gui.GDPCEvent;
import gov.usda.gdpc.gui.GDPCListener;

import java.util.ArrayList;
import java.util.List;

import javax.swing.ImageIcon;

import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.SimpleAlignment;
import net.maizegenetics.pal.alignment.SimplePhenotype;
import net.maizegenetics.pal.distance.DistanceMatrix;

import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;

/**
 *
 * @author terryc
 */
public class ConvertGDPCToPALPlugin extends AbstractPlugin implements GDPCListener {

    /**
     * Creates a new instance of ConvertGDPCToPALPlugin
     */
    public ConvertGDPCToPALPlugin() {
        super(null, false);
    }

    public DataSet performFunction(DataSet input) {

        List result = new ArrayList();

        try {

            List data = input.getDataOfType(GenotypeTable.class);
            if ((data != null) && (data.size() != 0)) {
                GDPCEvent event = new GDPCEvent(this);
                GenotypeTable table = (GenotypeTable) ((Datum) data.get(0)).getData();
                event.addData(table);
                DataSet tds = fireGenotypeTable(event);
                if (tds != null) {
                    result.add(tds);
                }
            }

            data = input.getDataOfType(GenotypeGroup.class);
            if ((data != null) && (data.size() != 0)) {
                GDPCEvent event = new GDPCEvent(this);
                GenotypeGroup group = (GenotypeGroup) ((Datum) data.get(0)).getData();
                GenotypeTable table = DefaultGenotypeTable.getInstance(group);
                event.addData(table);
                DataSet tds = fireGenotypeTable(event);
                if (tds != null) {
                    result.add(tds);
                }
            }

            data = input.getDataOfType(PhenotypeTable.class);
            if ((data != null) && (data.size() != 0)) {
                GDPCEvent event = new GDPCEvent(this);
                PhenotypeTable table = (PhenotypeTable) ((Datum) data.get(0)).getData();
                event.addData(table);
                DataSet tds = firePhenotypeTable(event);
                if (tds != null) {
                    result.add(tds);
                }
            }

            data = input.getDataOfType(TaxonParentGroup.class);
            if ((data != null) && (data.size() != 0)) {
                GDPCEvent event = new GDPCEvent(this);
                TaxonParentGroup group = (TaxonParentGroup) ((Datum) data.get(0)).getData();
                event.addData(group);
                DataSet tds = fireTaxonParentGroup(event);
                if (tds != null) {
                    result.add(tds);
                }
            }

        } catch (Exception ex) {
            ex.printStackTrace();
        }

        return DataSet.getDataSet(result, this);

    }

    public String getToolTipText() {
        return "";
    }

    public ImageIcon getIcon() {
        return null;
    }

    public String getButtonName() {
        return "";
    }

    public void taxonParentWorkingGroupChanged(GDPCEvent event) {
    }

    public void phenotypeTableChanged(GDPCEvent event) {
    }

    public DataSet getTaxonParentGroup(GDPCEvent event) {
        TaxonParentGroup group = (TaxonParentGroup) event.getData();
        DistanceMatrix matrix = PALUtil.createKinshipMatrix(group);
        Datum data = new Datum("KinshipMatrix", matrix, null);
        return new DataSet(data, this);
    }

    public DataSet fireTaxonParentGroup(GDPCEvent event) {
        DataSet result = getTaxonParentGroup(event);
        if (result != null) {
            fireDataSetReturned(result);
        }
        return result;
    }

    public void loadTaxonParentGroup(GDPCEvent event) {
        fireTaxonParentGroup(event);
    }

    public DataSet getPhenotypeTable(GDPCEvent event) {

        PhenotypeTable table = (PhenotypeTable) event.getData();

        if (table != null) {
            SimplePhenotype theCharacterAlignment = (SimplePhenotype) PALUtil.createCharacterAlignment(table);
            String str = table.numColumns() + " traits/environ";
            Datum data = new Datum(str, theCharacterAlignment, null);
            return new DataSet(data, this);
        }

        return null;

    }

    public DataSet firePhenotypeTable(GDPCEvent event) {
        DataSet result = getPhenotypeTable(event);
        if (result != null) {
            fireDataSetReturned(result);
        }
        return result;
    }

    public void loadPhenotypeTable(GDPCEvent event) {
        firePhenotypeTable(event);
    }

    public DataSet getGenotypeTable(GDPCEvent event) {

        GenotypeTable genoTable = (GenotypeTable) event.getData();

        if (genoTable == null) {
            return null;
        }

        GenotypeExperiment ge = (GenotypeExperiment) genoTable.getColumnHeading(0);
        String polyType = (String) ge.getProperty(GenotypeExperimentProperty.POLY_TYPE);
        if (polyType.equals(GenotypeExperimentProperty.POLY_TYPE_SEQUENCE)) {
            int numberOfGenes = genoTable.getColumnHeadings().size();
            for (int i = 0; i < numberOfGenes; i++) {
                ge = (GenotypeExperiment) genoTable.getColumnHeading(i);
                Alignment aa = (SimpleAlignment) PALUtil.createAnnotationAlignment(genoTable, i);
                Datum data = new Datum("Raw", aa, null);
                return new DataSet(data, this);
            }
        } else {

            Alignment[] theAnnotationAlignment2 = PALUtil.createMultiLocusAnnotatedAlignment(genoTable, null);
            Datum data = new Datum("GDPC", theAnnotationAlignment2[0], null);
            return new DataSet(data, this);
            
        }

        return null;

    }

    public DataSet fireGenotypeTable(GDPCEvent event) {
        DataSet result = getGenotypeTable(event);
        if (result != null) {
            fireDataSetReturned(result);
        }
        return result;
    }

    public void loadGenotypeTable(GDPCEvent event) {
        fireGenotypeTable(event);
    }

    public void loadPressed(GDPCEvent event) {
        // do nothing
    }
}
