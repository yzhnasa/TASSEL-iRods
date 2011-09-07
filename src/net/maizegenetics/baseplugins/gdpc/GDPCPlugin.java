/*
 * GDPCPlugin.java
 *
 */
package net.maizegenetics.baseplugins.gdpc;

import gov.usda.gdpc.*;

import gov.usda.gdpc.browser.Browser;
import gov.usda.gdpc.browser.BrowserSettings;
import gov.usda.gdpc.gui.GDPCListener;

import java.awt.Component;
import java.awt.Frame;

import java.net.URL;

import java.util.List;

import javax.swing.ImageIcon;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JPanel;

import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;

/**
 *
 * @author terryc
 */
public class GDPCPlugin extends AbstractPlugin implements DBGateway {

    final private DBGateway myDBGateway;
    final private Browser myBrowser;

    /** Creates a new instance of GDPCPlugin */
    public GDPCPlugin() {
        this(null, false);
    }

    public GDPCPlugin(Frame parentFrame, boolean isInteractive) {

        super(parentFrame, isInteractive);

        if (isInteractive) {
            Component glassPane = null;
            if ((parentFrame != null) && (parentFrame instanceof JFrame)) {
                glassPane = ((JFrame) parentFrame).getGlassPane();
            }
            myBrowser = Browser.createInstanceForTASSEL(glassPane);

            myDBGateway = myBrowser.getGateway();
        } else {
            myBrowser = null;
            myDBGateway = DefaultDBGateway.getInstance();
        }

    }

    public DataSet performFunction(DataSet input) {

        List data = input.getDataOfType(Filter.class);
        if ((data != null) && (data.size() != 0)) {
            Filter filter = (Filter) ((Datum) data.get(0)).getData();
            return performFunctionGivenFilter(filter);
        }

        data = input.getDataOfType(TaxonGroup.class);
        TaxonGroup taxonGroup = null;
        if ((data != null) && (data.size() != 0)) {
            taxonGroup = (TaxonGroup) ((Datum) data.get(0)).getData();
        }

        data = input.getDataOfType(GenotypeExperimentGroup.class);
        GenotypeExperimentGroup genotypeExpGroup = null;
        if ((data != null) && (data.size() != 0)) {
            genotypeExpGroup = (GenotypeExperimentGroup) ((Datum) data.get(0)).getData();
        }

        if ((taxonGroup != null) && (genotypeExpGroup != null)) {
            return getGenotypeTableDataSet(genotypeExpGroup, taxonGroup);
        }


        data = input.getDataOfType(EnvironmentExperimentGroup.class);
        EnvironmentExperimentGroup envExpGroup = null;
        if ((data != null) && (data.size() != 0)) {
            envExpGroup = (EnvironmentExperimentGroup) ((Datum) data.get(0)).getData();
        }

        data = input.getDataOfType(PhenotypeOntologyGroup.class);
        PhenotypeOntologyGroup ontologyGroup = null;
        if ((data != null) && (data.size() != 0)) {
            ontologyGroup = (PhenotypeOntologyGroup) ((Datum) data.get(0)).getData();
        }

        if ((taxonGroup != null) && (envExpGroup != null) && (ontologyGroup != null)) {
            return getPhenotypeGroupDataSet(envExpGroup, taxonGroup, ontologyGroup);
        }


        data = input.getDataOfType(Property[].class);
        if ((data != null) && (data.size() != 0)) {
            Property[] properties = (Property[]) ((Datum) data.get(0)).getData();
            return getDistinctPropertiesDataSet(properties);
        }


        throw new IllegalArgumentException("GDPCPlugin: performFunction: input does contain objects to perform any function.");

    }

    public DataSet performFunctionGivenFilter(Filter filter) {

        if (filter instanceof TaxonFilter) {
            return getTaxonGroupDataSet((TaxonFilter) filter);
        } else if (filter instanceof StudyFilter) {
            return getStudyGroupDataSet((StudyFilter) filter);
        } else if (filter instanceof LocusFilter) {
            return getLocusGroupDataSet((LocusFilter) filter);
        } else if (filter instanceof GenotypeExperimentFilter) {
            return getGenotypeExperimentGroupDataSet((GenotypeExperimentFilter) filter);
        } else if (filter instanceof EnvironmentExperimentFilter) {
            return getEnvironmentExperimentGroupDataSet((EnvironmentExperimentFilter) filter);
        } else if (filter instanceof TaxonParentFilter) {
            return getTaxonParentGroupDataSet((TaxonParentFilter) filter);
        } else if (filter instanceof PhenotypeOntologyFilter) {
            return getPhenotypeOntologyGroupDataSet((PhenotypeOntologyFilter) filter);
        } else if (filter instanceof LocalityFilter) {
            return getLocalityGroupDataSet((LocalityFilter) filter);
        } else {
            throw new IllegalArgumentException("GDPCPlugin: performFunctionGivenFilter: unknown filter type: " + filter.getClass().getName());
        }

    }

    public String getToolTipText() {
        return "Open Genomic Diversity and Phenotype Connection Data Browser";
    }

    public ImageIcon getIcon() {

        URL imageURL = GDPCPlugin.class.getResource("gdpc.gif");

        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }

    }

    public String getButtonName() {
        return "GDPC";
    }

    public JPanel getPanel() {
        return myBrowser;
    }

    public JMenu getMenu() {
        return myBrowser.getBrowserMenu();
    }

    public DataSet getTaxonGroupDataSet(TaxonFilter filter) {
        TaxonGroup result = myDBGateway.getTaxonGroup(filter);
        DataSet tds = new DataSet(new Datum("TaxonGroup", result, null), this);
        fireDataSetReturned(tds);
        return tds;
    }

    public TaxonGroup getTaxonGroup(TaxonFilter filter) {
        TaxonGroup result = myDBGateway.getTaxonGroup(filter);
        DataSet tds = new DataSet(new Datum("TaxonGroup", result, null), this);
        fireDataSetReturned(tds);
        return result;
    }

    public DataSet getPhenotypeGroupDataSet(EnvironmentExperimentGroup experimentGroup, TaxonGroup taxonGroup, PhenotypeOntologyGroup ontologyGroup) {
        PhenotypeGroup result = myDBGateway.getPhenotypeGroup(experimentGroup, taxonGroup, ontologyGroup);
        DataSet tds = new DataSet(new Datum("PhenotypeGroup", result, null), this);
        fireDataSetReturned(tds);
        return tds;
    }

    public PhenotypeGroup getPhenotypeGroup(EnvironmentExperimentGroup experimentGroup, TaxonGroup taxonGroup, PhenotypeOntologyGroup ontologyGroup) {
        PhenotypeGroup result = myDBGateway.getPhenotypeGroup(experimentGroup, taxonGroup, ontologyGroup);
        fireDataSetReturned(new DataSet(new Datum("PhenotypeGroup", result, null), this));
        return result;
    }

    public DataSet getStudyGroupDataSet(StudyFilter filter) {
        StudyGroup result = myDBGateway.getStudyGroup(filter);
        DataSet tds = new DataSet(new Datum("StudyGroup", result, null), this);
        fireDataSetReturned(tds);
        return tds;
    }

    public StudyGroup getStudyGroup(StudyFilter filter) {
        StudyGroup result = myDBGateway.getStudyGroup(filter);
        fireDataSetReturned(new DataSet(new Datum("StudyGroup", result, null), this));
        return result;
    }

    public DataSet getLocusGroupDataSet(LocusFilter filter) {
        LocusGroup result = myDBGateway.getLocusGroup(filter);
        DataSet tds = new DataSet(new Datum("LocusGroup", result, null), this);
        fireDataSetReturned(tds);
        return tds;
    }

    public LocusGroup getLocusGroup(LocusFilter filter) {
        LocusGroup result = myDBGateway.getLocusGroup(filter);
        fireDataSetReturned(new DataSet(new Datum("LocusGroup", result, null), this));
        return result;
    }

    public DataSet getGenotypeExperimentGroupDataSet(GenotypeExperimentFilter filter) {
        GenotypeExperimentGroup result = myDBGateway.getGenotypeExperimentGroup(filter);
        DataSet tds = new DataSet(new Datum("GenotypeExperimentGroup", result, null), this);
        fireDataSetReturned(tds);
        return tds;
    }

    public GenotypeExperimentGroup getGenotypeExperimentGroup(GenotypeExperimentFilter filter) {
        GenotypeExperimentGroup result = myDBGateway.getGenotypeExperimentGroup(filter);
        fireDataSetReturned(new DataSet(new Datum("GenotypeExperimentGroup", result, null), this));
        return result;
    }

    public DataSet getEnvironmentExperimentGroupDataSet(EnvironmentExperimentFilter filter) {
        EnvironmentExperimentGroup result = myDBGateway.getEnvironmentExperimentGroup(filter);
        DataSet tds = new DataSet(new Datum("EnvironmentExperimentGroup", result, null), this);
        fireDataSetReturned(tds);
        return tds;
    }

    public EnvironmentExperimentGroup getEnvironmentExperimentGroup(EnvironmentExperimentFilter filter) {
        EnvironmentExperimentGroup result = myDBGateway.getEnvironmentExperimentGroup(filter);
        fireDataSetReturned(new DataSet(new Datum("EnvironmentExperimentGroup", result, null), this));
        return result;
    }

    public boolean removeDBConnection(DBConnection connection) {
        return myDBGateway.removeDBConnection(connection);
    }

    public boolean addDBConnection(DBConnection connection) {
        return myDBGateway.addDBConnection(connection);
    }

    public GenotypeGroup getGenotypeGroup(GenotypeExperimentGroup experimentGroup, TaxonGroup taxonGroup) {
        GenotypeGroup result = myDBGateway.getGenotypeGroup(experimentGroup, taxonGroup);
        fireDataSetReturned(new DataSet(new Datum("GenotypeGroup", result, null), this));
        return result;
    }

    public DataSet getGenotypeTableDataSet(GenotypeExperimentGroup experimentGroup, TaxonGroup taxonGroup) {
        GenotypeTable result = myDBGateway.getGenotypeTable(experimentGroup, taxonGroup);
        DataSet tds = new DataSet(new Datum("GenotypeTable", result, null), this);
        fireDataSetReturned(tds);
        return tds;
    }

    public GenotypeTable getGenotypeTable(GenotypeExperimentGroup experimentGroup, TaxonGroup taxonGroup) {
        GenotypeTable result = myDBGateway.getGenotypeTable(experimentGroup, taxonGroup);
        fireDataSetReturned(new DataSet(new Datum("GenotypeTable", result, null), this));
        return result;
    }

    public DataSet getTaxonParentGroupDataSet(TaxonParentFilter filter) {
        TaxonParentGroup result = myDBGateway.getTaxonParentGroup(filter);
        DataSet tds = new DataSet(new Datum("TaxonParentGroup", result, null), this);
        fireDataSetReturned(tds);
        return tds;
    }

    public TaxonParentGroup getTaxonParentGroup(TaxonParentFilter filter) {
        TaxonParentGroup result = myDBGateway.getTaxonParentGroup(filter);
        fireDataSetReturned(new DataSet(new Datum("TaxonParentGroup", result, null), this));
        return result;
    }

    public DataSet getDistinctPropertiesDataSet(Property[] properties) {
        DistinctPropertyValues result = myDBGateway.getDistinctProperties(properties);
        DataSet tds = new DataSet(new Datum("DistinctPropertyValues", result, null), this);
        fireDataSetReturned(tds);
        return tds;
    }

    public DistinctPropertyValues getDistinctProperties(Property[] properties) {
        DistinctPropertyValues result = myDBGateway.getDistinctProperties(properties);
        fireDataSetReturned(new DataSet(new Datum("DistinctPropertyValues", result, null), this));
        return result;
    }

    public DataSet getLocalityGroupDataSet(LocalityFilter filter) {
        LocalityGroup result = myDBGateway.getLocalityGroup(filter);
        DataSet tds = new DataSet(new Datum("LocalityGroup", result, null), this);
        fireDataSetReturned(tds);
        return tds;
    }

    public LocalityGroup getLocalityGroup(LocalityFilter filter) {
        LocalityGroup result = myDBGateway.getLocalityGroup(filter);
        fireDataSetReturned(new DataSet(new Datum("LocalityGroup", result, null), this));
        return result;
    }

    public DataSet getPhenotypeOntologyGroupDataSet(PhenotypeOntologyFilter filter) {
        PhenotypeOntologyGroup result = myDBGateway.getPhenotypeOntologyGroup(filter);
        DataSet tds = new DataSet(new Datum("PheotypeOntologyGroup", result, null), this);
        fireDataSetReturned(tds);
        return tds;
    }

    public PhenotypeOntologyGroup getPhenotypeOntologyGroup(PhenotypeOntologyFilter filter) {
        PhenotypeOntologyGroup result = myDBGateway.getPhenotypeOntologyGroup(filter);
        fireDataSetReturned(new DataSet(new Datum("PheotypeOntologyGroup", result, null), this));
        return result;
    }

    public DBConnection[] getDBConnections() {
        return myDBGateway.getDBConnections();
    }

    public int getNumConnections() {
        return myDBGateway.getNumConnections();
    }

    public void addGDPCListener(GDPCListener listener) {
        myBrowser.addGDPCListener(listener);
    }

    public void saveSettings() {
        BrowserSettings.getInstance(myBrowser).saveSettings();
    }
}
