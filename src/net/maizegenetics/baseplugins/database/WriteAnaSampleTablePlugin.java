/*
 * WriteAnaSampleTablePlugin.java
 *
 * Created on June 7, 2007
 *
 */

package net.maizegenetics.baseplugins.database;


import gov.usda.gdpc.GenotypeExperiment;
import gov.usda.gdpc.database.DefaultJDBCConnection;
import gov.usda.gdpc.database.JDBCConnection;

import java.sql.SQLException;

import net.maizegenetics.pal.report.TableReport;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;

import javax.swing.*;

import java.util.Iterator;
import java.util.List;


/**
 *
 * @author terryc
 */
public class WriteAnaSampleTablePlugin extends AbstractPlugin {
    
    private static final String COLUMN_SITE_TYPE = "Site_Type";
    private static final String COLUMN_START_SITE = "StartSite";
    private static final String COLUMN_END_SITE = "EndSite";
    private static final String COLUMN_SITE = "Site";
    private static final String COLUMN_SITE_COUNT = "SiteCount";
    private static final String COLUMN_SEQ_SITES = "SegSites";
    private static final String COLUMN_PI = "Pi";
    private static final String COLUMN_THETA = "Theta";
    private static final String COLUMN_HAPLOTYPES = "Haplotypes";
    
    
    private static final String DIV_ALLELE_ASSAY_ID = "div_allele_assay_id";
    
    private static final String ANA_SAMPLE = "ana_sample";
    private static final String ANA_SAMPLE_ID = "ana_sample_id";
    private static final String ANA_SAMPLE_SAMPLE = "sample";
    private static final String ANA_SAMPLE_SITE = "site";
    private static final String ANA_SAMPLE_START_SITE = "start_site";
    private static final String ANA_SAMPLE_END_SITE = "end_site";
    private static final String ANA_SAMPLE_DIV_PI = "div_pi";
    private static final String ANA_SAMPLE_DIV_THETA = "div_theta";
    private static final String ANA_SAMPLE_LD_PNG = "ld_png";
    private static final String ANA_SAMPLE_TREE_PNG = "tree_png";
    
    public static final String LD_PNG_FILENAME_TAG = "LD_PNG";
    public static final String TREE_PNG_FILENAME_TAG = "TREE_PNG";
    public static final String GENOTYPE_EXP_TAG = "GENOTYPE_EXP";
    
    private final JDBCConnection myConnection;
    
    
    /**
     * Creates a new instance of WriteAnaSampleTablePlugin
     */
    public WriteAnaSampleTablePlugin(String driver, String url, String uid, String passwd) {
        super(null, false);
        myConnection = new DefaultJDBCConnection("NO_SOURCE", null, driver, url, uid, passwd);
        //myConnection = null;
    }
    
    
    public DataSet performFunction(DataSet input) {
        
        List <Datum> reports = input.getDataOfType(TableReport.class);
        if (reports.size() != 1) {
            throw new IllegalArgumentException("WriteAnaSampleTablePlugin: performFunction: input must have one TableReport.");
        }
        TableReport report = (TableReport) reports.get(0).getData();
        
        List <Datum> exps = input.getDataOfType(GenotypeExperiment.class);
        if (exps.size() != 1) {
            throw new IllegalArgumentException("WriteAnaSampleTablePlugin: performFunction: input must have one GenotypeExperiment.");
        }
        GenotypeExperiment exp = (GenotypeExperiment) exps.get(0).getData();
        
        List <Datum> filenames = input.getDataOfType(String.class);
        if (filenames.size() > 2) {
            throw new IllegalArgumentException("WriteAnaSampleTablePlugin: performFunction: input should have two or less filenames.");
        }
        String ldPngFilename = null;
        String treePngFilename = null;
        Iterator itr = filenames.iterator();
        while (itr.hasNext()) {
            Datum current = (Datum) itr.next();
            String name = current.getName();
            if (name.equals(LD_PNG_FILENAME_TAG)) {
                ldPngFilename = (String)current.getData();
            } else if (name.equals(TREE_PNG_FILENAME_TAG)) {
                treePngFilename = (String)current.getData();
            }
        }
        
        insertTableReport(report, exp, ldPngFilename, treePngFilename);

        return null;

    }
    
    
    private void insertTableReport(TableReport tr, GenotypeExperiment exp, String ldPng, String treePng) {
        
        Object [][] data = tr.getTableData();
        Object [] colNames = tr.getTableColumnNames();
        
        int siteIndex = findIndex(COLUMN_SITE, colNames);
        int startSiteIndex = findIndex(COLUMN_START_SITE, colNames);
        int endSiteIndex = findIndex(COLUMN_END_SITE, colNames);
        int divPIIndex = findIndex(COLUMN_PI, colNames);
        int divThetaIndex = findIndex(COLUMN_THETA, colNames);
        
        
        int rows = data.length;
        for(int i=0; i<rows; i++) {
            
            Object site = data[i][siteIndex];
            Object startSite = data[i][startSiteIndex];
            Object endSite = data[i][endSiteIndex];
            Object divPI = data[i][divPIIndex];
            Object divTheta = data[i][divThetaIndex];
            
            StringBuilder buffer = new StringBuilder();
            
            // insert into ana_sample (div_allele_assay_id, sample, site, start_site, end_site,
            // div_pi, div_theta, ld_png, tree_png) values (1234, 'maize', '344.44', ...)
            buffer.append("insert into ");
            buffer.append(ANA_SAMPLE);
            buffer.append(" (");
            buffer.append(DIV_ALLELE_ASSAY_ID);
            buffer.append(", ");
            buffer.append(ANA_SAMPLE_SAMPLE);
            buffer.append(", ");
            buffer.append(ANA_SAMPLE_SITE);
            buffer.append(", ");
            buffer.append(ANA_SAMPLE_START_SITE);
            buffer.append(", ");
            buffer.append(ANA_SAMPLE_END_SITE);
            buffer.append(", ");
            buffer.append(ANA_SAMPLE_DIV_PI);
            buffer.append(", ");
            buffer.append(ANA_SAMPLE_DIV_THETA);
            buffer.append(", ");
            buffer.append(ANA_SAMPLE_LD_PNG);
            buffer.append(", ");
            buffer.append(ANA_SAMPLE_TREE_PNG);
            buffer.append(") values (");
            buffer.append(exp.getID().getNumber());
            buffer.append(", ");
            buffer.append(" '");
            buffer.append("Maize");
            buffer.append("', '");
            buffer.append(site);
            buffer.append("', '");
            buffer.append(startSite);
            buffer.append("', '");
            buffer.append(endSite);
            buffer.append("', '");
            buffer.append(divPI);
            buffer.append("', '");
            buffer.append(divTheta);
            buffer.append("', '");
            buffer.append(ldPng);
            buffer.append("', '");
            buffer.append(treePng);
            buffer.append("')");
            
            String statement = buffer.toString();
            
            System.out.println("statement: " + statement);
            
            try {
                myConnection.executeUpdate(statement);
            } catch (SQLException e) {
                e.printStackTrace();
            }
            
        }
        
    }
    
    
    private int findIndex(String target, Object [] headings) {
        
        for (int i=0; i<headings.length; i++) {
            if (target.equals(headings[i].toString())) {
                return i;
            }
        }
        
        return -1;
        
    }
    
    
    private void printTableReport(TableReport tr) {
        
        Object [][] data = tr.getTableData();
        Object [] colNames = tr.getTableColumnNames();
        
        int cols = colNames.length;
        int rows = data.length;
        for(int j=0; j<cols; j++) {
            System.out.print(colNames[j]);
            if(j<(cols-1)) {
                System.out.print("\t");
            }
        }
        
        System.out.println();
        
        for(int i=0; i<rows; i++) {
            for(int j=0; j<cols; j++) {
                System.out.print(data[i][j]);
                if(j<(cols-1)) {
                    System.out.print("\t");
                }
            }
            System.out.println("\n");
        }
        
    }
    
    
    /**
     * Icon for this plugin to be used in buttons, etc.
     *
     * @return ImageIcon
     */
    public ImageIcon getIcon() {
        return null;
    }
    
    /**
     * Button name for this plugin to be used in buttons, etc.
     *
     * @return String
     */
    public String getButtonName() {
        return null;
    }
    
    /**
     * Tool Tip Text for this plugin
     *
     * @return String
     */
    public String getToolTipText() {
        return null;
    }
    
}
