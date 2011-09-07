/*
 * WriteAnaTraitAssocTablePlugin.java
 *
 * Created on June 11, 2007
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

import java.util.List;


/**
 *
 * @author terryc
 */
public class WriteAnaTraitAssocTablePlugin extends AbstractPlugin {
    
    private static final String COLUMN_TRAIT = "Trait";
    private static final String COLUMN_LOCUS = "Locus";
    private static final String COLUMN_SITE = "Site";
    private static final String COLUMN_DF = "df";
    private static final String COLUMN_F = "F";
    private static final String COLUMN_P = "p";
    private static final String COLUMN_MODEL_DF = "df_Model";
    private static final String COLUMN_ERROR_DF = "df_Error";
    private static final String COLUMN_ERROR_MS = "MS_Error";
    private static final String COLUMN_RSQ_MODEL = "Rsq_model";
    private static final String COLUMN_RSQ_MARKER = "Rsq_marker";
    
    
    private static final String DIV_ALLELE_ASSAY_ID = "div_allele_assay_id";
    
    private static final String ANA_TRAIT_ASSOC = "ana_trait_assoc";
    private static final String ANA_TRAIT_ASSOC_ID = "ana_trait_assoc_id";
    private static final String ANA_TRAIT_ASSOC_SITE = "site";
    private static final String ANA_TRAIT_ASSOC_METHOD = "method";
    private static final String ANA_TRAIT_ASSOC_DF = "df";
    private static final String ANA_TRAIT_ASSOC_F_STAT = "f_stat";
    private static final String ANA_TRAIT_ASSOC_P_VALUE = "p_value";
    private static final String ANA_TRAIT_ASSOC_MODEL_DF = "model_df";
    private static final String ANA_TRAIT_ASSOC_ERROR_MS = "error_ms";
    private static final String ANA_TRAIT_ASSOC_ERROR_DF = "error_df";
    private static final String ANA_TRAIT_ASSOC_RSQ_MODEL = "rsq_model";
    private static final String ANA_TRAIT_ASSOC_RSQ_MARKER = "rsq_marker";
    
    
    private static final String NULL = "NULL";
    
    
    public static final String METHOD_TAG = "Method";
    
    
    private final JDBCConnection myConnection;
    
    
    /**
     * Creates a new instance of WriteAnaTraitAssocTablePlugin
     */
    public WriteAnaTraitAssocTablePlugin(String driver, String url, String uid, String passwd) {
        super(null, false);
        myConnection = new DefaultJDBCConnection("NO_SOURCE", null, driver, url, uid, passwd);
    }
    
    
    public DataSet performFunction(DataSet input) {
        
        List <Datum> reports = input.getDataOfType(TableReport.class);
        if (reports.size() == 0) {
            throw new IllegalArgumentException("WriteAnaTraitAssocTablePlugin: performFunction: input must have one TableReport.");
        }
        TableReport report = (TableReport) reports.get(0).getData();
        
        List <Datum> exps = input.getDataOfType(GenotypeExperiment.class);
        if (exps.size() != 1) {
            throw new IllegalArgumentException("WriteAnaTraitAssocTablePlugin: performFunction: input must have one GenotypeExperiment.");
        }
        GenotypeExperiment exp = (GenotypeExperiment) exps.get(0).getData();
        
        List <Datum> methods = input.getDataOfType(String.class);
        if (methods.size() != 1) {
            throw new IllegalArgumentException("WriteAnaTraitAssocTablePlugin: performFunction: input must have one Method string.");
        }
        String method = (String) methods.get(0).getData();
        
        
        insertTableReport(report, exp, method);

        return null;
        
    }
    
    
    private void insertTableReport(TableReport tr, GenotypeExperiment exp, String method) {
        
        Object [][] data = tr.getTableData();
        Object [] colNames = tr.getTableColumnNames();
        
        int traitIndex = findIndex(COLUMN_TRAIT, colNames);
        int locusIndex = findIndex(COLUMN_LOCUS, colNames);
        int siteIndex = findIndex(COLUMN_SITE, colNames);
        int dfIndex = findIndex(COLUMN_DF, colNames);
        int fIndex = findIndex(COLUMN_F, colNames);
        int pIndex = findIndex(COLUMN_P, colNames);
        int modelDfIndex = findIndex(COLUMN_MODEL_DF, colNames);
        int errorDfIndex = findIndex(COLUMN_ERROR_DF, colNames);
        int errorMsIndex = findIndex(COLUMN_ERROR_MS, colNames);
        int rsqModelIndex = findIndex(COLUMN_RSQ_MODEL, colNames);
        int rsqMarkerIndex = findIndex(COLUMN_RSQ_MARKER, colNames);
        
        
        int rows = data.length;
        for(int i=0; i<rows; i++) {
            
            Object trait = data[i][traitIndex];
            Object locus = data[i][locusIndex];
            Object site = data[i][siteIndex];
            Object df = data[i][dfIndex];
            Object f = data[i][fIndex];
            Object p = data[i][pIndex];
            Object modelDf = data[i][modelDfIndex];
            Object errorDf = data[i][errorDfIndex];
            Object errorMs = data[i][errorMsIndex];
            Object rsqModel = data[i][rsqModelIndex];
            Object rsqMarker = data[i][rsqMarkerIndex];
            
            StringBuilder buffer = new StringBuilder();
            
            buffer.append("insert into ");
            buffer.append(ANA_TRAIT_ASSOC);
            buffer.append(" (");
            buffer.append(DIV_ALLELE_ASSAY_ID);
            buffer.append(", ");
            buffer.append(ANA_TRAIT_ASSOC_SITE);
            buffer.append(", ");
            buffer.append(ANA_TRAIT_ASSOC_METHOD);
            buffer.append(", ");
            buffer.append(ANA_TRAIT_ASSOC_DF);
            buffer.append(", ");
            buffer.append(ANA_TRAIT_ASSOC_F_STAT);
            buffer.append(", ");
            buffer.append(ANA_TRAIT_ASSOC_P_VALUE);
            buffer.append(", ");
            buffer.append(ANA_TRAIT_ASSOC_MODEL_DF);
            buffer.append(", ");
            buffer.append(ANA_TRAIT_ASSOC_ERROR_MS);
            buffer.append(", ");
            buffer.append(ANA_TRAIT_ASSOC_ERROR_DF);
            buffer.append(", ");
            buffer.append(ANA_TRAIT_ASSOC_RSQ_MODEL);
            buffer.append(", ");
            buffer.append(ANA_TRAIT_ASSOC_RSQ_MARKER);
            buffer.append(") values (");
            buffer.append(exp.getID().getNumber());
            buffer.append(", ");
            try {
                buffer.append(Integer.parseInt(site.toString()));
            } catch (Exception e) {
                buffer.append(NULL);
            }
            buffer.append(", '");
            buffer.append(method);
            buffer.append("', ");
            try {
                buffer.append(Double.parseDouble(df.toString()));
            } catch (Exception e) {
                buffer.append(NULL);
            }
            buffer.append(", ");
            try {
                buffer.append(Double.parseDouble(f.toString()));
            } catch (Exception e) {
                buffer.append(NULL);
            }
            buffer.append(", ");
            try {
                buffer.append(Double.parseDouble(p.toString()));
            } catch (Exception e) {
                buffer.append(NULL);
            }
            buffer.append(", ");
            try {
                buffer.append(Double.parseDouble(modelDf.toString()));
            } catch (Exception e) {
                buffer.append(NULL);
            }
            buffer.append(", ");
            try {
                buffer.append(Double.parseDouble(errorMs.toString()));
            } catch (Exception e) {
                buffer.append(NULL);
            }
            buffer.append(", ");
            try {
                buffer.append(Double.parseDouble(errorDf.toString()));
            } catch (Exception e) {
                buffer.append(NULL);
            }
            buffer.append(", ");
            try {
                buffer.append(Double.parseDouble(rsqModel.toString()));
            } catch (Exception e) {
                buffer.append(NULL);
            }
            buffer.append(", ");
            try {
                buffer.append(Double.parseDouble(rsqMarker.toString()));
            } catch (Exception e) {
                buffer.append(NULL);
            }
            buffer.append(")");
            
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
