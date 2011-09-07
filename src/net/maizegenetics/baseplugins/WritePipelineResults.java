/*
 * WritePipelineResults.java
 *
 * Created on October 30, 2007
 *
 */

package net.maizegenetics.baseplugins;


import gov.usda.gdpc.AlleleList;
import gov.usda.gdpc.EnvironmentExperiment;
import gov.usda.gdpc.EnvironmentExperimentGroup;
import gov.usda.gdpc.GenotypeExperiment;
import gov.usda.gdpc.GenotypeExperimentProperty;
import gov.usda.gdpc.GenotypeTable;
import gov.usda.gdpc.PhenotypeTable;
import java.io.FileWriter;

import net.maizegenetics.pal.report.TableReport;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;

import java.util.Arrays;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.HashMap;
import java.util.Map;

import java.io.File;

import javax.swing.ImageIcon;


/**
 *
 * @author terryc
 */
public class WritePipelineResults extends AbstractPlugin {
    
    private String myFilename = null;
    
    private List myGenotypeExpProps = new ArrayList();
    
    private Map myColumnHeadings = new HashMap();
    
    private int [] myColumnTypes = new int [10];
    
    private List myColumnHeadingsRow = new ArrayList();
    
    private short myHeadingIndex = 0;
    
    
    /**
     * Creates a new instance of WritePipelineResults
     */
    public WritePipelineResults(String filename) {
        super(null, false);
        myFilename = filename;
        myGenotypeExpProps.add(GenotypeExperimentProperty.NAME);
    }
    
    
    public DataSet performFunction(DataSet input) {
        
        TableReport report = null;
        List <Datum> reports = input.getDataOfType(TableReport.class);
        if (reports.size() > 0) {
            report = (TableReport) reports.get(0).getData();
        }
        
        GenotypeExperiment exp = null;
        List <Datum> exps = input.getDataOfType(GenotypeExperiment.class);
        if (exps.size() > 0) {
            exp = (GenotypeExperiment) exps.get(0).getData();
        }
        
        GenotypeTable genotypeTable = null;
        List <Datum> genotypes = input.getDataOfType(GenotypeTable.class);
        if (genotypes.size() > 0) {
            genotypeTable = (GenotypeTable) genotypes.get(0).getData();
        }
        
        PhenotypeTable phenotypeTable = null;
        List <Datum> phenotypes = input.getDataOfType(PhenotypeTable.class);
        if (phenotypes.size() > 0) {
            phenotypeTable = (PhenotypeTable) phenotypes.get(0).getData();
        }
        
        List <Datum> otherFields = input.getDataOfType(String.class);
        
        insertData(report, exp, genotypeTable, phenotypeTable, otherFields);

        return null;
        
    }
    
    
    private void insertData(TableReport tr, GenotypeExperiment exp, GenotypeTable genotypes, PhenotypeTable phenotypes, List otherFields) {
        
        
        Object [] row = new Object[100];
        
        try {
            
            
            Object [][] data = null;
            Object [] colNames = null;
            int numTableRows = 0;
            try {
                data = tr.getTableData();
                colNames = tr.getTableColumnNames();
                numTableRows = data.length;
            } catch (Exception e) {
            }
            
            boolean notFinished = true;
            int index = 0;
            
            while (notFinished) {
                
                try {
                    
                    if (exp != null) {
                        insertExperiment(exp, row);
                    }
                    
                    if (genotypes != null) {
                        insertGenotypes(genotypes, exp, row);
                    }
                    
                    if (phenotypes != null) {
                        insertPhenotypes(phenotypes, row);
                    }
                    
                    if (otherFields != null) {
                        insertOtherFields(otherFields, row);
                    }
                    
                    if (numTableRows == 0) {
                        notFinished = false;
                    } else {
                        insertTableReportRow(colNames, data[index++], row);
                    }
                    
                    if (index == numTableRows) {
                        notFinished = false;
                    }
                    
                } catch (Exception ex) {
                    ex.printStackTrace();
                    notFinished = false;
                }
                
            }
            
        } catch (Exception e) {
            e.printStackTrace();
        } finally {
            try {
                
                File tabFile = new File(myFilename);
                
                if (!tabFile.exists()) {
                    FileWriter writer = new FileWriter(tabFile, true);
                    
                    StringBuilder builder = new StringBuilder();
                    for (int i=0; i<myHeadingIndex; i++) {
                        if (myColumnHeadingsRow.get(i) != null) {
                            builder.append(myColumnHeadingsRow.get(i).toString());
                        }
                        builder.append("\t");
                    }
                    builder.append("\n");
                    
                    writer.write(builder.toString());
                    writer.close();
                }
                
                FileWriter writer = new FileWriter(tabFile, true);
                
                StringBuilder builder = new StringBuilder();
                for (int i=0; i<myHeadingIndex; i++) {
                    if (row[i] != null) {
                        builder.append(row[i].toString());
                    }
                    builder.append("\t");
                }
                builder.append("\n");
                
                writer.write(builder.toString());
                writer.close();
                
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
        
    }
    
    
    private int getColumnIndex(Class type, String name) {
        
        if ((type == null) || (name == null)) {
            return -1;
        }
        
        Map typeMap = (Map) myColumnHeadings.get(type);
        
        Integer cell = null;
        
        String nameUpperCase = name.toUpperCase();
        
        if (type.isAssignableFrom(TableReport.class)) {
            
            if (nameUpperCase.equals("P_MARKER")) {
                nameUpperCase = "P";
                name = "p";
            } else if (nameUpperCase.equals("DF_MARKER")) {
                nameUpperCase = "DF";
                name = "df";
            } else if (nameUpperCase.equals("F_MARKER")) {
                nameUpperCase = "F";
                name = "f";
            }
            
        }
        
        if (typeMap == null) {
            
            typeMap = new HashMap();
            myColumnHeadings.put(type, typeMap);
            cell = new Integer(getNewHeadingCell(name));
            typeMap.put(nameUpperCase, cell);
            
        } else {
            
            cell = (Integer) typeMap.get(nameUpperCase);
            
            if (cell == null) {
                cell = new Integer(getNewHeadingCell(name));
                typeMap.put(nameUpperCase, cell);
            }
            
        }
        
        return cell.intValue();
        
    }
    
    
    private void increaseColumnTypesArraySize(int hint) {
        int newSize = hint + 10;
        int [] result = new int[newSize];
        System.arraycopy(myColumnTypes, 0, result, 0, myColumnTypes.length);
        myColumnTypes = result;
    }
    
    
    private int getNewHeadingCell(String name) {
        
        /*
        try {
         
            StringBuilder builder = new StringBuilder();
            builder.append(name);
            builder.append("\t");
         
            boolean exists = (new File(myFilename)).exists();
            if (!exists) {
                builder.append("\n");
            }
         
            RandomAccessFile tabFile = new RandomAccessFile(myFilename, "rws");
            tabFile.readLine();
            tabFile.seek(tabFile.getFilePointer());
            tabFile.writeBytes(builder.toString());
            tabFile.close();
         
        } catch (Exception e) {
            e.printStackTrace();
        }
         */
        
        myColumnHeadingsRow.add(myHeadingIndex, name);
        
        return myHeadingIndex++;
        
    }
    
    
    private void insertTableReportRow(Object [] colNames, Object [] data, Object [] row) {
        
        for (int i=0; i<data.length; i++) {
            int index = getColumnIndex(TableReport.class, colNames[i].toString());
            row[index] = data[i];
        }
        
    }
    
    
    private void insertOtherFields(List otherFields, Object [] row) {
        
        Iterator itr = otherFields.iterator();
        while (itr.hasNext()) {
            Datum current = (Datum) itr.next();
            String name = current.getName();
            Object value = current.getData();
            int index = getColumnIndex(TableReport.class, name);
            row[index] = value;
        }
        
    }
    
    
    private void insertExperiment(GenotypeExperiment exp, Object [] row) {
        
        Iterator itr = myGenotypeExpProps.iterator();
        while (itr.hasNext()) {
            GenotypeExperimentProperty current = (GenotypeExperimentProperty) itr.next();
            int index = getColumnIndex(GenotypeExperiment.class, current.getName());
            row[index] = exp.getProperty(current);
        }
        
    }
    
    
    private void insertPhenotypes(PhenotypeTable phenotypes, Object [] row) {
        
        StringBuilder builder = new StringBuilder();
        
        EnvironmentExperimentGroup exps = phenotypes.getEnvironmentExperimentGroup();
        builder.append(exps.size());
        builder.append(" environments");
        
        for (int i=0; i<exps.size(); i++) {
            builder.append("; ");
            EnvironmentExperiment current = (EnvironmentExperiment) exps.get(i);
            if (current != null) {
                builder.append(" ");
                builder.append(current.toString());
            }
        }
        
        int index = getColumnIndex(PhenotypeTable.class, "Environments");
        row[index] = builder.toString();
        
    }
    
    
    private void insertGenotypes(GenotypeTable genotypes, GenotypeExperiment exp, Object [] row) {
        
        StringBuilder builder = new StringBuilder();
        
        Object [] snps = genotypes.getColumn(exp);
        builder.append(snps.length);
        builder.append(" taxa");
        
        for (int i=0; i<snps.length; i++) {
            builder.append("; ");
            AlleleList current = (AlleleList) snps[i];
            if (current != null) {
                builder.append(genotypes.getTaxon(i).getName());
                builder.append(" (");
                builder.append(current.toString());
                builder.append(")");
            }
        }
        
        int index = getColumnIndex(GenotypeTable.class, "Taxa");
        row[index] = builder.toString();
        
    }
    
    
    public void setGenotypeExperimentProperties(GenotypeExperimentProperty [] props) {
        myGenotypeExpProps.clear();
        myGenotypeExpProps.addAll(Arrays.asList(props));
    }
    
    
    public void setFilename(String filename) {
        myFilename = filename;
        myColumnHeadings = new HashMap();
        myColumnTypes = new int [10];
        myColumnHeadingsRow = new ArrayList();
        myHeadingIndex = 0;
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
