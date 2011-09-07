/*
 * WriteDelimitedPlugin.java
 *
 * Created on June 23, 2008
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

import net.maizegenetics.pal.report.TableReport;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;

import java.util.Arrays;
import java.util.ArrayList;
import java.util.List;
import java.util.HashMap;
import java.util.Map;

import java.io.File;
import java.io.FileWriter;
import java.io.Writer;

import javax.swing.ImageIcon;



/**
 *
 * @author terryc
 */
public class WriteDelimitedPlugin extends AbstractPlugin {
    
    private String myFilename = null;
    private String myFilenameBase = null;
    
    private List myGenotypeExpProps = new ArrayList();
    
    private Map myColumnHeadings = new HashMap();
    
    private List myHeadings = new ArrayList();
    
    private int mySaveCount = 0;
    
    private Writer myWriter = null;
    
    private boolean myHeadingsWritten = false;
    private int myNumHeadingsLastWritten = 0;
    
    private String myDelimiter = "\t";
    
    
    /**
     * Creates a new instance of WriteDelimitedPlugin
     */
    public WriteDelimitedPlugin(String filename) {
        super(null, false);
        myFilenameBase = filename;
        myGenotypeExpProps.add(GenotypeExperimentProperty.NAME);
    }
    
    
    public WriteDelimitedPlugin() {
        super(null, false);
        myGenotypeExpProps.add(GenotypeExperimentProperty.NAME);
    }
    
    
    private void init() {
        
        myColumnHeadings = new HashMap();
        myHeadings = new ArrayList();
        mySaveCount = 0;
        myHeadingsWritten = false;
        
        if (myFilename == null) {
            return;
        }
        
        try {
            
            File tabFile = new File(myFilename);
            myWriter = new FileWriter(tabFile, true);
            
        } catch (Exception e) {
            e.printStackTrace();
        }
        
    }
    
    
    public DataSet performFunction(DataSet input) {
        
        System.out.println("WriteDelimitedPlugin: performFunction: input start");
        System.out.println(input.toString());
        System.out.println("WriteDelimitedPlugin: performFunction: input end");
        
        List <Datum> reports = input.getDataOfType(TableReport.class);
        
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
        
        for (int i=0; i<reports.size(); i++) {
            
            TableReport report = (TableReport) reports.get(i).getData();
            String name = reports.get(i).getName();
            int periodIndex = myFilenameBase.lastIndexOf(".");
            String filename = null;
            if (periodIndex == -1) {
                filename = myFilenameBase + "_" + name + ".txt";
            } else {
                filename = myFilenameBase.substring(0, periodIndex) + "_" + name + myFilenameBase.substring(periodIndex);
            }
            setFullFilename(filename);
            
            insertData(report, exp, genotypeTable, phenotypeTable, otherFields);
            
            saveFile();
            
        }

        return null;
        
    }
    
    
    private void insertData(TableReport tr, GenotypeExperiment exp, GenotypeTable genotypes, PhenotypeTable phenotypes, List otherFields) {
        
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
                
                List row = new ArrayList();
                
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
                
                if (!myHeadingsWritten) {
                    writeHeadings();
                    myHeadingsWritten = true;
                }
                
                for (int i=0, m=row.size(); i<m; i++) {
                    if (i != 0) {
                        myWriter.write(myDelimiter);
                    }
                    myWriter.write(row.get(i).toString());
                }
                
                myWriter.write("\n");
                
            }
            
            if (mySaveCount++ == 100) {
                myWriter.flush();
                mySaveCount = 0;
            }
            
        } catch (Exception e) {
            e.printStackTrace();
        }
        
    }
    
    
    private short getColumnIndex(Class type, String name) {
        
        if ((type == null) || (name == null)) {
            return -1;
        }
        
        
        
        Map typeMap = (Map) myColumnHeadings.get(type);
        
        short cell = -1;
        
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
            myHeadings.add(name);
            cell = (short) (myHeadings.size() - 1);
            typeMap.put(nameUpperCase, cell);
            
        } else {
            
            Short temp = (Short) typeMap.get(nameUpperCase);
            if (temp == null) {
                myHeadings.add(name);
                cell = (short) (myHeadings.size() - 1);
                typeMap.put(nameUpperCase, cell);
            } else {
                cell = temp.shortValue();
            }
            
        }
        
        return cell;
        
    }
    
    
    private void insertTableReportRow(Object [] colNames, Object [] data, List row) {
        
        for (int i=0; i<data.length; i++) {
            short index = getColumnIndex(TableReport.class, colNames[i].toString());
            row.add(index, data[i]);
        }
        
    }
    
    
    private void insertOtherFields(List otherFields, List row) {
        
        int size = otherFields.size();
        for (int i=0; i<size; i++) {
            Datum current = (Datum) otherFields.get(i);
            String name = current.getName();
            Object value = current.getData();
            short index = getColumnIndex(TableReport.class, name);
            row.add(index, value);
        }
        
    }
    
    
    private void insertExperiment(GenotypeExperiment exp, List row) {
        
        int size = exp.size();
        for (int i=0; i<size; i++) {
            GenotypeExperimentProperty current = (GenotypeExperimentProperty) exp.get(i);
            short index = getColumnIndex(TableReport.class, current.getName());
            row.add(index, exp.getProperty(current));
        }
        
    }
    
    
    private void insertPhenotypes(PhenotypeTable phenotypes, List row) {
        
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
        
        short index = getColumnIndex(PhenotypeTable.class, "Environments");
        row.add(index, builder.toString());
        
    }
    
    
    private void insertGenotypes(GenotypeTable genotypes, GenotypeExperiment exp, List row) {
        
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
        
        short index = getColumnIndex(GenotypeTable.class, "Taxa");
        row.add(index, builder.toString());
        
    }
    
    
    public void setGenotypeExperimentProperties(GenotypeExperimentProperty [] props) {
        myGenotypeExpProps.clear();
        myGenotypeExpProps.addAll(Arrays.asList(props));
    }
    
    
    public void setFilename(String filename) {
        myFilenameBase = filename;
    }
    
    
    private void setFullFilename(String filename) {
        saveFile();
        myFilename = filename;
        init();
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
    
    
    private void writeHeadings() {
        
        if (myHeadings.size() == myNumHeadingsLastWritten) {
            return;
        }
        
        try {
            
            for (int i=0, m=myHeadings.size(); i<m; i++) {
                if (i != 0) {
                    myWriter.write(myDelimiter);
                }
                myWriter.write(myHeadings.get(i).toString());
            }
            
            myWriter.write("\n");
            
            myNumHeadingsLastWritten = myHeadings.size();
            
        } catch (Exception e) {
            // do nothing
        }
        
    }
    
    
    public void saveFile() {
        
        try {
            
            if (myWriter != null) {
                
                writeHeadings();
                
                myWriter.flush();
                myWriter.close();
                
                myWriter = null;
                
            }
            
        } catch (Exception e) {
            // do nothing
        }
        
    }
    
    
    /**
     * Sets the delimiter used between column values.   The
     * default in the tab character.
     *
     * @parm delimiter delimiter
     */
    public void setDelimiter(String delimiter) {
        myDelimiter = delimiter;
    }
    
    
    protected void finalize() throws Throwable {
        
        try {
            saveFile();
        } finally {
            super.finalize();
        }
        
    }
    
}
