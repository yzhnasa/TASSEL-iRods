/*
 * WriteXLSPlugin.java
 *
 * Created on September 5, 2007
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
import gov.usda.gdpc.util.Calendar;

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
import java.io.FileOutputStream;
import java.io.FileInputStream;

import javax.swing.ImageIcon;

import org.apache.poi.hssf.usermodel.HSSFCell;
import org.apache.poi.hssf.usermodel.HSSFCellStyle;
import org.apache.poi.hssf.usermodel.HSSFDataFormat;
import org.apache.poi.hssf.usermodel.HSSFRow;
import org.apache.poi.hssf.usermodel.HSSFSheet;
import org.apache.poi.hssf.usermodel.HSSFWorkbook;
import org.apache.poi.poifs.filesystem.POIFSFileSystem;


/**
 *
 * @author terryc
 */
public class WriteXLSPlugin extends AbstractPlugin {
    
    private String myFilename = null;
    private String myFilenameBase = null;
    
    private List myGenotypeExpProps = new ArrayList();
    
    private HSSFCellStyle myDateCellStyle = null;
    
    private Map myColumnHeadings = new HashMap();
    
    private int [] myColumnTypes = new int [10];
    
    private HSSFSheet mySheet = null;
    
    private HSSFRow myColumnHeadingsRow = null;
    
    private short myHeadingIndex = 0;
    
    private HSSFWorkbook myWorkbook = null;
    
    private int mySaveCount = 0;
    
    
    /**
     * Creates a new instance of WriteXLSPlugin
     */
    public WriteXLSPlugin(String filename) {
        super(null, false);
        myFilenameBase = filename;
        myGenotypeExpProps.add(GenotypeExperimentProperty.NAME);
    }
    
    
    public WriteXLSPlugin() {
        super(null, false);
        myGenotypeExpProps.add(GenotypeExperimentProperty.NAME);
    }
    
    
    private void init() {
        
        myColumnHeadings = new HashMap();
        myColumnTypes = new int [10];
        myColumnHeadingsRow = null;
        myDateCellStyle = null;
        myHeadingIndex = 0;
        myWorkbook = null;
        mySheet = null;
        mySaveCount = 0;
        
        if (myFilename == null) {
            return;
        }
        
        try {
            
            if (myWorkbook == null) {
                File xlsFile = new File(myFilename);
                
                if (!xlsFile.exists()) {
                    createXLSFile();
                } else {
                    
                    FileInputStream inputStream = new FileInputStream(myFilename);
                    POIFSFileSystem fs = new POIFSFileSystem(inputStream);
                    myWorkbook = new HSSFWorkbook(fs);
                    inputStream.close();
                    
                    mySheet = myWorkbook.getSheetAt(0);
                    
                    myColumnHeadingsRow = mySheet.getRow(mySheet.getFirstRowNum());
                    if (myColumnHeadingsRow == null) {
                        myColumnHeadingsRow = mySheet.createRow(mySheet.getFirstRowNum());
                    }
                    
                    if (myDateCellStyle == null) {
                        myDateCellStyle = myWorkbook.createCellStyle();
                        myDateCellStyle.setDataFormat(HSSFDataFormat.getBuiltinFormat("m/d/yy"));
                    }
                    
                }
                
            }
            
        } catch (Exception e) {
            e.printStackTrace();
        }
        
    }
    
    
    public DataSet performFunction(DataSet input) {
        
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
                filename = myFilenameBase + "_" + name + ".xls";
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
            HSSFRow row = null;
            while (notFinished) {
                
                try {
                    
                    row = mySheet.createRow(mySheet.getLastRowNum() + 1);
                    
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
            
            if (mySaveCount++ == 100) {
                saveFile();
                mySaveCount = 0;
            }
            
        } catch (Exception e) {
            e.printStackTrace();
        }
        
    }
    
    
    private void createXLSFile() {
        System.out.println("CREATING: " + myFilename);
        myWorkbook = new HSSFWorkbook();
        mySheet = myWorkbook.createSheet("Assoc. Analysis");
        myColumnHeadingsRow = mySheet.createRow(mySheet.getFirstRowNum());
        myDateCellStyle = myWorkbook.createCellStyle();
        myDateCellStyle.setDataFormat(HSSFDataFormat.getBuiltinFormat("m/d/yy"));
        saveFile();
    }
    
    
    private short getColumnIndex(Class type, String name) {
        
        if ((type == null) || (name == null)) {
            return -1;
        }
        
        
        
        Map typeMap = (Map) myColumnHeadings.get(type);
        
        HSSFCell cell = null;
        
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
            cell = getNewHeadingCell(name);
            typeMap.put(nameUpperCase, cell);
            setColumnType(cell, type, nameUpperCase);
            
        } else {
            
            cell = (HSSFCell) typeMap.get(nameUpperCase);
            
            if (cell == null) {
                cell = getNewHeadingCell(name);
                typeMap.put(nameUpperCase, cell);
                setColumnType(cell, type, nameUpperCase);
            }
            
        }
        
        return cell.getCellNum();
        
    }
    
    
    private void setColumnType(HSSFCell cell, Class type, String name) {
        
        if (cell.getCellNum() > myColumnTypes.length - 1) {
            increaseColumnTypesArraySize(cell.getCellNum());
        }
        
        if (type.isInstance(TableReport.class)) {
            if ((name.equals("P")) || (name.equals("DF")) || (name.equals("F"))
            || (name.equals("DF_MODEL")) || (name.equals("DF_ERROR"))
            || (name.equals("MS_ERROR")) || (name.equals("RSQ_MODEL"))
            || (name.equals("RSQ_MARKER")) || (name.equals("DF_MARKER"))
            || (name.equals("F_MARKER")) || (name.equals("P_MARKER"))) {
                myColumnTypes[cell.getCellNum()] = HSSFCell.CELL_TYPE_NUMERIC;
            }
        } else {
            myColumnTypes[cell.getCellNum()] = HSSFCell.CELL_TYPE_STRING;
        }
    }
    
    
    private void increaseColumnTypesArraySize(int hint) {
        int newSize = hint + 10;
        int [] result = new int[newSize];
        System.arraycopy(myColumnTypes, 0, result, 0, myColumnTypes.length);
        myColumnTypes = result;
    }
    
    
    private HSSFCell getNewHeadingCell(String name) {
        
        HSSFCell cell = myColumnHeadingsRow.createCell(myHeadingIndex++);
        cell.setCellValue(name);
        
        return cell;
        
    }
    
    
    private void insertTableReportRow(Object [] colNames, Object [] data, HSSFRow row) {
        
        for (int i=0; i<data.length; i++) {
            short index = getColumnIndex(TableReport.class, colNames[i].toString());
            HSSFCell cell = row.createCell(index);
            setCell(cell, data[i]);
        }
        
    }
    
    
    private void insertOtherFields(List otherFields, HSSFRow row) {
        
        Iterator itr = otherFields.iterator();
        while (itr.hasNext()) {
            Datum current = (Datum) itr.next();
            String name = current.getName();
            Object value = current.getData();
            short index = getColumnIndex(TableReport.class, name);
            HSSFCell cell = row.createCell(index);
            setCell(cell, value);
        }
        
    }
    
    
    private void insertExperiment(GenotypeExperiment exp, HSSFRow row) {
        
        Iterator itr = myGenotypeExpProps.iterator();
        while (itr.hasNext()) {
            GenotypeExperimentProperty current = (GenotypeExperimentProperty) itr.next();
            short index = getColumnIndex(GenotypeExperiment.class, current.getName());
            HSSFCell cell = row.createCell(index);
            setCell(cell, exp.getProperty(current));
        }
        
    }
    
    
    private void setCell(HSSFCell cell, Object data) {
        
        cell.setCellType(myColumnTypes[cell.getCellNum()]);
        
        if (data == null) {
            cell.setCellValue("");
        } else if (data instanceof Calendar) {
            cell.setCellStyle(myDateCellStyle);
            cell.setCellValue((Calendar) data);
        } else {
            try {
                cell.setCellValue(Double.parseDouble(data.toString()));
            } catch (Exception e) {
                cell.setCellValue(data.toString());
            }
        }
        
    }
    
    
    private void insertPhenotypes(PhenotypeTable phenotypes, HSSFRow row) {
        
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
        HSSFCell cell = row.createCell(index);
        setCell(cell, builder.toString());
        
    }
    
    
    private void insertGenotypes(GenotypeTable genotypes, GenotypeExperiment exp, HSSFRow row) {
        
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
        HSSFCell cell = row.createCell(index);
        setCell(cell, builder.toString());
        
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
    
    
    public void saveFile() {
        
        if (myFilename == null) {
            return;
        }
        
        try {
            FileOutputStream output = new FileOutputStream(myFilename);
            myWorkbook.write(output);
            output.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
        
    }
    
    
    protected void finalize() throws Throwable {
        
        try {
            saveFile();
        } finally {
            super.finalize();
        }
        
    }
    
}
