/*
 * GDPCImportPlugin.java
 *
 */

package net.maizegenetics.baseplugins.gdpc;


import gov.usda.gdpc.GenotypeTable;
import gov.usda.gdpc.Import;

import java.net.URL;

import javax.swing.ImageIcon;

import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;

/**
 *
 * @author terryc
 */
public class GDPCImportPlugin extends AbstractPlugin{
    
    public static enum IMPORT_TYPE {fasta};
    
    private IMPORT_TYPE myImportType = IMPORT_TYPE.fasta;
    
    private String myFilename = null;
    
    
    /** Creates a new instance of GDPCImportPlugin */
    public GDPCImportPlugin() {
        super(null, false);
    }
    
    
    /**
     * Performs import of file specified by setFileName() according
     * to IMPORT_TYPE specified by setImportType().
     *
     * @param input not used.
     */
    public DataSet performFunction(DataSet input) {
        return performFunction(myImportType, myFilename);
    }
    
    
    /**
     * Performs import of file specified according
     * to IMPORT_TYPE specified.
     *
     * @param type type of import.
     * @param file file to import.
     */
    public DataSet performFunction(IMPORT_TYPE type, String file) {
        
        try {
            
            if (type == IMPORT_TYPE.fasta) {
                GenotypeTable table = Import.readFasta(file);
                DataSet result = new DataSet(new Datum("GenotypeTable", table, null), this);
                fireDataSetReturned(result);
                return result;
            } else {
                return null;
            }
            
        } catch (Exception e) {
            e.printStackTrace();
            return null;
        }
        
    }
    
    
    public void setImportType(IMPORT_TYPE type) {
        myImportType = type;
    }
    
    
    public void setFileName(String filename) {
        myFilename = filename;
    }
    
    
    public ImageIcon getIcon() {
        
        URL imageURL = GDPCImportPlugin.class.getResource("gdpc.gif");
        
        if(imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
        
    }
    
    public String getToolTipText() {
        return "GDPC Import";
    }
    
    public String getButtonName() {
        return "GDPC Import";
    }
    
}
