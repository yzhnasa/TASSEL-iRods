package net.maizegenetics.analysis.gbs;

import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.util.ArgsEngine;
import net.maizegenetics.util.DirectoryCrawler;
import net.maizegenetics.util.MultiMemberGZIPInputStream;
import org.apache.log4j.Logger;

import javax.swing.*;
import java.awt.*;
import java.io.*;

/** 
 * Split large TagsByTaxaByteHDF5TagGroup file into small sub TBTs. Designed to submit genetic mapping jobs in cluster
 * 
 * @author Fei Lu
 */
public class PanAModelTrainingPlugin extends AbstractPlugin {

    static long timePoint1;
    private ArgsEngine engine = null;
    private Logger logger = Logger.getLogger(PanAModelTrainingPlugin.class);
    
    String trainingSetFileS = null;
    String wekaPath = null;
    String modelFileS = null;
    String reportDirS = null;

    public PanAModelTrainingPlugin() {
        super(null, false);
    }

    public PanAModelTrainingPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    private void printUsage() {
        logger.info(
                "\n\nUsage is as follows:\n"
                + " -t  training data set file, including unique reference tag\n"
                + " -w  path of weka library file\n"        
                + " -m  M5 model file\n"
                + " -r  directory of training report\n");
    }

    public DataSet performFunction(DataSet input) {
           
        return null;
    }

    @Override
    public void setParameters(String[] args) {
        if (args.length == 0) {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }
        
        if (engine == null) {
            engine = new ArgsEngine();
            engine.add("-t", "--training-data", true);
            engine.add("-w", "--weka-path", true);
            engine.add("-m", "--output-model", true);
            engine.add("-r", "--report-dir", true);
            engine.parse(args);
        }

        if (engine.getBoolean("-t")) {
            this.trainingSetFileS = engine.getString("-t");
        }
        else {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }

        if (engine.getBoolean("-w")) {
            this.wekaPath = engine.getString("-w");
        } 
        else {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }
        
        if (engine.getBoolean("-m")) {
            modelFileS = engine.getString("-m");
        } 
        else {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }
        
        if (engine.getBoolean("-r")) {
            reportDirS = engine.getString("-r");
        } 
        else {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }
    }

    @Override
    public ImageIcon getIcon() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String getButtonName() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String getToolTipText() {
        throw new UnsupportedOperationException("Not supported yet.");
    }
}
