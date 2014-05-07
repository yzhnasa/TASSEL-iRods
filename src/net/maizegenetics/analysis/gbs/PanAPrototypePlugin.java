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
public class PanAPrototypePlugin extends AbstractPlugin {

    static long timePoint1;
    private ArgsEngine engine = null;
    private Logger logger = Logger.getLogger(PanAPrototypePlugin.class);
    
    String tbtHDF5 = null;
    String outDirS = null;
    int chunkSize = 65536;

    public PanAPrototypePlugin() {
        super(null, false);
    }

    public PanAPrototypePlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    private void printUsage() {
        logger.info(
                "\n\nUsage is as follows:\n"
                + " -i  input TagsByTaxa(TBT) file, TagsByTaxaByteHDF5TagGroup format\n"
                + " -s  chunkSize, number of tags in a sub TBT. This determines the mapping calculation time usage in a node/computer. Default = 65536\n"        
                + " -o  output directory of sub TBTs\n");
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
            engine.add("-i", "--input-TBT", true);
            engine.add("-s", "--chunk-size", true);
            engine.add("-o", "--output-dir", true);
            engine.parse(args);
        }

        if (engine.getBoolean("-i")) {
            tbtHDF5 = engine.getString("-i");
        }
        else {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }

        if (engine.getBoolean("-s")) {
            chunkSize = Integer.valueOf(engine.getString("-s"));
        } 
        else {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }
        
        if (engine.getBoolean("-o")) {
            outDirS = engine.getString("-o");
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
