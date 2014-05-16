/*
 * UTagCountToTagPairPlugin
 */
package net.maizegenetics.analysis.gbs;

import java.awt.Frame;
import java.io.File;
import javax.swing.ImageIcon;
import net.maizegenetics.util.ArgsEngine;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import org.apache.log4j.Logger;

/**
 *
 * @author Fei Lu
 */
public class UTagCountToTagPairPlugin extends AbstractPlugin {

    private ArgsEngine engine = null;
    private Logger logger = Logger.getLogger(UTagCountToTagPairPlugin.class);
    private String parentDir = null;
    private String infile, outfile;
    private boolean useWorkDir = true;
    private double etr = 0.03;

    public UTagCountToTagPairPlugin() {
        super(null, false);
    }

    public UTagCountToTagPairPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    private void printUsage() {
        logger.info(
                "\n\nUsage is as follows:\n"
                + " -e  --error-tolerance-rate  Error tolerance rate in the network filter. (Default: " + etr + ")\n"
                + " -i  --input                 Input file of merged tag counts (required)\n"
                + " -o  --output                Output file of tag pairs (required)\n");
    }

    @Override
    public DataSet performFunction(DataSet input) {
        File pd;
        String mergedTagCountOfAllS, tagPairS;
        mergedTagCountOfAllS = new File(infile).getAbsolutePath();
        tagPairS = new File(outfile).getAbsolutePath();
        UNetworkFilter unf = new UNetworkFilter(mergedTagCountOfAllS, etr, tagPairS);
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
            engine.add("-e", "--error-tolerance-rate ", true);
            engine.add("-i", "--input", true);
            engine.add("-o", "--output", true);
            engine.parse(args);
        }

        //Check for input file
        if (engine.getBoolean("-i")) {
            infile = engine.getString("-i");
        } else {
            throw new IllegalArgumentException("\n\nMust supply input file name.\n\n");
        }

        //Check for output file
        if (engine.getBoolean("-o")) {
            outfile = engine.getString("-o");
        } else {
            throw new IllegalArgumentException("\n\nMust supply output file name.\n\n");
        }

        //Check for error tolerance
        if (engine.getBoolean("-e")) {
            etr = Double.parseDouble(engine.getString("-e"));
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
