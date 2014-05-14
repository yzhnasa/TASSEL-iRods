/*
 * ProductionPipeline
 */
package net.maizegenetics.analysis.gbs;

import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.PluginParameterTerry;
import net.maizegenetics.util.Utils;
import net.maizegenetics.util.DirectoryCrawler;

import org.apache.log4j.Logger;

import javax.swing.*;
import java.awt.*;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;

import java.text.SimpleDateFormat;
import java.util.Arrays;
import java.util.Date;
import org.apache.log4j.PropertyConfigurator;

/**
 *
 * @author Terry Casstevens
 */
public class ProductionPipeline extends ParameterConceptPlugin {

    private Logger myLogger = Logger.getLogger(ProductionPipeline.class);
    private static final SimpleDateFormat LOGGING_DATE_FORMAT = new SimpleDateFormat("yyyyMMdd HH:mm:ss");

    public enum PARAMETERS {

        inputDirectory, keyFile, enzyme, productionTOPM, outputGenotypeFile, archiveDirectory
    };

    protected PluginParameterTerry<String> myInputDirectory = new PluginParameterTerry.Builder<String>(PARAMETERS.inputDirectory, null, String.class).required(true).inFile()
            .description("Input directory containing fastq AND/OR qseq files").build();
    protected PluginParameterTerry<String> myKeyFile = new PluginParameterTerry.Builder<String>(PARAMETERS.keyFile, null, String.class).required(true).inFile()
            .description("Barcode Key File").build();
    protected PluginParameterTerry<String> myEnzyme = new PluginParameterTerry.Builder<String>(PARAMETERS.enzyme, null, String.class).required(true)
            .description("Enzyme used to create the GBS library").build();
    protected PluginParameterTerry<String> myProductionTOPM = new PluginParameterTerry.Builder<String>(PARAMETERS.productionTOPM, null, String.class).required(true).inFile()
            .description("Physical map file containing tags and corresponding variants (production TOPM)").build();
    protected PluginParameterTerry<String> myOutputGenotypeFile = new PluginParameterTerry.Builder<String>(PARAMETERS.outputGenotypeFile, null, String.class).required(true).outFile()
            .description("Output (target) HDF5 genotypes file to add new genotypes to (new file created if it doesn't exist)").build();
    protected PluginParameterTerry<String> myArchiveDirectory = new PluginParameterTerry.Builder<String>(PARAMETERS.archiveDirectory, null, String.class).required(true).outFile()
            .description("Archive directory where to move processed files").build();

    private String myOutputDirectory;
    private PrintStream myPrintStreamToLog;

    public ProductionPipeline(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    public DataSet performFunction(DataSet input) {
        myOutputDirectory = Utils.getDirectory(myOutputGenotypeFile.value());
        setupLogfile();
        myLogger.info(getTimeStamp());
        return super.performFunction(input);
    }

    @Override
    public DataSet processData(DataSet input) {

        try {
            String[] rawSeqFileNames = DirectoryCrawler.listFileNames(ProductionSNPCallerPlugin.rawSeqFileNameRegex, myInputDirectory.value());
            if ((rawSeqFileNames == null) || (rawSeqFileNames.length == 0)) {
                return null;
            }

            myLogger.info("Raw Sequence Files: " + Arrays.deepToString(rawSeqFileNames));
            myLogger.info("Parameters Passed to ProductionSNPCallerPlugin: " + Arrays.deepToString(getPluginArgs()));
            
            ProductionSNPCallerPlugin plugin = new ProductionSNPCallerPlugin();
            plugin.setParameters(getPluginArgs());

            printParameterValues();
            plugin.performFunction(null);

            for (int i = 0; i < rawSeqFileNames.length; i++) {
                File currentLocation = new File(rawSeqFileNames[i]);
                File newLocation = new File(myArchiveDirectory.value() + Utils.getFilename(rawSeqFileNames[i]));
                currentLocation.renameTo(newLocation);
                myLogger.info("Moved : " + currentLocation.getAbsolutePath() + " to: " + newLocation.getAbsolutePath());
            }

            return null;
        } finally {
            if (myPrintStreamToLog != null) {
                myPrintStreamToLog.close();
            }
        }

    }

    private String[] getPluginArgs() {
        String[] args = {
            "-i", myInputDirectory.value(),
            "-k", myKeyFile.value(),
            "-e", myEnzyme.value(),
            "-o", myOutputGenotypeFile.value(),
            "-m", myProductionTOPM.value()
        };
        return args;
    }

    private void setupLogfile() {

        String todayDate = new SimpleDateFormat("yyyyMMdd").format(new Date());
        String logFileName = todayDate + "_" + "ProductionPipeline" + ".log";

        File logFile = new File(myOutputDirectory + "/" + logFileName);
        myLogger.info("Log File: " + logFile.getAbsolutePath());
        try {
            logFile.createNewFile();
            myPrintStreamToLog = new PrintStream(new BufferedOutputStream(new FileOutputStream(logFile, true)));
            System.setOut(myPrintStreamToLog);
            System.setErr(myPrintStreamToLog);
        } catch (IOException ioe) {
            ioe.printStackTrace();
        }

        java.util.Properties props = new java.util.Properties();
        props.setProperty("log4j.logger.net.maizegenetics", "DEBUG, FILE");
        props.setProperty("log4j.appender.FILE", "org.apache.log4j.FileAppender");
        props.setProperty("log4j.appender.FILE.File", logFile.getAbsolutePath());
        props.setProperty("log4j.appender.FILE.ImmediateFlush", "true");
        props.setProperty("log4j.appender.FILE.Threshold", "debug");
        props.setProperty("log4j.appender.FILE.Append", "true");
        props.setProperty("log4j.appender.FILE.layout", "org.apache.log4j.TTCCLayout");
        PropertyConfigurator.configure(props);

    }

    /**
     * Convenience method to provide uniformly labeled timestamps
     */
    private static String getTimeStamp() {
        return "Timestamp: " + LOGGING_DATE_FORMAT.format(new Date()) + ": ";
    }

    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "Production Pipeline";
    }

    @Override
    public String getToolTipText() {
        return "Production Pipeline";
    }
}
