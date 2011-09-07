/*
 * SNPFileConverstionPipeline.java
 *
 * Created on January 15th, 2010
 */

package net.maizegenetics.pipeline;

import org.apache.log4j.PropertyConfigurator;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.ImportUtils;
import net.maizegenetics.pal.alignment.ExportUtils;

/**
 *
 * @author Jon
 */
public class AlleleFileConversionPipeline{

    protected String myMAPFile = null;
    protected String myGenotypeFile = null;

    protected String myOutputFileHapmap = null;
    protected String myOutputFilePlink = null;
    protected String myOutputFileFlapjack = null;
    protected String myOutputFileZip = null;

    protected Boolean myInputHapmap = false;
    protected Boolean myInputPlink = false;
    protected Boolean myInputFlapjack = false;
    protected Boolean myInputZip = false;

    protected Boolean myOutputHapmap = false;
    protected Boolean myOutputPlink = false;
    protected Boolean myOutputFlapjack = false;
    protected Boolean myOutputZip = false;


    /**
     * Creates a new instance of AlleleFileConversionPipeline
     */
    public AlleleFileConversionPipeline(String args[]) {
        java.util.Properties props = new java.util.Properties();
        props.setProperty("log4j.logger.net.maizegenetics", "INFO, stdout");
        props.setProperty("log4j.appender.stdout",
                "org.apache.log4j.ConsoleAppender");
        props.setProperty("log4j.appender.stdout.layout",
                "org.apache.log4j.TTCCLayout");
        PropertyConfigurator.configure(props);

        parseArgs(args);

        checkArgs();

        init();
    }

    public void parseArgs(String[] args) {

        int index = 0;

        while (index < args.length) {

            try {

                String current = args[index++];

                if (!current.startsWith("-")) {
                    printUsage();
                    System.exit(0);
                }

                if (current.equalsIgnoreCase("-ih")) {
                    if (!myInputHapmap && !myInputPlink && !myInputFlapjack && !myInputZip) {
                        myInputHapmap = true;
                        myMAPFile = args[index++].trim();
                    }
                    else {
                        printUsage();
                        System.out.println("Error: Can only specify one set of input files.  " + current + " is invalid.");
                        System.exit(0);
                    }
                }
                else if (current.equalsIgnoreCase("-ip")) {
                    if (!myInputHapmap && !myInputPlink && !myInputFlapjack && !myInputZip) {
                        myInputPlink = true;
                        myMAPFile = args[index++].trim();
                        myGenotypeFile = args[index++].trim();
                    }
                    else {
                        printUsage();
                        System.out.println("Error: Can only specify one set of input files.  " + current + " is invalid.");
                        System.exit(0);
                    }
                }
                else if (current.equalsIgnoreCase("-if")) {
                    if (!myInputHapmap && !myInputPlink && !myInputFlapjack && !myInputZip) {
                        myInputFlapjack = true;
                        myMAPFile = args[index++].trim();
                        myGenotypeFile = args[index++].trim();
                    }
                    else {
                        printUsage();
                        System.out.println("Error: Can only specify one set of input files.  " + current + " is invalid.");
                        System.exit(0);
                    }
                }
                else if (current.equalsIgnoreCase("-iz")) {
                    if (!myInputHapmap && !myInputPlink && !myInputFlapjack && !myInputZip) {
                        myInputZip = true;
                        myMAPFile = args[index++].trim();
                    }
                    else {
                        printUsage();
                        System.out.println("Error: Can only specify one set of input files.  " + current + " is invalid.");
                        System.exit(0);
                    }
                }
                else if (current.equalsIgnoreCase("-oh")) {
                    if (!myOutputHapmap) {
                        myOutputHapmap = true;
                        myOutputFileHapmap = args[index++].trim();
                    }
                    else {
                        printUsage();
                        System.out.println("Error: Already selected Hapmap format to output.");
                        System.exit(0);
                    }
                }
                else if (current.equalsIgnoreCase("-op")) {
                    if (!myOutputPlink) {
                        myOutputPlink = true;
                        myOutputFilePlink = args[index++].trim();
                    }
                    else {
                        printUsage();
                        System.out.println("Error: Already selected Plink format to output.");
                        System.exit(0);
                    }
                }
                else if (current.equalsIgnoreCase("-of")) {
                    if (!myOutputFlapjack) {
                        myOutputFlapjack = true;
                        myOutputFileFlapjack = args[index++].trim();
                    }
                    else {
                        printUsage();
                        System.out.println("Error: Already selected Flapjack format to output.");
                        System.exit(0);
                    }
                }
                else if (current.equalsIgnoreCase("-oz")) {
                    if (!myOutputZip) {
                        myOutputZip = true;
                        myOutputFileZip = args[index++].trim();
                    }
                    else {
                        printUsage();
                        System.out.println("Error: Already selected Zip as output");
                        System.exit(0);
                    }
                }
                else {
                    printUsage();
                    System.out.println("Error: Unknown parameter: " + current);
                    System.exit(0);
                }

            }

            catch (Exception e) {
                printUsage();
                System.exit(0);
            }
        }
    }

    private void init() {
        Alignment a = null;
        if (myInputHapmap) {
            a = ImportUtils.readFromHapmap(myMAPFile);
        }
        else if (myInputPlink) {
            a = ImportUtils.readFromPLINK(myGenotypeFile, myMAPFile);
        }
        else if (myInputFlapjack) {
            a = ImportUtils.readFromFlapjack(myGenotypeFile, myMAPFile);
        }
        else if (myInputZip) {
            a = ImportUtils.readFromZip(myMAPFile);
        }
        else {
            printUsage();
            System.out.println("Unknown Error.");
            System.exit(0);
        }

        if (myOutputHapmap) {
            ExportUtils.writeToHapmap(a, true, myOutputFileHapmap, '\t');
        }
        if (myOutputPlink) {
            ExportUtils.writeToPlink(a, myOutputFilePlink, '\t');
        }
        if (myOutputFlapjack) {
            ExportUtils.writeToFlapjack(a, myOutputFileFlapjack, '\t');
        }
        if (myOutputZip) {
            ExportUtils.writeToZip(a, myOutputFileZip);
        }
    }

    public static void main(String args[]) {
        AlleleFileConversionPipeline pipeline = new AlleleFileConversionPipeline(args);
    }

    protected void checkArgs() {
        if (!myInputHapmap && !myInputPlink && !myInputFlapjack && !myInputZip) {
            printUsage();
            System.out.println("Error: Input File must be specified.");
            System.exit(0);
        }
        else if (!myOutputHapmap && !myOutputPlink && !myOutputFlapjack && !myOutputZip) {
            printUsage();
            System.out.println("Error: Output File must be specified.");
            System.exit(0);
        }
    }

    protected void printUsage() {
        System.out.println("\nUsage: java net.maizegenetics.pipeline.AlleleFileConversionPipeline (-ih <hapmap file> | -ip <MAP file> <PED file> | -if <MAP file> <Genotype file> | -iz <Zip file>) [-oh <hapmap output file name>] [-op <plink output file names>] [-of <flapjack output file names>] [-oz <zip output file name>]");
        System.out.println("");
        System.out.println("  -ih <hapmap file>                 : Selects input file as hapmap");
        System.out.println("  -ip <MAP file> <PED file>         : Selects input files as plink");
        System.out.println("  -if <MAP file> <Genotype file>    : Selects input files as flapjack");
        System.out.println("  -iz <GZip file>                   : Selects input file as GZip");
        System.out.println("  -oh <hapmap output file name>     : Optionally specify file name for hapmap output");
        System.out.println("  -op <plink output file names>     : Optionally specify name for MAP and PED files");
        System.out.println("  -of <flapjack output file names>  : Optionally specify name for MAP and Genotype files");
        System.out.println("  -oz <GZip output file name>       : Optionally specify name for GZip file, one file created per chromosome");
        System.out.println("");
        System.out.println("Not required to specify file suffixes. At least one output file type must be selected.");
        System.out.println("");
    }
}