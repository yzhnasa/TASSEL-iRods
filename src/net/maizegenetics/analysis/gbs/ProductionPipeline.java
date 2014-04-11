/*
 * ProductionPipeline
 */
package net.maizegenetics.analysis.gbs;

import com.google.common.base.Splitter;
import com.google.common.collect.Iterables;
import com.google.common.io.Files;
import net.maizegenetics.analysis.imputation.FILLINImputationPlugin;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.prefs.TasselPrefs;
import net.maizegenetics.util.CheckSum;
import net.maizegenetics.util.SMTPClient;
import net.maizegenetics.util.Utils;
import org.apache.commons.io.FilenameUtils;
import org.apache.log4j.Logger;

import java.io.*;
import java.net.InetAddress;
import java.net.UnknownHostException;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.HashMap;
import java.util.Map;
import java.util.Properties;
import org.apache.log4j.PropertyConfigurator;

/**
 * This class is for running the GBS Production Pipeline. It is to be run from
 * within the sTASSEL.jar. The cron job should be set up to run the
 * run_pipeline.pl which has been modified to make this class the main() to be
 * run. The JVM memory settings within run_pipeline.pl should also be adjusted
 * upwards.
 *
 * cron example 20 3 * * * cd /workdir/tassel/tassel4-src && /usr/local/bin/perl
 * /workdir/tassel/tassel4-src/run_prod_cron.pl >>
 * /workdir/tassel/tassel4-src/20130808_cron.log 2>&1
 *
 * 20130718 Note from Jeff Glaubitz: A minor detail: ProductionSNPCallerPlugin
 * needs the key file name to end with "_key.txt". All the output files are
 * named after the key file but replacing "_key.txt" with the appropriate
 * extension.
 *
 * @author Dallas Kroon
 * @author Terry Casstevens
 *
 */
public class ProductionPipeline {

    private static final Logger myLogger = Logger.getLogger(ProductionPipeline.class);

    private static final String EMAIL_ADDRESS_DELIMITER = ";";

    private static final SimpleDateFormat LOGGING_DATE_FORMAT = new SimpleDateFormat("yyyyMMdd HH:mm:ss");

    private static final String RUN_FILE_SUFFIX = ".run";

    // host on which pipeline is running
    private String myApplicationHost = "unknown";

    // default server to be used to send email notifications
    private String myEmailHost = "appsmtp.mail.cornell.edu";

    // pipeline admin default email address
    private String[] myRecipientEmailAddresses = null;

    // directory into which the key and run files are placed
    // after being run a dated directory will be created
    // within this directory and artifacts such as the .log
    // and .xml file will be placed there
    private String myArchiveDirectory = "/SSD/prop_pipeline/arcvtmp/";

    // directory containing the haplotype files to be used
    // in imputation
    private String myHaplotypeDirectory = "/SSD/haplos/";

    // Input variable to ProductionSNPCallerPlugin
    private String myInputFolder = null;
    private String myEnzyme = null;
    private String myTopmFile = null;
    private String myOutputFolder = null;
    private String myKeyFile = null;

    private final boolean myRunImputation = true;

    private String myPropertiesFileContents = null;

    private final Map<String, String> myEmailSubjects = new HashMap<>();

    private static final String EXAMPLE_RUN_FILE
            = "emailHost=appsmtp.mail.cornell.edu\n"
            + "emailAddress=dek29@cornell.edu\n"
            + "archiveDirectory=/SSD/prop_pipeline/arcvtmp/\n"
            + "haplosDirectory=/SSD/haplos/\n"
            + "inputFolder=/workdir/tassel/tassel4-src/20130716test/raw_seq\n"
            + "enzyme=ApeKI\n"
            + "topmFile=/workdir/tassel/tassel4-src/20130716test/topm/AllZeaGBSv2.6ProdTOPM_20130605.topm.h5\n"
            + "outputFolder=/workdir/tassel/tassel4-src/20130716test/hap_maps\n"
            + "keyFile=/workdir/tassel/tassel4-src/20130716test/keyfile/MGP1_low_vol_2smallReps_key.txt";

    public ProductionPipeline(File runFile) {

        java.util.Properties props = new java.util.Properties();
        props.setProperty("log4j.logger.net.maizegenetics", "INFO, stdout");
        props.setProperty("log4j.appender.stdout",
                "org.apache.log4j.ConsoleAppender");
        props.setProperty("log4j.appender.stdout.layout",
                "org.apache.log4j.TTCCLayout");
        PropertyConfigurator.configure(props);

        try {
            myApplicationHost = InetAddress.getLocalHost().getHostName();
        } catch (UnknownHostException uhe) {
            // do nothing
        }

        executeRunFile(runFile);
    }

    private String getEmailSubjectRun(String runFile) {
        String result = myEmailSubjects.get(runFile);
        if (result != null) {
            return result;
        }

        StringBuilder builder = new StringBuilder();
        builder.append("Discussion: ");
        builder.append(getTimeStamp());
        builder.append("Run File: ").append(runFile).append(" ");
        result = builder.toString();
        myEmailSubjects.put(runFile, result);
        return result;
    }

    private String getEmailSubjectApp() {
        StringBuilder builder = new StringBuilder();
        builder.append("Discussion: ");
        builder.append(getTimeStamp());
        builder.append("Tassel GBS Production Pipeline on Host: ").append(myApplicationHost);
        return builder.toString();
    }

    /**
     *
     * @param runDirectoryIn Directory containing any number of .run files
     */
    private void executeRunFile(File runFile) {

        String currentRunFile = runFile.getAbsolutePath();
        String msgBody = "Starting to run " + currentRunFile + " on server " + myApplicationHost;
        sendAlertNotification(getEmailSubjectRun(currentRunFile), msgBody);
        String runFileContents = loadRunConfiguration(runFile);
        String todayDate = new SimpleDateFormat("yyyyMMdd").format(new Date());

        String fileName = FilenameUtils.removeExtension(runFile.getName());
        String fileNameBase = todayDate + "_" + fileName;
        String logFileName = fileNameBase + ".log";

        File logFile = new File(myOutputFolder + "/" + logFileName);
        try {
            if (!logFile.exists()) {
                logFile.createNewFile();
            }
            BufferedWriter bw = Utils.getBufferedWriter(logFile);
            bw.write("Contents of the .properties file:\n" + myPropertiesFileContents);
            bw.write(getTimeStamp() + "Contents of the .run file: " + "\n" + runFileContents);
            bw.write(getCurrentRunContext());
            bw.close();
        } catch (IOException ioe) {
            ioe.printStackTrace();
        }

        // redirect System.out and System.err to the log file
        PrintStream ps = null;
        try {
            ps = new PrintStream(new BufferedOutputStream(new FileOutputStream(logFile, true)));
        } catch (FileNotFoundException fnfe) {
            fnfe.printStackTrace();
        }
        System.setOut(ps);
        System.setErr(ps);

        System.out.println(getTimeStamp() + "Initializing ProductionSNPCallerPlugin \n");
        long start = System.nanoTime();

        String[] pluginArgs = getPipelinePluginArgs();
        StringBuilder builder = new StringBuilder();
        for (String s : pluginArgs) {
            builder.append(s + "\n");
        }
        System.out.println("Arguments passed to ProductionSNPCallerPlugin:\n" + builder.toString());
        ProductionSNPCallerPlugin pscp = new ProductionSNPCallerPlugin();
        System.out.println(getTimeStamp() + "Initialized ProductionSNPCallerPlugin \n");
        pscp.setParameters(pluginArgs);
        System.out.println(getTimeStamp() + "Done with ProductionSNPCallerPlugin.setParameters() \n");
        pscp.performFunction(null);
        System.out.println(getTimeStamp() + "Done with ProductionSNPCallerPlugin.performFunction() \n");

        double elapsedSeconds = (double) (System.nanoTime() - start) / 1_000_000_000.0;
        System.out.println(getTimeStamp() + "Time to run ProductionSNPCallerPlugin: " + elapsedSeconds + " sec.");

        if (myRunImputation) {
            start = System.nanoTime();
            String[] name = runFile.getName().split("\\.");
            String h5File = myOutputFolder + "/" + name[0] + ".hmp.h5";
            String haploDir = myHaplotypeDirectory + "/" + "AllZeaGBSv27.gX.hmp.txt.gz";
            String targetFile = myOutputFolder + "/" + name[0] + ".globalimp.hmp.h5";

            runImputation(h5File, haploDir, targetFile);
            elapsedSeconds = (double) (System.nanoTime() - start) / 1_000_000_000.0;
            System.out.println(getTimeStamp() + "Time to run Imputation: " + elapsedSeconds + " sec.");
        }

        String email = "Ran:\n " + myInputFolder
                + "\n\n  Tassel Pipeline Execution Time: " + elapsedSeconds + " seconds"
                + "\n\n Attachment:\n " + logFile.getAbsolutePath()
                + "\nRun on server: " + myApplicationHost;

        StringBuilder emailMsg = new StringBuilder(email);

        File toFile = new File(myArchiveDirectory + "/" + runFile.getName());

        boolean movedFile = false;
        try {
            Files.move(runFile, toFile);
            movedFile = true;
        } catch (IOException ioe) {
        }

        if (movedFile) {
            System.out.println("Moved file " + runFile.getAbsolutePath() + " to " + toFile.getAbsolutePath());
        } else {
            String msg = "******* COULD NOT MOVE FILE " + runFile.getAbsolutePath() + " TO " + toFile.getAbsolutePath()
                    + " on server: " + myApplicationHost;
            System.out.println(msg);
            sendAlertNotification(getEmailSubjectRun(currentRunFile), msg);
        }

        // send email notification that a .run file has been processed
        SMTPClient sc = new SMTPClient(myEmailHost, myRecipientEmailAddresses);

        try {
            sc.sendMessageWithAttachment(getEmailSubjectRun(currentRunFile), emailMsg.toString(), logFile.getAbsolutePath());
        } catch (javax.mail.MessagingException me) {
            // do nothing
        }
    }

    /**
     * @return Arguments to run ProductionSNPCallerPlugin
     */
    private String[] getPipelinePluginArgs() {
        String[] args = {
            "-i", myInputFolder,
            "-k", myKeyFile,
            "-e", myEnzyme,
            "-o", myOutputFolder,
            "-m", myTopmFile
        };
        return args;
    }

    private String getPipelinePluginArgsString() {
        StringBuilder builder = new StringBuilder();
        builder.append("-i ").append(myInputFolder);
        builder.append(" -k ").append(myKeyFile);
        builder.append(" -e ").append(myEnzyme);
        builder.append(" -o ").append(myOutputFolder);
        builder.append(" -m ").append(myTopmFile);
        return builder.toString();
    }

    /**
     * Load a file containing the information necessary to write out an XML
     * output file necessary for running the production pipeline
     *
     * @param aFileIn
     * @return
     */
    private String loadRunConfiguration(File aFileIn) {

        Properties props = new Properties();
        try {
            System.out.println(aFileIn.getAbsoluteFile());
            props.load(new FileInputStream(aFileIn));
        } catch (IOException ioe) {
            System.out.println("Problem loading run file configuration file:" + aFileIn);
            System.out.println("************** Example .properties file: ");
            System.out.println(EXAMPLE_RUN_FILE);
            ioe.printStackTrace();
            sendAlertNotification(getEmailSubjectApp(), "Properties file could not be loaded: "
                    + aFileIn + " on server " + myApplicationHost);
            System.exit(1);
        }

        // server used for sending email
        String configurationElement = "emailHost";
        myEmailHost = props.getProperty(configurationElement, myEmailHost);

        // to whom the email notifications should be sent
        configurationElement = "emailAddress";
        String address = props.getProperty(configurationElement);
        if (address != null) {
            Iterable<String> results = Splitter.on(EMAIL_ADDRESS_DELIMITER).split(address);
            myRecipientEmailAddresses = Iterables.toArray(results, String.class);
        }

        // directory into which the key and run files are placed after being run
        configurationElement = "archiveDirectory";
        myArchiveDirectory = props.getProperty(configurationElement);
        String response = testInputDirectory(aFileIn, myArchiveDirectory, configurationElement);
        if (response != null) {
            System.out.println(response);
        }

        // directory into which the key and run files are placed after being run
        configurationElement = "haplosDirectory";
        myHaplotypeDirectory = props.getProperty(configurationElement);
        response = testInputDirectory(aFileIn, myHaplotypeDirectory, configurationElement);
        if (response != null) {
            System.out.println(response);
        }

        configurationElement = "inputFolder";
        myInputFolder = props.getProperty(configurationElement);
        response = testInputDirectory(aFileIn, myInputFolder, configurationElement);
        if (response != null) {
            System.out.println(response);
        }

        configurationElement = "enzyme";
        myEnzyme = props.getProperty(configurationElement);
        response = testInputDirectory(aFileIn, myInputFolder, configurationElement);
        if (response != null) {
            System.out.println(response);
        }

        configurationElement = "topmFile";
        myTopmFile = props.getProperty(configurationElement);
        response = testInputDirectory(aFileIn, myInputFolder, configurationElement);
        if (response != null) {
            System.out.println(response);
        }

        configurationElement = "outputFolder";
        myOutputFolder = props.getProperty(configurationElement);
        response = testInputDirectory(aFileIn, myInputFolder, configurationElement);
        if (response != null) {
            System.out.println(response);
        }

        configurationElement = "keyFile";
        myKeyFile = props.getProperty(configurationElement);
        response = testInputDirectory(aFileIn, myInputFolder, configurationElement);
        if (response != null) {
            System.out.println(response);
        }

        BufferedReader br = null;
        StringBuilder sb = new StringBuilder();
        try {
            br = new BufferedReader(new FileReader(aFileIn));
            String line = null;
            while ((line = br.readLine()) != null) {
                sb.append(line + "\n");
            }
        } catch (IOException ioe) {

        }
        return sb.toString();
    }

    /**
     * Verify that element specifying a directory exists in a configuration or
     * properties file and that the directory exists on the file system.
     *
     * @param filename Name of configuration/properties/.run file
     * @param input String containing the fully qualified path of a directory
     * @param configurationElement
     *
     * @return Description of the problem with the input. If return is null then
     * there is no problem with the entered file or directory
     */
    private String testInputDirectory(File filename, String input, String configurationElement) {

        String response = null;
        if (input == null) {
            response = filename.getAbsolutePath() + " is missing a run configuration element:  " + configurationElement;
        } else {
            File aFileOrDir = new File(input);
            if (!aFileOrDir.exists()) {
                response = filename.getAbsolutePath() + "'s configuration element " + configurationElement + " does not exist.  Please confirm path and filename.";
            }
        }
        return response;
    }

    /**
     * Collect information about the context in which the pipeline is being run.
     * This information includes the following: 1) Current date and time. 2)
     * Name of the machine on which run occurs. 3) User account used 4) Files
     * used and their MD5sums 5) Contents of XML configuration file used to run
     * TasselPipeline
     *
     * @return Information about current run
     */
    private String getCurrentRunContext() {

        StringBuilder sb = new StringBuilder();

        // date
        sb.append(getTimeStamp()).append("\n");

        // user account name
        String user = System.getProperty("user.name");
        sb.append("User Account Name: ").append(user).append("\n");

        // hostname
        sb.append("Name of Machine on which JVM is Running: ").append(myApplicationHost).append("\n");

        // for each file in a directory, include the md5sum
        sb.append(getTimeStamp()).append("MD5: ").append(myKeyFile).append(": ").append(CheckSum.getMD5Checksum(myKeyFile)).append("\n");
        File inFldr = new File(myInputFolder);

        if (inFldr.isDirectory()) {
            File[] files = inFldr.listFiles();
            for (File f : files) {
                sb.append(getTimeStamp()).append("MD5: ").append(f.getPath()).append(": ").append(CheckSum.getMD5Checksum(f.getPath())).append("\n");
            }
        } else {
            sb.append(getTimeStamp()).append(CheckSum.getMD5Checksum(myInputFolder)).append("\n");
        }

        return sb.toString();
    }

    /**
     *
     * @param unImpTargetFile
     * @param donorFile "AllZeaGBSv27.gX.hmp.txt.gz"
     * @param impTargetFile
     * @return
     */
    private void runImputation(String unImpTargetFile, String donorFile, String impTargetFile) {
        String[] args2 = new String[]{
            "-hmp", unImpTargetFile,
            "-d", donorFile,
            "-o", impTargetFile,
            "-minMnCnt", "20",
            "-mxInbErr", "0.02",
            "-mxHybErr", "0.005",
            "-mnTestSite", "50",
            "-mxDonH", "10", // "-projA",
        };

        StringBuilder builder = new StringBuilder();
        for (String s : args2) {
            builder.append(s).append("\n");
        }
        System.out.println("Arguments passed to MinorWindowViterbiImputationPlugin:\n" + builder.toString());
        TasselPrefs.putAlignmentRetainRareAlleles(false);
        FILLINImputationPlugin plugin = new FILLINImputationPlugin();
        plugin.setParameters(args2);
        plugin.performFunction(null);
    }

    /**
     * Generates basic summary information about imputation such as change in
     * sites or taxa count change in missing data change in heterozygosity how
     * many sites changed their major allele
     *
     * @param originalFile
     * @param imputedFile
     * @return
     */
    public static String compareOriginalAgainstImputed(String originalFile, String imputedFile) {

        StringBuilder sb = new StringBuilder();
        GenotypeTable origAlignment = ImportUtils.readGuessFormat(originalFile);
        GenotypeTable impAlignment = ImportUtils.readGuessFormat(imputedFile);

        int siteCount = origAlignment.numberOfSites();
        int taxaCount = origAlignment.numberOfTaxa();
        int totalSiteCount = siteCount * taxaCount;

        int siteDelta = Math.abs(origAlignment.numberOfSites() - impAlignment.numberOfSites());
        int taxaDelta = Math.abs(origAlignment.numberOfTaxa() - impAlignment.numberOfTaxa());

        int origTotalSitesNotMissing = 0, impTotalSitesNotMissing = 0, totalSitesNotMissingDelta = 0;
        double origProportionNotMissing = 0, impProportionNotMissing = 0, proportionNotMissingDelta = 0;

        int origHetCount = 0, impHetCount = 0, hetCountDelta = 0;
        double origHetProportion = 0, impHetProportion = 0, hetProportionDelta = 0;
        int flipCount = 0;     // number of sites that have had a change in which allele is major

        int allelicChangeCount = 0;
        if (siteDelta == 0) {
            for (int i = 0; i < siteCount; i++) {

                // heterozygous counts
                origHetCount += origAlignment.heterozygousCount(i);
                impHetCount += impAlignment.heterozygousCount(i);
                hetCountDelta = impHetCount - origHetCount;

                //not missing
                origTotalSitesNotMissing += origAlignment.totalNonMissingForSite(i);
                impTotalSitesNotMissing += impAlignment.totalNonMissingForSite(i);

                // switching of major and minor allele
                byte origMajorAllele = origAlignment.majorAllele(i);
                byte impMajorAllele = impAlignment.majorAllele(i);
                if (origMajorAllele != impMajorAllele) {
                    flipCount++;
                    double diff = Math.abs(origAlignment.majorAlleleFrequency(i) - impAlignment.minorAlleleFrequency(i));
                    allelicChangeCount += (int) diff * taxaCount;
                } else {
                    double diff = Math.abs(origAlignment.majorAlleleFrequency(i) - impAlignment.majorAlleleFrequency(i));
                    allelicChangeCount += (int) diff * taxaCount;
                }
            }

            totalSitesNotMissingDelta = impTotalSitesNotMissing - origTotalSitesNotMissing;
            origProportionNotMissing = (double) origTotalSitesNotMissing / (double) totalSiteCount;
            impProportionNotMissing = (double) impTotalSitesNotMissing / (double) totalSiteCount;
            proportionNotMissingDelta = impProportionNotMissing - origProportionNotMissing;

            hetCountDelta = impHetCount - origHetCount;
            origHetProportion = (double) origHetCount / (double) totalSiteCount;
            impHetProportion = (double) impHetCount / (double) totalSiteCount;
            hetProportionDelta = impHetProportion - origHetProportion;

        }

        sb.append("\nSites: " + siteCount + "\tSite Delta: " + siteDelta);
        sb.append("\nTaxa: " + taxaCount + "\tTaxa Delta: " + taxaDelta);
        sb.append("\nTotal Sites: " + totalSiteCount);
        sb.append("\nSites Not Missing Original: " + origTotalSitesNotMissing
                + "\tSites Not Missing Imputed: " + impTotalSitesNotMissing
                + "\tSites Not Missing Delta: " + totalSitesNotMissingDelta);
        sb.append("\nProportion Not Missing Original: " + origProportionNotMissing
                + "\tProportion Not Missing Imputed: " + impProportionNotMissing
                + "\tProportion Not Missing Delta: " + proportionNotMissingDelta);
        sb.append("\nChange in Heterozygous Sites: " + hetCountDelta);
        sb.append("\nHeterozygous Sites Original: " + origHetCount
                + "\tHeterozygous Sites Imputed: " + impHetCount
                + "\tHet Delta: " + hetCountDelta);
        sb.append("\nHeterozygous Proportion Original: " + origHetProportion
                + "\tHeterozygous Proportion Imputed: " + impHetProportion
                + "\tHet Proportion Delta: " + hetProportionDelta);
        sb.append("\nTotal Alleles Changed: " + allelicChangeCount
                + "\tProportion of Alleles Changed: " + (double) allelicChangeCount / totalSiteCount);
        sb.append("\nNumber of Sites Changing Major Allele: " + flipCount
                + "\tMajor <-> Minor Proportion: " + (double) flipCount / (double) totalSiteCount);

        return sb.toString();
    }

    /**
     * Convenience method to provide uniformly labeled timestamps
     */
    private static String getTimeStamp() {
        return "Timestamp: " + LOGGING_DATE_FORMAT.format(new Date()) + ": ";
    }

    /**
     * Send email to pipeline administrator when issue is encountered that may
     * keep pipeline from successfully completing, e.g., cannot load .properties
     * file, or send notification of how many times the pipeline should run,
     * i.e., how many .run files are present.
     *
     * @param subject Subject line of email
     * @param message Message body of email
     */
    private void sendAlertNotification(String subject, String message) {
        System.out.println("myEmailHost: " + myEmailHost);
        System.out.println("myRecipientEmailAddresses: " + myRecipientEmailAddresses);
        SMTPClient sc = new SMTPClient(myEmailHost, myRecipientEmailAddresses);
        try {
            sc.sendMessage(subject, message);
        } catch (javax.mail.MessagingException me) {
            // do nothing
        }
    }

    public static void main(String[] args) {

        if (args.length != 1) {
            System.out.println("Usage: ProductionPipeline <run directory>");
            System.exit(1);
        }

        File inputDirectory = new File(args[0]);
        File[] runFiles = inputDirectory.listFiles(new FilenameFilter() {
            @Override
            public boolean accept(File dir, String name) {
                return name.toLowerCase().endsWith(RUN_FILE_SUFFIX);
            }
        });

        if ((runFiles == null) || (runFiles.length == 0)) {
            System.out.println("ProductionPipeline: Could not find a valid .run files in directory: " + args[0]);
            System.out.println("ProductionPipeline: Example .run file: ");
            System.out.println(EXAMPLE_RUN_FILE);
        }

        for (File current : runFiles) {
            System.out.println("ProductionPipeline: current run file: " + current.getAbsolutePath());
            new ProductionPipeline(current);
        }

    }
}
