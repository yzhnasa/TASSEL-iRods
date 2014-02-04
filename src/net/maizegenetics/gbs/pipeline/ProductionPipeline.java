/*
 * ProductionPipeline
 */
package net.maizegenetics.gbs.pipeline;

import com.google.common.base.Splitter;
import com.google.common.collect.Iterables;
import com.google.common.io.Files;

import net.maizegenetics.analysis.imputation.MinorWindowViterbiImputationPlugin;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.prefs.TasselPrefs;
import net.maizegenetics.util.CheckSum;
import net.maizegenetics.util.SMTPClient;

import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang.StringUtils;

import org.apache.log4j.Logger;

import java.io.*;

import java.net.InetAddress;
import java.net.UnknownHostException;

import java.text.SimpleDateFormat;

import java.util.Date;
import java.util.Properties;

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

    private static final String DEFAULT_APPLICATION_CONFIGURATION = "production_pipeline.properties";
    private static final String EMAIL_ADDRESS_DELIMITER = ";";
    // standard subject line GPP = GBS Production Pipeline
    private static final String EMAIL_SUBJECT_BASE = "GPP ";

    private static final SimpleDateFormat LOGGING_DATE_FORMAT = new SimpleDateFormat("yyyyMMdd HH:mm:ss");

    // host on which pipeline is running
    private String myApplicationHost = "unknown";

    // run file suffix
    private String myRunFileSuffix = ".run";

    // default server to be used to send email notifications
    private String myEmailHost = "appsmtp.mail.cornell.edu";

    // pipeline admin default email address
    private String[] myRecipientEmailAddresses = {"tmc46@cornell.edu"};

    // directory into which the key and run files are
    // placed for pickup by cron job
    private String myRunDirectory = "/SSD/prop_pipeline/run/";

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

    private String expectedCheckSum = "10e75e612ade7979f210958933da4de9";

    private String myPropertiesFileContents = null;

    private static final String EXAMPLE_CONFIG_FILE
            = "runFileSuffix=.run\n"
            + "emailHost=appsmtp.mail.cornell.edu\n"
            + "emailAddress=dek29@cornell.edu\n"
            + "runDirectory=/SSD/prop_pipeline/run/\n"
            + "archiveDirectory=/SSD/prop_pipeline/arcvtmp/\n"
            + "haplosDirectory=/SSD/haplos/\n";

    private static final String EXAMPLE_RUN_FILE
            = "inputFolder=/workdir/tassel/tassel4-src/20130716test/raw_seq\n"
            + "enzyme=ApeKI\n"
            + "topmFile=/workdir/tassel/tassel4-src/20130716test/topm/AllZeaGBSv2.6ProdTOPM_20130605.topm.h5\n"
            + "outputFolder=/workdir/tassel/tassel4-src/20130716test/hap_maps\n"
            + "keyFile=/workdir/tassel/tassel4-src/20130716test/keyfile/MGP1_low_vol_2smallReps_key.txt";

    public ProductionPipeline(String appPropertiesFile, boolean runCheckSum, boolean runImputation, String testingCheckSum) {

        try {
            myApplicationHost = InetAddress.getLocalHost().getHostName();
        } catch (UnknownHostException uhe) {
            // do nothing
        }

        // if no application property file is specified, try looking for something in the home directory
        if (appPropertiesFile == null) {
            appPropertiesFile = DEFAULT_APPLICATION_CONFIGURATION;
        }
        myPropertiesFileContents = loadApplicationConfiguration(appPropertiesFile);

        loadRunFiles(myRunDirectory, runCheckSum, runImputation, testingCheckSum);
    }

    /**
     *
     * @param runDirectoryIn Directory containing any number of .run files
     */
    private void loadRunFiles(String runDirectoryIn, boolean doCheckSum, boolean runImputation, String testingCheckSum) {

        File dir = new File(runDirectoryIn);

        if (!dir.exists()) {
            System.out.println("Could not find the directory containing .run files: " + dir.getPath());
            System.out.println("Exiting program.");
            sendAlertNotification(EMAIL_SUBJECT_BASE + "- Error", "Could not find directory: " + dir.getAbsolutePath()
                    + " on  server " + myApplicationHost);
            System.exit(1);
        }

        // get all property files in the directory
        File[] files = dir.listFiles(new FilenameFilter() {
            public boolean accept(File dir, String name) {
                return name.toLowerCase().endsWith(myRunFileSuffix);
            }
        });

        if (files == null) {
            System.out.println("************** Could not find a valid .run file ***************");
            System.out.println("************** Example .run file: ");
            System.out.println(EXAMPLE_RUN_FILE);   // display properly configured run file
            sendAlertNotification(EMAIL_SUBJECT_BASE + "- No Files", "No .run files found on " + myApplicationHost);
        } else {
            StringBuilder sb = new StringBuilder();
            for (File f : files) {
                sb.append(f + "\n");
            }
            sb.append("\nRunning on server: " + myApplicationHost + "\n");
            sendAlertNotification(EMAIL_SUBJECT_BASE + "- File Count: " + files.length, sb.toString());
        }

        for (File aFile : files) {

            String msgBody = "Starting to run " + aFile.getAbsolutePath() + " on server " + myApplicationHost;
            sendAlertNotification(EMAIL_SUBJECT_BASE + " File: " + aFile.getName(), msgBody);
            String runFileContents = loadRunConfiguration(aFile);
            SimpleDateFormat yyyyMMdd_format = new SimpleDateFormat("yyyyMMdd");
            String todayDate = yyyyMMdd_format.format(new Date());

            String fileName = FilenameUtils.removeExtension(aFile.getName());
            String fileNameBase = todayDate + "_" + fileName;
            String logFileName = fileNameBase + ".log";

            String contextLog = recordContext(doCheckSum);    // get file MD5sum

            String runFileMsg = getTimeStamp() + ": Contents of the .run file: ";

            File logFile = new File(myOutputFolder + "/" + logFileName);
            try {
                if (!logFile.exists()) {
                    logFile.createNewFile();
                }
                BufferedWriter bw = new BufferedWriter(new FileWriter(logFile.getAbsolutePath()));
                bw.write("Contents of the .properties file:\n" + myPropertiesFileContents);
                bw.write(runFileMsg + "\n" + runFileContents);
                bw.write(contextLog);
                bw.close();

            } catch (IOException ioe) {
                ioe.printStackTrace();
            }

            // redirect System.out and System.err to the log file
            PrintStream ps = null;
            try {
                ps = new PrintStream(new BufferedOutputStream(new FileOutputStream(logFile.getAbsolutePath(), true)));

            } catch (FileNotFoundException fnfe) {
                fnfe.printStackTrace();
            }
            System.setOut(ps);
            System.setErr(ps);

            System.out.println(getTimeStamp() + " Initializing ProductionSNPCallerPlugin \n");
            Date start = new Date();

            String[] pluginArgs = getPipelinePluginArgs();
            StringBuilder builder = new StringBuilder();
            for (String s : pluginArgs) {
                builder.append(s + "\n");
            }
            System.out.println("Arguments passed to ProductionSNPCallerPlugin:\n" + builder.toString());
            ProductionSNPCallerPlugin pscp = new ProductionSNPCallerPlugin();
            System.out.println(getTimeStamp() + " Initialized ProductionSNPCallerPlugin \n");
            pscp.setParameters(pluginArgs);
            System.out.println(getTimeStamp() + " Done with ProductionSNPCallerPlugin.setParameters() \n");
            pscp.performFunction(null);
            System.out.println(getTimeStamp() + " Done with ProductionSNPCallerPlugin.performFunction() \n");

            if (runImputation) {
                String[] name = aFile.getName().split("\\.");
                String h5File = myOutputFolder + "/" + name[0] + ".hmp.h5";
                String haploDir = myHaplotypeDirectory + "/" + "AllZeaGBSv27.gX.hmp.txt.gz";
                String targetFile = myOutputFolder + "/" + name[0] + ".globalimp.hmp.h5";

                runImputation(h5File, haploDir, targetFile);
            }
            Date stop = new Date();

            long startTime = start.getTime();
            long stopTime = stop.getTime();
            long diff = stopTime - startTime;
            long elapsedSeconds = diff / 1000;

            String emailSubject = EMAIL_SUBJECT_BASE + myInputFolder;
            String email = "Ran:\n " + myInputFolder
                    + "\n\n  Tassel Pipeline Execution Time: " + elapsedSeconds + " seconds"
                    + "\n\n Attachment:\n " + logFile.getAbsolutePath()
                    + "\nRun on server: " + myApplicationHost;

            StringBuilder emailMsg = new StringBuilder(email);

            // test section: get checksum on hapmap file(s), notify whether matching or not
            if (testingCheckSum != null) {

                // find hapmap file in output directory
                final String suffix = ".hmp.txt.gz";
                File outDir = new File(myOutputFolder);
                File[] hapMapFiles = outDir.listFiles(new FilenameFilter() {
                    @Override
                    public boolean accept(File dir, String name) {
                        return name.endsWith(suffix);
                    }
                });
                boolean passedTest = false;
                for (File f : hapMapFiles) {
                    // get checksum
                    String filename = f.getAbsolutePath();
                    String cksum = CheckSum.getChecksum(filename, "MD5");
                    if (cksum.equalsIgnoreCase(expectedCheckSum)) {

                        emailSubject = EMAIL_SUBJECT_BASE + "test passed";
                        passedTest = true;
                    }
                    emailMsg.append("\nFile: " + filename + "\tChecksum:" + cksum);
                }

                if (!passedTest) {
                    emailSubject = EMAIL_SUBJECT_BASE + "TEST FAILED!";
                }
            } else {

                emailSubject = EMAIL_SUBJECT_BASE + myInputFolder;
                File toFile = new File(myArchiveDirectory + "/" + aFile.getName());

                boolean movedFile = false;
                try {
                    Files.move(aFile, toFile);
                    movedFile = true;
                } catch (IOException ioe) {
                }

                if (movedFile) {
                    System.out.println("Moved file " + aFile.getAbsolutePath() + " to " + toFile.getAbsolutePath());
                } else {
                    String msg = "******* COULD NOT MOVE FILE " + aFile.getAbsolutePath() + " TO " + toFile.getAbsolutePath()
                            + " on server: " + myApplicationHost;
                    System.out.println(msg);
                    sendAlertNotification(EMAIL_SUBJECT_BASE + "- Error", msg);
                }
            }
            // send email notification that a .run file has been processed
            SMTPClient sc = new SMTPClient(myEmailHost, myRecipientEmailAddresses);

            try {
                sc.sendMessageWithAttachment(emailSubject, emailMsg.toString(), logFile.getAbsolutePath());
            } catch (javax.mail.MessagingException me) {   /* ignore */ }
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

    /**
     * Load application-wide properties file and initialize variables.
     *
     * @param aFileIn config file
     *
     * @return Contents of config file
     */
    private String loadApplicationConfiguration(String aFileIn) {

        Properties props = new Properties();
        try {
            File propsFile = new File(aFileIn);
            System.out.println(propsFile.getAbsoluteFile());

            props.load(new FileInputStream(aFileIn));
        } catch (IOException ioe) {
            System.out.println("Problem loading application configuration file:" + aFileIn);
            System.out.println("************** Example .properties file: ");
            System.out.println(EXAMPLE_CONFIG_FILE);
            ioe.printStackTrace();
            sendAlertNotification(EMAIL_SUBJECT_BASE + "- Error", "Properties file could not be loaded: "
                    + aFileIn + " on server " + myApplicationHost);
            System.exit(1);
        }

        String configurationElement = "runFileSuffix";
        myRunFileSuffix = props.getProperty(configurationElement, myRunFileSuffix);

        // server used for sending email
        configurationElement = "emailHost";
        myEmailHost = props.getProperty(configurationElement, myEmailHost);

        // to whom the email notifications should be sent
        configurationElement = "emailAddress";
        String address = props.getProperty(configurationElement);
        if (address != null) {
            Iterable<String> results = Splitter.on(EMAIL_ADDRESS_DELIMITER).split(address);
            myRecipientEmailAddresses = Iterables.toArray(results, String.class);
        }

        // directory where the .run files are expected to be
        configurationElement = "runDirectory";
        myRunDirectory = props.getProperty(configurationElement);
        String response = testInputDirectory(aFileIn, myRunDirectory, configurationElement);
        if (response != null) {
            System.out.println(response);
        }

        // directory into which the key and run files are placed after being run
        configurationElement = "archiveDirectory";
        myArchiveDirectory = props.getProperty(configurationElement);
        response = testInputDirectory(aFileIn, myArchiveDirectory, configurationElement);
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

        // read in the contents of the properties file so it can be placed into the log
        BufferedReader br = null;
        StringBuilder sb = new StringBuilder();
        try {
            br = new BufferedReader(new FileReader(aFileIn));
            String line = null;
            while ((line = br.readLine()) != null) {
                sb.append(line + "\n");
            }
        } catch (IOException ioe) {

        } finally {
            try {
                br.close();
            } catch (Exception e) {
                // do nothing
            }
        }
        return sb.toString();
    }

    /**
     * Load a file containing the information necessary to write out an XML
     * output file necessary for running the production pipeline
     *
     * @param aFileIn
     * @return
     */
    private String loadRunConfiguration(File aFileIn) {

        String usage = aFileIn.getName() + " is missing a run configuration element:  ";
        Properties props = new Properties();
        try {
            props.load(new FileInputStream(aFileIn));
        } catch (IOException ioe) {
            System.err.println("Issue loading run configuration file: " + aFileIn.getName());
            ioe.printStackTrace();
        }
        String configurationElement = "inputFolder";
        myInputFolder = props.getProperty(configurationElement);
        String response = testInputDirectory(aFileIn.getName(), myInputFolder, configurationElement);
        if (response != null) {
            System.out.println(response);
        }

        configurationElement = "enzyme";
        myEnzyme = props.getProperty(configurationElement);
        response = testInputDirectory(aFileIn.getName(), myInputFolder, configurationElement);
        if (response != null) {
            System.out.println(response);
        }

        configurationElement = "topmFile";
        myTopmFile = props.getProperty(configurationElement);
        response = testInputDirectory(aFileIn.getName(), myInputFolder, configurationElement);
        if (response != null) {
            System.out.println(response);
        }

        configurationElement = "outputFolder";
        myOutputFolder = props.getProperty(configurationElement);
        response = testInputDirectory(aFileIn.getName(), myInputFolder, configurationElement);
        if (response != null) {
            System.out.println(response);
        }

        configurationElement = "keyFile";
        myKeyFile = props.getProperty(configurationElement);
        response = testInputDirectory(aFileIn.getName(), myInputFolder, configurationElement);
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
     * properties file and that the directory exists on the filesystem.
     *
     * @param filename Name of configuration/properties/.run file
     * @param input String containing the fully qualified path of a directory
     * @param configurationElement
     * @return Description of the problem with the input. If return is null then
     * there is no problem with the entered file or directory
     */
    private String testInputDirectory(String filename, String input, String configurationElement) {

        String response = null;
        if (input == null) {
            response = filename + " is missing a run configuration element:  " + configurationElement;
        } else {
            File aFileOrDir = new File(input);
            if (!aFileOrDir.exists()) {
                response = filename + "'s configuration element " + configurationElement + " does not exist.  Please confirm path and filename.";
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
     * @param calculateChecksum Calculating MD5sum can be time consuming. This
     * allows it to be skipped when appropriate.
     * @return
     */
    private String recordContext(boolean calculateChecksum) {

        StringBuilder sb = new StringBuilder();

        // date
        sb.append(getTimeStamp()).append("\n");

        // user account name
        String user = System.getProperty("user.name");
        sb.append("User Account Name: ").append(user).append("\n");

        // hostname
        sb.append("Name of Machine on which JVM is Running: ").append(myApplicationHost).append("\n");

        String checkSumString = " md5sum: ";
        // for each file in a directory, include the md5sum
        if (calculateChecksum) {
            sb.append(getTimeStamp() + checkSumString + myKeyFile + " " + CheckSum.getMD5Checksum(myKeyFile) + "\n");
            File inFldr = new File(myInputFolder);

            if (inFldr.isDirectory()) {
                File[] files = inFldr.listFiles();
                for (File f : files) {
                    String fCheckSum = CheckSum.getMD5Checksum(f.getPath());
                    String msg = getTimeStamp() + checkSumString + f.getPath() + " " + fCheckSum + "\n";
                    sb.append(msg);
                }
            } else {
                sb.append(getTimeStamp() + CheckSum.getMD5Checksum(myInputFolder) + "\n");
            }

        } else {
            sb.append(getTimeStamp() + "MD5sum checking has been switched off using the --skipCheckSum argument");
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
            "-mxDonH", "10", //           "-projA",
        };

        StringBuilder builder = new StringBuilder();
        for (String s : args2) {
            builder.append(s).append("\n");
        }
        System.out.println("Arguments passed to MinorWindowViterbiImputationPlugin:\n" + builder.toString());
        System.out.println("TasselPrefs: " + TasselPrefs.getAlignmentRetainRareAlleles());
        TasselPrefs.putAlignmentRetainRareAlleles(false);
        MinorWindowViterbiImputationPlugin plugin = new MinorWindowViterbiImputationPlugin();
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
        return "Timestamp: " + LOGGING_DATE_FORMAT.format(new Date()) + " ";
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
        SMTPClient sc = new SMTPClient(myEmailHost, myRecipientEmailAddresses);
        try {
            sc.sendMessage(subject, message);
        } catch (javax.mail.MessagingException me) { /* ignore */        }
    }

    public static void main(String[] args) {

        String msg = "\n--skipCheckSum flag allows MD5sum checking to be skipped.\n"
                + "--skipImputation flag allows imputation to be skipped\n"
                + "--runTest flag is for use with a test data set and should be followed by the expected MD5Sum\n"
                + "--propsFile should be followed by a fully-qualified path name without spaces\n";

        boolean doCheckSum = true;
        boolean doImputation = true;
        String expectedChksm = null;
        String propsFile = "propsFile";
        String propsFilePath = null;
        String skipCheckSum = "skipCheckSum";
        String skipImputation = "skipImputation";
        String runTest = "runTest";
        if (args != null) {
            for (int i = 0; i < args.length; i++) {
                if (StringUtils.containsIgnoreCase(args[i], skipCheckSum)) {
                    doCheckSum = false;
                    System.out.println("Skipping Checksums");
                }
                if (StringUtils.containsIgnoreCase(args[i], skipImputation)) {
                    doImputation = false;
                    System.out.println("Skipping Imputation");
                }
                if (StringUtils.containsIgnoreCase(args[i], runTest)) {
                    if (args.length > i + 1) {
                        expectedChksm = args[i + 1];
                        i++;
                        System.out.println("Running testing and expecting the following MD5Sum: " + expectedChksm);
                    } else {
                        System.out.println("No checksum following --" + runTest + " flag.\n");
                        System.out.println(msg);
                    }
                }
                if (StringUtils.containsIgnoreCase(args[i], propsFile)) {
                    if (args.length > i + 1) {
                        propsFilePath = args[i + 1];
                        i++;
                        System.out.println("--" + propsFile + "\tUsing this properties file: " + propsFilePath);
                    } else {
                        System.out.println("No path following --" + propsFile + " flag.\n"
                                + "Will look for " + DEFAULT_APPLICATION_CONFIGURATION + " file in the current directory");
                        System.out.println(msg);
                    }
                }
            }
        } else {
            System.out.println(msg);
        }

        new ProductionPipeline(propsFilePath, doCheckSum, doImputation, expectedChksm);
    }
}
