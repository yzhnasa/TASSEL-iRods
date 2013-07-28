/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.gbs.pipeline;

import java.io.*;
import java.net.UnknownHostException;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.net.InetAddress;
import java.util.Properties;
import java.util.logging.Level;

import com.google.common.collect.Iterables;
import com.google.common.io.Files;
import net.maizegenetics.pipeline.TasselPipeline;
import net.maizegenetics.pipeline.TasselPipelineXMLUtil;
import net.maizegenetics.util.CheckSum;
import net.maizegenetics.util.SMTPClient;
import com.google.common.base.Splitter;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang.StringUtils;
import org.apache.log4j.Logger;


/**
 * This class is for running the GBS Production Pipeline.  It is to be run from within the sTASSEL.jar.  The cron job
 * should be set up to run the run_pipeline.pl which has been modified to make this class the main() to be run.  The
 * JVM memory settings within run_pipeline.pl should also be adjusted upwards.
 *
 * 20130718 Note from Jeff Glaubitz: A minor (and not yet documented) detail:  SeqToGenos needs the key file name to end with "_key.txt".
 * All the output files are named after the key file but replacing "_key.txt" with the appropriate extension.
 * User: dkroon
 * Date: 4/8/13
 */
public class ProductionPipelineMain {


    private final Logger myLogger = Logger.getLogger(ProductionPipelineMain.class);

//    myLogger.setLevel(Level.OFF);
    // these could reasonably be contents of a properties file
    private String applicationConfiguration = "production_pipeline.properties";

    private String runFileSuffix = ".run";                        // default run file suffix
    private String emailHost = "appsmtp.mail.cornell.edu";        // default server to be used to send email notifications
    private String[] emailAddresses = {"dek29@cornell.edu"};      // pipeline admin default email address
    private String emailAddressDelimiters = ";";                  // delimiter between email addresses
    private String emailSubjectBase = "GBS Production Pipeline "; // standard subject line

    private String runDirectory = "/workdir/prop_pipeline/run/";        // directory into which the key and run files are placed for pickup by cron job
    private String archiveDirectory = "/workdir/prop_pipeline/done/";   // directory into which the key and run files are placed after being run
                                                                        // a dated directory will be created within this directory and artifacts
                                                                        // such as the .log and .xml file will be placed there
    private boolean runCheckSum = true;             // switch for turning off checksum when appropriate

    private String todayDate = null;
    private String fileNameBase = null;     //  todayDate + property filename   to be used for XML and Log files
    private String anInputFolder= null;
    private String enzyme= null;
    private String topmFile= null;
    private String outputFolder= null;
    private String keyFile= null;
    private String hostName = "host unknown";

    private String[] tasselPipelineArgs;        // replacement for using an XML config file to pass relevant parameter to TasselPipeline

    private SimpleDateFormat dateFormat = new SimpleDateFormat("yyyyMMdd HH:mm:ss");      // make dateFormat uniform for all logging
    private String propertiesFileContents = null;

    //todo: log canonical paths and absolute paths of fastq files so as to resolve symlinks and document them

    private String exampleAppConfigFile =
            "runFileSuffix=.run\n" +
            "emailHost=appsmtp.mail.cornell.edu\n" +
            "emailAddress=dek29@cornell.edu\n" +
            "runDirectory=/workdir/prop_pipeline/run/\n" +
            "archiveDirectory=/workdir/prop_pipeline/arcv/\n";

    private String exampleRunFile =
            "inputFolder=/workdir/tassel/tassel4-src/20130716test/raw_seq\n" +
            "enzyme=ApeKI\n" +
            "topmFile=/workdir/tassel/tassel4-src/20130716test/topm/AllZeaGBSv2.6ProdTOPM_20130605.topm.h5\n" +
            "outputFolder=/workdir/tassel/tassel4-src/20130716test/hap_maps\n" +
            "keyFile=/workdir/tassel/tassel4-src/20130716test/keyfile/MGP1_low_vol_2smallReps_key.txt";


    public ProductionPipelineMain(boolean runCheckSumIn){
        runCheckSum = runCheckSumIn;
        init();
    }

    private void init(){

        propertiesFileContents = loadApplicationConfiguration(applicationConfiguration);

        loadRunFiles(runDirectory);
    }

    /**
     *
     * @param runDirectoryIn  Directory containing any number of .run files
     */
    private void loadRunFiles(String runDirectoryIn){

        File dir = new File(runDirectoryIn);

        if(!dir.exists()){
            System.out.println("Could not find the directory containing .run files: " + dir.getPath());
            System.out.println("Exiting program.");
            sendAlertNotification(emailSubjectBase + "- Error", "Could not find directory: " + dir.getAbsolutePath());
            System.exit(1);
        }

        // get all property files in the directory
        File[] files = dir.listFiles(new FilenameFilter() {
            public boolean accept(File dir, String name) {
                return name.toLowerCase().endsWith(runFileSuffix);
            }
        });

        if(files == null){
            System.out.println("************** Could not find a valid .run file ***************");
            System.out.println("************** Example .run file: ");
            System.out.println(exampleRunFile);   // display properly configured run file
            sendAlertNotification(emailSubjectBase + "- No Files", "No .run files found");
        }else{
            StringBuffer sb = new StringBuffer();
            for(File f: files){
                sb.append(f + "\n");
            }
            sendAlertNotification(emailSubjectBase + "- File Count: " + files.length, sb.toString());
        }

        for(File aFile: files){

            String runFileContents = loadRunConfiguration(aFile);
            SimpleDateFormat yyyyMMdd_format = new SimpleDateFormat("yyyyMMdd");
            todayDate = yyyyMMdd_format.format(new Date());

            String fileName = FilenameUtils.removeExtension(aFile.getName());
            fileNameBase =   todayDate + "_"  + fileName;
            String logFileName = fileNameBase + ".log";

            String contextLog = recordContext(new File(outputFolder), runCheckSum);    // get file MD5sum

            String runFileMsg = getTimeStamp() + ": Contents of the .run file: ";

            File logFile = new File(outputFolder + "/" + logFileName);
            try {
                if (!logFile.exists()) {
                    logFile.createNewFile();
                }
                BufferedWriter bw = new BufferedWriter(new FileWriter(logFile.getAbsolutePath()));
                bw.write("Contents of the .properties file:\n" + propertiesFileContents);
                bw.write(runFileMsg + "\n" + runFileContents);
                bw.write(contextLog);
                bw.close();

            } catch (IOException ioe) {
                ioe.printStackTrace();
            }

            // redirect System.out and System.err to the log file
            PrintStream ps = null;
            try{
                ps = new PrintStream(new BufferedOutputStream(new FileOutputStream(outputFolder + "/" + logFileName, true)));
            }catch(FileNotFoundException fnfe) {
                fnfe.printStackTrace();
            };
            System.setOut(ps);
            System.setErr(ps);


            System.out.println(getTimeStamp() + " ***********Now initializing ProductionSNPCallerPlugin*********");
            Date start = new Date();

            String[] pluginArgs = createPluginArgs();
            ProductionSNPCallerPlugin pscp = new ProductionSNPCallerPlugin();
            pscp.setParameters(pluginArgs);
            pscp.performFunction(null);

            Date stop = new Date();
            System.out.println(getTimeStamp() + " *********** ProductionSNPCallerPlugin has completed*********");

            long startTime = start.getTime();
            long stopTime = stop.getTime();
            long diff = stopTime - startTime;
            long elapsedSeconds = diff / 1000;


            File toFile = new File(archiveDirectory + "/" + aFile.getName());

            boolean movedFile = false;
            try{
                Files.copy(aFile, toFile);
                movedFile = true;
            }catch (IOException ioe){}

            if(movedFile){
                System.out.println("Moved file " + aFile.getAbsolutePath() + " to " + toFile.getAbsolutePath());
            }else{
                String msg = "******* COULD NOT MOVE FILE " + aFile.getAbsolutePath() + " TO " + toFile.getAbsolutePath();
                System.out.println(msg);
                sendAlertNotification(emailSubjectBase + "- Error", msg);
            }

            // send email notification that a .run file has been processed
            SMTPClient sc = new SMTPClient(emailHost, emailAddresses);
            String emailSubject = "GBS Production Pipeline " + this.anInputFolder;
            String emailMsg = "Ran:\n " + this.anInputFolder +
                               "\n\n  Tassel Pipeline Execution Time: " + elapsedSeconds + " seconds" +
                               "\n\n Attachment:\n " + logFile.getAbsolutePath();
            try{
                sc.sendMessageWithAttachment(emailSubject, emailMsg, logFileName);
            }catch (javax.mail.MessagingException me) {   /* ignore */ }
        }
    }



    /*
     *  Write out XML file for each instance of the TasselPipeline
     *
     * @param anXMLFile
     */
        private String[] createPluginArgs(){


            String[] args = {
                    "-i", anInputFolder,
                    "-k", keyFile,
                    "-e", enzyme,
                    "-o", outputFolder,
                    "-m", topmFile
            };
        return args;

    }

    /**
     * Load application-wide properties file and initialize variables.
     * @param aFileIn
     * @return
     */
    private String loadApplicationConfiguration(String aFileIn){
        boolean loaded = false;
        Properties props = new Properties();
        try{
            props.load(new FileInputStream(aFileIn));
            loaded = true;
        }catch(IOException ioe){
            System.out.println("Problem loading application configuration file:"  + aFileIn);
            System.out.println("************** Example .properties file: ");
            System.out.println(exampleAppConfigFile);
            ioe.printStackTrace();
        }

        // the properties file must load successfully or exit with email notification
        if(!loaded){
            sendAlertNotification(emailSubjectBase + "- Error", "Properties file could not be loaded: " + aFileIn);
            System.exit(1);
        }

        String configurationElement =  "runFileSuffix";
        runFileSuffix =    props.getProperty(configurationElement);

        configurationElement =  "emailHost";     // server used for sending email
        emailHost =    props.getProperty(configurationElement);

        configurationElement =  "emailAddress";     // to whom the email notifications should be sent
        String address =    props.getProperty(configurationElement);;
        Iterable<String> results = Splitter.on(emailAddressDelimiters).split(address);
        emailAddresses = Iterables.toArray(results, String.class);

        configurationElement =  "runDirectory";     // directory where the .run files are expected to be
        runDirectory =    props.getProperty(configurationElement);
        String response = testInputDirectory(aFileIn, runDirectory, configurationElement);
        if(response != null)  System.out.println(response);

        configurationElement =  "archiveDirectory";   // directory into which the key and run files are placed after being run
        archiveDirectory =    props.getProperty(configurationElement);
        response = testInputDirectory(aFileIn, archiveDirectory, configurationElement);
        if(response != null)  System.out.println(response);

        // read in the contents of the properties file so it can be placed into the log
        BufferedReader br = null;
        StringBuffer sb = new StringBuffer();
        try{
            br = new BufferedReader(new FileReader(aFileIn));
            String line = null;
            while( ( line = br.readLine() ) != null){
                sb.append(line + "\n");
            }
        }catch(IOException ioe){

        }
        return sb.toString();
    }

    /**
     * Load a file containing the information necessary to write out an
     * XML output file necessary for running the production pipeline
     *
     * @param aFileIn
     * @return
     */
    private String loadRunConfiguration(File aFileIn){

        String usage = aFileIn.getName() + " is missing a run configuration element:  ";
        Properties props = new Properties();
        try{
            props.load(new FileInputStream(aFileIn));
        }catch(IOException ioe){
            System.err.println("Issue loading run configuration file: " + aFileIn.getName());
            ioe.printStackTrace();
        }
        String configurationElement =  "inputFolder";
        anInputFolder =    props.getProperty(configurationElement);
        String response = testInputDirectory(aFileIn.getName(), anInputFolder, configurationElement);
        if(response != null)  System.out.println(response);

        configurationElement = "enzyme";
        enzyme =    props.getProperty(configurationElement);
        response = testInputDirectory(aFileIn.getName(), anInputFolder, configurationElement);
        if(response != null)  System.out.println(response);

        configurationElement = "topmFile";
        topmFile =    props.getProperty(configurationElement);
        response = testInputDirectory(aFileIn.getName(), anInputFolder, configurationElement);
        if(response != null)  System.out.println(response);

        configurationElement = "outputFolder";
        outputFolder =    props.getProperty(configurationElement);
        response = testInputDirectory(aFileIn.getName(), anInputFolder, configurationElement);
        if(response != null)  System.out.println(response);

        configurationElement = "keyFile";
        keyFile =     props.getProperty(configurationElement);
        response = testInputDirectory(aFileIn.getName(), anInputFolder, configurationElement);
        if(response != null)  System.out.println(response);
        
        BufferedReader br = null;
        StringBuffer sb = new StringBuffer();
        try{
            br = new BufferedReader(new FileReader(aFileIn));
            String line = null;
            while( ( line = br.readLine() ) != null){
                sb.append(line + "\n");
            }
        }catch(IOException ioe){
            
        }
        return sb.toString();
    }


    /**
     * Verify that element specifying a directory exists in a configuration or properties file
     * and that the directory exists on the filesystem.
     *
     * @param filename  Name of configuration/properties/.run file
     * @param input   String containing the fully qualified path of a directory
     * @param configurationElement
     * @return Description of the problem with the input.  If return is null then there is no problem with the entered file or directory
     */
    private String testInputDirectory(String filename, String input, String configurationElement){

        String response = null;
        if(input == null) {
            response = filename + " is missing a run configuration element:  " + configurationElement;
        }else{
            File aFileOrDir = new File(input);
            if(!aFileOrDir.exists()){
                response = filename + "'s configuration element " + configurationElement + " does not exist.  Please confirm path and filename.";
            }
        }
        return response;
    }


    /**
     * Collect information about the context in which the pipeline is being run.
     * This information includes the following:
     *      1) Current date and time.
     *      2) Name of the machine on which run occurs.
     *      3) User account used
     *      4) Files used and their MD5sums
     *      5) Contents of XML configuration file used to run TasselPipeline
     * @param outputFolder
     * @param calculateChecksum  Calculating MD5sum can be time consuming.  This allows it to be skipped when appropriate.
     * @return
     */

    private String recordContext(File outputFolder, boolean calculateChecksum){

        StringBuffer sb = new StringBuffer();
        
        // date
        sb.append(getTimeStamp()  + "\n");

        // user account name
        String userMsg = "User Account Name: ";
        String user = System.getProperty("user.name");
        sb.append(userMsg + user + "\n");

        // hostname
        String hostNameMsg = "Name of Machine on which JVM is Running: ";
        try {
            hostName = InetAddress.getLocalHost().getHostName();
        } catch (UnknownHostException e) {
            e.printStackTrace();
        }
        sb.append(hostNameMsg + hostName + "\n");


        // for each file in a directory, include the md5sum
        if(calculateChecksum){
            sb.append(getTimeStamp() + CheckSum.getMD5Checksum(keyFile) + "\n");
            sb.append(getTimeStamp() + CheckSum.getMD5Checksum(anInputFolder)+ "\n");
            sb.append(getTimeStamp() + CheckSum.getMD5Checksum(topmFile) + "\n") ;
        }else{
            sb.append(getTimeStamp() + "MD5sum checking has been switched off using the --skipCheckSum argument");
        }

        // todo: record plugin args
//        sb.append("XML file used to run the pipeline:\n" );
//        String[] xmlFileContents = createXMLFile(fileNameBase + ".xml");
//        tasselPipelineArgs = createXMLFile(fileNameBase + ".xml");
//        for(String arg: xmlFileContents) System.out.println("arg = " + arg);

//        sb.append(getTimeStamp() + "\n" + xmlFileContents + "\n");
        return sb.toString();
    }

    /**
     * Run the TasselPipeline by passing it the flags and values.
     * @param args  String[] containing flags and values
     */
    private void runTassel(String[] args){
        TasselPipeline tp = new TasselPipeline(args, null);
        tp.main(args);
    }

    private void runProductionSNPCallerPlugin(){

    }
    /**
     * Convenience method to provide uniformly labelled timestamps
     * @return
     */
    private String getTimeStamp(){
        String label = "Timestamp: ";
        Date now = new Date();
        return label + dateFormat.format(now) + " ";
    }

    /**
     * Send email to pipeline administrator when issue is encountered that
     * may keep pipeline from successfully completing, e.g., cannot load .properties file,
     * or send notification of how many times the pipeline should run, i.e., how many
     * .run files are present.
     * @param subject   Subject line of email
     * @param message   Message body of email
     */
    private void sendAlertNotification(String subject, String message){
        SMTPClient sc = new SMTPClient(emailHost, emailAddresses);
        try{
            sc.sendMessage(subject, message);
        }catch (javax.mail.MessagingException me) { /* ignore */  }
    }


    public static void main(String[] args){

        String msg = "--skipCheckSum flag allows MD5sum checking to be skipped. Use with discretion.";

        boolean doCheckSum = true;
        String skipCheckSum = "skipCheckSum";
        if(args != null){
            for(String arg: args){
                if(StringUtils.containsIgnoreCase(arg, skipCheckSum )){
                    doCheckSum = false;
                    break;
                }
            }
        }else{
            System.out.println(skipCheckSum);
        }

    	 new ProductionPipelineMain(doCheckSum);
    }
}
