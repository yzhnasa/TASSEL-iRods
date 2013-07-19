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
import java.util.Calendar;
import java.util.Properties;

import com.google.common.collect.Iterables;
import net.maizegenetics.pipeline.TasselPipeline;
import net.maizegenetics.util.CheckSum;
import net.maizegenetics.util.SMTPClient;
import com.google.common.base.Splitter;



/**
 * This class is for running the GBS Production Pipeline.  It is to be run from within the sTASSEL.jar.  The cron job
 * should be set up to run the run_pipeline.pl which has been modified to make this class the main() to be run.  The
 * JVM memory settings within run_pipeline.pl should also be adjusted upwards.
 *
 * 20130718 Note from Jeff Glaubitz: A minor (and not yet documented) detail:  SeqToGenos needs the key file name to end with "_key.txt".  All the output files are named after the key file but replacing "_key.txt" with the appropriate extension.
 * User: dkroon
 * Date: 4/8/13
 */
public class ProductionPipelineMain {


    // these could reasonably be contents of a properties file
    private String applicationConfiguration = "production_pipeline.properties";

    private String runFileSuffix = ".run";
    private String emailHost = "appsmtp.mail.cornell.edu";        // server to be used to send email notifications
    private String[] emailAddresses = {"dek29@cornell.edu"};            // to whom the email notifications should be sent.

    private String emailAddressDelimiters = "[.,;:]";           // any of these are acceptable delimiters between email addresses

    private String runDirectory = "/workdir/prop_pipeline/run/";  // directory into which the key and run files are placed for pickup by cron job
    private String archiveDirectory = "/workdir/prop_pipeline/done/";  // directory into which the key and run files are placed after being run
                                                                       // a dated directory will be created within this directory and artifacts
                                                                       // such as the .log and .xml file will be placed there
    private boolean runCheckSum = true;             // running the md5 checksum is time consuming so skip for testing

    private String todayDate = null;
    private String fileNameBase = null;  //  todayDate + property filename   to be used for XML and Log files//
    private String anInputFolder= null;
    private String enzyme= null;
    private String topmFile= null;
    private String outputFolder= null;
    private String keyFile= null;
    private String hostName = "host unknown";

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
        loadApplicationConfiguration(applicationConfiguration);
        // obtain directory in which the jar is being executed
        File baseLocation = new File(runDirectory);
            System.out.println("baseLocation = " + baseLocation);


        File dir = new File(baseLocation.getPath());

        if(!dir.exists()){
            System.out.println("Could not find the directory containing .run files: " + baseLocation.getPath());
            System.out.println("Exiting program.");
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
        }

        for(File aFile: files){

            String propFileContents = loadRunConfiguration(aFile);
            SimpleDateFormat yyyyMMdd_format = new SimpleDateFormat("yyyyMMdd");
            todayDate = yyyyMMdd_format.format(new Date());

            fileNameBase =   todayDate + "_"  + aFile.getName();

            String contextLog = recordContext(new File(outputFolder));    // get file checksums
            String propFileMsg = "Contents of the properties file: ";

            File logFile = new File(outputFolder + "/" + fileNameBase + ".log");
            try {
                if (!logFile.exists()) {
                    logFile.createNewFile();
                }
                BufferedWriter bw = new BufferedWriter(new FileWriter(logFile.getAbsolutePath()));
                bw.write(propFileMsg + "\n" + propFileContents);
                bw.write(contextLog);
                bw.close();

            } catch (IOException ioe) {
                ioe.printStackTrace();
            }

            PrintStream ps = null;
            try{
                ps = new PrintStream(new BufferedOutputStream(new FileOutputStream(fileNameBase + ".log", true)));
            }catch(FileNotFoundException fnfe) { /* ignore */ };
            System.setOut(ps);
            System.setErr(ps);

            String xmlFilename = fileNameBase + ".xml";

            for(String address: emailAddresses){
                SMTPClient sc = new SMTPClient(emailHost, address);
                String emailMsg = "Ran " + fileNameBase + ".xml";
                try{
                    sc.sendMessage("Ran job", emailMsg);
                }catch (javax.mail.MessagingException me) {
                    me.printStackTrace();
                }
            }
            String xmlFilePath = outputFolder  + "/" + xmlFilename;
            System.out.println("xmlFilePath = " + xmlFilePath);
           runTassel(xmlFilePath);
        }
    }

    /*
     *  write out XML file for each chromosome being run
     *
     * @param anXMLFile
     */
    private String createXMLFile( String anXMLFile){

        StringBuffer sb = new StringBuffer();
        sb.append("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n");
        sb.append(" <TasselPipeline>\n");
        sb.append("    <fork1>\n");
        sb.append("        <SeqToGenosPlugin>\n");
        sb.append("            <i>" + this.anInputFolder + "</i>\n");
        sb.append("            <k>" + this.keyFile + "</k>\n");
        sb.append("            <e>" + this.enzyme + "</e>\n");
        sb.append("            <o>" + this.outputFolder + "</o>\n");
        sb.append("            <m>" + this.topmFile + "</m>\n");
        sb.append("        </SeqToGenosPlugin>\n");
        sb.append("    </fork1>\n");
        sb.append("    <runfork1/>\n");
        sb.append("</TasselPipeline>\n");


        File xmlFile = new File(this.outputFolder + "/" + anXMLFile);
        try {
            if (!xmlFile.exists()) {
                xmlFile.createNewFile();
            }
            BufferedWriter bw = new BufferedWriter(new FileWriter(xmlFile.getAbsolutePath()));
            bw.write(sb.toString());
            bw.close();
            
        } catch (IOException ioe) {
            ioe.printStackTrace();
        }


        return sb.toString();
    }

    /**
     * Load application-wide properties file and initialize variables.
     * @param aFileIn
     * @return
     */
    private String loadApplicationConfiguration(String aFileIn){
        Properties props = new Properties();
        try{
            props.load(new FileInputStream(aFileIn));
        }catch(IOException ioe){
            System.err.println("Problem loading application configuration file:"  + aFileIn);
            ioe.printStackTrace();
        }

        String configurationElement =  "runFileSuffix";
        runFileSuffix =    props.getProperty(configurationElement);

        configurationElement =  "emailHost";     // to whom the email notifications should be sent
        emailHost =    props.getProperty(configurationElement);


        // directory into which the key and run files are placed for pickup by cron job
        configurationElement =  "emailAddress";     // to whom the email notifications should be sent
        String address =    props.getProperty(configurationElement);;
        Iterable<String> results = Splitter.on(emailAddressDelimiters).split(address);
        emailAddresses = Iterables.toArray(results, String.class);


        // directory into which the key and run files are placed after being run
        configurationElement =  "runDirectory";     // to whom the email notifications should be sent
        runDirectory =    props.getProperty(configurationElement);
        String response = testInputDirectory(aFileIn, runDirectory, configurationElement);
        if(response != null)  System.out.println(response);

        configurationElement =  "archiveDirectory";     // to whom the email notifications should be sent
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


    // record files used in run_pipeline
    // create XML file and records its contents

    private String recordContext(File outputFolder){

        StringBuffer sb = new StringBuffer();
        
        // date
        Calendar cal = Calendar.getInstance();
        SimpleDateFormat sdf = new SimpleDateFormat("yyyyMMdd HH:mm:ss");
        String dateMsg = "Date: ";
        sb.append(dateMsg + sdf.format(cal.getTime()) + "\n");

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
        if(runCheckSum){
            File prodFile = new File(topmFile);
            sb.append(CheckSum.getMD5Checksum(prodFile) + "\n") ;
            File inputFolder = new File(anInputFolder);
            if(inputFolder.isDirectory()){
                File[] fastqFile = inputFolder.listFiles();
                for(File aFile: fastqFile){
                    sb.append(CheckSum.getMD5Checksum(aFile) + "\n");
                }
            }
            sb.append(CheckSum.getMD5Checksum(inputFolder)+ "\n");
            File keyfileFile = new File(keyFile);
            sb.append(CheckSum.getMD5Checksum(keyfileFile) + "\n");
        }else{
            sb.append("Checksum calculation has been switched off.  To turn it back on, please edit the pipeline.properties file.\n");
        }

        sb.append("XML file used to run the pipeline:\n" );
        String xmlFileContents = createXMLFile(fileNameBase + ".xml");

        sb.append(xmlFileContents + "\n");
        return sb.toString();
    }

    private void runTassel(String xmlFile){

        System.out.println("The XML file for the TASSEL Pipeline to run: " + xmlFile);
        String[] args = { "-configFile", xmlFile};

        TasselPipeline tp = new TasselPipeline(args, null);
        tp.main(args);
    }
    
    public static void main(String[] args){

         boolean checkSumWanted = false;
    	 new ProductionPipelineMain(checkSumWanted);
    }
}
