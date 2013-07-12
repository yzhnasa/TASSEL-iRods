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
import net.maizegenetics.pipeline.TasselPipeline;
import net.maizegenetics.util.CheckSum;
import net.maizegenetics.util.SMTPClient;

/**
 * This class is for running the GBS Production Pipeline.  It is to be run from within the sTASSEL.jar.  The cron job
 * should be set up to run the run_pipeline.pl which has been modified to make this class the main() to be run.  The
 * JVM memory settings within run_pipeline.pl should also be adjusted upwards.
 * User: dkroon
 * Date: 4/8/13
 */
public class ProductionPipelineMain {


    private String projectName = "";  //  PI or Project Name
    private String pipelineName = "_ProdPipeline_";
    private String todayDate = null;
    private String fileNameBase = null;  //  projectName + pipelineName + todayDate  to be used for XML and Log files//    static String baseDir = null;
    private String anInputFolder= null;
    private String enzyme= null;
    private String topmFile= null;
    private static String outputFolder= null;
    private String keyFile= null;
    private String hostName = "host unknown";
    private static String emailHost = "appsmtp.mail.cornell.edu";
    private static String emailAddress = "dek29@cornell.edu";


    /*
     *  write out XML file for each chromosome being run
     *
     * @param outputFolder  Directory into which the XML file gets written
     * @param anXMLFile
     */
    private String createXMLFile(File outputFolder, File anXMLFile){

        StringBuffer sb = new StringBuffer();
        sb.append("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n");
        sb.append(" <TasselPipeline>\n");
        sb.append("    <fork1>\n");
        sb.append("        <RawReadsToHapMapPlugin>\n");
        sb.append("            <i>" + this.anInputFolder + "</i>\n");
        sb.append("            <k>" + this.keyFile + "</k>\n");
        sb.append("            <e>" + this.enzyme + "</e>\n");
        sb.append("            <o>" + this.outputFolder + "</o>\n");
        sb.append("            <m>" + this.topmFile + "</m>\n");
        sb.append("        </RawReadsToHapMapPlugin>\n");
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

    private String init(){
        
        String propFileContents = loadProperties("pipeline.properties");
        SimpleDateFormat yyyyMMdd_format = new SimpleDateFormat("yyyyMMdd");
        todayDate = yyyyMMdd_format.format(new Date());

        fileNameBase =  projectName + pipelineName + todayDate;

        String contextLog = recordContext(new File(outputFolder));
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
        return outputFolder + "/" + fileNameBase;
    }


    private String loadProperties(String filename){
        Properties props = new Properties();
        try{
            props.load(new FileInputStream(filename));
        }catch(IOException ioe){
            ioe.printStackTrace();
        }
        this.anInputFolder =    props.getProperty("inputFolder");
        this.enzyme =           props.getProperty("enzyme");
        this.topmFile =         props.getProperty("topmFile");
        this.outputFolder =     props.getProperty("outputFolder");
        this.keyFile =          props.getProperty("keyFile");
        
        BufferedReader br = null;
        StringBuffer sb = new StringBuffer();
        try{
            br = new BufferedReader(new FileReader(filename));
            String line = null;
            while( ( line = br.readLine() ) != null){
                sb.append(line + "\n");
            }
        }catch(IOException ioe){
            
        }
        return sb.toString();
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
//*      COMMENTED OUT FOR NOW BECAUSE IT ADDS TOO MUCH TIME TO TESTING
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
        

        sb.append("XML file used to run the pipeline:\n" );
        String xmlFileContents = createXMLFile(new File(this.outputFolder), new File(fileNameBase + ".xml"));

        sb.append(xmlFileContents + "\n");
        return sb.toString();
    }

    private void runTassel(String xmlFile){

        String[] args = { "-configFile", xmlFile};
        TasselPipeline tp = new TasselPipeline(args, null);
        tp.main(args);
    }
    
    public static void main(String[] args){
           
    	ProductionPipelineMain ppm = new ProductionPipelineMain();
    	String base = ppm.init();
        System.out.println("Completed init " + base);
       
        
        PrintStream ps = null;
        try{
            ps = new PrintStream(new BufferedOutputStream(new FileOutputStream(base + ".log", true)));
        }catch(FileNotFoundException fnfe) { /* ignore */ };
        System.setOut(ps);
        System.setErr(ps);
        
        String xmlFilename = base + ".xml";
        SMTPClient sc = new SMTPClient(emailHost, emailAddress);
        String emailMsg = "Ran " + base + ".xml";
        try{
            sc.sendMessage("Ran job", emailMsg);
        }catch (javax.mail.MessagingException me) {
            me.printStackTrace();
        }
    	ppm.runTassel(xmlFilename);
        
    }
}
