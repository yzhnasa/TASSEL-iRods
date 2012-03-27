/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.pipeline;

import java.io.*;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;
import net.maizegenetics.pal.alignment.*;
import net.maizegenetics.pal.ids.*;
import net.maizegenetics.pal.popgen.BitNeighborFinder;
import net.maizegenetics.util.ProgressListener;
import net.maizegenetics.util.Utils;

/**
 *
 * @author edbuckler
 */
public class EdTests {
 //   String gbsFile="/Users/edbuckler/SolexaAnal/GBS/build111217/imputed/allZea20111217_scv10mF8maf002_mgs_E1pLD5.imp95_1024.c10.hmp.txt";
    String gbsFile="/Users/edbuckler/SolexaAnal/GBS/build20120110/imp/Zea20120110_scv10mF8maf002_mgs_E1pLD5kpUn.imp95_1024.c10.hmp.txt";
    //String hapFileAGP1="/Users/edbuckler/SolexaAnal/GBS/build111217/imputed/chr10.CSHLALLBGI.h90_f12.Q87Q87Union.hmp.txt";
    String hapFileAGP1="/Users/edbuckler/SolexaAnal/HapMapV2/HapMapV2RefgenV2_NAMfounders_teosinte_20110927_chr10.hmp.txt";
    String hapFileAGP2="/Users/edbuckler/SolexaAnal/GBS/build111217/imputed/chr10.CSHLALLBGI.h90_f12.Q87Q87Union.hmp.txt";
    TBitAlignment gbsMap=null;
    TBitAlignment hapMap=null;

    

    public EdTests() {
//        convertFilesToFast();
        gbsMap=(TBitAlignment)readGZOfSBit(gbsFile, false);
        hapMap=(TBitAlignment)readGZOfSBit(hapFileAGP1, false);
//        for(int i=0; i<hapMap.getSequenceCount(); i++) {
//            System.out.print("\""+hapMap.getTaxaName(i)+"\",");
//        }
//        System.out.println("");
       String[] taxa={"B73","B97","CML103","CML228","CML247","CML277","CML322","CML333","CML52","CML69","HP301",
           "IL14H","Ki11","Ki3","KY21","M162W","M37W","MO17","MO18W","MS71","NC350","NC358","OH43","OH7B","P39",
           "TIL01","TIL02","TIL03","TIL04","TIL05","TIL06","TIL07","TIL08","TIL09","TIL10","TIL11","TIL12",
           "TIL15","TIL16","TIL17","TIL25","TIL14","TX303","TZI8","W22","W64A"};
       IdGroup testHD=new SimpleIdGroup(taxa);
       //I am using this testHD as the taxa names from hapMap do not match gbsMap.  
       //We will probably need to load a TreeMap to redirect between the two.
       BitNeighborFinder bnf=new BitNeighborFinder(testHD, gbsMap, hapMap);
       
        
        
//        compareSitesInFiles();
 //       hapMap2=new MutableTBitNucleotideAlignment(hapMap, hapMap.getMaxNumAlleles(), hapMap.retainsRareAlleles());
//        convertHapMapSites();
//        hapMap=hapMap2;
//        compareSitesInFiles();
    }
  
    
    private void convertFilesToFast() {
        SBitAlignment gbs=(SBitAlignment)ImportUtils.readFromHapmap(gbsFile, (ProgressListener)null);
        TBitAlignment gbsOut=TBitAlignment.getInstance(gbs);
        writeAlignmentToSerialGZ(gbsOut, gbsFile);
        SBitAlignment hapmapIn=(SBitAlignment)ImportUtils.readFromHapmap(hapFileAGP1, (ProgressListener)null);
        TBitAlignment hapmapOut=TBitAlignment.getInstance(hapmapIn);
        writeAlignmentToSerialGZ(hapmapOut, hapFileAGP1);
    }
    
    private void compareSitesInFiles() {
        int count=0;
        Locus lx=hapMap.getLocus(0);
        for (int i = 0; i < gbsMap.getSiteCount(); i++) {
            int position=gbsMap.getPositionInLocus(i);
            int hapSite=hapMap.getSiteOfPhysicalPosition(position, lx);
            if(hapSite>0) System.out.printf("%d %d %d %n",count++, position, hapSite);
        }
    }
    
     
    public static void main(String[] args) {
        EdTests et=new EdTests();
        
    }
    
    
    
    static Alignment readGZOfSBit(String inFile, boolean isSBit) {
        SBitAlignment sba=null;
        TBitAlignment tba=null;
        long time=System.currentTimeMillis();
        try {
            File theFile = new File(Utils.addSuffixIfNeeded(inFile, ".gz"));
            System.out.println("Reading:"+theFile);
            FileInputStream fis = new FileInputStream(theFile);
            GZIPInputStream gs = new GZIPInputStream(fis);
            ObjectInputStream ois = new ObjectInputStream(gs);
            if(isSBit) {sba=(SBitAlignment)ois.readObject();
            System.out.println(sba.getSiteCount());
            System.out.println(sba.getSequenceCount());}
            else {tba=(TBitAlignment)ois.readObject();
            System.out.println(tba.getSiteCount());
            System.out.println(tba.getSequenceCount());}
            ois.close();
            fis.close();
            
        } catch (Exception ee) {
            //sendErrorMessage("Data could not be saved: " + ee);
            ee.printStackTrace();
        }
        System.out.println("Time:"+(System.currentTimeMillis()-time));
        if(isSBit) return sba;
        return tba;
    }
    
    static void writeAlignmentToSerialGZ(Alignment sba, String outFile) {
        //String inFile="/Users/edbuckler/SolexaAnal/GBS/build111217/imputed/allZea20111217_scv10mF8maf002_mgs_E1pLD5.imp95_1024.c10.hmp.txt";
        long time=System.currentTimeMillis();
        try {

            File theFile = new File(Utils.addSuffixIfNeeded(outFile, ".gz"));
            FileOutputStream fos = new FileOutputStream(theFile);
            GZIPOutputStream gz = new GZIPOutputStream(fos);
            ObjectOutputStream oos = new ObjectOutputStream(gz);
            oos.writeObject(sba);
            oos.flush();
            oos.close();
            fos.close();
            System.out.println("Wrote:"+theFile.toString());
        } catch (Exception ee) {
            //sendErrorMessage("Data could not be saved: " + ee);
            ee.printStackTrace();
        }
        System.out.println("Time:"+(System.currentTimeMillis()-time));
    }
}
