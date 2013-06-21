package net.maizegenetics.pal.ids;

import java.io.BufferedReader;
import java.util.ArrayList;
import java.util.HashMap;
import net.maizegenetics.util.Utils;

/**
 * Utilities for reading and writing IdGroup and PedigreeIdGroups.
 * 
 * @author Ed Buckler
 */
public class IdGroupIOUtils {

    private IdGroupIOUtils() {
    }
    
    public static IdGroup readPedigree(String fileName) {
        BufferedReader fileIn = null;
        try {
            fileIn = Utils.getBufferedReader(fileName, 1000000);
            fileIn.mark(1<<16);
            String line;
            ArrayList<Identifier> taxaNames=new ArrayList<Identifier>();
            line=fileIn.readLine();
            if(line.contains("<Name>")) {
                //parse headers
                String[] s=line.split("\\t");
                int countCat=0, countNum=0;
                for (String si : s) {
                    if(si.startsWith("<#")) {countNum++;}
                    else if(si.startsWith("<")) {countCat++;}
                }
                HashMap<String, String>[] catHash=new HashMap[countCat];
                HashMap<String, String>[] numHash=new HashMap[countNum];
                
                
            } else {
               fileIn.reset();
            }
            while((line=fileIn.readLine())!=null) {
                String[] s=line.split("\\t");
                taxaNames.add(new Identifier(s[0]));
            }
            return new SimpleIdGroup(taxaNames);
        } catch(Exception e) {
            System.err.println("Error in Reading Pedigree File:"+fileName);
            e.printStackTrace();
        }    
        return null;
        
    }
    
    
    
}
