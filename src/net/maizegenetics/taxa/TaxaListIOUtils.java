package net.maizegenetics.taxa;

/**
 * Utilities for reading and writing IdGroup and PedigreeIdGroups.
 * 
 * @author Ed Buckler
 */
public class TaxaListIOUtils {

    private TaxaListIOUtils() {
    }
   /*
    public static TreeMultimap<String,String> getMapOfTextAnnotatedIds(TaxaList annoIdGroup, String annoName) {
        TreeMultimap<String,String> annoMap=TreeMultimap.create();
        for (int i = 0; i < annoIdGroup.numberOfTaxa(); i++) {
            if(annoIdGroup.getTaxon(i) instanceof AnnotatedIdentifier) {
                AnnotatedIdentifier ai=(AnnotatedIdentifier)annoIdGroup.getTaxon(i);
                for (String key : ai.getTextAnnotation(annoName)) {
                    annoMap.put(key,ai.getFullName());
                } 
            }
        }
        return annoMap;
    }
    
    public static TreeMultimap<String,String> getMapOfTextAnnotatedIds(TaxaList annoIdGroup, String keyAnnoName, String valueAnnoName) {
        TreeMultimap<String,String> annoMap=TreeMultimap.create();
        for (int i = 0; i < annoIdGroup.numberOfTaxa(); i++) {
            if(annoIdGroup.getTaxon(i) instanceof AnnotatedIdentifier) {
                AnnotatedIdentifier ai=(AnnotatedIdentifier)annoIdGroup.getTaxon(i);
                for (String key : ai.getTextAnnotation(keyAnnoName)) {
                    String[] values=ai.getTextAnnotation(valueAnnoName);
                    for (String value : values) {
                        annoMap.put(key,value);
                    }
                    
                } 
            }
        }
        return annoMap;
    }
    
    public static TaxaList readPedigree(String fileName) {
        BufferedReader fileIn = null;
        try {
            fileIn = Utils.getBufferedReader(fileName, 1000000);
            fileIn.mark(1<<16);
            String line;
            ArrayList<Taxon> taxaNames=new ArrayList<Taxon>();
            line=fileIn.readLine();
            int indexOfName=0;
            String[] headers=null;
            boolean[] isQuant=null;
            if(line.contains("<Name>")) {
                //parse headers
                headers=line.split("\\t");
                isQuant=new boolean[headers.length];
                for (int i = 0; i < headers.length; i++) {
                    if(headers[i].equals("<Name>")) {indexOfName=i; continue;}
                    headers[i]=headers[i].replace(">", "");
                    headers[i]=headers[i].replace("<", "");
                    if(headers[i].startsWith("#")) {
                        isQuant[i]=true;
                        headers[i]=headers[i].replace("#", "");
                    } else {
                        isQuant[i]=false;
                    }
                    
                }
            } else {
               fileIn.reset();
            }
            while((line=fileIn.readLine())!=null) {
                String[] s=line.split("\\t");
                if(headers!=null) {
                    AnnotatedIdentifier anID=new AnnotatedIdentifier(s[indexOfName]);
                    for (int i = 0; i < s.length; i++) {
                        if(i==indexOfName) continue;
                        String[] cs=s[i].split(";");
                        for(String ta: cs) {
                            if(isQuant[i]) {
                                if(ta.equals("NA")) {anID.addAnnotation(headers[i], Double.NaN);}
                                else {anID.addAnnotation(headers[i], Double.parseDouble(ta));}
                            }
                            else {anID.addAnnotation(headers[i], ta);}
                        }
                    }
                    taxaNames.add(anID);
                } else {
                    taxaNames.add(new Taxon(s[indexOfName]));
                }
            }
            return new SimpleIdGroup(taxaNames);
        } catch(Exception e) {
            System.err.println("Error in Reading Pedigree File:"+fileName);
            e.printStackTrace();
        }    
        return null;
        
    }
    
   */
    
}
