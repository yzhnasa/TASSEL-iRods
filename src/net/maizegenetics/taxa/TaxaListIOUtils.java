package net.maizegenetics.taxa;

import com.google.common.collect.TreeMultimap;
import net.maizegenetics.util.Utils;

import java.io.BufferedReader;
import java.util.Set;

/**
 * Utilities for reading and writing IdGroup and PedigreeIdGroups.
 * 
 * @author Ed Buckler
 */
public class TaxaListIOUtils {

    private TaxaListIOUtils() {
    }

    public static TreeMultimap<String,Taxon> getMapOfTaxonByAnnotation(TaxaList taxaList, String annotation) {
        TreeMultimap<String,Taxon> annoMap=TreeMultimap.create();
        for (Taxon taxon : taxaList) {
            for (String value : taxon.getTextAnnotation(annotation)) {
                annoMap.put(value,taxon);
            }
        }
        return annoMap;
    }

    /**
     * Returns a subsetted taxa list based on annotation value.  For example, return all taxa where
     * {@literal GermType=Inbred}.
     * @param baseTaxaList base annotated taxa list
     * @param annotation annotation name (key)
     * @param annoValue  annotation value being tested for
     * @return TaxaList equal to the annotation value
     */
    public static TaxaList subsetTaxaListByAnnotation(TaxaList baseTaxaList, String annotation, String annoValue) {
        TaxaListBuilder tlb=new TaxaListBuilder();
        for (Taxon taxon : baseTaxaList) {
            for (String value : taxon.getTextAnnotation(annotation)) {
                if(value.equals(annoValue)) {
                    tlb.add(taxon);
                    break;
                }
            }
        }
        return tlb.build();
    }

    public static TaxaList retainSpecificAnnotations(TaxaList tl, String[] annotationsToKeep) {
      return null;
    }

    public static TaxaList removeSpecificAnnotations(TaxaList tl, String[] annotationsToRemove) {
        return null;
    }

    public static Set<String> allUniqueAnnotation(TaxaList tl) {
        return null;
    }


    /**
     * Returns an annotated TaxaList from a text annotation file in matrix format.  This is a tab delimited file.
     * First row in the file with the field {@literal taxaNameField} is the header row.  {@literal taxaNameField} indicated the taxon
     * name, all other fields are user defined.  The fields become the keys for the taxa annotation.
     * Quantitative fields should be tagged with "#" sign, e.g. {@literal <#INBREEDF>}.  Multiple values are supported per key, and
     * additional values can be either described with an additional column or ";" to delimit values with the same key.
     * <p></p>
     *{@literal <Name>	<GermType>	<Set>	<#InbreedF>	<Set>}<br></br>
     *{@literal B73	Inbred	Goodman282	0.98    ISU;IBMFounder}<br></br>
     *{@literal MO17	Inbred	Goodman282	0.98    UMC;IBMFounder}<br></br>
     * <p></p>
     * Produces:<br></br>
     * {@literal B73<GermType=Inbred;Set=Goodman282;Set=ISU;Set=IBMFounder;f=0.98>}<br></br>
     * {@literal MO17<GermType=Inbred;Set=Goodman282;Set=UMC;Set=IBMFounder;f=0.98>}<br></br>
     * The standardized keys are described in the {@link net.maizegenetics.taxa.Taxon}, and these constant fields are all upper
     * case.
     * @param fileName with complete path
     * @param taxaNameField field name with the taxon name
     * @return TaxaList with annotations
     */
    public static TaxaList readTaxaAnnotationFile(String fileName, String taxaNameField) {
        BufferedReader fileIn;
        try {
            fileIn = Utils.getBufferedReader(fileName, 1000000);
            fileIn.mark(1<<16);
            String line;
            TaxaListBuilder tlb=new TaxaListBuilder();
            line=fileIn.readLine();
            int indexOfName=0;
            String[] headers=null;
            boolean[] isQuant=null;
            if(line.contains(taxaNameField)) {
                //parse headers
                headers=line.split("\\t");
                isQuant=new boolean[headers.length];
                for (int i = 0; i < headers.length; i++) {
                    if(headers[i].equals(taxaNameField)) {indexOfName=i; continue;}
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
                    Taxon.Builder anID=new Taxon.Builder(s[indexOfName]);
                    for (int i = 0; i < s.length; i++) {
                        if(i==indexOfName) continue;
                        String[] cs=s[i].split(";");
                        for(String ta: cs) {
                            if(ta==null || ta.isEmpty()) continue;
                            if(isQuant[i]) {
                                if(ta.equals("NA")) {anID.addAnno(headers[i], Double.NaN);}
                                else {anID.addAnno(headers[i], Double.parseDouble(ta));}
                            }
                            else {anID.addAnno(headers[i], ta);}
                        }
                    }
                    tlb.add(anID.build());
                } else {
                    tlb.add(new Taxon(s[indexOfName]));
                }
            }
            return tlb.build();
        } catch(Exception e) {
            System.err.println("Error in Reading Pedigree File:"+fileName);
            e.printStackTrace();
        }    
        return null;
    }
    
}
