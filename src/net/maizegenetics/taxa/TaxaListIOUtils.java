package net.maizegenetics.taxa;

import com.google.common.collect.ImmutableMultimap;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Multimap;
import com.google.common.collect.Ordering;
import net.maizegenetics.util.Utils;

import java.io.BufferedReader;
import java.util.Arrays;
import java.util.Map;
import java.util.Set;

/**
 * Utilities for reading and writing IdGroup and PedigreeIdGroups.
 * 
 * @author Ed Buckler
 */
public class TaxaListIOUtils {

    private TaxaListIOUtils() {
    }


    /**
     * Create a Multimap of all the taxa associated with a particular annotation value.
     * @param taxaList input taxa list with annotation associated with
     * @param annotation annotation key used to create the multimap, the values of these keys become the key of the
     *                   resulting Multimap
     * @return Map of AnnotationValues -> Taxon
     */
    public static Multimap<String,Taxon> getMapOfTaxonByAnnotation(TaxaList taxaList, String annotation) {
        ImmutableMultimap.Builder<String,Taxon> annoMap=new ImmutableMultimap.Builder<String,Taxon>().orderKeysBy(Ordering.natural());
        for (Taxon taxon : taxaList) {
            for (String value : taxon.getTextAnnotation(annotation)) {
                annoMap.put(value,taxon);
            }
        }
        return annoMap.build();
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

    /**
     * Creates a new taxa list with the taxa only retaining annotations within a specified list.  All taxa are retained,
     * only the annotations are changed.
     * @param baseTaxaList
     * @param annotationsToKeep the retained keys annotation
     * @return new TaxaList with a subset of the annotations
     */
    public static TaxaList retainSpecificAnnotations(TaxaList baseTaxaList, String[] annotationsToKeep) {
        Set<String> keepers=new ImmutableSet.Builder<String>().addAll(Arrays.asList(annotationsToKeep)).build();
        TaxaListBuilder tlb=new TaxaListBuilder();
        for (Taxon taxon : baseTaxaList) {
            Taxon.Builder tb=new Taxon.Builder(taxon.getName());
            for (Map.Entry<String,String> entry : taxon.getAllAnnotationEntries()) {
                if(keepers.contains(entry.getKey())) {tb.addAnno(entry.getKey(), entry.getValue());}
            }
            tlb.add(tb.build());
        }
        return tlb.build();
    }

    /**
     * Creates a new taxa list with the taxa retaining annotations EXCEPT those specified by the list.
     * All taxa are retained,
     * only the annotations are changed.
     * @param baseTaxaList
     * @param annotationsToRemove the retained keys annotation
     * @return new TaxaList with a subset of the annotations
     */
    public static TaxaList removeSpecificAnnotations(TaxaList baseTaxaList, String[] annotationsToRemove) {
        Set<String> keepers=new ImmutableSet.Builder<String>().addAll(Arrays.asList(annotationsToRemove)).build();
        TaxaListBuilder tlb=new TaxaListBuilder();
        for (Taxon taxon : baseTaxaList) {
            Taxon.Builder tb=new Taxon.Builder(taxon.getName());
            for (Map.Entry<String,String> entry : taxon.getAllAnnotationEntries()) {
                if(!keepers.contains(entry.getKey())) {tb.addAnno(entry.getKey(), entry.getValue());}
            }
            tlb.add(tb.build());
        }
        return tlb.build();
    }

    public static Set<String> allUniqueAnnotation(TaxaList tl) {

        return null;
    }

    public static void exportAnnotatedTaxaListTable() {
        //TODO
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
