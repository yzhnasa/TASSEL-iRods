/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.pal.site;

import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.Locus;
import java.util.HashMap;

/**
 * Container class for holding information about the site.  This includes information
 * on position, MAF, coverage,
 *
 * @author Ed Buckler
 */
public class AnnotatedSite {
    //Required 
    public final Locus myLocus;
    public final int myPosition;
    public final byte myStrand;
    public String mySNPIDs=null;
    //equals should be based on these
    
    //ideas for locus
    //try to parse to number,  not set to number Int.Max.  Sort by number, text.  
    
    
    public byte[] myAlleleFreqOrder=null;
    public int[] myAlleleCnt=null;
    public float myMAF= Float.NaN;;
    public float mySiteCoverage= Float.NaN;;
    
    public byte myReferenceAllele= Alignment.UNKNOWN_ALLELE;
    public byte myAncestralAllele=Alignment.UNKNOWN_ALLELE;
    public byte myMajorAllele=Alignment.UNKNOWN_ALLELE;;
    public byte myHighDepthAllele=Alignment.UNKNOWN_ALLELE;
    private HashMap myAnnoMap=null;


    public AnnotatedSite(Locus locus, int position, byte strand) {
        this.myLocus=locus;
        this.myPosition = position;
        this.myStrand=strand;
    }

    public AnnotatedSite(int position, byte[] myAlleleFreqOrder, int[] myAlleleCnt, float maf, float siteCov, String mySNPIDs) {
        this(null,position,(byte)1);
        this.myAlleleFreqOrder = myAlleleFreqOrder;
        this.myAlleleCnt = myAlleleCnt;
        this.myMAF = maf;
        this.mySiteCoverage = siteCov;
        this.mySNPIDs = mySNPIDs;
    }

    public void addAnnotation(String annoName, String value) {
        if(myAnnoMap==null) myAnnoMap=new HashMap<String,Object>(10);
        myAnnoMap.put(annoName, value);
    }

    public Object getAnnotation(String annoName) {
        if(myAnnoMap==null) return null;
        return myAnnoMap.get(annoName);
    }

    public double getAnnotationAsNumber(String annoName) {
        if(myAnnoMap==null) return Float.NaN;
        Object r=myAnnoMap.get(annoName);
        if(r instanceof Number) return ((Number)r).doubleValue();
        return Float.NaN;
    }

    public String getAnnotationAsString(String annoName) {
        if(myAnnoMap==null) return "";
        Object r=myAnnoMap.get(annoName);
        return r.toString();
    }
    
    public int[][] getAllelesSortedByFrequency() {
        int result[][] = new int[2][myAlleleCnt.length];
            for (int i = 0; i < myAlleleFreqOrder.length; i++) {
               result[0][i]=myAlleleFreqOrder[i];
               result[1][i]=myAlleleCnt[i];
           }
         return result;
    }
    
    public int getAlleleTotal() {
        int total=0;
        for (int b : myAlleleCnt) {total+=b;}
        return total;
    }
    
    
}
