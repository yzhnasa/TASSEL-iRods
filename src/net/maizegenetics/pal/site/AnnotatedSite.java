/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.pal.site;

import com.google.common.collect.ComparisonChain;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.Locus;
import java.util.HashMap;
import java.util.Map;
import java.util.Objects;

/**
 * Container class for holding information about the site.  This includes information
 * on position, MAF, coverage.  This class is immutable.
 * 
 * Locus, position, cM, strand, and name (ID) are all set on instantiation.  They are sorted
 * based on the order listed above.  
 *
 * @author Ed Buckler
 */
public class AnnotatedSite implements Comparable<AnnotatedSite>{
    /**Locus of the site (required)*/
    private final Locus myLocus;
    /**Physical position of the site (unknown = Float.NaN)*/
    private final int myPosition;
    /**Strand of the site (unknown = Byte.MIN_VALUE)*/
    private final byte myStrand;
    /**Genetic position in centiMorgans (unknown = Float.NaN)*/
    private final float cM;
    /**Name of the site (default = SLocus_Position)*/
    private final String mySNPID;
    //equals should be based on these


    //These are used to support alignments and their analyses
    private final byte[] myAlleleFreqOrder;
    private final int[] myAlleleCnt;
    private final float myMAF;
    private final float mySiteCoverage;
    private final byte myMajorAllele;
    
    
    private byte myReferenceAllele= Alignment.UNKNOWN_ALLELE;
    private byte myAncestralAllele=Alignment.UNKNOWN_ALLELE;
    private byte myHighDepthAllele=Alignment.UNKNOWN_ALLELE;
    
    //Custom annotation are stored in the map
    private HashMap myAnnoMap=null;


    public AnnotatedSite(Locus locus, int position, float cM, byte strand, String snpID) {
        this(locus, position, cM, strand, snpID, null, null, Float.NaN, Float.NaN, 
            Alignment.UNKNOWN_ALLELE, Alignment.UNKNOWN_ALLELE, Alignment.UNKNOWN_ALLELE, Alignment.UNKNOWN_ALLELE,
            null);
    }

    public AnnotatedSite(Locus locus, int position, float cM, byte strand, String snpID,
            byte[] alleleFreqOrder, int[] alleleCnt, float maf, float siteCov, byte majorAllele, 
            byte referenceAllele, byte ancestralAllele, byte highDepthAllele, Map annoMap) {
        this.myLocus=locus;
        this.myPosition = position;
        this.cM=cM;
        this.myStrand=strand;
        this.mySNPID=snpID;
  
        this.myAlleleFreqOrder = alleleFreqOrder;
        this.myAlleleCnt = alleleCnt;
        this.myMAF = maf;
        this.mySiteCoverage = siteCov;
        
        this.myMajorAllele=majorAllele;
        this.myReferenceAllele=referenceAllele;
        this.myAncestralAllele=ancestralAllele;
        this.myHighDepthAllele=highDepthAllele;
        
        if(annoMap==null) {this.myAnnoMap=null;}
        else {myAnnoMap=new HashMap(annoMap);}
    }

    public void addAnnotation(String annoName, String value) {
        if(myAnnoMap==null) myAnnoMap=new HashMap<String,Object>(10);
        myAnnoMap.put(annoName, value);
    }

    public Object getAnnotation(String annoName) {
        if(myAnnoMap==null) return null;
//        switch (annoName) {  //TODO: uncomment once in Java 7
//            case "locus":return myLocus;
//            case "position":return myPosition;
//            case "cM":return cM;
//            case "strand":return myStrand;
//            case "snpID":return mySNPID;
//        }
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

    @Override
    public int hashCode() {
        int hash = 7;
        hash = 37 * hash + Objects.hashCode(this.myLocus);
        hash = 37 * hash + this.myPosition;
        hash = 37 * hash + this.myStrand;
        hash = 37 * hash + Float.floatToIntBits(this.cM);
        hash = 37 * hash + Objects.hashCode(this.mySNPID);
        return hash;
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        return (compareTo((AnnotatedSite) obj)==0);
    }
    
    
    @Override
    public int compareTo(AnnotatedSite o) {
        return ComparisonChain.start()
                .compare(myLocus,o.getLocus())
                .compare(myPosition,o.getPosition())
                .compare(cM,o.getcM())
                .compare(myStrand,o.getStrand())
                .compare(mySNPID, o.getSNPID())
                .result();
    }

    /**Return the locus (generally a chromosome) of a site*/
    public Locus getLocus() {
        return myLocus;
    }

    /**Return the physical position of a site*/
    public int getPosition() {
        return myPosition;
    }

    /**Return the strand for a site definition*/
    public byte getStrand() {
        return myStrand;
    }

    /**Return the strand for a site definition*/
    public float getcM() {
        return cM;
    }

    /**Return the ID (name) for a site*/
    public String getSNPID() {
        if (mySNPID == null) {
            return "S" + getLocus().getChromosomeName() + "_" + myPosition;
        } else {
            return mySNPID;
        }
    }

    /**Return the alleles ordered by allele frequency*/
    public byte[] getAlleleFreqOrder() {
        return myAlleleFreqOrder;
    }

    /**Return the alleles frequency counts ordered from most frequent to least frequent*/
    public int[] getAlleleCnt() {
        return myAlleleCnt;
    }

    /**Return the minor allele frequency*/
    public float getMAF() {
        return myMAF;
    }

    /**Returns the proportion of genotypes scored at a given site*/
    public float getSiteCoverage() {
        return mySiteCoverage;
    }

    /**Returns the reference allele*/
    public byte getReferenceAllele() {
        return myReferenceAllele;
    }

     /**Returns the ancestral allele*/
    public byte getAncestralAllele() {
        return myAncestralAllele;
    }

    /**Returns the major allele*/
    public byte getMajorAllele() {
        return myMajorAllele;
    }

    /**Returns the high depth allele*/
    public byte getHighDepthAllele() {
        return myHighDepthAllele;
    }
    
}
