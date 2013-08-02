/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.pal.site;

import com.google.common.collect.ComparisonChain;
import net.maizegenetics.pal.alignment.Alignment;
import java.util.HashMap;
import java.util.Map;
import net.maizegenetics.pal.alignment.Alignment.ALLELE_SCOPE_TYPE;
//import java.util.Objects;

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
    private final Chromosome myChromosome;
    /**Physical position of the site (unknown = Float.NaN)*/
    private final int myPosition;
    /**Strand of the site (unknown = Byte.MIN_VALUE)*/
    private final byte myStrand;
    /**Genetic position in centiMorgans (unknown = Float.NaN)*/
    private final float cM;
    /**Name of the site (default = SLocus_Position)*/
    private final String mySNPID;
    /**Is type Nucleotide or Text*/
    private final boolean isNucleotide;
    /**Whether the variant define the nature of the indel*/
    private final boolean isIndel;
    /**Define the nature of the polymorphism {"ACTAT","-"} or {"A","C","G"} or {"100","103","106"}
     */
    private final String[] myKnownVariants;

    //Field below here should probably stay in annotatedSite, while those above go to CoreSite

    //These need to be in caches of genotype
//    private final byte[] myAlleleFreqOrder;
//    private final int[] myAlleleCnt;
//    private final byte myMajorAllele;
    
    //These could perhaps be an array of 
    private final ALLELE_SCOPE_TYPE myScopeForMAF;
    private final float myMAF;
    private final float mySiteCoverage;
    private byte myReferenceAllele= Alignment.UNKNOWN_ALLELE;
    private byte myGlobalMajorAllele= Alignment.UNKNOWN_ALLELE;
    private byte myAncestralAllele=Alignment.UNKNOWN_ALLELE;
    private byte myHighDepthAllele=Alignment.UNKNOWN_ALLELE;
    
    //Custom annotation are stored in the map
    private Map myAnnoMap=null;


    public AnnotatedSite(Chromosome locus, int position, float cM, byte strand, String snpID) {
        this(locus, position, cM, strand, snpID, true, false, null, ALLELE_SCOPE_TYPE.Global_Frequency, Float.NaN, Float.NaN, 
            Alignment.UNKNOWN_ALLELE, Alignment.UNKNOWN_ALLELE, Alignment.UNKNOWN_ALLELE, Alignment.UNKNOWN_ALLELE,
            null);
    }

    public AnnotatedSite(Chromosome locus, int position, float cM, byte strand, String snpID,
            boolean nucleotide, boolean indel, String[] knownVariants, ALLELE_SCOPE_TYPE scopeForMAF ,
            float maf, float siteCov, byte majorGlobalAllele, 
            byte referenceAllele, byte ancestralAllele, byte highDepthAllele, Map annoMap) {
        this.myChromosome=locus;
        this.myPosition = position;
        this.cM=cM;
        this.myStrand=strand;
        this.mySNPID=snpID;
        this.isNucleotide=nucleotide;
        this.isIndel=indel;
        this.myKnownVariants=knownVariants;
        
  
        this.myScopeForMAF=scopeForMAF;
        this.myMAF = maf;
        this.mySiteCoverage = siteCov;
        
        this.myGlobalMajorAllele=majorGlobalAllele;
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

    @Override
    public int hashCode() {
        //TODO:  this hash code should be stored
        int hash = 7;
        hash = 37 * hash + this.myChromosome.hashCode();
        hash = 37 * hash + this.myPosition;
        hash = 37 * hash + this.myStrand;
        hash = 37 * hash + Float.floatToIntBits(this.cM);
        hash = 37 * hash + this.getSNPID().hashCode();
        return hash;
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == this) {return true;}
        if (!(obj instanceof AnnotatedSite)) {return false;}
        AnnotatedSite o=(AnnotatedSite)obj;
        int result=ComparisonChain.start()
                .compare(myPosition,o.getPosition())  //position is most discriminating for speed
                .compare(myChromosome,o.getLocus())
                .compare(cM,o.getcM())
                .compare(myStrand,o.getStrand())
                .compare(mySNPID, o.getSNPID())
                .result();
        return (result==0);
    }
    
    
    @Override
    public int compareTo(AnnotatedSite o) {
        return ComparisonChain.start()
                .compare(myChromosome,o.getLocus())
                .compare(myPosition,o.getPosition())
                .compare(cM,o.getcM())
                .compare(myStrand,o.getStrand())
                .compare(mySNPID, o.getSNPID())
                .result();
    }

    /**Return the locus (generally a chromosome) of a site*/
    public Chromosome getLocus() {
        return myChromosome;
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
            return "S" + getLocus().getName() + "_" + myPosition;
        } else {
            return mySNPID;
        }
    }

    public boolean isIsNucleotide() {
        return isNucleotide;
    }

    public boolean isIsIndel() {
        return isIndel;
    }

    public String[] getMyKnownVariants() {
        return myKnownVariants;
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
    public byte getGlobalMajorAllele() {
        return myGlobalMajorAllele;
    }

    /**Returns the high depth allele*/
    public byte getHighDepthAllele() {
        return myHighDepthAllele;
    }
    
}
