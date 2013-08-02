/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.pal.site;

import com.google.common.collect.ComparisonChain;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.Alignment.ALLELE_SCOPE_TYPE;

import java.util.HashMap;
import java.util.Map;
//import java.util.Objects;

/**
 * Container class for holding information about the site.  This includes information
 * on position, MAF, coverage.  This class is immutable.
 * 
 * Locus, position, myCM, strand, and name (ID) are all set on instantiation.  They are sorted
 * based on the order listed above.  
 *
 * @author Ed Buckler
 */
public class AnnotatedPosition implements Position {
    /**Locus of the site (required)*/
    protected final Chromosome myChromosome;
    /**Physical position of the site (unknown = Float.NaN)*/
    protected final int myPosition;
    /**Strand of the site (unknown = Byte.MIN_VALUE)*/
    protected final byte myStrand;
    /**Genetic position in centiMorgans (unknown = Float.NaN)*/
    protected final float myCM;
    /**Name of the site (default = SLocus_Position)*/
    protected final String mySNPID;
    /**Is type Nucleotide or Text*/
    protected final boolean isNucleotide;
    /**Whether the variant define the nature of the indel*/
    protected final boolean isIndel;
    /**Define the nature of the polymorphism {"ACTAT","-"} or {"A","C","G"} or {"100","103","106"}
     */
    protected final String[] myKnownVariants;
    
    //These could perhaps be an array of 
    private final ALLELE_SCOPE_TYPE myScopeForMAF;
    private final float myMAF;
    private final float mySiteCoverage;
    private byte myReferenceAllele;
    private byte myGlobalMajorAllele;
    private byte myAncestralAllele;
    private byte myHighDepthAllele;
    
    //Custom annotation are stored in the map
    private Map myAnnoMap=null;

    public static class Builder {
        // Required parameters
        private final Chromosome myChromosome;
        private final int myPosition;
        // Optional parameters - initialized to default values
        private byte myStrand=1;
        private float myCM=Float.NaN;
        private String mySNPID=null;
        private boolean isNucleotide=true;
        private boolean isIndel=false;
        private String[] myKnownVariants=null;

        private ALLELE_SCOPE_TYPE myScopeForMAF=ALLELE_SCOPE_TYPE.Global_Frequency;
        private float myMAF = Float.NaN;
        private float mySiteCoverage = Float.NaN;
        private byte myGlobalMajorAllele=Alignment.UNKNOWN_ALLELE;
        private byte myReferenceAllele=Alignment.UNKNOWN_ALLELE;
        private byte myAncestralAllele=Alignment.UNKNOWN_ALLELE;
        private byte myHighDepthAllele=Alignment.UNKNOWN_ALLELE;
        private Map myAnnoMap=null;
        public Builder(Chromosome myChromosome, int myPosition) {
            this.myChromosome = myChromosome;
            this.myPosition = myPosition;
        }
        public Builder strand(byte val) {myStrand = val; return this;}
        public Builder cM(float val) {myCM = val; return this;}
        public Builder snpName(String val) {mySNPID = val; return this;}
        public Builder nucleotide(boolean val) {isNucleotide = val; return this; }
        public Builder indel(boolean val) {isIndel = val; return this;}
        public Builder knownVariants(String[] val) {myKnownVariants = val; return this;}
        public Builder maf(float val) {myMAF = val; return this;}
        public Builder siteCoverage(float val) {mySiteCoverage = val; return this;}
        //TODO still need to write the code for the rest

        public AnnotatedPosition build() {
            return new AnnotatedPosition(this);
        }
    }
    private AnnotatedPosition(Builder builder) {
        this.myChromosome = builder.myChromosome;
        myPosition = builder.myPosition;
        myStrand = builder.myStrand;
        myCM = builder.myCM;
        mySNPID = builder.mySNPID;
        isNucleotide = builder.isNucleotide;
        isIndel = builder.isIndel;
        myKnownVariants = builder.myKnownVariants;

        myScopeForMAF= builder.myScopeForMAF;
        myMAF = builder. myMAF;
        mySiteCoverage = builder. mySiteCoverage;

        myGlobalMajorAllele= builder.myGlobalMajorAllele;
        myReferenceAllele= builder.myReferenceAllele;
        myAncestralAllele= builder.myAncestralAllele;
        myHighDepthAllele= builder.myHighDepthAllele;
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
//            case "myCM":return myCM;
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

    @Override
    public int hashCode() {
        //TODO:  this hash code should be stored
        int hash = 7;
        hash = 37 * hash + this.myChromosome.hashCode();
        hash = 37 * hash + this.myPosition;
        hash = 37 * hash + this.myStrand;
        hash = 37 * hash + Float.floatToIntBits(this.myCM);
        hash = 37 * hash + this.getSNPID().hashCode();
        return hash;
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == this) {return true;}
        if (!(obj instanceof Position)) {return false;}
        Position o=(Position)obj;
        int result= ComparisonChain.start()
                .compare(myPosition,o.getPosition())  //position is most discriminating for speed
                .compare(myChromosome,o.getLocus())
                .compare(myCM, o.getCM())
                .compare(myStrand,o.getStrand())
                .compare(getSNPID(), o.getSNPID())
                .result();
        return (result==0);
    }

    @Override
    public int compareTo(Position o) {
        return ComparisonChain.start()
                .compare(myChromosome,o.getLocus())
                .compare(myPosition,o.getPosition())
                .compare(myCM, o.getCM())
                .compare(myStrand,o.getStrand())
                .compare(getSNPID(), o.getSNPID())
                .result();
    }

    /**Return the locus (generally a chromosome) of a site*/
    @Override
    public Chromosome getLocus() {
        return myChromosome;
    }

    /**Return the physical position of a site*/
    @Override
    public int getPosition() {
        return myPosition;
    }

    /**Return the strand for a site definition*/
    @Override
    public byte getStrand() {
        return myStrand;
    }

    /**Return the strand for a site definition*/
    @Override
    public float getCM() {
        return myCM;
    }

    /**Return the ID (name) for a site*/
    @Override
    public String getSNPID() {
        if (mySNPID == null) {
            return "S" + getLocus().getName() + "_" + myPosition;
        } else {
            return mySNPID;
        }
    }

    @Override
    public boolean isNucleotide() {
        return isNucleotide;
    }

    @Override
    public boolean isIndel() {
        return isIndel;
    }

    @Override
    public String[] getMyKnownVariants() {
        return myKnownVariants;
    }
    
}
