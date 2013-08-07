/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.pal.site;

import com.google.common.collect.HashMultimap;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.NucleotideAlignmentConstants;

import java.util.Map;
//import java.util.Objects;

/**
 * Provide information on a site and its annotations.  This includes information
 * on position, MAF, coverage.  This class is immutable.
 * <p></p>
 * The annotations are all set using the builder.
 *
 * @author Ed Buckler
 */
public final class CoreAnnotatedPosition implements AnnotatedPosition {
    private final Position myCorePosition;

    //These could perhaps be an array of
    private final float myMAF;
    private final float mySiteCoverage;
    private final byte myReferenceAllele;
    private final byte myGlobalMajorAllele;
    private final byte myAncestralAllele;
    private final byte myHighDepthAllele;

    //Custom annotation are stored in the map
    private final HashMultimap<String, Object> myAnnoMap=null;

    /**
     * A builder for creating immutable CoreAnnotatedPosition instances. AnnotatedPositions are
     * built off a base of a CorePosition, so build it first.
     *<p> Example:
     * <pre>   {@code
     * Position cp= new CorePosition.Builder(new Chromosome("1"),1232).build();
     * CoreAnnotatedPosition ap= new CoreAnnotatedPosition.Builder(cp)
     *    .maf(0.05f)
     *    .ancAllele(NucleotideAlignmentConstants.C_ALLELE)
     *    .build();}</pre>
     * <p>This would create nucleotide position on chromosome 1 at position 1232.  The MAF is 0.05 and the ancestral allele
     * is C.
     */
    public static class Builder {
        // Required parameters
        private final Position myCorePosition;

        //in an allele annotation objects
        private float myMAF = Float.NaN;
        private float mySiteCoverage = Float.NaN;
        private byte myGlobalMajorAllele=Alignment.UNKNOWN_ALLELE;
        private byte myReferenceAllele=Alignment.UNKNOWN_ALLELE;
        private byte myAncestralAllele=Alignment.UNKNOWN_ALLELE;
        private byte myHighDepthAllele=Alignment.UNKNOWN_ALLELE;

        //in an general annotation object
        private final Map myAnnoMap=null;

        /**Constructor requires a Position before annotation of the position*/
        public Builder(Position aCorePosition) {
            this.myCorePosition = aCorePosition;
        }
        /**Constructor requires a Position before annotation of the position*/
        public Builder(Chromosome chr, int position) {
            this.myCorePosition = new CorePosition.Builder(chr, position).build();
        }
        /**Set Minor Allele Frequency annotation (default=Float.NaN)*/
        public Builder maf(float val) {myMAF = val; return this;}
        /**Set site coverage annotation (default=Float.NaN)*/
        public Builder siteCoverage(float val) {mySiteCoverage = val; return this;}
        /**Set major allele annotation (default=Alignment.UNKNOWN_ALLELE)*/
        public Builder majAllele(byte val) {myGlobalMajorAllele = val; return this;}
        /**Set reference allele annotation (default=Alignment.UNKNOWN_ALLELE)*/
        public Builder refAllele(byte val) {myReferenceAllele = val; return this;}
        /**Set ancestral allele annotation (default=Alignment.UNKNOWN_ALLELE)*/
        public Builder ancAllele(byte val) {myAncestralAllele = val; return this;}
        /**Set high depth allele annotation (default=Alignment.UNKNOWN_ALLELE)*/
        public Builder hiDepthAllele(byte val) {myHighDepthAllele = val; return this;}


        public CoreAnnotatedPosition build() {
            return new CoreAnnotatedPosition(this);
        }
    }
    private CoreAnnotatedPosition(Builder builder) {
        this.myCorePosition = builder.myCorePosition;

        myMAF = builder.myMAF;
        mySiteCoverage = builder.mySiteCoverage;
        myGlobalMajorAllele= builder.myGlobalMajorAllele;
        myReferenceAllele= builder.myReferenceAllele;
        myAncestralAllele= builder.myAncestralAllele;
        myHighDepthAllele= builder.myHighDepthAllele;
    }

    @Override
    public String toString() {
        StringBuilder sb=new StringBuilder("Position");
        sb.append("\tChr:").append(getChromosome().getName());
        sb.append("\tPos:").append(getPosition());
        sb.append("\tName:").append(getSNPID());
        sb.append("\tMAF:").append(getGlobalMAF());
        sb.append("\tRef:").append(NucleotideAlignmentConstants.getHaplotypeNucleotide(myReferenceAllele));
        return sb.toString();
    }

    @Override
    public Object[] getAnnotation(String annoName) {
        if(myAnnoMap==null) return null;
//        switch (annoName) {  //TODO: uncomment once in Java 7
//            case "locus":return myLocus;
//            case "position":return myPosition;
//            case "myCM":return myCM;
//            case "strand":return myStrand;
//            case "snpID":return mySNPID;
//        }
        return myAnnoMap.get(annoName).toArray();
    }

    @Override
    public String[] getTextAnnotation(String annoName) {
        try{return myAnnoMap.get(annoName).toArray(new String[0]);}
        catch(Exception e) {
            return null;
        }
    }

    @Override
    public Double[] getQuantAnnotation(String annoName) {
        try{return myAnnoMap.get(annoName).toArray(new Double[0]);}
        catch(Exception e) {
            return null;
        }
    }

    @Override
    public float getGlobalMAF() {
        return myMAF;
    }

    @Override
    public float getGlobalSiteCoverage() {
        return mySiteCoverage;
    }

    @Override
    public byte getReferenceAllele() {
        return myReferenceAllele;
    }

    @Override
    public byte getAncestralAllele() {
        return myAncestralAllele;
    }

    @Override
    public byte getGlobalMajorAllele() {
        return myGlobalMajorAllele;
    }

    @Override
    public byte getHighDepthAllele() {
        return myHighDepthAllele;
    }

    @Override
    public int hashCode() {
        return myCorePosition.hashCode();
    }

    @Override
    public boolean equals(Object obj) {
        return myCorePosition.equals(obj);
    }

    @Override
    public int compareTo(Position o) {
        return myCorePosition.compareTo(o);
    }

    @Override
    public Chromosome getChromosome() {
        return myCorePosition.getChromosome();
    }

    @Override
    public int getPosition() {
        return myCorePosition.getPosition();
    }

    @Override
    public byte getStrand() {
        return myCorePosition.getStrand();
    }

    @Override
    public float getCM() {
        return myCorePosition.getCM();
    }

    @Override
    public String getSNPID() {
        return myCorePosition.getSNPID();
    }

    @Override
    public boolean isNucleotide() {
        return myCorePosition.isNucleotide();
    }

    @Override
    public boolean isIndel() {
        return myCorePosition.isIndel();
    }

    @Override
    public String[] getKnownVariants() {
        return myCorePosition.getKnownVariants();
    }
    
}
