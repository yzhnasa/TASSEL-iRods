/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.pal.site;

import com.google.common.collect.ImmutableMultimap;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.NucleotideAlignmentConstants;

import java.util.Arrays;
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
    private final GeneralAnnotation myGA;
    private final float myMAF;
    private final float mySiteCoverage;
    private final byte[] myAlleleValue;

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
        private byte[] myAlleles=new byte[Allele.COUNT];

        //in an general annotation object
        private ImmutableMultimap.Builder<String, Object> myAnnoMapBld=null;
        private ImmutableMultimap<String, Object> myAnnoMap=null;

        /**Constructor requires a Position before annotation of the position*/
        public Builder(Position aCorePosition) {
            this.myCorePosition = aCorePosition;
            Arrays.fill(myAlleles,Alignment.UNKNOWN_ALLELE);
        }
        /**Constructor requires a Position before annotation of the position*/
        public Builder(Chromosome chr, int position) {
            this(new CorePosition.Builder(chr, position).build());
        }
        /**Set Minor Allele Frequency annotation (default=Float.NaN)*/
        public Builder maf(float val) {myMAF = val; return this;}
        /**Set site coverage annotation (default=Float.NaN)*/
        public Builder siteCoverage(float val) {mySiteCoverage = val; return this;}
        /**Set allele annotation by Allele type (default=Alignment.UNKNOWN_ALLELE)*/
        public Builder allele(Allele aT, byte val) {myAlleles[aT.index()] = val; return this;}
        /**Add non-standard annotation*/
        public Builder addAnno(String key, String value) {
            if(myAnnoMapBld==null) {
                myAnnoMapBld=new ImmutableMultimap.Builder();
            }
            myAnnoMapBld.put(key, value);
            return this;
        }
        /**Add non-standard annotation*/
        public Builder addAnno(String key, Number value) {
            if(myAnnoMapBld==null) {
                myAnnoMapBld=new ImmutableMultimap.Builder();
            }
            myAnnoMapBld.put(key, value);
            return this;
        }

        public CoreAnnotatedPosition build() {
            if(myAnnoMapBld!=null) myAnnoMap=myAnnoMapBld.build();
            return new CoreAnnotatedPosition(this);
        }
    }
    private CoreAnnotatedPosition(Builder builder) {
        this.myCorePosition = builder.myCorePosition;

        myMAF = builder.myMAF;
        mySiteCoverage = builder.mySiteCoverage;
        myAlleleValue=builder.myAlleles;
        myGA = new AbstractAnnotation(builder.myAnnoMap);
    }

    @Override
    public String toString() {
        StringBuilder sb=new StringBuilder("Position");
        sb.append("\tChr:").append(getChromosome().getName());
        sb.append("\tPos:").append(getPosition());
        sb.append("\tName:").append(getSNPID());
        sb.append("\tMAF:").append(getGlobalMAF());
        sb.append("\tRef:").append(NucleotideAlignmentConstants.getHaplotypeNucleotide(getAllele(Allele.REF)));
        return sb.toString();
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
    public byte getAllele(Allele alleleType) {
        return myAlleleValue[alleleType.index()];
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

    @Override
    public Object[] getAnnotation(String annoName) {
        return myGA.getAnnotation(annoName);
//        switch (annoName) {  //TODO: uncomment once in Java 7
//            case "locus":return myLocus;
//            case "position":return myPosition;
//            case "myCM":return myCM;
//            case "strand":return myStrand;
//            case "snpID":return mySNPID;
//        }
//        return myGA.getAnnotation(annoName);
    }

    @Override
    public String[] getTextAnnotation(String annoName) {
        return myGA.getTextAnnotation(annoName);
    }

    @Override
    public double[] getQuantAnnotation(String annoName) {
        return myGA.getQuantAnnotation(annoName);
    }


    @Override
    public String getConsensusAnnotation(String annoName) {
        return myGA.getConsensusAnnotation(annoName);
    }

    @Override
    public double getAverageAnnotation(String annoName) {
        return myGA.getAverageAnnotation(annoName);
    }
    
}
