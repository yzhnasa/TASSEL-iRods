package net.maizegenetics.pal.site;

import com.google.common.collect.ComparisonChain;

/**
 * Defines the central attributes of a position in the genome.
 */
public final class CorePosition implements Position {
    /**Locus of the site (required)*/
    private final Chromosome myChromosome;
    /**Physical position of the site (unknown = Float.NaN)*/
    private final int myPosition;
    /**Strand of the site (unknown = Byte.MIN_VALUE)*/
    private final byte myStrand;
    /**Genetic position in centiMorgans (unknown = Float.NaN)*/
    private final float myCM;
    /**Name of the site (default = SLocus_Position)*/
    private final String mySNPID;
    /**Is type Nucleotide or Text*/
    private final boolean isNucleotide;
    /**Whether the variant define the nature of the indel*/
    private final boolean isIndel;
    /**Define the nature of the polymorphism {"ACTAT","-"} or {"A","C","G"} or {"100","103","106"}
     */
    private final String[] myKnownVariants;
    private final int hashCode;

    /**
     * A builder for creating immutable CorePosition instances.
     *<p> Example:
     * <pre>   {@code
     * Position cp= new CorePosition.Builder(new Chromosome("1"),1232)
     *   .nucleotide(false)
     *   .knownVariants(new String[]{"A","C"})
     *   .cM(1.4f)
     *   .snpName("PZA123")
     *   .build();}</pre>
     * <p>This would create nucleotide position on chromosome 1 at position 1232 (1.4cM).  The position is named PZA123.
     */
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


        /**
         * Constructor for Builder, requires chromosome and position
         * @param myChromosome
         * @param myPosition
         */
        public Builder(Chromosome myChromosome, int myPosition) {
            this.myChromosome = myChromosome;
            this.myPosition = myPosition;
        }
        /**Set strand (default=1)*/
        public Builder strand(byte val) {myStrand = val; return this;}
        /**Set strand (default=Float.NaN)*/
        public Builder cM(float val) {myCM = val; return this;}
        /**Set SNP name (default="S"+Chromosome+"_"+position)*/
        public Builder snpName(String val) {mySNPID = val; return this;}
        /**Set whether position is nucleotide (default=true)*/
        public Builder nucleotide(boolean val) {isNucleotide = val; return this; }
        /**Set whether position is indel (default=false)*/
        public Builder indel(boolean val) {isIndel = val; return this;}
        /**Set text definition of variants (default=null)*/
        public Builder knownVariants(String[] val) {myKnownVariants = val; return this;}
        public CorePosition build() {
            return new CorePosition(this);
        }
    }
    private CorePosition(Builder builder) {
        myChromosome = builder.myChromosome;
        myPosition = builder.myPosition;
        myStrand = builder.myStrand;
        myCM = builder.myCM;
        mySNPID = builder.mySNPID;
        isNucleotide = builder.isNucleotide;
        isIndel = builder.isIndel;
        myKnownVariants = builder.myKnownVariants;
        hashCode=hashCode();
    }

    @Override
    public String toString() {
        StringBuilder sb=new StringBuilder("Position");
        sb.append("\tChr:").append(getChromosome().getName());
        sb.append("\tPos:").append(getPosition());
        sb.append("\tName:").append(getSNPID());
        return sb.toString();
    }

    @Override
    public int hashCode() {
        return this.hashCode;    //To change body of overridden methods use File | Settings | File Templates.
    }

    public int calcHashCode() {
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
                .compare(myPosition, o.getPosition())  //position is most discriminating for speed
                .compare(myChromosome, o.getChromosome())
                .compare(myCM, o.getCM())
                .compare(myStrand,o.getStrand())
                .result();
        if(result!=0) return false;
        return getSNPID().equals(o.getSNPID()); //This is done last as the string comparison is expensive
    }

    @Override
    public int compareTo(Position o) {
        int result=ComparisonChain.start()
                .compare(myChromosome,o.getChromosome())
                .compare(myPosition,o.getPosition())
//                .compare(myCM, o.getCM())
//                .compare(myStrand,o.getStrand())
                .result();
        if(result!=0) return result;
        return getSNPID().compareTo(o.getSNPID()); //This is done last as the string comparison is expensive
    }

    @Override
    public Chromosome getChromosome() {
        return myChromosome;
    }


    @Override
    public int getPosition() {
        return myPosition;
    }


    @Override
    public byte getStrand() {
        return myStrand;
    }


    @Override
    public float getCM() {
        return myCM;
    }


    @Override
    public String getSNPID() {
        if (mySNPID == null) {
            return (new StringBuilder("S").append(getChromosome().getName()).append("_").append(myPosition)).toString();
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
    public String[] getKnownVariants() {
        return myKnownVariants;
    }
}
