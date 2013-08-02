package net.maizegenetics.pal.site;

import com.google.common.collect.ComparisonChain;

/**
 * Defines the central attributes of a position in the genome.
 */
public class CorePosition implements Position {
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


    @Override
    public Chromosome getLocus() {
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
