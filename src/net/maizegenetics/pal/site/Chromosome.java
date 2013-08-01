/*
 * Locus
 */
package net.maizegenetics.pal.site;

import net.maizegenetics.pal.alignment.*;
import com.google.common.collect.ComparisonChain;
import java.io.Serializable;

import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

/**
 * Defines the chromosome structure and length.  The name and length recorded for 
 * each chromosome. 
 * 
 * @author terry
 */
public class Chromosome implements Serializable, Comparable<Chromosome>{

    private static final long serialVersionUID = -5197800047652332969L;
    public static Locus UNKNOWN = new Locus("Unknown");
    private final String myName;
    private final int myChromosomeNumber;
    private final int myLength;
    private final Map<String, Integer> myFeatures;

    /**
     *
     * @param name Name of the chromosome 
     * @param length Length of chromosome in base pairs
     * @param features Map of features about the chromosome
     */
    public Chromosome(String name, int length, Map<String, Integer> features) {
        myName = name;
        myLength = length;
        int convChr=Integer.MAX_VALUE;
        try{convChr=Integer.parseInt(name);}
        catch(NumberFormatException ne) {}
        myChromosomeNumber=convChr;
        if(features==null) {myFeatures=null;}
        else {myFeatures = new HashMap<String, Integer>(features);}
    }
    
    public Chromosome(String name) {
        this(name, -1, null);
    }

    public static Chromosome getMergedInstance(Chromosome chr1, Chromosome chr2) {
        String chromosome = chr1.getName();
        if (!chromosome.equals(chr2.getName())) {
            throw new IllegalArgumentException("Locus: getInstance: Chromosome Names must be the same.  Chromosome1: " + chromosome + "  Chromosome2: " + chr2.getName());
        }
        int length = Math.max(chr1.getLength(), chr2.getLength());
        Map features1 = chr1.getFeatures();
        Map features2 = chr2.getFeatures();
        Map newfeatures = new HashMap<String, Integer>(features1);
        newfeatures.putAll(features2);
        return new Chromosome(chromosome, length, newfeatures);
    }

    public String getName() {
        return myName;
    }
    
    /**
     * Returns the interger value of the chromosome (if name is not a number then
     * Integer.MAX_VALUE is returned)
     */
    public int getChromosomeNumber() {
        return myChromosomeNumber;
    }

    public int getLength() {
        return myLength;
    }

    public int getFeaturePosition(String key) {
        if(myFeatures==null) return Integer.MIN_VALUE;
        return myFeatures.get(key);
    }
    
    protected Map getFeatures() {
        return Collections.unmodifiableMap(myFeatures);
    }

    @Override
    public String toString() {
        return getName();
    }

    @Override
    public int hashCode() {
        int hash = 7;
        hash = 79 * hash + (this.myName != null ? this.myName.hashCode() : 0);
        return hash;
    }
    
    @Override
    public boolean equals(Object obj) {
        return (compareTo((Chromosome)obj)==0);
    }

    @Override
    public int compareTo(Chromosome o) {
        return ComparisonChain.start()
                .compare(myChromosomeNumber,o.getChromosomeNumber())
                .compare(myName,o.getName())
                .result();
    }
}
