package net.maizegenetics.pal.site;

/**
 * Created with IntelliJ IDEA.
 * User: edbuckler
 * Date: 8/2/13
 * Time: 5:46 PM
 * To change this template use File | Settings | File Templates.
 */
public interface Position extends Comparable<Position> {

    /**Return the locus (generally a chromosome) of a site*/
    Chromosome getLocus();

    /**Return the physical position of a site*/
    int getPosition();

    /**Return the strand for a site definition*/
    byte getStrand();

    /**Return the strand for a site definition*/
    float getCM();

    /**Return the ID (name) for a site*/
    String getSNPID();

    boolean isNucleotide();

    boolean isIndel();

    /**Returns the nature of the polymorphism {"ACTAT","-"} or {"A","C","G"} or {"100","103","106"}
     */
    String[] getMyKnownVariants();
}
