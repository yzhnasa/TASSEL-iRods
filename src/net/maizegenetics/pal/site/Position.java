package net.maizegenetics.pal.site;

/**
 * Defines a genomic positions and its known variants.  Includes attributes of chromosome, position, strand,
 * centiMorgans, name (or SNP ID), whether this position is a nucleotide, or includes an indel.
 *
 * @author Ed Buckler
 */
public interface Position extends Comparable<Position> {

    /**Return the locus (generally a chromosome) of a site*/
    Chromosome getChromosome();

    /**Return the physical position of a site*/
    int getPosition();

    /**Return the strand for a site definition*/
    byte getStrand();

    /**Return the strand for a site definition*/
    float getCM();

    /**Return the ID (name) for a site*/
    String getSNPID();

    /**Whether the position is a nucleotide position or another marker type (SSR, AFLP, RAPD, CNV, which are recoded
     * with text states)*/
    boolean isNucleotide();

    /**Whether the position includes indels, which would be defined in the variants*/
    boolean isIndel();

    /**Returns the nature of the polymorphism {"ACTAT","-"} or {"A","C","G"} or {"100","103","106"}
     */
    String[] getKnownVariants();
}
