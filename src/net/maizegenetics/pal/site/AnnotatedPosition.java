package net.maizegenetics.pal.site;

/**
 * Annotations about the allele frequency and state of a polymorphism for various allele scopes.
 *
 */
public interface AnnotatedPosition extends GeneralAnnotation, Position {
    /**Return the minor allele frequency in a global scope*/
    public float getGlobalMAF();

    /**Returns the proportion of genotypes scored at a given site*/
    public float getGlobalSiteCoverage();

    /**Returns the reference allele*/
    public byte getReferenceAllele();

    /**Returns the ancestral allele*/
    public byte getAncestralAllele();

    /**Returns the major allele*/
    public byte getGlobalMajorAllele();

    /**Returns the high depth allele*/
    public byte getHighDepthAllele();
}
