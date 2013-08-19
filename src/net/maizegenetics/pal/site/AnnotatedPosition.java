package net.maizegenetics.pal.site;

/**
 * Annotations about the allele frequency and state of a polymorphism for various allele scopes.
 *
 */
public interface AnnotatedPosition extends GeneralAnnotation, Position {
    /**
     * Allele types recorded in an annotated position.  If unknown,
     * Alignment.UNKNOWN_ALLELE is returned.
     */
    public enum Allele {  //The indices are used in effectively as map (EnumMap is not used as it requires 4X more memory)
        /**Reference Allele*/
        REF(0),
        /**Major (most frequent) allele from the globally defined alignment*/
        GLBMAJ(1),
        /**Minor (second most frequent) allele from the globally defined alignment*/
        GLBMIN(2),
        /**Ancestral allele defined by evolutionary comparison*/
        ANC(3),
        /**High depth allele as defined from DNA sequencing analysis*/
        HIDEP(4);
        private final int index;
        /**Count of the number of allele types*/
        public final static int COUNT=Allele.values().length;
        Allele(int index) {this.index=index;}
        /**Sequential index that can be use for primitive arrays*/
        public int index() {return index;}
    }

    /**Return the minor allele frequency in a global scope*/
    public float getGlobalMAF();

    /**Returns the proportion of genotypes scored at a given site*/
    public float getGlobalSiteCoverage();

    /**Return the allele specified by alleleType, if unknown Alignment.Unknown is return*/
    public byte getAllele(Allele alleleType);

}
