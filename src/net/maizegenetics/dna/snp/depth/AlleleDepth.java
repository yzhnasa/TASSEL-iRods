/*
 *  AlleleDepth
 */
package net.maizegenetics.dna.snp.depth;

/**
 * Provides DNA read depth data for a genotype table.
 * @author Terry Casstevens
 */
public interface AlleleDepth {

    /**
     * Returns depth count for each diploid allele at the given taxon and site.
     * The array of depths is sized as determined by NUMBER_NUCLEOTIDE_ALLELES (6), and it is ordered as in
     * NUCLEOTIDE_ALLELES (A, C, G, T, +(insertion), - (deletion))
     *
     * The depths from 0 to 127 are recorded exactly.  Depths above 127 are stored as approximate logs (equation), and
     * they are represented by negative numbers.
     *
     * @param taxon taxon
     * @param site site
     *
     * @return array of counts
     */
    public byte[] getDepthForAlleles(int taxon, int site);
}
