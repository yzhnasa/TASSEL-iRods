/*
 * BitStorageNew
 */
package net.maizegenetics.dna.snp.bit;

import net.maizegenetics.util.BitSet;

/**
 * Interface provides genotypes in a binary fashion to be used in rapid
 * computation. See the package descriptions
 * {@link net.maizegenetics.dna.snp.bit} for more information on how bits
 * encoding is used throughout TASSEL.
 *
 * @author Terry Casstevens
 * @author Ed Buckler
 */
public interface BitStorageNew {

    /**
     * Returns sequence of true/false values indicating whether taxon at each
     * site matches a specific allele.
     *
     * @param taxon taxon
     *
     * @return sequence of true/false values.
     */
    public BitSet allelePresenceForAllSites(int taxon);

    /**
     * Returns sequence of true/false values indicating whether site at each
     * taxon matches a specific allele.
     *
     * @param site site
     *
     * @return sequence of true/false values.
     */
    public BitSet allelePresenceForAllTaxa(int site);

    /**
     * Returns sequence of true/false values indicating whether taxon at sites
     * (in given blocks, 64 sites per block including start block but excluding
     * end block) matches a specific allele.
     *
     * @param taxon taxon
     * @param startBlock starting block
     * @param endBlock end block
     *
     * @return sequence of true/false values.
     */
    public long[] allelePresenceForSitesBlock(int taxon, int startBlock, int endBlock);

}
