/*
 * BitStorage
 */
package net.maizegenetics.pal.alignment.bit;

import net.maizegenetics.util.BitSet;

/**
 *
 * @author terry
 */
public interface BitStorage {

    /**
     * Returns sequence of true/false values indicating whether taxon at each
     * site matches a specific allele (based on frequency). Allele number of
     * value 0 would be the major allele. Allele number of value 1 would be the
     * minor allele. Allele number of value 2 would be the third most frequent
     * allele value and so on.
     *
     * @param taxon taxon
     * @param alleleNumber allele number
     *
     * @return sequence of true/false values.
     */
    public BitSet getAllelePresenceForAllSites(int taxon, int alleleNumber);

    /**
     * Returns sequence of true/false values indicating whether site at each
     * taxon matches a specific allele (based on frequency). Allele number of
     * value 0 would be the major allele. Allele number of value 1 would be the
     * minor allele. Allele number of value 2 would be the third most frequent
     * allele value and so on.
     *
     * @param site site
     * @param alleleNumber allele number
     *
     * @return sequence of true/false values.
     */
    public BitSet getAllelePresenceForAllTaxa(int site, int alleleNumber);

    /**
     * Returns sequence of true/false values indicating whether taxon at sites
     * (in given blocks, 64 sites per block including start block but excluding
     * end block) matches a specific allele (based on frequency). Allele number
     * of value 0 would be the major allele. Allele number of value 1 would be
     * the minor allele. Allele number of value 2 would be the third most
     * frequent allele value and so on.
     *
     * @param taxon taxon
     * @param alleleNumber allele number
     * @param startBlock starting block
     * @param endBlock end block
     *
     * @return sequence of true/false values.
     */
    public long[] getAllelePresenceForSitesBlock(int taxon, int alleleNumber, int startBlock, int endBlock);

    /**
     * Returns sequence of true/false values indicating whether taxon at each
     * site for given parent matches a specific allele (based on frequency).
     * Allele number of value 0 would be the major allele. Allele number of
     * value 1 would be the minor allele. Allele number of value 2 would be the
     * third most frequent allele value and so on.
     *
     * @param taxon taxon
     * @param firstParent true for first parent (false for second parent)
     * @param alleleNumber allele number
     *
     * @return sequence of true/false values.
     */
    public BitSet getPhasedAllelePresenceForAllSites(int taxon, boolean firstParent, int alleleNumber);

    /**
     * Returns sequence of true/false values indicating whether site at each
     * taxon for given parent matches a specific allele (based on frequency).
     * Allele number of value 0 would be the major allele. Allele number of
     * value 1 would be the minor allele. Allele number of value 2 would be the
     * third most frequent allele value and so on.
     *
     * @param site site
     * @param firstParent true for first parent (false for second parent)
     * @param alleleNumber allele number
     *
     * @return sequence of true/false values.
     */
    public BitSet getPhasedAllelePresenceForAllTaxa(int site, boolean firstParent, int alleleNumber);

    /**
     * Returns sequence of true/false values indicating whether taxon at sites
     * (in given blocks, 64 sites per block including start block but excluding
     * end block) for given parent matches a specific allele (based on
     * frequency). Allele number of value 0 would be the major allele. Allele
     * number of value 1 would be the minor allele. Allele number of value 2
     * would be the third most frequent allele value and so on.
     *
     * @param taxon taxon
     * @param firstParent true for first parent (false for second parent)
     * @param alleleNumber allele number
     * @param startBlock starting block
     * @param endBlock end block
     *
     * @return sequence of true/false values.
     */
    public long[] getPhasedAllelePresenceForSitesBlock(int taxon, boolean firstParent, int alleleNumber, int startBlock, int endBlock);
}
