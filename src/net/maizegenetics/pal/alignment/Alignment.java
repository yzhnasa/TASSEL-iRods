/*
 * Alignment
 */
package net.maizegenetics.pal.alignment;

import net.maizegenetics.pal.alignment.bit.BitStorage;
import net.maizegenetics.pal.alignment.genotype.Genotype;
import net.maizegenetics.pal.position.Chromosome;
import net.maizegenetics.pal.position.PositionList;
import net.maizegenetics.pal.taxa.TaxaList;
import net.maizegenetics.util.BitSet;

/**
 * This supports heterozygous diploid alignments.
 *
 * @author terry
 */
public interface Alignment extends Genotype {

    /**
     * This encoding is used to lump together allele values with frequencies too
     * low to be retained as one of the maximum number of alleles.
     */
    public static byte RARE_ALLELE = 0xE;
    public static byte RARE_DIPLOID_ALLELE = (byte) 0xEE;
    public static String RARE_ALLELE_STR = "Z";
    public static byte UNKNOWN_ALLELE = 0xF;
    public static byte UNKNOWN_DIPLOID_ALLELE = (byte) 0xFF;
    public static String UNKNOWN_ALLELE_STR = "N";
    public static String UNKNOWN_DIPLOID_ALLELE_STR = "N:N";
    public static char UNKNOWN_ALLELE_CHAR = 'N';

    public static enum SITE_SCORE_TYPE {

        None, MixedScoreTypes, QualityScore, ImputedProbablity, Dosage
    };

    /**
     * This defines the possible allele scope types.
     */
    public static enum ALLELE_SCOPE_TYPE {

        /**
         * This is the default where alleles are sorted by frequency. Same as
         * getAlleles().
         */
        Frequency,
        /**
         * This sorts alleles based on there depth value.
         */
        Depth,
        /**
         * This uses the allele frequency of a base/global Alignment determine
         * sort order of alleles. That Alignment is usually a superset.
         */
        Global_Frequency,
        /**
         * This sorts alleles based on the reference sequence.
         */
        Reference
    };

    /**
     * Returns the immutable Genotype matrix. Taxa and Positions are not part of
     * the matrix. This method is used for copying Alignments, when either the
     * Taxa or Positions have changed.
     *
     * @return genotype matrix
     */
    public Genotype getGenotypeMatrix();

    /**
     * Returns diploid values for given taxon, chromosome, and physical
     * position. The chromosome and physical position should map to an unique
     * site.
     *
     * @param taxon taxon
     * @param chromosome chromosome
     * @param physicalPosition physical position
     *
     * @return first four bits are the first allele value and the second four
     * bits are the second allele value.
     */
    public byte getBase(int taxon, Chromosome chromosome, int physicalPosition);

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

    /**
     * Return reference diploid allele values at given site.
     *
     * @param site site
     *
     * @return first four bits are the first allele value and the second four
     * bits are the second allele value.
     */
    public byte getReferenceAllele(int site);

    /**
     * Returns reference sequence of diploid allele values for given taxon in
     * specified range (end site not included). Each value in array contains
     * both diploid values. First four bits holds the first allele, and the
     * second four bits holds the second allele.
     *
     * @param startSite start site
     * @param endSite end site
     *
     * @return reference sequence of diploid allele values.
     */
    public byte[] getReference(int startSite, int endSite);

    /**
     * Returns reference sequence of diploid allele values. Each value in array
     * contains both diploid values. First four bits holds the first allele, and
     * the second four bits holds the second allele.
     *
     * @return reference sequence of diploid allele values.
     */
    public byte[] getReference();

    /**
     * Return whether this alignment has defined reference sequence.
     *
     * @return true if this alignment has reference sequence.
     */
    public boolean hasReference();

    /**
     * Get SNP ID for specified site.
     *
     * @param site site
     * @return site name
     */
    public String getSNPID(int site);

    /**
     * Return number of sites for given chromosome.
     *
     * @param chromosome chromosome
     *
     * @return number of sites
     */
    public int getChromosomeSiteCount(Chromosome chromosome);

    /**
     * Get the first (inclusive) and last (exclusive) site of the specified
     * chromosome in this alignment.
     *
     * @param chromosome chromosome
     *
     * @return first and last site
     */
    public int[] getStartAndEndOfChromosome(Chromosome chromosome);

    /**
     * Returns number of sequences (taxa).
     *
     * @return number of sequences
     */
    public int getSequenceCount();

    /**
     * Return the position list for the alignment.
     *
     * @return PositionList for all sites.
     */
    public PositionList getPositionList();

    /**
     * Returns the physical position at given site.
     *
     * @param site site
     *
     * @return physical position
     */
    public int getPositionInChromosome(int site);

    /**
     * Return site of given physical position in chromosome. If the physical
     * position doesn't exist, (-(insertion point) - 1) is returned. If
     * chromosome is not found, an exception is thrown.
     *
     * @param physicalPosition physical position
     * @param chromosome chromosome. if null, the first chromosome is used.
     *
     * @return index
     */
    public int getSiteOfPhysicalPosition(int physicalPosition, Chromosome chromosome);

    /**
     * Return site of given physical position / SNP ID in chromosome. If the
     * physical position doesn't exist, (-(insertion point) - 1) is returned. If
     * chromosome is not found, an exception is thrown. This is to support
     * multiple sites with the same physical position but different SNP IDs.
     *
     * @param physicalPosition physical position
     * @param chromosome chromosome. if null, the first chromosome is used.
     * @param snpID SNP ID
     *
     * @return index
     */
    public int getSiteOfPhysicalPosition(int physicalPosition, Chromosome chromosome, String snpID);

    /**
     * Returns all physical positions.
     *
     * @return physical positions.
     */
    public int[] getPhysicalPositions();

    /**
     * Return Chromosome Name for given site.
     *
     * @param site site
     *
     * @return Chromosome Name
     */
    public String getChromosomeName(int site);

    /**
     * Return Chromosome for given site.
     *
     * @param site site
     *
     * @return Chromosome
     */
    public Chromosome getChromosome(int site);

    /**
     * Return Chromosome with matching name. First to match will be returned.
     *
     * @param name name
     *
     * @return Chromosome
     */
    public Chromosome getChromosome(String name);

    /**
     * Return all chromosomes.
     *
     * @return chromosomes
     */
    public Chromosome[] getChromosomes();

    /**
     * Return number of chromosomes.
     *
     * @return number of chromosomes
     */
    public int getNumChromosomes();

    /**
     * Returns starting site for each chromosome.
     *
     * @return starting site for each chromosome.
     */
    public int[] getChromosomesOffsets();

    /**
     * Returns the site score of the given sequence and site.
     *
     * @param seq sequence index
     * @param site site
     *
     * @return site score.
     */
    public float getSiteScore(int seq, int site);

    /**
     * Returns the site scores.
     *
     * @return site scores.
     */
    public float[][] getSiteScores();

    /**
     * Returns true if this alignment has site scores.
     *
     * @return true if this alignment has site scores.
     */
    public boolean hasSiteScores();

    /**
     * Return what type of site scores this alignment has.
     *
     * @return site score type.
     */
    public Alignment.SITE_SCORE_TYPE getSiteScoreType();

    /**
     * Return size of indel at given site.
     *
     * @param site site
     *
     * @return indel size
     */
    public int getIndelSize(int site);

    /**
     * Returns whether give site is an indel.
     *
     * @param site site
     *
     * @return true if indel
     */
    public boolean isIndel(int site);

    /**
     * Return taxa list of this alignment.
     *
     * @return taxa list.
     */
    public TaxaList getTaxaList();

    /**
     * Return taxa name at given index.
     *
     * @param index
     *
     * @return taxa name
     */
    public String getTaxaName(int index);

    /**
     * Return full taxa name at given index.
     *
     * @param index
     * @return full taxa name
     */
    public String getFullTaxaName(int index);

    /**
     * Gets the Genome Assembly.
     *
     * @return the genome assembly.
     */
    public String getGenomeAssembly();

    /**
     * Return whether is positive strand at given site.
     *
     * @param site site
     *
     * @return whether is positive strand.
     */
    public boolean isPositiveStrand(int site);

    /**
     * Returns individual alignments within this alignment.
     *
     * @return list of alignments.
     */
    public Alignment[] getAlignments();

    /**
     * Returns Genetic Map of this Alignment if available.
     *
     * @return Genetic Map
     */
    public GeneticMap getGeneticMap();

    /**
     * Returns depth count for each diploid allele at the given taxon and site.
     *
     * @param taxon taxon
     * @param site site
     *
     * @return two counts
     */
    public byte[] getDepthForAlleles(int taxon, int site);

    /**
     * Returns all alleles at given site in order defined by scope.
     *
     * @param scope scope
     * @param site site
     *
     * @return alleles
     */
    public byte[] getAllelesByScope(Alignment.ALLELE_SCOPE_TYPE scope, int site);

    /**
     * Returns sequence of true/false values indicating whether site at each
     * taxon matches a specific allele (based on scope).
     *
     * @param scope scope
     * @param site site
     * @param alleleNumber allele number
     *
     * @return sequence of true/false values.
     */
    public BitSet getAllelePresenceForAllTaxaByScope(Alignment.ALLELE_SCOPE_TYPE scope, int site, int alleleNumber);

    /**
     * Returns BitStorage for this Genotype
     *
     * @param scopeType type
     *
     * @return BitStorage
     */
    public BitStorage getBitStorage(Alignment.ALLELE_SCOPE_TYPE scopeType);
}
