package net.maizegenetics.pal.site;

import net.maizegenetics.pal.alignment.Locus;

import java.util.List;

/**
 * List of sites in the genome.
 *
 * @author Terry Casstevens and Ed Buckler
 */
public interface SiteList extends List<AnnotatedSite> {

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
     * Get SNP IDs.
     *
     * @return site names.
     */
    public String[] getSNPIDs();

    /**
     * Get SNP ID for specified site.
     *
     * @param site site
     * @return site name
     */
    public String getSNPID(int site);

    /**
     * Returns total number of sites of this alignment.
     *
     * @return number of sites
     */
    public int getSiteCount();

    /**
     * Return number of sites for given locus.
     *
     * @param locus locus
     *
     * @return number of sites
     */
    public int getLocusSiteCount(Locus locus);

    /**
     * Get the first (inclusive) and last (exclusive) site of the specified
     * locus in this alignment.
     *
     * @param locus locus
     *
     * @return first and last site
     */
    public int[] getStartAndEndOfLocus(Locus locus);

    /**
     * Returns the physical position at given site.
     *
     * @param site site
     *
     * @return physical position
     */
    public int getPositionInLocus(int site);

    /**
     * Return site of given physical position in locus. If the physical position
     * doesn't exist, (-(insertion point) - 1) is returned. If locus is not
     * found, an exception is thrown.
     *
     * @param physicalPosition physical position
     * @param locus locus. if null, the first locus is used.
     *
     * @return index
     */
    public int getSiteOfPhysicalPosition(int physicalPosition, Locus locus);

    /**
     * Return site of given physical position / SNP ID in locus. If the physical
     * position doesn't exist, (-(insertion point) - 1) is returned. If locus is
     * not found, an exception is thrown. This is to support multiple sites with
     * the same physical position but different SNP IDs.
     *
     * @param physicalPosition physical position
     * @param locus locus. if null, the first locus is used.
     * @param snpID SNP ID
     *
     * @return index
     */
    public int getSiteOfPhysicalPosition(int physicalPosition, Locus locus, String snpID);

    /**
     * Returns all physical positions.
     *
     * @return physical positions.
     */
    public int[] getPhysicalPositions();

    /**
     * Returns position type for given site. (eg. I=intron, E=exon, P=promoter,
     * 1=first, 2=second, 3=third, etc.)
     *
     * @param site site
     *
     * @return position type
     */
    public byte getPositionType(int site);

    /**
     * Returns position types. (eg. I=intron, E=exon, P=promoter, 1=first,
     * 2=second, 3=third, etc.)
     *
     * @return position types
     */
    public byte[] getPositionTypes();

    /**
     * Return Locus Name for given site.
     *
     * @param site site
     *
     * @return Locus Name
     */
    public String getLocusName(int site);

    /**
     * Return Locus for given site.
     *
     * @param site site
     *
     * @return Locus
     */
    public Locus getLocus(int site);

    /**
     * Return Locus with matching name. First to match will be returned.
     *
     * @param name name
     *
     * @return Locus
     */
    public Locus getLocus(String name);

    /**
     * Return all loci.
     *
     * @return loci
     */
    public Locus[] getLoci();

    /**
     * Return number of loci.
     *
     * @return number of loci
     */
    public int getNumLoci();

    /**
     * Returns starting site for each locus.
     *
     * @return starting site for each locus.
     */
    public int[] getLociOffsets();

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
     * Return most common allele at given site. Gap is included as state.
     * Heterozygous count one for each allele value. Homozygous counts two for
     * the allele value.
     *
     * @param site site
     *
     * @return most common allele
     */
    public byte getMajorAllele(int site);

    /**
     * Return most common allele at given site. Gap is included as state.
     * Heterozygous count one for each allele value. Homozygous counts two for
     * the allele value.
     *
     * @param site site
     *
     * @return most common allele as String
     */
    public String getMajorAlleleAsString(int site);

    /**
     * Return most common minor allele at given site. Gap is included as state.
     * Heterozygous count one for each allele value. Homozygous counts two for
     * the allele value.
     *
     * @param site site
     *
     * @return most common minor allele
     */
    public byte getMinorAllele(int site);

    /**
     * Return most common minor allele at given site. Gap is included as state.
     * Heterozygous count one for each allele value. Homozygous counts two for
     * the allele value.
     *
     * @param site site
     *
     * @return most common minor allele as String
     */
    public String getMinorAlleleAsString(int site);

    /**
     * Return all minor alleles at given site. Gap is included as state.
     * Heterozygous count one for each allele value. Homozygous counts two for
     * the allele value.
     *
     * @param site site
     *
     * @return all minor alleles
     */
    public byte[] getMinorAlleles(int site);

    /**
     * Returns all alleles at given site in order of frequency. Gap is included
     * as state. Heterozygous count one for each allele value. Homozygous counts
     * two for the allele value.
     *
     * @param site site
     *
     * @return all alleles
     */
    public byte[] getAlleles(int site);

    /**
     * Return frequency for most common minor allele at given site. Gap is
     * included as state. Heterozygous count one for each allele value.
     * Homozygous counts two for the allele value.
     *
     * @param site site
     *
     * @return frequency
     */
    public double getMinorAlleleFrequency(int site);

    /**
     * Return frequency for major allele at given site. Gap is included as
     * state. Heterozygous count one for each allele value. Homozygous counts
     * two for the allele value.
     *
     * @param site site
     *
     * @return frequency
     */
    public double getMajorAlleleFrequency(int site);

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
}
