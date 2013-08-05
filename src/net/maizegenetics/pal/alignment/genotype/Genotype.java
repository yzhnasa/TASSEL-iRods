/*
 * Genotype
 */
package net.maizegenetics.pal.alignment.genotype;

import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.Locus;
import net.maizegenetics.pal.alignment.Locus;
import net.maizegenetics.util.ProgressListener;

/**
 *
 * @author terry
 */
public interface Genotype {

    /**
     * Returns diploid value (genotype) for a given taxon and site.
     *
     * @param taxon taxon
     * @param site site
     *
     * @return high four bits generally encode the more frequent allele and the
     * lower four bits encode the less frequent allele.
     */
    public byte getBase(int taxon, int site);

    /**
     * Returns diploid values for given taxon and site. Same values as
     * getBase(), except two values are already separated into two bytes.
     *
     * @param taxon taxon
     * @param site site
     *
     * @return first byte (index 0) holds first allele value in right-most four
     * bits. second byte (index 1) holds second allele value in right-most four
     * bits.
     */
    public byte[] getBaseArray(int taxon, int site);

    /**
     * Returns diploid values for given taxon, locus, and physical position. The
     * locus and physical position should map to an unique site.
     *
     * @param taxon taxon
     * @param locus locus
     * @param physicalPosition physical position
     *
     * @return first four bits are the first allele value and the second four
     * bits are the second allele value.
     */
    public byte getBase(int taxon, Locus locus, int physicalPosition);

    /**
     * Returns sequence of diploid allele values for given taxon in specified
     * range (end site excluded). Each value in array is what would be returned
     * by getBase().
     *
     * @param taxon taxon
     * @param startSite start site
     * @param endSite end site
     *
     * @return sequence of diploid allele values.
     */
    public byte[] getBaseRange(int taxon, int startSite, int endSite);

    /**
     * Returns sequence of diploid allele values for all sites for given taxon.
     * Each value in array is what would be returned by getBase().
     *
     * @param taxon taxon
     *
     * @return sequence of diploid allele values.
     */
    public byte[] getBaseRow(int taxon);

    /**
     * Returns string representation of diploid values returned by getBase() for
     * given taxon and site. The two allele values will be separated by a colon
     * (:) delimiter.
     *
     * @param taxon taxon
     * @param site site
     *
     * @return string representation of diploid values.
     */
    public String getBaseAsString(int taxon, int site);

    /**
     * Returns string representation of diploid alleles for given taxon in
     * specified range (end site excluded). Each value in string is what would
     * be returned by getBaseAsString().
     *
     * @param taxon taxon
     * @param startSite start site
     * @param endSite end site
     *
     * @return string representation of alleles in range
     */
    public String getBaseAsStringRange(int taxon, int startSite, int endSite);

    /**
     * Returns string representation of diploid alleles for given taxon for all
     * sites. Each value in string is what would be returned by
     * getBaseAsString().
     *
     * @param taxon taxon
     *
     * @return string representation of alleles
     */
    public String getBaseAsStringRow(int taxon);

    /**
     * Returns string representation of diploid values returned by
     * getBaseArray() for given taxon and site. Same two allele values as
     * getBaseAsString(), except already separated into two Strings.
     *
     * @param taxon taxon
     * @param site site
     *
     * @return string representations of diploid values.
     */
    public String[] getBaseAsStringArray(int taxon, int site);
    
    /**
     * Returns whether allele values at given taxon and site are heterozygous.
     * If two values returned by getBase() are different, this will return
     * false.
     *
     * @param taxon taxon
     * @param site site
     *
     * @return whether heterozygous
     */
    public boolean isHeterozygous(int taxon, int site);

    /**
     * Returns number of heterozygous taxa at given site.
     *
     * @param site site
     *
     * @return number of heterozygous taxa
     */
    public int getHeterozygousCount(int site);
    
    /**
     * Returns whether all sites are polymorphic.
     *
     * @return true if all sites are polymorphic.
     */
    public boolean isAllPolymorphic();

    /**
     * Return whether given site is polymorphic.
     *
     * @param site site
     *
     * @return true if given site is polymorphic.
     */
    public boolean isPolymorphic(int site);
    
    /**
     * Returns individual alignments within this alignment.
     *
     * @return list of alignments.
     */
    public Alignment[] getAlignments();
    
    /**
     * Returns whether this alignment is phased.
     *
     * @return true if phased.
     */
    public boolean isPhased();
    
    /**
     * Returns true if this Alignment retains rare alleles. If false, rare
     * alleles are recorded as unknown.
     *
     * @return whether rare alleles are retained.
     */
    public boolean retainsRareAlleles();

    /**
     * Returns allele values as strings for all sites. The first dimension of
     * the array indexes the sites. The second dimension indexes the allele
     * values for given site. The indices for the allele values are used as the
     * codes to store data. These codes (indices) are returned by the getBase()
     * methods. If only one array of allele values is returned, that is the
     * encoding for all sites.
     *
     * @return allele values for all sites.
     */
    public String[][] getAlleleEncodings();

    /**
     * Same as getAlleleEncodings() for only one site.
     *
     * @param site site
     *
     * @return allele values for given site.
     */
    public String[] getAlleleEncodings(int site);

    /**
     * Returns String representation of allele value at site.
     *
     * @param site site
     * @param value allele value
     *
     * @return String representation
     */
    public String getBaseAsString(int site, byte value);

    /**
     * Returns String representation of diploid allele value at site.
     *
     * @param site site
     * @param value diploid allele value
     *
     * @return String representation
     */
    public String getDiploidAsString(int site, byte value);

    /**
     * Return max number of alleles retained by this alignment.
     *
     * @return max number of alleles.
     */
    public int getMaxNumAlleles();

    /**
     * Returns max number of alleles retained by this alignment including rare
     * (if retained) and unknown (if phased).
     *
     * @return total number alleles
     */
    public int getTotalNumAlleles();

    /**
     * Returns total number of non-missing allele values for given site. This
     * can be twice the number of taxa, as diploid values are supported.
     *
     * @param site site
     * @return number of non-missing allele values.
     */
    public int getTotalGametesNotMissing(int site);

    /**
     * Returns total number of non-missing taxa for given site. Taxa are
     * considered missing only if both allele values are Unknown (N).
     *
     * @param site site
     *
     * @return number of non-missing taxa..
     */
    public int getTotalNotMissing(int site);

    /**
     * Returns the minor allele count for given site.
     *
     * @param site site
     * @return minor allele count
     */
    public int getMinorAlleleCount(int site);

    /**
     * Returns the major allele count for given site.
     *
     * @param site site
     * @return major allele count
     */
    public int getMajorAlleleCount(int site);

    /**
     * Returns counts of all diploid combinations from highest frequency to
     * lowest for whole alignment. Resulting double dimension array holds
     * diploids (Strings) in result[0]. And the counts are in result[1] (Longs).
     *
     * @return diploid counts.
     */
    public Object[][] getDiploidCounts();

    /**
     * Returns counts of all major/minor allele combinations from highest
     * frequency to lowest for whole alignment. Resulting double dimension array
     * holds major/minor allele (Strings) in result[0]. And the counts are in
     * result[1] (Longs).
     *
     * @return diploid counts.
     */
    public Object[][] getMajorMinorCounts();

    /**
     * Returns total number of non-missing allele values for given taxon. This
     * can be twice the number of sites, as diploid values are supported.
     *
     * @param taxon taxon
     *
     * @return number of non-missing allele values.
     */
    public int getTotalGametesNotMissingForTaxon(int taxon);

    /**
     * Returns number of heterozygous sites at given taxon.
     *
     * @param taxon taxon
     *
     * @return number of heterozygous sites
     */
    public int getHeterozygousCountForTaxon(int taxon);

    /**
     * Returns total number of non-missing sites for given taxon. Sites are
     * considered missing only if both allele values are Unknown (N).
     *
     * @param taxon taxon
     *
     * @return number of non-missing sites.
     */
    public int getTotalNotMissingForTaxon(int taxon);

    /**
     * Return whether alignment will execute quickly for site optimized
     * operations. SBitAlignment is obviously friendly. But so would
     * FilterAlignment is only sites have been filtered for example.
     *
     * @return whether optimized for site operations.
     */
    public boolean isSBitFriendly();

    /**
     * Return whether alignment will execute quickly for taxa optimized
     * operations. TBitAlignment is obviously friendly. But so would
     * FilterAlignment is only taxa have been filtered for example.
     *
     * @return whether optimized for taxa operations.
     */
    public boolean isTBitFriendly();

    /**
     * Optimizes this Alignment for Taxa based operations.
     */
    public void optimizeForTaxa(ProgressListener listener);

    /**
     * Optimizes this Alignment for Site based operations.
     */
    public void optimizeForSites(ProgressListener listener);

}
