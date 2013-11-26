/*
 * FilterAlignment
 */
package net.maizegenetics.dna.snp;

import net.maizegenetics.dna.snp.bit.BitStorage;
import net.maizegenetics.dna.snp.bit.DynamicBitStorage;
import net.maizegenetics.dna.snp.genotype.Genotype;
import net.maizegenetics.dna.snp.genotype.GenotypeBuilder;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.map.PositionListBuilder;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.taxa.Taxon;
import net.maizegenetics.util.BitSet;
import org.apache.log4j.Logger;

import java.util.*;

/**
 * All taxa and site filtering should be controlled through this class. It
 * essentially creates views of the baseAlignment
 *
 * @author terry
 */
public class FilterAlignment implements Alignment {

    private static final long serialVersionUID = -5197800047652332969L;
    private static final Logger myLogger = Logger.getLogger(FilterAlignment.class);
    private final boolean myIsTaxaFilter;
    private final boolean myIsSiteFilter;
    private final boolean myIsSiteFilterByRange;
    private final Alignment myBaseAlignment;
    private final TaxaList myTaxaList;
    private final int[] myTaxaRedirect;
    private final int[] mySiteRedirect;
    private final int myRangeStart;
    private final int myRangeEnd;
    private Chromosome[] myChromosomes;
    private int[] myChromosomeOffsets;
    private PositionList myPositionList;
    private final Genotype myGenotype;
    private final Map<ALLELE_SORT_TYPE, BitStorage> myBitStorage = new EnumMap<ALLELE_SORT_TYPE, BitStorage>(ALLELE_SORT_TYPE.class);

    private FilterAlignment(Alignment a, TaxaList subList, int[] taxaRedirect, FilterAlignment original) {

        myTaxaList = subList;

        if (myTaxaList.getTaxaCount() != taxaRedirect.length) {
            throw new IllegalArgumentException("FilterAlignment: init: subList should be same size as taxaRedirect.");
        }

        myIsTaxaFilter = true;
        myBaseAlignment = a;
        myTaxaRedirect = taxaRedirect;

        if (original == null) {
            myIsSiteFilter = false;
            myIsSiteFilterByRange = false;
            mySiteRedirect = null;
            myRangeStart = -1;
            myRangeEnd = -1;
            myChromosomes = myBaseAlignment.chromosomes();
            myChromosomeOffsets = myBaseAlignment.chromosomesOffsets();
        } else {
            myIsSiteFilter = original.isSiteFilter();
            myIsSiteFilterByRange = original.isSiteFilterByRange();
            mySiteRedirect = original.getSiteRedirect();
            myRangeStart = original.getRangeStart();
            myRangeEnd = original.getRangeEnd();
            myChromosomes = original.chromosomes();
            myChromosomeOffsets = original.chromosomesOffsets();
        }

        if (myIsSiteFilter) {
            myGenotype = GenotypeBuilder.getFilteredInstance(myBaseAlignment.genotypeMatrix(), numberOfTaxa(), myTaxaRedirect, numberOfSites(), mySiteRedirect);
        } else {
            myGenotype = GenotypeBuilder.getFilteredInstance(myBaseAlignment.genotypeMatrix(), numberOfTaxa(), myTaxaRedirect, numberOfSites(), myRangeStart, myRangeEnd);
        }

    }

    /**
     * This returns FilterAlignment with only specified subTaxaList. Defaults to
     * retain unknown taxa.
     *
     * @param a alignment
     * @param subTaxaList subset id group
     *
     * @return filter alignment
     */
    public static Alignment getInstance(Alignment a, TaxaList subTaxaList) {
        return getInstance(a, subTaxaList, true);
    }

    /**
     * This returns FilterAlignment with only specified subTaxaList. If
     * retainUnknownTaxa is true then Alignment will return unknown values for
     * missing taxa.
     *
     * @param a alignment
     * @param subTaxaList subset id group
     * @param retainUnknownTaxa whether to retain unknown taxa
     *
     * @return filter alignment
     */
    public static Alignment getInstance(Alignment a, TaxaList subTaxaList, boolean retainUnknownTaxa) {

        Alignment baseAlignment = a;
        FilterAlignment original = null;
        if (baseAlignment instanceof FilterAlignment) {
            original = (FilterAlignment) a;
            baseAlignment = ((FilterAlignment) a).getBaseAlignment();
        }

        List<Integer> taxaRedirectList = new ArrayList<Integer>();
        List<Taxon> idList = new ArrayList<Taxon>();
        boolean noNeedToFilter = true;
        if (subTaxaList.getTaxaCount() != a.numberOfTaxa()) {
            noNeedToFilter = false;
        }
        for (int i = 0, n = subTaxaList.getTaxaCount(); i < n; i++) {
            List<Integer> ion = a.taxa().getIndicesMatchingTaxon(subTaxaList.get(i));

            if ((ion.size() != 1) || (ion.get(0) != i)) {
                noNeedToFilter = false;
            }

            if (ion.isEmpty()) {
                if (retainUnknownTaxa) {
                    taxaRedirectList.add(-1);
                    idList.add(subTaxaList.get(i));
                }
            } else {
                for (int x = 0; x < ion.size(); x++) {
                    if (a instanceof FilterAlignment) {
                        taxaRedirectList.add(((FilterAlignment) a).translateTaxon(ion.get(x)));
                    } else {
                        taxaRedirectList.add(ion.get(x));
                    }
                    idList.add(a.taxa().get(ion.get(x)));
                }
            }
        }

        if (noNeedToFilter) {
            return a;
        }

        int[] taxaRedirect = new int[taxaRedirectList.size()];
        for (int j = 0, n = taxaRedirectList.size(); j < n; j++) {
            taxaRedirect[j] = (int) taxaRedirectList.get(j);
        }

        TaxaList resultTaxaList = new TaxaListBuilder().addAll(idList).build();

        return new FilterAlignment(baseAlignment, resultTaxaList, taxaRedirect, original);

    }

    /**
     * Removes specified IDs.
     *
     * @param a alignment to filter
     * @param subTaxaList specified IDs
     *
     * @return Filtered Alignment
     */
    public static Alignment getInstanceRemoveIDs(Alignment a, TaxaList subTaxaList) {

        TaxaListBuilder result = new TaxaListBuilder();
        TaxaList current = a.taxa();
        for (int i = 0, n = current.getTaxaCount(); i < n; i++) {
            if (subTaxaList.getIndicesMatchingTaxon(current.get(i)).isEmpty()) {
                result.add(current.get(i));
            }
        }
        return FilterAlignment.getInstance(a, result.build());

    }

    /**
     * Constructor
     *
     * @param a base alignment
     * @param startSite start site (included)
     * @param endSite end site (included)
     */
    private FilterAlignment(Alignment a, int startSite, int endSite, FilterAlignment original) {

        myTaxaList = original == null ? a.taxa() : original.taxa();

        if (startSite > endSite) {
            throw new IllegalArgumentException("FilterAlignment: init: start site: " + startSite + " is larger than end site: " + endSite);
        }

        if ((startSite < 0) || (startSite > a.numberOfSites() - 1)) {
            throw new IllegalArgumentException("FilterAlignment: init: start site: " + startSite + " is out of range.");
        }

        if ((endSite < 0) || (endSite > a.numberOfSites() - 1)) {
            throw new IllegalArgumentException("FilterAlignment: init: end site: " + endSite + " is out of range.");
        }

        myBaseAlignment = a;
        myIsSiteFilterByRange = true;
        myIsSiteFilter = false;
        myRangeStart = startSite;
        myRangeEnd = endSite;
        mySiteRedirect = null;
        getLociFromBase();

        if (original == null) {
            myIsTaxaFilter = false;
            myTaxaRedirect = null;
        } else {
            myIsTaxaFilter = original.isTaxaFilter();
            myTaxaRedirect = original.getTaxaRedirect();
        }

        if (myIsSiteFilter) {
            myGenotype = GenotypeBuilder.getFilteredInstance(myBaseAlignment.genotypeMatrix(), numberOfTaxa(), myTaxaRedirect, numberOfSites(), mySiteRedirect);
        } else {
            myGenotype = GenotypeBuilder.getFilteredInstance(myBaseAlignment.genotypeMatrix(), numberOfTaxa(), myTaxaRedirect, numberOfSites(), myRangeStart, myRangeEnd);
        }

    }

    /**
     * Constructor
     *
     * @param a base alignment
     * @param subSites site to include
     */
    private FilterAlignment(Alignment a, int[] subSites, FilterAlignment original) {

        myTaxaList = original == null ? a.taxa() : original.taxa();

        myBaseAlignment = a;
        myIsSiteFilter = true;
        myIsSiteFilterByRange = false;
        if ((subSites == null) || (subSites.length == 0)) {
            mySiteRedirect = new int[0];
        } else {
            mySiteRedirect = new int[subSites.length];
            Arrays.sort(subSites);
            System.arraycopy(subSites, 0, mySiteRedirect, 0, subSites.length);
        }
        myRangeStart = -1;
        myRangeEnd = -1;
        getLociFromBase();

        if (original == null) {
            myIsTaxaFilter = false;
            myTaxaRedirect = null;
        } else {
            myIsTaxaFilter = original.isTaxaFilter();
            myTaxaRedirect = original.getTaxaRedirect();
        }

        if (myIsSiteFilter) {
            myGenotype = GenotypeBuilder.getFilteredInstance(myBaseAlignment.genotypeMatrix(), numberOfTaxa(), myTaxaRedirect, numberOfSites(), mySiteRedirect);
        } else {
            myGenotype = GenotypeBuilder.getFilteredInstance(myBaseAlignment.genotypeMatrix(), numberOfTaxa(), myTaxaRedirect, numberOfSites(), myRangeStart, myRangeEnd);
        }

    }

    public static FilterAlignment getInstance(Alignment a, int[] subSites) {

        if (a instanceof FilterAlignment) {
            FilterAlignment original = (FilterAlignment) a;
            Alignment baseAlignment = ((FilterAlignment) a).getBaseAlignment();
            if (original.isSiteFilter()) {
                int[] newSubSites = new int[subSites.length];
                for (int i = 0; i < subSites.length; i++) {
                    newSubSites[i] = original.translateSite(subSites[i]);
                }
                return new FilterAlignment(baseAlignment, newSubSites, original);
            } else if (original.isSiteFilterByRange()) {
                int[] newSubSites = new int[subSites.length];
                for (int i = 0; i < subSites.length; i++) {
                    newSubSites[i] = original.translateSite(subSites[i]);
                }
                return new FilterAlignment(baseAlignment, newSubSites, original);
            } else if (original.isTaxaFilter()) {
                return new FilterAlignment(baseAlignment, subSites, original);
            } else {
                throw new IllegalStateException("FilterAlignment: getInstance: original not in known state.");
            }
        } else {
            return new FilterAlignment(a, subSites, null);
        }

    }

    public static FilterAlignment getInstance(Alignment a, String[] siteNamesToKeep) {

        Arrays.sort(siteNamesToKeep);
        int[] temp = new int[siteNamesToKeep.length];
        int count = 0;
        for (int i = 0, n = a.numberOfSites(); i < n; i++) {
            if (Arrays.binarySearch(siteNamesToKeep, a.siteName(i)) >= 0) {
                temp[count++] = i;
                if (count == siteNamesToKeep.length) {
                    break;
                }
            }
        }

        int[] result = null;
        if (count == siteNamesToKeep.length) {
            result = temp;
        } else {
            result = new int[count];
            System.arraycopy(temp, 0, result, 0, count);
        }
        return getInstance(a, result);

    }

    public static FilterAlignment getInstanceRemoveSiteNames(Alignment a, String[] siteNamesToRemove) {

        Arrays.sort(siteNamesToRemove);
        int[] temp = new int[a.numberOfSites()];
        int count = 0;
        for (int i = 0, n = a.numberOfSites(); i < n; i++) {
            if (Arrays.binarySearch(siteNamesToRemove, a.siteName(i)) < 0) {
                temp[count++] = i;
            }
        }

        int[] result = null;
        if (count == temp.length) {
            result = temp;
        } else {
            result = new int[count];
            System.arraycopy(temp, 0, result, 0, count);
        }
        return getInstance(a, result);

    }

    public static FilterAlignment getInstance(Alignment a, String chromosome, int startPhysicalPos, int endPhysicalPos) {
        return getInstance(a, a.chromosome(chromosome), startPhysicalPos, endPhysicalPos);
    }

    public static FilterAlignment getInstance(Alignment a, Chromosome chromosome, int startPhysicalPos, int endPhysicalPos) {

        int startSite = a.siteOfPhysicalPosition(startPhysicalPos, chromosome);
        if (startSite < 0) {
            startSite = -(startSite + 1);
        }

        int endSite = a.siteOfPhysicalPosition(endPhysicalPos, chromosome);
        if (endSite < 0) {
            endSite = -(endSite + 2);
        }

        if (startSite > endSite) {
            myLogger.warn("getInstance: start site: " + startSite + " from physical pos: " + startPhysicalPos + " is larger than end site: " + endSite + " from physical pos: " + endPhysicalPos);
            return null;
        }

        return getInstance(a, startSite, endSite);

    }

    public static FilterAlignment getInstance(Alignment a, Chromosome chromosome) {
        int[] endStart = a.startAndEndOfChromosome(chromosome);
        return getInstance(a, endStart[0], endStart[1] - 1);
    }

    /**
     * Factory method that returns a FilterAlignment viewing sites between start
     * site and end site inclusive.
     *
     * @param a alignment
     * @param startSite start site
     * @param endSite end site
     *
     * @return Filter Alignment
     */
    public static FilterAlignment getInstance(Alignment a, int startSite, int endSite) {

        if (a instanceof FilterAlignment) {
            FilterAlignment original = (FilterAlignment) a;
            Alignment baseAlignment = ((FilterAlignment) a).getBaseAlignment();
            if (original.isSiteFilter()) {
                int[] subSites = new int[endSite - startSite + 1];
                int[] originalSites = original.getSiteRedirect();
                for (int i = startSite; i <= endSite; i++) {
                    subSites[i - startSite] = originalSites[i];
                }
                return new FilterAlignment(baseAlignment, subSites, original);
            } else if (original.isSiteFilterByRange()) {
                return new FilterAlignment(baseAlignment, original.translateSite(startSite), original.translateSite(endSite), original);
            } else if (original.isTaxaFilter()) {
                return new FilterAlignment(baseAlignment, startSite, endSite, original);
            } else {
                throw new IllegalStateException("FilterAlignment: getInstance: original not in known state.");
            }
        } else {
            return new FilterAlignment(a, startSite, endSite, null);
        }

    }

    @Override
    public byte genotype(int taxon, int site) {
        return myGenotype.genotype(taxon, site);
    }

    @Override
    public byte[] genotypeRange(int taxon, int startSite, int endSite) {
        return myGenotype.genotypeRange(taxon, startSite, endSite);
    }

    @Override
    public byte genotype(int taxon, Chromosome chromosome, int physicalPosition) {
        return myGenotype.genotype(taxon, myPositionList.siteOfPhysicalPosition(physicalPosition, chromosome));
    }

    public int translateSite(int site) {

        if (myIsSiteFilterByRange) {
            return site + myRangeStart;
        } else if (myIsSiteFilter) {
            return mySiteRedirect[site];
        } else {
            return site;
        }

    }

    /**
     * Returns site of this FilterAlignment based on given site from embedded
     * Alignment.
     *
     * @param site site
     * @return site in this alignment
     */
    public int reverseTranslateSite(int site) {

        if (myIsSiteFilterByRange) {
            return site - myRangeStart;
        } else if (myIsSiteFilter) {
            return Arrays.binarySearch(mySiteRedirect, site);
        } else {
            return site;
        }

    }

    /**
     * Returns sites from original alignment that are viewable (not filtered) by
     * this filter alignment.
     *
     * @return list of sites
     */
    public int[] getBaseSitesShown() {
        int numSites = numberOfSites();
        int[] result = new int[numSites];
        for (int i = 0; i < numSites; i++) {
            result[i] = translateSite(i);
        }
        return result;
    }

    public int translateTaxon(int taxon) {

        if (myIsTaxaFilter) {
            return myTaxaRedirect[taxon];
        } else {
            return taxon;
        }

    }

    private void getLociFromBase() {

        if ((!myIsSiteFilter) && (!myIsSiteFilterByRange)) {
            myChromosomes = myBaseAlignment.chromosomes();
            myChromosomeOffsets = myBaseAlignment.chromosomesOffsets();
            return;
        }

        int numSites = numberOfSites();
        List<Chromosome> chromosomes = new ArrayList<Chromosome>();
        List<Integer> offsets = new ArrayList<Integer>();
        for (int i = 0; i < numSites; i++) {
            Chromosome current = chromosome(i);
            if (!chromosomes.contains(current)) {
                chromosomes.add(current);
                offsets.add(i);
            }
        }

        myChromosomes = new Chromosome[chromosomes.size()];
        chromosomes.toArray(myChromosomes);

        myChromosomeOffsets = new int[offsets.size()];
        for (int i = 0; i < offsets.size(); i++) {
            myChromosomeOffsets[i] = (Integer) offsets.get(i);
        }

    }

    @Override
    public int indelSize(int site) {
        return myBaseAlignment.indelSize(translateSite(site));
    }

    @Override
    public Chromosome chromosome(int site) {
        return myBaseAlignment.chromosome(translateSite(site));
    }

    @Override
    public int chromosomalPosition(int site) {
        return myBaseAlignment.chromosomalPosition(translateSite(site));
    }

    @Override
    public String chromosomeName(int site) {
        return myBaseAlignment.chromosomeName(translateSite(site));
    }

    @Override
    public Chromosome chromosome(String name) {
        return myBaseAlignment.chromosome(name);
    }

    @Override
    public Chromosome[] chromosomes() {
        return myChromosomes;
    }

    @Override
    public int numChromosomes() {
        return myChromosomes.length;
    }

    @Override
    public int[] chromosomesOffsets() {
        return myChromosomeOffsets;
    }

    @Override
    public int[] startAndEndOfChromosome(Chromosome chromosome) {
        for (int i = 0; i < numChromosomes(); i++) {
            if (chromosome.equals(myChromosomes[i])) {
                int end = 0;
                if (i == numChromosomes() - 1) {
                    end = numberOfSites();
                } else {
                    end = myChromosomeOffsets[i + 1];
                }
                return new int[]{myChromosomeOffsets[i], end};
            }
        }
        throw new IllegalArgumentException("FilterAlignment: getStartAndEndOfLocus: this locus not defined: " + chromosome.getName());
    }

    @Override
    public float[][] siteScores() {

        if (!myBaseAlignment.hasSiteScores()) {
            return null;
        }

        int numSites = numberOfSites();
        int numSeqs = numberOfTaxa();
        float[][] result = new float[numSeqs][numSites];
        for (int i = 0; i < numSites; i++) {
            for (int j = 0; j < numSeqs; j++) {
                int taxaIndex = translateTaxon(j);
                if (taxaIndex == -1) {
                    result[j][i] = -9;
                } else {
                    result[j][i] = myBaseAlignment.siteScore(taxaIndex, translateSite(i));
                }
            }
        }

        return result;

    }

    @Override
    public byte referenceGenotype(int site) {
        return myBaseAlignment.referenceGenotype(translateSite(site));
    }

    @Override
    public int numberOfSites() {

        if (myIsSiteFilterByRange) {
            return myRangeEnd - myRangeStart + 1;
        } else if (myIsSiteFilter) {
            return mySiteRedirect.length;
        } else {
            return myBaseAlignment.numberOfSites();
        }

    }

    @Override
    public String siteName(int site) {
        return myBaseAlignment.siteName(translateSite(site));
    }

    @Override
    public boolean hasReference() {
        return myBaseAlignment.hasReference();
    }

    @Override
    public boolean isIndel(int site) {
        return myBaseAlignment.isIndel(translateSite(site));
    }

    @Override
    public int siteOfPhysicalPosition(int physicalPosition, Chromosome chromosome) {
        int temp = myBaseAlignment.siteOfPhysicalPosition(physicalPosition, chromosome);
        if (temp < 0) {
            temp = -(temp + 1);
            return -(reverseTranslateSite(temp) + 1);
        }
        return reverseTranslateSite(temp);
    }

    @Override
    public int siteOfPhysicalPosition(int physicalPosition, Chromosome chromosome, String snpName) {
        int temp = myBaseAlignment.siteOfPhysicalPosition(physicalPosition, chromosome, snpName);
        if (temp < 0) {
            temp = -(temp + 1);
            return -(reverseTranslateSite(temp) + 1);
        }
        return reverseTranslateSite(temp);
    }

    public Alignment getBaseAlignment() {
        return myBaseAlignment;
    }

    public boolean isTaxaFilter() {
        return myIsTaxaFilter;
    }

    public boolean isSiteFilter() {
        return myIsSiteFilter;
    }

    public boolean isSiteFilterByRange() {
        return myIsSiteFilterByRange;
    }

    protected int[] getTaxaRedirect() {
        return myTaxaRedirect;
    }

    protected int[] getSiteRedirect() {
        return mySiteRedirect;
    }

    protected int getRangeStart() {
        return myRangeStart;
    }

    protected int getRangeEnd() {
        return myRangeEnd;
    }

    @Override
    public byte[] genotypeArray(int taxon, int site) {
        return myGenotype.genotypeArray(taxon, site);
    }

    @Override
    public byte[] genotypeAllSites(int taxon) {
        return myGenotype.genotypeAllSites(taxon);
    }

    @Override
    public BitSet allelePresenceForAllSites(int taxon, int alleleNumber) {
        return bitStorage(ALLELE_SORT_TYPE.Frequency).allelePresenceForAllSites(taxon, alleleNumber);
    }

    @Override
    public BitSet allelePresenceForAllTaxa(int site, int alleleNumber) {
        return bitStorage(ALLELE_SORT_TYPE.Frequency).allelePresenceForAllTaxa(site, alleleNumber);
    }

    @Override
    public long[] allelePresenceForSitesBlock(int taxon, int alleleNumber, int startBlock, int endBlock) {
        return bitStorage(ALLELE_SORT_TYPE.Frequency).allelePresenceForSitesBlock(taxon, alleleNumber, startBlock, endBlock);
    }

    @Override
    public String genotypeAsString(int taxon, int site) {
        return myGenotype.genotypeAsString(taxon, site);
    }

    @Override
    public String[] genotypeAsStringArray(int taxon, int site) {
        return myGenotype.genotypeAsStringArray(taxon, site);
    }

    @Override
    public byte[] referenceGenotypeForAllSites() {
        if ((myIsSiteFilterByRange) || (myIsSiteFilter)) {
            byte[] result = new byte[numberOfSites()];
            for (int i = 0, n = numberOfSites(); i < n; i++) {
                result[i] = referenceGenotype(i);
            }
            return result;
        } else {
            return myBaseAlignment.referenceGenotypeForAllSites();
        }
    }

    @Override
    public boolean isHeterozygous(int taxon, int site) {
        return myGenotype.isHeterozygous(taxon, site);
    }

    @Override
    public int[] physicalPositions() {
        if ((myIsSiteFilterByRange) || (myIsSiteFilter)) {
            int numSites = numberOfSites();
            int[] result = new int[numSites];
            for (int i = 0; i < numSites; i++) {
                result[i] = chromosomalPosition(i);
            }
            return result;
        } else {
            return myBaseAlignment.physicalPositions();
        }
    }

    @Override
    public TaxaList taxa() {
        return myTaxaList;
    }

    @Override
    public int numberOfTaxa() {
        return myTaxaList.getTaxaCount();
    }

    @Override
    public float siteScore(int taxon, int site) {
        int taxaIndex = translateTaxon(taxon);
        if (taxaIndex == -1) {
            return Float.NaN;
        } else {
            return myBaseAlignment.siteScore(taxaIndex, translateSite(site));
        }
    }

    @Override
    public boolean hasSiteScores() {
        return myBaseAlignment.hasSiteScores();
    }

    @Override
    public SITE_SCORE_TYPE siteScoreType() {
        return myBaseAlignment.siteScoreType();
    }

    @Override
    public String genomeVersion() {
        return myBaseAlignment.genomeVersion();
    }

    @Override
    public boolean isPositiveStrand(int site) {
        return myBaseAlignment.isPositiveStrand(translateSite(site));
    }

    @Override
    public boolean isPhased() {
        return myGenotype.isPhased();
    }

    @Override
    public boolean retainsRareAlleles() {
        return myGenotype.retainsRareAlleles();
    }

    @Override
    public String[][] alleleDefinitions() {
        return alleleDefinitions();
    }

    @Override
    public String[] alleleDefinitions(int site) {
        return alleleDefinitions(site);
    }

    @Override
    public String genotypeAsString(int site, byte value) {
        return myGenotype.genotypeAsString(site, value);
    }

    @Override
    public int maxNumAlleles() {
        return myGenotype.maxNumAlleles();
    }

    @Override
    public byte[] alleles(int site) {
        return myGenotype.alleles(site);
    }

    @Override
    public int[][] allelesSortedByFrequency(int site) {
        return myGenotype.allelesSortedByFrequency(site);
    }

    @Override
    public double majorAlleleFrequency(int site) {
        return myGenotype.majorAlleleFrequency(site);
    }

    @Override
    public int heterozygousCount(int site) {
        return myGenotype.heterozygousCount(site);
    }

    @Override
    public boolean isPolymorphic(int site) {
        return myGenotype.isPolymorphic(site);
    }

    @Override
    public int totalGametesNonMissingForSite(int site) {
        return myGenotype.totalGametesNonMissingForSite(site);
    }

    @Override
    public int totalNonMissingForSite(int site) {
        return myGenotype.totalNonMissingForSite(site);
    }

    @Override
    public int minorAlleleCount(int site) {
        return myGenotype.minorAlleleCount(site);
    }

    @Override
    public double minorAlleleFrequency(int site) {
        return myGenotype.minorAlleleFrequency(site);
    }

    @Override
    public int majorAlleleCount(int site) {
        return myGenotype.majorAlleleCount(site);
    }

    @Override
    public byte majorAllele(int site) {
        return myGenotype.majorAllele(site);
    }

    @Override
    public byte minorAllele(int site) {
        return myGenotype.minorAllele(site);
    }

    @Override
    public Object[][] genosSortedByFrequency(int site) {
        return myGenotype.genosSortedByFrequency(site);
    }

    @Override
    public int totalGametesNonMissingForTaxon(int taxon) {
        return myGenotype.totalGametesNonMissingForTaxon(taxon);
    }

    @Override
    public int totalNonMissingForTaxon(int taxon) {
        return myGenotype.totalNonMissingForTaxon(taxon);
    }

    @Override
    public int heterozygousCountForTaxon(int taxon) {
        return myGenotype.heterozygousCountForTaxon(taxon);
    }

    @Override
    public byte[] depthForAlleles(int taxon, int site) {
        int taxaIndex = translateTaxon(taxon);
        if (taxaIndex == -1) {
            return null;
        } else {
            return myBaseAlignment.depthForAlleles(taxaIndex, translateSite(site));
        }
    }

    @Override
    public byte[] allelesBySortType(ALLELE_SORT_TYPE scope, int site) {
        if (scope == ALLELE_SORT_TYPE.Frequency) {
            return alleles(site);
        } else {
            return myBaseAlignment.allelesBySortType(scope, translateSite(site));
        }
    }

    @Override
    public BitSet allelePresenceForAllTaxaBySortType(ALLELE_SORT_TYPE type, int site, int alleleNumber) {
        return bitStorage(type).allelePresenceForAllTaxa(site, alleleNumber);
    }

    @Override
    public BitSet haplotypeAllelePresenceForAllSites(int taxon, boolean firstParent, int alleleNumber) {
        return bitStorage(ALLELE_SORT_TYPE.Frequency).haplotypeAllelePresenceForAllSites(taxon, firstParent, alleleNumber);
    }

    @Override
    public BitSet haplotypeAllelePresenceForAllTaxa(int site, boolean firstParent, int alleleNumber) {
        return bitStorage(ALLELE_SORT_TYPE.Frequency).haplotypeAllelePresenceForAllTaxa(site, firstParent, alleleNumber);
    }

    @Override
    public long[] haplotypeAllelePresenceForSitesBlock(int taxon, boolean firstParent, int alleleNumber, int startBlock, int endBlock) {
        return bitStorage(ALLELE_SORT_TYPE.Frequency).haplotypeAllelePresenceForSitesBlock(taxon, firstParent, alleleNumber, startBlock, endBlock);
    }

    @Override
    public String genotypeAsStringRange(int taxon, int startSite, int endSite) {
        return myGenotype.genotypeAsStringRange(taxon, startSite, endSite);
    }

    @Override
    public String genotypeAsStringRow(int taxon) {
        return myGenotype.genotypeAsStringRow(taxon);
    }

    @Override
    public byte[] referenceGenotypes(int startSite, int endSite) {
        if (!hasReference()) {
            return null;
        }

        byte[] result = new byte[endSite - startSite];
        for (int i = startSite; i < endSite; i++) {
            result[i] = referenceGenotype(i);
        }
        return result;
    }

    @Override
    public int chromosomeSiteCount(Chromosome chromosome) {
        int[] startEnd = startAndEndOfChromosome(chromosome);
        return startEnd[1] - startEnd[0];
    }

    @Override
    public boolean isAllPolymorphic() {
        return myGenotype.isAllPolymorphic();
    }

    @Override
    public String majorAlleleAsString(int site) {
        return myGenotype.majorAlleleAsString(site);
    }

    @Override
    public String minorAlleleAsString(int site) {
        return myGenotype.minorAlleleAsString(site);
    }

    @Override
    public byte[] minorAlleles(int site) {
        return myGenotype.minorAlleles(site);
    }

    @Override
    public String taxaName(int index) {
        return myTaxaList.getTaxaName(index);
    }

    @Override
    public Alignment[] compositeAlignments() {
        return new Alignment[]{this};
    }

    @Override
    public String diploidAsString(int site, byte value) {
        return myGenotype.diploidAsString(site, value);
    }

    @Override
    public Object[][] genoCounts() {
        return myGenotype.genoCounts();
    }

    @Override
    public Object[][] majorMinorCounts() {
        return myGenotype.majorMinorCounts();
    }

    @Override
    public BitStorage bitStorage(ALLELE_SORT_TYPE scopeType) {

        BitStorage result = myBitStorage.get(scopeType);
        if (result != null) {
            return result;
        }

        switch (scopeType) {
            case Frequency:
                result = new DynamicBitStorage(myGenotype, scopeType, myGenotype.majorAlleleForAllSites(), myGenotype.minorAlleleForAllSites());
                break;
            case Reference:
                result = DynamicBitStorage.getInstance(myGenotype, scopeType, referenceGenotypeForAllSites());
                break;
            default:
                myLogger.warn("getBitStorage: Unsupported type: " + scopeType);
                return null;
        }

        myBitStorage.put(scopeType, result);
        return result;
    }

    @Override
    public PositionList positionList() {
        if (myPositionList == null) {
            PositionListBuilder pLB = new PositionListBuilder();
            PositionList basePL = getBaseAlignment().positionList();
            for (int i = 0; i < numberOfSites(); i++) {
                pLB.add(basePL.get(translateSite(i)));
            }
            myPositionList = pLB.build();
        }
        return myPositionList;
    }

    @Override
    public Genotype genotypeMatrix() {
        return myGenotype;
    }
}
