/*
 * FilterAlignment
 */
package net.maizegenetics.pal.alignment;

import net.maizegenetics.pal.alignment.bit.BitStorage;
import net.maizegenetics.pal.alignment.genotype.Genotype;
import net.maizegenetics.pal.position.Chromosome;
import net.maizegenetics.pal.position.PositionList;
import net.maizegenetics.pal.position.PositionListBuilder;
import net.maizegenetics.pal.taxa.TaxaList;
import net.maizegenetics.pal.taxa.TaxaListBuilder;
import net.maizegenetics.pal.taxa.Taxon;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;
import net.maizegenetics.util.UnmodifiableBitSet;
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
            myChromosomes = myBaseAlignment.getChromosomes();
            myChromosomeOffsets = myBaseAlignment.getChromosomesOffsets();
        } else {
            myIsSiteFilter = original.isSiteFilter();
            myIsSiteFilterByRange = original.isSiteFilterByRange();
            mySiteRedirect = original.getSiteRedirect();
            myRangeStart = original.getRangeStart();
            myRangeEnd = original.getRangeEnd();
            myChromosomes = original.getChromosomes();
            myChromosomeOffsets = original.getChromosomesOffsets();
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
        if (subTaxaList.getTaxaCount() != a.getTaxaCount()) {
            noNeedToFilter = false;
        }
        for (int i = 0, n = subTaxaList.getTaxaCount(); i < n; i++) {
            List<Integer> ion = a.getTaxaList().getIndicesMatchingTaxon(subTaxaList.get(i));

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
                    idList.add(a.getTaxaList().get(ion.get(x)));
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
        TaxaList current = a.getTaxaList();
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

        myTaxaList = original == null ? a.getTaxaList() : original.getTaxaList();

        if (startSite > endSite) {
            throw new IllegalArgumentException("FilterAlignment: init: start site: " + startSite + " is larger than end site: " + endSite);
        }

        if ((startSite < 0) || (startSite > a.getSiteCount() - 1)) {
            throw new IllegalArgumentException("FilterAlignment: init: start site: " + startSite + " is out of range.");
        }

        if ((endSite < 0) || (endSite > a.getSiteCount() - 1)) {
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

    }

    /**
     * Constructor
     *
     * @param a base alignment
     * @param subSites site to include
     */
    private FilterAlignment(Alignment a, int[] subSites, FilterAlignment original) {

        myTaxaList = original == null ? a.getTaxaList() : original.getTaxaList();

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
        for (int i = 0, n = a.getSiteCount(); i < n; i++) {
            if (Arrays.binarySearch(siteNamesToKeep, a.getSNPID(i)) >= 0) {
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
        int[] temp = new int[a.getSiteCount()];
        int count = 0;
        for (int i = 0, n = a.getSiteCount(); i < n; i++) {
            if (Arrays.binarySearch(siteNamesToRemove, a.getSNPID(i)) < 0) {
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
        return getInstance(a, a.getChromosome(chromosome), startPhysicalPos, endPhysicalPos);
    }

    public static FilterAlignment getInstance(Alignment a, Chromosome chromosome, int startPhysicalPos, int endPhysicalPos) {

        int startSite = a.getSiteOfPhysicalPosition(startPhysicalPos, chromosome);
        if (startSite < 0) {
            startSite = -(startSite + 1);
        }

        int endSite = a.getSiteOfPhysicalPosition(endPhysicalPos, chromosome);
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
        int[] endStart = a.getStartAndEndOfChromosome(chromosome);
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
    public byte getBase(int taxon, int site) {
        int taxaIndex = translateTaxon(taxon);
        if (taxaIndex == -1) {
            return Alignment.UNKNOWN_ALLELE;
        } else {
            return myBaseAlignment.getBase(taxaIndex, translateSite(site));
        }
    }

    @Override
    public byte[] getBaseRange(int taxon, int startSite, int endSite) {

        int siteCount = endSite - startSite;
        byte[] result = new byte[siteCount];
        int taxaIndex = translateTaxon(taxon);
        if (taxaIndex == -1) {
            for (int i = 0; i < siteCount; i++) {
                result[i] = Alignment.UNKNOWN_ALLELE;
            }
            return result;
        } else {
            if (myIsSiteFilterByRange) {
                for (int i = 0; i < siteCount; i++) {
                    result[i] = myBaseAlignment.getBase(taxaIndex, startSite + myRangeStart + i);
                }
                return result;
            } else if (myIsSiteFilter) {
                for (int i = 0; i < siteCount; i++) {
                    result[i] = myBaseAlignment.getBase(taxaIndex, mySiteRedirect[startSite + i]);
                }
                return result;
            } else {
                for (int i = 0; i < siteCount; i++) {
                    result[i] = myBaseAlignment.getBase(taxaIndex, startSite + i);
                }
                return result;
            }
        }

    }

    @Override
    public byte getBase(int taxon, Chromosome locus, int physicalPosition) {
        int site = getSiteOfPhysicalPosition(physicalPosition, locus);
        int taxaIndex = translateTaxon(taxon);
        if (taxaIndex == -1) {
            return Alignment.UNKNOWN_ALLELE;
        } else {
            return myBaseAlignment.getBase(taxaIndex, site);
        }
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
        int numSites = getSiteCount();
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
            myChromosomes = myBaseAlignment.getChromosomes();
            myChromosomeOffsets = myBaseAlignment.getChromosomesOffsets();
            return;
        }

        int numSites = getSiteCount();
        List<Chromosome> chromosomes = new ArrayList<Chromosome>();
        List<Integer> offsets = new ArrayList<Integer>();
        for (int i = 0; i < numSites; i++) {
            Chromosome current = getChromosome(i);
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
    public int getIndelSize(int site) {
        return myBaseAlignment.getIndelSize(translateSite(site));
    }

    @Override
    public Chromosome getChromosome(int site) {
        return myBaseAlignment.getChromosome(translateSite(site));
    }

    @Override
    public int getPositionInChromosome(int site) {
        return myBaseAlignment.getPositionInChromosome(translateSite(site));
    }

    @Override
    public String getChromosomeName(int site) {
        return myBaseAlignment.getChromosomeName(translateSite(site));
    }

    @Override
    public Chromosome getChromosome(String name) {
        return myBaseAlignment.getChromosome(name);
    }

    @Override
    public Chromosome[] getChromosomes() {
        return myChromosomes;
    }

    @Override
    public int getNumChromosomes() {
        return myChromosomes.length;
    }

    @Override
    public int[] getChromosomesOffsets() {
        return myChromosomeOffsets;
    }

    @Override
    public int[] getStartAndEndOfChromosome(Chromosome chromosome) {
        for (int i = 0; i < getNumChromosomes(); i++) {
            if (chromosome.equals(myChromosomes[i])) {
                int end = 0;
                if (i == getNumChromosomes() - 1) {
                    end = getSiteCount();
                } else {
                    end = myChromosomeOffsets[i + 1];
                }
                return new int[]{myChromosomeOffsets[i], end};
            }
        }
        throw new IllegalArgumentException("FilterAlignment: getStartAndEndOfLocus: this locus not defined: " + chromosome.getName());
    }

    @Override
    public float[][] getSiteScores() {

        if (!myBaseAlignment.hasSiteScores()) {
            return null;
        }

        int numSites = getSiteCount();
        int numSeqs = getTaxaCount();
        float[][] result = new float[numSeqs][numSites];
        for (int i = 0; i < numSites; i++) {
            for (int j = 0; j < numSeqs; j++) {
                int taxaIndex = translateTaxon(j);
                if (taxaIndex == -1) {
                    result[j][i] = -9;
                } else {
                    result[j][i] = myBaseAlignment.getSiteScore(taxaIndex, translateSite(i));
                }
            }
        }

        return result;

    }

    @Override
    public byte getReferenceAllele(int site) {
        return myBaseAlignment.getReferenceAllele(translateSite(site));
    }

    @Override
    public int getSiteCount() {

        if (myIsSiteFilterByRange) {
            return myRangeEnd - myRangeStart + 1;
        } else if (myIsSiteFilter) {
            return mySiteRedirect.length;
        } else {
            return myBaseAlignment.getSiteCount();
        }

    }

    @Override
    public String getSNPID(int site) {
        return myBaseAlignment.getSNPID(translateSite(site));
    }

    @Override
    public String[] getSNPIDs() {

        if ((myIsSiteFilterByRange) || (myIsSiteFilter)) {
            int numSites = getSiteCount();
            String[] result = new String[numSites];
            for (int i = 0; i < numSites; i++) {
                result[i] = myBaseAlignment.getSNPID(translateSite(i));
            }

            return result;
        } else {
            return myBaseAlignment.getSNPIDs();
        }

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
    public int getSiteOfPhysicalPosition(int physicalPosition, Chromosome chromosome) {
        int temp = myBaseAlignment.getSiteOfPhysicalPosition(physicalPosition, chromosome);
        if (temp < 0) {
            temp = -(temp + 1);
            return -(reverseTranslateSite(temp) + 1);
        }
        return reverseTranslateSite(temp);
    }

    @Override
    public int getSiteOfPhysicalPosition(int physicalPosition, Chromosome chromosome, String snpID) {
        int temp = myBaseAlignment.getSiteOfPhysicalPosition(physicalPosition, chromosome, snpID);
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
    public byte[] getBaseArray(int taxon, int site) {
        int taxaIndex = translateTaxon(taxon);
        if (taxaIndex == -1) {
            return new byte[]{Alignment.UNKNOWN_ALLELE, Alignment.UNKNOWN_ALLELE};
        } else {
            return myBaseAlignment.getBaseArray(taxaIndex, translateSite(site));
        }
    }

    @Override
    public byte[] getBaseRow(int taxon) {
        int siteCount = getSiteCount();
        byte[] result = new byte[siteCount];
        int taxaIndex = translateTaxon(taxon);
        if (taxaIndex == -1) {
            for (int i = 0; i < siteCount; i++) {
                result[i] = Alignment.UNKNOWN_ALLELE;
            }
            return result;
        } else {
            if (myIsSiteFilterByRange) {
                for (int i = 0; i < siteCount; i++) {
                    result[i] = myBaseAlignment.getBase(taxaIndex, myRangeStart + i);
                }
                return result;
            } else if (myIsSiteFilter) {
                for (int i = 0; i < siteCount; i++) {
                    result[i] = myBaseAlignment.getBase(taxaIndex, mySiteRedirect[i]);
                }
                return result;
            } else {
                return myBaseAlignment.getBaseRow(taxaIndex);
            }
        }
    }

    @Override
    public BitSet getAllelePresenceForAllSites(int taxon, int alleleNumber) {
        if (myIsSiteFilter || myIsSiteFilterByRange) {
            throw new IllegalStateException("FilterAlignment: getAllelePresenceForAllSites: This Filter Alignment has had Sites removed.  You need to optimize for taxa before calling this.");
        } else {
            return myBaseAlignment.getAllelePresenceForAllSites(translateTaxon(taxon), alleleNumber);
        }
    }

    @Override
    public BitSet getAllelePresenceForAllTaxa(int site, int alleleNumber) {
        if (myIsTaxaFilter) {
            throw new IllegalStateException("FilterAlignment: getAllelePresenceForAllTaxa: This Filter Alignment has had Taxa removed.  You need to optimize for sites before calling this.");
        } else {
            return myBaseAlignment.getAllelePresenceForAllTaxa(translateSite(site), alleleNumber);
        }
    }

    @Override
    public long[] getAllelePresenceForSitesBlock(int taxon, int alleleNumber, int startBlock, int endBlock) {
        if (myIsSiteFilter || myIsSiteFilterByRange) {
            throw new IllegalStateException("FilterAlignment: getAllelePresenceForSitesBlock: This Filter Alignment has had Sites removed.  You need to optimize for taxa before calling this.");
        } else {
            return myBaseAlignment.getAllelePresenceForSitesBlock(translateTaxon(taxon), alleleNumber, startBlock, endBlock);
        }
    }

    @Override
    public String getBaseAsString(int taxon, int site) {
        int taxaIndex = translateTaxon(taxon);
        if (taxaIndex == -1) {
            return Alignment.UNKNOWN_ALLELE_STR;
        } else {
            return myBaseAlignment.getBaseAsString(taxaIndex, translateSite(site));
        }
    }

    @Override
    public String[] getBaseAsStringArray(int taxon, int site) {
        int taxaIndex = translateTaxon(taxon);
        if (taxaIndex == -1) {
            return new String[]{Alignment.UNKNOWN_ALLELE_STR, Alignment.UNKNOWN_ALLELE_STR};
        } else {
            return myBaseAlignment.getBaseAsStringArray(taxaIndex, translateSite(site));
        }
    }

    @Override
    public byte[] getReference() {
        if ((myIsSiteFilterByRange) || (myIsSiteFilter)) {
            byte[] result = new byte[getSiteCount()];
            for (int i = 0, n = getSiteCount(); i < n; i++) {
                result[i] = getReferenceAllele(i);
            }
            return result;
        } else {
            return myBaseAlignment.getReference();
        }
    }

    @Override
    public int getTotalNumAlleles() {
        return myBaseAlignment.getTotalNumAlleles();
    }

    @Override
    public boolean isHeterozygous(int taxon, int site) {
        int taxaIndex = translateTaxon(taxon);
        if (taxaIndex == -1) {
            return false;
        } else {
            return myBaseAlignment.isHeterozygous(taxaIndex, translateSite(site));
        }
    }

    @Override
    public int[] getPhysicalPositions() {
        if ((myIsSiteFilterByRange) || (myIsSiteFilter)) {
            int numSites = getSiteCount();
            int[] result = new int[numSites];
            for (int i = 0; i < numSites; i++) {
                result[i] = getPositionInChromosome(i);
            }
            return result;
        } else {
            return myBaseAlignment.getPhysicalPositions();
        }
    }

    @Override
    public TaxaList getTaxaList() {
        return myTaxaList;
    }

    @Override
    public int getSequenceCount() {
        return myTaxaList.getTaxaCount();
    }

    @Override
    public int getTaxaCount() {
        return myTaxaList.getTaxaCount();
    }

    @Override
    public float getSiteScore(int seq, int site) {
        int taxaIndex = translateTaxon(seq);
        if (taxaIndex == -1) {
            return Float.NaN;
        } else {
            return myBaseAlignment.getSiteScore(taxaIndex, translateSite(site));
        }
    }

    @Override
    public boolean hasSiteScores() {
        return myBaseAlignment.hasSiteScores();
    }

    @Override
    public SITE_SCORE_TYPE getSiteScoreType() {
        return myBaseAlignment.getSiteScoreType();
    }

    @Override
    public String getGenomeAssembly() {
        return myBaseAlignment.getGenomeAssembly();
    }

    @Override
    public boolean isPositiveStrand(int site) {
        return myBaseAlignment.isPositiveStrand(translateSite(site));
    }

    @Override
    public boolean isPhased() {
        return myBaseAlignment.isPhased();
    }

    @Override
    public GeneticMap getGeneticMap() {
        return myBaseAlignment.getGeneticMap();
    }

    @Override
    public boolean retainsRareAlleles() {
        return myBaseAlignment.retainsRareAlleles();
    }

    @Override
    public String[][] getAlleleEncodings() {
        String[][] encodings = myBaseAlignment.getAlleleEncodings();
        if (encodings.length == 1) {
            return encodings;
        } else if ((myIsSiteFilterByRange) || (myIsSiteFilter)) {
            int numSites = getSiteCount();
            String[][] result = new String[numSites][];
            for (int i = 0; i < numSites; i++) {
                result[i] = getAlleleEncodings(i);
            }
            return result;
        } else {
            return encodings;
        }
    }

    @Override
    public String[] getAlleleEncodings(int site) {
        return myBaseAlignment.getAlleleEncodings(translateSite(site));
    }

    @Override
    public String getBaseAsString(int site, byte value) {
        return myBaseAlignment.getBaseAsString(translateSite(site), value);
    }

    @Override
    public int getMaxNumAlleles() {
        return myBaseAlignment.getMaxNumAlleles();
    }

    @Override
    public byte[] getAlleles(int site) {
        if (myIsTaxaFilter) {
            int[][] alleles = getAllelesSortedByFrequency(site);
            int resultSize = alleles[0].length;
            byte[] result = new byte[resultSize];
            for (int i = 0; i < resultSize; i++) {
                result[i] = (byte) alleles[0][i];
            }
            return result;
        } else {
            return myBaseAlignment.getAlleles(translateSite(site));
        }
    }

    @Override
    public int[][] getAllelesSortedByFrequency(int site) {
        if (myIsTaxaFilter) {
            return AlignmentUtils.getAllelesSortedByFrequency(this, site);
        } else {
            return myBaseAlignment.getAllelesSortedByFrequency(translateSite(site));
        }
    }

    @Override
    public double getMajorAlleleFrequency(int site) {
        if (myIsTaxaFilter) {
            int[][] alleles = getAllelesSortedByFrequency(site);

            int numAlleles = alleles[0].length;
            if (numAlleles >= 1) {
                int totalNonMissing = 0;
                for (int i = 0; i < numAlleles; i++) {
                    totalNonMissing += alleles[1][i];
                }
                return (double) alleles[1][0] / (double) totalNonMissing;
            } else {
                return 0.0;
            }
        } else {
            return myBaseAlignment.getMajorAlleleFrequency(translateSite(site));
        }
    }

    @Override
    public int getHeterozygousCount(int site) {
        if (myIsTaxaFilter) {
            int result = 0;
            for (int i = 0, n = getTaxaCount(); i < n; i++) {
                if (isHeterozygous(i, site)) {
                    result++;
                }
            }
            return result;
        } else {
            return myBaseAlignment.getHeterozygousCount(translateSite(site));
        }
    }

    @Override
    public boolean isPolymorphic(int site) {
        if (myIsTaxaFilter) {
            byte first = Alignment.UNKNOWN_ALLELE;
            for (int i = 0, n = getTaxaCount(); i < n; i++) {
                byte[] current = getBaseArray(i, site);
                if (current[0] != Alignment.UNKNOWN_ALLELE) {
                    if (first == Alignment.UNKNOWN_ALLELE) {
                        first = current[0];
                    } else if (first != current[0]) {
                        return true;
                    }
                }
                if (current[1] != Alignment.UNKNOWN_ALLELE) {
                    if (first == Alignment.UNKNOWN_ALLELE) {
                        first = current[1];
                    } else if (first != current[1]) {
                        return true;
                    }
                }
            }

            return false;
        } else {
            return myBaseAlignment.isPolymorphic(translateSite(site));
        }
    }

    @Override
    public int getTotalGametesNotMissing(int site) {
        if (myIsTaxaFilter) {
            int result = 0;
            for (int i = 0, n = getTaxaCount(); i < n; i++) {
                byte[] current = getBaseArray(i, site);
                if (current[0] != Alignment.UNKNOWN_ALLELE) {
                    result++;
                }
                if (current[1] != Alignment.UNKNOWN_ALLELE) {
                    result++;
                }
            }
            return result;
        } else {
            return myBaseAlignment.getTotalGametesNotMissing(translateSite(site));
        }
    }

    @Override
    public int getTotalNotMissing(int site) {
        if (myIsTaxaFilter) {
            int result = 0;
            for (int i = 0, n = getTaxaCount(); i < n; i++) {
                byte[] current = getBaseArray(i, site);
                if (current[0] != Alignment.UNKNOWN_ALLELE) {
                    result++;
                }
                if (current[1] != Alignment.UNKNOWN_ALLELE) {
                    result++;
                }
            }
            return result;
        } else {
            return myBaseAlignment.getTotalNotMissing(translateSite(site));
        }
    }

    @Override
    public int getMinorAlleleCount(int site) {
        if (myIsTaxaFilter) {
            int[][] alleles = getAllelesSortedByFrequency(site);

            if (alleles[0].length >= 2) {
                return alleles[1][1];
            } else {
                return 0;
            }
        } else {
            return myBaseAlignment.getMinorAlleleCount(translateSite(site));
        }
    }

    @Override
    public double getMinorAlleleFrequency(int site) {
        if (myIsTaxaFilter) {
            int[][] alleles = getAllelesSortedByFrequency(site);

            int numAlleles = alleles[0].length;
            if (numAlleles >= 2) {
                int totalNonMissing = 0;
                for (int i = 0; i < numAlleles; i++) {
                    totalNonMissing += alleles[1][i];
                }
                return (double) alleles[1][1] / (double) totalNonMissing;
            } else {
                return 0.0;
            }
        } else {
            return myBaseAlignment.getMinorAlleleFrequency(translateSite(site));
        }
    }

    @Override
    public int getMajorAlleleCount(int site) {
        if (myIsTaxaFilter) {
            int[][] alleles = getAllelesSortedByFrequency(site);

            if (alleles[0].length >= 1) {
                return alleles[1][0];
            } else {
                return 0;
            }
        } else {
            return myBaseAlignment.getMajorAlleleCount(translateSite(site));
        }
    }

    @Override
    public byte getMajorAllele(int site) {
        if (myIsTaxaFilter) {
            int[][] alleles = getAllelesSortedByFrequency(site);

            if (alleles[0].length >= 1) {
                return (byte) alleles[0][0];
            } else {
                return Alignment.UNKNOWN_ALLELE;
            }
        } else {
            return myBaseAlignment.getMajorAllele(translateSite(site));
        }
    }

    @Override
    public byte getMinorAllele(int site) {
        if (myIsTaxaFilter) {
            int[][] alleles = getAllelesSortedByFrequency(site);

            if (alleles[0].length >= 1) {
                return (byte) alleles[0][0];
            } else {
                return Alignment.UNKNOWN_ALLELE;
            }
        } else {
            return myBaseAlignment.getMinorAllele(translateSite(site));
        }
    }

    @Override
    public Object[][] getDiploidsSortedByFrequency(int site) {
        if (myIsTaxaFilter) {
            return AlignmentUtils.getDiploidsSortedByFrequency(this, site);
        } else {
            return myBaseAlignment.getDiploidsSortedByFrequency(translateSite(site));
        }
    }

    @Override
    public int getTotalGametesNotMissingForTaxon(int taxon) {
        if (myIsSiteFilter || myIsSiteFilterByRange) {
            int result = 0;
            for (int i = 0, n = getSiteCount(); i < n; i++) {
                byte[] current = getBaseArray(taxon, i);
                if (current[0] != Alignment.UNKNOWN_ALLELE) {
                    result++;
                }
                if (current[1] != Alignment.UNKNOWN_ALLELE) {
                    result++;
                }
            }
            return result;
        } else {
            return myBaseAlignment.getTotalGametesNotMissingForTaxon(translateTaxon(taxon));
        }
    }

    @Override
    public int getTotalNotMissingForTaxon(int taxon) {
        if (myIsSiteFilter || myIsSiteFilterByRange) {
            int result = 0;
            for (int i = 0, n = getSiteCount(); i < n; i++) {
                byte[] current = getBaseArray(taxon, i);
                if (current[0] != Alignment.UNKNOWN_ALLELE) {
                    result++;
                }
                if (current[1] != Alignment.UNKNOWN_ALLELE) {
                    result++;
                }
            }
            return result;
        } else {
            return myBaseAlignment.getTotalNotMissingForTaxon(translateTaxon(taxon));
        }
    }

    @Override
    public int getHeterozygousCountForTaxon(int taxon) {
        if (myIsSiteFilter || myIsSiteFilterByRange) {
            int result = 0;
            for (int i = 0, n = getSiteCount(); i < n; i++) {
                if (isHeterozygous(taxon, i)) {
                    result++;
                }
            }
            return result;
        } else {
            return myBaseAlignment.getHeterozygousCountForTaxon(translateTaxon(taxon));
        }
    }

    @Override
    public byte[] getDepthForAlleles(int taxon, int site) {
        int taxaIndex = translateTaxon(taxon);
        if (taxaIndex == -1) {
            return null;
        } else {
            return myBaseAlignment.getDepthForAlleles(taxaIndex, translateSite(site));
        }
    }

    @Override
    public byte[] getAllelesByScope(ALLELE_SCOPE_TYPE scope, int site) {
        if (scope == ALLELE_SCOPE_TYPE.Frequency) {
            return getAlleles(site);
        } else {
            return myBaseAlignment.getAllelesByScope(scope, translateSite(site));
        }
    }

    @Override
    public BitSet getAllelePresenceForAllTaxaByScope(ALLELE_SCOPE_TYPE scope, int site, int alleleNumber) {
        if (scope == ALLELE_SCOPE_TYPE.Frequency) {
            return getAllelePresenceForAllTaxa(site, alleleNumber);
        } else {
            int numTaxa = getSequenceCount();
            BitSet result = new OpenBitSet(numTaxa);
            BitSet baseBitSet = myBaseAlignment.getAllelePresenceForAllTaxaByScope(scope, translateSite(site), alleleNumber);
            for (int i = 0; i < numTaxa; i++) {
                if (baseBitSet.fastGet(translateTaxon(i))) {
                    result.fastSet(i);
                }
            }
            return UnmodifiableBitSet.getInstance(result);
        }

    }

    @Override
    public BitSet getPhasedAllelePresenceForAllSites(int taxon, boolean firstParent, int alleleNumber) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public BitSet getPhasedAllelePresenceForAllTaxa(int site, boolean firstParent, int alleleNumber) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public long[] getPhasedAllelePresenceForSitesBlock(int taxon, boolean firstParent, int alleleNumber, int startBlock, int endBlock) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String getBaseAsStringRange(int taxon, int startSite, int endSite) {
        StringBuilder builder = new StringBuilder();
        for (int i = startSite; i < endSite; i++) {
            if (i != startSite) {
                builder.append(";");
            }
            builder.append(getBaseAsString(taxon, i));
        }
        return builder.toString();
    }

    @Override
    public String getBaseAsStringRow(int taxon) {
        return getBaseAsStringRange(taxon, 0, getSiteCount());
    }

    @Override
    public byte[] getReference(int startSite, int endSite) {
        if (!hasReference()) {
            return null;
        }

        byte[] result = new byte[endSite - startSite];
        for (int i = startSite; i < endSite; i++) {
            result[i] = getReferenceAllele(i);
        }
        return result;
    }

    @Override
    public int getChromosomeSiteCount(Chromosome chromosome) {
        int[] startEnd = getStartAndEndOfChromosome(chromosome);
        return startEnd[1] - startEnd[0];
    }

    @Override
    public boolean isAllPolymorphic() {
        for (int i = 0, n = getSiteCount(); i < n; i++) {
            if (!isPolymorphic(i)) {
                return false;
            }
        }

        return true;
    }

    @Override
    public String getMajorAlleleAsString(int site) {
        return getBaseAsString(site, getMajorAllele(site));
    }

    @Override
    public String getMinorAlleleAsString(int site) {
        return getBaseAsString(site, getMinorAllele(site));
    }

    @Override
    public byte[] getMinorAlleles(int site) {
        int[][] alleles = getAllelesSortedByFrequency(site);
        int resultSize = alleles[0].length - 1;
        byte[] result = new byte[resultSize];
        for (int i = 0; i < resultSize; i++) {
            result[i] = (byte) alleles[0][i + 1];
        }
        return result;
    }

    @Override
    public String getTaxaName(int index) {
        return myTaxaList.getTaxaName(index);
    }

    @Override
    public String getFullTaxaName(int index) {
        return myTaxaList.getFullTaxaName(index);
    }

    @Override
    public Alignment[] getAlignments() {
        return new Alignment[]{this};
    }

    @Override
    public String getDiploidAsString(int site, byte value) {
        String[] alleleStates = getAlleleEncodings(site);
        return alleleStates[(value >>> 4) & 0xf] + ":" + alleleStates[value & 0xf];
    }

    @Override
    public Object[][] getDiploidCounts() {

        int numSites = getSiteCount();
        int numTaxa = getSequenceCount();

        Map<String, Long> diploidValueCounts = new HashMap<String, Long>();
        for (int c = 0; c < numSites; c++) {
            Object[][] diploids = getDiploidsSortedByFrequency(c);
            for (int i = 0; i < diploids[0].length; i++) {
                String current = (String) diploids[0][i];
                Long count = (long) ((Integer) diploids[1][i]).intValue();
                Long num = diploidValueCounts.get(current);
                if (num == null) {
                    diploidValueCounts.put(current, count);
                } else {
                    diploidValueCounts.put(current, (num + count));
                }
            }
        }

        Object[][] result = new Object[2][diploidValueCounts.size()];

        int i = 0;
        Iterator itr = diploidValueCounts.keySet().iterator();
        while (itr.hasNext()) {
            String key = (String) itr.next();
            Long count = diploidValueCounts.get(key);
            result[0][i] = key;
            result[1][i++] = count;
        }

        boolean change = true;
        while (change) {

            change = false;

            for (int k = 0, n = diploidValueCounts.size() - 1; k < n; k++) {

                if ((Long) result[1][k] < (Long) result[1][k + 1]) {

                    Object temp = result[0][k];
                    result[0][k] = result[0][k + 1];
                    result[0][k + 1] = temp;

                    Object tempCount = result[1][k];
                    result[1][k] = result[1][k + 1];
                    result[1][k + 1] = tempCount;

                    change = true;
                }
            }

        }

        return result;
    }

    @Override
    public Object[][] getMajorMinorCounts() {

        String[][] alleleStates = getAlleleEncodings();

        if (alleleStates.length != 1) {
            return new Object[0][0];
        }

        int numSites = getSiteCount();
        long[][] counts = new long[16][16];

        if (getMaxNumAlleles() >= 2) {
            for (int site = 0; site < numSites; site++) {
                byte[] alleles = getAlleles(site);
                byte indexI = alleles[0];
                byte indexJ = alleles[1];
                if (indexJ == UNKNOWN_ALLELE) {
                    indexJ = indexI;
                }
                counts[indexI][indexJ]++;
            }
        } else {
            for (int site = 0; site < numSites; site++) {
                byte[] alleles = getAlleles(site);
                byte indexI = alleles[0];
                counts[indexI][indexI]++;
            }
        }

        int numAlleles = 0;
        for (byte x = 0; x < 16; x++) {
            for (byte y = 0; y < 16; y++) {
                if (counts[x][y] != 0) {
                    numAlleles++;
                }
            }
        }

        Object[][] result = new Object[2][numAlleles];
        int nextResult = 0;
        for (byte x = 0; x < 16; x++) {
            for (byte y = 0; y < 16; y++) {
                if (counts[x][y] != 0) {
                    result[0][nextResult] = getBaseAsString(0, x) + ":" + getBaseAsString(0, y);
                    result[1][nextResult++] = counts[x][y];
                }
            }
        }

        boolean change = true;
        while (change) {

            change = false;

            for (int k = 0; k < numAlleles - 1; k++) {

                if ((Long) result[1][k] < (Long) result[1][k + 1]) {

                    Object temp = result[0][k];
                    result[0][k] = result[0][k + 1];
                    result[0][k + 1] = temp;

                    Object tempCount = result[1][k];
                    result[1][k] = result[1][k + 1];
                    result[1][k + 1] = tempCount;

                    change = true;
                }
            }

        }

        return result;
    }

    @Override
    public BitStorage getBitStorage(ALLELE_SCOPE_TYPE scopeType) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public PositionList getPositionList() {
        PositionListBuilder pLB=new PositionListBuilder();
        PositionList basePL=getBaseAlignment().getPositionList();
        for (int i=0; i<getSiteCount(); i++) {
            pLB.add(basePL.get(translateSite(i)));
        }
        return pLB.build();
    }

    @Override
    public Genotype getGenotypeMatrix() {
        throw new UnsupportedOperationException("Not supported yet.");
    }
}
