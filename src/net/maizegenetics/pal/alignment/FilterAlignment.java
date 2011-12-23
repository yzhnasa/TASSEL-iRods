/*
 * FilterAlignment
 */
package net.maizegenetics.pal.alignment;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.Identifier;
import net.maizegenetics.pal.ids.SimpleIdGroup;
import net.maizegenetics.util.BitSet;

/**
 * All taxa and site filtering should be controlled through this class.
 * It essentially creates views of the baseAlignment
 * @author terry
 */
public class FilterAlignment extends AbstractAlignment {

    private static final long serialVersionUID = -5197800047652332969L;
    private final boolean myIsTaxaFilter;
    private final boolean myIsSiteFilter;
    private final boolean myIsSiteFilterByRange;
    private final Alignment myBaseAlignment;
    private final int[] myTaxaRedirect;
    private final int[] mySiteRedirect;
    private final int myRangeStart;
    private final int myRangeEnd;
    private final Locus[] myLoci;

    private FilterAlignment(Alignment a, IdGroup subIdGroup, int[] taxaRedirect, FilterAlignment original) {

        super(subIdGroup);

        if (subIdGroup.getIdCount() != taxaRedirect.length) {
            throw new IllegalArgumentException("FilterAlignment: init: subIdGroup should be same size as taxaRedirect.");
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
            myLoci = myBaseAlignment.getLoci();
        } else {
            myIsSiteFilter = original.isSiteFilter();
            myIsSiteFilterByRange = original.isSiteFilterByRange();
            mySiteRedirect = original.getSiteRedirect();
            myRangeStart = original.getRangeStart();
            myRangeEnd = original.getRangeEnd();
            myLoci = original.getLoci();
        }

    }

    public static FilterAlignment getInstance(Alignment a, IdGroup subIdGroup) {

        Alignment baseAlignment = a;
        FilterAlignment original = null;
        if (baseAlignment instanceof FilterAlignment) {
            original = (FilterAlignment) a;
            baseAlignment = ((FilterAlignment) a).getBaseAlignment();
        }

        int[] taxaRedirect = new int[subIdGroup.getIdCount()];
        for (int i = 0; i < subIdGroup.getIdCount(); i++) {
            int ion = a.getIdGroup().whichIdNumber(subIdGroup.getIdentifier(i));
            if (ion < 0) {
                taxaRedirect[i] = -1;
            } else {
                int idn = baseAlignment.getIdGroup().whichIdNumber(subIdGroup.getIdentifier(i));
                if (idn > -1) {
                    taxaRedirect[i] = idn;
                } else {
                    taxaRedirect[i] = -1;
                }
            }
        }

        return new FilterAlignment(baseAlignment, subIdGroup, taxaRedirect, original);

    }

    /**
     * Removes specified IDs.
     *
     * @param a alignment to filter
     * @param subIdGroup specified IDs
     *
     * @return Filtered Alignment
     */
    public static FilterAlignment getInstanceRemoveIDs(Alignment a, IdGroup subIdGroup) {

        List result = new ArrayList();
        IdGroup current = a.getIdGroup();
        for (int i = 0, n = current.getIdCount(); i < n; i++) {
            if (subIdGroup.whichIdNumber(current.getIdentifier(i)) == -1) {
                result.add(current.getIdentifier(i));
            }
        }
        Identifier[] ids = new Identifier[result.size()];
        result.toArray(ids);
        return FilterAlignment.getInstance(a, new SimpleIdGroup(ids));

    }

    /**
     * Constructor
     *
     * @param a base alignment
     * @param startSite start site (included)
     * @param endSite end site (included)
     */
    private FilterAlignment(Alignment a, int startSite, int endSite, FilterAlignment original) {

        super(original == null ? a.getIdGroup() : original.getIdGroup());

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
        myLoci = getLociFromBase();

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

        super(original == null ? a.getIdGroup() : original.getIdGroup());

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
        myLoci = getLociFromBase();

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

    public byte getBase(int taxon, int site) {
        int taxaIndex = translateTaxon(taxon);
        if (taxaIndex == -1) {
            return Alignment.UNKNOWN_ALLELE;
        } else {
            return myBaseAlignment.getBase(taxaIndex, translateSite(site));
        }
    }

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

    public byte getBase(int taxon, Locus locus, int physicalPosition) {
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
     * Returns site of this FilterAlignment based on give site from
     * embedded Alignment.
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

    public int translateTaxon(int taxon) {

        if (myIsTaxaFilter) {
            return myTaxaRedirect[taxon];
        } else {
            return taxon;
        }

    }

    private Locus[] getLociFromBase() {

        if ((!myIsSiteFilter) && (!myIsSiteFilterByRange)) {
            return myBaseAlignment.getLoci();
        }

        int numSites = getSiteCount();
        List loci = new ArrayList();
        for (int i = 0; i < numSites; i++) {
            Locus current = getLocus(i);
            if (!loci.contains(current)) {
                loci.add(current);
            }
        }

        Locus[] result = new Locus[loci.size()];
        loci.toArray(result);
        return result;

    }

    public int getIndelSize(int site) {
        return myBaseAlignment.getIndelSize(translateSite(site));
    }

    public Locus[] getLoci() {
        return myLoci;
    }

    public Locus getLocus(int site) {
        return myBaseAlignment.getLocus(translateSite(site));
    }

    public int getLocusSiteCount(Locus locus) {

        int numSites = getSiteCount();
        int result = 0;
        for (int i = 0; i < numSites; i++) {
            if (getLocus(translateSite(i)) == locus) {
                result++;
            }
        }

        return result;

    }

    public int getNumLoci() {
        return myLoci.length;
    }

    public int getPositionInLocus(int site) {
        return myBaseAlignment.getPositionInLocus(translateSite(site));
    }

    public byte getPositionType(int site) {
        return myBaseAlignment.getPositionType(translateSite(site));
    }

    public byte[] getPositionTypes() {

        int numSites = getSiteCount();
        byte[] result = new byte[numSites];
        for (int i = 0; i < numSites; i++) {
            result[i] = myBaseAlignment.getPositionType(translateSite(i));
        }

        return result;

    }

    public float[][] getSiteScores() {

        if (!myBaseAlignment.hasSiteScores()) {
            return null;
        }

        int numSites = getSiteCount();
        int numSeqs = getSequenceCount();
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

    public byte getReferenceAllele(int site) {
        return myBaseAlignment.getReferenceAllele(translateSite(site));
    }

    public int getSiteCount() {

        if (myIsSiteFilterByRange) {
            return myRangeEnd - myRangeStart + 1;
        } else if (myIsSiteFilter) {
            return mySiteRedirect.length;
        } else {
            return myBaseAlignment.getSiteCount();
        }

    }

    public String getSNPID(int site) {
        return myBaseAlignment.getSNPID(translateSite(site));
    }

    public String[] getSNPIDs() {

        int numSites = getSiteCount();
        String[] result = new String[numSites];
        for (int i = 0; i < numSites; i++) {
            result[i] = myBaseAlignment.getSNPID(translateSite(i));
        }

        return result;

    }

    public boolean hasReference() {
        return myBaseAlignment.hasReference();
    }

    public boolean isIndel(int site) {
        return myBaseAlignment.isIndel(translateSite(site));
    }

    public int getSiteOfPhysicalPosition(int physicalPosition, Locus locus) {
        int temp = myBaseAlignment.getSiteOfPhysicalPosition(physicalPosition, locus);
        if (temp < 0) {
            temp = -temp;
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

    public byte getMajorAllele(int site) {
        return myBaseAlignment.getMajorAllele(translateSite(site));
    }

    public byte getMinorAllele(int site) {
        return myBaseAlignment.getMinorAllele(translateSite(site));
    }

    public byte[] getMinorAlleles(int site) {
        return myBaseAlignment.getMinorAlleles(translateSite(site));
    }

    public byte[] getAlleles(int site) {
        return myBaseAlignment.getAlleles(translateSite(site));
    }

    public double getMinorAlleleFrequency(int site) {
        return myBaseAlignment.getMinorAlleleFrequency(translateSite(site));
    }

    public int[][] getAllelesSortedByFrequency(int site) {
        return myBaseAlignment.getAllelesSortedByFrequency(translateSite(site));
    }

    public byte[] getBaseArray(int taxon, int site) {
        int taxaIndex = translateTaxon(taxon);
        if (taxaIndex == -1) {
            return new byte[]{Alignment.UNKNOWN_ALLELE, Alignment.UNKNOWN_ALLELE};
        } else {
            return myBaseAlignment.getBaseArray(taxaIndex, translateSite(site));
        }
    }

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
                for (int i = 0; i < siteCount; i++) {
                    result[i] = myBaseAlignment.getBase(taxaIndex, i);
                }
                return result;
            }
        }
    }

    public BitSet getAllelePresensceForAllSites(int taxon, int alleleNumber) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public BitSet getAllelePresenceForAllTaxa(int site, int alleleNumber) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public long[] getAllelePresensceForSitesBlock(int taxon, int alleleNumber, int startBlock, int endBlock) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public String getBaseAsString(int taxon, int site) {
        int taxaIndex = translateTaxon(taxon);
        if (taxaIndex == -1) {
            return Alignment.UNKNOWN_ALLELE_STR;
        } else {
            return myBaseAlignment.getBaseAsString(taxaIndex, translateSite(site));
        }
    }

    public String[] getBaseAsStringArray(int taxon, int site) {
        int taxaIndex = translateTaxon(taxon);
        if (taxaIndex == -1) {
            return new String[]{Alignment.UNKNOWN_ALLELE_STR, Alignment.UNKNOWN_ALLELE_STR};
        } else {
            return myBaseAlignment.getBaseAsStringArray(taxaIndex, translateSite(site));
        }
    }

    public byte[] getReference() {
        if ((myIsSiteFilterByRange) || (myIsSiteFilter)) {
            return myBaseAlignment.getReference(0, getSiteCount());
        } else {
            return myBaseAlignment.getReference();
        }
    }

    public boolean isHeterozygous(int taxon, int site) {
        int taxaIndex = translateTaxon(taxon);
        if (taxaIndex == -1) {
            return false;
        } else {
            return myBaseAlignment.isHeterozygous(taxaIndex, translateSite(site));
        }
    }

    public int[] getPhysicalPositions() {
        if ((myIsSiteFilterByRange) || (myIsSiteFilter)) {
            int numSites = getSiteCount();
            int[] result = new int[numSites];
            for (int i = 0; i < numSites; i++) {
                result[i] = getPositionInLocus(i);
            }
            return result;
        } else {
            return myBaseAlignment.getPhysicalPositions();
        }
    }

    public String getLocusName(int site) {
        return myBaseAlignment.getLocusName(translateSite(site));
    }

    public int[] getLociOffsets() {

        int[] orgOffsets = myBaseAlignment.getLociOffsets();

        if ((!myIsSiteFilterByRange) && (!myIsSiteFilter)) {
            return orgOffsets;
        }

        List newOffsets = new ArrayList();
        for (int i = 0; i < orgOffsets.length; i++) {
            int current = reverseTranslateSite(i);
            if (current >= 0) {
                newOffsets.add(i);
            }
        }
        int[] result = new int[newOffsets.size()];
        for (int i = 0; i < newOffsets.size(); i++) {
            result[i] = (Integer) newOffsets.get(i);
        }
        return result;

    }

    public float getSiteScore(int seq, int site) {
        int taxaIndex = translateTaxon(seq);
        if (taxaIndex == -1) {
            return Float.NaN;
        } else {
            return myBaseAlignment.getSiteScore(taxaIndex, translateSite(site));
        }
    }

    public boolean hasSiteScores() {
        return myBaseAlignment.hasSiteScores();
    }

    public SITE_SCORE_TYPE getSiteScoreType() {
        return myBaseAlignment.getSiteScoreType();
    }

    public boolean isPolymorphic(int site) {
        return myBaseAlignment.isPolymorphic(translateSite(site));
    }

    public double getMajorAlleleFrequency(int site) {
        return myBaseAlignment.getMajorAlleleFrequency(translateSite(site));
    }

    public String getGenomeAssembly() {
        return myBaseAlignment.getGenomeAssembly();
    }

    public boolean isPositiveStrand(int site) {
        return myBaseAlignment.isPositiveStrand(translateSite(site));
    }

    public boolean isPhased() {
        return myBaseAlignment.isPhased();
    }

    public GeneticMap getGeneticMap() {
        return myBaseAlignment.getGeneticMap();
    }

    public boolean retainsRareAlleles() {
        return myBaseAlignment.retainsRareAlleles();
    }

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

    public String[] getAlleleEncodings(int site) {
        return myBaseAlignment.getAlleleEncodings(translateSite(site));
    }

    public String getBaseAsString(int site, byte value) {
        return myBaseAlignment.getBaseAsString(translateSite(site), value);
    }

    public int getMaxNumAlleles() {
        return myBaseAlignment.getMaxNumAlleles();
    }
}
