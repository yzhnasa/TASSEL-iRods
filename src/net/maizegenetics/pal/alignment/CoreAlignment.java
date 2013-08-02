/*
 * CoreAlignment
 */
package net.maizegenetics.pal.alignment;

import net.maizegenetics.pal.alignment.score.SiteScore;
import net.maizegenetics.pal.alignment.genotype.Genotype;
import net.maizegenetics.pal.alignment.bit.BitStorage;
import net.maizegenetics.pal.alignment.depth.AlleleDepth;
import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.site.SiteList;
import net.maizegenetics.pal.taxa.TaxaList;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.ProgressListener;

/**
 *
 * @author terry
 */
public class CoreAlignment implements Alignment {

    private final Genotype myGenotype;
    private final BitStorage myBitStorage;
    private final SiteList mySiteList;
    private final TaxaList myTaxaList;
    private final SiteScore mySiteScore;
    private final AlleleDepth myAlleleDepth;

    private CoreAlignment(Genotype genotype, BitStorage bitStorage, SiteList siteList, TaxaList taxaList, SiteScore siteScore, AlleleDepth alleleDepth) {
        myGenotype = genotype;
        myBitStorage = bitStorage;
        mySiteList = siteList;
        myTaxaList = taxaList;
        mySiteScore = siteScore;
        myAlleleDepth = alleleDepth;
    }

    @Override
    public byte getBase(int taxon, int site) {
        return myGenotype.getBase(taxon, site);
    }

    @Override
    public byte[] getBaseArray(int taxon, int site) {
        return myGenotype.getBaseArray(taxon, site);
    }

    @Override
    public byte getBase(int taxon, Locus locus, int physicalPosition) {
        return myGenotype.getBase(taxon, locus, physicalPosition);
    }

    @Override
    public byte[] getBaseRange(int taxon, int startSite, int endSite) {
        return myGenotype.getBaseRange(taxon, startSite, endSite);
    }

    @Override
    public byte[] getBaseRow(int taxon) {
        return myGenotype.getBaseRow(taxon);
    }

    @Override
    public BitSet getAllelePresenceForAllSites(int taxon, int alleleNumber) {
        return myBitStorage.getAllelePresenceForAllSites(taxon, alleleNumber);
    }

    @Override
    public BitSet getAllelePresenceForAllTaxa(int site, int alleleNumber) {
        return myBitStorage.getAllelePresenceForAllTaxa(site, alleleNumber);
    }

    @Override
    public long[] getAllelePresenceForSitesBlock(int taxon, int alleleNumber, int startBlock, int endBlock) {
        return myBitStorage.getAllelePresenceForSitesBlock(taxon, alleleNumber, startBlock, endBlock);
    }

    @Override
    public BitSet getPhasedAllelePresenceForAllSites(int taxon, boolean firstParent, int alleleNumber) {
        return myBitStorage.getPhasedAllelePresenceForAllSites(taxon, firstParent, alleleNumber);
    }

    @Override
    public BitSet getPhasedAllelePresenceForAllTaxa(int site, boolean firstParent, int alleleNumber) {
        return myBitStorage.getPhasedAllelePresenceForAllTaxa(site, firstParent, alleleNumber);
    }

    @Override
    public long[] getPhasedAllelePresenceForSitesBlock(int taxon, boolean firstParent, int alleleNumber, int startBlock, int endBlock) {
        return myBitStorage.getPhasedAllelePresenceForSitesBlock(taxon, firstParent, alleleNumber, startBlock, endBlock);
    }

    @Override
    public String getBaseAsString(int taxon, int site) {
        return myGenotype.getBaseAsString(taxon, site);
    }

    @Override
    public String getBaseAsStringRange(int taxon, int startSite, int endSite) {
        return myGenotype.getBaseAsStringRange(taxon, startSite, endSite);
    }

    @Override
    public String getBaseAsStringRow(int taxon) {
        return myGenotype.getBaseAsStringRow(taxon);
    }

    @Override
    public String[] getBaseAsStringArray(int taxon, int site) {
        return myGenotype.getBaseAsStringArray(taxon, site);
    }

    @Override
    public byte getReferenceAllele(int site) {
        return mySiteList.getReferenceAllele(site);
    }

    @Override
    public byte[] getReference(int startSite, int endSite) {
        return mySiteList.getReference(startSite, endSite);
    }

    @Override
    public byte[] getReference() {
        return mySiteList.getReference();
    }

    @Override
    public boolean hasReference() {
        return mySiteList.hasReference();
    }

    @Override
    public boolean isHeterozygous(int taxon, int site) {
        return myGenotype.isHeterozygous(taxon, site);
    }

    @Override
    public int getHeterozygousCount(int site) {
        return myGenotype.getHeterozygousCount(site);
    }

    @Override
    public String[] getSNPIDs() {
        return mySiteList.getSNPIDs();
    }

    @Override
    public String getSNPID(int site) {
        return mySiteList.getSNPID(site);
    }

    @Override
    public int getSiteCount() {
        return mySiteList.getSiteCount();
    }

    @Override
    public int getLocusSiteCount(Locus locus) {
        return mySiteList.getLocusSiteCount(locus);
    }

    @Override
    public int[] getStartAndEndOfLocus(Locus locus) {
        return mySiteList.getStartAndEndOfLocus(locus);
    }

    @Override
    public int getSequenceCount() {
        return myTaxaList.getSequenceCount();
    }

    @Override
    public int getTaxaCount() {
        return myTaxaList.getTaxaCount();
    }

    @Override
    public int getPositionInLocus(int site) {
        return mySiteList.getPositionInLocus(site);
    }

    @Override
    public int getSiteOfPhysicalPosition(int physicalPosition, Locus locus) {
        return mySiteList.getSiteOfPhysicalPosition(physicalPosition, locus);
    }

    @Override
    public int getSiteOfPhysicalPosition(int physicalPosition, Locus locus, String snpID) {
        return mySiteList.getSiteOfPhysicalPosition(physicalPosition, locus, snpID);
    }

    @Override
    public int[] getPhysicalPositions() {
        return mySiteList.getPhysicalPositions();
    }

    @Override
    public byte getPositionType(int site) {
        return mySiteList.getPositionType(site);
    }

    @Override
    public byte[] getPositionTypes() {
        return mySiteList.getPositionTypes();
    }

    @Override
    public String getLocusName(int site) {
        return mySiteList.getLocusName(site);
    }

    @Override
    public Locus getLocus(int site) {
        return mySiteList.getLocus(site);
    }

    @Override
    public Locus getLocus(String name) {
        return mySiteList.getLocus(name);
    }

    @Override
    public Locus[] getLoci() {
        return mySiteList.getLoci();
    }

    @Override
    public int getNumLoci() {
        return mySiteList.getNumLoci();
    }

    @Override
    public int[] getLociOffsets() {
        return mySiteList.getLociOffsets();
    }

    @Override
    public float getSiteScore(int seq, int site) {
        if (mySiteScore == null) {
            throw new IllegalStateException("CoreAlignment: getSiteScore: This Alignment has no Site Scores.");
        }
        return mySiteScore.getSiteScore(seq, site);
    }

    @Override
    public float[][] getSiteScores() {
        if (mySiteScore == null) {
            throw new IllegalStateException("CoreAlignment: getSiteScores: This Alignment has no Site Scores.");
        }
        return mySiteScore.getSiteScores();
    }

    @Override
    public boolean hasSiteScores() {
        if (mySiteScore == null) {
            return false;
        } else {
            return true;
        }
    }

    @Override
    public SITE_SCORE_TYPE getSiteScoreType() {
        return mySiteScore.getSiteScoreType();
    }

    @Override
    public int getIndelSize(int site) {
        return mySiteList.getIndelSize(site);
    }

    @Override
    public boolean isIndel(int site) {
        return mySiteList.isIndel(site);
    }

    @Override
    public boolean isAllPolymorphic() {
        return myGenotype.isAllPolymorphic();
    }

    @Override
    public boolean isPolymorphic(int site) {
        return myGenotype.isPolymorphic(site);
    }

    @Override
    public byte getMajorAllele(int site) {
        return mySiteList.getMajorAllele(site);
    }

    @Override
    public String getMajorAlleleAsString(int site) {
        return mySiteList.getMajorAlleleAsString(site);
    }

    @Override
    public byte getMinorAllele(int site) {
        return mySiteList.getMinorAllele(site);
    }

    @Override
    public String getMinorAlleleAsString(int site) {
        return mySiteList.getMinorAlleleAsString(site);
    }

    @Override
    public byte[] getMinorAlleles(int site) {
        return mySiteList.getMinorAlleles(site);
    }

    @Override
    public byte[] getAlleles(int site) {
        return mySiteList.getAlleles(site);
    }

    @Override
    public double getMinorAlleleFrequency(int site) {
        return mySiteList.getMinorAlleleFrequency(site);
    }

    @Override
    public double getMajorAlleleFrequency(int site) {
        return mySiteList.getMajorAlleleFrequency(site);
    }

    @Override
    public IdGroup getIdGroup() {
        // This Method signature needs to be changed.
        // return myTaxaList;
        return null;
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
    public String getGenomeAssembly() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public boolean isPositiveStrand(int site) {
        return mySiteList.isPositiveStrand(site);
    }

    @Override
    public Alignment[] getAlignments() {
        return myGenotype.getAlignments();
    }

    @Override
    public int[][] getAllelesSortedByFrequency(int site) {
        return mySiteList.getAllelesSortedByFrequency(site);
    }

    @Override
    public Object[][] getDiploidssSortedByFrequency(int site) {
        return mySiteList.getDiploidssSortedByFrequency(site);
    }

    @Override
    public boolean isPhased() {
        return myGenotype.isPhased();
    }

    @Override
    public GeneticMap getGeneticMap() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public boolean retainsRareAlleles() {
        return myGenotype.retainsRareAlleles();
    }

    @Override
    public String[][] getAlleleEncodings() {
        return myGenotype.getAlleleEncodings();
    }

    @Override
    public String[] getAlleleEncodings(int site) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public String getBaseAsString(int site, byte value) {
        return myGenotype.getBaseAsString(site, value);
    }

    @Override
    public String getDiploidAsString(int site, byte value) {
        return myGenotype.getDiploidAsString(site, value);
    }

    @Override
    public int getMaxNumAlleles() {
        return myGenotype.getMaxNumAlleles();
    }

    @Override
    public int getTotalNumAlleles() {
        return myGenotype.getTotalNumAlleles();
    }

    @Override
    public int getTotalGametesNotMissing(int site) {
        return myGenotype.getTotalGametesNotMissing(site);
    }

    @Override
    public int getTotalNotMissing(int site) {
        return myGenotype.getTotalNotMissing(site);
    }

    @Override
    public int getMinorAlleleCount(int site) {
        return myGenotype.getMinorAlleleCount(site);
    }

    @Override
    public int getMajorAlleleCount(int site) {
        return myGenotype.getMajorAlleleCount(site);
    }

    @Override
    public Object[][] getDiploidCounts() {
        return myGenotype.getDiploidCounts();
    }

    @Override
    public Object[][] getMajorMinorCounts() {
        return myGenotype.getMajorMinorCounts();
    }

    @Override
    public int getTotalGametesNotMissingForTaxon(int taxon) {
        return myGenotype.getTotalGametesNotMissingForTaxon(taxon);
    }

    @Override
    public int getHeterozygousCountForTaxon(int taxon) {
        return myGenotype.getHeterozygousCountForTaxon(taxon);
    }

    @Override
    public int getTotalNotMissingForTaxon(int taxon) {
        return myGenotype.getTotalNotMissingForTaxon(taxon);
    }

    @Override
    public boolean isSBitFriendly() {
        return myGenotype.isSBitFriendly();
    }

    @Override
    public boolean isTBitFriendly() {
        return myGenotype.isTBitFriendly();
    }

    @Override
    public void optimizeForTaxa(ProgressListener listener) {
        myGenotype.optimizeForTaxa(listener);
    }

    @Override
    public void optimizeForSites(ProgressListener listener) {
        optimizeForSites(listener);
    }

    @Override
    public byte[] getDepthForAlleles(int taxon, int site) {
        return myAlleleDepth.getDepthForAlleles(taxon, site);
    }

    @Override
    public byte[] getAllelesByScope(ALLELE_SCOPE_TYPE scope, int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public BitSet getAllelePresenceForAllTaxaByScope(ALLELE_SCOPE_TYPE scope, int site, int alleleNumber) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
}
