/*
 * CoreAlignment
 */
package net.maizegenetics.dna.snp;

import java.util.HashMap;
import java.util.Map;
import static net.maizegenetics.dna.snp.Alignment.ALLELE_SORT_TYPE.Frequency;
import static net.maizegenetics.dna.snp.Alignment.ALLELE_SORT_TYPE.Reference;
import net.maizegenetics.dna.snp.score.SiteScore;
import net.maizegenetics.dna.snp.genotype.Genotype;
import net.maizegenetics.dna.snp.bit.BitStorage;
import net.maizegenetics.dna.snp.bit.DynamicBitStorage;
import net.maizegenetics.dna.snp.depth.AlleleDepth;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.util.BitSet;
import org.apache.log4j.Logger;

/**
 *
 * @author terry
 */
public class CoreAlignment implements Alignment {

    private static final Logger myLogger = Logger.getLogger(CoreAlignment.class);
    private final Genotype myGenotype;
    private final Map<ALLELE_SORT_TYPE, BitStorage> myBitStorage = new HashMap<ALLELE_SORT_TYPE, BitStorage>();
    private final PositionList myPositionList;
    private final TaxaList myTaxaList;
    private final SiteScore mySiteScore;
    private final AlleleDepth myAlleleDepth;
    private final int mySiteCount;
    private final int myTaxaCount;

    CoreAlignment(Genotype genotype, PositionList positionList, TaxaList taxaList, SiteScore siteScore, AlleleDepth alleleDepth) {
        //todo need check dimensions
        myGenotype = genotype;
        myPositionList = positionList;
        myTaxaList = taxaList;
        mySiteScore = siteScore;
        myAlleleDepth = alleleDepth;
        mySiteCount = myPositionList.getSiteCount();
        myTaxaCount = myTaxaList.getTaxaCount();
    }

    CoreAlignment(Genotype genotype, PositionList positionList, TaxaList taxaList) {
        this(genotype, positionList, taxaList, null, null);
    }

    @Override
    public Genotype getGenotypeMatrix() {
        return myGenotype;
    }

    @Override
    public byte genotype(int taxon, int site) {
        return myGenotype.genotype(taxon, site);
    }

    @Override
    public byte[] genotypeArray(int taxon, int site) {
        return myGenotype.genotypeArray(taxon, site);
    }

    @Override
    public byte genotype(int taxon, Chromosome chromosome, int physicalPosition) {
        return myGenotype.genotype(taxon, myPositionList.siteOfPhysicalPosition(physicalPosition, chromosome));
    }

    @Override
    public byte[] genotypeRange(int taxon, int startSite, int endSite) {
        return myGenotype.genotypeRange(taxon, startSite, endSite);
    }

    @Override
    public byte[] genotypeAllSites(int taxon) {
        return myGenotype.genotypeAllSites(taxon);
    }

    @Override
    public BitSet allelePresenceForAllSites(int taxon, int alleleNumber) {
        return getBitStorage(ALLELE_SORT_TYPE.Frequency).allelePresenceForAllSites(taxon, alleleNumber);
    }

    @Override
    public BitSet allelePresenceForAllTaxa(int site, int alleleNumber) {
        return getBitStorage(ALLELE_SORT_TYPE.Frequency).allelePresenceForAllTaxa(site, alleleNumber);
    }

    @Override
    public long[] allelePresenceForSitesBlock(int taxon, int alleleNumber, int startBlock, int endBlock) {
        return getBitStorage(ALLELE_SORT_TYPE.Frequency).allelePresenceForSitesBlock(taxon, alleleNumber, startBlock, endBlock);
    }

    @Override
    public BitSet haplotypeAllelePresenceForAllSites(int taxon, boolean firstParent, int alleleNumber) {
        return getBitStorage(ALLELE_SORT_TYPE.Frequency).haplotypeAllelePresenceForAllSites(taxon, firstParent, alleleNumber);
    }

    @Override
    public BitSet haplotypeAllelePresenceForAllTaxa(int site, boolean firstParent, int alleleNumber) {
        return getBitStorage(ALLELE_SORT_TYPE.Frequency).haplotypeAllelePresenceForAllTaxa(site, firstParent, alleleNumber);
    }

    @Override
    public long[] haplotypeAllelePresenceForSitesBlock(int taxon, boolean firstParent, int alleleNumber, int startBlock, int endBlock) {
        return getBitStorage(ALLELE_SORT_TYPE.Frequency).haplotypeAllelePresenceForSitesBlock(taxon, firstParent, alleleNumber, startBlock, endBlock);
    }

    @Override
    public String genotypeAsString(int taxon, int site) {
        return myGenotype.genotypeAsString(taxon, site);
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
    public String[] genotypeAsStringArray(int taxon, int site) {
        return myGenotype.genotypeAsStringArray(taxon, site);
    }

    @Override
    public byte referenceGenotype(int site) {
        return myPositionList.referenceGenotype(site);
    }

    @Override
    public byte[] referenceGenotypes(int startSite, int endSite) {
        return myPositionList.referenceGenotypes(startSite, endSite);
    }

    @Override
    public byte[] referenceGenotypeForAllSites() {
        return myPositionList.referenceGenotypeForAllSites();
    }

    @Override
    public boolean hasReference() {
        return myPositionList.hasReference();
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
    public PositionList getPositionList() {
        return myPositionList;
    }

    @Override
    public String siteName(int site) {
        return myPositionList.siteName(site);
    }

    @Override
    public int numberOfSites() {
        return mySiteCount;
    }

    @Override
    public int getChromosomeSiteCount(Chromosome chromosome) {
        return myPositionList.getChromosomeSiteCount(chromosome);
    }

    @Override
    public int[] getStartAndEndOfChromosome(Chromosome chromosome) {
        return myPositionList.getStartAndEndOfChromosome(chromosome);
    }

    @Override
    public int numberOfTaxa() {
        return myTaxaCount;
    }

    @Override
    public int chromosomalPosition(int site) {
        return myPositionList.chromosomalPosition(site);
    }

    @Override
    public int siteOfPhysicalPosition(int physicalPosition, Chromosome chromosome) {
        return myPositionList.siteOfPhysicalPosition(physicalPosition, chromosome);
    }

    @Override
    public int siteOfPhysicalPosition(int physicalPosition, Chromosome chromosome, String snpName) {
        return myPositionList.siteOfPhysicalPosition(physicalPosition, chromosome, snpName);
    }

    @Override
    public int[] getPhysicalPositions() {
        return myPositionList.getPhysicalPositions();
    }

    @Override
    public String getChromosomeName(int site) {
        return myPositionList.getChromosomeName(site);
    }

    @Override
    public Chromosome getChromosome(int site) {
        return myPositionList.getChromosome(site);
    }

    @Override
    public Chromosome getChromosome(String name) {
        return myPositionList.getChromosome(name);
    }

    @Override
    public Chromosome[] getChromosomes() {
        return myPositionList.getChromosomes();
    }

    @Override
    public int getNumChromosomes() {
        return myPositionList.getNumChromosomes();
    }

    @Override
    public int[] getChromosomesOffsets() {
        return myPositionList.getChromosomesOffsets();
    }

    @Override
    public float siteScore(int taxon, int site) {
        if (mySiteScore == null) {
            throw new IllegalStateException("CoreAlignment: getSiteScore: This Alignment has no Site Scores.");
        }
        return mySiteScore.siteScore(taxon, site);
    }

    @Override
    public float[][] siteScores() {
        if (mySiteScore == null) {
            throw new IllegalStateException("CoreAlignment: getSiteScores: This Alignment has no Site Scores.");
        }
        return mySiteScore.siteScores();
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
    public SITE_SCORE_TYPE siteScoreType() {
        return mySiteScore.siteScoreType();
    }

    @Override
    public int getIndelSize(int site) {
        return myPositionList.getIndelSize(site);
    }

    @Override
    public boolean isIndel(int site) {
        return myPositionList.isIndel(site);
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
        return myGenotype.getMajorAllele(site);
    }

    @Override
    public String getMajorAlleleAsString(int site) {
        return myGenotype.getMajorAlleleAsString(site);
    }

    @Override
    public byte getMinorAllele(int site) {
        return myGenotype.getMinorAllele(site);
    }

    @Override
    public String getMinorAlleleAsString(int site) {
        return myGenotype.getMinorAlleleAsString(site);
    }

    @Override
    public byte[] getMinorAlleles(int site) {
        return myGenotype.getMinorAlleles(site);
    }

    @Override
    public byte[] getAlleles(int site) {
        return myGenotype.getAlleles(site);
    }

    @Override
    public double getMinorAlleleFrequency(int site) {
        return myGenotype.getMinorAlleleFrequency(site);
    }

    @Override
    public double getMajorAlleleFrequency(int site) {
        return myGenotype.getMajorAlleleFrequency(site);
    }

    @Override
    public TaxaList taxa() {
        return myTaxaList;
    }

    @Override
    public String taxaName(int index) {
        return myTaxaList.getTaxaName(index);
    }

    @Override
    public String genomeVersion() {
        return myPositionList.genomeVersion();
    }

    @Override
    public boolean isPositiveStrand(int site) {
        return myPositionList.isPositiveStrand(site);
    }

    @Override
    public Alignment[] getAlignments() {
        return new Alignment[]{this};
    }

    @Override
    public int[][] allelesSortedByFrequency(int site) {
        return myGenotype.allelesSortedByFrequency(site);
    }

    @Override
    public Object[][] genosSortedByFrequency(int site) {
        return myGenotype.genosSortedByFrequency(site);
    }

    @Override
    public boolean isPhased() {
        return myGenotype.isPhased();
    }

    @Override
    public GeneticMap getGeneticMap() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public boolean retainsRareAlleles() {
        return myGenotype.retainsRareAlleles();
    }

    @Override
    public String[][] alleleDefinitions() {
        return myGenotype.alleleDefinitions();
    }

    @Override
    public String[] alleleDefinitions(int site) {
        return myGenotype.alleleDefinitions(site);
    }

    @Override
    public String genotypeAsString(int site, byte value) {
        return myGenotype.genotypeAsString(site, value);
    }

    @Override
    public String getDiploidAsString(int site, byte value) {
        return myGenotype.getDiploidAsString(site, value);
    }

    @Override
    public int maxNumAlleles() {
        return myGenotype.maxNumAlleles();
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
    public int getMinorAlleleCount(int site) {
        return myGenotype.getMinorAlleleCount(site);
    }

    @Override
    public int getMajorAlleleCount(int site) {
        return myGenotype.getMajorAlleleCount(site);
    }

    @Override
    public Object[][] genoCounts() {
        return myGenotype.genoCounts();
    }

    @Override
    public Object[][] getMajorMinorCounts() {
        return myGenotype.getMajorMinorCounts();
    }

    @Override
    public int totalGametesNonMissingForTaxon(int taxon) {
        return myGenotype.totalGametesNonMissingForTaxon(taxon);
    }

    @Override
    public int getHeterozygousCountForTaxon(int taxon) {
        return myGenotype.getHeterozygousCountForTaxon(taxon);
    }

    @Override
    public int totalNonMissingForTaxon(int taxon) {
        return myGenotype.totalNonMissingForTaxon(taxon);
    }

    @Override
    public byte[] getDepthForAlleles(int taxon, int site) {
        return myAlleleDepth.getDepthForAlleles(taxon, site);
    }

    @Override
    public byte[] getAllelesByScope(ALLELE_SORT_TYPE scope, int site) {
        switch (scope) {
            case Frequency:
                return getAlleles(site);
            default:
                myLogger.warn("getAllelesByScope: Unsupported type: " + scope);
                return null;
        }
    }

    @Override
    public BitSet allelePresenceForAllTaxaBySortType(ALLELE_SORT_TYPE type, int site, int alleleNumber) {
        return getBitStorage(type).allelePresenceForAllTaxa(site, alleleNumber);
    }

    @Override
    public BitStorage getBitStorage(ALLELE_SORT_TYPE scopeType) {

        BitStorage result = myBitStorage.get(scopeType);
        if (result != null) {
            return result;
        }

        switch (scopeType) {
            case Frequency:
                result = new DynamicBitStorage(myGenotype, scopeType, myGenotype.getMajorAlleleForAllSites(), myGenotype.getMinorAlleleForAllSites());
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
}
