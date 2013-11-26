/*
 * CombineAlignment
 */
package net.maizegenetics.dna.snp;

import net.maizegenetics.dna.snp.bit.BitStorage;
import net.maizegenetics.dna.snp.genotype.Genotype;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListUtils;
import net.maizegenetics.util.BitSet;

import java.util.*;

/**
 *
 * @author terry
 */
public class CombineAlignment implements Alignment {

    private static final long serialVersionUID = -5197800047652332969L;
    private final Alignment[] myAlignments;
    private final int[] mySiteOffsets;
    private final Map<Chromosome, Alignment> myChromosomes = new HashMap<>();
    private Chromosome[] myChromosomesList;
    private int[] myChromosomesOffsets;
    private final TaxaList myTaxaList;
    private String[][] myAlleleStates;

    private CombineAlignment(TaxaList taxaList, Alignment[] alignments) {

        myTaxaList = taxaList;
        myAlignments = alignments;
        mySiteOffsets = new int[alignments.length + 1];

        mySiteOffsets[0] = 0;
        int count = 0;
        for (int i = 0; i < alignments.length; i++) {
            count = alignments[i].numberOfSites() + count;
            mySiteOffsets[i + 1] = count;

            Chromosome[] chromosomes = alignments[i].getChromosomes();
            for (int j = 0; j < chromosomes.length; j++) {
                myChromosomes.put(chromosomes[j], alignments[i]);
            }
        }

        initChromosomes();
    }

    /**
     * This factory method combines given alignments. If only one alignment,
     * then it is returned unchanged. Otherwise, this requires that each
     * alignment has the same Identifiers in the same order.
     *
     * @param alignments
     * @return
     */
    public static Alignment getInstance(Alignment[] alignments) {

        if ((alignments == null) || (alignments.length == 0)) {
            throw new IllegalArgumentException("CombineAlignment: getInstance: must provide alignments.");
        }

        if (alignments.length == 1) {
            return alignments[0];
        }

        TaxaList firstGroup = alignments[0].taxa();
        for (int i = 1; i < alignments.length; i++) {
            if (!areTaxaListsEqual(firstGroup, alignments[i].taxa())) {
                throw new IllegalArgumentException("CombineAlignment: getInstance: TaxaLists do not match.");
            }
        }

        return new CombineAlignment(firstGroup, alignments);

    }

    /**
     * This factory method combines given alignments. If only one alignment,
     * then it is returned unchanged. If isUnion equals true, a union join of
     * the Identifiers will be used to construct the combination. Any alignment
     * not containing one of the Identifiers will return unknown value for those
     * locations. If isUnion equals false, a intersect join of the Identifiers
     * will be used.
     *
     * @param alignments alignments to combine
     * @param isUnion whether to union or intersect join
     * @return
     */
    public static Alignment getInstance(Alignment[] alignments, boolean isUnion) {

        if ((alignments == null) || (alignments.length == 0)) {
            throw new IllegalArgumentException("CombineAlignment: getInstance: must provide alignments.");
        }

        if (alignments.length == 1) {
            return alignments[0];
        }

        TaxaList[] groups = new TaxaList[alignments.length];
        for (int i = 0; i < alignments.length; i++) {
            groups[i] = alignments[i].taxa();
        }
        TaxaList newTaxa = null;
        if (isUnion) {
            newTaxa = TaxaListUtils.getAllTaxa(groups);
        } else {
            newTaxa = TaxaListUtils.getCommonTaxa(groups);
        }

        Alignment[] newAlignmentNews = new Alignment[alignments.length];
        for (int i = 0; i < alignments.length; i++) {
            newAlignmentNews[i] = FilterAlignment.getInstance(alignments[i], newTaxa);
        }

        return new CombineAlignment(newTaxa, newAlignmentNews);

    }

    private static boolean areTaxaListsEqual(TaxaList first, TaxaList second) {

        if (first.getTaxaCount() != second.getTaxaCount()) {
            return false;
        }

        for (int i = 0, n = first.getTaxaCount(); i < n; i++) {
            if (!first.get(i).equals(second.get(i))) {
                return false;
            }
        }

        return true;

    }

    private void initChromosomes() {

        List<Integer> offsets = new ArrayList<>();
        List<Chromosome> chromosomes = new ArrayList<>();
        for (int i = 0; i < myAlignments.length; i++) {
            chromosomes.addAll(Arrays.asList(myAlignments[i].getChromosomes()));
            int[] tempOffsets = myAlignments[i].getChromosomesOffsets();
            for (int j = 0; j < tempOffsets.length; j++) {
                offsets.add(tempOffsets[j] + mySiteOffsets[i]);
            }
        }

        myChromosomesList = new Chromosome[chromosomes.size()];
        myChromosomesList = chromosomes.toArray(myChromosomesList);

        myChromosomesOffsets = new int[offsets.size()];
        for (int i = 0; i < offsets.size(); i++) {
            myChromosomesOffsets[i] = (Integer) offsets.get(i);
        }

        if (myChromosomesOffsets.length != myChromosomesList.length) {
            throw new IllegalStateException("CombineAlignment: initChromosomes: number chromosomes offsets should equal number of chromosomes.");
        }

    }

    public byte genotype(int taxon, int site) {
        int translate = translateSite(site);
        return myAlignments[translate].genotype(taxon, site - mySiteOffsets[translate]);
    }

    @Override
    public byte[] genotypeRange(int taxon, int startSite, int endSite) {

        byte[] result = new byte[endSite - startSite];
        int count = 0;
        int firstAlign = translateSite(startSite);
        int secondAlign = translateSite(endSite);
        for (int i = firstAlign; i <= secondAlign; i++) {
            int firstSite = 0;
            if (i == firstAlign) {
                firstSite = startSite - mySiteOffsets[firstAlign];
            }
            int secondSite = 0;
            if (firstAlign == secondAlign) {
                secondSite = endSite - mySiteOffsets[firstAlign];
            } else if (i != secondAlign) {
                secondSite = myAlignments[i].numberOfSites();
            } else {
                secondSite = endSite - mySiteOffsets[secondAlign];
            }
            for (int s = firstSite; s < secondSite; s++) {
                result[count++] = myAlignments[i].genotype(taxon, s);
            }
        }
        return result;

    }

    @Override
    public byte genotype(int taxon, Chromosome locus, int physicalPosition) {
        int site = siteOfPhysicalPosition(physicalPosition, locus);
        int translate = translateSite(site);
        return myAlignments[translate].genotype(taxon, site - mySiteOffsets[translate]);
    }

    /**
     * Returns which alignment to use.
     *
     * @param site
     * @return alignment index.
     */
    public int translateSite(int site) {

        for (int i = 1; i < mySiteOffsets.length; i++) {
            if (mySiteOffsets[i] > site) {
                return i - 1;
            }
        }
        throw new IndexOutOfBoundsException("CombineAlignment: translateSite: index out of range: " + site);

    }

    @Override
    public boolean hasReference() {

        for (int i = 0; i < myAlignments.length; i++) {
            if (!myAlignments[i].hasReference()) {
                return false;
            }
        }

        return true;
    }

    @Override
    public String siteName(int site) {
        int translate = translateSite(site);
        return myAlignments[translate].siteName(site - mySiteOffsets[translate]);
    }

    @Override
    public int numberOfSites() {
        return mySiteOffsets[mySiteOffsets.length - 1];
    }

    @Override
    public int getChromosomeSiteCount(Chromosome locus) {
        return ((Alignment) myChromosomes.get(locus)).getChromosomeSiteCount(locus);
    }

    @Override
    public int chromosomalPosition(int site) {
        int translate = translateSite(site);
        return myAlignments[translate].chromosomalPosition(site - mySiteOffsets[translate]);
    }

    @Override
    public int siteOfPhysicalPosition(int physicalPosition, Chromosome locus) {
        Alignment align = ((Alignment) myChromosomes.get(locus));
        int i = -1;
        for (int j = 0; j < myAlignments.length; j++) {
            if (myAlignments[j] == align) {
                i = j;
                break;
            }
        }
        if (i == -1) {
            return -1;
        }
        return mySiteOffsets[i] + align.siteOfPhysicalPosition(physicalPosition, locus);
    }

    @Override
    public int siteOfPhysicalPosition(int physicalPosition, Chromosome locus, String snpName) {
        Alignment align = ((Alignment) myChromosomes.get(locus));
        int i = -1;
        for (int j = 0; j < myAlignments.length; j++) {
            if (myAlignments[j] == align) {
                i = j;
                break;
            }
        }
        if (i == -1) {
            return -1;
        }
        return mySiteOffsets[i] + align.siteOfPhysicalPosition(physicalPosition, locus, snpName);
    }

    @Override
    public Chromosome getChromosome(int site) {
        int translate = translateSite(site);
        return myAlignments[translate].getChromosome(site - mySiteOffsets[translate]);
    }

    @Override
    public Chromosome[] getChromosomes() {
        return myChromosomesList;
    }

    @Override
    public int numChromosomes() {
        if (myChromosomesList == null) {
            return 0;
        } else {
            return myChromosomesList.length;
        }
    }

    @Override
    public float[][] siteScores() {

        if (!hasSiteScores()) {
            return null;
        }

        int numSeqs = numberOfTaxa();
        float[][] result = new float[numSeqs][numberOfSites()];
        for (int a = 0, n = myAlignments.length; a < n; a++) {
            if (myAlignments[a].hasSiteScores()) {
                for (int s = 0, m = myAlignments[a].numberOfSites(); s < m; s++) {
                    for (int t = 0; t < numSeqs; t++) {
                        result[t][mySiteOffsets[a] + s] = myAlignments[a].siteScore(t, s);
                    }
                }
            }
        }

        return result;

    }

    @Override
    public float siteScore(int taxon, int site) {
        int translate = translateSite(site);
        return myAlignments[translate].siteScore(taxon, site - mySiteOffsets[translate]);
    }

    @Override
    public boolean hasSiteScores() {
        for (Alignment align : myAlignments) {
            if (align.hasSiteScores()) {
                return true;
            }
        }
        return false;
    }

    @Override
    public int getIndelSize(int site) {
        int translate = translateSite(site);
        return myAlignments[translate].getIndelSize(site - mySiteOffsets[translate]);
    }

    @Override
    public boolean isIndel(int site) {
        int translate = translateSite(site);
        return myAlignments[translate].isIndel(site - mySiteOffsets[translate]);
    }

    @Override
    public byte referenceGenotype(int site) {
        int translate = translateSite(site);
        return myAlignments[translate].referenceGenotype(site - mySiteOffsets[translate]);
    }

    @Override
    public Alignment[] getAlignments() {
        return myAlignments;
    }

    @Override
    public byte getMajorAllele(int site) {
        int translate = translateSite(site);
        return myAlignments[translate].getMajorAllele(site - mySiteOffsets[translate]);
    }

    @Override
    public byte getMinorAllele(int site) {
        int translate = translateSite(site);
        return myAlignments[translate].getMinorAllele(site - mySiteOffsets[translate]);
    }

    @Override
    public byte[] getMinorAlleles(int site) {
        int translate = translateSite(site);
        return myAlignments[translate].getMinorAlleles(site - mySiteOffsets[translate]);
    }

    @Override
    public byte[] getAlleles(int site) {
        int translate = translateSite(site);
        return myAlignments[translate].getAlleles(site - mySiteOffsets[translate]);
    }

    @Override
    public double getMinorAlleleFrequency(int site) {
        int translate = translateSite(site);
        return myAlignments[translate].getMinorAlleleFrequency(site - mySiteOffsets[translate]);
    }

    @Override
    public int[][] allelesSortedByFrequency(int site) {
        int translate = translateSite(site);
        return myAlignments[translate].allelesSortedByFrequency(site - mySiteOffsets[translate]);
    }

    @Override
    public byte[] genotypeArray(int taxon, int site) {
        int translate = translateSite(site);
        return myAlignments[translate].genotypeArray(taxon, site - mySiteOffsets[translate]);
    }

    @Override
    public byte[] genotypeAllSites(int taxon) {
        byte[] result = new byte[numberOfSites()];
        for (int i = 0; i < myAlignments.length; i++) {
            byte[] current = myAlignments[i].genotypeAllSites(taxon);
            System.arraycopy(current, 0, result, myChromosomesOffsets[i], current.length);
        }
        return result;
    }

    @Override
    public BitSet allelePresenceForAllSites(int taxon, int alleleNumber) {
        throw new UnsupportedOperationException("CombineAlignment: getAllelePresenceForAllSites: This operation isn't possible as it spans multiple AlignmentNews. It needs to be optimized for taxa first.");
    }

    @Override
    public BitSet allelePresenceForAllTaxa(int site, int alleleNumber) {
        int translate = translateSite(site);
        return myAlignments[translate].allelePresenceForAllTaxa(site - mySiteOffsets[translate], alleleNumber);
    }

    @Override
    public long[] allelePresenceForSitesBlock(int taxon, int alleleNumber, int startBlock, int endBlock) {
        throw new UnsupportedOperationException("CombineAlignment: getAllelePresenceForSitesBlock: This operation isn't possible as it spans multiple AlignmentNews. It needs to be optimized for taxa first.");
    }

    @Override
    public String genotypeAsString(int taxon, int site) {
        int translate = translateSite(site);
        return myAlignments[translate].genotypeAsString(taxon, site - mySiteOffsets[translate]);
    }

    @Override
    public String[] genotypeAsStringArray(int taxon, int site) {
        int translate = translateSite(site);
        return myAlignments[translate].genotypeAsStringArray(taxon, site - mySiteOffsets[translate]);
    }

    @Override
    public byte[] referenceGenotypes(int startSite, int endSite) {
        int numSites = endSite - startSite;
        byte[] result = new byte[numSites];
        for (int i = 0; i < numSites; i++) {
            result[i] = referenceGenotype(startSite + i);
        }
        return result;
    }

    @Override
    public byte[] referenceGenotypeForAllSites() {

        for (int i = 0; i < myAlignments.length; i++) {
            if (!myAlignments[i].hasReference()) {
                return null;
            }
        }

        byte[] result = new byte[numberOfSites()];
        int count = 0;
        for (int i = 0; i < myAlignments.length; i++) {
            byte[] current = myAlignments[i].referenceGenotypeForAllSites();
            for (int j = 0; j < current.length; j++) {
                result[count++] = current[j];
            }
        }
        return result;

    }

    @Override
    public boolean isHeterozygous(int taxon, int site) {
        int translate = translateSite(site);
        return myAlignments[translate].isHeterozygous(taxon, site - mySiteOffsets[translate]);
    }

    @Override
    public int[] physicalPositions() {

        boolean allNull = true;
        for (int i = 0; i < myAlignments.length; i++) {
            int[] current = myAlignments[0].physicalPositions();
            if ((current != null) && (current.length != 0)) {
                allNull = false;
                break;
            }
        }

        if (allNull) {
            return null;
        } else {
            int[] result = new int[numberOfSites()];
            int count = 0;
            for (int i = 0; i < myAlignments.length; i++) {
                int[] current = myAlignments[i].physicalPositions();
                for (int j = 0; j < current.length; j++) {
                    result[count++] = current[j];
                }
            }
            return result;
        }
    }

    @Override
    public String getChromosomeName(int site) {
        int translate = translateSite(site);
        return myAlignments[translate].getChromosomeName(site - mySiteOffsets[translate]);
    }

    @Override
    public int[] getChromosomesOffsets() {
        return myChromosomesOffsets;
    }

    @Override
    public SITE_SCORE_TYPE siteScoreType() {
        SITE_SCORE_TYPE first = myAlignments[0].siteScoreType();
        for (int i = 1; i < myAlignments.length; i++) {
            if (first != myAlignments[i].siteScoreType()) {
                return SITE_SCORE_TYPE.MixedScoreTypes;
            }
        }
        return first;
    }

    @Override
    public boolean isAllPolymorphic() {
        for (int i = 0; i < myAlignments.length; i++) {
            if (!myAlignments[i].isAllPolymorphic()) {
                return false;
            }
        }
        return true;
    }

    @Override
    public boolean isPolymorphic(int site) {
        int translate = translateSite(site);
        return myAlignments[translate].isPolymorphic(site - mySiteOffsets[translate]);
    }

    @Override
    public double getMajorAlleleFrequency(int site) {
        int translate = translateSite(site);
        return myAlignments[translate].getMajorAlleleFrequency(site - mySiteOffsets[translate]);
    }

    @Override
    public String genomeVersion() {
        String first = myAlignments[0].genomeVersion();
        if (first == null) {
            return null;
        }
        for (int i = 1; i < myAlignments.length; i++) {
            String current = myAlignments[i].genomeVersion();
            if ((current != null) && (!first.equals(current))) {
                return null;
            }
        }
        return first;
    }

    @Override
    public boolean isPositiveStrand(int site) {
        int translate = translateSite(site);
        return myAlignments[translate].isPositiveStrand(site - mySiteOffsets[translate]);
    }

    @Override
    public boolean isPhased() {
        for (int i = 0; i < myAlignments.length; i++) {
            if (myAlignments[i].isPhased() == false) {
                return false;
            }
        }
        return true;
    }

    @Override
    public GeneticMap getGeneticMap() {
        GeneticMap result = myAlignments[0].getGeneticMap();
        for (int i = 1; i < myAlignments.length; i++) {
            GeneticMap current = myAlignments[i].getGeneticMap();
            if ((current == null) || (!current.equals(result))) {
                return null;
            }
        }
        return result;
    }

    @Override
    public boolean retainsRareAlleles() {
        for (int i = 0; i < myAlignments.length; i++) {
            if (myAlignments[i].retainsRareAlleles() == false) {
                return false;
            }
        }
        return true;
    }

    @Override
    public String[][] alleleDefinitions() {

        if (myAlleleStates != null) {
            return myAlleleStates;
        }

        boolean allTheSame = true;
        String[][] encodings = myAlignments[0].alleleDefinitions();
        if (encodings.length == 1) {
            for (int i = 1; i < myAlignments.length; i++) {
                String[][] current = myAlignments[i].alleleDefinitions();
                if ((current.length == 1) && (encodings[0].length == current[0].length)) {
                    for (int j = 0; j < encodings[0].length; j++) {
                        if (!current[0][j].equals(encodings[0][j])) {
                            allTheSame = false;
                            break;
                        }
                    }
                } else {
                    allTheSame = false;
                    break;
                }

                if (!allTheSame) {
                    break;
                }
            }
        } else {
            allTheSame = false;
        }

        if (allTheSame) {
            myAlleleStates = encodings;
        } else {
            String[][] result = new String[numberOfSites()][];
            int count = 0;
            for (int i = 0; i < myAlignments.length; i++) {
                for (int j = 0, n = myAlignments[i].numberOfSites(); j < n; j++) {
                    result[count++] = myAlignments[i].alleleDefinitions(j);
                }
            }
            myAlleleStates = result;
        }

        return myAlleleStates;

    }

    @Override
    public String[] alleleDefinitions(int site) {
        int translate = translateSite(site);
        return myAlignments[translate].alleleDefinitions(site - mySiteOffsets[translate]);
    }

    @Override
    public String genotypeAsString(int site, byte value) {
        int translate = translateSite(site);
        return myAlignments[translate].genotypeAsString(site - mySiteOffsets[translate], value);
    }

    @Override
    public int maxNumAlleles() {
        int result = 999999;
        for (int i = 0; i < myAlignments.length; i++) {
            if (myAlignments[i].maxNumAlleles() < result) {
                result = myAlignments[i].maxNumAlleles();
            }
        }
        return result;
    }

    @Override
    public int totalGametesNonMissingForSite(int site) {
        int translate = translateSite(site);
        return myAlignments[translate].totalGametesNonMissingForSite(site - mySiteOffsets[translate]);
    }

    @Override
    public int getHeterozygousCount(int site) {
        int translate = translateSite(site);
        return myAlignments[translate].getHeterozygousCount(site - mySiteOffsets[translate]);
    }

    @Override
    public int getMinorAlleleCount(int site) {
        int translate = translateSite(site);
        return myAlignments[translate].getMinorAlleleCount(site - mySiteOffsets[translate]);
    }

    @Override
    public int getMajorAlleleCount(int site) {
        int translate = translateSite(site);
        return myAlignments[translate].getMajorAlleleCount(site - mySiteOffsets[translate]);
    }

    @Override
    public Object[][] genosSortedByFrequency(int site) {
        int translate = translateSite(site);
        return myAlignments[translate].genosSortedByFrequency(site - mySiteOffsets[translate]);
    }

    @Override
    public byte[] getAllelesByScope(ALLELE_SORT_TYPE scope, int site) {
        int translate = translateSite(site);
        return myAlignments[translate].getAllelesByScope(scope, site - mySiteOffsets[translate]);
    }

    @Override
    public BitSet allelePresenceForAllTaxaBySortType(ALLELE_SORT_TYPE type, int site, int alleleNumber) {
        int translate = translateSite(site);
        return myAlignments[translate].allelePresenceForAllTaxaBySortType(type, site - mySiteOffsets[translate], alleleNumber);
    }

    @Override
    public BitSet haplotypeAllelePresenceForAllSites(int taxon, boolean firstParent, int alleleNumber) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public BitSet haplotypeAllelePresenceForAllTaxa(int site, boolean firstParent, int alleleNumber) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public long[] haplotypeAllelePresenceForSitesBlock(int taxon, boolean firstParent, int alleleNumber, int startBlock, int endBlock) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String genotypeAsStringRange(int taxon, int startSite, int endSite) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String genotypeAsStringRow(int taxon) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int[] getStartAndEndOfChromosome(Chromosome chromosome) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int numberOfTaxa() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public Chromosome getChromosome(String name) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String getMajorAlleleAsString(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String getMinorAlleleAsString(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public TaxaList taxa() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String taxaName(int index) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String getDiploidAsString(int site, byte value) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int totalNonMissingForSite(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public Object[][] genoCounts() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public Object[][] getMajorMinorCounts() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int totalGametesNonMissingForTaxon(int taxon) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int getHeterozygousCountForTaxon(int taxon) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int totalNonMissingForTaxon(int taxon) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public byte[] getDepthForAlleles(int taxon, int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public BitStorage getBitStorage(ALLELE_SORT_TYPE scopeType) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public PositionList positionList() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public Genotype getGenotypeMatrix() {
        throw new UnsupportedOperationException("Not supported yet.");
    }
}
