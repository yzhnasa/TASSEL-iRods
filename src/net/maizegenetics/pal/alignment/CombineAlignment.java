/*
 * CombineAlignment
 */
package net.maizegenetics.pal.alignment;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import net.maizegenetics.pal.alignment.bit.BitStorage;
import net.maizegenetics.pal.alignment.genotype.Genotype;
import net.maizegenetics.pal.position.Chromosome;
import net.maizegenetics.pal.position.PositionList;
import net.maizegenetics.pal.taxa.TaxaList;
import net.maizegenetics.pal.taxa.TaxaListUtils;
import net.maizegenetics.util.BitSet;

/**
 *
 * @author terry
 */
public class CombineAlignment implements Alignment {

    private static final long serialVersionUID = -5197800047652332969L;
    private final Alignment[] myAlignments;
    private final int[] mySiteOffsets;
    private final Map myChromosomes = new HashMap();
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
            count = alignments[i].getSiteCount() + count;
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

        TaxaList firstGroup = alignments[0].getTaxaList();
        for (int i = 1; i < alignments.length; i++) {
            if (!areTaxaListsEqual(firstGroup, alignments[i].getTaxaList())) {
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
            groups[i] = alignments[i].getTaxaList();
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

        List offsets = new ArrayList();
        List<Chromosome> chromosomes = new ArrayList();
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

    public byte getBase(int taxon, int site) {
        int translate = translateSite(site);
        return myAlignments[translate].getBase(taxon, site - mySiteOffsets[translate]);
    }

    @Override
    public byte[] getBaseRange(int taxon, int startSite, int endSite) {

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
                secondSite = myAlignments[i].getSiteCount();
            } else {
                secondSite = endSite - mySiteOffsets[secondAlign];
            }
            for (int s = firstSite; s < secondSite; s++) {
                result[count++] = myAlignments[i].getBase(taxon, s);
            }
        }
        return result;

    }

    @Override
    public byte getBase(int taxon, Chromosome locus, int physicalPosition) {
        int site = getSiteOfPhysicalPosition(physicalPosition, locus);
        int translate = translateSite(site);
        return myAlignments[translate].getBase(taxon, site - mySiteOffsets[translate]);
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
    public String getSNPID(int site) {
        int translate = translateSite(site);
        return myAlignments[translate].getSNPID(site - mySiteOffsets[translate]);
    }

    @Override
    public int getSiteCount() {
        return mySiteOffsets[mySiteOffsets.length - 1];
    }

    @Override
    public int getChromosomeSiteCount(Chromosome locus) {
        return ((Alignment) myChromosomes.get(locus)).getChromosomeSiteCount(locus);
    }

    @Override
    public int getPositionInChromosome(int site) {
        int translate = translateSite(site);
        return myAlignments[translate].getPositionInChromosome(site - mySiteOffsets[translate]);
    }

    @Override
    public int getSiteOfPhysicalPosition(int physicalPosition, Chromosome locus) {
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
        return mySiteOffsets[i] + align.getSiteOfPhysicalPosition(physicalPosition, locus);
    }

    @Override
    public int getSiteOfPhysicalPosition(int physicalPosition, Chromosome locus, String snpID) {
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
        return mySiteOffsets[i] + align.getSiteOfPhysicalPosition(physicalPosition, locus, snpID);
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
    public int getNumChromosomes() {
        if (myChromosomesList == null) {
            return 0;
        } else {
            return myChromosomesList.length;
        }
    }

    @Override
    public float[][] getSiteScores() {

        if (!hasSiteScores()) {
            return null;
        }

        int numSeqs = getTaxaCount();
        float[][] result = new float[numSeqs][getSiteCount()];
        for (int a = 0, n = myAlignments.length; a < n; a++) {
            if (myAlignments[a].hasSiteScores()) {
                for (int s = 0, m = myAlignments[a].getSiteCount(); s < m; s++) {
                    for (int t = 0; t < numSeqs; t++) {
                        result[t][mySiteOffsets[a] + s] = myAlignments[a].getSiteScore(t, s);
                    }
                }
            }
        }

        return result;

    }

    @Override
    public float getSiteScore(int seq, int site) {
        int translate = translateSite(site);
        return myAlignments[translate].getSiteScore(seq, site - mySiteOffsets[translate]);
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
    public byte getReferenceAllele(int site) {
        int translate = translateSite(site);
        return myAlignments[translate].getReferenceAllele(site - mySiteOffsets[translate]);
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
    public int[][] getAllelesSortedByFrequency(int site) {
        int translate = translateSite(site);
        return myAlignments[translate].getAllelesSortedByFrequency(site - mySiteOffsets[translate]);
    }

    @Override
    public byte[] getBaseArray(int taxon, int site) {
        int translate = translateSite(site);
        return myAlignments[translate].getBaseArray(taxon, site - mySiteOffsets[translate]);
    }

    @Override
    public byte[] getBaseRow(int taxon) {
        byte[] result = new byte[getSiteCount()];
        for (int i = 0; i < myAlignments.length; i++) {
            byte[] current = myAlignments[i].getBaseRow(taxon);
            System.arraycopy(current, 0, result, myChromosomesOffsets[i], current.length);
        }
        return result;
    }

    @Override
    public BitSet getAllelePresenceForAllSites(int taxon, int alleleNumber) {
        throw new UnsupportedOperationException("CombineAlignment: getAllelePresenceForAllSites: This operation isn't possible as it spans multiple AlignmentNews. It needs to be optimized for taxa first.");
    }

    @Override
    public BitSet getAllelePresenceForAllTaxa(int site, int alleleNumber) {
        int translate = translateSite(site);
        return myAlignments[translate].getAllelePresenceForAllTaxa(site - mySiteOffsets[translate], alleleNumber);
    }

    @Override
    public long[] getAllelePresenceForSitesBlock(int taxon, int alleleNumber, int startBlock, int endBlock) {
        throw new UnsupportedOperationException("CombineAlignment: getAllelePresenceForSitesBlock: This operation isn't possible as it spans multiple AlignmentNews. It needs to be optimized for taxa first.");
    }

    @Override
    public String getBaseAsString(int taxon, int site) {
        int translate = translateSite(site);
        return myAlignments[translate].getBaseAsString(taxon, site - mySiteOffsets[translate]);
    }

    @Override
    public String[] getBaseAsStringArray(int taxon, int site) {
        int translate = translateSite(site);
        return myAlignments[translate].getBaseAsStringArray(taxon, site - mySiteOffsets[translate]);
    }

    @Override
    public byte[] getReference(int startSite, int endSite) {
        int numSites = endSite - startSite;
        byte[] result = new byte[numSites];
        for (int i = 0; i < numSites; i++) {
            result[i] = getReferenceAllele(startSite + i);
        }
        return result;
    }

    @Override
    public byte[] getReference() {

        for (int i = 0; i < myAlignments.length; i++) {
            if (!myAlignments[i].hasReference()) {
                return null;
            }
        }

        byte[] result = new byte[getSiteCount()];
        int count = 0;
        for (int i = 0; i < myAlignments.length; i++) {
            byte[] current = myAlignments[i].getReference();
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
    public int[] getPhysicalPositions() {

        boolean allNull = true;
        for (int i = 0; i < myAlignments.length; i++) {
            int[] current = myAlignments[0].getPhysicalPositions();
            if ((current != null) && (current.length != 0)) {
                allNull = false;
                break;
            }
        }

        if (allNull) {
            return null;
        } else {
            int[] result = new int[getSiteCount()];
            int count = 0;
            for (int i = 0; i < myAlignments.length; i++) {
                int[] current = myAlignments[i].getPhysicalPositions();
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
    public SITE_SCORE_TYPE getSiteScoreType() {
        SITE_SCORE_TYPE first = myAlignments[0].getSiteScoreType();
        for (int i = 1; i < myAlignments.length; i++) {
            if (first != myAlignments[i].getSiteScoreType()) {
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
    public String getGenomeAssembly() {
        String first = myAlignments[0].getGenomeAssembly();
        if (first == null) {
            return null;
        }
        for (int i = 1; i < myAlignments.length; i++) {
            String current = myAlignments[i].getGenomeAssembly();
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
    public String[][] getAlleleEncodings() {

        if (myAlleleStates != null) {
            return myAlleleStates;
        }

        boolean allTheSame = true;
        String[][] encodings = myAlignments[0].getAlleleEncodings();
        if (encodings.length == 1) {
            for (int i = 1; i < myAlignments.length; i++) {
                String[][] current = myAlignments[i].getAlleleEncodings();
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
            String[][] result = new String[getSiteCount()][];
            int count = 0;
            for (int i = 0; i < myAlignments.length; i++) {
                for (int j = 0, n = myAlignments[i].getSiteCount(); j < n; j++) {
                    result[count++] = myAlignments[i].getAlleleEncodings(j);
                }
            }
            myAlleleStates = result;
        }

        return myAlleleStates;

    }

    @Override
    public String[] getAlleleEncodings(int site) {
        int translate = translateSite(site);
        return myAlignments[translate].getAlleleEncodings(site - mySiteOffsets[translate]);
    }

    @Override
    public String getBaseAsString(int site, byte value) {
        int translate = translateSite(site);
        return myAlignments[translate].getBaseAsString(site - mySiteOffsets[translate], value);
    }

    @Override
    public int getMaxNumAlleles() {
        int result = 999999;
        for (int i = 0; i < myAlignments.length; i++) {
            if (myAlignments[i].getMaxNumAlleles() < result) {
                result = myAlignments[i].getMaxNumAlleles();
            }
        }
        return result;
    }

    @Override
    public int getTotalGametesNotMissing(int site) {
        int translate = translateSite(site);
        return myAlignments[translate].getTotalGametesNotMissing(site - mySiteOffsets[translate]);
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
    public Object[][] getDiploidsSortedByFrequency(int site) {
        int translate = translateSite(site);
        return myAlignments[translate].getDiploidsSortedByFrequency(site - mySiteOffsets[translate]);
    }

    @Override
    public byte[] getAllelesByScope(ALLELE_SCOPE_TYPE scope, int site) {
        int translate = translateSite(site);
        return myAlignments[translate].getAllelesByScope(scope, site - mySiteOffsets[translate]);
    }

    @Override
    public BitSet getAllelePresenceForAllTaxaByScope(ALLELE_SCOPE_TYPE scope, int site, int alleleNumber) {
        int translate = translateSite(site);
        return myAlignments[translate].getAllelePresenceForAllTaxaByScope(scope, site - mySiteOffsets[translate], alleleNumber);
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
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String getBaseAsStringRow(int taxon) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int[] getStartAndEndOfChromosome(Chromosome chromosome) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int getSequenceCount() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int getTaxaCount() {
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
    public TaxaList getTaxaList() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String getTaxaName(int index) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String getFullTaxaName(int index) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String getDiploidAsString(int site, byte value) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int getTotalNotMissing(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public Object[][] getDiploidCounts() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public Object[][] getMajorMinorCounts() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int getTotalGametesNotMissingForTaxon(int taxon) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int getHeterozygousCountForTaxon(int taxon) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int getTotalNotMissingForTaxon(int taxon) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public byte[] getDepthForAlleles(int taxon, int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public BitStorage getBitStorage(ALLELE_SCOPE_TYPE scopeType) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public PositionList getPositionList() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public Genotype getGenotypeMatrix() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public byte[] getMajorAlleleForAllSites() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public byte[] getMinorAlleleForAllSites() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public byte[] getGenotypeForAllSites(int taxon) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public byte[] getGenotypeForSiteRange(int taxon, int start, int end) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public byte[] getGenotypeForAllTaxa(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public void transposeData(boolean siteInnerLoop) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
}
