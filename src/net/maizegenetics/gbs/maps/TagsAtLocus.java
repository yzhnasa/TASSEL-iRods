/*
 * TagsAtLocus
 */
package net.maizegenetics.gbs.maps;

import org.apache.commons.math.distribution.BinomialDistributionImpl;

import java.util.ArrayList;
import java.util.Arrays;

import net.maizegenetics.gbs.tagdist.TagsByTaxa;
import net.maizegenetics.gbs.util.BaseEncoder;
import net.maizegenetics.pal.alignment.Alignment;

import java.util.List;

import org.biojava3.core.sequence.DNASequence;
import org.biojava3.core.sequence.compound.AmbiguityDNACompoundSet;
import org.biojava3.alignment.template.Profile;
import org.biojava3.alignment.Alignments;
import org.biojava3.alignment.SimpleGapPenalty;
import org.biojava3.alignment.SubstitutionMatrixHelper;
import org.biojava3.alignment.template.SequencePair;
import org.biojava3.alignment.template.SubstitutionMatrix;
import org.biojava3.alignment.Alignments.PairwiseSequenceAlignerType;
import org.biojava3.core.sequence.compound.NucleotideCompound;

import net.maizegenetics.pal.ids.SimpleIdGroup;
import net.maizegenetics.pal.alignment.Locus;
import net.maizegenetics.pal.alignment.BitAlignment;
import net.maizegenetics.pal.alignment.AlignmentUtils;
import net.maizegenetics.pal.alignment.NucleotideAlignmentConstants;
import net.maizegenetics.prefs.TasselPrefs;

/**
 *
 * @author jcg233
 */
public class TagsAtLocus {

    ArrayList<SingleTagByTaxa> theTags = new ArrayList<SingleTagByTaxa>();
    private int minStartPosition;
    private int maxStartPosition;
    private int chromosome;
    private byte strand;
    private int indexOfRef;
    private int[] tagIndices = null;  // redirect from aligned tag indices to index in theTags
    private int[] positionsOfVariableSites;
    private byte[][] callsAtVariableSitesByTag;
    private final static int maxSNPsPerLocus = 64;
    private final static int maxAlignmentSize = 150;
    private final static double errorRate = 0.01;
    private final static int maxCountAtGeno = 500;
    private SubstitutionMatrix<NucleotideCompound> subMatrix = SubstitutionMatrixHelper.getNuc4_4();
    private SimpleGapPenalty gapPen = new SimpleGapPenalty((short) 5, (short) 2);
    private final static int[] likelihoodRatioThreshAlleleCnt = new int[maxCountAtGeno];  // index = sample size; value = min count of less tagged allele for likelihood ratio > 1
    // if less tagged allele has counts < likelihoodRatioThreshAlleleCnt[totalCount], call it a homozygote
    // where likelihood ratio = (binomial likelihood het) / (binomial likelihood all less tagged alleles are errors)

    static {   // initialize the likelihood ratio cutoffs for quantitative SNP calling
        System.out.println("\n\nInitializing the cutoffs for quantitative SNP calling likelihood ratio (pHet/pErr) >1\n");
        System.out.println("totalReadsForSNPInIndiv\tminLessTaggedAlleleCountForHet");
        for (int trials = 0; trials < 2; ++trials) {
            likelihoodRatioThreshAlleleCnt[trials] = 1;
        }
        int lastThresh = 1;
        for (int trials = 2; trials < likelihoodRatioThreshAlleleCnt.length; ++trials) {
            BinomialDistributionImpl binomHet = new BinomialDistributionImpl(trials, 0.5);
            BinomialDistributionImpl binomErr = new BinomialDistributionImpl(trials, errorRate);
            double LikeRatio;
            try {
                LikeRatio = binomHet.cumulativeProbability(lastThresh) / (1 - binomErr.cumulativeProbability(lastThresh) + binomErr.probability(lastThresh));
                while (LikeRatio <= 1.0) {
                    ++lastThresh;
                    LikeRatio = binomHet.cumulativeProbability(lastThresh) / (1 - binomErr.cumulativeProbability(lastThresh) + binomErr.probability(lastThresh));
                }
                likelihoodRatioThreshAlleleCnt[trials] = lastThresh;
                System.out.println(trials + "\t" + lastThresh);
            } catch (Exception e) {
                System.err.println("Error in the TagsAtLocus.BinomialDistributionImpl");
            }
        }
        System.out.println("\n");
    }

    public TagsAtLocus(int chromosome, byte strand, int startPosition, boolean includeRefGenome) {
        this.chromosome = chromosome;
        this.strand = includeRefGenome ? 1 : strand;
        this.minStartPosition = startPosition;
        this.maxStartPosition = startPosition;
        positionsOfVariableSites = null;
        callsAtVariableSitesByTag = null;
    }

    public void addTag(int tagTOPMIndex, TagsOnPhysicalMap theTOPM, TagsByTaxa theTBT, boolean includeRefGenome) {
        SingleTagByTaxa singleTBT = new SingleTagByTaxa(tagTOPMIndex, theTOPM, theTBT, includeRefGenome);
        if (singleTBT.taxaWithTag > 0) {
            theTags.add(singleTBT);
            if (singleTBT.startPosition > minStartPosition) {
                maxStartPosition = singleTBT.startPosition;
            }
        }
    }

    public int getSize() {
        return theTags.size();
    }

    public int getChromosome() {
        return chromosome;
    }

    public byte getStrand() {
        return strand;
    }

    public int getMinStartPosition() {
        return minStartPosition;
    }

    public int getMaxStartPosition() {
        return maxStartPosition;
    }

    public void setMinStartPosition(int newMinStartPosition) {
        minStartPosition = newMinStartPosition;
    }

    public int getTOPMIndexOfTag(int tagIndex) {
        return theTags.get(tagIndex).tagTOPMIndex;
    }

    public byte getCallAtVariableSiteForTag(int site, int tagIndex) {
        return callsAtVariableSitesByTag[site][tagIndex];
    }

    public int getNumberTaxaCovered() {
        if (theTags.size() < 1) {
            return 0;
        }
        int nTaxaCovered = 0;
        boolean[] covered = new boolean[theTags.get(0).tagDist.length];  // initializes to false
        for (SingleTagByTaxa sTBT : theTags) {
            for (int tx = 0; tx < covered.length; ++tx) {
                if (sTBT.tagDist[tx] > 0) {
                    covered[tx] = true;
                }
            }
        }
        for (int tx = 0; tx < covered.length; ++tx) {
            if (covered[tx]) {
                ++nTaxaCovered;
            }
        }
        return nTaxaCovered;
    }

    private void assignRefTag() {
        indexOfRef = Integer.MIN_VALUE;
        int lengthOfRef = Integer.MIN_VALUE;
        int tagIndex = 0;
        for (SingleTagByTaxa sTBT : theTags) {
            if (sTBT.divergence == 0 && sTBT.tagLength > lengthOfRef) {
                indexOfRef = tagIndex;
                lengthOfRef = sTBT.tagLength;
            }
            ++tagIndex;
        }
    }

    public byte[][] getSNPCallsQuant(boolean callBiallelicSNPsWithGap) {
        if (theTags.size() < 2) {
            return null;
        }
        Alignment tagAlignment = this.getVariableSites();
        if (tagAlignment == null || tagAlignment.getSiteCount() < 1) {
            return null;
        }
        int nSites = tagAlignment.getSiteCount();
        int nTaxa = theTags.get(0).tagDist.length;
        if (nTaxa < 1) {
            return null;
        }
        byte[][] callsBySite = new byte[nSites][nTaxa];
        final int nAlignedTags = tagAlignment.getSequenceCount();
        tagIndices = new int[nAlignedTags];
        callsAtVariableSitesByTag = new byte[nSites][theTags.size()];
        for (int tg = 0; tg < nAlignedTags; tg++) {
            tagIndices[tg] = Integer.parseInt(tagAlignment.getTaxaName(tg).split("_")[0]);  // taxaName in tagAlignment is set to indexInTheTags_"refSeq"|"no"
            for (int s = 0; s < nSites; s++) {
                callsAtVariableSitesByTag[s][tagIndices[tg]] = tagAlignment.getBase(tg, s);
            }
        }
        positionsOfVariableSites = new int[nSites];
        for (int s = 0; s < nSites; s++) {
            positionsOfVariableSites[s] = tagAlignment.getPositionInLocus(s);
            for (int tx = 0; tx < nTaxa; tx++) {
                int[] alleleCounts = new int[Byte.MAX_VALUE];
                for (int tg = 0; tg < nAlignedTags; tg++) {
                    int tagIndex = tagIndices[tg];
                    byte baseToAdd = callsAtVariableSitesByTag[s][tagIndex];
                    if (baseToAdd == Alignment.UNKNOWN_DIPLOID_ALLELE && callBiallelicSNPsWithGap) {
                        baseToAdd = NucleotideAlignmentConstants.GAP_DIPLOID_ALLELE;
                    }
                    alleleCounts[baseToAdd] += theTags.get(tagIndex).tagDist[tx];
                }
                callsBySite[s][tx] = resolveQuantGeno(alleleCounts);
            }
        }
        return callsBySite;
    }

    public byte[][] getSNPCallsQuant(String refSeq, boolean callBiallelicSNPsWithGap) {
        if (theTags.size() < 2) {
            return null;
        }
        Alignment tagAlignment = this.getVariableSites(refSeq);
        if (tagAlignment == null || tagAlignment.getSiteCount() < 1) {
            return null;
        }
        int nSites = tagAlignment.getSiteCount();
        int nTaxa = theTags.get(0).tagDist.length;  // the number of taxa is the same for all tags
        if (nTaxa < 1) {
            return null;
        }
        byte[][] callsBySite = new byte[nSites][nTaxa];
        final int nAlignedTags = tagAlignment.getSequenceCount();
        tagIndices = new int[nAlignedTags];  // the reference sequence is not included
        callsAtVariableSitesByTag = new byte[nSites][theTags.size()];
        for (int tg = 0; tg < nAlignedTags; tg++) {
            int indexInTheTags = Integer.parseInt(tagAlignment.getTaxaName(tg)); // taxaName in tagAlignment is set to indexInTheTags
            tagIndices[tg] = indexInTheTags;
            for (int s = 0; s < nSites; s++) {
                callsAtVariableSitesByTag[s][tagIndices[tg]] = tagAlignment.getBase(tg, s);
            }
        }
        positionsOfVariableSites = new int[nSites];
        for (int s = 0; s < nSites; s++) {
            positionsOfVariableSites[s] = tagAlignment.getPositionInLocus(s);
            for (int tx = 0; tx < nTaxa; tx++) {
                int[] alleleCounts = new int[Byte.MAX_VALUE];
                for (int tg = 0; tg < nAlignedTags; tg++) {
                    int tagIndex = tagIndices[tg];
                    byte baseToAdd = callsAtVariableSitesByTag[s][tagIndex];
                    if (baseToAdd == Alignment.UNKNOWN_DIPLOID_ALLELE && callBiallelicSNPsWithGap && maxStartPosition == minStartPosition) {
                        baseToAdd = NucleotideAlignmentConstants.GAP_DIPLOID_ALLELE;
                    }
                    alleleCounts[baseToAdd] += theTags.get(tagIndex).tagDist[tx];
                }
                callsBySite[s][tx] = resolveQuantGeno(alleleCounts);
            }
        }
        return callsBySite;
    }

    public int[] getPositionsOfVariableSites() {
        return positionsOfVariableSites;
    }

    private Alignment getVariableSites() {
        if (theTags.size() < 2) {
            return null;
        }
        if (theTags.size() > maxAlignmentSize) {
            return null;   // should we use a maxAlignmentSize (upper limit of tags) here?
        }
        this.assignRefTag();
        boolean checkReplicateProfiles = false;
        boolean printOutRefWithGaps = false;
        List<DNASequence> lst = new ArrayList<DNASequence>();
        int tagIndex = 0;
        for (SingleTagByTaxa sTBT : theTags) {
            DNASequence ds = new DNASequence(sTBT.tagTrimmed);
            String refMark = (tagIndex == indexOfRef) ? "refSeq" : "no";
            ds.setOriginalHeader(tagIndex + "_" + refMark);    // OriginalHeader set to indexInTheTags_'refSeq'|'no'
            ds.setCompoundSet(AmbiguityDNACompoundSet.getDNACompoundSet());
            lst.add(ds);
            ++tagIndex;
        }
        Profile<DNASequence, NucleotideCompound> profile = Alignments.getMultipleSequenceAlignment(lst);
        if (checkReplicateProfiles) {
            System.out.printf("Clustal1:%d%n%s%n", minStartPosition, profile);
            Profile<DNASequence, NucleotideCompound> profile2 = Alignments.getMultipleSequenceAlignment(lst);
            System.out.printf("Clustal2:%d%n%s%n", minStartPosition, profile2);
        }
        String[] aseqs = new String[theTags.size()];
        String[] names = new String[theTags.size()];
        boolean refTagWithGaps = false;
        int[] positions = null;
        for (int i = 0; i < aseqs.length; i++) {
            aseqs[i] = profile.getAlignedSequence(i + 1).getSequenceAsString();
            names[i] = profile.getAlignedSequence(i + 1).getOriginalSequence().getOriginalHeader();
            if (names[i].split("_")[1].equals("refSeq")) {  // names were set to indexInTheTags_"refSeq"|"no"
                if (aseqs[i].contains("-")) {
                    refTagWithGaps = true;
                    positions = new int[aseqs[i].length()];
                    positions[0] = 0;
                    for (int site = 1; site < aseqs[i].length(); site++) {
                        positions[site] = (aseqs[i].charAt(site) == '-') ? (positions[site - 1]) : (positions[site - 1] + 1);
                    }
                }
            }
        }
        profile = null;
        Alignment aa = null;
        if (refTagWithGaps) {
            aa = BitAlignment.getNucleotideInstance(new SimpleIdGroup(names), aseqs, null, null, positions, TasselPrefs.getAlignmentMaxAllelesToRetain(), new Locus[]{Locus.UNKNOWN}, new int[]{0}, null, TasselPrefs.getAlignmentRetainRareAlleles(), true);
        } else {
            aa = BitAlignment.getNucleotideInstance(new SimpleIdGroup(names), aseqs, null, null, null, TasselPrefs.getAlignmentMaxAllelesToRetain(), new Locus[]{Locus.UNKNOWN}, new int[]{0}, null, TasselPrefs.getAlignmentRetainRareAlleles(), true);
        }
        Alignment faa = AlignmentUtils.removeSitesBasedOnFreqIgnoreMissing(aa, 0.000001, 1.0, 2);
        if (printOutRefWithGaps && refTagWithGaps) {
            System.out.println("chr" + chromosome + "  pos:" + minStartPosition + "  strand:" + strand + "  AA\n" + aa.toString().trim());
            System.out.println("chr" + chromosome + "  pos:" + minStartPosition + "  strand:" + strand + "  FAA\n" + faa.toString());
        }
        if (faa.getSiteCount() > maxSNPsPerLocus || faa.getSiteCount() < 1 || faa.getSequenceCount() < 2) {
            return null;
        }
        return faa;
    }

    private Alignment getVariableSites(String refSeq) {
        if (theTags.size() < 2) {
            return null;
        }
        boolean printOutAlignments = true;
        int startRefGenIndex = 0, endRefGenIndex = 1, startTagIndex = 2, endTagIndex = 3; // relevant indices in alignStats[]  (alignedTagLen=4)
        DNASequence dsRefSeq = new DNASequence(refSeq);
        dsRefSeq.setCompoundSet(AmbiguityDNACompoundSet.getDNACompoundSet());
        int minRefGenIndex = Integer.MAX_VALUE, maxRefGenIndex = Integer.MIN_VALUE;
        ArrayList<SequencePair<DNASequence, NucleotideCompound>> pairwiseAligns = new ArrayList<SequencePair<DNASequence, NucleotideCompound>>();
        int tagIndex = 0;
        for (SingleTagByTaxa sTBT : theTags) {
            DNASequence ds = new DNASequence(sTBT.tagTrimmed);
            ds.setCompoundSet(AmbiguityDNACompoundSet.getDNACompoundSet());
            SequencePair<DNASequence, NucleotideCompound> psa = Alignments.getPairwiseAlignment(ds, dsRefSeq, PairwiseSequenceAlignerType.LOCAL, gapPen, subMatrix);
            int[] alignStats = getAlignStats(psa, printOutAlignments, tagIndex, sTBT.tagLength, sTBT.tagStrand);
            minRefGenIndex = adjustMinRefGenIndex(minRefGenIndex, alignStats[startRefGenIndex], alignStats[startTagIndex]);
            maxRefGenIndex = adjustMaxRefGenIndex(maxRefGenIndex, alignStats[endRefGenIndex], alignStats[endTagIndex], sTBT.tagLength, refSeq);
            pairwiseAligns.add(psa);
            ++tagIndex;
        }
        minStartPosition += minRefGenIndex - 1;
        if (printOutAlignments && minStartPosition > 10000000 && minStartPosition < 10100000) {
            System.out.println("minRefGenIndex:" + minRefGenIndex + "  maxRefGenIndex:" + maxRefGenIndex + "  ChrPositionAtMinRefGenIndex:" + minStartPosition + "\n");
        }
        String[] aseqs = new String[theTags.size()];  // omit the reference genome sequence
        String[] names = new String[theTags.size()];
        char[][] myAlign = getAlignment(pairwiseAligns, refSeq, minRefGenIndex, maxRefGenIndex, aseqs, names);
        if (printOutAlignments && minStartPosition > 10000000 && minStartPosition < 10100000) {
            writeAlignment(refSeq, myAlign, minRefGenIndex, maxRefGenIndex);
        }
        Alignment a = null;
        a = BitAlignment.getNucleotideInstance(new SimpleIdGroup(names), aseqs, null, null, null, TasselPrefs.getAlignmentMaxAllelesToRetain(), new Locus[]{Locus.UNKNOWN}, new int[]{0}, null, TasselPrefs.getAlignmentRetainRareAlleles(), true);
        Alignment fa = AlignmentUtils.removeSitesBasedOnFreqIgnoreMissing(a, 0.000001, 1.0, 2);
        if (printOutAlignments && minStartPosition > 10000000 && minStartPosition < 10100000) {
            System.out.println("chr" + chromosome + "  pos:" + minStartPosition + "  strand:" + strand + "  FA (alignment filtered for polymorphic sites):\n" + fa.toString());
        }
        if (fa.getSiteCount() > maxSNPsPerLocus * 5 || fa.getSiteCount() < 1 || fa.getSequenceCount() < 2) {
            return null;
        }
        return fa;
    }

    private int[] getAlignStats(SequencePair<DNASequence, NucleotideCompound> psa, boolean printOutAlignments, int tagIndex, int tagLength, byte tagStrand) {
        int[] alignStats = new int[5];
        int startRefGenIndex = 0, endRefGenIndex = 1, startTagIndex = 2, endTagIndex = 3, alignedTagLen = 4; // indices in alignStats[]
        alignStats[startRefGenIndex] = psa.getIndexInTargetAt(1);
        alignStats[endRefGenIndex] = psa.getIndexInTargetAt(psa.getLength());
        alignStats[startTagIndex] = psa.getIndexInQueryAt(1);
        alignStats[endTagIndex] = psa.getIndexInQueryAt(psa.getLength());
        alignStats[alignedTagLen] = alignStats[endTagIndex] - alignStats[startTagIndex] + 1;
        if (printOutAlignments && minStartPosition > 10000000 && minStartPosition < 10100000) {
            System.out.println("tagIndex:" + tagIndex
                    + "  startRefGenIndex:" + alignStats[startRefGenIndex]
                    + "  endRefGenIndex:" + alignStats[endRefGenIndex]
                    + "  tagLength:" + tagLength
                    + "  startTagIndex:" + alignStats[startTagIndex]
                    + "  endTagIndex:" + alignStats[endTagIndex]
                    + "  alignedTagLen:" + alignStats[alignedTagLen]
                    + "  originalStrand: " + tagStrand + "\n"
                    + psa);
        }
        return alignStats;
    }

    private int adjustMinRefGenIndex(int minRefGenIndex, int startRefGenIndex, int startTagIndex) {
        if (startRefGenIndex < minRefGenIndex) {
            minRefGenIndex = startRefGenIndex;
        }
        if (startTagIndex > 1 && startTagIndex < 4 && startRefGenIndex - startTagIndex + 1 < minRefGenIndex && startRefGenIndex - startTagIndex + 1 > 0) {
            // extend regional alignment if there was soft clipping of 1 or 2 bases and there is sufficient 5' refSeq available
            minRefGenIndex = startRefGenIndex - startTagIndex + 1;
        }
        return minRefGenIndex;
    }

    private int adjustMaxRefGenIndex(int maxRefGenIndex, int endRefGenIndex, int endTagIndex, int tagLength, String refSeq) {
        if (endRefGenIndex > maxRefGenIndex) {
            maxRefGenIndex = endRefGenIndex;
        }
        if (endTagIndex < tagLength
                && tagLength - endTagIndex < 3
                && endRefGenIndex + tagLength - endTagIndex > maxRefGenIndex
                && endRefGenIndex + tagLength - endTagIndex <= refSeq.length()) {
            // extend regional alignment if there was soft clipping of 1 or 2 bases and there is sufficient 3' refSeq available
            maxRefGenIndex = endRefGenIndex + tagLength - endTagIndex;
        }
        return maxRefGenIndex;
    }

    private char[][] getAlignment(ArrayList<SequencePair<DNASequence, NucleotideCompound>> pairwiseAligns,
            String refSeq, int minRefGenIndex, int maxRefGenIndex, String[] aseqs, String[] names) {
        int totAlignedLen = maxRefGenIndex - minRefGenIndex + 1;
        char[][] myAlign = new char[theTags.size()][totAlignedLen];  // omit the reference genome sequence
        for (int t = 0; t < myAlign.length; t++) {
            Arrays.fill(myAlign[t], 'N');
        }
        int tagIndex = 0;
        for (SingleTagByTaxa sTBT : theTags) {
            SequencePair<DNASequence, NucleotideCompound> psa = pairwiseAligns.get(tagIndex);
            int tagStart = psa.getIndexInQueryAt(1);
            int tagEnd = psa.getIndexInQueryAt(psa.getLength());
            int refSeqStart = psa.getIndexInTargetAt(1);
            int refSeqEnd = psa.getIndexInTargetAt(psa.getLength());
            if (tagStart > 1 && tagStart < 4 && refSeqStart - minRefGenIndex - tagStart + 1 > -1) {
                // extend tag start if there was soft clipping of 1 or 2 bases and there is sufficient 5' refSeq available
                for (int offset = tagStart - 1; offset > 0; offset--) {
                    myAlign[tagIndex][refSeqStart - minRefGenIndex - offset] = sTBT.tagTrimmed.charAt(tagStart - offset - 1);
                }
            }
            for (int i = 1; i <= psa.getLength(); i++) {
                char refBase = psa.getCompoundInTargetAt(i).getBase().charAt(0);
                if (refBase != '-') {
                    myAlign[tagIndex][psa.getIndexInTargetAt(i) - minRefGenIndex] = psa.getCompoundInQueryAt(i).getBase().charAt(0);
                }
            }
            int extension = sTBT.tagLength - tagEnd;
            if (extension > 0 && extension < 3 && refSeqEnd - minRefGenIndex + extension < refSeq.length()) {
                // extend tag end if there was soft clipping of 1 or 2 bases and there is sufficient 3' refSeq available
                for (int offset = sTBT.tagLength - tagEnd; offset > 0; offset--) {
                    myAlign[tagIndex][refSeqEnd - minRefGenIndex + offset] = sTBT.tagTrimmed.charAt(sTBT.tagLength - offset);
                }
            }
            aseqs[tagIndex] = new String(myAlign[tagIndex]);
            names[tagIndex] = tagIndex + "";
            ++tagIndex;
        }
        return myAlign;
    }

    private void writeAlignment(String refSeq, char[][] myAlign, int minRefGenIndex, int maxRefGenIndex) {
        System.out.println("All tags in the region aligned to the reference sequence (first line) (insertions relative to the reference excluded):");
        System.out.println(refSeq.substring(minRefGenIndex - 1, maxRefGenIndex));
        for (int tagIndex = 0; tagIndex < myAlign.length; tagIndex++) {
            for (int b = 0; b < myAlign[tagIndex].length; b++) {
                System.out.print(myAlign[tagIndex][b]);
            }
            System.out.print("\n");
        }
        System.out.print("\n");
    }

    private byte resolveQuantGeno(int[] alleleCounts) {
        int[][] sortedAlleleCounts = sortAllelesByCount(alleleCounts);
        int a1Count = sortedAlleleCounts[1][0];
        if (a1Count == 0) {
            return Alignment.UNKNOWN_DIPLOID_ALLELE;
        }
        int a2Count = sortedAlleleCounts[1][1];  // What if a3Count = a2Count? -- this situation is not dealt with
        byte a1 = (byte) sortedAlleleCounts[0][0];
        if (a2Count == 0) {
            return a1;
        }
        byte a2 = (byte) sortedAlleleCounts[0][1];
        int totCount = a1Count + a2Count;
        if (totCount < maxCountAtGeno) {
            if (a2Count < likelihoodRatioThreshAlleleCnt[totCount]) {
                return a1;  // call it a homozygote
            } else {
                return AlignmentUtils.getDiploidValue(a1, a2);  // call it a het
                //return IUPACNucleotides.getDegerateSNPByteFromTwoSNPs(a1, a2); // call it a het
            }
        } else {
            if (a2Count / totCount < 0.1) {
                return a1;  // call it a homozygote
            } else {
                return AlignmentUtils.getDiploidValue(a1, a2);  // call it a het
                //return IUPACNucleotides.getDegerateSNPByteFromTwoSNPs(a1, a2); // call it a het
            }
        }
    }

    private int[][] sortAllelesByCount(int[] alleleCounts) {
        byte[] alleles = {'A', 'C', 'G', 'T', '-'};  // note that 'N' is not included as an allele
        int[][] result = new int[2][alleles.length]; // result[0][i]=allele; result[1][i]=count
        for (int i = 0; i < alleles.length; i++) {
            result[0][i] = alleles[i];
            result[1][i] = alleleCounts[alleles[i]];
        }
        boolean change = true;
        while (change) { // sort the alleles by descending frequency
            change = false;
            for (int k = 0; k < alleles.length - 1; k++) {
                if (result[1][k] < result[1][k + 1]) {
                    int temp = result[0][k];
                    result[0][k] = result[0][k + 1];
                    result[0][k + 1] = temp;
                    int tempCount = result[1][k];
                    result[1][k] = result[1][k + 1];
                    result[1][k + 1] = tempCount;
                    change = true;
                }
            }
        }
        return result;
    }

    private String padTagWithNs(SingleTagByTaxa tag, String refSeq) {
        StringBuilder sb = new StringBuilder();
        char[] nullBases;
        if (tag.tagStrand == -1) {
            nullBases = new char[tag.startPosition - minStartPosition - tag.tagTrimmed.length() + 1];
        } else {
            nullBases = new char[tag.startPosition - minStartPosition];
        }
        Arrays.fill(nullBases, 'N');
        sb.append(nullBases);
        sb.append(tag.tagTrimmed);
        if (tag.tagStrand == -1) {
            nullBases = new char[refSeq.length() - sb.length()];
        } else {
            nullBases = new char[refSeq.length() - sb.length()];
        }

        Arrays.fill(nullBases, 'N');
        sb.append(nullBases);
        return sb.toString();
    }
}

class SingleTagByTaxa {

    int tagTOPMIndex;
    long[] tag;
    byte tagLength;
    int startPosition;
    byte tagStrand;
    int divergence;
    String tagTrimmed;
    int tagTBTIndex; //index in the TBT
    int taxaWithTag;
    byte[] tagDist;  // observed count of the tag for each taxon

    SingleTagByTaxa(int tagTOPMIndex, TagsOnPhysicalMap theTOPM, TagsByTaxa theTBT, boolean includeRefGenome) {
        tagStrand = Byte.MIN_VALUE;
        this.tagTOPMIndex = tagTOPMIndex;
        tag = theTOPM.getTag(tagTOPMIndex);
        tagTBTIndex = theTBT.getTagIndex(tag);
        taxaWithTag = (tagTBTIndex > -1) ? theTBT.getNumberOfTaxaWithTag(tagTBTIndex) : 0;
        if (taxaWithTag > 0) {  // tags with 0 taxaWithTag will not be added to TagsAtLocus
            startPosition = theTOPM.getStartPosition(tagTOPMIndex);
            tagLength = (byte) theTOPM.getTagLength(tagTOPMIndex);
            divergence = theTOPM.getDivergence(tagTOPMIndex);
            tagTrimmed = BaseEncoder.getSequenceFromLong(tag).substring(0, tagLength);
            tagStrand = theTOPM.getStrand(tagTOPMIndex);
            if (includeRefGenome) {
                if (tagStrand == -1) {
                    tagTrimmed = BaseEncoder.getReverseComplement(tagTrimmed);
                }
            } else if (tagLength < theTOPM.getTagSizeInLong() * BaseEncoder.chunkSize) {
                tagTrimmed = tagTrimmed
                        + theTOPM.getNullTag().substring(0, theTOPM.getTagSizeInLong() * BaseEncoder.chunkSize - tagLength).replace("A", "N");
            }
            tagDist = theTBT.getTaxaReadCountsForTag(tagTBTIndex);
        }
    }
}
