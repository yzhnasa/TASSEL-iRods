/*
 * TagsToSNPByAlignmentPlugin
 */
package net.maizegenetics.gbs.pipeline;

import java.awt.Frame;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;

import java.util.Arrays;
import java.util.HashMap;

import javax.swing.ImageIcon;

import net.maizegenetics.gbs.maps.TagsAtLocus;
import net.maizegenetics.gbs.maps.TagsOnPhysicalMap;
import net.maizegenetics.gbs.tagdist.TagsByTaxa;
import net.maizegenetics.gbs.tagdist.TagsByTaxaBitFileMap;
import net.maizegenetics.gbs.tagdist.TagsByTaxaByteFileMap;
import net.maizegenetics.gbs.tagdist.TagsByTaxaByteHDF5TagGroups;
import net.maizegenetics.util.ArgsEngine;
import net.maizegenetics.gbs.util.BaseEncoder;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.AlignmentUtils;
import net.maizegenetics.pal.alignment.ExportUtils;
import net.maizegenetics.pal.alignment.Locus;
import net.maizegenetics.pal.alignment.MutableNucleotideAlignment;
import net.maizegenetics.pal.alignment.NucleotideAlignmentConstants;
import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.SimpleIdGroup;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;

import org.apache.log4j.Logger;

import org.biojava3.core.util.ConcurrencyTools;

/**
 * This class aligns tags at the same physical location against one another, calls SNPs,
 * and then outputs the SNPs to a HapMap file.
 *
 * It is multi-threaded, as there are substantial speed increases with it.
 * @author edbuckler
 */
public class TagsToSNPByAlignmentPlugin extends AbstractPlugin {

    static int maxSize = 200000;  //normally 200K;
    private double minF = -2.0, minMAF = 0.01;
    private int minMAC = 10;
    //    static boolean ignoreTriallelic=false;
    private boolean inclRare = false;  // false = only call the two most common alleles at a site
    private boolean inclGaps = false;  // false = ignore sites where the major or the 1st minor alleles are gaps
    private boolean callBiallelicSNPsWithGap = false;  // true = call sites with a biallelic SNP plus a gap (e.g., A/C/-)
    private boolean isUpdateTOPM = false;
    private boolean useTBTByte = false;
    static double defaultMinPropTaxaWithLocus = 0.1;
    private static Logger myLogger = Logger.getLogger(TagsToSNPByAlignmentPlugin.class);
    TagsOnPhysicalMap theTOPM = null;
    TagsByTaxa theTBT = null;
    File inputFile = null;
    private String inTOPMFile = null;
    private String outTOPMFile = null;
    private boolean usePedigree = false;
    HashMap<String, Double> taxaFs = null;
    boolean[] useTaxaForMinF = null;
    int nInbredTaxa = Integer.MIN_VALUE;
    String outputFilePrefix = null;
    String outHapMap = null;
    int startChr = Integer.MAX_VALUE;
    int endChr = Integer.MIN_VALUE;
    private static ArgsEngine myArgsEngine = null;
    int minTaxaWithLocus;
    private boolean includeReferenceGenome = false;
    private String refGenomeFileStr = null;
    private long[] refGenomeChr = null;
    private boolean fuzzyStartPositions = false;
    int locusBorder = 0;
    final static int chr = 0, str = 1, startPosit = 2;  // indices of these position attributes in array returned by theTOPM.getPositionArray(i)

    public TagsToSNPByAlignmentPlugin() {
        super(null, false);
    }

    public TagsToSNPByAlignmentPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    @Override
    public DataSet performFunction(DataSet input) {
        myLogger.info("Finding SNPs in " + inputFile.getAbsolutePath() + ".");
        myLogger.info(String.format("StartChr:%d EndChr:%d %n", startChr, endChr));
        theTOPM.sortTable(true);
        myLogger.info("\nHere are the first 5 tags in the TOPM (sorted by position):");
        theTOPM.printRows(5, true, true);
        for (int i = startChr; i <= endChr; i++) {
            myLogger.info("\n\nProcessing chromosome " + i + "...");
            String out = outHapMap + ".c" + i;
            myLogger.info("Creating Mutable Alignment to hold genotypes for chr" + i + " (maximum number of sites = " + maxSize + ")");
            MutableNucleotideAlignment theMSA = createMutableAlignment(theTBT, i, i, maxSize + 100);
            if (includeReferenceGenome) {
                refGenomeChr = readReferenceGenomeChr(refGenomeFileStr, i);
            }
            runTagsToSNPByAlignment(theMSA, out, i, false);
            myLogger.info("Finished processing chromosome " + i + "\n\n");
        }
        if (this.isUpdateTOPM) {
            if (outTOPMFile.endsWith(".txt")) {
                theTOPM.writeTextFile(new File(outTOPMFile));
            } else {
                theTOPM.writeBinaryFile(new File(outTOPMFile));
            }
        }
        ConcurrencyTools.shutdown();
        return null;
    }

    private void printUsage() {
        myLogger.info(
                "\n\n\nThe available options for the TagsToSNPByAlignmentPlugin are as follows:\n"
                + "-i       Input .tbt file\n"
                + "-y       Use byte-formatted TBT file (*.tbt.byte)\n"
                + "-m       TagsOnPhysicalMap file containing genomic position of tags\n"
                + "-mUpd    Update TagsOnPhysicalMap file with allele calls, saved to specified file\n"
                + "-o       Output directory (default current directory)\n"
                + "-mxSites Maximum number of sites (SNPs) output per chromosome (default: " + maxSize + ")\n"
                + "-mnF     Minimum F (inbreeding coefficient) (default: " + minF + "  = no filter)\n"
                + "-p       Pedigree file containing full sample names (or expected names after merging) & expected inbreeding\n"
                + "         coefficient (F) for each.  Only taxa with expected F >= mnF used to calculate F = 1-Ho/He.\n"
                + "         (default: use ALL taxa to calculate F)\n"
                + "-mnMAF   Minimum minor allele frequency (default: " + minMAF + ")\n"
                + "-mnMAC   Minimum minor allele count (default: " + minMAC + ")\n"
                + "-mnLCov  Minimum locus coverage (proportion of Taxa with a genotype) (default: " + defaultMinPropTaxaWithLocus + ")\n"
                + "-ref     Path to reference genome in fasta format. DEFAULT: Don't use reference\n"
                + "         genome (instead, align tags with identical starting postions against each other)\n"
                + "-LocusBorder  All tags on either strand with start postions that differ by less than the specified\n"
                + "              integer (LocusBorder) are aligned to the reference genome to call SNPs at a locus.\n"
                + "              By default (without the -LocusBorder option), only tags with identical start postions and\n"
                + "              strand are grouped as a locus.\n"
                + "              Use of the -LocusBorder option requires that the -ref option is also invoked.\n"
                + "-inclRare  Include the rare alleles at site (3 or 4th states) (default: " + inclRare + ")\n"
                + "-inclGaps  Include sites where major or minor allele is a GAP (default: " + inclGaps + ")\n"
                + "-callBiSNPsWGap  Include sites where the third allele is a GAP (default: " + callBiallelicSNPsWithGap + ")\n"
                + "-s       Start chromosome\n"
                + "-e       End chromosome\n\n\n");
    }

    @Override
    public void setParameters(String[] args) {
        if (args.length == 0) {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }

        if (myArgsEngine == null) {
            myArgsEngine = new ArgsEngine();
            myArgsEngine.add("-i", "--input-file", true);
            myArgsEngine.add("-y", "--useTBTByte", false);
            myArgsEngine.add("-m", "--physical-map", true);
            myArgsEngine.add("-mUpd", "--update-physical-map", true);
            myArgsEngine.add("-o", "--output-directory", true);
            myArgsEngine.add("-mxSites", "--max-sites-per-chr", true);
            myArgsEngine.add("-mnF", "--minFInbreeding", true);
            myArgsEngine.add("-p", "--pedigree-file", true);
            myArgsEngine.add("-mnMAF", "--minMinorAlleleFreq", true);
            myArgsEngine.add("-mnMAC", "--minMinorAlleleCount", true);
            myArgsEngine.add("-mnLCov", "--minLocusCov", true);
            myArgsEngine.add("-ref", "--referenceGenome", true);
            myArgsEngine.add("-LocusBorder", "--locus-border", true);
            myArgsEngine.add("-inclRare", "--includeRare", false);
            myArgsEngine.add("-inclGaps", "--includeGaps", false);
            myArgsEngine.add("-callBiSNPsWGap", "--callBiSNPsWGap", false);
            myArgsEngine.add("-s", "--start-chromosome", true);
            myArgsEngine.add("-e", "--end-chromosome", true);
        }
        myArgsEngine.parse(args);

        if (myArgsEngine.getBoolean("-y")) {
            useTBTByte = true;
        }
        if (myArgsEngine.getBoolean("-i")) {
            String inputFileName = myArgsEngine.getString("-i");
            inputFile = new File(inputFileName);
            if (!inputFile.exists() || !inputFile.isFile()) {
                printUsage();
                throw new IllegalArgumentException("Can't find the TagsByTaxa input file (-i option: " + myArgsEngine.getString("-i") + ").");
            }
            outputFilePrefix = inputFile.getParentFile().getName();
            if (inputFileName.endsWith(".hdf") || inputFileName.endsWith(".h5")) {
                theTBT = new TagsByTaxaByteHDF5TagGroups(inputFileName);
            } else if (useTBTByte) {
                theTBT = new TagsByTaxaByteFileMap(inputFileName);
            } else {
                theTBT = new TagsByTaxaBitFileMap(inputFileName);
            }
        } else {
            printUsage();
            throw new IllegalArgumentException("Please specify a TagsByTaxa input file (-i option).");
        }
        if (myArgsEngine.getBoolean("-m")) {
            inTOPMFile = myArgsEngine.getString("-m");
            File inTOPMFileTest = new File(inTOPMFile);
            if (!inTOPMFileTest.exists() || !inTOPMFileTest.isFile()) {
                printUsage();
                throw new IllegalArgumentException("Can't find the TOPM input file (-m option: " + inTOPMFile + ").");
            }
            inTOPMFileTest = null;
            boolean loadBinary = (inTOPMFile.endsWith(".txt")) ? false : true;
            theTOPM = new TagsOnPhysicalMap(inTOPMFile, loadBinary);
        } else {
            printUsage();
            throw new IllegalArgumentException("Please specify a physical map file.");
        }
        if (myArgsEngine.getBoolean("-mUpd")) {
            this.isUpdateTOPM = true;
            this.outTOPMFile = myArgsEngine.getString("-mUpd");
        }

        //Set output directory and use the input TagsByTaxa filename for the output file prefix
        if (myArgsEngine.getBoolean("-o")) {
            outHapMap = myArgsEngine.getString("-o") + File.separator + outputFilePrefix;
        } else {
            outHapMap = inputFile.getParent() + File.separator + outputFilePrefix;
        }
        if (myArgsEngine.getBoolean("-mxSites")) {
            maxSize = Integer.parseInt(myArgsEngine.getString("-mxSites"));
        }
        if (myArgsEngine.getBoolean("-mnF")) {
            minF = Double.parseDouble(myArgsEngine.getString("-mnF"));
        }
        if (myArgsEngine.getBoolean("-p")) {
            String pedigreeFileStr = myArgsEngine.getString("-p");
            File pedigreeFile = new File(pedigreeFileStr);
            if (!pedigreeFile.exists() || !pedigreeFile.isFile()) {
                printUsage();
                throw new IllegalArgumentException("Can't find the pedigree input file (-p option: " + pedigreeFileStr + ").");
            }
            taxaFs = readTaxaFsFromFile(pedigreeFile);
            if (taxaFs == null) {
                throw new IllegalArgumentException("Problem reading the pedigree file. Progam aborted.");
            }
            if (!maskNonInbredTaxa()) {
                throw new IllegalArgumentException("Mismatch between taxa names in the pedigree file and TBT. Progam aborted.");
            }
            usePedigree = true;
        }
        if (myArgsEngine.getBoolean("-mnMAF")) {
            minMAF = Double.parseDouble(myArgsEngine.getString("-mnMAF"));
        }
        if (myArgsEngine.getBoolean("-mnMAC")) {
            minMAC = Integer.parseInt(myArgsEngine.getString("-mnMAC"));
        }
        minTaxaWithLocus = (int) Math.round(theTBT.getTaxaCount() * defaultMinPropTaxaWithLocus);
        if (myArgsEngine.getBoolean("-mnLCov")) {
            double minPropTaxaWithLocus = Double.parseDouble(myArgsEngine.getString("-mnLCov"));
            minTaxaWithLocus = (int) Math.round(theTBT.getTaxaCount() * minPropTaxaWithLocus);
        }
        if (myArgsEngine.getBoolean("-ref")) {
            refGenomeFileStr = myArgsEngine.getString("-ref");
            File refGenomeFile = new File(refGenomeFileStr);
            if (!refGenomeFile.exists() || !refGenomeFile.isFile()) {
                printUsage();
                throw new IllegalArgumentException("Can't find the reference genome fasta file (-ref option: " + refGenomeFileStr + ").");
            }
            includeReferenceGenome = true;
            refGenomeFile = null;
            System.gc();
        }
        if (myArgsEngine.getBoolean("-LocusBorder")) {
            if (!includeReferenceGenome) {
                printUsage();
                throw new IllegalArgumentException("The -LocusBorder option requires that the -ref option (referenceGenome) is also invoked.");
            }
            locusBorder = Integer.parseInt(myArgsEngine.getString("-LocusBorder"));
            fuzzyStartPositions = true;
        }
        if (myArgsEngine.getBoolean("-inclRare")) {
            inclRare = true;
        }
        if (myArgsEngine.getBoolean("-inclGaps")) {
            inclGaps = true;
        }
        if (myArgsEngine.getBoolean("-callBiSNPsWGap")) {
            callBiallelicSNPsWithGap = true;
        }
        if (myArgsEngine.getBoolean("-s")) {
            startChr = Integer.parseInt(myArgsEngine.getString("-s"));
        } else {
            printUsage();
            throw new IllegalArgumentException("Please specify start and end chromosome numbers.");
        }
        if (myArgsEngine.getBoolean("-e")) {
            endChr = Integer.parseInt(myArgsEngine.getString("-e"));
        } else {
            printUsage();
            throw new IllegalArgumentException("Please specify start and end chromosome numbers.");
        }
        if (endChr - startChr < 0) {
            printUsage();
            throw new IllegalArgumentException("Error: The start chromosome is larger than the end chromosome.");
        }
    }

    @Override
    public ImageIcon getIcon() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String getButtonName() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String getToolTipText() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public void runTagsToSNPByAlignment(MutableNucleotideAlignment theMSA, String outHapMap, int targetChromo, boolean requireGeneticSupport) {
        myLogger.info(String.format("minTaxaWithLocus:%d MinF:%g MinMAF:%g MinMAC:%d %n", minTaxaWithLocus, minF, minMAF, minMAC));
        myLogger.info(String.format("includeRare:%s includeGaps:%s %n", inclRare, inclGaps));
        long time = System.currentTimeMillis();
        TagsAtLocus currTAL = new TagsAtLocus(Integer.MIN_VALUE, Byte.MIN_VALUE, Integer.MIN_VALUE, includeReferenceGenome);
        int[] currPos = null;
        int countLoci = 0;
        for (int i = 0; (i < theTOPM.getSize()) && (theMSA.getSiteCount() < (maxSize - 1000)); i++) {
            int ri = theTOPM.getReadIndexForPositionIndex(i);  // process tags in order of physical position
            int[] newPos = theTOPM.getPositionArray(ri);
            if (newPos[chr] != targetChromo) {
                continue;    //Skip tags from other chromosomes
            }
            if (requireGeneticSupport && (theTOPM.getMapP(ri) < 2)) {
                continue; //Skip tags with low mapP scores
            }
            if ((fuzzyStartPositions && nearbyTag(newPos, currPos)) || Arrays.equals(newPos, currPos)) {
                currTAL.addTag(ri, theTOPM, theTBT, includeReferenceGenome);
            } else {
                if ((currTAL.getSize() > 1) && (currTAL.getNumberTaxaCovered() > minTaxaWithLocus)) {  // finish the current TAL
                    addSitesToMutableAlignment(currTAL, theMSA);  // note that with fuzzyStartPositions there may be no overlapping tags!!
                    countLoci++;
                    if (theMSA.getSiteCount() % 100 == 0) {
                        double rate = (double) theMSA.getSiteCount() / (double) (System.currentTimeMillis() - time);
                        myLogger.info(String.format(
                                "Chr:%d Pos:%d Loci=%d SNPs=%d rate=%g SNP/millisec %n", currPos[chr], currPos[startPosit], countLoci, theMSA.getSiteCount(), rate));
                    }
                }
                currPos = newPos; // start a new TAL with the current tag
                if ((currPos[str] != TagsOnPhysicalMap.byteMissing) && (currPos[startPosit] != TagsOnPhysicalMap.intMissing)) {  // we already know that currPos[chr]==targetChromo
                    currTAL = new TagsAtLocus(currPos[chr], (byte) currPos[str], currPos[startPosit], includeReferenceGenome);
                    currTAL.addTag(ri, theTOPM, theTBT, includeReferenceGenome);
                } else {
                    currPos = null;  // invalid position
                }
            }
        }
        if ((currTAL.getSize() > 1) && (currTAL.getNumberTaxaCovered() > minTaxaWithLocus)) { // then finish the final TAL for the targetChromo
            addSitesToMutableAlignment(currTAL, theMSA);
        }
        if (theMSA.getSiteCount() > 0) {
            theMSA.clean();
            //theMSA.sortSiteByPhysicalPosition();
            ExportUtils.writeToHapmap(theMSA, false, outHapMap, '\t', null);
            //ExportUtils.writeToHapmap(theMSA, false, outHapMap, '\t');
        }
        myLogger.info("Number of marker sites recorded for chr" + targetChromo + ": " + theMSA.getSiteCount());
    }

    /**
     * Fills an array of Locus objects with one locus for each supplied chromosome.  Creates a MutableNucleotideAlignment using
     * the locus list and TBT profile as input.  Returns the MSA object.
     */
    private static MutableNucleotideAlignment createMutableAlignment(TagsByTaxa theTBT, int startChr, int endChr, int maxSites) {
        Locus[] theL = new Locus[endChr - startChr + 1];
        for (int i = 0; i < theL.length; i++) {
            theL[i] = new Locus("" + (startChr + i), "" + (startChr + i), -1, -1, null, null);
        }
        IdGroup taxa = new SimpleIdGroup(theTBT.getTaxaNames());
        MutableNucleotideAlignment theMSA = MutableNucleotideAlignment.getInstance(taxa, maxSites, taxa.getIdCount(), maxSites);
        //MutableNucleotideAlignment theMSA = new MutableNucleotideAlignment(theTBT.getTaxaNames(), maxSites, theL);
        return theMSA;
    }

    boolean nearbyTag(int[] newTagPos, int[] currTagPos) {
        if (newTagPos == null || currTagPos == null) {
            return false;
        }
        // because we move through the TOPM in positional order, the newTag startPosition is guaranteed to be >= that of the current tag
        if (newTagPos[chr] == currTagPos[chr] && newTagPos[startPosit] - currTagPos[startPosit] < locusBorder) {  // &&newTagPos[str]==currTagPos[str]
            // grab all of the tags that align to a local region (until a gap > tolerance is reached)
            currTagPos[startPosit] = newTagPos[startPosit];
            return true;
        }
        return false;
    }

    private synchronized void addSitesToMutableAlignment(TagsAtLocus theTAL, MutableNucleotideAlignment theMSA) {
        if (theTAL.getSize() < 2) {
            return;  // need at least two (overlapping!) sequences to make an alignment
        }
        byte[][] callsBySite;
        if (includeReferenceGenome) {
            String refSeqInRegion = getRefSeqInRegion(theTAL);
            callsBySite = theTAL.getSNPCallsQuant(refSeqInRegion, callBiallelicSNPsWithGap);
        } else {
            callsBySite = theTAL.getSNPCallsQuant(callBiallelicSNPsWithGap);
        }
        if (callsBySite == null) {
            return;
        }
        int[] positionsInLocus = theTAL.getPositionsOfVariableSites();
        int strand = theTAL.getStrand();
        for (int s = 0; s < callsBySite.length; s++) {
            byte[] calls = callsBySite[s];
            byte[] alleles = null;
            if ((alleles = isSiteGood(calls)) == null) {
                continue;
            }
            int currSite = theMSA.getSiteCount();
            //theMSA.setLocusOfSite(currSite, "" + theTAL.getChromosome());
            String chromosome = String.valueOf(theTAL.getChromosome());
            theMSA.setLocusOfSite(currSite, new Locus(chromosome, chromosome, -1, -1, null, null));
            int position = (strand == -1) ? theTAL.getMinStartPosition() - positionsInLocus[s] : theTAL.getMinStartPosition() + positionsInLocus[s];
            theMSA.setPositionOfSite(currSite, position);
            //theMSA.setStrandOfSite(currSite, (byte) '+');  // all minus strand genotypes will be complemented to plus strand
            for (int tx = 0; tx < theTBT.getTaxaCount(); tx++) {
                if (strand == -1) {
                    theMSA.setBase(tx, currSite, complement(calls[tx]));  // complement to plus strand
                } else {
                    theMSA.setBase(tx, currSite, calls[tx]);
                }
            }
            if (this.isUpdateTOPM) {
                updateTOPM(theTAL, s, position, strand, alleles);
            }
            if (currSite % 100 == 0) {
                System.out.printf("Site:%d Position:%d %n", currSite, position);
            }
        }
    }

    private void updateTOPM(TagsAtLocus myTAL, int variableSite, int position, int strand, byte[] alleles) {
        for (int tg = 0; tg < myTAL.getSize(); tg++) {
            byte baseToAdd = myTAL.getCallAtVariableSiteForTag(variableSite, tg);
            if (baseToAdd == 0 || baseToAdd == Alignment.UNKNOWN_DIPLOID_ALLELE) {
                continue;
            }
            boolean matched = false;
            for (byte cb : alleles) {
                if (baseToAdd == cb) {
                    matched = true;
                    break;
                }
            }
            // so that all tags in the tagAlignment have the same corresponding variants in the TOPM, add a variant no matter what (set to N if needed)
            int topmTagIndex = myTAL.getTOPMIndexOfTag(tg);
            byte offset = (byte) (position - myTAL.getMinStartPosition());
            if (!matched) {
                baseToAdd = Alignment.UNKNOWN_DIPLOID_ALLELE;
            }
            if (strand == -1) {
                baseToAdd = complement(baseToAdd);  // record everything relative to the plus strand
            }
            theTOPM.addVariant(topmTagIndex, offset, baseToAdd);
        }
    }

    /**
     *
     * @param calls
     * @return
     */
    private byte[] isSiteGood(byte[] calls) {
        int[][] alleles = AlignmentUtils.getAllelesSortedByFrequency(calls);
        //int[][] alleles = getSortedAlleleCounts(calls);
        int aCnt = alleles[1][0] + alleles[1][1];
        double theMAF = (double) alleles[1][1] / (double) aCnt;
        if ((theMAF < minMAF) && (alleles[1][1] < minMAC)) {
            return null;  // note that a site only needs to pass one of the criteria, minMAF &/or minMAC
        }
        byte majAllele = (byte) alleles[0][0];
        byte minAllele = (byte) alleles[0][1];
        if (!inclGaps && ((majAllele == NucleotideAlignmentConstants.GAP_DIPLOID_ALLELE) || (minAllele == NucleotideAlignmentConstants.GAP_DIPLOID_ALLELE))) {
            return null;
        }
        //byte hetG = IUPACNucleotides.getDegerateSNPByteFromTwoSNPs(majAllele, minAllele);
        byte hetG = AlignmentUtils.getDiploidValue(majAllele, minAllele);
        double obsF = calculateF(calls, alleles, hetG, theMAF);
        if (obsF < minF) {
            return null;
        }
        if (!inclRare) {
            if (callBiallelicSNPsWithGap) {
                for (int i = 0; i < calls.length; i++) {
                    if ((calls[i] != majAllele) && (calls[i] != minAllele) && (calls[i] != hetG) && (calls[i] != NucleotideAlignmentConstants.GAP_DIPLOID_ALLELE) && (calls[i] != '0')) {
                        calls[i] = Alignment.UNKNOWN_DIPLOID_ALLELE;
                    }
                }
            } else {
                for (int i = 0; i < calls.length; i++) {
                    if ((calls[i] != majAllele) && (calls[i] != minAllele) && (calls[i] != hetG)) {
                        calls[i] = Alignment.UNKNOWN_DIPLOID_ALLELE;
                    }
                }
            }
        }
        byte[] majMinAlleles = {majAllele, minAllele};
        return majMinAlleles;
    }

    /*
    private int[][] getSortedAlleleCounts(byte[] calls) {
        byte[] nuc = {'A', 'C', 'G', 'T', '-', '+', 'N'};
        int[] nucIndex = new int[Byte.MAX_VALUE];
        for (int i = 0; i < nuc.length; i++) {
            nucIndex[nuc[i]] = i;
        }
        int[] cnt = new int[nuc.length];
        for (byte dc : calls) {
            byte[] cc = IUPACNucleotides.getDiploidValueFromIUPACCode(dc);
            cnt[nucIndex[cc[0]]]++;
            cnt[nucIndex[cc[1]]]++;
        }
        int[][] alleles = new int[2][nuc.length - 1];  // alleles[0]=allele; alleles[1]=count
        for (int i = 0; i < nuc.length - 1; i++) {  // "i<nuc.length-1" stops N from being included
            alleles[0][i] = nuc[i];
            alleles[1][i] = cnt[i];
        }
        boolean change = true;  // sort the alleles by descending frequency
        while (change) {
            change = false;
            for (int k = 0; k < nuc.length - 2; k++) {
                if (alleles[1][k] < alleles[1][k + 1]) {
                    int temp = alleles[0][k];
                    alleles[0][k] = alleles[0][k + 1];
                    alleles[0][k + 1] = temp;
                    int tempCount = alleles[1][k];
                    alleles[1][k] = alleles[1][k + 1];
                    alleles[1][k + 1] = tempCount;
                    change = true;
                }
            }
        }
        return alleles;
    }
     */

    private double calculateF(byte[] calls, int[][] alleles, byte hetG, double theMAF) {
        boolean report = false;
        double obsF;
        int hetGCnt = 0;
        if (usePedigree) {
            byte[] callsToUse = filterCallsForInbreds(calls);
            //int[][] allelesToUse = getSortedAlleleCounts(callsToUse);
            int[][] allelesToUse = AlignmentUtils.getAllelesSortedByFrequency(callsToUse);
            int aCnt = allelesToUse[1][0] + allelesToUse[1][1];
            double newMAF = (double) allelesToUse[1][1] / (double) aCnt;
            if (newMAF <= 0.0) {
                return 1.0;  // lack of variation in the known inbreds will NOT reject a SNP
            }
            byte majAllele = (byte) allelesToUse[0][0];
            byte minAllele = (byte) allelesToUse[0][1];
            //byte newHetG = IUPACNucleotides.getDegerateSNPByteFromTwoSNPs(majAllele, minAllele);
            byte newHetG = AlignmentUtils.getDiploidValue(majAllele, minAllele);
            for (byte i : callsToUse) {
                if (i == newHetG) {
                    hetGCnt++;
                }
            }
            int majGCnt = (allelesToUse[1][0] - hetGCnt) / 2; // number of homozygous major allele genotypes
            int minGCnt = (allelesToUse[1][1] - hetGCnt) / 2; // number of homozygous minor allele genotypes
            double propHets = (double) hetGCnt / (double) (hetGCnt + majGCnt + minGCnt);
            double expHets = 2.0 * newMAF * (1 - newMAF);
            obsF = 1.0 - (propHets / expHets);
            if (report) {
                System.out.printf("%d %d %d propHets:%g expHets:%g obsF:%g %n", majGCnt, minGCnt, hetGCnt, propHets, expHets, obsF);
            }
            return obsF;
        } else {
            for (byte i : calls) {
                if (i == hetG) {
                    hetGCnt++;
                }
            }
            int majGCnt = (alleles[1][0] - hetGCnt) / 2; // number of homozygous major allele genotypes
            int minGCnt = (alleles[1][1] - hetGCnt) / 2; // number of homozygous minor allele genotypes
            double propHets = (double) hetGCnt / (double) (hetGCnt + majGCnt + minGCnt);
            double expHets = 2.0 * theMAF * (1 - theMAF);
            obsF = 1.0 - (propHets / expHets);
            if (report) {
                System.out.printf("%d %d %d propHets:%g expHets:%g obsF:%g %n", majGCnt, minGCnt, hetGCnt, propHets, expHets, obsF);
            }
            return obsF;
        }
    }

    private byte[] filterCallsForInbreds(byte[] calls) {
        byte[] callsForInbredsOnly = new byte[nInbredTaxa];
        int inbred = 0;
        for (int taxon = 0; taxon < calls.length; taxon++) {
            if (useTaxaForMinF[taxon]) {
                callsForInbredsOnly[inbred] = calls[taxon];
                inbred++;
            }
        }
        return callsForInbredsOnly;
    }

    public static HashMap<String, Double> readTaxaFsFromFile(File pedigreeFile) {
        HashMap<String, Double> taxaFs = new HashMap<String, Double>();
        String inputLine = "Nothing has been read from the pedigree input file yet";
        int nameCol = -1, fCol = -1, nTaxa = 0;
        try {
            BufferedReader br = new BufferedReader(new FileReader(pedigreeFile), 65536);
            inputLine = br.readLine();  // header line
            String[] cells = inputLine.split("\t");  // headers
            for (int col = 0; col < cells.length; col++) {
                if (cells[col].equalsIgnoreCase("Name")) {
                    nameCol = col;
                }
                if (cells[col].equalsIgnoreCase("F")) {
                    fCol = col;
                }
            }
            if (nameCol > -1 && fCol > -1) {
                while ((inputLine = br.readLine()) != null) {
                    cells = inputLine.split("\t");
                    if (cells[fCol].equals("NA")) {
                        taxaFs.put(cells[nameCol], -2.0);
                    } else {
                        taxaFs.put(cells[nameCol], Double.parseDouble(cells[fCol]));
                    }
                    ++nTaxa;
                }
            } else {
                throw new Exception("Name and/or F column not found in header");
            }
        } catch (Exception e) {
            myLogger.error("Catch in reading pedigree file e=" + e);
            e.printStackTrace();
            System.out.println(inputLine);
            return null;
        }
        myLogger.info(nTaxa + " taxa read from the pedigree file");
        return taxaFs;
    }

    private boolean maskNonInbredTaxa() {
        useTaxaForMinF = new boolean[theTBT.getTaxaCount()];  // initialized to false
        nInbredTaxa = 0;
        try {
            for (int taxon = 0; taxon < theTBT.getTaxaCount(); taxon++) {
                if (taxaFs.containsKey(theTBT.getTaxaName(taxon))) {
                    if (taxaFs.get(theTBT.getTaxaName(taxon)) >= minF) {
                        useTaxaForMinF[taxon] = true;
                        nInbredTaxa++;
                    }
                } else {
                    throw new Exception("Taxon " + theTBT.getTaxaName(taxon) + " not found in the pedigree file");
                }
            }
            myLogger.info(nInbredTaxa + " taxa with an Expected F >= the mnF of " + minF + " were found in the input TBT");
            return true;
        } catch (Exception e) {
            myLogger.error("Mismatch between TBT and pedigree file e=" + e);
            e.printStackTrace();
            return false;
        }
    }

    private long[] readReferenceGenomeChr(String inFileStr, int targetChr) {
        int nBases = getLengthOfReferenceGenomeChr(inFileStr, targetChr);
        int basesPerLong = BaseEncoder.chunkSize;
        int nLongs = (nBases % basesPerLong == 0) ? nBases / basesPerLong : (nBases / basesPerLong) + 1;
        long[] refGenomeChrAsLongs = new long[nLongs];
        myLogger.info("\n\nReading in the target chromosome " + targetChr + " from the reference genome fasta file: " + inFileStr);
        String temp = "Nothing has been read yet from the reference genome fasta file";
        try {
            BufferedReader br = new BufferedReader(new FileReader(new File(inFileStr)));
            StringBuilder currStrB = new StringBuilder();
            int currChr = Integer.MIN_VALUE, chunk = 0;
            while (br.ready()) {
                temp = br.readLine().trim();
                if (temp.startsWith(">")) {
                    if (chunk > 0) {
                        break;  // finished reading the targetChr (no need to read the rest of the file)
                    }
                    String chrS = temp.replace(">", "");
                    chrS = chrS.replace("chr", "");
                    currChr = Integer.parseInt(chrS);  // don't need to catch exception because getLengthOfReferenceGenomeChr() would have caught it already
                    myLogger.info("Currently reading chromosome " + currChr + " (target chromosome = " + targetChr + ")");
                } else if (currChr == targetChr) {
                    currStrB.append(temp.replace("N", "A")); // BaseEncoder encodes sequences with N's as (long) -1
                    while (currStrB.length() >= basesPerLong) {
                        refGenomeChrAsLongs[chunk] = BaseEncoder.getLongFromSeq(currStrB.substring(0, basesPerLong));
                        currStrB = (currStrB.length() > basesPerLong) ? new StringBuilder(currStrB.substring(basesPerLong)) : new StringBuilder();
                        chunk++;
                        if (chunk % 1000000 == 0) {
                            myLogger.info(chunk + " chunks of " + basesPerLong + " bases read from the reference genome fasta file for chromosome " + targetChr);
                        }
                    }
                }
            }
            if (currStrB.length() > 0) {
                refGenomeChrAsLongs[chunk] = BaseEncoder.getLongFromSeq(currStrB.toString());
                chunk++;
            }
            myLogger.info("\n\nFinished reading target chromosome " + targetChr + " into a total of " + chunk + " " + basesPerLong + "bp chunks\n\n");
            if (chunk != nLongs) {
                throw new Exception("The number of 32 base chunks read (" + chunk + ") was not equal to the expected number (" + nLongs + ")");
            }
            br.close();
        } catch (Exception e) {
            myLogger.error("Exception caught while reading the reference genome fasta file at line.  Error=" + e);
            e.printStackTrace();
            System.exit(1);
        }
        return refGenomeChrAsLongs;
    }

    private int getLengthOfReferenceGenomeChr(String inFileStr, int targetChr) {
        myLogger.info("\n\nDetermining the length (in bases) of target chromosome " + targetChr + " in the reference genome fasta file: " + inFileStr);
        String temp = "Nothing has been read yet from the reference genome fasta file";
        int line = 0, nBases = 0;
        try {
            BufferedReader br = new BufferedReader(new FileReader(new File(inFileStr)));
            int currChr = Integer.MIN_VALUE;
            while (br.ready()) {
                temp = br.readLine().trim();
                line++;
                if (line % 1000000 == 0) {
                    myLogger.info(line + " lines read from the reference genome fasta file");
                }
                if (temp.startsWith(">")) {
                    if (nBases > 0) {
                        break;  // finished reading the targetChr (no need to read the rest of the file)
                    }
                    String chrS = temp.replace(">", "");
                    chrS = chrS.replace("chr", "");
                    try {
                        currChr = Integer.parseInt(chrS);
                    } catch (NumberFormatException e) {
                        myLogger.error("\n\nTagsToSNPByAlignment detected a non-numeric chromosome name in the reference genome sequence fasta file: " + chrS
                                + "\n\nPlease change the FASTA headers in your reference genome sequence to integers "
                                + "(>1, >2, >3, etc.) OR to 'chr' followed by an integer (>chr1, >chr2, >chr3, etc.)\n\n");
                        System.exit(1);
                    }
                    myLogger.info("Currently reading chromosome " + currChr + " (target chromosome = " + targetChr + ")");
                } else if (currChr == targetChr) {
                    nBases += temp.length();
                }
            }
            myLogger.info("The target chromosome " + targetChr + " is " + nBases + " bases long");
            br.close();
        } catch (Exception e) {
            myLogger.error("Exception caught while reading the reference genome fasta file at line " + line + "  e=" + e);
            e.printStackTrace();
            System.exit(1);
        }
        return nBases;
    }

    private String getRefSeqInRegion(TagsAtLocus theTAL) {
        int basesPerLong = BaseEncoder.chunkSize;
        int refSeqStartPosition = theTAL.getMinStartPosition() - 128;
        int startIndex = Math.max((refSeqStartPosition / basesPerLong) - 1, 0);
        int refSeqEndPosition = theTAL.getMaxStartPosition() + 128;
        int endIndex = Math.min((refSeqEndPosition / basesPerLong) + 1, refGenomeChr.length - 1);
        StringBuilder sb = new StringBuilder();
        for (int i = startIndex; i <= endIndex; ++i) {
            sb.append(BaseEncoder.getSequenceFromLong(refGenomeChr[i]));
        }
        theTAL.setMinStartPosition(startIndex * basesPerLong + 1);
        return sb.toString();
    }

    public static byte complement(byte geno) {
        byte comp = Byte.MIN_VALUE;
        switch (geno) {
            case 'A':
                comp = 'T';
                break;
            case 'C':
                comp = 'G';
                break;
            case 'G':
                comp = 'C';
                break;
            case 'T':
                comp = 'A';
                break;
            case 'K':
                comp = 'M';
                break;
            case 'M':
                comp = 'K';
                break;
            case 'R':
                comp = 'Y';
                break;
            case 'S':
                comp = 'S';
                break;
            case 'W':
                comp = 'W';
                break;
            case 'Y':
                comp = 'R';
                break;
            case '-':
                comp = '-';
                break;  // both strands have the deletion
            case '+':
                comp = '+';
                break;  // both strands have the insertion
            case '0':
                comp = '0';
                break;
            case 'N':
                comp = 'N';
                break;
            default:
                comp = 'N';
                break;
        }
        return comp;
    }
}
