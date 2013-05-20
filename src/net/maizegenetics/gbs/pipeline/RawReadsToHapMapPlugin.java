/*
 * RawReadsToHapMapPlugin
 */
package net.maizegenetics.gbs.pipeline;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import net.maizegenetics.util.MultiMemberGZIPInputStream;
import net.maizegenetics.gbs.homology.ParseBarcodeRead;
import net.maizegenetics.gbs.homology.ReadBarcodeResult;
import net.maizegenetics.util.ArgsEngine;
import net.maizegenetics.util.DirectoryCrawler;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.gbs.homology.TagMatchFinder;
import net.maizegenetics.pal.alignment.MutableNucleotideAlignment;
import net.maizegenetics.pal.alignment.ExportUtils;
import net.maizegenetics.pal.alignment.Locus;
import java.awt.Frame;
import java.util.ArrayList;
import java.util.TreeMap;
import java.util.TreeSet;
import javax.swing.ImageIcon;
import net.maizegenetics.gbs.maps.TOPMInterface;
import net.maizegenetics.gbs.maps.TOPMUtils;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.AlignmentUtils;
import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.SimpleIdGroup;
import org.apache.log4j.Logger;

/**
 * This pipeline converts all of the fastq AND/OR qseq files in the input folder
 * to hapmap genotype files (nOutputHapMapFiles = nChrs X nRawSeqFiles). We
 * refer to this step as the "Production Pipeline".
 *
 * It requires a TOPM with variants added from a previous "Discovery Pipeline"
 * run.
 *
 * @author jcg233
 */
public class RawReadsToHapMapPlugin extends AbstractPlugin {

    private final Logger myLogger = Logger.getLogger(RawReadsToHapMapPlugin.class);
    private ArgsEngine myArgsEngine = null;
    private String[] myRawSeqFileNames = null;
    private String myKeyFile = null;
    private String myEnzyme = null;
    private String myOutputDir = null;
    private TOPMInterface topm = null;
    private int maxDivergence = 0;
    private int[] chromosomes = null;
    private boolean fastq = true;
    private IdGroup taxaIDGroup = null;

    public RawReadsToHapMapPlugin() {
        super(null, false);
    }

    public RawReadsToHapMapPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    private void printUsage() {
        myLogger.info(
                "\nThe options for the RawReadsToHapMapPlugin TASSEL plugin are as follows:\n"
                + "-i  Input directory containing fastq AND/OR qseq files\n"
                + "-k  Barcode key file\n"
                + "-e  Enzyme used to create the GBS library\n"
                + "-o  Output directory\n"
                + "-m  Physical map file containing alignments and variants (production TOPM)\n");
//                + "-d  Maximum divergence (edit distance) between new read and previously mapped read (Default: 0 = perfect matches only)\n");  // NOT IMPLEMENTED YET
    }

    @Override
    public void setParameters(String[] args) {
        if (args.length == 0) {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }
        if (myArgsEngine == null) {
            myArgsEngine = new ArgsEngine();
            myArgsEngine.add("-i", "--input-directory", true);
            myArgsEngine.add("-k", "--key-file", true);
            myArgsEngine.add("-e", "--enzyme", true);
            myArgsEngine.add("-o", "--output-directory", true);
            myArgsEngine.add("-m", "--physical-map", true);
            myArgsEngine.add("-d", "--divergence", true);
        }
        myArgsEngine.parse(args);
        String tempDirectory = myArgsEngine.getString("-i");
        if (tempDirectory != null) {
            File rawSeqDirectory = new File(tempDirectory);
            if (!rawSeqDirectory.isDirectory()) {
                printUsage();
                throw new IllegalArgumentException("setParameters: The input name you supplied is not a directory: " + tempDirectory);
            }
            myRawSeqFileNames = DirectoryCrawler.listFileNames(rawSeqFileNameRegex, rawSeqDirectory.getAbsolutePath());
            if (myRawSeqFileNames.length == 0 || myRawSeqFileNames == null) {
                printUsage();
                throw new IllegalArgumentException(noMatchingRawSeqFileNamesMessage + tempDirectory);
            } else {
                myLogger.info("RawReadsToHapMapPlugin: setParameters: Using the following GBS raw sequence data files:");
                for (String filename : myRawSeqFileNames) {
                    myLogger.info(filename);
                }
            }
        }
        if (myArgsEngine.getBoolean("-k")) {
            myKeyFile = myArgsEngine.getString("-k");
        } else {
            printUsage();
            throw new IllegalArgumentException("Please specify a key file (option -k).");
        }
        if (myArgsEngine.getBoolean("-e")) {
            myEnzyme = myArgsEngine.getString("-e");
        } else {
            printUsage();
            throw new IllegalArgumentException("Please specify the enzyme used to create the GBS library.");
        }
        if (myArgsEngine.getBoolean("-o")) {
            myOutputDir = myArgsEngine.getString("-o");
            File outDirectory = new File(myOutputDir);
            if (!outDirectory.isDirectory()) {
                printUsage();
                throw new IllegalArgumentException("The output name you supplied (option -o) is not a directory: " + myOutputDir);
            }
            outDirectory = null;
        } else {
            printUsage();
            throw new IllegalArgumentException("Please specify an output directory (option -o).");
        }
        if (myArgsEngine.getBoolean("-m")) {
            topm = TOPMUtils.readTOPM(myArgsEngine.getString("-m"));
        } else {
            printUsage();
            throw new IllegalArgumentException("Please specify a TagsOnPhysicalMap file (-m)");
        }
        if (myArgsEngine.getBoolean("-d")) {
            maxDivergence = Integer.parseInt(myArgsEngine.getString("-d"));
        }
    }

    @Override
    public DataSet performFunction(DataSet input) {
//        myLogger.addAppender(new ConsoleAppender(new SimpleLayout()));
        if ((myEnzyme == null) || (myEnzyme.length() == 0)) {
            printUsage();
            throw new IllegalStateException("performFunction: enzyme must be set.");
        }
        // TODO - More checks to validate parameters...
        translateRawReadsToHapmap();
        return null;
    }

    private void translateRawReadsToHapmap() {
        for (int laneNum = 0; laneNum < myRawSeqFileNames.length; laneNum++) {
            int[] counters = {0, 0, 0, 0, 0, 0}; // 0:allReads 1:goodBarcodedReads 2:goodMatched 3:perfectMatches 4:imperfectMatches 5:singleImperfectMatches
            ParseBarcodeRead thePBR = setUpBarcodes(laneNum);
            if (thePBR == null || thePBR.getBarCodeCount() == 0) {
                System.out.println("No barcodes found. Skipping this flowcell lane.");
                continue;
            }
            MutableNucleotideAlignment[] outMSA = setUpMutableNucleotideAlignments(thePBR);
            myLogger.info("Looking for known SNPs in sequence reads...");
            String temp = "";
            BufferedReader br = getBufferedReader(laneNum);
            try {
                while ((temp = br.readLine()) != null) {
                    if (counters[0] % 1000000 == 0) {
                        reportProgress(counters);
                    }
                    ReadBarcodeResult rr = readSequenceRead(br, temp, thePBR, counters);
                    if (rr != null) {
                        counters[1]++;  // goodBarcodedReads
                        int tagIndex = topm.getTagIndex(rr.getRead());
                        if (tagIndex >= 0) {
                            counters[3]++;  // perfectMatches
                        }
                        if (tagIndex < 0 && maxDivergence > 0) {
                            tagIndex = findBestImperfectMatch(rr.getRead(), counters);
                        }
                        if (tagIndex < 0) {
                            continue;
                        }
                        counters[2]++;  // goodMatched++;
                        recordVariantsFromTag(outMSA, tagIndex, taxaIDGroup.whichIdNumber(rr.getTaxonName()));
                    }
                }
                br.close();
            } catch (Exception e) {
                System.out.println("Catch in translateRawReadsToHapmap() at nReads=" + counters[0] + " e=" + e);
                System.out.println(temp);
                e.printStackTrace();
            }
            writeHapMapFiles(outMSA, laneNum, counters);
        }
    }

    private ParseBarcodeRead setUpBarcodes(int laneNum) {
        System.gc();
        System.out.println("\nWorking on GBS raw sequence file: " + myRawSeqFileNames[laneNum]);
        fastq = true;
        if (myRawSeqFileNames[laneNum].substring(myRawSeqFileNames[laneNum].lastIndexOf(File.separator)).contains("qseq")) {
            fastq = false;
        }
        if (fastq) {
            System.out.println("\tThis file is assumed to be in fastq format");
        } else {
            System.out.println("\tThis file contains 'qseq' in its name so is assumed to be in qseq format");
        }
        File rawSeqFile = new File(myRawSeqFileNames[laneNum]);
        String[] np = rawSeqFile.getName().split("_");

        ParseBarcodeRead thePBR = null;
        if (np.length == 3) {
            thePBR = new ParseBarcodeRead(myKeyFile, myEnzyme, np[0], np[1]);
        } else if (np.length == 4) {
            thePBR = new ParseBarcodeRead(myKeyFile, myEnzyme, np[0], np[2]);
        } else if (np.length == 5) {
            thePBR = new ParseBarcodeRead(myKeyFile, myEnzyme, np[1], np[3]);
        } else {
            printFileNameConventions(myRawSeqFileNames[laneNum]);
            return thePBR;
        }
        System.out.println("Total barcodes found in key file for this lane:" + thePBR.getBarCodeCount());
        return thePBR;
    }

    private MutableNucleotideAlignment[] setUpMutableNucleotideAlignments(ParseBarcodeRead thePBR) {
        myLogger.info("\nCounting sites in TOPM file.");
        ArrayList<int[]> uniquePositions = getUniquePositions();  // need this method in TOPMInterface
        myLogger.info("Creating alignment objects to hold the genotypic data (one per chromosome in the TOPM).");
        MutableNucleotideAlignment[] outMSA = new MutableNucleotideAlignment[chromosomes.length];
        for (int i = 0; i < outMSA.length; i++) {
            outMSA[i] = MutableNucleotideAlignment.getInstance(new SimpleIdGroup(thePBR.getTaxaNames()), uniquePositions.get(i).length);
        }
        taxaIDGroup = outMSA[0].getIdGroup();
        myLogger.info("Adding sites from the TOPM file to the alignment objects.");
        for (int i = 0; i < outMSA.length; i++) {
            int currSite = 0;
            for (int j = 0; j < uniquePositions.get(i).length; j++) {
                String chromosome = Integer.toString(chromosomes[i]);
                outMSA[i].addSite(currSite);
                outMSA[i].setLocusOfSite(currSite, new Locus(chromosome, chromosome, -1, -1, null, null));
                outMSA[i].setPositionOfSite(currSite, uniquePositions.get(i)[j]);
                currSite++;
            }
            outMSA[i].clean();
        }
        return outMSA;
    }

    private ArrayList<int[]> getUniquePositions() {
        ArrayList<int[]> uniquePositions = new ArrayList<int[]>();
        chromosomes = topm.getChromosomes();
        for (int i = 0; i < chromosomes.length; i++) {
            uniquePositions.add(topm.getUniquePositions(chromosomes[i]));
        }
        return uniquePositions;
    }

    private BufferedReader getBufferedReader(int laneNum) {
        BufferedReader br = null;
        try {
            if (myRawSeqFileNames[laneNum].endsWith(".gz")) {
                br = new BufferedReader(new InputStreamReader(new MultiMemberGZIPInputStream(new FileInputStream(myRawSeqFileNames[laneNum]))));
            } else {
                br = new BufferedReader(new FileReader(myRawSeqFileNames[laneNum]), 65536);
            }
        } catch (Exception e) {
            System.out.println("Catch in getBufferedReader(): e=" + e);
            e.printStackTrace();
        }
        return br;
    }

    private void reportProgress(int[] counters) {
        System.out.println(
                "totalReads:" + counters[0]
                + " goodBarcodedReads:" + counters[1]
                + " goodMatchedToTOPM:" + counters[2]
                + " perfectMatches:" + counters[3]
                + " nearMatches:" + counters[4]
                + " uniqueNearMatches:" + counters[5]);
    }

    private ReadBarcodeResult readSequenceRead(BufferedReader br, String temp, ParseBarcodeRead thePBR, int[] counters) {
        ReadBarcodeResult rr = null;
        String sl = "";
        try {
            if (fastq) {
                sl = br.readLine();    // read the 2nd line in the set of 4 lines = sequence
                temp = br.readLine();  // skip the 3rd line
                temp = br.readLine();  // skip the 4th in the set of 4 lines = quality score (note that the QS is scaled differently in Cassava 1.8 - we don't use it so it is not corrected here)
                rr = thePBR.parseReadIntoTagAndTaxa(sl, null, true, 0);
            } else {  // qseq
                String[] jj = temp.split("\\s");
                sl = jj[8];
                rr = thePBR.parseReadIntoTagAndTaxa(sl, null, false, 0);
            }
        } catch (Exception e) {
            System.out.println("Catch in readSequenceRead() at nReads=" + counters[0] + " e=" + e);
            System.out.println(temp);
            e.printStackTrace();
        }
        counters[0]++;  // allReads
        return rr;
    }

    private int findBestImperfectMatch(long[] read, int[] counters) {
        // this method is not ready for prime time -- to resolve a tie, it currently chooses a random tag out of the tied tags
        int tagIndex = -1;
        TagMatchFinder tmf = new TagMatchFinder(topm);
        TreeMap<Integer, Integer> bestHitsAndDiv = tmf.findMatchesWithIntLengthWords(read, maxDivergence, true);
        if (bestHitsAndDiv.size() > 0) {
            counters[4]++; // imperfectMatches
            if (bestHitsAndDiv.size() == 1) {
                counters[5]++; // singleImperfectMatches
            }
            tagIndex = bestHitsAndDiv.firstKey();  // a random tag (firstKey) chosen to resolve the tie = suboptimal behavior
        }
        return tagIndex;
    }

    private void recordVariantsFromTag(MutableNucleotideAlignment[] outMSA, int tagIndex, int taxonIndex) {
        int chromosome = topm.getChromosome(tagIndex);
        if (chromosome == Integer.MIN_VALUE) {
            return;
        }
        int chrIndex = 0;
        for (int i = 0; i < chromosomes.length; i++) {
            if (chromosomes[i] == chromosome) {
                chrIndex = i;
                break;
            }
        }
        Locus locus = topm.getLocus(tagIndex);
        int startPos = topm.getStartPosition(tagIndex);
        for (int variant = 0; variant < topm.getMaxNumVariants(); variant++) {
            byte newBase = topm.getVariantDef(tagIndex, variant); // Nb: this should return Tassel4 allele encodings
            if ((newBase == topm.getByteMissing()) || (newBase == Alignment.UNKNOWN_ALLELE)) {
                continue;
            }
            int offset = topm.getVariantPosOff(tagIndex, variant);
            int pos = startPos + offset;
            int currSite = outMSA[chrIndex].getSiteOfPhysicalPosition(pos, locus);
            if (currSite < 0) {
                continue;
            }
            byte newGeno = AlignmentUtils.getDiploidValue(newBase, newBase);
            byte prevGeno = outMSA[chrIndex].getBase(taxonIndex, currSite);
            if (prevGeno == Alignment.UNKNOWN_DIPLOID_ALLELE) {
                outMSA[chrIndex].setBase(taxonIndex, currSite, newGeno);
            } else if (newGeno != prevGeno) {
                outMSA[chrIndex].setBase(taxonIndex, currSite, resolveGenoFromCallPair(prevGeno, newBase));
            }
        }
    }

    private void writeHapMapFiles(MutableNucleotideAlignment[] outMSA, int laneNum, int[] counters) {
        for (int i = 0; i < outMSA.length; i++) {
            outMSA[i].clean();
            AlignmentFilterByGBSUtils.getCoverage_MAF_F_Dist(outMSA[i], false);
            String outFileS = myOutputDir + myRawSeqFileNames[laneNum].substring(myRawSeqFileNames[laneNum].lastIndexOf(File.separator));
            outFileS = outFileS.replaceAll(rawSeqFileNameReplaceRegex, "_c" + chromosomes[i]); // ".hmp.txt" gets added by ExportUtils.writeToHapmap
            ExportUtils.writeToHapmap(outMSA[i], false, outFileS, '\t', this);
        }
        System.out.println("Total number of reads in lane=" + counters[0]);
        System.out.println("Total number of good, barcoded reads=" + counters[1]);
        int filesDone = laneNum + 1;
        System.out.println("Finished reading " + filesDone + " of " + myRawSeqFileNames.length + " sequence files: " + myRawSeqFileNames[laneNum] + "\n");
    }

    private void printFileNameConventions(String actualFileName) {
        System.out.println("Error in parsing file name:");
        System.out.println("   The raw sequence filename does not contain either 3, 4, or 5 underscore-delimited values.");
        System.out.println("   Acceptable file naming conventions include the following (where FLOWCELL indicates the flowcell name and LANE is an integer):");
        System.out.println("       FLOWCELL_LANE_fastq.gz");
        System.out.println("       FLOWCELL_s_LANE_fastq.gz");
        System.out.println("       code_FLOWCELL_s_LANE_fastq.gz");
        System.out.println("       FLOWCELL_LANE_fastq.txt.gz");
        System.out.println("       FLOWCELL_s_LANE_fastq.txt.gz");
        System.out.println("       code_FLOWCELL_s_LANE_fastq.txt.gz");
        System.out.println("       FLOWCELL_LANE_qseq.txt.gz");
        System.out.println("       FLOWCELL_s_LANE_qseq.txt.gz");
        System.out.println("       code_FLOWCELL_s_LANE_qseq.txt.gz");
        System.out.println("");
        System.out.println("   Actual Filename: " + actualFileName);
    }
    private String rawSeqFileNameRegex =
            "(?i)" + // case insensitve
            ".*\\.fq" + "$|"
            + ".*\\.fq\\.gz" + "$|"
            + ".*\\.fastq" + "$|"
            + ".*_fastq\\.txt" + "$|"
            + ".*_fastq\\.gz" + "$|"
            + ".*_fastq\\.txt\\.gz" + "$|"
            + ".*_sequence\\.txt" + "$|"
            + ".*_sequence\\.txt\\.gz" + "$|"
            + ".*_qseq\\.txt" + "$|"
            + ".*_qseq\\.txt\\.gz" + "$";
    //            \\. denotes escape . so it doesn't mean 'any char'
    // NOTE: If you add addtional file naming conventions here, you must also
    //       add them to rawSeqFileNameReplaceRegex immediately below
    private String rawSeqFileNameReplaceRegex =
            "(?i)" + // case insensitve
            "\\.fq" + "$|"
            + "\\.fq\\.gz" + "$|"
            + "\\.fastq" + "$|"
            + "_fastq\\.txt" + "$|"
            + "_fastq\\.gz" + "$|"
            + "_fastq\\.txt\\.gz" + "$|"
            + "_sequence\\.txt" + "$|"
            + "_sequence\\.txt\\.gz" + "$|"
            + "_qseq\\.txt" + "$|"
            + "_qseq\\.txt\\.gz" + "$";
    private String noMatchingRawSeqFileNamesMessage =
            "Couldn't find any files that end with "
            + "\".fq\", "
            + "\".fq.gz\", "
            + "\".fastq\", "
            + "\"_fastq.txt\", "
            + "\"_fastq.gz\", "
            + "\"_fastq.txt.gz\", "
            + "\"_sequence.txt\", "
            + "\"_sequence.txt.gz\", "
            + "\"_qseq.txt\", or "
            + "\"_qseq.txt.gz\" "
            + "in the supplied directory: ";

    private byte resolveGenoFromCallPair(byte currGeno, byte newAllele) {
        TreeSet<Byte> alleles = new TreeSet<Byte>();
        byte[] currAlleles = AlignmentUtils.getDiploidValues(currGeno);
        alleles.add(currAlleles[0]);
        alleles.add(currAlleles[1]);
        alleles.add(newAllele);
        if (alleles.size() > 2) {
            return Alignment.UNKNOWN_DIPLOID_ALLELE;
        }
        byte[] alleleArray = new byte[alleles.size()];
        int i = 0;
        for (Byte allele : alleles) {
            alleleArray[i++] = allele.byteValue();
        }
        if (alleleArray.length == 2) {
            return AlignmentUtils.getDiploidValue(alleleArray[0], alleleArray[1]);
        }
        return AlignmentUtils.getDiploidValue(alleleArray[0], alleleArray[0]);
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
}
