/*
 * SeqToGenosPlugin
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
import net.maizegenetics.pal.alignment.ExportUtils;
import net.maizegenetics.pal.alignment.Locus;
import java.awt.Frame;
import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.TreeMap;
import java.util.TreeSet;
import javax.swing.ImageIcon;
import net.maizegenetics.gbs.maps.TOPMInterface;
import net.maizegenetics.gbs.maps.TOPMUtils;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.MutableNucleotideDepthAlignment;
import net.maizegenetics.pal.ids.SimpleIdGroup;
import net.maizegenetics.util.VCFUtil;
import org.apache.commons.math.distribution.BinomialDistributionImpl;
import org.apache.log4j.Logger;

/**
 * This plugin converts all of the fastq (and/or qseq) files in the input
 * folder and keyfile to genotypes and adds these to a genotype file in HDF5 format.
 * We refer to this step as the "Production Pipeline".
 * 
 * The output format is HDF5 genotypes with allelic depth stored. SNP calling 
 * is quantitative with the option of using either the Glaubitz/Buckler binomial
 * method (pHet/pErr > 1 = het), or the VCF/Stacks method.
 * 
 * Samples on multiple lanes with the same LibraryPrepID are merged prior to 
 * SNP calling (so that SNP calling is based upon all available reads).
 *
 * It requires a TOPM with variants added from a previous "Discovery Pipeline"
 * run.  In binary topm or HDF5 format (TOPMInterface).
 *
 * @author jcg233
 */
public class SeqToGenosPlugin extends AbstractPlugin {

    private final Logger myLogger = Logger.getLogger(SeqToGenosPlugin.class);
    private ArgsEngine myArgsEngine = null;
    private String[] myRawSeqFileNames = null;
    private String myKeyFile = null;
    private String myEnzyme = null;
    private String myOutputDir = null;
    private TOPMInterface topm = null;
    private int maxDivergence = 0;
    private int[] chromosomes = null;
    private boolean fastq = true;
    private HashMap<String,Integer> KeyFileColumns = new HashMap<String,Integer>();
    private TreeMap<String,Boolean> FlowcellLanes = new TreeMap<String,Boolean>();  // true = corresponding fastq (or qseq) file is present in input directory
    private TreeMap<String,String> FullNameToFinalName = new TreeMap<String,String>();
    private TreeMap<String,ArrayList<String>> LibraryPrepIDToFlowCellLanes = new TreeMap<String,ArrayList<String>>();
    private TreeMap<String,String> LibraryPrepIDToSampleName = new TreeMap<String,String>();
    private HashMap<String,Integer> FinalNameToTaxonIndex = new HashMap<String,Integer>();
    private MutableNucleotideDepthAlignment[] genos = null;
    private TreeMap<String,Integer> RawReadCountsForFullSampleName = new TreeMap<String,Integer>();
    private TreeMap<String,Integer> RawReadCountsForFinalSampleName = new TreeMap<String,Integer>();
    private boolean vcfL = false;  // if true, use vcf likelihood method for calling hets
    private double errorRate = 0.01;
    private final static int maxCountAtGeno = 500;  // maximum value for likelihoodRatioThreshAlleleCnt[] lookup table
    private static int[] likelihoodRatioThreshAlleleCnt = null;  // index = sample size; value = min count of less tagged allele for likelihood ratio > 1
    // if less tagged allele has counts < likelihoodRatioThreshAlleleCnt[totalCount], call it a homozygote
    // where likelihood ratio = (binomial likelihood het) / (binomial likelihood all less tagged alleles are errors)

    public SeqToGenosPlugin() {
        super(null, false);
    }

    public SeqToGenosPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    private void printUsage() {
        myLogger.info(
            "\nThe options for the TASSEL SeqToGenosPlugin are as follows:\n"
            + "  -i   Input directory containing fastq AND/OR qseq files\n"
            + "  -k   Barcode key file\n"
            + "  -m   Physical map file containing alignments and variants (production TOPM)\n"
            + "  -e   Enzyme used to create the GBS library\n"
            + "  -vL  Use VCF likelihood method to call heterozygotes (default: use tasselGBS likelihood ratio method)"
            + "  -o   Output directory\n"
//            + "  -d  Maximum divergence (edit distance) between new read and previously mapped read (Default: 0 = perfect matches only)\n"  // NOT IMPLEMENTED YET
        );
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
            myArgsEngine.add("-m", "--physical-map", true);
            myArgsEngine.add("-e", "--enzyme", true);
            myArgsEngine.add("-vL", "--VCF-likelihood", false);
            myArgsEngine.add("-o", "--output-directory", true);
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
                Arrays.sort(myRawSeqFileNames);
                myLogger.info("RawReadsToHapMapPlugin: setParameters: \n\nThe following GBS raw sequence data files were found in the input folder (and sub-folders):");
                for (String filename : myRawSeqFileNames) {
                    System.out.println("   "+filename);
                }
                System.out.println("\n");
            }
        } else {
            printUsage();
            throw new IllegalArgumentException("Please specify an input directory containing fastq (or qseq) files (option -i).");
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
        if (myArgsEngine.getBoolean("-vL")) {
            vcfL = true;
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
            if (topm.getSize()==0) {
                throw new IllegalStateException("TagsOnPhysicalMap file not available or is empty");
            }
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
        if ((myEnzyme == null) || (myEnzyme.length() == 0)) {
            printUsage();
            throw new IllegalStateException("performFunction: enzyme must be set.");
        }
        if (!vcfL) {
            setLikelihoodThresh(errorRate);
        }
        readKeyFile();
        matchKeyFileToAvailableRawSeqFiles();
        setUpMutableNucleotideAlignmentsWithDepth();
        readRawSequencesAndRecordDepth();
        callGenotypes();
        writeHapMapFiles();
        writeReadsPerSampleReports();
        return null;
    }

    private void readRawSequencesAndRecordDepth() {
        for (int fileNum = 0; fileNum < myRawSeqFileNames.length; fileNum++) {
            int[] counters = {0, 0, 0, 0, 0, 0}; // 0:allReads 1:goodBarcodedReads 2:goodMatched 3:perfectMatches 4:imperfectMatches 5:singleImperfectMatches
            ParseBarcodeRead thePBR = setUpBarcodes(fileNum);
            if (thePBR == null || thePBR.getBarCodeCount() == 0) {
                System.out.println("No barcodes found. Skipping this flowcell lane.");
                continue;
            }
            myLogger.info("Looking for known SNPs in sequence reads...");
            String temp = "Nothing has been read from the raw sequence file yet";
            BufferedReader br = getBufferedReaderForRawSeqFile(fileNum);
            try {
                while ((temp = br.readLine()) != null) {
                    if (counters[0] % 1000000 == 0) {
                        reportProgress(counters);
                    }
                    ReadBarcodeResult rr = readSequenceRead(br, temp, thePBR, counters);
                    if (rr != null) {
                        counters[1]++;  // goodBarcodedReads
                        RawReadCountsForFullSampleName.put(rr.getTaxonName(),RawReadCountsForFullSampleName.get(rr.getTaxonName())+1);
                        RawReadCountsForFinalSampleName.put(FullNameToFinalName.get(rr.getTaxonName()),RawReadCountsForFinalSampleName.get(FullNameToFinalName.get(rr.getTaxonName()))+1);
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
                        int taxonIndex = FinalNameToTaxonIndex.get(FullNameToFinalName.get(rr.getTaxonName()));
                        incrementDepthForTagVariants(tagIndex, taxonIndex);
                    }
                }
                br.close();
            } catch (Exception e) {
                System.out.println("Catch in readRawSequencesAndRecordDepth() at nReads=" + counters[0] + " e=" + e);
                System.out.println("Last line read: "+temp);
                e.printStackTrace();
            }
            System.out.println("Total number of reads in lane=" + counters[0]);
            System.out.println("Total number of good, barcoded reads=" + counters[1]);
            System.out.println("Finished reading " + (fileNum+1) + " of " + myRawSeqFileNames.length + " sequence files: " + myRawSeqFileNames[fileNum] + "\n");
        }
    }
    
    private void readKeyFile() {
        FlowcellLanes.clear();
        FullNameToFinalName.clear();
        LibraryPrepIDToFlowCellLanes.clear();
        LibraryPrepIDToSampleName.clear();
        String inputLine = "Nothing has been read from the keyfile yet";
        try {
            BufferedReader br = new BufferedReader(new FileReader(myKeyFile), 65536);
            int currLine = 0;
            while ((inputLine = br.readLine()) != null) {
                if (currLine == 0) {
                    parseKeyFileHeader(inputLine);
                } else {
                    populateKeyFileFields(inputLine);
                }
                currLine++;
            }
        } catch (Exception e) {
            System.out.println("Couldn't read key file: " + e);
            System.out.println("Last line read from key file: " + inputLine);
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    private void parseKeyFileHeader(String headerLine) {
        headerLine.trim();
        String[] header = headerLine.split("\\t");
        KeyFileColumns.clear();
        for (int col = 0; col < header.length; col++) {
            if (header[col].equalsIgnoreCase("Flowcell")) {
                KeyFileColumns.put("Flowcell", col);
            } else if (header[col].equalsIgnoreCase("Lane")) {
                KeyFileColumns.put("Lane", col);
            } else if (header[col].equalsIgnoreCase("Barcode")) {
                KeyFileColumns.put("Barcode", col);
            } else if (header[col].equalsIgnoreCase("DNASample") || header[col].equalsIgnoreCase("Sample")) {
                KeyFileColumns.put("Sample", col);
            } else if (header[col].equalsIgnoreCase("LibraryPrepID")) {
                KeyFileColumns.put("LibPrepID", col);
            } else if (header[col].equalsIgnoreCase("Enzyme")) {
                KeyFileColumns.put("Enzyme", col);
            }
        }
        if (!confirmKeyFileHeader()) {
            throwBadKeyFileError();
        }
    }
    
    private boolean confirmKeyFileHeader() {
        if (!KeyFileColumns.containsKey("Flowcell"))
            return false;
        if (!KeyFileColumns.containsKey("Lane"))
            return false;
        if (!KeyFileColumns.containsKey("Barcode"))
            return false;
        if (!KeyFileColumns.containsKey("Sample"))
            return false;
        if (!KeyFileColumns.containsKey("LibPrepID"))
            return false;
        if (!KeyFileColumns.containsKey("Enzyme"))
            return false;
        return true;
    }
    
    private void throwBadKeyFileError() {
        String badKeyFileMessage =
            "The keyfile does not conform to expections.\n" +
            "It must contain columns with the following (exact) headers:\n"+
            "   Flowcell\n"+
            "   Lane\n" +
            "   Barcode\n" +
            "   DNASample (or \"Sample\")\n" +
            "   LibraryPrepID\n" +
            "   Enzyme\n" +
            "\n";
        try {
            throw new IllegalStateException(badKeyFileMessage);
        } catch (Exception e) {
            System.out.println("Couldn't read key file: " + e);
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    private void populateKeyFileFields(String keyFileLine) {
        keyFileLine.trim();
        String[] cells = keyFileLine.split("\\t");
 
        String sample = cells[KeyFileColumns.get("Sample")];
        String libPrepID = cells[KeyFileColumns.get("LibPrepID")];
        String fullName=sample+":"+cells[KeyFileColumns.get("Flowcell")]+":"+cells[KeyFileColumns.get("Lane")]+":"+libPrepID;
        FullNameToFinalName.put(fullName, fullName);

        String flowCellLane = cells[KeyFileColumns.get("Flowcell")]+":"+cells[KeyFileColumns.get("Lane")];
        FlowcellLanes.put(flowCellLane, false);
        ArrayList<String> flowcellLanesForLibPrep = LibraryPrepIDToFlowCellLanes.get(libPrepID);
        if (flowcellLanesForLibPrep == null) {
            LibraryPrepIDToFlowCellLanes.put(libPrepID, flowcellLanesForLibPrep = new ArrayList<String>());
        }
        flowcellLanesForLibPrep.add(flowCellLane);

        String prevSample = LibraryPrepIDToSampleName.get(libPrepID);
        if (prevSample == null) {
            LibraryPrepIDToSampleName.put(libPrepID, sample);
        } else if (!prevSample.contentEquals(sample)) {
            try {
                throw new IllegalStateException("\nThe key file contains different Sample names (\""+prevSample+"\" and \""+sample+"\") for the sample LibraryPrepID ("+libPrepID+")\n\n");
            } catch (Exception e) {
                System.out.println("Error in key file: " + e);
                e.printStackTrace();
                System.exit(1);
            }
        }
    }
    
    private void matchKeyFileToAvailableRawSeqFiles() {
        System.out.println("\nThe following raw sequence files in the input directory conform to one of our file naming conventions and have corresponding samples in the barcode key file:");
        for (int fileNum = 0; fileNum < myRawSeqFileNames.length; fileNum++) {
            String[] flowcellLane = parseRawSeqFileName(myRawSeqFileNames[fileNum]);
            if (flowcellLane != null && FlowcellLanes.containsKey(flowcellLane[0]+":"+flowcellLane[1])) {
                FlowcellLanes.put(flowcellLane[0]+":"+flowcellLane[1], true);  // change from false to true
                System.out.println("  "+myRawSeqFileNames[fileNum]);
            }
        }
        System.out.println("\n");
    }

    private ParseBarcodeRead setUpBarcodes(int fileNum) {
        System.gc();
        System.out.println("\nWorking on GBS raw sequence file: " + myRawSeqFileNames[fileNum]);
        fastq = true;
        if (myRawSeqFileNames[fileNum].substring(myRawSeqFileNames[fileNum].lastIndexOf(File.separator)).contains("qseq")) {
            fastq = false;
        }
        if (fastq) {
            System.out.println("\tThis file is assumed to be in fastq format");
        } else {
            System.out.println("\tThis file contains 'qseq' in its name so is assumed to be in qseq format");
        }
        String[] flowcellLane = parseRawSeqFileName(myRawSeqFileNames[fileNum]);
        if (flowcellLane == null) {
            return null;
        } else {
            ParseBarcodeRead thePBR = new ParseBarcodeRead(myKeyFile, myEnzyme, flowcellLane[0], flowcellLane[1]);
            System.out.println("Total barcodes found in key file for this lane:" + thePBR.getBarCodeCount());
            return thePBR;
        }
    }
    
    /**
     * Parses out the flowcell and lane from the raw GBS sequence filename (fastq or qseq file)
     * @param rawSeqFileName
     * @return String[2] where element[0]=flowcell and element[1]=lane
     */
    private String[] parseRawSeqFileName(String rawSeqFileName) {
        File rawSeqFile = new File(rawSeqFileName);
        String[] FileNameParts = rawSeqFile.getName().split("_");
        if (FileNameParts.length == 3) {
            return new String[] {FileNameParts[0], FileNameParts[1]};
        } else if (FileNameParts.length == 4) {
            return new String[] {FileNameParts[0], FileNameParts[2]};
        } else if (FileNameParts.length == 5) {
            return new String[] {FileNameParts[1], FileNameParts[3]};
        } else {
            printFileNameConventions(rawSeqFileName);
            return null;
        }
    }

    private void setUpMutableNucleotideAlignmentsWithDepth() {
        String[] finalSampleNames = getFinalSampleNames();
        myLogger.info("\nCounting sites in TOPM file.");
        ArrayList<int[]> uniquePositions = getUniquePositions();
        myLogger.info("\nCreating alignment objects to hold the genotypic data (one per chromosome in the TOPM).");
        genos = new MutableNucleotideDepthAlignment[chromosomes.length];
        for (int i = 0; i < genos.length; i++) {
            genos[i] = MutableNucleotideDepthAlignment.getInstance(new SimpleIdGroup(finalSampleNames), uniquePositions.get(i).length); // keeps depth for all 6 alleles
        }
        generateQuickTaxaLookup(finalSampleNames);
        myLogger.info("   Adding sites from the TOPM file to the alignment objects.");
        for (int i = 0; i < genos.length; i++) {
            int currSite = 0;
            for (int j = 0; j < uniquePositions.get(i).length; j++) {
                String chromosome = Integer.toString(chromosomes[i]);
                genos[i].addSite(currSite);
                genos[i].setLocusOfSite(currSite, new Locus(chromosome, chromosome, -1, -1, null, null));
                genos[i].setPositionOfSite(currSite, uniquePositions.get(i)[j]);
                currSite++;
            }
            genos[i].clean();
        }
    }
    
    private void generateQuickTaxaLookup(String[] finalSampleNames) {
        for (int taxonIndex=0; taxonIndex<finalSampleNames.length; taxonIndex++) {
            FinalNameToTaxonIndex.put(finalSampleNames[taxonIndex], taxonIndex);
        }
    }
    
    private String[] getFinalSampleNames() {
        TreeSet<String> finalSampleNamesTS = new TreeSet<String>(); // this will keep the names sorted
        TreeSet<String> samplesInKeyFileWithNoRawSeqFile = new TreeSet<String>();
        for (String LibPrepID : LibraryPrepIDToSampleName.keySet()) {
            String sample = LibraryPrepIDToSampleName.get(LibPrepID);
            ArrayList<String> flowcellLanesForLibPrep = LibraryPrepIDToFlowCellLanes.get(LibPrepID);
            String tempFullName="(NoCorrespondingRawSeqFileForLibPrepID):"+LibPrepID;
            int nRepSamplesWithRawSeqFile = 0;
            for (String flowcellLane : flowcellLanesForLibPrep) {
                if (FlowcellLanes.get(flowcellLane)) { // is fastq (or qseq) file available?
                    nRepSamplesWithRawSeqFile++;
                    tempFullName = sample+":"+flowcellLane+":"+LibPrepID;
                    RawReadCountsForFullSampleName.put(tempFullName, 0);
                } else {
                    samplesInKeyFileWithNoRawSeqFile.add(sample+":"+flowcellLane+":"+LibPrepID);
                }
            }
            if (nRepSamplesWithRawSeqFile == 1) {
                finalSampleNamesTS.add(tempFullName);
                RawReadCountsForFinalSampleName.put(tempFullName, 0);
            } else if (nRepSamplesWithRawSeqFile > 1) {
                String finalName = sample+":MRG:"+nRepSamplesWithRawSeqFile+":"+LibPrepID;
                finalSampleNamesTS.add(finalName);
                RawReadCountsForFinalSampleName.put(finalName, 0);
                for (String flowcellLane : flowcellLanesForLibPrep) {
                    if (FlowcellLanes.get(flowcellLane)) {
                        FullNameToFinalName.put(sample+":"+flowcellLane+":"+LibPrepID, finalName);
                    }
                }
            }
            if (samplesInKeyFileWithNoRawSeqFile.size() > 0) {
                reportOnMissingSamples(samplesInKeyFileWithNoRawSeqFile);
            }
        }
        return finalSampleNamesTS.toArray(new String[0]);
    }
    
    private void reportOnMissingSamples(TreeSet<String> samplesInKeyFileWithNoRawSeqFile) {
        System.out.println("\nThe follow samples in the barcode key file will be absent from the results because there is no corresponding raw sequence (fastq or qseq) file in the input directory:");
        for (String missingSample : samplesInKeyFileWithNoRawSeqFile) {
            System.out.println("   "+missingSample);
        }
        System.out.println("\n");
    }

    private ArrayList<int[]> getUniquePositions() {
        ArrayList<int[]> uniquePositions = new ArrayList<int[]>();
        chromosomes = topm.getChromosomes();
        for (int i = 0; i < chromosomes.length; i++) {
            uniquePositions.add(topm.getUniquePositions(chromosomes[i]));
        }
        return uniquePositions;
    }

    private BufferedReader getBufferedReaderForRawSeqFile(int fileNum) {
        BufferedReader br = null;
        try {
            if (myRawSeqFileNames[fileNum].endsWith(".gz")) {
                br = new BufferedReader(new InputStreamReader(new MultiMemberGZIPInputStream(new FileInputStream(myRawSeqFileNames[fileNum]))));
            } else {
                br = new BufferedReader(new FileReader(myRawSeqFileNames[fileNum]), 65536);
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
            + "  goodBarcodedReads:" + counters[1]
            + "  goodMatchedToTOPM:" + counters[2]
//            + "  perfectMatches:" + counters[3]
//            + "  nearMatches:" + counters[4]
//            + "  uniqueNearMatches:" + counters[5]
        );
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

    private void incrementDepthForTagVariants(int tagIndex, int taxonIndex) {
        int chromosome = topm.getChromosome(tagIndex);
        if (chromosome == TOPMInterface.INT_MISSING) {
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
            if ((newBase == TOPMInterface.BYTE_MISSING) || (newBase == Alignment.UNKNOWN_ALLELE)) {
                continue;
            }
            int offset = topm.getVariantPosOff(tagIndex, variant);
            int pos = startPos + offset;
            int currSite = genos[chrIndex].getSiteOfPhysicalPosition(pos, locus);
            if (currSite < 0) {
                continue;
            }
            genos[chrIndex].incrementDepthForAllele(taxonIndex, currSite, newBase);
        }
    }
    
    private void callGenotypes() {
        System.out.print("\nCalling genotypes...");
        for (int taxon = 0, nTaxa = genos[0].getSequenceCount(); taxon < nTaxa; taxon++) {
            for (int chrIndex = 0, nChrs = genos.length; chrIndex < nChrs; chrIndex++) {
                byte[] callsForTaxon = resolveGenosForTaxon(genos[chrIndex].getDepthsForAllSites(taxon));
                genos[chrIndex].setBaseRange(taxon, 0, callsForTaxon);
            }
        }
        System.out.print("   ...done\n");
    }
    
    private byte[] resolveGenosForTaxon(byte[][] depthsForTaxon) {
        int nAlleles = depthsForTaxon.length;
        byte[] depthsAtSite = new byte[nAlleles];
        int nSites = depthsForTaxon[0].length;
        byte[] genos = new byte[nSites];
        for (int site = 0; site < nSites; site++) {
            for (int allele = 0; allele < nAlleles; allele++) {
                depthsAtSite[allele] = depthsForTaxon[allele][site];
            }
            genos[site] = resolveGeno(depthsAtSite);
        }
        return genos;
    }
    
    private byte resolveGeno(byte[] depths) {
        if (vcfL) {
            int nAlleles = depths.length;
            byte[] alleles = new byte[nAlleles];
            int[] intDepths = new int[nAlleles];
            for (byte a = 0; a < nAlleles; a++) {
                alleles[a] = a;
                intDepths[a] = depths[a];
            }
            return VCFUtil.resolveVCFGeno(alleles, intDepths);
        }
        int count = 0;
        for (int a = 0; a < depths.length; a++) {
            count += depths[a];
        }
        if (count == 0) {
            return Alignment.UNKNOWN_DIPLOID_ALLELE;
        }
        // check for each possible homozygote
        for (int a = 0; a < depths.length; a++) {
            if ((count - depths[a]) == 0) {
                byte byteA = (byte) a;
                return (byte) ((byteA << 4) | byteA);
            }
        }
        return resolveHetGeno(depths);
    }
    
    private byte resolveHetGeno(byte[] depths) {
        int max = 0;
        byte maxAllele = Alignment.UNKNOWN_ALLELE;
        int nextMax = 0;
        byte nextMaxAllele = Alignment.UNKNOWN_ALLELE;
        for (int a = 0; a < depths.length; a++) {
            if (depths[a] > max) {
                nextMax = max;
                nextMaxAllele = maxAllele;
                max = depths[a];
                maxAllele = (byte) a;
            } else if (depths[a] > nextMax) {
                nextMax = depths[a];
                nextMaxAllele = (byte) a;
            }
        }
        // use the Glaubitz/Buckler LR method (if binomialPHet/binomialPErr > 1, call it a het)
        int totCount = max + nextMax;
        if (totCount < maxCountAtGeno) {
            if (nextMax < likelihoodRatioThreshAlleleCnt[totCount]) {
                return (byte) ((maxAllele << 4) | maxAllele); // call it a homozygote
            } else {
                return (byte) ((maxAllele << 4) | nextMaxAllele); // call it a het
            }
        } else {
            if (nextMax / totCount < 0.1) {
                return (byte) ((maxAllele << 4) | maxAllele); // call it a homozygote
            } else {
                return (byte) ((maxAllele << 4) | nextMaxAllele); // call it a het
            }
        }
    }

    private void writeHapMapFiles() {
        for (int i = 0; i < genos.length; i++) {
            genos[i].clean();
            System.out.println("\n\nDistribution of site summary stats for chromosome "+chromosomes[i]+":\n");
            AlignmentFilterByGBSUtils.getCoverage_MAF_F_Dist(genos[i], false);
            String outFileS = myOutputDir + myKeyFile.substring(myKeyFile.lastIndexOf(File.separator));
            outFileS = outFileS.replaceAll(".txt", "_chr" + chromosomes[i] + ".hmp.txt.gz");
            outFileS = outFileS.replaceAll("_key", "");
            ExportUtils.writeToHapmap(genos[i], false, outFileS, '\t', this);
            System.out.println("Genotypes written to:\n"+outFileS+"\n\n");
        }
    }
    
    private void writeReadsPerSampleReports() {
        System.out.print("\nWriting ReadsPerSample log files...");
        String outFileS = myOutputDir + myKeyFile.substring(myKeyFile.lastIndexOf(File.separator));;
        outFileS = outFileS.replaceAll(".txt", "_ReadsPerSample.log");
        outFileS = outFileS.replaceAll("_key", "");
        try {
            BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new File(outFileS))), 65536);
            bw.write("FullSampleName\tgoodBarcodedReads\n");
            for (String fullSampleName : RawReadCountsForFullSampleName.keySet()) {
                bw.write(fullSampleName+"\t"+RawReadCountsForFullSampleName.get(fullSampleName)+"\n");
            }
            bw.close();
        } catch (Exception e) {
            System.out.println("Couldn't write to ReadsPerSample log file: " + e);
            e.printStackTrace();
            System.exit(1);
        }
        outFileS = outFileS.replaceAll("_ReadsPerSample.log", "_ReadsPerLibPrepID.log");
        try {
            BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new File(outFileS))), 65536);
            bw.write("FinalSampleName\tgoodBarcodedReads\n");
            for (String finalSampleName : RawReadCountsForFinalSampleName.keySet()) {
                bw.write(finalSampleName+"\t"+RawReadCountsForFinalSampleName.get(finalSampleName)+"\n");
            }
            bw.close();
        } catch (Exception e) {
            System.out.println("Couldn't write to ReadsPerLibPrepID log file: " + e);
            e.printStackTrace();
            System.exit(1);
        }
        System.out.print("   ...done\n");
    }

    static void setLikelihoodThresh(double errorRate) {   // initialize the likelihood ratio cutoffs for quantitative SNP calling
        likelihoodRatioThreshAlleleCnt = new int[maxCountAtGeno];
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
