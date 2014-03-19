/*
 * ProductionSNPCallerPlugin
 */
package net.maizegenetics.analysis.gbs;

import cern.colt.list.IntArrayList;
import com.google.common.collect.ImmutableTable;
import com.google.common.collect.Multimap;
import com.google.common.collect.Table;
import com.google.common.collect.TreeMultimap;
import net.maizegenetics.dna.map.*;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;
import net.maizegenetics.util.ArgsEngine;
import net.maizegenetics.util.DirectoryCrawler;
import net.maizegenetics.util.MultiMemberGZIPInputStream;
import org.apache.log4j.Logger;
import javax.swing.*;
import java.awt.*;
import java.io.*;
import java.util.*;
import net.maizegenetics.dna.snp.depth.AlleleDepthUtil;
import net.maizegenetics.dna.snp.genotypecall.BasicGenotypeMergeRule;

/**
 * This plugin converts all of the fastq (and/or qseq) files in the input
 * folder and keyfile to genotypes and adds these to a genotype file in HDF5 format.
 * We refer to this step as the "Production Pipeline".
 * 
 * The output format is HDF5 genotypes with allelic depth stored. SNP calling 
 * is quantitative with the option of using either the Glaubitz/Buckler binomial
 * method (pHet/pErr > 1 = het) (=default), or the VCF/Stacks method.
 * 
 * Merging of samples with the same LibraryPrepID is handled by GenotypeTableBuilder.addTaxon(),
 * with the genotypes re-called based upon the new depths.  Therefore the input 
 * and output GenotypeTableBuilder must be mutable, using closeUnfinished() rather
 * than build()
 * 
 * It requires a TOPM with variants added from a previous "Discovery Pipeline"
 * run.  In binary topm or HDF5 format (TOPMInterface).
 * 
 * TODO add the Stacks likelihood method to BasicGenotypeMergeRule
 *
 * @author Jeff Glaubitz
 * @author Ed Buckler
 */
public class ProductionSNPCallerPlugin extends AbstractPlugin {

    private final Logger myLogger = Logger.getLogger(ProductionSNPCallerPlugin.class);
    private ArgsEngine myArgsEngine = null;
    private String[] myRawSeqFileNames = null;
    private String myKeyFile = null;
    private String myEnzyme = null;
    private String myOutputDir = null;
    private String myTargetHDF5file = null;
    private TOPMInterface topm = null;
    private int maxDivergence = 0;
    private int[] chromosomes = null;
    private boolean fastq = true;
    private Map<String,Integer> keyFileColumns = new HashMap<>();
    private Map<String,Boolean> flowcellLanesInKey = new TreeMap<>();  // flowcellLanes present in key file (stored here as "Flowcell:Lane"); True if corresponding fastq (or qseq) available in input directory
    private Map<String,String> seqFileNameToFlowcellLane = new HashMap<>(); // map from name of fastq or qseq file to "Flowcell:Lane"
    private Set<String> flowcellLanesInKeyAndDir = new TreeSet<>();  // fastq (or qseq) file names present in input directory that have a "Flowcell:Lane" in the key file
    private Map<String,String> fullNameToHDF5Name = new TreeMap<>();
    private Multimap<String, String> libraryPrepIDToFlowCellLanes = TreeMultimap.create();
    private Map<String,String> libraryPrepIDToSampleName = new TreeMap<String,String>();

    private GenotypeTableBuilder genos = null; //output genotype table
    private TaxaList taxaList=null;
    private PositionList myPositionList=null;
    private IntArrayList[] obsTagsForEachTaxon=null;
    private Table<Integer,Integer,Integer> positionToSite = null;  // indices = chrIndices.  For a given position (key), each HashMap provides the site in the MutableNucleotideDepthAlignment (value)
    private int totalNSites = 0;

    //Documentation of read depth per sample (one recored per replicate)
    private Map<String,Integer> rawReadCountsForFullSampleName = new TreeMap<>();
    private Map<String,Integer> matchedReadCountsForFullSampleName = new TreeMap<>();

    private boolean stacksL = false;  // if true, use STACKS likelihood method for calling hets
    private boolean keepOpen = false; // if true, keep the HDF5 genotypes file open for future edits ( i.e., final close is: genos.closeUnfinished() )
    private double errorRate = 0.01;
    private BasicGenotypeMergeRule genoMergeRule = null;
    
    public ProductionSNPCallerPlugin() {
        super(null, false);
    }

    public ProductionSNPCallerPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    private void printUsage() {
        myLogger.info(
            "\n\n\nThe options for the TASSEL ProductionSNPCallerPlugin are as follows:\n"
            + "  -i   Input directory containing fastq AND/OR qseq files\n"
            + "  -k   Barcode key file\n"
            + "  -e   Enzyme used to create the GBS library\n"
            + "  -m   Physical map file containing tags and corresponding variants (production TOPM)\n"
            + "  -o   Output (target) HDF5 genotypes file to add new genotypes to (new file created if it doesn't exist)\n"
            + "  -eR  Average sequencing error rate per base (used to decide between heterozygous and homozygous calls) (default: "+errorRate+")\n"
            + "  -ko  Keep hdf5 genotypes open for future runs that add more taxa or more depth\n (default: finalize hdf5 file)"
//            + "  -sL  Use STACKS likelihood method to call heterozygotes (default: use tasselGBS likelihood ratio method)\n\n\n"
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
            myArgsEngine.add("-i",  "--input-directory", true);
            myArgsEngine.add("-k",  "--key-file", true);
            myArgsEngine.add("-e",  "--enzyme", true);
            myArgsEngine.add("-m",  "--physical-map", true);
            myArgsEngine.add("-o",  "--target-HDF5", true);
            myArgsEngine.add("-eR", "--seqErrRate", true);
            myArgsEngine.add("-ko", "--keep-open", false);
            myArgsEngine.add("-sL", "--STACKS-likelihood", false);
            myArgsEngine.add("-d",  "--divergence", true);
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
                myLogger.info("ProductionSNPCallerPlugin:\n\nThe following GBS raw sequence data files were found in the input folder (and sub-folders):");
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
        if (myArgsEngine.getBoolean("-m")) {
            topm = TOPMUtils.readTOPM(myArgsEngine.getString("-m"));
            if (topm.getSize()==0) {
                throw new IllegalStateException("TagsOnPhysicalMap file not available or is empty");
            }
        } else {
            printUsage();
            throw new IllegalArgumentException("Please specify a TagsOnPhysicalMap file (-m)");
        }
        if (myArgsEngine.getBoolean("-eR")) {
            errorRate = Double.parseDouble(myArgsEngine.getString("-eR"));
        }
        if (myArgsEngine.getBoolean("-o")) {
            myTargetHDF5file = myArgsEngine.getString("-o");
            myOutputDir = myTargetHDF5file.substring(0, myTargetHDF5file.lastIndexOf(File.separator));
            File outDir = new File(myOutputDir);
            if (!outDir.isDirectory()) {
                throw new IllegalArgumentException("The directory containing (or to contain) the target HDF5 genotypes file (option -t) does not exist"); 
            }
        } else {
            printUsage();
            throw new IllegalArgumentException("Please specify a target HDF5 genotypes file (option -t)");
        }
        if (myArgsEngine.getBoolean("-sL")) {
            stacksL = true;
        }
        if (myArgsEngine.getBoolean("-ko")) {
            keepOpen = true;
        }
        if (myArgsEngine.getBoolean("-d")) {
            maxDivergence = Integer.parseInt(myArgsEngine.getString("-d"));
        }
    }

    @Override
    public DataSet performFunction(DataSet input) {
        readKeyFile();  // TODO: read/write full set of metadata
        matchKeyFileToAvailableRawSeqFiles();
        myPositionList = getUniquePositions();
        generateFastSiteLookup(myPositionList);
        int nFilesProcessed = 0;
        for (int fileNum = 0; fileNum < myRawSeqFileNames.length; fileNum++) {
            if ( !flowcellLanesInKeyAndDir.contains(myRawSeqFileNames[fileNum]) ) continue;  // skip fastq/qseq files that are not in the key file
            int[] counters = {0, 0, 0, 0, 0, 0}; // 0:allReads 1:goodBarcodedReads 2:goodMatched 3:perfectMatches 4:imperfectMatches 5:singleImperfectMatches
            System.out.println("\nLooking for known SNPs in sequence reads from file:");
            System.out.println("  "+myRawSeqFileNames[fileNum]+"\n");
            readRawSequencesAndRecordDepth(fileNum, counters);  // TODO: read the machine name from the fastq/qseq file
            callGenotypes();
            ++nFilesProcessed;
            reportTotals(fileNum, counters, nFilesProcessed);
            if(fileNum == myRawSeqFileNames.length-1 && !keepOpen) {
                genos.build();
            } else {
                genos.closeUnfinished();
            }
        }
        writeReadsPerSampleReports();
        return null;
    }

    private void readRawSequencesAndRecordDepth(int fileNum, int[] counters) {
        ParseBarcodeRead thePBR = setUpBarcodes(fileNum);
        if (thePBR == null || thePBR.getBarCodeCount() == 0) {
            System.out.println("No barcodes found. Skipping this raw sequence file.");
            return;
        }
        setUpGenotypeTableBuilder(fileNum);
        obsTagsForEachTaxon = new IntArrayList[taxaList.numberOfTaxa()];
        for (int t = 0; t < obsTagsForEachTaxon.length; t++) {
            obsTagsForEachTaxon[t]=new IntArrayList(750_000); // initial capacity
        }
        String temp = "Nothing has been read from the raw sequence file yet";
        BufferedReader br = getBufferedReaderForRawSeqFile(fileNum);
        try {
            while ((temp = br.readLine()) != null) {
                if (counters[0] % 1000000 == 0)  reportProgress(counters);
                ReadBarcodeResult rr = readSequenceRead(br, temp, thePBR, counters);
                if (rr != null) {
                    counters[1]++;  // goodBarcodedReads
                    rawReadCountsForFullSampleName.put(rr.getTaxonName(),rawReadCountsForFullSampleName.get(rr.getTaxonName())+1);
                    int tagIndex = topm.getTagIndex(rr.getRead());
                    if (tagIndex >= 0)  counters[3]++;  // perfectMatches
                    if (tagIndex < 0 && maxDivergence > 0)  tagIndex = findBestImperfectMatch(rr.getRead(), counters);
                    if (tagIndex < 0)  continue;
                    counters[2]++;  // goodMatched++;
                    matchedReadCountsForFullSampleName.put(rr.getTaxonName(),matchedReadCountsForFullSampleName.get(rr.getTaxonName())+1);
                    int taxonIndex = taxaList.indexOf(fullNameToHDF5Name.get(rr.getTaxonName()));
                    obsTagsForEachTaxon[taxonIndex].add(tagIndex);
                }
            }
            br.close();
        } catch (Exception e) {
            System.out.println("Catch in readRawSequencesAndCallGenos() at nReads=" + counters[0] + " e=" + e);
            System.out.println("Last line read: "+temp);
            e.printStackTrace();
        }
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
   
    private void reportTotals(int fileNum, int[] counters, int nFilesProcessed) {
        System.out.println("Total number of reads in lane=" + counters[0]);
        System.out.println("Total number of good, barcoded reads=" + counters[1]);
        System.out.println("Total number of good, barcoded reads matched to the TOPM=" + counters[2]);
        System.out.println("Finished reading "+nFilesProcessed+" of "+flowcellLanesInKeyAndDir.size()+" sequence files: "+myRawSeqFileNames[fileNum]+"\n");
    }
    
    private void readKeyFile() {
        flowcellLanesInKey.clear();
        flowcellLanesInKeyAndDir.clear();
        fullNameToHDF5Name.clear();
        libraryPrepIDToFlowCellLanes.clear();
        libraryPrepIDToSampleName.clear();
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
        keyFileColumns.clear();
        for (int col = 0; col < header.length; col++) {
            if (header[col].equalsIgnoreCase("Flowcell")) {
                keyFileColumns.put("Flowcell", col);
            } else if (header[col].equalsIgnoreCase("Lane")) {
                keyFileColumns.put("Lane", col);
            } else if (header[col].equalsIgnoreCase("Barcode")) {
                keyFileColumns.put("Barcode", col);
            } else if (header[col].equalsIgnoreCase("DNASample") || header[col].equalsIgnoreCase("Sample")) {
                keyFileColumns.put("Sample", col);
            } else if (header[col].equalsIgnoreCase("LibraryPrepID")) {
                keyFileColumns.put("LibPrepID", col);
            } else if (header[col].equalsIgnoreCase("Enzyme")) {
                keyFileColumns.put("Enzyme", col);
            }
        }
        if (!confirmKeyFileHeader()) {
            throwBadKeyFileError();
        }
    }
    
    private boolean confirmKeyFileHeader() {
        if (!keyFileColumns.containsKey("Flowcell"))
            return false;
        if (!keyFileColumns.containsKey("Lane"))
            return false;
        if (!keyFileColumns.containsKey("Barcode"))
            return false;
        if (!keyFileColumns.containsKey("Sample"))
            return false;
        if (!keyFileColumns.containsKey("LibPrepID"))
            return false;
        if (!keyFileColumns.containsKey("Enzyme"))
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
 
        String sample = cells[keyFileColumns.get("Sample")];
        String libPrepID = cells[keyFileColumns.get("LibPrepID")];
        String fullName=sample+":"+cells[keyFileColumns.get("Flowcell")]+":"+cells[keyFileColumns.get("Lane")]+":"+libPrepID;
        rawReadCountsForFullSampleName.put(fullName, 0);
        matchedReadCountsForFullSampleName.put(fullName, 0);
        fullNameToHDF5Name.put(fullName, sample+":"+libPrepID);

        String flowCellLane = cells[keyFileColumns.get("Flowcell")]+":"+cells[keyFileColumns.get("Lane")];
        flowcellLanesInKey.put(flowCellLane,false);  // false = no corresponding fastq or qseq file (default)
        libraryPrepIDToFlowCellLanes.put(libPrepID,flowCellLane);

        String prevSample = libraryPrepIDToSampleName.get(libPrepID);
        if (prevSample == null) {
            libraryPrepIDToSampleName.put(libPrepID, sample);
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
            if (flowcellLane != null && flowcellLanesInKey.containsKey(flowcellLane[0]+":"+flowcellLane[1])) {
                flowcellLanesInKey.put(flowcellLane[0]+":"+flowcellLane[1], true);  // true = fastq (or qseq) file is available for this flowcell:lane
                seqFileNameToFlowcellLane.put(myRawSeqFileNames[fileNum],flowcellLane[0]+":"+flowcellLane[1]);
                flowcellLanesInKeyAndDir.add(myRawSeqFileNames[fileNum]);
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
    
    private void setUpGenotypeTableBuilder(int fileNum) {
        genoMergeRule = new BasicGenotypeMergeRule(errorRate);
        File hdf5File = new File(myTargetHDF5file);
        if (hdf5File.exists()) {
            System.out.println("\nGenotypes from the raw sequence file:\n  "+myRawSeqFileNames[fileNum]);
            System.out.println("will be added to existing HDF5 file:\n  "+myTargetHDF5file+"\n");
            genos = GenotypeTableBuilder.mergeTaxaIncremental(myTargetHDF5file, genoMergeRule);
        } else {
            System.out.println("\nThe target HDF5 file:\n  "+myTargetHDF5file);
            System.out.println("does not exist. A new HDF5 file of that name will be created \nto hold the genotypes from the raw sequence file:");
            System.out.println("  "+myRawSeqFileNames[fileNum]+"\n");
            genos = GenotypeTableBuilder.getTaxaIncremental(myPositionList, myTargetHDF5file);
        }
        taxaList= new TaxaListBuilder().addAll(getHDF5TaxaNames(fileNum)).sortTaxaAlphabetically().build();
    }

    /**
     * Gets the list of desired HDF5 taxa names ("Sample:LibraryPrepID") for the LibraryPrepIDs in the
     * key file for the corresponding fastq file. 
     * @return String[]
     */
    private String[] getHDF5TaxaNames(int fileNum) {
        String currFlowcellLane = seqFileNameToFlowcellLane.get(myRawSeqFileNames[fileNum]);
        TreeSet<String> hdf5TaxaNamesTS = new TreeSet<String>(); // this will keep the names sorted
        for (String LibPrepID : libraryPrepIDToSampleName.keySet()) {
            Collection<String> flowcellLanesForLibPrep = libraryPrepIDToFlowCellLanes.get(LibPrepID);
            for (String flowcellLane : flowcellLanesForLibPrep) {
                if (flowcellLane.equals(currFlowcellLane)) { // libraryPrepID was run on the current flowcellLane
                    String hdf5Name = libraryPrepIDToSampleName.get(LibPrepID)+":"+LibPrepID;
//                    hdf5Name = hdf5Name.replaceAll(":", " "); // for sorting of taxa based on the short name (" " sorts before any acceptable chars) (matches HDF5 sorting)
                    hdf5TaxaNamesTS.add(hdf5Name);
                }
            }
        }
        String[] hdf5TaxaNames = hdf5TaxaNamesTS.toArray(new String[0]);
//        for (int i = 0; i < finalNames.length; i++) {
//            finalNames[i] = finalNames[i].replaceAll(" ", ":");
//        }
        return hdf5TaxaNames;
    }

    private PositionList getUniquePositions() {
        //todo Move this code to TOPM
        myLogger.info("\nCounting sites in TOPM file");
        PositionListBuilder plb=new PositionListBuilder();
        chromosomes = topm.getChromosomes();
        for (Integer chrNum : chromosomes) {
            Chromosome chr=new Chromosome(chrNum.toString());
            for (int pos : topm.getUniquePositions(chrNum)) {
                Position p=new GeneralPosition.Builder(chr,pos).build();
                plb.add(p);
            }
        }
        PositionList pl=plb.sortPositions().build();
        totalNSites = pl.numberOfSites();
        System.out.println("In total, the TOPM contains "+chromosomes.length+" chromosomes and "+totalNSites+" sites.");
        return pl;
    }
    
    private void generateFastSiteLookup(PositionList pl) {
        ImmutableTable.Builder ptsB=new ImmutableTable.Builder<Chromosome, Integer, Integer>();
        for (int i = 0; i < pl.numberOfSites(); i++) {
            Position position=pl.get(i);
            ptsB.put(position.getChromosome().getChromosomeNumber(),position.getPosition(),i);
        }
        positionToSite=ptsB.build();
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

    private void incrementDepthForTagVariants(int tagIndex, int[][] alleleDepths, int increment) {
        int chromosome = topm.getChromosome(tagIndex);
        if (chromosome == TOPMInterface.INT_MISSING) return;
        int startPos = topm.getStartPosition(tagIndex);
        for (int variant = 0; variant < topm.getMaxNumVariants(); variant++) {
            byte newBase = topm.getVariantDef(tagIndex, variant);
            if ((newBase == TOPMInterface.BYTE_MISSING) || (newBase == GenotypeTable.UNKNOWN_ALLELE)) continue;
            int offset = topm.getVariantPosOff(tagIndex, variant);
            int pos = startPos + offset;
//            int currSite = genos.getSiteOfPhysicalPosition(pos, locus);
            int currSite = positionToSite.get(chromosome,pos);
            if (currSite < 0) continue;
            alleleDepths[newBase][currSite]+=increment;
        }
    }
    
    private void callGenotypes() {
        System.out.println("\nCalling genotypes...");
        for (int currTaxonIndex = 0; currTaxonIndex < obsTagsForEachTaxon.length; currTaxonIndex++) {
            IntArrayList currTagList=obsTagsForEachTaxon[currTaxonIndex];
            currTagList.sort();
            int[][] alleleDepths=new int[NucleotideAlignmentConstants.NUMBER_NUCLEOTIDE_ALLELES][myPositionList.numberOfSites()];
            int prevTag=currTagList.getQuick(0);
            int currInc=0;
            for (int t = 0; t < currTagList.size(); t++) {
                int tag=currTagList.getQuick(t);
                if(tag==prevTag) {currInc++;}
                else {
                    incrementDepthForTagVariants(prevTag,alleleDepths,currInc);
                    prevTag=tag;
                    currInc=1;
                }
            }
            incrementDepthForTagVariants(prevTag,alleleDepths,currInc);
            byte[][] byteDepths = AlleleDepthUtil.depthIntToByte(alleleDepths);
            byte[] taxonGenos = resolveGenosForTaxon(byteDepths);
            genos.addTaxon(taxaList.get(currTaxonIndex),taxonGenos,byteDepths);
            System.out.println("  finished calling genotypes for "+taxaList.get(currTaxonIndex).getName());
        }
        System.out.println("Finished calling genotypes for "+obsTagsForEachTaxon.length+" taxa\n");
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
            genos[site] = genoMergeRule.callBasedOnDepth(depthsAtSite);
        }
        return genos;
    }
    
    private void writeReadsPerSampleReports() {
        System.out.print("\nWriting ReadsPerSample log file...");
        String outFileS = myOutputDir + myKeyFile.substring(myKeyFile.lastIndexOf(File.separator));
        outFileS = outFileS.replaceAll(".txt", "_ReadsPerSample.log");
        outFileS = outFileS.replaceAll("_key", "");
        try {
            BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new File(outFileS))), 65536);
            bw.write("FullSampleName\tgoodBarcodedReads\tgoodReadsMatchedToTOPM\n");
            for (String fullSampleName : rawReadCountsForFullSampleName.keySet()) {
                bw.write(fullSampleName+"\t"+rawReadCountsForFullSampleName.get(fullSampleName)+"\t"+matchedReadCountsForFullSampleName.get(fullSampleName)+"\n");
            }
            bw.close();
        } catch (Exception e) {
            System.out.println("Couldn't write to ReadsPerSample log file: " + e);
            e.printStackTrace();
            System.exit(1);
        }
        System.out.print("   ...done\n");
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
