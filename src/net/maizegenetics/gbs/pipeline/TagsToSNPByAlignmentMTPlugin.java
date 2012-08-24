/*
 * TagsToSNPByAlignmentMTPlugin
 */
package net.maizegenetics.gbs.pipeline;

import java.awt.Frame;

import java.io.File;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map.Entry;
import java.util.TreeMap;


import javax.swing.ImageIcon;

import net.maizegenetics.gbs.maps.TagsOnPhysicalMap;
import net.maizegenetics.pal.alignment.MutableNucleotideAlignment;
import net.maizegenetics.gbs.tagdist.TagsByTaxa;
import net.maizegenetics.gbs.tagdist.TagsByTaxaBitFileMap;
import net.maizegenetics.util.ArgsEngine;
import net.maizegenetics.gbs.util.BaseEncoder;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.AlignmentUtils;
import net.maizegenetics.pal.alignment.ExportUtils;
import net.maizegenetics.pal.alignment.Locus;
import net.maizegenetics.pal.alignment.NucleotideAlignmentConstants;
import net.maizegenetics.pal.alignment.BitAlignment;
import net.maizegenetics.pal.ids.SimpleIdGroup;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;

import org.apache.log4j.Logger;

import org.biojava3.alignment.Alignments;
import org.biojava3.alignment.template.Profile;
import org.biojava3.core.sequence.DNASequence;
import org.biojava3.core.sequence.compound.AmbiguityDNACompoundSet;
import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.biojava3.core.util.ConcurrencyTools;

/**
 * This class aligns tags at the same physical location against one another, calls SNPs,
 * and then outputs the SNPs to a HapMap file.  
 * 
 * It is multi-threaded, as there are substantial speed increases with it.
 * @author edbuckler
 */
public class TagsToSNPByAlignmentMTPlugin extends AbstractPlugin {

    static int minTaxaCnt = 0;  // a tag must show up in GREATER THAN minTaxaCnt taxa to be included in a sequence alignment
    static int maxSize = 200000;  //normally 200K;
    private double minF = -2.0, minMAF = 0.01;
    private int minMAC = 10;
    //    static boolean ignoreTriallelic=false;
    private boolean inclRare = false;  // false = only call the two most common alleles at a site
    private boolean inclGaps = false;  // false = ignore sites where the major or the 1st minor alleles are gaps
    private boolean isUpdateTOPM = false;
    private final static int maxSNPsPerLocus = 64;
    private final static int maxAlignmentSize = 150;
    static double defaultMinPropTaxaWithLocus = 0.1;
    private static Logger myLogger = Logger.getLogger(TagsToSNPByAlignmentMTPlugin.class);
    TagsOnPhysicalMap theTOPM = null;
    TagsByTaxaBitFileMap theTBT = null;
    File inputFile = null;
    private String inTOPMFile = null;
    private String outTOPMFile = null;
    String outputFilePrefix = null;
    String outHapMap = null;
    int startChr = 0;
    int endChr = 0;
    private static ArgsEngine myArgsEngine = null;
    int minTaxaWithLocus;
    private static String polyN = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";

    public TagsToSNPByAlignmentMTPlugin() {
        super(null, false);
    }

    public TagsToSNPByAlignmentMTPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    public DataSet performFunction(DataSet input) {
        myLogger.info("Finding SNPs in " + inputFile.getAbsolutePath() + ".");
        theTOPM.sortTable(true);
        theTOPM.printRows(5, true, true);
        myLogger.info(String.format("StartChr:%d EndChr:%d %n", startChr, endChr));
        for (int i = startChr; i <= endChr; i++) {
            myLogger.info("\n\nProcessing chromosome " + i + "...");
            String out = outHapMap + ".c" + i;
            myLogger.info("Creating Mutable Alignment");
            MutableNucleotideAlignment theMSA = createMutableAlignment(theTBT, i, i, maxSize + 100);
            runTagsToSNPByAlignment(theTOPM, theTBT, theMSA, out, i, false, minTaxaWithLocus);
            myLogger.info("Finished processing chromosome " + i + "\n\n");
        }
        if (this.isUpdateTOPM) {
            if (outTOPMFile.endsWith(".txt")) {
                theTOPM.writeTextFile(new File(outTOPMFile));
            } else {
                theTOPM.writeBinaryFile(new File(outTOPMFile));
            }
        }
        return null;
    }

    private void printUsage() {
        myLogger.info(
                "\nUsage is as follows:\n"
                + "-i       Input .tbt file\n"
                + "-o       Output directory (default current directory)\n"
                + "-m       TagsOnPhysicalMap file containing genomic position of tags\n"
                + "-mUpd    Update TagsOnPhysicalMap file with allele calls, saved to specified file\n"
                + "-mnF     Minimum F (inbreeding coefficient) (default: " + minF + "  = no filter)\n"
                + "-mnMAF   Minimum minor allele frequency (default: " + minMAF + ")\n"
                + "-mnMAC   Minimum minor allele count (default: " + minMAC + ")\n"
                + "-mnLCov  Minimum locus coverage (proportion of Taxa with a genotype) (default: " + defaultMinPropTaxaWithLocus + ")\n"
                + "-inclRare  Include the rare alleles at site (3 or 4th states) (default: " + inclRare + ")\n"
                + "-inclGaps  Include sites where major or minor allele is a GAP (default: " + inclGaps + ")\n"
                + "-s       Start chromosome\n"
                + "-e       End chromosome");
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
            myArgsEngine.add("-o", "--output-directory", true);
            myArgsEngine.add("-m", "--physical-map", true);
            myArgsEngine.add("-mUpd", "--update-physical-map", true);
            myArgsEngine.add("-mnF", "--minFInbreeding", true);
            myArgsEngine.add("-mnMAF", "--minMinorAlleleFreq", true);
            myArgsEngine.add("-mnMAC", "--minMinorAlleleCount", true);
            myArgsEngine.add("-mnLCov", "--minLocusCov", true);
            myArgsEngine.add("-inclRare", "--includeRare", false);
            myArgsEngine.add("-inclGaps", "--includeGaps", false);
            myArgsEngine.add("-s", "--start-chromosome", true);
            myArgsEngine.add("-e", "--end-chromosome", true);
        }
        myArgsEngine.parse(args);

        if (myArgsEngine.getBoolean("-i")) {
            inputFile = new File(myArgsEngine.getString("-i"));
            if (!inputFile.exists() || !inputFile.isFile()) {
                printUsage();
                throw new IllegalArgumentException("Can't find the TagsByTaxa input file (-i option: " + myArgsEngine.getString("-i") + ").");
            }
            outputFilePrefix = inputFile.getParentFile().getName();
            theTBT = new TagsByTaxaBitFileMap(myArgsEngine.getString("-i"));
        } else {
            printUsage();
            throw new IllegalArgumentException("Please specify a TagsByTaxa input file (-i option).");
        }

        //Set output directory and use the input TagsByTaxa filename for the output file prefix
        if (myArgsEngine.getBoolean("-o")) {
            outHapMap = myArgsEngine.getString("-o") + File.separator + outputFilePrefix;
        } else {
            outHapMap = inputFile.getParent() + File.separator + outputFilePrefix;
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
        if (myArgsEngine.getBoolean("-mnF")) {
            minF = Double.parseDouble(myArgsEngine.getString("-mnF"));
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
        if (myArgsEngine.getBoolean("-inclRare")) {
            inclRare = true;
        }
        if (myArgsEngine.getBoolean("-inclGaps")) {
            inclGaps = true;
        }
        if (myArgsEngine.getBoolean("-s")) {
            startChr = Integer.parseInt(myArgsEngine.getString("-s"));
            if (startChr == 0) {
                printUsage();
                throw new IllegalArgumentException("Error: start chromosome is 0.");
            }
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
            throw new IllegalArgumentException("Error: The start chromosome is higher than the end chromosome.");
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

    public void runTagsToSNPByAlignment(TagsOnPhysicalMap theTOPM, TagsByTaxa theTBT, MutableNucleotideAlignment theMSA,
            String outHapMap, int targetChromo, boolean requireGeneticSupport, int minTaxaWithLocus) {

        System.out.printf("minTaxaWithLocus:%d MinF:%g MinMAF:%g MinMAC:%d %n", minTaxaWithLocus, minF, minMAF, minMAC);
        System.out.printf("includeRare:%s includeGaps:%s %n", inclRare, inclGaps);

        final int cores = Runtime.getRuntime().availableProcessors();
        System.out.println("runTagsToSNPByAlignment processors available:" + cores);
//        ExecutorService pool = Executors.newCachedThreadPool();
//        ThreadPoolExecutor tpe=(ThreadPoolExecutor)pool;
//        ConcurrencyTools.setThreadPool(tpe);
//        System.out.println("Max Pool Size "+tpe.getMaximumPoolSize());
//        System.out.println("Core Pool Size "+tpe.getCorePoolSize());
        long time = System.currentTimeMillis();
        int pauses = 0;

        TreeMap<Integer, String> ts = new TreeMap<Integer, String>();  //TBT index, and the sequence of the tag
        TreeMap<Integer, String> refMarks = new TreeMap<Integer, String>();
        int[] currPos = null;
        int countLoci = 0, totalTaxaCnt = 0;
        int[] longestSeqZeroDiv = {Integer.MIN_VALUE, Integer.MIN_VALUE};  // longestSeqZeroDiv[0] = tbtID, longestSeqZeroDiv[1] = length in bases

        //Loop over all tags in order of physical position
        for (int i = 0; (i < theTOPM.getSize()) && (theMSA.getSiteCount() < (maxSize - 1000)); i++) {
            int ri = theTOPM.getReadIndexForPositionIndex(i);
            if (theTOPM.getChromosome(ri) != targetChromo) {
                continue;    //Skip tags from other chromosomes
            }
            if (requireGeneticSupport && (theTOPM.getMapP(ri) < 2)) {
                continue; //Skip tags with low mapP scores
            }
            if (Arrays.equals(currPos, theTOPM.getPositionArray(ri))) {  // add tag to current ts
                totalTaxaCnt += addTagToTreeMap(ts, ri, longestSeqZeroDiv); // ts<key=tagIndexInTBT,value=TagAsAString>;
            } else {  // finish the current ts and then start a new one with the current tag
                if ((ts.size() > 1) && (totalTaxaCnt > minTaxaWithLocus)) {  // make a sequence alignment and call SNPs from the current ts
                    markRefSeqTag(ts, longestSeqZeroDiv, refMarks);  // mark the longest tag with zero divergence from the reference as "refSeq" (others as "no")
                    SNPsFromSequenceAlignment a = new SNPsFromSequenceAlignment(ts, refMarks, theTBT, currPos[0], currPos[1], currPos[2], theMSA);
                    a.setName("Locus" + countLoci);
//                    while(tpe.getActiveCount()>(cores*10)) {  //4 seems like too little
//                        try{Thread.sleep(0,100000); pauses++;} catch (Exception e) {e.printStackTrace();}
//                    }
//                    tpe.submit(a);
                    a.run();
                    countLoci++;
                    if (theMSA.getSiteCount() % 100 == 0) {
                        double rate = (double) theMSA.getSiteCount() / (double) (System.currentTimeMillis() - time);
                        System.out.printf("Chr:%d Pos:%d Loci=%d SNPs=%d rate=%g SNP/millisec ActiveThreads=%d %n",
                                currPos[0], currPos[2], countLoci, theMSA.getSiteCount(), rate, -1);//, tpe.getActiveCount());
                    }
                }
                ts = new TreeMap<Integer, String>();  // start a new ts for the next locus
                refMarks = new TreeMap<Integer, String>();  // ditto for refMarks
                longestSeqZeroDiv[0] = Integer.MIN_VALUE;
                longestSeqZeroDiv[1] = Integer.MIN_VALUE;
                totalTaxaCnt = 0;
                currPos = theTOPM.getPositionArray(ri);  // currPos[0]=chr (int)  // currPos[1]=strand (byte)  // currPos[2]=minPosition (int)  
                if ((currPos[1] != TagsOnPhysicalMap.byteMissing) && (currPos[2] != TagsOnPhysicalMap.intMissing)) {  // we already know that currPos[0]==targetChromo
                    totalTaxaCnt += addTagToTreeMap(ts, ri, longestSeqZeroDiv);
                } else {
                    currPos = null;  // invalid position
                }
            }
        }
//        System.out.println("Total Pauses or Yields:"+pauses);
//        System.out.println("Main ThreadsCnt:"+Thread.activeCount()+" AlignmentThreadsCnt:"+ConcurrencyTools.getThreadPool().getActiveCount());
//        try {
//            // Wait a while for existing tasks to terminate
//            if (!pool.awaitTermination(5, TimeUnit.SECONDS)) {
//                pool.shutdownNow(); // Cancel currently executing tasks
//                // Wait a while for tasks to respond to being cancelled
//                if (!pool.awaitTermination(5, TimeUnit.SECONDS)) {System.err.println("Pool did not terminate");}
//                else {System.out.println("Pool did terminate");}
//            }
//        } catch (InterruptedException ie) {
//            System.err.println("Pool did not terminate");
//            // (Re-)Cancel if current thread also interrupted
//            pool.shutdownNow();
//            // Preserve interrupt status
//            Thread.currentThread().interrupt();
//        }
//        System.out.println("TC:"+Thread.activeCount()+" BJC"+ConcurrencyTools.getThreadPool().getActiveCount());
        ConcurrencyTools.shutdown();
        System.out.println("Number of marker sites recorded for chr" + targetChromo + ": " + theMSA.getSiteCount());
        if (theMSA.getSiteCount() > 0) {
            //theMSA.sortSiteByPhysicalPosition();
            theMSA.clean();
            ExportUtils.writeToHapmap(theMSA, false, outHapMap, '\t', this);
        }
        System.out.printf("done.\n\n");
        System.out.println("Active thread count:" + Thread.activeCount());
    }

    private int addTagToTreeMap(TreeMap<Integer, String> ts, int TagIndex, int[] longestSeqZeroDiv) {
        long[] tag = theTOPM.getTag(TagIndex);
        int hit = theTBT.getTagIndex(tag);
        int cnt = (hit > -1) ? theTBT.getNumberOfTaxaWithTag(hit) : 0;
        if (cnt > minTaxaCnt) {
            String trimmedSeq = BaseEncoder.getSequenceFromLong(tag).substring(0, theTOPM.getTagLength(TagIndex));
            if (trimmedSeq.length() < 64) {
                trimmedSeq = trimmedSeq + polyN.substring(0, 64 - trimmedSeq.length());
            }
            if (trimmedSeq != null) {
                ts.put(hit, trimmedSeq);  // put sequence into a TreeMap value indexed by its index in the TBT
                if ((theTOPM.getDivergence(TagIndex) == 0) && (theTOPM.getTagLength(TagIndex) > longestSeqZeroDiv[1])) {
                    longestSeqZeroDiv[0] = hit;
                    longestSeqZeroDiv[1] = theTOPM.getTagLength(TagIndex);
                }
                return cnt;
            }
            return 0;
        }
        return 0;
    }

    private void markRefSeqTag(TreeMap<Integer, String> ts, int[] longestSeqZeroDiv, TreeMap<Integer, String> refMarks) {
        for (Entry<Integer, String> id : ts.entrySet()) {
            if (id.getKey() == longestSeqZeroDiv[0]) {
                refMarks.put(id.getKey(), "refSeq");
            } else {
                refMarks.put(id.getKey(), "no");
            }
        }
    }

    /** Fills an array of Locus objects with one locus for each supplied chromosome.  Creates a MutableNucleotideAlignment using
    the locus list and TBT profile as input.  Returns the MSA object. */
    private static MutableNucleotideAlignment createMutableAlignment(TagsByTaxa theTBT, int startChr, int endChr, int maxSites) {
        Locus[] theL = new Locus[endChr - startChr + 1];
        for (int i = 0; i < theL.length; i++) {
            theL[i] = new Locus("" + (startChr + i), "" + (startChr + i), -1, -1, null, null);
        }
        //MutableNucleotideAlignment theMSA = new MutableNucleotideAlignment(theTBT.getTaxaNames(), maxSites, theL);
        MutableNucleotideAlignment theMSA = MutableNucleotideAlignment.getInstance(new SimpleIdGroup(theTBT.getTaxaNames()), maxSites);
        return theMSA;
    }

    private synchronized void addSiteToMutableAlignment(TagsByTaxa theTBT, int chromosome, int strand,
            int startPos, Alignment tagAlignment, MutableNucleotideAlignment theMSA) {
//        int refGenTaxa=theMSA.getIdGroup().whichIdNumber("B73RefGenV2:Sanger:1:A1");
//        System.out.println("TaxaRef:"+theMSA.getTaxaName(refGenTaxa)+" pos:"+startPos);
        for (int s = 0; s < tagAlignment.getSiteCount(); s++) {
            byte[] calls = makeSNPCalls(theTBT, tagAlignment, s);
            byte[] genotypes = null;
            if ((genotypes = isSiteGood(calls, this.inclRare, this.inclGaps, this.minMAF, this.minMAC, this.minF)) == null) {
                continue;
            }
            int currSite = theMSA.getSiteCount();
            Locus chrom = new Locus(String.valueOf(chromosome), String.valueOf(chromosome), -1, -1, null, null);
            theMSA.setLocusOfSite(currSite, chrom);
            int position = (strand == -1) ? (startPos - tagAlignment.getPositionInLocus(s)) : (startPos + tagAlignment.getPositionInLocus(s));
            theMSA.setPositionOfSite(currSite, position);
            //theMSA.setStrandOfSite(currSite, (byte) '+');  // all minus strand genotypes will be complemented to plus strand  
            for (int tx = 0; tx < theTBT.getTaxaCount(); tx++) {
                if (strand == -1) {
                    theMSA.setBase(tx, currSite, complement(calls[tx]));  // complement to plus strand
                    //                   if(tx==refGenTaxa) System.out.printf("setcom pos: %d tx:%d s:%d b:%d %n", startPos, tx, currSite, complement(calls[tx]));
                } else {
                    theMSA.setBase(tx, currSite, calls[tx]);
                    //                   if(tx==refGenTaxa) System.out.printf("setnor pos: %d tx:%d s:%d b:%d %n", startPos, tx, currSite, complement(calls[tx]));
                }
            }
            if (this.isUpdateTOPM) {
                updateTOPM(tagAlignment, s, position, strand, genotypes);
            }
            if (currSite % 100 == 0) {
                System.out.printf("Site:%d Position:%d %n", currSite, position);
            }
        }
    }

    private void updateTOPM(Alignment tagAlignment, int tagAlignSite, int position, int strand, byte[] genotypes) {
        for (int tg = 0; tg < tagAlignment.getSequenceCount(); tg++) {
            byte baseToAdd = tagAlignment.getBase(tg, tagAlignSite);
            boolean matched = false;
            for (byte cb : genotypes) {
                if (baseToAdd == cb) {
                    matched = true;
                    break;
                }
            }
            // so that all tags in the tagAlignment have the same corresponding variants in the TOPM, add a variant no matter what (set to N if needed)
            int tbtTagIndex = Integer.parseInt(tagAlignment.getTaxaName(tg).split("_")[0]);  // taxaName in tagAlignment is set to indexInTheTBT_"refSeq"|"no"
            long[] tag = theTBT.getTag(tbtTagIndex);
            int topmTagIndex = theTOPM.getTagIndex(tag);
            byte offset = (byte) (position - theTOPM.getStartPosition(topmTagIndex));
            if (!matched) {
                baseToAdd = Alignment.UNKNOWN_DIPLOID_ALLELE;
            }
            if (strand == -1) {
                baseToAdd = complement(baseToAdd);  // record everything relative to the plus strand 
            }
            int e = theTOPM.addVariant(topmTagIndex, offset, baseToAdd);
//            System.out.printf("%s %d %d %s %d %n",BaseEncoder.getSequenceFromLong(tag), position, offset, (char)baseToAdd, e);
//            System.out.println(topmTagIndex+":"+tg+":"+theTOPM.printRow(topmTagIndex));
        }
    }

    private byte[] makeSNPCalls(TagsByTaxa theTBT, Alignment tagAlignment, int tagAlignSite) {
        byte[] calls = new byte[theTBT.getTaxaCount()];
        Arrays.fill(calls, Byte.MIN_VALUE);
        for (int tg = 0; tg < tagAlignment.getSequenceCount(); tg++) {
            int tagIndex = Integer.parseInt(tagAlignment.getTaxaName(tg).split("_")[0]);  // taxaName in tagAlignment is set to indexInTheTBT_"refSeq"|"no"
            byte baseToAdd = tagAlignment.getBase(tg, tagAlignSite);
            for (int tx = 0; tx < theTBT.getTaxaCount(); tx++) {
                if (theTBT.getReadCountForTagTaxon(tagIndex, tx) > 0) {
                    byte currentBase = calls[tx];
                    if (currentBase == Byte.MIN_VALUE) {
                        calls[tx] = baseToAdd;
                    } else if (baseToAdd != currentBase) {
                        calls[tx] = resolveSNPByteFromCallPair(currentBase, baseToAdd); // More than 2 alleles in a taxon at a site (eg, A/C/G) --> N
                    }
                }
            }
        }
        for (int tx = 0; tx < theTBT.getTaxaCount(); tx++) {
            if (calls[tx] == Byte.MIN_VALUE) {
                calls[tx] = Alignment.UNKNOWN_DIPLOID_ALLELE;
            }
        }
        return calls;
    }

    /**
     * 
     * @param calls
     * @param includeRare
     * @param includeGap
     * @param minMAF
     * @param minMAC
     * @param minF
     * @return 
     */
    private byte[] isSiteGood(byte[] calls, boolean includeRare, boolean includeGap, double minMAF, int minMAC, double minF) {
        byte[] nuc = {'A', 'C', 'G', 'T', '-', '+', 'N'};
        int[] nucIndex = new int[Byte.MAX_VALUE];
        for (int i = 0; i < nuc.length; i++) {
            nucIndex[nuc[i]] = i;
            //    System.out.println(nuc[i]+":"+(char)nuc[i]);
        }
        int[] cnt = new int[nuc.length];
        for (byte dc : calls) {
            //byte[] cc = IUPACNucleotides.getDiploidValueFromIUPACCode(dc);
            byte[] cc = AlignmentUtils.getDiploidValues(dc);
            cnt[nucIndex[cc[0]]]++;
            cnt[nucIndex[cc[1]]]++;
        }
        int[][] alleles = new int[2][nuc.length - 1];  // alleles[0]=allele; alleles[1]=count   
        for (int i = 0; i < nuc.length - 1; i++) {  // "i<nuc.length-1" stops N from being included
            alleles[0][i] = nuc[i];
            alleles[1][i] = cnt[i];
        }
        boolean change = true;
        while (change) { // sort the alleles by descending frequency
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
        int aCnt = alleles[1][0] + alleles[1][1];
        double theMAF = (double) alleles[1][1] / (double) aCnt;
        if ((theMAF < minMAF) && (alleles[1][1] < minMAC)) {
            return null;
        }
        byte majAllele = (byte) alleles[0][0];
        byte minAllele = (byte) alleles[0][1];
        if (!includeGap && ((majAllele == NucleotideAlignmentConstants.GAP_DIPLOID_ALLELE) || (minAllele == NucleotideAlignmentConstants.GAP_DIPLOID_ALLELE))) {
            return null;
        }
        //byte hetG = IUPACNucleotides.getDegerateSNPByteFromTwoSNPs(majAllele, minAllele);
        byte hetG = AlignmentUtils.getDiploidValue(majAllele, minAllele);
        int hetGCnt = 0;
        for (byte i : calls) {
            if (i == hetG) {
                hetGCnt++;
            }
        }
        int majGCnt = (alleles[1][0] - hetGCnt) / 2; // number of homozygous major allele genotypes
        int minGCnt = (alleles[1][1] - hetGCnt) / 2; // number of homozygous minor allele genotypes
        double propHets = (double) hetGCnt / (double) (hetGCnt + majGCnt + minGCnt);
        double expHets = 2.0 * theMAF * (1 - theMAF);
        double obsF = 1.0 - (propHets / expHets);
        if (obsF < minF) {
            return null;
        }
//        System.out.printf("%d %d %d propHets:%g expHets:%g obsF:%g %n", majGCnt, minGCnt, hetGCnt, propHets, expHets, obsF);
        if (!includeRare) {
            for (int i = 0; i < calls.length; i++) {
                if ((calls[i] != majAllele) && (calls[i] != minAllele) && (calls[i] != hetG)) {
                    calls[i] = Alignment.UNKNOWN_DIPLOID_ALLELE;
                }
            }
        }
        byte[] genotypes = {majAllele, minAllele};
        return genotypes;
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

    private class SNPsFromSequenceAlignment extends Thread {

        TreeMap<Integer, String> tagSeqMap;
        TreeMap<Integer, String> refMarks;
        TagsByTaxa theTBT;
        int chromosome, strand, startPos;
        Alignment tagAlignment = null;
        MutableNucleotideAlignment theMSA;

        public SNPsFromSequenceAlignment(TreeMap<Integer, String> tagSeqMap, TreeMap<Integer, String> refMarks, TagsByTaxa theTBT, int chromosome, int strand,
                int startPos, MutableNucleotideAlignment theMSA) {
            this.tagSeqMap = tagSeqMap;
            this.refMarks = refMarks;
            this.theTBT = theTBT;
            this.chromosome = chromosome;
            this.strand = strand;
            this.startPos = startPos;
            this.theMSA = theMSA;
        }

        public void run() {  // Leave in the comments below
            List<DNASequence> lst = new ArrayList<DNASequence>();
            if (tagSeqMap.size() > maxAlignmentSize) {
                return;   // upper limit of tags
            }
            for (Entry<Integer, String> id : tagSeqMap.entrySet()) {
                DNASequence ds = new DNASequence(id.getValue());
                ds.setOriginalHeader(id.getKey().toString() + "_" + refMarks.get(id.getKey()));  // OriginalHeader will be set to indexInTheTBT_'refSeq'|'no'  
                ds.setCompoundSet(AmbiguityDNACompoundSet.getDNACompoundSet());
                lst.add(ds);
            }
            Profile<DNASequence, NucleotideCompound> profile = Alignments.getMultipleSequenceAlignment(lst);
//            System.out.printf("Clustal1:%d%n%s%n", startPos,profile);
//            Profile<DNASequence, NucleotideCompound> profile2 = Alignments.getMultipleSequenceAlignment(lst);
//            System.out.printf("Clustal2:%d%n%s%n", startPos, profile2);
            String[] aseqs = new String[tagSeqMap.size()];
            String[] names = new String[tagSeqMap.size()];
            boolean refTagWithGaps = false;
            int[] positions = null;
            for (int i = 0; i < aseqs.length; i++) {
                aseqs[i] = profile.getAlignedSequence(i + 1).getSequenceAsString();
                names[i] = profile.getAlignedSequence(i + 1).getOriginalSequence().getOriginalHeader();
                if (names[i].split("_")[1].equals("refSeq")) {  // names were set to indexInTheTBT_"refSeq"|"no"
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
                //aa = SimpleAlignment.getInstance(new SimpleIdGroup(names), aseqs, new IUPACNucleotides(), positions);
                aa = BitAlignment.getNucleotideInstance(new SimpleIdGroup(names), aseqs, null, null, positions, Alignment.DEFAULT_MAX_NUM_ALLELES, null, null, null, false, true);
                //public static SBitAlignment getNucleotideInstance(IdGroup idGroup, String[] data, GeneticMap map, byte[] reference, int[] variableSites,
                //int maxNumAlleles, Locus[] loci, int[] lociOffsets, String[] snpIDs, boolean retainRareAlleles) {
            } else {
                //aa = SimpleAlignment.getInstance(new SimpleIdGroup(names), aseqs, new IUPACNucleotides());
                aa = BitAlignment.getNucleotideInstance(new SimpleIdGroup(names), aseqs, null, null, null, Alignment.DEFAULT_MAX_NUM_ALLELES, null, null, null, false, true);
            }
            //Alignment faa = AnnotatedAlignmentUtils.removeConstantSitesIgnoreMissing(aa);
            Alignment faa = AlignmentUtils.removeSitesBasedOnFreqIgnoreMissing(aa, 0.000001, 1.0, 2);
//            if (refTagWithGaps) { 
//                System.out.println("chr"+chromosome+"  pos:"+startPos+"  strand:"+strand+"  AA\n"+aa.toString().trim());
//                System.out.println("chr"+chromosome+"  pos:"+startPos+"  strand:"+strand+"  FAA\n"+faa.toString()); 
//            }
            //        Alignment iaa=AnnotatedAlignmentUtils.extractIndels(aa, true);
//            System.out.printf("FilteredAlign:%d%n%s%n", startPos,faa.toString());
            //System.out.println(iaa.toString());
            if (faa.getSiteCount() > maxSNPsPerLocus) {
                return;
            }
            addSiteToMutableAlignment(theTBT, chromosome, strand, startPos, faa, theMSA);
        }
    }

    /**
     * Resolves the appropriate IUPACNucleotide from the given callPair (currCall, newCall)
     * 
     * CurrCall is any valid IUPACNucleotide (except '+') while newCall is restricted to A,C,G,T,-,N
     * 
     * @param   currCall, the current genotypic call from previous tag(s) at the locus
     * @param   newCall,  the new allele from the current tag to be combined with currCall
     *                    to make a new genotype
     * @return  resolved byte (valid IUPACNucleotide) 
     */
    public static byte resolveSNPByteFromCallPair(byte currCall, byte newCall) {
        byte snpByte;
        if (newCall == 'A') {
            switch (currCall) {  // conflicts (more than 2 alleles) get set to N
                case 'A':
                    snpByte = 'A';
                    break;
                case 'C':
                    snpByte = 'M';
                    break;
                case 'G':
                    snpByte = 'R';
                    break;
                case 'T':
                    snpByte = 'W';
                    break;
                case 'K':
                    snpByte = 'N';
                    break;
                case 'M':
                    snpByte = 'M';
                    break;
                case 'R':
                    snpByte = 'R';
                    break;
                case 'S':
                    snpByte = 'N';
                    break;
                case 'W':
                    snpByte = 'W';
                    break;
                case 'Y':
                    snpByte = 'N';
                    break;
                case '-':
                    snpByte = '0';
                    break;
                case '+':
                    snpByte = 'N';
                    break; // it should not be possible for currCall to be '+'
                case '0':
                    snpByte = '0';
                    break;
                case 'N':
                    snpByte = 'N';
                    break; // was set to N because of a previous conflict, so should stay as N 
                default:
                    snpByte = 'N';
                    break;
            }
        } else if (newCall == 'C') {
            switch (currCall) {  // conflicts (more than 2 alleles) get set to N
                case 'A':
                    snpByte = 'M';
                    break;
                case 'C':
                    snpByte = 'C';
                    break;
                case 'G':
                    snpByte = 'S';
                    break;
                case 'T':
                    snpByte = 'Y';
                    break;
                case 'K':
                    snpByte = 'N';
                    break;
                case 'M':
                    snpByte = 'M';
                    break;
                case 'R':
                    snpByte = 'N';
                    break;
                case 'S':
                    snpByte = 'S';
                    break;
                case 'W':
                    snpByte = 'N';
                    break;
                case 'Y':
                    snpByte = 'Y';
                    break;
                case '-':
                    snpByte = '0';
                    break;
                case '+':
                    snpByte = 'N';
                    break; // it should not be possible for currCall to be '+'
                case '0':
                    snpByte = '0';
                    break;
                case 'N':
                    snpByte = 'N';
                    break; // was set to N because of a previous conflict, so should stay as N 
                default:
                    snpByte = 'N';
                    break;
            }
        } else if (newCall == 'G') {
            switch (currCall) {  // conflicts (more than 2 alleles) get set to N
                case 'A':
                    snpByte = 'R';
                    break;
                case 'C':
                    snpByte = 'S';
                    break;
                case 'G':
                    snpByte = 'G';
                    break;
                case 'T':
                    snpByte = 'K';
                    break;
                case 'K':
                    snpByte = 'K';
                    break;
                case 'M':
                    snpByte = 'N';
                    break;
                case 'R':
                    snpByte = 'R';
                    break;
                case 'S':
                    snpByte = 'S';
                    break;
                case 'W':
                    snpByte = 'N';
                    break;
                case 'Y':
                    snpByte = 'N';
                    break;
                case '-':
                    snpByte = '0';
                    break;
                case '+':
                    snpByte = 'N';
                    break; // it should not be possible for currCall to be '+'
                case '0':
                    snpByte = '0';
                    break;
                case 'N':
                    snpByte = 'N';
                    break; // was set to N because of a previous conflict, so should stay as N 
                default:
                    snpByte = 'N';
                    break;
            }
        } else if (newCall == 'T') {
            switch (currCall) {  // conflicts (more than 2 alleles) get set to N
                case 'A':
                    snpByte = 'W';
                    break;
                case 'C':
                    snpByte = 'Y';
                    break;
                case 'G':
                    snpByte = 'K';
                    break;
                case 'T':
                    snpByte = 'T';
                    break;
                case 'K':
                    snpByte = 'K';
                    break;
                case 'M':
                    snpByte = 'N';
                    break;
                case 'R':
                    snpByte = 'N';
                    break;
                case 'S':
                    snpByte = 'N';
                    break;
                case 'W':
                    snpByte = 'W';
                    break;
                case 'Y':
                    snpByte = 'Y';
                    break;
                case '-':
                    snpByte = '0';
                    break;
                case '+':
                    snpByte = 'N';
                    break; // it should not be possible for currCall to be '+'
                case '0':
                    snpByte = '0';
                    break;
                case 'N':
                    snpByte = 'N';
                    break; // was set to N because of a previous conflict, so should stay as N 
                default:
                    snpByte = 'N';
                    break;
            }
        } else if (newCall == '-') {  // conflicts (more than 2 alleles) get set to N
            switch (currCall) {
                case 'A':
                    snpByte = '0';
                    break;
                case 'C':
                    snpByte = '0';
                    break;
                case 'G':
                    snpByte = '0';
                    break;
                case 'T':
                    snpByte = '0';
                    break;
                case 'K':
                    snpByte = 'N';
                    break;
                case 'M':
                    snpByte = 'N';
                    break;
                case 'R':
                    snpByte = 'N';
                    break;
                case 'S':
                    snpByte = 'N';
                    break;
                case 'W':
                    snpByte = 'N';
                    break;
                case 'Y':
                    snpByte = 'N';
                    break;
                case '-':
                    snpByte = '-';
                    break;
                case '+':
                    snpByte = 'N';
                    break; // it should not be possible for currCall to be '+'
                case '0':
                    snpByte = '0';
                    break;
                case 'N':
                    snpByte = 'N';
                    break; // was set to N because of a previous conflict, so should stay as N 
                default:
                    snpByte = 'N';
                    break;
            }
        } else if (newCall == 'N') {
            switch (currCall) {
                case 'A':
                    snpByte = 'A';
                    break;
                case 'C':
                    snpByte = 'C';
                    break;
                case 'G':
                    snpByte = 'G';
                    break;
                case 'T':
                    snpByte = 'T';
                    break;
                case 'K':
                    snpByte = 'N';
                    break;
                case 'M':
                    snpByte = 'N';
                    break;
                case 'R':
                    snpByte = 'N';
                    break;
                case 'S':
                    snpByte = 'N';
                    break;
                case 'W':
                    snpByte = 'N';
                    break;
                case 'Y':
                    snpByte = 'N';
                    break;
                case '-':
                    snpByte = '-';
                    break;
                case '+':
                    snpByte = 'N';
                    break; // it should not be possible for currCall to be '+'
                case '0':
                    snpByte = '0';
                    break;
                case 'N':
                    snpByte = 'N';
                    break; // was set to N because of a previous conflict, so should stay as N 
                default:
                    snpByte = 'N';
                    break;
            }
        } else {
            snpByte = 'N';
        }
        return snpByte;
    }
}
