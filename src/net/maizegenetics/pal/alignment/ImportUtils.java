/*
 * ImportUtils
 */
package net.maizegenetics.pal.alignment;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.logging.Level;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;
import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.SimpleIdGroup;
import net.maizegenetics.prefs.TasselPrefs;
import net.maizegenetics.util.ExceptionUtils;
import net.maizegenetics.util.OpenBitSet;
import net.maizegenetics.util.ProgressListener;
import net.maizegenetics.util.Utils;
import org.apache.log4j.Logger;

/**
 * The class imports Alignment from various file formats.
 *
 * @author terry
 */
public class ImportUtils {

    private static final Logger myLogger = Logger.getLogger(ImportUtils.class);
    private static final Pattern WHITESPACE_PATTERN = Pattern.compile("\\s");
    public static final int NUM_HAPMAP_NON_TAXA_HEADERS = 11;
    public static final int HAPMAP_SNPID_COLUMN_INDEX = 0;
    public static final int HAPMAP_CHROMOSOME_COLUMN_INDEX = 2;
    public static final int HAPMAP_POSITION_COLUMN_INDEX = 3;
    public static final int NUM_VCF_NON_TAXA_COLUMNS = 9;
    public static final int VCF_CHROMOSOME_COLUMN_INDEX = 0;
    public static final int VCF_POSITION_COLUMN_INDEX = 1;
    public static final int VCF_SNPID_COLUMN_INDEX = 2;
    public static final int VCF_REF_COLUMN_INDEX = 3;
    public static final int VCF_ALT_COLUMN_INDEX = 4;
    public static final int VCF_FORMAT_COLUMN_INDEX = 8;

    private ImportUtils() {
        // Utility Class - do not instantiate.
    }
    
    //Main method
    static public void main(String[] args) {
        String testVCFFile = "C:\\Users\\yz79\\Documents\\Tassel\\Tassel4\\Data\\testVCF.txt";
        String outVCFFile = "C:\\Users\\yz79\\Documents\\Tassel\\Tassel4\\Data\\testVCFout.vcf";
        String outHapmap = "C:\\Users\\yz79\\Documents\\Tassel\\Tassel4\\Data\\testHapmap.txt";
        String outHapmap2 = "C:\\Users\\yz79\\Documents\\Tassel\\Tassel4\\Data\\testHapmap2.txt";
        String outHapmap3 = "C:\\Users\\yz79\\Documents\\Tassel\\Tassel4\\Data\\testHapmap3.txt";
        String inHapmap = "C:\\Users\\yz79\\Documents\\Tassel\\TASSELTutorialData3\\TASSELTutorialData\\data\\mdp_genotype.hmp.txt";
        Alignment a = readFromVCF(testVCFFile, null);
        ExportUtils.writeToHapmap(a, true, outHapmap, '\t', null);
        Alignment b = readFromHapmap(inHapmap, null);
        ExportUtils.writeToHapmap(b, true, outHapmap3, '\t', null);
        ExportUtils.writeToVCF(b, outVCFFile, '\t');
        Alignment c = readFromVCF(outVCFFile, null);
        ExportUtils.writeToHapmap(c, true, outHapmap2, '\t', null);
    }
    
    /*
     * Counts number of Header rows in a VCF files
     * Use in conjunction with util.getNumberLines to count numSites
     */
    private static int getNumHeaderRowsVCF(String filename) {
        BufferedReader fileIn = null;
        try {
            int numHeader = 0;
            fileIn = Utils.getBufferedReader(filename, 1000000);
            
            String currLine = fileIn.readLine();
            while(currLine != null) {
                if (currLine.substring(0, 1).equals("#")) {
                    numHeader++;
                } else {
                    break;
                }
                currLine = fileIn.readLine();
            }
            
            return numHeader;
        } catch (Exception e) {
            e.printStackTrace();
            throw new IllegalArgumentException("Error in getNumSiteVCF, unable to read VCF file: " + ExceptionUtils.getExceptionCauses(e));
        } finally {
            try {
                fileIn.close();
            } catch (Exception e) {
                // do nothing
            }
        }
    }
    
    // add support for SNP ID?
    public static Alignment readFromVCF(final String filename, ProgressListener listener) {
        
        int maxAlleles = 3;
        int minPosition = Integer.MAX_VALUE;
        String currLocus = null;
        List<Locus> loci = new ArrayList<Locus>();
        List<Integer> lociOffsets = new ArrayList<Integer>();
        
        Pattern colonPattern = Pattern.compile(":");
        Pattern commaPattern = Pattern.compile(",");
        
        long currentTime = System.currentTimeMillis();
        int numHeader = getNumHeaderRowsVCF(filename);
        int numSites = Utils.getNumberLines(filename)- numHeader;
        myLogger.info("readFromVCF: Number of Header Rows: " + numHeader);
        myLogger.info("readFromVCF: Number of Sites: " + numSites);
        
        long prevTime = currentTime;
        currentTime = System.currentTimeMillis();
        myLogger.info("readFromVCF: Time to count lines: " + ((currentTime - prevTime) / 1000));
        
        
        BufferedReader fileIn = null;
        try {
            fileIn = Utils.getBufferedReader(filename, 1000000);
            String currLine = null;
            for(int i = 0; i < numHeader; i++) {
                currLine = fileIn.readLine();
            }
            
            // get taxa names
            String[] header = WHITESPACE_PATTERN.split(currLine);
            int numTaxa = header.length - NUM_VCF_NON_TAXA_COLUMNS;
            String[] taxaNames = new String[numTaxa];
            System.arraycopy(header, NUM_VCF_NON_TAXA_COLUMNS, taxaNames, 0, numTaxa);
            IdGroup idGroup = new SimpleIdGroup(taxaNames);
            
            String[] snpID = new String[numSites];
            
            MutableVCFAlignment result = MutableVCFAlignment.getInstance(idGroup, numSites, numTaxa, numSites);
            int prevPos = -1;
            int currPos = -1;
            int locusStart = 0;
            
            for(int site = 0; site < numSites; site++) {
                currLine = fileIn.readLine();
                String[] currSite = WHITESPACE_PATTERN.split(currLine);
//                1	6407	S1000_6407	.	A,C	20	PASS	NS=929;DP=2210;AF=0.01,0.01	GT:AD:DP:GQ:PL	./.	0/0:2,0,0:2:79:0,6,72	1/2:0,1,1:2:54:0,5,34
                String temp = currSite[VCF_CHROMOSOME_COLUMN_INDEX];
                currPos = Integer.parseInt(currSite[VCF_POSITION_COLUMN_INDEX]);
                snpID[site] = currSite[VCF_SNPID_COLUMN_INDEX];
                String ref = currSite[VCF_REF_COLUMN_INDEX];
                String alt = currSite[VCF_ALT_COLUMN_INDEX].replaceAll(",", "");
                String format = currSite[VCF_FORMAT_COLUMN_INDEX];
                String[] dataByTaxa = new String[numTaxa];
                System.arraycopy(currSite, NUM_VCF_NON_TAXA_COLUMNS, dataByTaxa, 0, numTaxa);
                
                // find alleles for current site, check to see if number of alleles is supported
                int numAlleles = alt.length() + 1;
                if (numAlleles > maxAlleles) {
                    throw new IllegalStateException("ImportUtils: readFromVCF: number of Alleles is larger than allowed currently in TASSEL: " + numAlleles + " alleles found, " + maxAlleles + " alleles allowed, in line " + (numHeader + site + 1));
                }
                
                String alleleString = ref + alt;
                byte[] alleles = new byte[numAlleles];
                for(int allele = 0; allele < numAlleles; allele++) {
                    String currAllele = alleleString.substring(allele, allele + 1);
                    if (currAllele.equals(".")) {
                        alleles[allele] = (byte) -1;
                    } else if (currAllele.equalsIgnoreCase("A")) {
                        alleles[allele] = NucleotideAlignmentConstants.A_ALLELE;
                    } else if (currAllele.equalsIgnoreCase("C")) {
                        alleles[allele] = NucleotideAlignmentConstants.C_ALLELE;
                    } else if (currAllele.equalsIgnoreCase("G")) {
                        alleles[allele] = NucleotideAlignmentConstants.G_ALLELE;
                    } else if (currAllele.equalsIgnoreCase("T")) {
                        alleles[allele] = NucleotideAlignmentConstants.T_ALLELE;
                    } else if (currAllele.equals("+")) {
                        alleles[allele] = NucleotideAlignmentConstants.INSERT_ALLELE;
                    } else if (currAllele.equals("-")) {
                        alleles[allele] = NucleotideAlignmentConstants.GAP_ALLELE;
                    } else {
                        throw new IllegalStateException("ImportUtils: readFromVCF: a unsupported allele detected in this VCF file: " + currAllele + " in line " + (numHeader + site + 1));
                    }
                }
                result.setCommonAlleles(site, alleles);
                
                // get the possible alleles for each site in to an byte array
//                result.setCommonAlleles(site, values);
                
                // figure out order of format
                int genoIndex = -1;
                int alleleDepthIndex = -1;
                
                String[] formatSplit = colonPattern.split(format);
                int numDataFields = formatSplit.length;
                for(int i = 0; i < formatSplit.length; i++) {
                    if (formatSplit[i].equalsIgnoreCase("GT")) {
                        genoIndex = i;
                    } else if (formatSplit[i].equalsIgnoreCase("AD")) {
                        alleleDepthIndex = i;
                    }
                }
                
                if (genoIndex == -1) {
                    throw new IllegalStateException("ImportUtils: readFromVCF: no genotype data found in this VCF file at line: " + (numHeader + site + 1));
                }
                
                if (alleleDepthIndex == -1) {
                    throw new IllegalStateException("ImportUtils: readFromVCF: no allele depth data found in this VCF file at line: " + (numHeader + site + 1));
                }
                
                for(int taxa = 0; taxa < numTaxa; taxa++) {
                    String[] dataSplit = colonPattern.split(dataByTaxa[taxa]);
                    
                    // for whatever reason if the actual data fields do not match up with the format column
                    // assume the data is unknown
                    byte value = (byte) 0xFF;
                    if (dataSplit.length != numDataFields) {
                        result.setBase(taxa, site, value);
                    } else {
                        String[] genotypes = Pattern.compile("[/|]").split(dataSplit[genoIndex]);
//                        String genotypes = dataSplit[genoIndex].replaceAll("/", "").replaceAll("|", "");
                        if (genotypes.length > 2) {
                            throw new IllegalStateException("ImportUtils: readFromVCF: number of genotypes larger than supported by TASSEL at line: " + (numHeader + site + 1));
                        }
                        for(int i = 0; i < genotypes.length; i++) {
                            value <<= 4;
                            String currGenotype = genotypes[i];
                            if (currGenotype.equals(".")) {
                                value |= 0x0F;
                            } else {
                                int currGenoInt = Integer.parseInt(currGenotype);
                                if (currGenoInt == 14) {
                                    value |= 0x0E;
                                } else {
                                    currGenotype = alleleString.substring(currGenoInt, currGenoInt + 1);
                                    if (currGenotype.equalsIgnoreCase("A")) {
                                        value |= 0x00;
                                    } else if (currGenotype.equalsIgnoreCase("C")) {
                                        value |= 0x01;
                                    } else if (currGenotype.equalsIgnoreCase("G")) {
                                        value |= 0x02;
                                    } else if (currGenotype.equalsIgnoreCase("T")) {
                                        value |= 0x03;
                                    } else if (currGenotype.equals("+")) {
                                        value |= 0x04;
                                    } else if (currGenotype.equals("-")) {
                                        value |= 0x05;
                                    } else {
                                        throw new IllegalStateException("ImportUtils: readFromVCF: a unsupported allele detected in this VCF file: " + currGenotype + " in line " + (numHeader + site + 1));
                                    }
                                }
                            }
                        }
                        result.setBase(taxa, site, value);
                        
                        String alleleDepths = dataSplit[alleleDepthIndex];
                        String[] stringDepths = commaPattern.split(alleleDepths);
                        
                        if (stringDepths.length != alleleString.length()) {
                            throw new IllegalStateException("ImportUtils: readFromVCF: number of allele depth values does not match number of alleles in line: " + (numHeader + site + 1) + " taxa number: " + taxa);
                        }
                        
                        byte[] depths = new byte[stringDepths.length];
                        for(int i = 0; i < stringDepths.length; i++) {
                            int depth = Integer.parseInt(stringDepths[i]);
                            if (depth > 127) {
                                myLogger.info("Depth value for genotype " + i + " had an original value of " + depth + ". Converted to the maximum of 127. In line: " + (numHeader + site + 1) + " taxa number: " + taxa);
                                depth = 127;
                                depths[i] = (byte) depth;
                            }
                        }
                        result.setDepthForAlleles(taxa, site, alleles);
                    }
                }
                        
                result.setPositionOfSite(site, currPos);
                if (currLocus == null) {
                    currLocus = temp;
                    minPosition = currPos;
                } else if (!temp.equals(currLocus)) {
                    Locus newLocus = new Locus(currLocus, currLocus, minPosition, prevPos, null, null);
                    for(int i = locusStart; i < site; i++) {
                        result.setLocusOfSite(i, newLocus);
                    }
                    currLocus = temp;
                    minPosition = currPos;
                    locusStart = site;
                    prevPos = -1;
                }
                
                if (currPos < prevPos) {
                    throw new IllegalStateException("ImportUtils: readFromVCF: Sites are not properly sorted for chromosome: " + currLocus + " at " + currPos + " and " + prevPos);
                }
                
                prevPos = currPos;
            }
            
            if (currLocus != null) {
                Locus newLocus = new Locus(currLocus, currLocus, minPosition, prevPos, null, null);
                for(int i = locusStart; i < numSites; i++) {
                    result.setLocusOfSite(i, newLocus);
                }
            }
            
            return result;
        } catch (Exception e) {
            e.printStackTrace();
            throw new IllegalArgumentException("ImportUtils: readFromVCF: Problem creating Alignment: " + filename + ": " + ExceptionUtils.getExceptionCauses(e));
        } finally {
            try {
                fileIn.close();
            } catch (Exception ex) {
                // do nothing
            }
        }
    }

    public static Alignment readFromHapmap(final String filename, ProgressListener listener) {
        return readFromHapmap(filename, true, listener);
    }

    public static Alignment readFromHapmap(final String filename, boolean isSBit, ProgressListener listener) {

        int minPosition = Integer.MAX_VALUE;
        String currLocus = null;
        List<Locus> loci = new ArrayList<Locus>();
        List<Integer> lociOffsets = new ArrayList<Integer>();

        long currentTime = System.currentTimeMillis();
        int numSites = Utils.getNumberLines(filename) - 1;
        myLogger.info("readFromHapmap: Number of Sites: " + numSites);

        long prevTime = currentTime;
        currentTime = System.currentTimeMillis();
        myLogger.info("readFromHapmap: Time to count lines: " + ((currentTime - prevTime) / 1000));


        BufferedReader fileIn = null;
        try {
            int numThreads = Runtime.getRuntime().availableProcessors();
            ExecutorService pool = Executors.newFixedThreadPool(numThreads);

            fileIn = Utils.getBufferedReader(filename, 1000000);
            String[] header = WHITESPACE_PATTERN.split(fileIn.readLine());
            int lineInFile = 1;
            int numTaxa = header.length - NUM_HAPMAP_NON_TAXA_HEADERS;
            String[] snpIDs = new String[numSites];
            int prevPosition = -1;

            OpenBitSet[][] theData;
            byte[][] alleles = new byte[numSites][TasselPrefs.getAlignmentMaxAllelesToRetain()];
            int numDataRows = TasselPrefs.getAlignmentMaxAllelesToRetain();
            if (TasselPrefs.getAlignmentRetainRareAlleles()) {
                numDataRows++;
            }
            int numSitesToProcess = 1;
            if (isSBit) {
                theData = new OpenBitSet[numDataRows][numSites];
                numSitesToProcess = 1;
            } else {
                theData = new OpenBitSet[numDataRows][numTaxa];
                for (int al = 0; al < numDataRows; al++) {
                    for (int t = 0; t < numTaxa; t++) {
                        theData[al][t] = new OpenBitSet(numSites);
                    }
                }
                numSitesToProcess = 64;
            }

            int[] physicalPositions = new int[numSites];
            int count = 0;
            String[][] tokens = new String[numSitesToProcess][];
            int currentSite = 0;
            for (int site = 0; site < numSites; site++) {

                lineInFile++;

                String input = fileIn.readLine();
                tokens[count] = WHITESPACE_PATTERN.split(input);

                snpIDs[site] = new String(tokens[count][HAPMAP_SNPID_COLUMN_INDEX]);
                int position = Integer.parseInt(tokens[count][HAPMAP_POSITION_COLUMN_INDEX]);
                String temp = new String(tokens[count][HAPMAP_CHROMOSOME_COLUMN_INDEX]);
                if (currLocus == null) {
                    lociOffsets.add(site);
                    currLocus = temp;
                    minPosition = position;
                    prevPosition = -1;
                } else if (!temp.equals(currLocus)) {
                    loci.add(new Locus(currLocus, currLocus, minPosition, prevPosition, null, null));
                    lociOffsets.add(site);
                    currLocus = temp;
                    minPosition = position;
                    prevPosition = -1;
                }

                if (position < prevPosition) {
                    throw new IllegalStateException("ImportUtils: readFromHapmap: Sites are not properly sorted for chromosome: " + currLocus + " at " + position + " and " + prevPosition);
                }

                count++;

                if (count == numSitesToProcess) {
                    pool.execute(ProcessLineOfHapmap.getInstance(alleles, theData, TasselPrefs.getAlignmentRetainRareAlleles(), tokens, count, currentSite, numTaxa, lineInFile, isSBit));
                    count = 0;
                    currentSite += numSitesToProcess;
                    tokens = new String[numSitesToProcess][];
                }

                physicalPositions[site] = position;
                prevPosition = position;

                if (listener != null) {
                    listener.progress((int) (((double) (site + 1) / (double) numSites) * 100.0), null);
                }
            }

            if (count != 0) {
                pool.submit(ProcessLineOfHapmap.getInstance(alleles, theData, TasselPrefs.getAlignmentRetainRareAlleles(), tokens, count, currentSite, numTaxa, lineInFile, isSBit));
            }


            pool.shutdown();
            if (!pool.awaitTermination(6000, TimeUnit.SECONDS)) {
                throw new IllegalStateException("ImportUtils: readFromHapmap: processing threads timed out.");
            }

            if (currLocus != null) {
                loci.add(new Locus(currLocus, currLocus, minPosition, prevPosition, null, null));
            }

            prevTime = currentTime;
            currentTime = System.currentTimeMillis();
            myLogger.info("readFromHapmap: Time to read file: " + ((currentTime - prevTime) / 1000));

            String[] taxaNames = new String[numTaxa];
            System.arraycopy(header, NUM_HAPMAP_NON_TAXA_HEADERS, taxaNames, 0, numTaxa);
            IdGroup idGroup = new SimpleIdGroup(taxaNames);

            Locus[] lociFinal = new Locus[loci.size()];
            loci.toArray(lociFinal);
            int[] offsetsFinal = new int[lociOffsets.size()];
            for (int i = 0; i < lociOffsets.size(); i++) {
                offsetsFinal[i] = ((Integer) lociOffsets.get(i)).intValue();
            }

            Alignment result = BitAlignment.getNucleotideInstance(idGroup, alleles, theData, null, null, physicalPositions, TasselPrefs.getAlignmentMaxAllelesToRetain(), lociFinal, offsetsFinal, snpIDs, TasselPrefs.getAlignmentRetainRareAlleles(), isSBit);

            prevTime = currentTime;
            currentTime = System.currentTimeMillis();
            myLogger.info("readFromHapmap: Time to create Alignment: " + ((currentTime - prevTime) / 1000));

            return result;
        } catch (Exception e) {
            e.printStackTrace();
            throw new IllegalArgumentException("ImportUtils: readFromHapmap: Problem creating Alignment: " + filename + ": " + ExceptionUtils.getExceptionCauses(e));
        } finally {
            try {
                fileIn.close();
            } catch (Exception ex) {
                // do nothing
            }
        }

    }

    public static Alignment readFasta(String filename, boolean isSBit) throws FileNotFoundException, IOException {

        BufferedReader reader = Utils.getBufferedReader(filename);

        List taxa = new ArrayList();
        List sequences = new ArrayList();

        String line = null;
        line = reader.readLine();
        boolean sequence = false;
        int sequenceLength = -1;
        int count = 1;
        while (line != null) {

            line = line.trim();

            if (line.startsWith(";")) {
                line = reader.readLine();
            } else if (line.startsWith(">")) {
                StringTokenizer tokens = new StringTokenizer(line);
                String taxaName = tokens.nextToken();
                if (taxaName.length() == 1) {
                    taxaName = tokens.nextToken();
                } else {
                    taxaName = taxaName.substring(1).trim();
                }
                taxa.add(taxaName);
                sequence = true;
                line = reader.readLine();
            } else if (sequence) {
                StringBuilder builder = new StringBuilder();
                while ((line != null) && (!line.startsWith(">")) && (!line.startsWith(";"))) {
                    line = line.trim().toUpperCase();
                    builder.append(line);
                    line = reader.readLine();
                }
                String temp = builder.toString();
                if (sequenceLength == -1) {
                    sequenceLength = temp.length();
                } else if (sequenceLength != temp.length()) {
                    throw new IllegalStateException("ImportUtils: readFasta: Sequence: " + count + " Differs in Length.");
                }
                sequences.add(temp);
                sequence = false;
                count++;
            } else {
                myLogger.error("readFasta: file: " + filename + " invalid format.");
                throw new IllegalArgumentException("Import: readFasta: invalid format.");
            }

        }

        String[] taxaNames = new String[taxa.size()];
        taxa.toArray(taxaNames);
        IdGroup idGroup = new SimpleIdGroup(taxaNames);

        String[] sequenceArray = new String[sequences.size()];
        sequences.toArray(sequenceArray);

        Locus unknown = new Locus("Unknown", "0", 0, sequenceArray[0].length(), null, null);
        return BitAlignment.getNucleotideInstance(idGroup, sequenceArray, null, null, null, TasselPrefs.getAlignmentMaxAllelesToRetain(), new Locus[]{unknown}, new int[]{0}, null, TasselPrefs.getAlignmentRetainRareAlleles(), isSBit);

    }

    public static Alignment readAlignmentFromSerialGZ(String inFile) {

        Alignment alignment = null;
        long time = System.currentTimeMillis();
        FileInputStream fis = null;
        GZIPInputStream gs = null;
        ObjectInputStream ois = null;
        try {
            File theFile = new File(Utils.addSuffixIfNeeded(inFile, ".serial.gz"));
            myLogger.info("readAlignmentFromSerialGZ: Reading:" + theFile);
            fis = new FileInputStream(theFile);
            gs = new GZIPInputStream(fis);
            ois = new ObjectInputStream(gs);
            alignment = (Alignment) ois.readObject();

        } catch (Exception ee) {
            ee.printStackTrace();
        } finally {
            try {
                ois.close();
                gs.close();
                fis.close();
            } catch (Exception e) {
                // do nothing
            }
        }
        myLogger.info("readAlignmentFromSerialGZ: Time: " + (System.currentTimeMillis() - time) + "  Sites: " + alignment.getSiteCount() + "  Taxa: " + alignment.getSequenceCount());
        return alignment;
    }
}
