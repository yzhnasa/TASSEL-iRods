/*
 * ImportUtils
 */
package net.maizegenetics.dna.snp;

import net.maizegenetics.dna.snp.io.BuilderFromHapMap;
import net.maizegenetics.dna.snp.io.BuilderFromVCF;
import net.maizegenetics.util.ProgressListener;
import net.maizegenetics.util.Utils;
import org.apache.log4j.Logger;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;
import java.util.regex.Pattern;
import net.maizegenetics.dna.map.PositionListBuilder;
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTable;
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTableBuilder;
import net.maizegenetics.taxa.TaxaList;
import net.maizegenetics.taxa.TaxaListBuilder;

/**
 * Methods for importing GenotypeTables from various file formats.
 *
 * @author Terry Casstevens
 * @author Ed Buckler
 */
public class ImportUtils {

    private static final Logger myLogger = Logger.getLogger(ImportUtils.class);
    private static final Pattern WHITESPACE_PATTERN = Pattern.compile("\\s");
    public static final int NUM_HAPMAP_NON_TAXA_HEADERS = 11;

    private ImportUtils() {
        // Utility Class - do not instantiate.
    }

    public static GenotypeTable readGuessFormat(String fileName) {
        try {
            if (fileName.endsWith(".h5")) {
                //return BuilderFromGenotypeHDF5.getBuilder(fileName).build();
                return GenotypeTableBuilder.getInstance(fileName);
            } else if (fileName.endsWith("hmp.txt.gz") || fileName.endsWith("hmp.txt")) {
                return readFromHapmap(fileName, null);
            } else if (fileName.endsWith(".vcf") || fileName.endsWith(".vcf.gz")) {
                return readFromVCF(fileName, null);
            }
            return null;
        } catch (Exception e) {
            System.err.println("Error reading:" + fileName);
            e.printStackTrace();
            return null;
        }

    }

    public static GenotypeTable readFromVCF(final String filename, ProgressListener listener, boolean ignoreDepth) {
        if (ignoreDepth) {
            return BuilderFromVCF.getBuilder(filename).keepDepth().build();
        }
        return BuilderFromVCF.getBuilder(filename).build();
    }

    public static GenotypeTable readFromVCF(final String filename, ProgressListener listener) {
        return readFromVCF(filename, listener, true);
    }

    /**
     * Read GenotypeTable from HapMap file
     *
     * @param filename input HapMap file name
     * @return a genotype table
     */
    public static GenotypeTable readFromHapmap(final String filename) {
        return BuilderFromHapMap.getBuilder(filename).build();
    }

    /**
     * Read GenotypeTable from HapMap file
     *
     * @param filename input HapMap file name
     * @param listener progress listener to track reading rate
     * @return a genotype table
     */
    public static GenotypeTable readFromHapmap(final String filename, ProgressListener listener) {
        return BuilderFromHapMap.getBuilder(filename).build();
    }

    public static GenotypeTable readFromPLink(final String pedFilename, final String mapFilename, ProgressListener listener) {
//        return readFromPLink(pedFilename, mapFilename, true, listener);
        return null; //TODO restore PLINK
    }

//    public static Alignment readFromPLink(final String pedFilename, final String mapFilename, boolean isSBit, ProgressListener listener) {
//
//        int minPosition = Integer.MAX_VALUE;
//        String currLocus = null;
//        List<Chromosome> loci = new ArrayList<Chromosome>();
//        List<Integer> lociOffsets = new ArrayList<Integer>();
//
//        long currentTime = System.currentTimeMillis();
//        int numSites = Utils.getNumberLines(pedFilename) - 1;
//        myLogger.info("readFromHapmap: Number of Sites: " + numSites);
//
//        long prevTime = currentTime;
//        currentTime = System.currentTimeMillis();
//        myLogger.info("readFromHapmap: Time to count lines: " + ((currentTime - prevTime) / 1000));
//
//
//        BufferedReader fileIn = null;
//        try {
//            int numThreads = Runtime.getRuntime().availableProcessors();
//            ExecutorService pool = Executors.newFixedThreadPool(numThreads);
//
//            fileIn = Utils.getBufferedReader(pedFilename, 1000000);
//            String[] header = WHITESPACE_PATTERN.split(fileIn.readLine());
//            int lineInFile = 1;
//            int numTaxa = header.length - NUM_HAPMAP_NON_TAXA_HEADERS;
//            String[] snpIDs = new String[numSites];
//            int prevPosition = -1;
//
//            OpenBitSet[][] theData;
//            byte[][] alleles = new byte[numSites][TasselPrefs.getAlignmentMaxAllelesToRetain()];
//            int numDataRows = TasselPrefs.getAlignmentMaxAllelesToRetain();
//            if (TasselPrefs.getAlignmentRetainRareAlleles()) {
//                numDataRows++;
//            }
//            int numSitesToProcess = 1;
//            if (isSBit) {
//                theData = new OpenBitSet[numDataRows][numSites];
//                numSitesToProcess = 1;
//            } else {
//                theData = new OpenBitSet[numDataRows][numTaxa];
//                for (int al = 0; al < numDataRows; al++) {
//                    for (int t = 0; t < numTaxa; t++) {
//                        theData[al][t] = new OpenBitSet(numSites);
//                    }
//                }
//                numSitesToProcess = 64;
//            }
//
//            int[] physicalPositions = new int[numSites];
//            int count = 0;
//            String[][] tokens = new String[numSitesToProcess][];
//            int currentSite = 0;
//            for (int site = 0; site < numSites; site++) {
//
//                lineInFile++;
//
//                String input = fileIn.readLine();
//                tokens[count] = WHITESPACE_PATTERN.split(input);
//
//                snpIDs[site] = new String(tokens[count][HAPMAP_SNPID_COLUMN_INDEX]);
//                int position = Integer.parseInt(tokens[count][HAPMAP_POSITION_COLUMN_INDEX]);
//                String temp = new String(tokens[count][HAPMAP_CHROMOSOME_COLUMN_INDEX]);
//                if (currLocus == null) {
//                    lociOffsets.add(site);
//                    currLocus = temp;
//                    minPosition = position;
//                    prevPosition = -1;
//                } else if (!temp.equals(currLocus)) {
//                    loci.add(new Chromosome(currLocus, currLocus, minPosition, prevPosition, null, null));
//                    lociOffsets.add(site);
//                    currLocus = temp;
//                    minPosition = position;
//                    prevPosition = -1;
//                }
//
//                if (position < prevPosition) {
//                    throw new IllegalStateException("ImportUtils: readFromHapmap: Sites are not properly sorted for chromosome: " + currLocus + " at " + position + " and " + prevPosition);
//                }
//
//                count++;
//
//                if (count == numSitesToProcess) {
//                    pool.execute(ProcessLineOfHapmap.getInstance(alleles, theData, TasselPrefs.getAlignmentRetainRareAlleles(), tokens, count, currentSite, numTaxa, lineInFile, isSBit));
//                    count = 0;
//                    currentSite += numSitesToProcess;
//                    tokens = new String[numSitesToProcess][];
//                }
//
//                physicalPositions[site] = position;
//                prevPosition = position;
//
//                if (listener != null) {
//                    listener.progress((int) (((double) (site + 1) / (double) numSites) * 100.0), null);
//                }
//            }
//
//            if (count != 0) {
//                pool.execute(ProcessLineOfHapmap.getInstance(alleles, theData, TasselPrefs.getAlignmentRetainRareAlleles(), tokens, count, currentSite, numTaxa, lineInFile, isSBit));
//            }
//
//
//            pool.shutdown();
//            if (!pool.awaitTermination(6000, TimeUnit.SECONDS)) {
//                throw new IllegalStateException("ImportUtils: readFromHapmap: processing threads timed out.");
//            }
//
//            if (currLocus != null) {
//                loci.add(new Chromosome(currLocus, currLocus, minPosition, prevPosition, null, null));
//            }
//
//            prevTime = currentTime;
//            currentTime = System.currentTimeMillis();
//            myLogger.info("readFromHapmap: Time to read file: " + ((currentTime - prevTime) / 1000));
//
//            String[] taxaNames = new String[numTaxa];
//            System.arraycopy(header, NUM_HAPMAP_NON_TAXA_HEADERS, taxaNames, 0, numTaxa);
//            IdGroup idGroup = new SimpleIdGroup(taxaNames);
//
//            Chromosome[] lociFinal = new Chromosome[loci.size()];
//            loci.toArray(lociFinal);
//            int[] offsetsFinal = new int[lociOffsets.size()];
//            for (int i = 0; i < lociOffsets.size(); i++) {
//                offsetsFinal[i] = ((Integer) lociOffsets.get(i)).intValue();
//            }
//
//            Alignment result = BitAlignment.getNucleotideInstance(idGroup, alleles, theData, null, null, physicalPositions, TasselPrefs.getAlignmentMaxAllelesToRetain(), lociFinal, offsetsFinal, snpIDs, TasselPrefs.getAlignmentRetainRareAlleles(), isSBit);
//
//            prevTime = currentTime;
//            currentTime = System.currentTimeMillis();
//            myLogger.info("readFromHapmap: Time to create Alignment: " + ((currentTime - prevTime) / 1000));
//
//            return result;
//        } catch (Exception e) {
//            e.printStackTrace();
//            throw new IllegalArgumentException("ImportUtils: readFromHapmap: Problem creating Alignment: " + pedFilename + ": " + ExceptionUtils.getExceptionCauses(e));
//        } finally {
//            try {
//                fileIn.close();
//            } catch (Exception ex) {
//                // do nothing
//            }
//        }
//
//    }
    
    public static GenotypeTable readFasta(String filename) throws FileNotFoundException, IOException {

        BufferedReader reader = Utils.getBufferedReader(filename);

        List<String> taxa = new ArrayList<>();
        List<String> sequences = new ArrayList<>();

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
        TaxaList taxaList = (new TaxaListBuilder()).addAll(taxaNames).build();

        String[] sequenceArray = new String[sequences.size()];
        sequences.toArray(sequenceArray);

        GenotypeCallTable genotype = GenotypeCallTableBuilder.getUnphasedNucleotideGenotypeBuilder(sequences.size(), sequenceLength)
                .setBases(sequenceArray)
                .build();

        return GenotypeTableBuilder.getInstance(genotype, PositionListBuilder.getInstance(sequenceLength), taxaList);

    }
}
