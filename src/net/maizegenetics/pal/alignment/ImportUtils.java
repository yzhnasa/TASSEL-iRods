/*
 * ImportUtils
 */
package net.maizegenetics.pal.alignment;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;

import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.SimpleIdGroup;

import net.maizegenetics.util.ExceptionUtils;
import net.maizegenetics.util.Utils;

import org.apache.log4j.Logger;

/**
 * The class imports Alignment from
 * various file formats.
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

    private ImportUtils() {
        // Utility Class - do not instantiate.
    }

    /**
     * This read a hapmap formatted file and creates an Alignment.
     *
     * @param filename file name
     * @param chrom chromosome to get.  null to get everything.
     *
     * @return Alignment
     */
    public static Alignment readFromHapmap(String filename, String chrom) {

        String[][] chromsAvailable = getFileInfo(filename, HAPMAP_CHROMOSOME_COLUMN_INDEX, HAPMAP_SNPID_COLUMN_INDEX);

        int numSites = 0;
        int depthInFile = 1;
        if (chrom == null) {
            for (int i = 1; i < chromsAvailable.length; i++) {
                numSites += Integer.parseInt(chromsAvailable[i][1]);
            }
            depthInFile = 1;
        } else {
            int[] chromInfo = getChromInFileInfo(chromsAvailable, chrom);
            numSites = chromInfo[0];
            if (numSites == 0) {
                throw new IllegalStateException("Desired chromosome not available in " + filename + "."); //chrom not available
            }
            depthInFile = chromInfo[1]; //how many lines into .txt, data for desired chromosome begins
        }

        int minPosition = Integer.MAX_VALUE;
        String currLocus = null;
        List<Locus> loci = new ArrayList<Locus>();
        List<Integer> lociOffsets = new ArrayList<Integer>();

        try {

            BufferedReader fileIn = Utils.getBufferedReader(filename, 1000000);
            String[] header = WHITESPACE_PATTERN.split(fileIn.readLine());
            int lineInFile = 1;
            int numTaxa = header.length - NUM_HAPMAP_NON_TAXA_HEADERS;
            String[] snpIDs = new String[numSites];
            int prevPosition = -1;
            for (int depth = 0; depth < depthInFile - 1; depth++) {
                fileIn.readLine(); //shift .txt reader to desired position in .txt file
                lineInFile++;
            }
            byte[][] theData = new byte[numTaxa][numSites];
            int[] physicalPositions = new int[numSites];
            for (int site = 0; site < numSites; site++) {
                String[] s = WHITESPACE_PATTERN.split(fileIn.readLine());
                lineInFile++;

                int position = Integer.parseInt(s[HAPMAP_POSITION_COLUMN_INDEX]);

                if (currLocus == null) {
                    lociOffsets.add(site);
                    currLocus = s[HAPMAP_CHROMOSOME_COLUMN_INDEX];
                    minPosition = position;
                    prevPosition = -1;
                } else if (!s[HAPMAP_CHROMOSOME_COLUMN_INDEX].equals(currLocus)) {
                    loci.add(new Locus(currLocus, currLocus, minPosition, prevPosition, null, null));
                    lociOffsets.add(site);
                    currLocus = s[HAPMAP_CHROMOSOME_COLUMN_INDEX];
                    minPosition = position;
                    prevPosition = -1;
                }

                if (position < prevPosition) {
                    throw new IllegalStateException("Sites are not properly sorted for chromosome: " + currLocus + " at " + position + " and " + prevPosition);
                }

                if (numTaxa + NUM_HAPMAP_NON_TAXA_HEADERS != s.length) {
                    throw new IllegalStateException("Number of Taxa: " + numTaxa + " does not match number of values: " + (s.length - NUM_HAPMAP_NON_TAXA_HEADERS) + " at line in file: " + lineInFile + " site: " + site);
                }
                for (int i = 0; i < numTaxa; i++) {
                    theData[i][site] = NucleotideAlignmentConstants.getNucleotideDiploidByte(s[NUM_HAPMAP_NON_TAXA_HEADERS + i]);
                }
                snpIDs[site] = s[HAPMAP_SNPID_COLUMN_INDEX];
                physicalPositions[site] = position;
                prevPosition = position;
            }

            if (currLocus != null) {
                loci.add(new Locus(currLocus, currLocus, minPosition, prevPosition, null, null));
            }

            String[] taxaNames = new String[numTaxa];
            System.arraycopy(header, NUM_HAPMAP_NON_TAXA_HEADERS, taxaNames, 0, numTaxa);
            IdGroup idGroup = new SimpleIdGroup(taxaNames);

            Locus[] lociFinal = new Locus[loci.size()];
            loci.toArray(lociFinal);
            int[] offsetsFinal = new int[lociOffsets.size()];
            for (int i = 0; i < lociOffsets.size(); i++) {
                offsetsFinal[i] = ((Integer) lociOffsets.get(i)).intValue();
            }

            return SBitAlignment.getNucleotideInstance(idGroup, theData, null, null, physicalPositions, Alignment.DEFAULT_MAX_NUM_ALLELES, lociFinal, offsetsFinal, snpIDs, true);
        } catch (Exception e) {
            e.printStackTrace();
            throw new IllegalArgumentException("Problem creating Alignment: " + filename + ": " + ExceptionUtils.getExceptionCauses(e));
        }

    }

    public static Alignment readFromHapmap(String filename) {

        int minPosition = Integer.MAX_VALUE;
        String currLocus = null;
        List<Locus> loci = new ArrayList<Locus>();
        List<Integer> lociOffsets = new ArrayList<Integer>();

        long currentTime = System.currentTimeMillis();
        int numSites = -1;
        BufferedReader reader = null;
        try {
            reader = Utils.getBufferedReader(filename, 1000000);
            while (reader.readLine() != null) {
                numSites++;
            }
        } catch (Exception e) {
            e.printStackTrace();
            throw new IllegalArgumentException("Problem creating Alignment: " + filename + ": " + ExceptionUtils.getExceptionCauses(e));
        } finally {
            try {
                reader.close();
            } catch (Exception ex) {
                // do nothing
            }
        }
        long prevTime = currentTime;
        currentTime = System.currentTimeMillis();
        System.out.println("Time to count lines: " + ((currentTime - prevTime) / 1000));


        BufferedReader fileIn = null;
        try {
            ExecutorService pool = Executors.newFixedThreadPool(10);

            fileIn = Utils.getBufferedReader(filename, 1000000);
            String[] header = WHITESPACE_PATTERN.split(fileIn.readLine());
            int lineInFile = 1;
            int numTaxa = header.length - NUM_HAPMAP_NON_TAXA_HEADERS;
            String[] snpIDs = new String[numSites];
            int prevPosition = -1;
            byte[][] theData = new byte[numTaxa][numSites];
            int[] physicalPositions = new int[numSites];
            for (int site = 0; site < numSites; site++) {

                lineInFile++;

                String input = fileIn.readLine();
                Matcher matcher = WHITESPACE_PATTERN.matcher(input);


                int position = 0;

                int index = 0;
                int columnIndex = -1;
                while (matcher.find()) {
                    columnIndex++;

                    if (columnIndex == HAPMAP_SNPID_COLUMN_INDEX) {
                        snpIDs[site] = input.substring(index, matcher.start());
                    } else if (columnIndex == HAPMAP_POSITION_COLUMN_INDEX) {
                        position = Integer.parseInt(input.substring(index, matcher.start()));
                    } else if (columnIndex == HAPMAP_CHROMOSOME_COLUMN_INDEX) {
                        String temp = input.substring(index, matcher.start());
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
                    }

                    index = matcher.end();

                    if (columnIndex == (NUM_HAPMAP_NON_TAXA_HEADERS - 1)) {
                        break;
                    }

                }

                if (position < prevPosition) {
                    throw new IllegalStateException("Sites are not properly sorted for chromosome: " + currLocus + " at " + position + " and " + prevPosition);
                }

                pool.execute(new ProcessLine(theData, matcher, input, site, index, numTaxa, lineInFile));

                physicalPositions[site] = position;
                prevPosition = position;
            }

            pool.shutdown();
            if (!pool.awaitTermination(120, TimeUnit.SECONDS)) {
                throw new IllegalStateException("ImportUtils: readFromHapmap: processing threads timed out.");
            }

            if (currLocus != null) {
                loci.add(new Locus(currLocus, currLocus, minPosition, prevPosition, null, null));
            }

            prevTime = currentTime;
            currentTime = System.currentTimeMillis();
            System.out.println("Time to read file: " + ((currentTime - prevTime) / 1000));

            String[] taxaNames = new String[numTaxa];
            System.arraycopy(header, NUM_HAPMAP_NON_TAXA_HEADERS, taxaNames, 0, numTaxa);
            IdGroup idGroup = new SimpleIdGroup(taxaNames);

            Locus[] lociFinal = new Locus[loci.size()];
            loci.toArray(lociFinal);
            int[] offsetsFinal = new int[lociOffsets.size()];
            for (int i = 0; i < lociOffsets.size(); i++) {
                offsetsFinal[i] = ((Integer) lociOffsets.get(i)).intValue();
            }

            Alignment result = SBitAlignment.getNucleotideInstance(idGroup, theData, null, null, physicalPositions, Alignment.DEFAULT_MAX_NUM_ALLELES, lociFinal, offsetsFinal, snpIDs, true);

            prevTime = currentTime;
            currentTime = System.currentTimeMillis();
            System.out.println("Time to create Alignment: " + ((currentTime - prevTime) / 1000));

            return result;
        } catch (Exception e) {
            e.printStackTrace();
            throw new IllegalArgumentException("Problem creating Alignment: " + filename + ": " + ExceptionUtils.getExceptionCauses(e));
        } finally {
            try {
                fileIn.close();
            } catch (Exception ex) {
                // do nothing
            }
        }

    }

    private static class ProcessLine implements Runnable {

        private byte[][] myData;
        private Matcher myMatcher;
        private String myLine;
        private int mySite;
        private int myIndex;
        private int myNumTaxa;
        private int myLineInFile;

        public ProcessLine(byte[][] data, Matcher matcher, String line, int site, int index, int numTaxa, int lineInFile) {
            myData = data;
            myMatcher = matcher;
            myLine = line;
            mySite = site;
            myIndex = index;
            myNumTaxa = numTaxa;
            myLineInFile = lineInFile;
        }

        public void run() {
            int columnIndex = -1;
            while (myMatcher.find()) {
                columnIndex++;
                try {
                    myData[columnIndex][mySite] = NucleotideAlignmentConstants.getNucleotideDiploidByte(myLine.substring(myIndex, myMatcher.start()));
                } catch (IndexOutOfBoundsException ex) {
                    throw new IllegalStateException("Number of Taxa: " + myNumTaxa + " does not match number of values at line in file: " + myLineInFile + " site: " + mySite);
                }
                myIndex = myMatcher.end();
            }
            columnIndex++;
            myData[columnIndex][mySite] = NucleotideAlignmentConstants.getNucleotideDiploidByte(myLine.substring(myIndex));
        }
    }

    /**
     * Get chromosome IDs for all represented chromosomes in a given .map file
     * along with number of sites on each chromosome
     * currently unique to Plink files
     *
     * @param infileName
     * @param chromIndex, column number where chromosome ID are stored, varies based on format
     * Hapmap index = 2
     * Plink index = 0
     * FlapJack index = 1
     * @param snpIDIndex, column number where SNP ID are stored, varies based on format
     * Hapmap index = 0
     * Plink index = 1
     * Flapjack index = 0
     *
     * @return stringArray
     */
    public static String[][] getFileInfo(String infileName, int chromIndex, int snpIDIndex) {
        try {
            ArrayList<String> stringList = new ArrayList<String>();
            ArrayList<Integer> countList = new ArrayList<Integer>();
            ArrayList<Integer> idLengthList = new ArrayList<Integer>();
            int count = 0;
            String previousChrom = "";
            int SNPidLength = 0;
            BufferedReader fileIn = Utils.getBufferedReader(infileName, 1000000);
            String fileInLine = fileIn.readLine();
            while (fileInLine != null) {
                String[] fileLine = WHITESPACE_PATTERN.split(fileInLine);
                String currentChrom = fileLine[chromIndex].trim();
                String SNPID = fileLine[snpIDIndex].trim();
                if (!currentChrom.equals(previousChrom)) {
                    if (stringList.size() != 0) {
                        countList.add(new Integer(count));
                        idLengthList.add(new Integer(SNPidLength));
                        count = 0;
                        SNPidLength = 0;
                    }
                    stringList.add(currentChrom);
                }
                count++;
                if (SNPID.length() > SNPidLength) {
                    SNPidLength = SNPID.length();
                }
                previousChrom = currentChrom;
                fileInLine = fileIn.readLine();
            }
            countList.add(new Integer(count));
            idLengthList.add(new Integer(SNPidLength));
            String[][] stringArray = new String[stringList.size()][3];
            for (int i = 0; i < stringList.size(); i++) {
                stringArray[i][0] = stringList.get(i);
                stringArray[i][1] = countList.get(i).toString();
                stringArray[i][2] = idLengthList.get(i).toString();
            }
            return stringArray;
        } catch (Exception e) {
            System.err.println("File IO in getFileInfo: " + e);
        }
        return null;
    }

    /**
     * Returns the number of sites on a given chromosome as well as how far into a .map file that the
     * sites for a given chromosome begins, works in conjunction with getChromsAvailableCounts
     *
     * @param chromCounts
     * @param chrom
     *
     * @return chromInfo[0] = numSites, chromInfo[1] = depth in .map file the sites begin
     * chromInfo[2] = length of longest SNP id within the chromosome
     */
    public static int[] getChromInFileInfo(String[][] chromCounts, String chrom) {
        int[] chromInfo = {0, 0, 0};
        for (int i = 0; i < chromCounts.length; i++) {
            if (chromCounts[i][0].equals(chrom)) {
                chromInfo[0] = Integer.parseInt(chromCounts[i][1]);
                chromInfo[2] = Integer.parseInt(chromCounts[i][2]);
                break;
            }
            chromInfo[1] += Integer.parseInt(chromCounts[i][1]);
        }
        return chromInfo;
    }

    public static Alignment readFasta(String filename) throws FileNotFoundException, IOException {

        BufferedReader reader = Utils.getBufferedReader(filename);

        List taxa = new ArrayList();
        List sequences = new ArrayList();

        String line = null;
        line = reader.readLine();
        boolean sequence = false;
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
                sequences.add(builder.toString());
                sequence = false;
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
        return SBitAlignment.getNucleotideInstance(idGroup, sequenceArray, null, null, null, Alignment.DEFAULT_MAX_NUM_ALLELES, new Locus[]{unknown}, new int[]{0}, null, true);

    }
}
