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
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import java.util.regex.Pattern;

import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.SimpleIdGroup;

import net.maizegenetics.util.ExceptionUtils;
import net.maizegenetics.util.ProgressListener;
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

    public static Alignment readFromHapmap(final String filename, ProgressListener listener) {

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
            int numThreads = 30;
            int currentFuture = 0;
            Future[] futures = new Future[numThreads];
            ExecutorService pool = Executors.newFixedThreadPool(numThreads);

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
                String[] tokens = WHITESPACE_PATTERN.split(input);

                int position = 0;

                snpIDs[site] = new String(tokens[HAPMAP_SNPID_COLUMN_INDEX]);
                position = Integer.parseInt(tokens[HAPMAP_POSITION_COLUMN_INDEX]);
                String temp = new String(tokens[HAPMAP_CHROMOSOME_COLUMN_INDEX]);
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
                    throw new IllegalStateException("Sites are not properly sorted for chromosome: " + currLocus + " at " + position + " and " + prevPosition);
                }

                futures[currentFuture++] = pool.submit(ProcessLineFromHapmap.getInstance(theData, tokens, site, numTaxa, lineInFile));
                if (currentFuture == numThreads) {
                    for (int i = 0; i < numThreads; i++) {
                        futures[i].get();
                    }
                    currentFuture = 0;
                }

                physicalPositions[site] = position;
                prevPosition = position;

                if (listener != null) {
                    listener.progress((int) (((double) (site + 1) / (double) numSites) * 80.0), null);
                }
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
