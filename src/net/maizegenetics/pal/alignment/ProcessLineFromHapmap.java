/*
 * ProcessLineFromHapmap
 */
package net.maizegenetics.pal.alignment;

import java.util.LinkedList;
import java.util.Queue;

/**
 *
 * @author terry
 */
public class ProcessLineFromHapmap implements Runnable {

    private static int myNumInstances = 0;
    private static Queue myInstances = new LinkedList();
    private byte[][] myData;
    private String[] myTokens;
    private int mySite;
    private int myNumTaxa;
    private int myLineInFile;
    private boolean myComplete = false;

    private ProcessLineFromHapmap(byte[][] data, String[] tokens, int site, int numTaxa, int lineInFile) {
        setVariables(data, tokens, site, numTaxa, lineInFile);
    }

    public static ProcessLineFromHapmap getInstance(byte[][] data, String[] tokens, int site, int numTaxa, int lineInFile) {
        if (myNumInstances <= 35) {
            return new ProcessLineFromHapmap(data, tokens, site, numTaxa, lineInFile);
        } else {
            ProcessLineFromHapmap result;
            while ((result = (ProcessLineFromHapmap) myInstances.poll()) == null) {
            }
            result.setVariables(data, tokens, site, numTaxa, lineInFile);
            return result;
        }
    }

    private void setVariables(byte[][] data, String[] tokens, int site, int numTaxa, int lineInFile) {
        myData = data;
        myTokens = tokens;
        mySite = site;
        myNumTaxa = numTaxa;
        myLineInFile = lineInFile;
        myComplete = false;
    }

    private void clearVariables() {
        myData = null;
        myTokens = null;
    }

    public void run() {
        if (myComplete) {
            throw new IllegalStateException("ImportUtils: ProcessLineFromHapmap: run: trying to run completed instance.");
        }
        myComplete = true;

        for (int i = 0; i < myNumTaxa; i++) {
            try {
                myData[i][mySite] = NucleotideAlignmentConstants.getNucleotideDiploidByte(myTokens[ImportUtils.NUM_HAPMAP_NON_TAXA_HEADERS + i]);
            } catch (IndexOutOfBoundsException ex) {
                throw new IllegalStateException("Number of Taxa: " + myNumTaxa + " does not match number of values at line in file: " + myLineInFile + " site: " + mySite);
            } catch (Exception e) {
                throw new IllegalStateException("Line in File: " + myLineInFile, e);
            }
        }

        clearVariables();
        myInstances.offer(this);

    }
}
