/*
 * ProcessLineOfHapmap
 */
package net.maizegenetics.pal.alignment;

import java.util.LinkedList;
import java.util.Queue;
import net.maizegenetics.util.OpenBitSet;

/**
 *
 * @author terry
 */
public class ProcessLineOfHapmap implements Runnable {

    private static int myNumInstances = 0;
    private static Queue myInstances = new LinkedList();
    private String[] myTokens;
    private int mySite;
    private int myNumTaxa;
    private int myLineInFile;
    private boolean myComplete = false;
    private OpenBitSet[][] myData;
    private byte[][] myAlleles;
    private int myNumAlleles;
    private int myMaxNumAlleles;
    private boolean myRetainRareAlleles;
    private boolean myIsSBit;

    private ProcessLineOfHapmap(byte[][] alleles, OpenBitSet[][] data, boolean retainRareAlleles, String[] tokens, int site, int numTaxa, int lineInFile, boolean isSBit) {
        setVariables(alleles, data, retainRareAlleles, tokens, site, numTaxa, lineInFile, isSBit);
    }

    public static ProcessLineOfHapmap getInstance(byte[][] alleles, OpenBitSet[][] data, boolean retainRareAlleles, String[] tokens, int site, int numTaxa, int lineInFile, boolean isSBit) {
        if (myNumInstances <= 35) {
            return new ProcessLineOfHapmap(alleles, data, retainRareAlleles, tokens, site, numTaxa, lineInFile, isSBit);
        } else {
            ProcessLineOfHapmap result;
            while ((result = (ProcessLineOfHapmap) myInstances.poll()) == null) {
            }
            result.setVariables(alleles, data, retainRareAlleles, tokens, site, numTaxa, lineInFile, isSBit);
            return result;
        }
    }

    private void setVariables(byte[][] alleles, OpenBitSet[][] data, boolean retainRareAlleles, String[] tokens, int site, int numTaxa, int lineInFile, boolean isSBit) {
        myData = data;
        myTokens = tokens;
        mySite = site;
        myNumTaxa = numTaxa;
        myLineInFile = lineInFile;
        myComplete = false;
        myAlleles = alleles;
        myNumAlleles = myAlleles[0].length;
        myRetainRareAlleles = retainRareAlleles;
        myMaxNumAlleles = myData.length;
        if (myRetainRareAlleles) {
            myMaxNumAlleles--;
        }
        myIsSBit = isSBit;
    }

    private void clearVariables() {
        myData = null;
        myTokens = null;
        myAlleles = null;
    }

    public void run() {
        try {
            if (myComplete) {
                throw new IllegalStateException("ImportUtils: ProcessLineOfHapmap: run: trying to run completed instance.");
            }
            myComplete = true;

            byte[] data = new byte[myNumTaxa];
            for (int i = 0; i < myNumTaxa; i++) {
                try {
                    data[i] = NucleotideAlignmentConstants.getNucleotideDiploidByte(myTokens[ImportUtils.NUM_HAPMAP_NON_TAXA_HEADERS + i]);
                } catch (IndexOutOfBoundsException ex) {
                    throw new IllegalStateException("Number of Taxa: " + myNumTaxa + " does not match number of values at line in file: " + myLineInFile + " site: " + mySite);
                } catch (Exception e) {
                    throw new IllegalStateException("Line in File: " + myLineInFile, e);
                }
            }

            setAlleles(data);
            setBits(data);

            clearVariables();
            myInstances.offer(this);
        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    private void setAlleles(byte[] data) {
        int[][] alleles = AlignmentUtils.getAllelesSortedByFrequency(data);
        int resultSize = alleles[0].length;
        for (int i = 0; i < myNumAlleles; i++) {
            myAlleles[mySite][i] = (i < resultSize) ? (byte) alleles[0][i] : Alignment.UNKNOWN_ALLELE;
        }
    }

    private void setBits(byte[] data) {
        byte[] cb = new byte[2];
        for (int t = 0; t < myNumTaxa; t++) {
            cb[0] = (byte) ((data[t] >>> 4) & 0xf);
            cb[1] = (byte) (data[t] & 0xf);
            for (int i = 0; i < 2; i++) {
                if (cb[i] != Alignment.UNKNOWN_ALLELE) {
                    boolean isRare = true;
                    for (int j = 0; j < myMaxNumAlleles; j++) {
                        if (cb[i] == myAlleles[mySite][j]) {
                            if (myIsSBit) {
                                myData[j][mySite].fastSet(t);
                            } else {
                                myData[j][t].fastSet(mySite);
                            }
                            isRare = false;
                            break;
                        }
                    }
                    if (isRare && myRetainRareAlleles) {
                        if (myIsSBit) {
                            myData[myMaxNumAlleles][mySite].fastSet(t);
                        } else {
                            myData[myMaxNumAlleles][t].fastSet(mySite);
                        }
                    }
                }
            }
        }
    }
}
