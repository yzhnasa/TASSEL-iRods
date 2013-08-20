/*
 *  ByteGenotype
 */
package net.maizegenetics.pal.alignment.genotype;

import org.apache.log4j.Logger;

/**
 *
 * @author Terry Casstevens
 */
public class ByteGenotype extends AbstractGenotype {

    private static final Logger myLogger = Logger.getLogger(ByteGenotype.class);
    private final byte[][] myGenotype;
    private byte[][] myGenotypeSiteTaxa = null;

    ByteGenotype(byte[][] genotype, boolean phased, String[][] alleleEncodings) {
        super(genotype.length, genotype[0].length, phased, alleleEncodings);
        myGenotype = new byte[myTaxaCount][mySiteCount];
        for (int t = 0; t < myTaxaCount; t++) {
            System.arraycopy(genotype[t], 0, myGenotype[t], 0, mySiteCount);
        }
    }

    @Override
    public byte getBase(int taxon, int site) {
        return myGenotype[taxon][site];
    }

    @Override
    byte[][] getBasesSiteTaxa() {
        if (myGenotypeSiteTaxa == null) {
            initTaxaSite();
        }
        return myGenotypeSiteTaxa;
    }

    @Override
    byte[][] getBasesTaxaSite() {
        return myGenotype;
    }

    private void initTaxaSite() {
        myGenotypeSiteTaxa = new byte[mySiteCount][myTaxaCount];
        for (int bigS = 0; bigS < mySiteCount; bigS += 64) {
            int length = (mySiteCount - bigS < 64) ? mySiteCount - bigS : 64;
            for (int t = 0; t < myTaxaCount; t++) {
                for (int s = 0; s < length; s++) {
                    myGenotypeSiteTaxa[bigS + s][t] = myGenotype[t][bigS + s];
                }
            }
        }
    }
}
