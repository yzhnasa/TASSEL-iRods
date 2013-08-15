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
}
