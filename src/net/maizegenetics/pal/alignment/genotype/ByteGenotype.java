/*
 *  ByteGenotype
 */
package net.maizegenetics.pal.alignment.genotype;

import net.maizegenetics.util.SuperByteMatrix;
import net.maizegenetics.util.SuperByteMatrixBuilder;
import org.apache.log4j.Logger;

/**
 *
 * @author Terry Casstevens
 */
public class ByteGenotype extends AbstractGenotype {

    private static final Logger myLogger = Logger.getLogger(ByteGenotype.class);
    private final SuperByteMatrix myGenotype;

    ByteGenotype(byte[][] genotype, boolean phased, String[][] alleleEncodings) {
        super(genotype.length, genotype[0].length, phased, alleleEncodings);
        myGenotype = SuperByteMatrixBuilder.getInstance(myTaxaCount, mySiteCount);
        System.out.println("myGenotype class: " + myGenotype.getClass().getName());
        for (int t = 0; t < myTaxaCount; t++) {
            for (int s = 0; s < mySiteCount; s++) {
                myGenotype.set(t, s, genotype[t][s]);
            }
        }
    }

    ByteGenotype(SuperByteMatrix genotype, boolean phased, String[][] alleleEncodings) {
        super(genotype.getNumRows(), genotype.getNumColumns(), phased, alleleEncodings);
        myGenotype = genotype;
        System.out.println("myGenotype class: " + myGenotype.getClass().getName());
    }

    @Override
    public byte getBase(int taxon, int site) {
        return myGenotype.get(taxon, site);
    }

    @Override
    public byte[] getGenotypeForAllSites(int taxon) {
        return myGenotype.getAllColumns(taxon);
    }

    @Override
    public byte[] getGenotypeForSiteRange(int taxon, int start, int end) {
        return myGenotype.getColumnRange(taxon, start, end);
    }

    @Override
    public byte[] getGenotypeForAllTaxa(int site) {
        return myGenotype.getAllRows(site);
    }
}
