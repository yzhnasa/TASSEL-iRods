/*
 *  GenotypeBuilder
 */
package net.maizegenetics.pal.alignment.genotype;

import net.maizegenetics.pal.alignment.NucleotideAlignmentConstants;
import net.maizegenetics.util.SuperByteMatrix;
import net.maizegenetics.util.SuperByteMatrixBuilder;

/**
 *
 * @author Terry Casstevens
 */
public class GenotypeBuilder {

    private SuperByteMatrix myGenotype;
    private final boolean myIsPhased;
    private final String[][] myAlleleEncodings;

    private GenotypeBuilder(int numTaxa, int numSites, boolean phased, String[][] alleleEncodings) {
        myGenotype = SuperByteMatrixBuilder.getInstance(numTaxa, numSites);
        myIsPhased = phased;
        myAlleleEncodings = alleleEncodings;
    }

    public static GenotypeBuilder getUnphasedNucleotideGenotypeBuilder(int numTaxa, int numSites) {
        return new GenotypeBuilder(numTaxa, numSites, false, NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES);
    }

    public void setBase(int taxon, int site, byte value) {
        myGenotype.set(taxon, site, value);
    }

    public void setBaseRangeForTaxon(int taxon, int startSite, byte[] value) {
        for (int i = 0; i < value.length; i++) {
            myGenotype.set(taxon, i, value[i]);
        }
    }

    public Genotype build() {
        SuperByteMatrix temp = myGenotype;
        myGenotype = null;
        return new ByteGenotype(temp, myIsPhased, myAlleleEncodings);
    }

    public Genotype buildHDF5(String filename) {
        SuperByteMatrix temp = myGenotype;
        myGenotype = null;
        throw new UnsupportedOperationException();
    }
}
