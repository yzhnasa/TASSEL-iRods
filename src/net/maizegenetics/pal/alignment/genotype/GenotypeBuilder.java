/*
 *  GenotypeBuilder
 */
package net.maizegenetics.pal.alignment.genotype;

import net.maizegenetics.pal.alignment.NucleotideAlignmentConstants;

/**
 *
 * @author Terry Casstevens
 */
public class GenotypeBuilder {

    private final byte[][] myGenotype;
    private final boolean myIsPhased;
    private final String[][] myAlleleEncodings;

    private GenotypeBuilder(int numTaxa, int numSites, boolean phased, String[][] alleleEncodings) {
        myGenotype = new byte[numTaxa][numSites];
        myIsPhased = phased;
        myAlleleEncodings = alleleEncodings;
    }

    public static GenotypeBuilder getUnphasedNucleotideGenotypeBuilder(int numTaxa, int numSites) {
        return new GenotypeBuilder(numTaxa, numSites, false, NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES);
    }

    public void setBase(int taxon, int site, byte value) {
        myGenotype[taxon][site] = value;
    }

    public Genotype build() {
        return new ByteGenotype(myGenotype, myIsPhased, myAlleleEncodings);
    }

    public Genotype buildHDF5(String filename) {
        throw new UnsupportedOperationException();
    }
}
