/*
 *  GenotypeBuilder
 */
package net.maizegenetics.dna.snp.genotypecall;

import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Reader;
import net.maizegenetics.dna.snp.Alignment;
import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.util.SuperByteMatrix;
import net.maizegenetics.util.SuperByteMatrixBuilder;

import java.util.regex.Pattern;

/**
 *
 * @author Terry Casstevens
 */
public class GenotypeBuilder {

    private SuperByteMatrix myGenotype;
    private boolean myIsPhased = false;
    private String[][] myAlleleEncodings = NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES;

    private GenotypeBuilder(SuperByteMatrix genotype) {
        myGenotype = genotype;
    }

    /**
     * Get Genotype Builder given number of taxa and sites. Performance
     * optimized for site loop inside taxon loop. Default is unphased and
     * NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES encoding.
     *
     * @param numTaxa number of taxa
     * @param numSites number of sites.
     *
     * @return Genotype Builder
     */
    public static GenotypeBuilder getInstance(int numTaxa, int numSites) {
        return getUnphasedNucleotideGenotypeBuilder(numTaxa, numSites);
    }

    /**
     * Get Genotype Builder given number of taxa and sites. Performance
     * optimized for taxon loop inside site loop. Default is unphased and
     * NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES encoding.
     *
     * @param numTaxa number of taxa
     * @param numSites number of sites.
     *
     * @return Genotype Builder
     */
    public static GenotypeBuilder getInstanceTranspose(int numTaxa, int numSites) {
        SuperByteMatrix matrix = SuperByteMatrixBuilder.getInstanceTranspose(numTaxa, numSites);
        matrix.setAll(Alignment.UNKNOWN_DIPLOID_ALLELE);
        return new GenotypeBuilder(matrix);
    }

    /**
     * Get Genotype Builder given number of taxa and sites. Performance
     * optimized for site loop inside taxon loop. Default is unphased and
     * NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES encoding.
     *
     * @param numTaxa number of taxa
     * @param numSites number of sites.
     *
     * @return Genotype Builder
     */
    public static GenotypeBuilder getUnphasedNucleotideGenotypeBuilder(int numTaxa, int numSites) {
        SuperByteMatrix matrix = SuperByteMatrixBuilder.getInstance(numTaxa, numSites);
        matrix.setAll(Alignment.UNKNOWN_DIPLOID_ALLELE);
        return new GenotypeBuilder(matrix);
    }

    public static Genotype getFilteredInstance(Genotype genotype, int numTaxa, int[] taxaRedirect, int numSites, int rangeStart, int rangeEnd) {
        return new FilterGenotype(genotype, numTaxa, taxaRedirect, numSites, rangeStart, rangeEnd);
    }

    public static Genotype getFilteredInstance(Genotype genotype, int numTaxa, int[] taxaRedirect, int numSites, int[] siteRedirect) {
        return new FilterGenotype(genotype, numTaxa, taxaRedirect, numSites, siteRedirect);
    }

    public GenotypeBuilder setBase(int taxon, int site, byte value) {
        myGenotype.set(taxon, site, value);
        return this;
    }

    public GenotypeBuilder setBaseRangeForTaxon(int taxon, int startSite, byte[] value) {
        //TODO this needs an array copy method, startSite was eliminated
        for (int i = 0; i < value.length; i++) {
            myGenotype.set(taxon, i + startSite, value[i]);
        }
        return this;
    }

    public GenotypeBuilder setBases(String[] data) {

        int numTaxa = data.length;

        int numSites = data[0].length();

        for (int site = 0; site < numSites; site++) {
            for (int taxon = 0; taxon < numTaxa; taxon++) {
                setBase(taxon, site, NucleotideAlignmentConstants.getNucleotideDiploidByte(data[taxon].charAt(site)));
            }
        }

        return this;

    }

    public GenotypeBuilder setBases(String[][] data) {

        if ((data == null) || (data.length == 0)) {
            throw new IllegalArgumentException("BitAlignment: getInstance: data can not be empty.");
        }

        int numTaxa = data.length;
        int numSites = data[0].length;

        for (int site = 0; site < numSites; site++) {
            if (data[0][0].contains(":")) {
                Pattern colon = Pattern.compile(":");
                for (int taxon = 0; taxon < numTaxa; taxon++) {
                    if (data[taxon][site].equalsIgnoreCase(Alignment.UNKNOWN_DIPLOID_ALLELE_STR)) {
                        setBase(taxon, site, Alignment.UNKNOWN_DIPLOID_ALLELE);
                    } else if (data[taxon][site].equals("?") || data[taxon][site].equals("?:?")) {
                        setBase(taxon, site, Alignment.UNKNOWN_DIPLOID_ALLELE);
                    } else {
                        String[] siteval = colon.split(data[taxon][site]);
                        byte first = NucleotideAlignmentConstants.getNucleotideAlleleByte(siteval[0]);
                        byte second = NucleotideAlignmentConstants.getNucleotideAlleleByte(siteval[1]);
                        setBase(taxon, site, (byte) ((first << 4) | second));
                    }
                }
            } else {
                for (int taxon = 0; taxon < numTaxa; taxon++) {
                    if (data[taxon][site].equalsIgnoreCase(Alignment.UNKNOWN_ALLELE_STR)) {
                        setBase(taxon, site, Alignment.UNKNOWN_DIPLOID_ALLELE);
                    } else if (data[taxon][site].equals("?")) {
                        setBase(taxon, site, Alignment.UNKNOWN_DIPLOID_ALLELE);
                    } else {
                        setBase(taxon, site, NucleotideAlignmentConstants.getNucleotideAlleleByte(data[taxon][site]));
                    }
                }
            }
        }
        return this;
    }

    public GenotypeBuilder isPhased(boolean isPhased) {
        myIsPhased = isPhased;
        return this;
    }

    public GenotypeBuilder alleleEncodings(String[][] alleleEncodings) {
        myAlleleEncodings = alleleEncodings;
        return this;
    }

    public int getTaxaCount() {
        return myGenotype.getNumRows();
    }

    public int getSiteCount() {
        return myGenotype.getNumColumns();
    }

    public void reorderTaxa(int[] newIndices) {
        myGenotype.reorderRows(newIndices);
    }

    public void reorderPositions(int[] newIndices) {
        myGenotype.reorderColumns(newIndices);
    }

    public Genotype build() {
        SuperByteMatrix temp = myGenotype;
        myGenotype = null;
        if (NucleotideAlignmentConstants.isNucleotideEncodings(myAlleleEncodings)) {
            return new NucleotideGenotype(temp, myIsPhased);
        } else {
            return new ByteGenotype(temp, myIsPhased, myAlleleEncodings);
        }
    }

    public static Genotype buildHDF5(String filename) {
        return HDF5ByteGenotype.getInstance(HDF5Factory.openForReading(filename));
    }

    public static Genotype buildHDF5(IHDF5Reader reader) {
        return HDF5ByteGenotype.getInstance(reader);
    }
}
