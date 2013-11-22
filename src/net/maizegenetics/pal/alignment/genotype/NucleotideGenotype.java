/*
 *  NucleotideGenotype
 */
package net.maizegenetics.pal.alignment.genotype;

import net.maizegenetics.dna.snp.NucleotideAlignmentConstants;
import net.maizegenetics.util.SuperByteMatrix;

/**
 *
 * @author Terry Casstevens
 */
class NucleotideGenotype extends ByteGenotype {

    NucleotideGenotype(SuperByteMatrix genotype, boolean phased) {
        super(genotype, phased, NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES);
    }

    @Override
    public String getBaseAsString(int taxon, int site) {
        return NucleotideAlignmentConstants.getNucleotideIUPAC(getBase(taxon, site));
    }

    @Override
    public String getDiploidAsString(int site, byte value) {
        return NucleotideAlignmentConstants.getNucleotideIUPAC(value);
    }

    @Override
    public int getMaxNumAlleles() {
        return NucleotideAlignmentConstants.NUMBER_NUCLEOTIDE_ALLELES;
    }

    @Override
    public boolean retainsRareAlleles() {
        return false;
    }
}
