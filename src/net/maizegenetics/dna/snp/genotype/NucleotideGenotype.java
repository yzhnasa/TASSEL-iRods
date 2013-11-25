/*
 *  NucleotideGenotype
 */
package net.maizegenetics.dna.snp.genotype;

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
    public String genotypeAsString(int taxon, int site) {
        return NucleotideAlignmentConstants.getNucleotideIUPAC(genotype(taxon, site));
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
