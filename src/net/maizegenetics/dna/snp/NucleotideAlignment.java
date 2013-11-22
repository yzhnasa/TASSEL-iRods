/*
 *  NucleotideAlignment
 */
package net.maizegenetics.dna.snp;

import net.maizegenetics.dna.snp.depth.AlleleDepth;
import net.maizegenetics.dna.snp.genotype.Genotype;
import net.maizegenetics.dna.snp.score.SiteScore;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.pal.taxa.TaxaList;

/**
 *
 * @author Terry Casstevens
 */
public class NucleotideAlignment extends CoreAlignment {

    public NucleotideAlignment(Genotype genotype, PositionList positionList, TaxaList taxaList, SiteScore siteScore, AlleleDepth alleleDepth) {
        super(genotype, positionList, taxaList, siteScore, alleleDepth);
    }

    @Override
    public String getBaseAsStringRange(int taxon, int startSite, int endSite) {
        StringBuilder builder = new StringBuilder();
        for (int i = startSite; i < endSite; i++) {
            builder.append(getBaseAsString(taxon, i));
        }
        return builder.toString();
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
    public boolean isIndel(int site) {
        int[][] alleles = getAllelesSortedByFrequency(site);
        int numAlleles = Math.min(alleles[0].length, 2);
        for (int i = 0; i < numAlleles; i++) {
            if ((alleles[0][i] == NucleotideAlignmentConstants.INSERT_ALLELE) || (alleles[0][i] == NucleotideAlignmentConstants.GAP_ALLELE)) {
                return true;
            }
        }
        return false;
    }
}
