/*
 *  AlignmentBuilder
 */
package net.maizegenetics.pal.alignment;

import net.maizegenetics.pal.alignment.depth.AlleleDepth;
import net.maizegenetics.pal.alignment.genotype.Genotype;
import net.maizegenetics.pal.alignment.genotype.GenotypeBuilder;
import net.maizegenetics.pal.alignment.score.SiteScore;
import net.maizegenetics.pal.site.PositionList;
import net.maizegenetics.pal.taxa.TaxaList;

/**
 *
 * @author Terry Casstevens
 */
public class AlignmentBuilder {

    private AlignmentBuilder() {
    }

    public static Alignment getInstance(Genotype genotype, PositionList positionList, TaxaList taxaList, SiteScore siteScore, AlleleDepth alleleDepth) {
        return new CoreAlignment(genotype, positionList, taxaList, siteScore, alleleDepth);
    }

    public static Alignment getInstance(Genotype genotype, PositionList positionList, TaxaList taxaList) {
        return new CoreAlignment(genotype, positionList, taxaList);
    }

    public static Alignment getInstanceOnlyMajorMinor(Alignment alignment) {
        int numTaxa = alignment.getTaxaCount();
        int numSites = alignment.getSiteCount();
        GenotypeBuilder builder = GenotypeBuilder.getInstance(numTaxa, numSites);
        byte[] majorAllele = new byte[64];
        byte[] minorAllele = new byte[64];
        for (int bigS = 0; bigS < numSites; bigS += 64) {
            int blockSize = Math.min(64, numSites - bigS);

            for (int s = 0; s < blockSize; s++) {
                majorAllele[s] = alignment.getMajorAllele(s + bigS);
                minorAllele[s] = alignment.getMinorAllele(s + bigS);
            }

            for (int t = 0; t < numTaxa; t++) {
                for (int s = 0; s < blockSize; s++) {
                    byte[] currentAlleles = alignment.getBaseArray(t, s + bigS);
                    if ((currentAlleles[0] != majorAllele[s]) && (currentAlleles[0] != minorAllele[s])) {
                        currentAlleles[0] = Alignment.UNKNOWN_ALLELE;
                    }
                    if ((currentAlleles[1] != majorAllele[s]) && (currentAlleles[1] != minorAllele[s])) {
                        currentAlleles[1] = Alignment.UNKNOWN_ALLELE;
                    }
                    builder.setBase(t, s, AlignmentUtils.getDiploidValue(currentAlleles[0], currentAlleles[1]));
                }
            }
        }
        return new CoreAlignment(builder.build(), alignment.getPositionList(), alignment.getTaxaList());
    }
}
