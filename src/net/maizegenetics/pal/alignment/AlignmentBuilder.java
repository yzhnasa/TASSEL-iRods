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

    public static Alignment getHomozygousInstance(Alignment alignment) {
        int numTaxa = alignment.getTaxaCount();
        int numSites = alignment.getSiteCount();
        GenotypeBuilder builder = GenotypeBuilder.getInstance(numTaxa, numSites);
        //TODO this would be even faster to work through the SuperByteMatrix, as knowledge of site or taxa is not needed.
        for (int t = 0; t < numTaxa; t++) {
            for (int s = 0; s < numSites; s++) {
                byte currGeno=alignment.getBase(t, s);
                if(AlignmentUtils.isHeterozygous(currGeno)) {
                    builder.setBase(t, s, Alignment.UNKNOWN_DIPLOID_ALLELE);
                } else {
                    builder.setBase(t, s, currGeno);
                }
            }

        }
        return new CoreAlignment(builder.build(), alignment.getPositionList(), alignment.getTaxaList());
    }

    /**
     * Returns a taxa optimized version of a filtered alignment.  Only needed in performance critical situations
     * like imputation.
     * @param alignment
     * @return alignment backed by a single SuperByteMatrix
     */
    public static Alignment getGenotypeCopyInstance(FilterAlignment alignment) {
        return getGenotypeCopyInstance(alignment);
    }

    /**
     * Returns a taxa optimized version of a combine alignment.  Only needed in performance critical situations
     * like imputation.
     * @param alignment
     * @return alignment backed by a single SuperByteMatrix
     */
    public static Alignment getGenotypeCopyInstance(CombineAlignment alignment) {
        return getGenotypeCopyInstance(alignment);
    }

    /**
     * This is private as only want the method used by combine and filter alignment.  Since other data structures
     * are immutable and optimized that is unneeded except for these alignment types.
     * @param alignment
     * @return alignment backed by a single SuperByteMatrix
     */
    private static Alignment getGenotypeCopyInstance(Alignment alignment) {
        int numTaxa = alignment.getTaxaCount();
        int numSites = alignment.getSiteCount();
        GenotypeBuilder builder = GenotypeBuilder.getInstance(numTaxa, numSites);
        for (int t = 0; t < numTaxa; t++) {
            for (int s = 0; s < numSites; s++) { builder.setBase(t, s, alignment.getBase(t, s));}
        }
        return new CoreAlignment(builder.build(), alignment.getPositionList(), alignment.getTaxaList());
    }


}
