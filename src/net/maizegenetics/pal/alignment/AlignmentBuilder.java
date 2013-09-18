/*
 *  AlignmentBuilder
 */
package net.maizegenetics.pal.alignment;

import net.maizegenetics.pal.alignment.depth.AlleleDepth;
import net.maizegenetics.pal.alignment.genotype.Genotype;
import net.maizegenetics.pal.alignment.genotype.GenotypeBuilder;
import net.maizegenetics.pal.alignment.score.SiteScore;
import net.maizegenetics.pal.site.Position;
import net.maizegenetics.pal.site.PositionArrayList;
import net.maizegenetics.pal.site.PositionList;
import net.maizegenetics.pal.taxa.TaxaList;
import net.maizegenetics.pal.taxa.TaxaListBuilder;
import net.maizegenetics.pal.taxa.Taxon;

import java.util.ArrayList;

/**
 *
 * @author Terry Casstevens
 * @author Ed Buckler
 */
public class AlignmentBuilder {

    private Genotype genotype=null;

    //Fields for incremental taxa
    private PositionList positionList=null;
    private TaxaListBuilder taxaListBuilder=null;
    private ArrayList<byte[]> incGeno=null;

    //Fields for incremental sites
    private TaxaList taxaList=null;
    private PositionArrayList.Builder posListBuilder=null;

    private enum BuildType{TAXA_INC, SITE_INC, GENO_EDIT};
    private BuildType myBuildType;

    private AlignmentBuilder(PositionList positionList) {
        this.positionList=positionList;
        this.myBuildType=BuildType.TAXA_INC;
        incGeno=new ArrayList<>();
        taxaListBuilder=new TaxaListBuilder();
    }

    private AlignmentBuilder(TaxaList taxaList) {
        this.taxaList=taxaList;
        this.myBuildType=BuildType.SITE_INC;
        incGeno=new ArrayList<>();
        taxaListBuilder=new TaxaListBuilder();
    }

    public static AlignmentBuilder getTaxaIncremental(PositionList positionList) {
        return new AlignmentBuilder(positionList);
    }

    public static AlignmentBuilder getSiteIncremental(TaxaList taxaList) {
        return new AlignmentBuilder(taxaList);
    }

    public AlignmentBuilder addSite(Position pos, byte[] genos) {
        if(posListBuilder==null) throw new IllegalArgumentException("addSite only be used with AlignmentBuilder.getSiteIncremental");
        if(genos.length!=taxaList.getTaxaCount()) throw new IndexOutOfBoundsException("Number of taxa and genotypes do not agree");
        posListBuilder.add(pos);
        incGeno.add(genos);
        return this;
    }

    public AlignmentBuilder addTaxon(Taxon taxon, byte[] genos) {
        if(taxaListBuilder==null) throw new IllegalArgumentException("addSite only be used with AlignmentBuilder.getTaxaIncremental");
        if(genos.length!=positionList.getSiteCount()) throw new IndexOutOfBoundsException("Number of sites and genotypes do not agree");
        taxaListBuilder.add(taxon);
        incGeno.add(genos);
        return this;
    }

    public Alignment build(){
        switch (myBuildType) {
            case TAXA_INC: {
                //TODO optional sort
                TaxaList tl=taxaListBuilder.build();
                GenotypeBuilder gB=GenotypeBuilder.getInstance(tl.getTaxaCount(),positionList.getSiteCount());
                for (int i=0; i<incGeno.size(); i++) {
                    gB.setBaseRangeForTaxon(i, 0, incGeno.get(i));
                }
                return new CoreAlignment(gB.build(), positionList, tl);
            }
            case SITE_INC: {
                //TODO validate sort order, sort if needed
                PositionList pl=posListBuilder.build();
                GenotypeBuilder gB=GenotypeBuilder.getInstance(taxaList.getTaxaCount(),pl.getSiteCount());
                for (int i=0; i<taxaList.getTaxaCount(); i++) {
                    byte[] b=incGeno.get(i);
                    for (int s=0; s<b.length; s++) {
                        gB.setBase(i,s,b[s]);
                    }

                }
                return new CoreAlignment(gB.build(), pl, taxaList);
            }
        }
        return null;
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
        return copyGenotypeInstance(alignment);
    }

    /**
     * Returns a taxa optimized version of a combine alignment.  Only needed in performance critical situations
     * like imputation.
     * @param alignment
     * @return alignment backed by a single SuperByteMatrix
     */
    public static Alignment getGenotypeCopyInstance(CombineAlignment alignment) {
        return copyGenotypeInstance(alignment);
    }

    /**
     * This is private as the method used by combine and filter alignment, as other data structures
     * are immutable and optimized that is unneeded except for these alignment types.
     * @param alignment
     * @return alignment backed by a single SuperByteMatrix
     */
    private static Alignment copyGenotypeInstance(Alignment alignment) {
        int numTaxa = alignment.getTaxaCount();
        int numSites = alignment.getSiteCount();
        GenotypeBuilder builder = GenotypeBuilder.getInstance(numTaxa, numSites);
        for (int t = 0; t < numTaxa; t++) {
            for (int s = 0; s < numSites; s++) { builder.setBase(t, s, alignment.getBase(t, s));}
        }
        return new CoreAlignment(builder.build(), alignment.getPositionList(), alignment.getTaxaList());
    }


}
