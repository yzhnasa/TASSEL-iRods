/*
 *  AlignmentBuilder
 */
package net.maizegenetics.pal.alignment;

import ch.systemsx.cisd.base.mdarray.MDArray;
import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Reader;
import ch.systemsx.cisd.hdf5.IHDF5Writer;
import ch.systemsx.cisd.hdf5.IHDF5WriterConfigurator;
import net.maizegenetics.pal.alignment.depth.AlleleDepth;
import net.maizegenetics.pal.alignment.genotype.Genotype;
import net.maizegenetics.pal.alignment.genotype.GenotypeBuilder;
import net.maizegenetics.pal.alignment.score.SiteScore;
import net.maizegenetics.pal.position.Position;
import net.maizegenetics.pal.position.PositionArrayList;
import net.maizegenetics.pal.position.PositionHDF5List;
import net.maizegenetics.pal.position.PositionList;
import net.maizegenetics.pal.taxa.TaxaList;
import net.maizegenetics.pal.taxa.TaxaListBuilder;
import net.maizegenetics.pal.taxa.Taxon;
import net.maizegenetics.util.HDF5Utils;

import java.util.ArrayList;
import java.util.Arrays;

/**
 * Create new alignments.  New alignments are built from a minimum of TaxaList, PositionList, and Genotypes.  Depth and Scores are optional
 * features of alignments.
 * <p></p>
 * If you know the taxa,position, and genotypes are known from the beginning use:
 * Alignment a=AlignmentBuilder.getInstance(genotype, positionList, taxaList);
 *
 * In many situations only
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
    private boolean isHDF5=false;
    private IHDF5Writer writer=null;
    private BuildType myBuildType;

    /*
    Builder for in memory taxa incremental
     */
    private AlignmentBuilder(PositionList positionList) {
        this.positionList=positionList;
        this.myBuildType=BuildType.TAXA_INC;
        incGeno=new ArrayList<>();
        taxaListBuilder=new TaxaListBuilder();
    }

    /*
        Builder for in memory site incremental
     */
    private AlignmentBuilder(TaxaList taxaList) {
        this.taxaList=taxaList;
        this.myBuildType=BuildType.SITE_INC;
        incGeno=new ArrayList<>();
        taxaListBuilder=new TaxaListBuilder();
    }

    public static AlignmentBuilder getTaxaIncremental(PositionList positionList) {
        return new AlignmentBuilder(positionList);
    }

    public static AlignmentBuilder getTaxaIncremental(PositionList positionList, String newHDF5File) {
        return new AlignmentBuilder(positionList,newHDF5File);
    }

    public static AlignmentBuilder getSiteIncremental(TaxaList taxaList) {
        return new AlignmentBuilder(taxaList);
    }

    public AlignmentBuilder addSite(Position pos, byte[] genos) {
        if((myBuildType!=BuildType.SITE_INC)||isHDF5) throw new IllegalArgumentException("addSite only be used with AlignmentBuilder.getSiteIncremental and without HDF5");
        if(genos.length!=taxaList.getTaxaCount()) throw new IndexOutOfBoundsException("Number of taxa and genotypes do not agree");
        posListBuilder.add(pos);
        incGeno.add(genos);
        return this;
    }

    public AlignmentBuilder addTaxon(Taxon taxon, byte[] genos) {
        return addTaxon(taxon, genos, null);
    }

    public AlignmentBuilder addTaxon(Taxon taxon, byte[] genos, byte[] depth) {
        if(myBuildType!=BuildType.TAXA_INC) throw new IllegalArgumentException("addTaxon only be used with AlignmentBuilder.getTaxaIncremental");
        if(genos.length!=positionList.getSiteCount()) throw new IndexOutOfBoundsException("Number of sites and genotypes do not agree");
        if(isHDF5) {
            addTaxon(writer, taxon, genos, null);

        } else {
            taxaListBuilder.add(taxon);
            incGeno.add(genos);
        }
        return this;
    }

    public boolean isHDF5() {
        return isHDF5;
    }


    public Alignment build(){
        if(isHDF5) {
            String name=writer.getFile().getAbsolutePath();
            annotateHDF5File(writer);
            writer.close();
            return getInstance(name);
        }
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

    /**
     * Standard approach for creating a new Alignment
     * @param genotype
     * @param positionList
     * @param taxaList
     * @return new alignment
     */
    public static Alignment getInstance(Genotype genotype, PositionList positionList, TaxaList taxaList) {
        return new CoreAlignment(genotype, positionList, taxaList);
    }

    /**
     * Creates a new HDF5 file alignment based on existing Genotype, PositionList, and TaxaList.
     * @param genotype
     * @param positionList
     * @param taxaList
     * @param hdf5File name of the file
     * @return alignment backed by new HDF5 file
     */
    public static Alignment getInstance(Genotype genotype, PositionList positionList, TaxaList taxaList, String hdf5File) {
        AlignmentBuilder aB=AlignmentBuilder.getTaxaIncremental(positionList,hdf5File);
        for (int i=0; i<taxaList.getTaxaCount(); i++) {
            aB.addTaxon(taxaList.get(i),genotype.getBaseRow(i));
        }
        return aB.build();
    }

    /**
     * Creates a new HDF5 file alignment based on an existing alignment.
     * @param a existing alignment
     * @param hdf5File name of the file
     * @return alignment backed by new HDF5 file
     */
    public static Alignment getInstance(Alignment a, String hdf5File) {
        return getInstance(a.getGenotypeMatrix(),a.getPositionList(),a.getTaxaList(),hdf5File);
    }

    public static Alignment getInstance(String hdf5File) {
        IHDF5Reader reader=HDF5Factory.openForReading(hdf5File);
        TaxaList tL=new TaxaListBuilder().buildFromHDF5(reader);
        PositionList pL=new PositionHDF5List.Builder(reader).build();
        Genotype geno=GenotypeBuilder.buildHDF5(reader);
        return AlignmentBuilder.getInstance(geno, pL, tL);
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

    /*
    HDF5 Alignment section.
     */

    /**
     * Creates a new HDF5 file that can used with TaxaIncremental addition.
     * @param positionList
     * @param hdf5File
     */
    private AlignmentBuilder(PositionList positionList, String hdf5File) {
        IHDF5WriterConfigurator config = HDF5Factory.configure(hdf5File);
        //config.overwrite();
        config.dontUseExtendableDataTypes();
        writer = config.writer();
        this.positionList=new PositionHDF5List.Builder(writer,positionList).build();  //create a new position list

        setupGenotypeTaxaInHDF5(writer);
        this.myBuildType=BuildType.TAXA_INC;
        isHDF5=true;
        taxaListBuilder=new TaxaListBuilder();
    }



    private synchronized void setupGenotypeTaxaInHDF5(IHDF5Writer writer) {
        writer.setIntAttribute(HapMapHDF5Constants.DEFAULT_ATTRIBUTES_PATH, HapMapHDF5Constants.MAX_NUM_ALLELES,
                NucleotideAlignmentConstants.NUMBER_NUCLEOTIDE_ALLELES);
        writer.setIntAttribute(HapMapHDF5Constants.DEFAULT_ATTRIBUTES_PATH, HapMapHDF5Constants.BLOCK_SIZE,
                1<<16);
        writer.setBooleanAttribute(HapMapHDF5Constants.DEFAULT_ATTRIBUTES_PATH, HapMapHDF5Constants.RETAIN_RARE_ALLELES, false);
        writer.setIntAttribute(HapMapHDF5Constants.DEFAULT_ATTRIBUTES_PATH, HapMapHDF5Constants.NUM_TAXA, 0);
        String[][] aEncodings = NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES;
        int numEncodings = aEncodings.length;
        int numStates = aEncodings[0].length;
        MDArray<String> alleleEncodings = new MDArray<>(String.class, new int[]{numEncodings, numStates});
        for (int s = 0; s < numEncodings; s++) {
            for (int x = 0; x < numStates; x++) {
                alleleEncodings.set(aEncodings[s][x], s, x);
            }
        }
        writer.createStringMDArray(HapMapHDF5Constants.ALLELE_STATES, 100, new int[]{numEncodings, numStates});
        writer.writeStringMDArray(HapMapHDF5Constants.ALLELE_STATES, alleleEncodings);
        MDArray<String> alleleEncodingReadAgain = writer.readStringMDArray(HapMapHDF5Constants.ALLELE_STATES);
        if (alleleEncodings.equals(alleleEncodingReadAgain) == false) {
            throw new IllegalStateException("ExportUtils: writeToMutableHDF5: Mismatch Allele States, expected '" + alleleEncodings + "', found '" + alleleEncodingReadAgain + "'!");
        }
        writer.createGroup(HapMapHDF5Constants.GENOTYPES);
    }

    /**
     * Code needed to add a Taxon to HDF5, potentially split into functions in TaxaListBuilder & GenotypeBuilder
     */
    private synchronized void addTaxon(IHDF5Writer myWriter, Taxon id, byte[] genotype, byte[][] depth) {
        int chunk=1<<16;
        int myNumSites=positionList.getSiteCount();
        if(myNumSites<chunk) chunk=myNumSites;
        String basesPath = HapMapHDF5Constants.GENOTYPES + "/" + id.getFullName();
        if(myWriter.exists(basesPath)) throw new IllegalStateException("Taxa Name Already Exists:"+basesPath);
        if(genotype.length!=myNumSites) throw new IllegalStateException("Setting all genotypes in addTaxon.  Wrong number of sites");
        myWriter.createByteArray(basesPath, myNumSites, chunk, HapMapHDF5Constants.intDeflation);
        HDF5Utils.writeHDF5EntireArray(basesPath, myWriter, genotype.length, 1<<16, genotype);
        if(depth!=null) {
            if(depth.length!=6) throw new IllegalStateException("Just set A, C, G, T, -, + all at once");
            if(depth[0].length!=myNumSites) throw new IllegalStateException("Setting all depth in addTaxon.  Wrong number of sites");
            //           myWriter.writeByteMatrix(getTaxaDepthPath(taxonIndex), depth, HapMapHDF5Constants.intDeflation);
        }
    }

    /**
     * Annotates the HDF5 HapMap file after it is built.
     * Currently, placed in the AlignmentBuilder as it still above genotypes, taxa, and sites.
     * @param writer
     */
    private void annotateHDF5File(IHDF5Writer writer) {
        int hdf5GenoBlock=writer.getIntAttribute(HapMapHDF5Constants.DEFAULT_ATTRIBUTES_PATH, HapMapHDF5Constants.BLOCK_SIZE);
        int sites=writer.getIntAttribute(HapMapHDF5Constants.DEFAULT_ATTRIBUTES_PATH, HapMapHDF5Constants.NUM_SITES);
        TaxaList tL=new TaxaListBuilder().buildFromHDF5(writer);
        int taxa=tL.getTaxaCount();
        writer.setIntAttribute(HapMapHDF5Constants.DEFAULT_ATTRIBUTES_PATH, HapMapHDF5Constants.NUM_TAXA,taxa);
        int[][] af=new int[NucleotideAlignmentConstants.NUMBER_NUCLEOTIDE_ALLELES][sites];
        byte[][] afOrder=new byte[NucleotideAlignmentConstants.NUMBER_NUCLEOTIDE_ALLELES][sites];
        float[] coverage=new float[taxa];
        float[] hets=new float[taxa];
        for (int taxon = 0; taxon < taxa; taxon++) {
            String basesPath = HapMapHDF5Constants.GENOTYPES + "/" + tL.getFullTaxaName(taxon);
            byte[] genotype=writer.readByteArray(basesPath);
            int covSum=0;  //coverage of the taxon
            int hetSum=0;
            for (int s = 0; s < sites; s++) {
                byte[] b = AlignmentUtils.getDiploidValues(genotype[s]);
                if(b[0]<6) af[b[0]][s]++;
                if(b[1]<6) af[b[1]][s]++;
                if(AlignmentUtils.isHeterozygous(genotype[s])) hetSum++;
                if(genotype[s]!=Alignment.UNKNOWN_DIPLOID_ALLELE) covSum++;
            }
            coverage[taxon]=(float)covSum/(float)sites;
            hets[taxon]=(float)hetSum/(float)covSum;
        }
        float[] maf=new float[sites];
        float[] paf=new float[sites];
        int baseMask=0xF;
        for (int s = 0; s < sites; s++) {
            int sum=0;
            int[] cntAndAllele=new int[NucleotideAlignmentConstants.NUMBER_NUCLEOTIDE_ALLELES];
            for (byte i = 0; i < 6; i++) {
                cntAndAllele[i]=(af[i][s]<<4)|(5-i);  //size | allele (the 5-i is to get the sort right, so if case of ties A is first)
                sum+=af[i][s];
            }
            Arrays.sort(cntAndAllele);  //ascending quick sort, there are faster ways
            //http://stackoverflow.com/questions/2786899/fastest-sort-of-fixed-length-6-int-array
            for (byte i = 0; i < 6; i++) {
                afOrder[5-i][s]=(cntAndAllele[i]>0xF)?((byte)(5-(baseMask&cntAndAllele[i]))):Alignment.UNKNOWN_ALLELE;
            }
            if(afOrder[1][s]!=Alignment.UNKNOWN_ALLELE) maf[s]=(float)af[afOrder[1][s]][s]/(float)sum;
            paf[s]=(float)sum/(float)(2*taxa);
        }
        writer.createGroup(HapMapHDF5Constants.SITE_DESC);
        int chunk=(sites<hdf5GenoBlock)?sites:hdf5GenoBlock;
        writer.createIntMatrix(HapMapHDF5Constants.ALLELE_CNT, 6, sites, 1, chunk, HapMapHDF5Constants.intDeflation);
        writer.createByteMatrix(HapMapHDF5Constants.ALLELE_FREQ_ORD, 6, sites, 1, chunk, HapMapHDF5Constants.intDeflation);
        writer.createFloatArray(HapMapHDF5Constants.MAF, sites, chunk, HapMapHDF5Constants.floatDeflation);
        writer.createFloatArray(HapMapHDF5Constants.SITECOV, sites, chunk, HapMapHDF5Constants.floatDeflation);
        writer.createGroup(HapMapHDF5Constants.TAXA_DESC);
        chunk=(tL.getTaxaCount()<hdf5GenoBlock)?tL.getTaxaCount():hdf5GenoBlock;
        writer.createFloatArray(HapMapHDF5Constants.TAXACOV, tL.getTaxaCount(), chunk, HapMapHDF5Constants.floatDeflation);
        writer.createFloatArray(HapMapHDF5Constants.TAXAHET, tL.getTaxaCount(), chunk, HapMapHDF5Constants.floatDeflation);
        if(af[0].length>0) HDF5Utils.writeHDF5EntireArray(HapMapHDF5Constants.ALLELE_CNT, writer, af[0].length, 1<<16, af);
        if(afOrder[0].length>0) HDF5Utils.writeHDF5EntireArray(HapMapHDF5Constants.ALLELE_FREQ_ORD, writer, afOrder[0].length, 1<<16, afOrder);
        if(maf.length>0) HDF5Utils.writeHDF5EntireArray(HapMapHDF5Constants.MAF, writer, maf.length, 1<<16, maf);
        if(paf.length>0) HDF5Utils.writeHDF5EntireArray(HapMapHDF5Constants.SITECOV, writer, paf.length, 1<<16, paf);

        if(coverage.length>0) HDF5Utils.writeHDF5EntireArray(HapMapHDF5Constants.TAXACOV, writer, coverage.length, 1<<16, coverage);
        if(hets.length>0) {
            System.out.println("Hets Length:"+hets.length);
            HDF5Utils.writeHDF5EntireArray(HapMapHDF5Constants.TAXAHET, writer, hets.length, 1<<16, hets);
        }
    }
}
