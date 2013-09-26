package net.maizegenetics.pal.alignment.io;

import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Reader;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.AlignmentBuilder;
import net.maizegenetics.pal.alignment.genotype.Genotype;
import net.maizegenetics.pal.alignment.genotype.GenotypeBuilder;
import net.maizegenetics.pal.position.PositionList;
import net.maizegenetics.pal.position.PositionListBuilder;
import net.maizegenetics.pal.taxa.TaxaList;
import net.maizegenetics.pal.taxa.TaxaListBuilder;
import net.maizegenetics.pal.taxa.Taxon;
import org.apache.log4j.Logger;

import java.util.ArrayList;
import java.util.List;
import java.util.TreeSet;

/**
 * Defines xxxx
 *
 * @author Ed Buckler
 */
public class BuilderFromGenotypeHDF5 {
    private static final Logger myLogger=Logger.getLogger(BuilderFromGenotypeHDF5.class);
    private final String infile;

    private BuilderFromGenotypeHDF5(String infile) {
        this.infile=infile;
    }

    public static BuilderFromGenotypeHDF5 getBuilder(String infile) {
        return new BuilderFromGenotypeHDF5(infile);
    }

    //TODO provide options on caching to use, read only some sites, etc.
    //TODO update to the newest version
    //TODO subset??
    public Alignment build() {
        IHDF5Reader reader=HDF5Factory.openForReading(infile);
        TaxaList tL=new TaxaListBuilder().buildFromHDF5(reader);
        PositionList pL=new PositionListBuilder(reader).build();
        Genotype geno=GenotypeBuilder.buildHDF5(reader);
        return AlignmentBuilder.getInstance(geno,pL, tL);
    }

    /**
     * This merge multiple alignment together into one ByteNucleotideHDF5 File.
     * This is designed for putting multiple chromosomes together into one whole genome file.
     * @param infiles array of input alignment names
     * @param newMerge name of ByteNucleotideHDF5
     */
    public static void mergeToMutableHDF5(String[] infiles, String newMerge) {
        if ((infiles == null) || (infiles.length == 0)) {
            return ;
        }
        System.out.println("Opening Existing Position List");
        List<PositionList> inPL=new ArrayList<>();
        PositionListBuilder palBuild=new PositionListBuilder();
        System.out.println("Combining Position List");
        for (String infile : infiles) {
            PositionList pl=new PositionListBuilder(infile).build();
            inPL.add(pl);
            palBuild.addAll(pl);
        }
        System.out.println("Sorting Position List");
        PositionList pal=palBuild.build(); //In memory position list
//        for (int i=0; i<pal.getSiteCount(); i+=10000) {
//            System.out.println(inPL.get(0).get(i).toString());
//            System.out.println(pal.get(i).toString());
//        }
        System.out.println("Creating Position List Look Up");
        int[][] oldSiteToNewSite=new int[inPL.size()][];
        int misses=0;
        int hits=0;
        for (int i=0; i<inPL.size(); i++) {
            PositionList aPL=inPL.get(i);
            oldSiteToNewSite[i]=new int[aPL.getSiteCount()];
            for (int j=0; j<aPL.size(); j++) {
                oldSiteToNewSite[i][j]=pal.indexOf(aPL.get(j));
                if(oldSiteToNewSite[i][j]<0) {
                    misses++;
                    System.out.println(oldSiteToNewSite[i][j]);
                    System.out.println(aPL.get(j).toString());
                    System.out.println(pal.get(-oldSiteToNewSite[i][j]).toString());
                    System.out.printf("misses:%d hits:%d %n", misses, hits);
                } else {
                    hits++;
                }
            }
        }
        System.out.println(misses);
        inPL=null;
        int numberOfSites=pal.getSiteCount();
        //Get taxa List
        List<TaxaList> inTL=new ArrayList<>();
        TreeSet<Taxon> taxa = new TreeSet<>();
        for (String infile : infiles) {
            IHDF5Reader reader=HDF5Factory.openForReading(infile);
            TaxaList aTL=new TaxaListBuilder().buildFromHDF5(reader);
            reader.close();
            taxa.addAll(aTL);
            inTL.add(aTL);
        }
        TaxaList newTaxaList=new TaxaListBuilder().addAll(taxa).build();
        int[][] oldTaxaToNewTaxa=new int[inTL.size()][];
        for (int i=0; i<inTL.size(); i++) {
            TaxaList aTL=inTL.get(i);
            oldTaxaToNewTaxa[i]=new int[aTL.getTaxaCount()];
            for (int j=0; j<aTL.size(); j++) {
                oldTaxaToNewTaxa[i][j]=newTaxaList.indexOf(aTL.get(j));
                if(oldTaxaToNewTaxa[i][j]<0) {
                    System.out.println(oldTaxaToNewTaxa[i][j]);
                    System.out.println(aTL.get(j).toString());
                }
            }
        }
        AlignmentBuilder ab=AlignmentBuilder.getTaxaIncremental(pal,newMerge);
        System.gc();

//        System.out.println(Arrays.deepToString(oldTaxaToNewTaxa));

        //Transfer the genotypes
        System.out.println("Opening alignments");
        List<Alignment> sourceA=new ArrayList<>();
        for (String infile : infiles) {
            sourceA.add(AlignmentBuilder.getInstance(infile));
        }


        for (Taxon aT : newTaxaList) {
            System.out.println("Write taxon:"+aT.getFullName());
            byte[] geno=new byte[numberOfSites];
            for (int i=0; i<sourceA.size(); i++) {
                int taxonIndex=sourceA.get(i).getTaxaList().getIndicesMatchingTaxon(aT).get(0);
                byte[] r=sourceA.get(i).getBaseRow(taxonIndex);
                for (int j=0; j<oldSiteToNewSite[i].length; j++) {
                    geno[oldSiteToNewSite[i][j]]=r[j];
                }
            }
            ab.addTaxon(aT,geno);
        }

    }


}
