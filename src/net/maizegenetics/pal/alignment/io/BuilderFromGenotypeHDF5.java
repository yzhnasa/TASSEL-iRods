package net.maizegenetics.pal.alignment.io;

import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Reader;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.AlignmentBuilder;
import net.maizegenetics.pal.alignment.genotype.Genotype;
import net.maizegenetics.pal.alignment.genotype.GenotypeBuilder;
import net.maizegenetics.pal.site.PositionHDF5List;
import net.maizegenetics.pal.site.PositionList;
import net.maizegenetics.pal.taxa.TaxaList;
import net.maizegenetics.pal.taxa.TaxaListBuilder;
import org.apache.log4j.Logger;

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
        TaxaList tL=new TaxaListBuilder().buildFromHDF5(infile);
        IHDF5Reader reader=HDF5Factory.openForReading(infile);
        PositionList pL=new PositionHDF5List.Builder(reader).build();
        Genotype geno=GenotypeBuilder.buildHDF5(reader);
        return AlignmentBuilder.getInstance(geno,pL, tL);
    }


}
