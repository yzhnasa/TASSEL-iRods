/*
 *  HDF5ByteGenotype
 */
package net.maizegenetics.pal.alignment.genotype;

import ch.systemsx.cisd.hdf5.IHDF5Reader;

/**
 *
 * @author Terry Casstevens
 */
public class HDF5ByteGenotype extends AbstractGenotype {

    private final IHDF5Reader myHDF5Reader;

    private HDF5ByteGenotype(IHDF5Reader reader, int numTaxa, int numSites, boolean phased, String[][] alleleEncodings) {
        super(numTaxa, numSites, phased, alleleEncodings);
        myHDF5Reader = reader;
    }

    HDF5ByteGenotype getInstance() {
        return new HDF5ByteGenotype(null, 0, 0, false, null);
    }

    @Override
    public byte getBase(int taxon, int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public byte[] getGenotypeForAllSites(int taxon) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public byte[] getGenotypeForSiteRange(int taxon, int start, int end) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public byte[] getGenotypeForAllTaxa(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
}
