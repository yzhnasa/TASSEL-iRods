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

    private HDF5ByteGenotype(IHDF5Reader reader, int taxaCount, int siteCount, boolean phased, String[][] alleleEncodings) {
        super(taxaCount, siteCount, phased, alleleEncodings);
        myHDF5Reader = reader;
    }

    HDF5ByteGenotype getInstance() {
        return new HDF5ByteGenotype(null, 0, 0, false, null);
    }

    @Override
    public byte getBase(int taxon, int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
}
