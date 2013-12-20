/*
 *  HDF5AlleleDepth
 */
package net.maizegenetics.dna.snp.depth;

import ch.systemsx.cisd.hdf5.IHDF5Reader;

/**
 *
 * @author Terry Casstevens
 */
public class HDF5AlleleDepth extends AbstractAlleleDepth {

    private final IHDF5Reader myReader;

    HDF5AlleleDepth(IHDF5Reader reader) {
        super(6);
        myReader = reader;
    }

    @Override
    public int[] depthForAlleles(int taxon, int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int depthForAllele(int taxon, int site, int allele) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public byte depthForAlleleByte(int taxon, int site, int allele) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

}
