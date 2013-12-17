/*
 *  HDF5AlleleDepth
 */
package net.maizegenetics.dna.snp.depth;

import ch.systemsx.cisd.hdf5.IHDF5Reader;

/**
 *
 * @author Terry Casstevens
 */
public class HDF5AlleleDepth implements AlleleDepth {

    private final IHDF5Reader myReader;

    public HDF5AlleleDepth(IHDF5Reader reader) {
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
    public int depth(int taxon, int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

}
