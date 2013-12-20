/*
 *  HDF5AlleleDepth
 */
package net.maizegenetics.dna.snp.depth;

import ch.systemsx.cisd.hdf5.IHDF5Reader;

import java.util.LinkedHashMap;
import java.util.Map;
import net.maizegenetics.dna.snp.HapMapHDF5Constants;
import net.maizegenetics.taxa.TaxaList;

/**
 *
 * @author Terry Casstevens
 */
public class HDF5AlleleDepth extends AbstractAlleleDepth {
    
    private static int MAX_CACHE_SIZE = 1 << 16;
    private static final int HDF5_BLOCK = 1 << 16;
    private final Map<Long, byte[][]> myDepthCache = new LinkedHashMap<Long, byte[][]>((3 * MAX_CACHE_SIZE) / 2) {
        @Override
        protected boolean removeEldestEntry(Map.Entry eldest) {
            return size() > MAX_CACHE_SIZE;
        }
    };
    
    private final IHDF5Reader myReader;
    private final int myNumSites;
    private final TaxaList myTaxa;
    
    HDF5AlleleDepth(IHDF5Reader reader) {
        super(6);
        myReader = reader;
        myNumSites = reader.getIntAttribute(HapMapHDF5Constants.DEFAULT_ATTRIBUTES_PATH, HapMapHDF5Constants.NUM_SITES);
        // TODO: Some how we need to get taxa names.  Either passed
        // in as a parameter or read from the HDF5 file.
        myTaxa = null;
    }
    
    private static long getCacheKey(int taxon, int site) {
        return ((long) taxon << 33) + (site / HDF5_BLOCK);
    }
    
    public byte[] depthForAllelesBytes(int taxon, int site) {
        long key = getCacheKey(taxon, site);
        byte[][] data = myDepthCache.get(key);
        if (data == null) {
            data = cacheDepthBlock(taxon, site, key);
        }
        byte[] result = new byte[6];
        for (int i = 0; i < 6; i++) {
            result[i] = data[i][site % MAX_CACHE_SIZE];
        }
        return result;
    }
    
    public byte[][] depthForAllelesBytes(int taxon) {
        byte[][] result = new byte[6][myNumSites];
        for (int site = 0; site < myNumSites; site++) {
            long key = getCacheKey(taxon, site);
            byte[][] data = myDepthCache.get(key);
            if (data == null) {
                data = cacheDepthBlock(taxon, site, key);
            }
            for (int i = 0; i < 6; i++) {
                result[i][site] = data[i][site % MAX_CACHE_SIZE];
            }
        }
        return result;
    }
    
    private String getTaxaDepthPath(int taxon) {
        return HapMapHDF5Constants.DEPTH + "/" + myTaxa.taxaName(taxon);
    }
    
    private byte[][] cacheDepthBlock(int taxon, int site, long key) {
        int start = (site / MAX_CACHE_SIZE) * MAX_CACHE_SIZE;
        int realSiteCache = (myNumSites - start < MAX_CACHE_SIZE) ? myNumSites - start : MAX_CACHE_SIZE;
        byte[][] data = myReader.readByteMatrixBlockWithOffset(getTaxaDepthPath(taxon), 6, realSiteCache, 0, start);
        if (data == null) {
            return null;
        }
        myDepthCache.put(key, data);
        return data;
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
