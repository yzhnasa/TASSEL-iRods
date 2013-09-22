/*
 *  HDF5ByteGenotype
 */
package net.maizegenetics.pal.alignment.genotype;

import ch.systemsx.cisd.hdf5.IHDF5Reader;
import com.google.common.cache.CacheBuilder;
import com.google.common.cache.CacheLoader;
import com.google.common.cache.LoadingCache;
import net.maizegenetics.pal.alignment.HapMapHDF5Constants;
import net.maizegenetics.pal.alignment.NucleotideAlignmentConstants;
import net.maizegenetics.pal.taxa.TaxaList;
import net.maizegenetics.pal.taxa.TaxaListBuilder;

/**
 *
 * @author Ed Buckler
 * @author Terry Casstevens
 */
public class HDF5ByteGenotype extends AbstractGenotype {

    private static final int SHIFT_AMOUNT = 16;
    private final String[] genotypePaths;
    /**
     * Byte representations of DNA sequences are stored in blocks of 65536 sites
     */
    private static final int HDF5_GENOTYPE_BLOCK_SIZE = 1 << SHIFT_AMOUNT;
    private final IHDF5Reader myHDF5Reader;
    private final CacheLoader<Long, byte[]> myGenoLoader = new CacheLoader<Long, byte[]>() {
        public byte[] load(Long key) {
            long offset = getSiteStartFromKey(key) << SHIFT_AMOUNT;
            byte[] data;
            synchronized (myHDF5Reader) {
                data = myHDF5Reader.readAsByteArrayBlockWithOffset(getTaxaGenoPath(getTaxonFromKey(key)), HDF5_GENOTYPE_BLOCK_SIZE, offset);
            }
            return data;
        }
    };
    private final LoadingCache<Long, byte[]> myGenoCache;



    private static long getCacheKey(int taxon, int site) {
        return ((long) taxon << 33) + (site / HDF5_GENOTYPE_BLOCK_SIZE);
    }

    private static int getTaxonFromKey(long key) {
        return (int) (key >>> 33);
    }

    private static int getSiteStartFromKey(long key) {
        return (int) ((key << 33) >>> 33);
    }

    private String getTaxaGenoPath(int taxon) {
        return genotypePaths[taxon];
     //   return HapMapHDF5Constants.GENOTYPES + "/taxon" + taxon;
    }

    private HDF5ByteGenotype(IHDF5Reader reader, int numTaxa, int numSites, boolean phased, String[][] alleleEncodings) {
        super(numTaxa, numSites, phased, alleleEncodings);
        genotypePaths=new String[numTaxa];
        TaxaList tL=new TaxaListBuilder().buildFromHDF5(reader);  //not the most efficient thing to do, but ensures sort is the same.
        for (int i=0; i<numTaxa; i++) {
            genotypePaths[i]=HapMapHDF5Constants.GENOTYPES + "/" + tL.getFullTaxaName(i);
        }
        myHDF5Reader = reader;
        myGenoCache = CacheBuilder.newBuilder()
                .maximumSize((3 * getTaxaCount()) / 2)
                .build(myGenoLoader);
    }

    static HDF5ByteGenotype getInstance(IHDF5Reader reader) {
        int numTaxa=reader.getIntAttribute(HapMapHDF5Constants.DEFAULT_ATTRIBUTES_PATH, HapMapHDF5Constants.NUM_TAXA);
        int numSites=reader.getIntAttribute(HapMapHDF5Constants.DEFAULT_ATTRIBUTES_PATH, HapMapHDF5Constants.NUM_SITES);
        String[][] alleleEncodings=NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES;
        return new HDF5ByteGenotype(reader, numTaxa, numSites, false, alleleEncodings);
    }

    @Override
    public byte getBase(int taxon, int site) {
        long key = getCacheKey(taxon, site);
        try {
            byte[] data = myGenoCache.get(key);
            return data[site % HDF5_GENOTYPE_BLOCK_SIZE];
        } catch (Exception e) {
            e.printStackTrace();
            throw new IllegalStateException("HDF5ByteGenotype: getBase: Error getting base from cache.");
        }
    }

    @Override
    public void transposeData(boolean siteInnerLoop) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
}
