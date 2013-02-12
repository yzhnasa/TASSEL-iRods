/*
 * HapMapHDF5Constants
 */
package net.maizegenetics.pal.alignment;

/**
 *
 * @author terry
 */
public final class HapMapHDF5Constants {
    // Paths
    public static final String ROOT = "/";
    public static final String TAXA = "Taxa";
    public static final String ALLELE_STATES = "AlleleStates";
    public static final String POSITIONS = "Positions";
    public static final String ALLELES = "Alleles";
    public static final String TBIT = "TBit";
    public static final String SBIT = "SBit";
    public static final String LOCI = "SeqRegion";
    public static final String LOCUS_OFFSETS = "SeqRegionOffsets";
    public static final String SNP_IDS = "SnpIds";
    
    // Attributes
    public static final String DEFAULT_ATTRIBUTES_PATH = ROOT;
    public static final String NUM_SBIT_WORDS = "numSBitWords";
    public static final String NUM_TBIT_WORDS = "numTBitWords";
    public static final String NUM_TAXA = "numTaxa";
    public static final String NUM_SITES = "numSites";
    public static final String MAX_NUM_ALLELES = "maxNumAlleles";
    public static final String RETAIN_RARE_ALLELES = "retainRareAlleles";
    
    private HapMapHDF5Constants() {
        // do not instantiate
    }
    
}
