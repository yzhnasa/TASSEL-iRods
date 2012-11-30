/*
 * HDF5Constants
 */
package net.maizegenetics.pal.alignment;

/**
 *
 * @author terry
 */
public final class HDF5Constants {
    // Paths
    public static final String ROOT = "/";
    public static final String TAXA = "Taxa";
    public static final String ALLELE_STATES = "AlleleStates";
    public static final String POSITIONS = "Positions";
    public static final String ALLELES = "Alleles";
    public static final String TBIT = "TBit";
    public static final String SBIT = "SBit";
    public static final String LOCI = "Loci";
    public static final String LOCUS_OFFSETS = "LocusOffsets";
    
    // Attributes
    public static final String NUM_WORDS_PATH = ROOT;
    public static final String NUM_SBIT_WORDS = "numSBitWords";
    public static final String NUM_TBIT_WORDS = "numTBitWords";
    public static final String NUM_TAXA_PATH = ROOT;
    public static final String NUM_TAXA = "numTaxa";
    public static final String NUM_SITES_PATH = ROOT;
    public static final String NUM_SITES = "numSites";
    
    private HDF5Constants() {
        // do not instantiate
    }
    
}
