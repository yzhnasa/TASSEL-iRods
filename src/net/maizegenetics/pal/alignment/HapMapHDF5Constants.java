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
    public static final String SITE_DESC = "SiteDesc";
    public static final String TAXA_DESC = "TaxaDesc";
    public static final String MAF_DESC = SITE_DESC+"/MAF";
    public static final String SITECOV_DESC = SITE_DESC+"/SiteCoverage";
    public static final String HET_DESC = SITE_DESC+"/HET";
    public static final String LD_DESC = SITE_DESC+"/LD";
    public static final String LDR2_DESC = LD_DESC+"/R2";
    public static final String LDP_DESC = LD_DESC+"/P";
    public static final String ERROR_DESC = SITE_DESC+"/ERROR";
    public static final String BPECERROR_DESC = ERROR_DESC+"/BPECERROR";
    public static final String BPECAVGR2_DESC = ERROR_DESC+"/BPECR2";
    public static final String IBSERROR_DESC = ERROR_DESC+"/IBSError";
    public static final String IBSERRORTESTS_DESC = ERROR_DESC+"/IBSErrorTest";
    public static final String IBSMINORERROR_DESC = ERROR_DESC+"/IBSMinorError";
    public static final String IBSMINORERRORTESTS_DESC = ERROR_DESC+"/IBSMinorErrorTests";
    
    // Attributes
    public static final String DEFAULT_ATTRIBUTES_PATH = ROOT;
    public static final String NUM_SBIT_WORDS = "numSBitWords";
    public static final String NUM_TBIT_WORDS = "numTBitWords";
    public static final String NUM_TAXA = "numTaxa";
    public static final String NUM_SITES = "numSites";
    public static final String MAX_NUM_ALLELES = "maxNumAlleles";
    public static final String RETAIN_RARE_ALLELES = "retainRareAlleles";
    public static final String NUM_LD_BINS = LD_DESC+"/numLDBins";
    public static final String LD_BINS = LD_DESC+"/binsLD";
    
    private HapMapHDF5Constants() {
        // do not instantiate
    }
    
}
