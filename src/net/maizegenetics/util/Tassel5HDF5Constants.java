/*
 * HDF5Constants
 */
package net.maizegenetics.util;

import ch.systemsx.cisd.hdf5.HDF5FloatStorageFeatures;
import ch.systemsx.cisd.hdf5.HDF5GenericStorageFeatures;
import ch.systemsx.cisd.hdf5.HDF5IntStorageFeatures;

/**
 * Definition of attributes and paths for Tassel HDF5 file format.
 *
 * @author Terry Casstevens
 * @author Ed Buckler
 */
public final class Tassel5HDF5Constants {

    public static final String ROOT = "/";

    // Genotypes Module
    public static final String GENOTYPES_MODULE = "Genotypes";
    public static final String GENOTYPES_ATTRIBUTES_PATH = GENOTYPES_MODULE + "/";
    public static final String GENOTYPES_MAX_NUM_ALLELES = "maxNumAlleles";
    public static final String GENOTYPES_RETAIN_RARE_ALLELES = "retainRareAlleles";
    public static final String GENOTYPES_NUM_TAXA = "numTaxa";
    public static final String GENOTYPES_SCORE_TYPE = "scoreType";
    public static final String GENOTYPES_ALLELE_STATES = GENOTYPES_MODULE + "/AlleleStates";
    public static final String GENO_DESC= GENOTYPES_MODULE + "/_Descriptors/";
    public static final String ALLELE_CNT = GENO_DESC+"AlleleCnt";
    public static final String MAF = GENO_DESC+"MAF";
    public static final String SITECOV = GENO_DESC+"SiteCoverage";
    public static final String ALLELE_FREQ_ORD = GENO_DESC+"AlleleFreqOrder";
    public static final String TAXACOV = GENO_DESC+"TaxaCoverage";
    public static final String TAXAHET = GENO_DESC+"TaxaHet";

    public static final int BLOCK_SIZE=1<<16;


    public static final String getGenotypesPedigreePath(String taxon) {
        return GENOTYPES_MODULE + "/" + taxon + "/pedigree";
    }

    public static final String getGenotypesCallsPath(String taxon) {
        return GENOTYPES_MODULE + "/" + taxon + "/calls";
    }

    public static final String getGenotypesDepthPath(String taxon) {
        return GENOTYPES_MODULE + "/" + taxon + "/depth";
    }

    public static final String getGenotypesScorePath(String taxon) {
        return GENOTYPES_MODULE + "/" + taxon + "/score";
    }

    //Taxa Module
    public static final String TAXA_MODULE = "Taxa";
    public static final String TAXA_ATTRIBUTES_PATH = TAXA_MODULE + "/";
    public static final String TAXA_NUM_TAXA = "numTaxa";

    public static final String getTaxonPath(String taxon) {
        return TAXA_MODULE + "/" + taxon;
    }

    //Position Module
    public static final String POSITION_MODULE = "Positions";
    public static final String POSITION_ATTRIBUTES_PATH = POSITION_MODULE + "/";
    public static final String POSITION_NUM_SITES = "numSites";
    public static final String POSITIONS = POSITION_ATTRIBUTES_PATH + "Positions";
    public static final String CHROMOSOMES = POSITION_ATTRIBUTES_PATH + "Chromosomes";
    public static final String CHROMOSOME_INDICES = POSITION_ATTRIBUTES_PATH + "ChromosomeIndices";
    public static final String SNP_IDS = POSITION_ATTRIBUTES_PATH + "SnpIds";

    //Standard Compression (deflation) levels
    public static final HDF5IntStorageFeatures intDeflation = HDF5IntStorageFeatures.createDeflation(2);
    public static final HDF5GenericStorageFeatures genDeflation = HDF5GenericStorageFeatures.createDeflation(2);
    public static final HDF5FloatStorageFeatures floatDeflation = HDF5FloatStorageFeatures.createDeflation(2);


    private Tassel5HDF5Constants() {
        // do not instantiate
    }

}
