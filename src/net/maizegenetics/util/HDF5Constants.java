/*
 * HDF5Constants
 */
package net.maizegenetics.util;

/**
 * Definition of attributes and paths for Tassel HDF5 file format.
 *
 * @author Terry Casstevens
 */
public final class HDF5Constants {

    public static final String ROOT = "/";

    // Genotypes Module
    public static final String GENOTYPES_MODULE = "Genotypes";
    public static final String GENOTYPES_ATTRIBUTES_PATH = GENOTYPES_MODULE + "/";
    public static final String GENOTYPES_MAX_NUM_ALLELES = "maxNumAlleles";
    public static final String GENOTYPES_RETAIN_RARE_ALLELES = "retainRareAlleles";
    public static final String GENOTYPES_NUM_TAXA = "numTaxa";
    public static final String GENOTYPES_SCORE_TYPE = "scoreType";
    public static final String GENOTYPES_ALLELE_STATES = GENOTYPES_MODULE + "/AlleleStates";

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

    private HDF5Constants() {
        // do not instantiate
    }

}
