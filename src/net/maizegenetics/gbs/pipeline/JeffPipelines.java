/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.gbs.pipeline;

/**
 *
 * @author terry
 */
public class JeffPipelines {

    public static void main(String[] args) {
        runTagsToSNPByAlignmentPlugin();
    }

    public static void runTagsToSNPByAlignmentPlugin() {

        String baseDirMDPLowVol = "/Users/jcg233/Documents/GBS/MDP1_low_vol/";
        String[] MDPLowVolArgs = new String[]{
            "-i", baseDirMDPLowVol + "C08L7ACXX_6_min2.tbt.byte",
            "-y", // use TagsByTaxaByte
            "-o", baseDirMDPLowVol + "hapmap/testLocusLog",
            "-m", baseDirMDPLowVol + "MGP1_low_vol_min2_wPosit.topm.bin",
            //            "-mUpd", baseDir+"",
            //"-ref", "maize_agp_v2.fasta",
            //"-LocusBorder", "150",
            "-mnF", "0.8",
            "-mnMAF", "0.005",
            "-mnMAC", "99999", // this will never be satified: this way -mnMAF overrides it
            "-mnLCov", "0.1", // Minimum locus coverage (proportion of Taxa)
            //            "-inclGaps",  // Include sites where major or minor allele is a GAP
            //            "-callBiSNPsWGap",  //call sites with a biallelic SNP plus a gap (e.g., A/C/-)
            "-s", "10", // Start chromosome
            "-e", "10" // End chromosome
        };

        String baseDirJulyBuild = "/Volumes/nextgen/Zea/build20120701/";
        String[] JulyBuildArgs = new String[]{
            "-i", baseDirJulyBuild + "05_TBT/04_PivotMergedTaxaTBT/mergedTBTHDF5_mergedtaxa_pivot_20120628.h5",
//            "-y", // use TagsByTaxaByte
            "-o", "/Users/jcg233/Documents/GBS/JulyBuild/hapmap",
            "-m", baseDirJulyBuild + "04_TOPM/AllZeaMasterTags_c10_20120703.topm",
            //            "-mUpd", baseDir+"",
            //"-ref", "maize_agp_v2.fasta",
            //"-LocusBorder", "150",
            "-p", "50_KeyFiles/AllZeaPedigree20120730.txt",
            "-mnF", "0.8",
            "-mnMAF", "0.001",
            "-mnMAC", "99999", // this will never be satified: this way -mnMAF overrides it
            "-mnLCov", "0.1", // Minimum locus coverage (proportion of Taxa)
            //            "-inclGaps",  // Include sites where major or minor allele is a GAP
            //            "-callBiSNPsWGap",  //call sites with a biallelic SNP plus a gap (e.g., A/C/-)
            "-s", "10", // Start chromosome
            "-e", "10" // End chromosome
        };

        TagsToSNPByAlignmentPlugin plugin = new TagsToSNPByAlignmentPlugin();
        plugin.setParameters(MDPLowVolArgs);
        plugin.performFunction(null);

    }
}
