/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.gbs.pipeline;

import net.maizegenetics.baseplugins.ExtractHapmapSubsetPlugin;

/**
 *
 * @author terry
 */
public class JeffPipelines {

    public static void main(String[] args) {
        runTagsToSNPByAlignmentPlugin();
//        runExtractHapmapSubsetPlugin();
//        runCompareGenosBetweenHapMapFilesPlugin();
    }

    public static void runTagsToSNPByAlignmentPlugin() {

        String baseDirMDPLowVol = "/Users/jcg233/Documents/GBS/MDP1_low_vol/";
        String[] MDPLowVolArgs = new String[]{
            "-i", baseDirMDPLowVol + "C08L7ACXX_6_min2.tbt.byte",
            "-y", // use TagsByTaxaByte
            "-o", baseDirMDPLowVol + "tassel4/hapmap/vcfRef/MDP1_low_vol_test9.chr+.hmp.txt",
            "-vcf",
            "-m", baseDirMDPLowVol + "MGP1_low_vol_min2_wPosit.topm.bin",
//            "-mUpd", baseDir+"",
//            "-ref", baseDirMDPLowVol + "maize_agp_v2_chr10.fasta",
            "-ref", "/Users/jcg233/Documents/GBS/refGenome/ZmB73_RefGen_v2.fa",
//            "-LocusBorder", "150",
            "-p", baseDirMDPLowVol + "MDP1_betterFake_ped.txt", 
            "-mnF", "0.8",
            "-mnMAF", "0.005",
            "-mnMAC", "99999", // this will never be satified: this way -mnMAF overrides it
            "-mnLCov", "0.1", // Minimum locus coverage (proportion of Taxa)
//            "-inclRare",
//            "-inclGaps",  // Include sites where major or minor allele is a GAP
//            "-callBiSNPsWGap",  //call sites with a biallelic SNP plus a gap (e.g., A/C/-)
            "-sC", "10", // Start chromosome
            "-eC", "10" // End chromosome
        };

        String baseDirJulyBuild = "/Users/jcg233/Documents/GBS/";
        String[] JulyBuildArgs = new String[]{
            "-i", baseDirJulyBuild + "20120701Build/tbt/mergeTBTHDF5_mergedtaxa_pivot_20120921.h5",
            "-m", baseDirJulyBuild + "20120701Build/topm/AllZeaMasterTags_c10_20120703.topm",
            "-p", baseDirJulyBuild + "20120701Build/AllZeaPedigree2012oct01C.txt",
            "-o", baseDirJulyBuild + "20120701Build/hapmap",
            "-ref", baseDirJulyBuild + "refGenome/ZmB73_RefGen_v2.fa",
            "-mnF", "0.8",
            "-mnMAF", "0.001",
            "-mnMAC", "10",
            "-mnLCov", "0.1", // Minimum locus coverage (proportion of Taxa)
            "-mxSites", "2000",
//            "-y", // use TagsByTaxaByte
//            "-mUpd", baseDir+"",
//            "-LocusBorder", "150",
//            "-inclRare",
//            "-inclGaps",  // Include sites where major or minor allele is a GAP
//            "-callBiSNPsWGap",  //call sites with a biallelic SNP plus a gap (e.g., A/C/-)
            "-sC", "7", // Start chromosome
            "-eC", "7" // End chromosome
        };

        TagsToSNPByAlignmentPlugin plugin = new TagsToSNPByAlignmentPlugin();
        plugin.setParameters(MDPLowVolArgs);
        plugin.performFunction(null);

    }
    
    public static void runExtractHapmapSubsetPlugin() {
        String baseDir = "/Volumes/nextgen/Zea/build20120701/06_HapMap/RC2/04_BPECFilteredSNPs/";
        String outDir =  "/Users/jcg233/Documents/GBS/ShilpaNIL28FMJuly2012BuildRC2BPEC/";
//        for (int chr=1; chr<11; chr++) {
//            String[] args = new String[]{
//                "-h", baseDir+"rje22_BPEC_AllZea_GBS_Build_July_2012_RC-2_chr"+chr+".hmp.txt.gz",
//                "-o", outDir+"ShilpaNIL28FMJuly2012BuildRC2BPECOrig_chr"+chr+".hmp.txt.gz",
//                "-p", outDir+"ShilpaNIL28FMSamples20120701buildUTF8.txt"
//            };
//            ExtractHapmapSubsetPlugin plugin = new ExtractHapmapSubsetPlugin(null);
//            plugin.setParameters(args);
//            plugin.performFunction(null);
//        }

//        baseDir =  "/Users/jcg233/Documents/GBS/ShilpaNIL28FMJuly2012BuildRC2BPEC/";
//        for (int chr=1; chr<11; chr++) {
//            String[] args = new String[]{
//                "-h", baseDir+"ShilpaNIL28FMJuly2012BuildRC2BPECOrig_chr"+chr+".hmp.txt.gz",
//                "-o", baseDir+"ShilpaNIL28FMJuly2012BuildRC2BPECOrig216taxa_chr"+chr+".hmp.txt.gz",
//                "-p", baseDir+"ShilpaNIL28FMSamples20120701build216taxa.txt"
//            };
//            ExtractHapmapSubsetPlugin plugin = new ExtractHapmapSubsetPlugin(null);
//            plugin.setParameters(args);
//            plugin.performFunction(null);
//        }

//        baseDir = "/Volumes/nextgen/Zea/build20120701/06_HapMap/RC2.1/04_BPECFilteredSNPs/";
//        outDir =  "/Users/jcg233/Documents/GBS/ZakFMJuly2012BuildRC2BPEC/";
//        for (int chr=1; chr<11; chr++) {
//            String[] args = new String[]{
//                "-h", baseDir+"AllTaxa_BPEC_AllZea_GBS_Build_July_2012_RC-2.1_chr"+chr+".hmp.txt.gz",
//                "-o", outDir+"ZakKRNCulmFMJuly2012BuildRC2-1BPEC_chr"+chr+".hmp.txt.gz",
//                "-p", outDir+"ZakKRNCulmFMSamples01072012Build.txt"
//            };
//            ExtractHapmapSubsetPlugin plugin = new ExtractHapmapSubsetPlugin(null);
//            plugin.setParameters(args);
//            plugin.performFunction(null);
//        }

        baseDir = "/Users/jcg233/Documents/GBS/20120701BuildRC2-1BPEC/";
        outDir =  "/Users/jcg233/Documents/GBS/ZakFMJuly2012BuildRC2BPEC/";
        for (int chr=5; chr<11; chr++) {
            String[] args = new String[]{
                "-h", baseDir+"AllTaxa_BPEC_AllZea_GBS_Build_July_2012_RC-2.1_chr"+chr+".hmp.txt.gz",
                "-o", outDir+"ZakKRNCulmFMHighCovBC2S3July2012BuildRC2-1BPEC_chr"+chr+".hmp.txt.gz",
                "-p", outDir+"ZakKRNCulmFMHighCovBC2S3Samples01072012Build.txt",
                "-a", "2" // at least 2 "alleles" (actually, genotypes) = polymorphic
            };
            ExtractHapmapSubsetPlugin plugin = new ExtractHapmapSubsetPlugin(null);
            plugin.setParameters(args);
            plugin.performFunction(null);
        }
    }
    
    public static void runCompareGenosBetweenHapMapFilesPlugin() {
        String baseDir = "/Users/jcg233/Documents/GBS/MDP1_low_vol/";
        String[] tassel4WRefVsTassel3Args = new String[]{
            "-hmp1", baseDir+"tassel4/hapmap/withRef/MDP1_low_vol_mergeSNPs.c+.hmp.txt",
//            "-hmp2", baseDir+"hapmap/MDP1_low_vol_noRefOption_mergeSNPs.c+.hmp.txt",
            "-hmp2", baseDir+"hapmap/MDP1_low_vol_RefOptionWOutput2_mergeSNPs.c+.hmp.txt",
            "-sC",   "10",
            "-eC",   "10",
            "-syn",  baseDir+"282synsWRef.txt",
            "-o",    baseDir+"tassel4WRefVsTassel3WRefGenoCompare.txt",
        };

        CompareGenosBetweenHapMapFilesPlugin plugin = new CompareGenosBetweenHapMapFilesPlugin();
        plugin.setParameters(tassel4WRefVsTassel3Args);
        plugin.performFunction(null);
    }
}
