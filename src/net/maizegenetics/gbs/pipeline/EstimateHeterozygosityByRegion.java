/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.gbs.pipeline;

import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.BitAlignmentHDF5;
import net.maizegenetics.pal.alignment.ExportUtils;
import net.maizegenetics.pal.alignment.ImportUtils;
import net.maizegenetics.pal.popgen.HeterozygosityProfiler;
import net.maizegenetics.util.ProgressListener;

/**
 * Estimate Heterozygosity for a region in genome for a combinations of listed taxa
 * @author edbuckler
 */
public class EstimateHeterozygosityByRegion {
    private Alignment tA;

    public EstimateHeterozygosityByRegion(String unImpTargetFile, String exportFile) {
        tA=BitAlignmentHDF5.getInstance(unImpTargetFile); 
        System.out.println(tA.toString());
    }
    
    
    
    
    /**
     *
     * @param args
     */
    public static void main(String[] args) {
      
      String root="/Users/edbuckler/SolexaAnal/GBS/build20120701/06_HapMap/";
//        String root="/Volumes/LaCie/build20120110/imp/";
       // String origHMFile=root+"SNP55K_maize282_AGPv2_chr10_20100513.hmp.txt.gz";
        String origHMFile=root+"rje22_BPEC_AllZea_GBS_Build_July_2012_RC-2_chr10.hmp.txt.gz";
        String origH5File=root+"rje22_BPEC_AllZea_GBS_Build_July_2012_RC-2_chr10.hmp.h5";
 //       String donorFile=root+"DTMAfounder20120110.imp.hmp.txt";
        String unImpTargetFile=root+"ZeaSyn20120110.hmp.txt";
        String impTargetFile=root+"ZeaSyn20120110.imp.hmp.txt";

        boolean buildInput=false;
        boolean filterTrue=true;
        if(buildInput) {
            System.out.println("Reading:"+origHMFile);
            Alignment hA=ImportUtils.readFromHapmap(origHMFile, false, (ProgressListener)null);
            System.out.println("Completed Reading:"+origHMFile);
            System.out.println("Writing:"+origH5File);
            ExportUtils.writeToHDF5(hA, origH5File);
            System.out.println("Completed Writing:"+origH5File);
        }

  //      HeterozygosityProfiler hp=new HeterozygosityProfiler(unImpTargetFile, impTargetFile);

        HeterozygosityProfiler hp=new HeterozygosityProfiler(origH5File, impTargetFile);

    }
    
}
