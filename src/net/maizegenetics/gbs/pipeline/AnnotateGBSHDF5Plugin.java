/*
 * BiParentalErrorCorrectionPlugin
 */
package net.maizegenetics.gbs.pipeline;


import cern.colt.list.IntArrayList;
import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Writer;
import ch.systemsx.cisd.hdf5.IHDF5WriterConfigurator;
import java.awt.Frame;


import java.util.ArrayList;
import java.util.Arrays;

import javax.swing.ImageIcon;

import net.maizegenetics.util.ArgsEngine;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.AlignmentUtils;
import net.maizegenetics.pal.alignment.BitAlignment;
import net.maizegenetics.pal.alignment.BitAlignmentHDF5;
import net.maizegenetics.pal.alignment.ExportUtils;
import net.maizegenetics.pal.alignment.FilterAlignment;
import net.maizegenetics.pal.alignment.HapMapHDF5Constants;
import net.maizegenetics.pal.alignment.ImportUtils;
import net.maizegenetics.pal.alignment.Locus;
import net.maizegenetics.pal.popgen.IBDErrorByTaxon;
import net.maizegenetics.pal.popgen.LinkageDisequilibrium;
import net.maizegenetics.pal.statistics.FisherExact;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.BitUtil;
import net.maizegenetics.util.OpenBitSet;

import org.apache.log4j.Logger;

/**
 * Tools for characterizing all the SNPs with MAF, LD patterns, mapping distribution, etc.
 *
 * Error rates are bounded away from zero, but adding 0.5 error to all error
 * rates that that were observed to be zero.
 *
 * @author edbuckler
 */
public class AnnotateGBSHDF5Plugin extends AbstractPlugin {
    private int minSite=500;
    private double IBSthreshold=0.01;
    private float[] errorRate, errorsMinorRate;
    //should we focus on minor allele agreement??
    private int[] errorCnt, mjCorrCnt, mnCorrCnt;
    private double maxErrorRate = 0.05;
    private boolean myRemoveUntestedError = true;
    private String outHapMap = null;
    private static ArgsEngine engine = new ArgsEngine();
    int start = 1, end = 1;
    private String infile;
    private static final Logger myLogger = Logger.getLogger(AnnotateGBSHDF5Plugin.class);
    int homoweight = 1;  //weight of homozygous genotypes to heterozgyous genotypes
    //1 assumes low coverage and each homozygote is only counted once, while hets get a count for each allele
    //2 assumes better coverage, and homozygotes are counted twice.
    private String snpLogFileName;
    private SNPLogging snpLogging = null;
    private byte[] minorGenotype;

    public AnnotateGBSHDF5Plugin() {
        super(null, false);
    }

    public AnnotateGBSHDF5Plugin(Frame parentFrame) {
        super(parentFrame, false);
    }
  
    public void reportPercentilesOfErrorRates() {
        //TODO make it saveable to file or stdout
        float[] le = Arrays.copyOf(errorRate, errorRate.length);
        Arrays.sort(le);
        System.out.println("Percentile\tErrorRate");
        for (int i = 0; i < le.length; i += (le.length / 20)) {
            System.out.printf("%.2g\t%.3g%n", ((double) i / (double) le.length), le[i]);
        }
    }
    
   public static void reportPercentilesOfErrorRates(double[] arr, int intervals) {
        //TODO make it saveable to file or stdout
        double[] le = Arrays.copyOf(arr, arr.length);
        Arrays.sort(le);
        System.out.println("Percentile\tErrorRate");
        for (int i = 0; i < le.length; i += (le.length / intervals)) {
            System.out.printf("%.2g\t%.3g%n", ((double) i / (double) le.length), le[i]);
        }
    }

    @Override
    public void setParameters(String[] args) {
        if (args.length == 0) {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }

        engine.add("-hmp", "-hmpFile", true);
        engine.add("-o", "--outFile", true);
        engine.add("-oE", "--outErrorFile", true);
        engine.add("-oB", "--outBinDistFile", true);
        engine.add("-sC", "--startChrom", true);
        engine.add("-eC", "--endChrom", true);
        engine.add("-mxE", "--maxError", true);
        engine.add("-kpUT", "--keepUntested", false);
        engine.add("-snpLog", "", true);
        engine.parse(args);
        if (engine.getBoolean("-sC")) {
            start = Integer.parseInt(engine.getString("-sC"));
        }
        if (engine.getBoolean("-eC")) {
            end = Integer.parseInt(engine.getString("-eC"));
        }
        infile = engine.getString("-hmp");
        if (engine.getBoolean("-snpLog")) {
            snpLogFileName = engine.getString("-snpLog");
        }
        snpLogging = new SNPLogging(snpLogFileName, this.getClass());
        performFunction(null);
    }

    public void setOutHapMap(String outHapMap) {
        this.outHapMap = outHapMap;
    }

    public void setMaxErrorRate(double maxErrorRate) {
        this.maxErrorRate = maxErrorRate;
    }


    private void printUsage() {
        myLogger.info(
                "\n\n\nAvailable options for the BiParentalErrorCorrectionPlugin are as follows:\n"
                + "-hmp   Input HapMap file\n"
                + "-o     Output HapMap file\n"
                + "-sC    Start chromosome\n"
                + "-eC    End chromosome\n"
                + "-mxE    Maximum error\n"
                + "-kpUT   Keep untested SNPs for error (default remove)\n"
                + "-snpLog  SNPs Removed Log file name\n\n");
    }

    @Override
    public DataSet performFunction(DataSet input) {
        while (start <= end) {
            int chr = start;
            ArrayList<Datum> dList = new ArrayList<Datum>();
            String currFile = infile.replace("+", "" + chr);
            System.out.println("Reading: " + currFile);
            BitAlignmentHDF5 align;
            Alignment a;
            try {
             //   align = BitAlignmentHDF5.getInstance(currFile);
                a=ImportUtils.readFromHapmap(infile, null);

            } catch (Exception e) {
                myLogger.info("Could not read input hapmap file for chr" + chr + ":\n\t" + currFile + "\n\tSkipping...");
                continue;
            }
            System.out.println("Finished Reading: " + currFile);
            String root="/Users/edbuckler/SolexaAnal/GBS/build20120701/06_HapMap/";
            String outfile=root+"anno2.hmp.h5";
            annotateMAF(a,  outfile);
            annotateIBSError(a,  outfile);
            annotateLD(a,outfile,3);
            start++;
        }
        snpLogging.close();
        return null;
    }
    
    private void annotateMAF(Alignment ah5, String hdf5File) {
        float[] maf=new float[ah5.getSequenceCount()];
        float[] hets=new float[ah5.getSequenceCount()];
        float[] scov=new float[ah5.getSequenceCount()];
        float taxa=(float)ah5.getSequenceCount();
        for (int i = 0; i < maf.length; i++) {
            maf[i]=(float)ah5.getMinorAlleleFrequency(i);
            int siteCovCnt=ah5.getTotalGametesNotMissing(i)/2;
            scov[i]=(float)siteCovCnt/taxa;
            hets[i]=(float)ah5.getHeterozygousCount(i)/siteCovCnt;  //terry's proportion hets is off all sites
        }
        IHDF5WriterConfigurator config = HDF5Factory.configure(hdf5File);
        myLogger.info("Annotating HDF5 file: " + hdf5File);
//        config.overwrite();
//        config.dontUseExtendableDataTypes();
        IHDF5Writer h5w = config.writer();
        if(!h5w.exists(HapMapHDF5Constants.SITE_DESC)) h5w.createGroup(HapMapHDF5Constants.SITE_DESC);
        if(!h5w.exists(HapMapHDF5Constants.MAF_DESC)) h5w.createFloatArray(HapMapHDF5Constants.MAF_DESC, maf.length);
        h5w.writeFloatArray(HapMapHDF5Constants.MAF_DESC, maf);
        if(!h5w.exists(HapMapHDF5Constants.HET_DESC)) h5w.createFloatArray(HapMapHDF5Constants.HET_DESC, hets.length);
        h5w.writeFloatArray(HapMapHDF5Constants.HET_DESC, hets);
        if(!h5w.exists(HapMapHDF5Constants.SITECOV_DESC)) h5w.createFloatArray(HapMapHDF5Constants.SITECOV_DESC, scov.length);
        h5w.writeFloatArray(HapMapHDF5Constants.SITECOV_DESC, scov);
    }
    
    private void annotateLD(Alignment a, String hdf5File, int minMinorCnt) {
        FisherExact myFisherExact=new FisherExact(a.getSequenceCount() + 10);
        a.optimizeForSites(null);
        int minCnt=20;
        int sites=a.getSiteCount();
        int binOf10=8;
        float[][] maxR2=new float[sites][binOf10];
        float[][] minP=new float[sites][binOf10];
        int[][] testsLD=new int[sites][binOf10];
        for (int i = 0; i < a.getSiteCount(); i++) {
            Arrays.fill(maxR2[i], Float.NaN);
            Arrays.fill(minP[i], Float.NaN);
            int[][] contig = new int[2][2];
            Locus myLocus=a.getLocus(i);
            int leftSite=Math.abs(a.getSiteOfPhysicalPosition(a.getPositionInLocus(i)-1000000,myLocus));
            int rightSite=Math.abs(a.getSiteOfPhysicalPosition(a.getPositionInLocus(i)+1000000,myLocus));
            BitSet rMj = a.getAllelePresenceForAllTaxa(i, 0);
            BitSet rMn = a.getAllelePresenceForAllTaxa(i, 1); 
            if(rMn.cardinality()<minMinorCnt) continue;
            for (int j = leftSite; j < rightSite; j++) {
                if((j<0)||(j>=a.getSiteCount())||(i==j)) continue;
                int dist=Math.abs(a.getPositionInLocus(i) - a.getPositionInLocus(j));
                int bin=(Integer.numberOfTrailingZeros(Integer.highestOneBit(dist))-7)/2;
                if(bin<0) bin=0;
               // System.out.printf("Dist: %d %d bin: %d %n",dist, Integer.highestOneBit(dist), bin);
                BitSet cMj = a.getAllelePresenceForAllTaxa(j, 0);
                BitSet cMn = a.getAllelePresenceForAllTaxa(j, 1);
                int n = 0;
                n += contig[1][1] = (int) OpenBitSet.intersectionCount(rMn, cMn);
                n += contig[1][0] = (int) OpenBitSet.intersectionCount(rMn, cMj);
                if(contig[1][0]+contig[1][1]<minMinorCnt) continue;
                n += contig[0][1] = (int) OpenBitSet.intersectionCount(rMj, cMn);
                if(contig[0][1]+contig[1][1]<minMinorCnt) continue;
                n += contig[0][0] = (int) OpenBitSet.intersectionCount(rMj, cMj);
                if(n<minCnt) continue;
                double rValue = LinkageDisequilibrium.calculateRSqr(contig[0][0], contig[1][0], contig[0][1], contig[1][1], minCnt);
                if (Double.isNaN(rValue)) {
                    continue;
                }
                double pValue=myFisherExact.getTwoTailedP(contig[0][0], contig[1][0], contig[0][1], contig[1][1]);
                testsLD[i][bin]++;
                if (Float.isNaN(maxR2[i][bin])||(rValue > maxR2[i][bin])) {
                    maxR2[i][bin] = (float)rValue;
                }
                if (Float.isNaN(minP[i][bin])||(pValue < minP[i][bin])) {
                    minP[i][bin] = (float)pValue;
                }
            }
            System.out.printf("s:%d %s %s %s %n",i,Arrays.toString(testsLD[i]),Arrays.toString(maxR2[i]),Arrays.toString(minP[i]));
        }
    }
    
    private void initDiploids(Alignment a) {
        minorGenotype=new byte[a.getSiteCount()];
        for (int i = 0; i < minorGenotype.length; i++) {
            minorGenotype[i]=AlignmentUtils.getDiploidValue(a.getMinorAllele(i), a.getMinorAllele(i));
        }
    }

    private void annotateIBSError(Alignment a, String hdf5File) {
        a.optimizeForTaxa(null);
        int sites=a.getSiteCount();
        errorCnt=new int[sites];
        mjCorrCnt=new int[sites];
        mnCorrCnt=new int[sites];
        for (int bt = 0; bt < a.getSequenceCount(); bt++) {
            IBDErrorByTaxon iebt=new IBDErrorByTaxon(bt,a,50, 1000, 50,0.02);
            mjCorrCnt=addTwoVector(mjCorrCnt,iebt.getMajorCorrectCounts());
            mnCorrCnt=addTwoVector(mnCorrCnt,iebt.getMinorCorrectCounts());
            errorCnt=addTwoVector(errorCnt,iebt.getErrorCounts());
        }
        IHDF5WriterConfigurator config = HDF5Factory.configure(hdf5File);
        myLogger.info("Annotating HDF5 file: " + hdf5File);
//        config.overwrite();
//        config.dontUseExtendableDataTypes();
        IHDF5Writer h5w = config.writer();
        if(!h5w.exists(HapMapHDF5Constants.SITE_DESC)) h5w.createGroup(HapMapHDF5Constants.SITE_DESC);
        if(!h5w.exists(HapMapHDF5Constants.IBSMAJORCORR_DESC)) h5w.createIntArray(HapMapHDF5Constants.IBSMAJORCORR_DESC, mjCorrCnt.length);
        h5w.writeIntArray(HapMapHDF5Constants.IBSMAJORCORR_DESC, mjCorrCnt);
        if(!h5w.exists(HapMapHDF5Constants.IBSMINORCORR_DESC)) h5w.createIntArray(HapMapHDF5Constants.IBSMINORCORR_DESC, mnCorrCnt.length);
        h5w.writeIntArray(HapMapHDF5Constants.IBSMINORCORR_DESC, mnCorrCnt);
        if(!h5w.exists(HapMapHDF5Constants.IBSERROR_DESC)) h5w.createIntArray(HapMapHDF5Constants.IBSERROR_DESC, errorCnt.length);
        h5w.writeIntArray(HapMapHDF5Constants.IBSERROR_DESC, errorCnt);
        float[] errorRate=new float[sites]; 
        float[] mnErrorRate=new float[sites];
        for (int i = 0; i < errorCnt.length; i++) {
            System.out.printf("%d %d %d %d %n",i,mjCorrCnt[i],mnCorrCnt[i],errorCnt[i]);
            int tests=mjCorrCnt[i]+mnCorrCnt[i]+errorCnt[i];
            if(tests==0) {errorRate[i]=Float.NaN;}
            else {errorRate[i]=(float)errorCnt[i]/(float)tests;}
            if((mnCorrCnt[i]+errorCnt[i])==0) {mnErrorRate[i]=Float.NaN;}
            else {mnErrorRate[i]=(float)errorCnt[i]/(float)((2*mnCorrCnt[i])+errorCnt[i]);}
            System.out.printf("%d %d %d %d %d %g %g %n",i,mjCorrCnt[i],mnCorrCnt[i],errorCnt[i], tests, errorRate[i], mnErrorRate[i]);
        }
        if(!h5w.exists(HapMapHDF5Constants.IBSERRORRATE_DESC)) h5w.createFloatArray(HapMapHDF5Constants.IBSERRORRATE_DESC, errorRate.length);
        h5w.writeFloatArray(HapMapHDF5Constants.IBSERRORRATE_DESC, errorRate);
        if(!h5w.exists(HapMapHDF5Constants.IBSMINORERRORRATE_DESC)) h5w.createFloatArray(HapMapHDF5Constants.IBSMINORERRORRATE_DESC, mnErrorRate.length);
        h5w.writeFloatArray(HapMapHDF5Constants.IBSMINORERRORRATE_DESC, mnErrorRate);
    }
    
    private int sum(int[] v) {
        int s=0;
        for (int i : v) {
            s+=i;
        }
        return s;
    }
    
    private int[] addTwoVector(int[] src, int[] addendum) {
        for (int i = 0; i < addendum.length; i++) {
            src[i]+=addendum[i];
        }
        return src;
    }

    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "FilterErrorForBiparental";
    }

    @Override
    public String getToolTipText() {
        return "Filters and estimates error rates from biparental populations";
    }
    
    public static void main(String[] args) {
        String root="/Users/edbuckler/SolexaAnal/GBS/build20120701/06_HapMap/";
        String infile=root+"Z0NE00N_chr10S.hmp.txt.gz";
 //       String infile=root+"All_chr10S.hmp.txt.gz";
        String outfile=root+"anno.hmp.h5";
//        Alignment a = ImportUtils.readFromHapmap(infile, null);
//        ExportUtils.writeToHDF5(a, outfile);
        
 //       infile=outfile;

        String[] args2 = new String[]{
            "-hmp", infile,
            "-o", outfile,
            "-mxE", "0.01",
            
//            "-sC", "10", // Start chromosome
//            "-eC", "10" // End chromosome
        };

        AnnotateGBSHDF5Plugin plugin = new AnnotateGBSHDF5Plugin();
        plugin.setParameters(args2);
        plugin.performFunction(null);
    }
}
