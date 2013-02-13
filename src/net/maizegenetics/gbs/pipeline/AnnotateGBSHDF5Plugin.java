/*
 * BiParentalErrorCorrectionPlugin
 */
package net.maizegenetics.gbs.pipeline;


import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Writer;
import ch.systemsx.cisd.hdf5.IHDF5WriterConfigurator;
import java.awt.Frame;


import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.swing.ImageIcon;

import net.maizegenetics.util.ArgsEngine;
import net.maizegenetics.pal.alignment.MutableNucleotideAlignment;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.AlignmentUtils;
import net.maizegenetics.pal.alignment.BitAlignment;
import net.maizegenetics.pal.alignment.BitAlignmentHDF5;
import net.maizegenetics.pal.alignment.ExportUtils;
import net.maizegenetics.pal.alignment.HapMapHDF5Constants;
import net.maizegenetics.pal.alignment.ImportUtils;
import net.maizegenetics.pal.distance.IBSDistanceMatrix;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.util.BitUtil;

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
    private int[] errorCnt, correctCnt, minorCorrectCnt, testCnt;
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

    private void calculateErrorByIBS(Alignment a,MutableNucleotideAlignment mna, String hdf5File) {
        int taxaNum=a.getSequenceCount();
        int siteCnt=a.getSiteCount();
        int totalWords=a.getAllelePresenceForAllSites(0, 0).getNumWords();
        double[] cvg=new double[a.getSequenceCount()];
        for (int i = 0; i < cvg.length; i++) {
            cvg[i]=(double)a.getTotalNotMissingForTaxon(i)/(double)siteCnt;
        }
        errorCnt=new int[siteCnt];
        correctCnt=new int[siteCnt];
        testCnt=new int[siteCnt];
        minorCorrectCnt=new int[siteCnt];
        errorRate=new float[siteCnt];
        errorsMinorRate=new float[siteCnt];
        int maxErrorKickOut=(int)((double)minSite*0.01*1.5);
        System.out.println("maxErrorKickOut:"+maxErrorKickOut);
        int IBShits=0;
        for (int i = 0; i < taxaNum; i++) {
            long[] iMj = a.getAllelePresenceForAllSites(i, 0).getBits();
            long[] iMn = a.getAllelePresenceForAllSites(i, 1).getBits();
            if(cvg[i]*siteCnt<minSite) continue;
            for (int j = 0; j < i; j++) {
                if(cvg[j]*siteCnt<minSite) continue;
                int estimatedBin=(int)Math.round((double)minSite/(cvg[i]*cvg[j]));
                if(estimatedBin>=siteCnt) continue;
                int binInWords=estimatedBin/64;
                long[] jMj = a.getAllelePresenceForAllSites(j, 0).getBits();
                long[] jMn = a.getAllelePresenceForAllSites(j, 1).getBits();
                for (int bin = 0; bin < totalWords; bin+=binInWords) {
                    int endBlock=bin+binInWords-1;
                    if(endBlock>=totalWords) endBlock=totalWords-1;
                    double dist=computeHetBitDistances(iMj, iMn, jMj, jMn, bin, endBlock, minSite, maxErrorKickOut);
                  //  System.out.printf("%s %g %s %g bin:%d %g %n",a.getTaxaName(i), cvg[i], a.getTaxaName(j), cvg[j], bin, dist);
                    if(dist<IBSthreshold) {//found IBS, score 
                        //System.out.printf("%s %s bin:%d dist:%g %n",a.getTaxaName(i), a.getTaxaName(j), bin, dist);
                        int endSite=(endBlock*64)+63;
                        if(endSite>=siteCnt) endSite=siteCnt-1;
                        recordErrorsForIBSRegion(mna,bin*64,endSite,i,j);
                        IBShits++;
                    }
                }
             
            }
            if(i%100==0) System.out.printf("%s IBSCnt: %d %n",a.getTaxaName(i), IBShits);
        }
        IHDF5WriterConfigurator config = HDF5Factory.configure(hdf5File);
        myLogger.info("Annotating HDF5 file with IBSError: " + hdf5File);
        IHDF5Writer h5w = config.writer();
        int totalErrors=0, totalCorrect=0, totalCnt=0;
        for (int s = 0; s < siteCnt; s++) {
            totalErrors+=errorCnt[s];
            totalCnt+=testCnt[s];
            totalCorrect+=correctCnt[s];
            if(testCnt[s]==0) {
                errorRate[s]=Float.NaN;
                errorsMinorRate[s]=Float.NaN;
            }
            else if((errorCnt[s]+minorCorrectCnt[s])==0) {
                errorRate[s]=0;
                errorsMinorRate[s]=Float.NaN;
            }
            else{
                errorRate[s]=(float)errorCnt[s]/(float)testCnt[s];
                errorsMinorRate[s]=(float)errorCnt[s]/(float)(errorCnt[s]+(2*minorCorrectCnt[s]));
               // if(errorsMinorRate[s]>(maxErrorRate*10)) errorRate[s]=errorsMinorRate[s];
            }
            System.out.printf("Site:%d %d %d %d %d Error: %g %g %n",s,testCnt[s], correctCnt[s], 
                    minorCorrectCnt[s], errorCnt[s], errorRate[s],errorsMinorRate[s]);
            
        }
        if(!h5w.exists(HapMapHDF5Constants.ERROR_DESC)) h5w.createGroup(HapMapHDF5Constants.ERROR_DESC);
        if(!h5w.exists(HapMapHDF5Constants.IBSERROR_DESC)) h5w.createFloatArray(HapMapHDF5Constants.IBSERROR_DESC, errorRate.length);
        h5w.writeFloatArray(HapMapHDF5Constants.IBSERROR_DESC, errorRate);
        if(!h5w.exists(HapMapHDF5Constants.IBSMINORERROR_DESC)) h5w.createFloatArray(HapMapHDF5Constants.IBSMINORERROR_DESC, errorRate.length);
        h5w.writeFloatArray(HapMapHDF5Constants.IBSMINORERROR_DESC, errorsMinorRate);
        System.out.printf("Error:%d Correct:%d Total:%d %n",totalErrors, totalCorrect, totalCnt);
    }
    
    
    private void recordErrorsForIBSRegion(Alignment a, int startBase, int endBase, int taxon1, int taxon2) {
//        System.out.printf("%d %d %n",startBase, endBase);
//        byte[] b1r=a.getBaseRange(taxon1, startBase, endBase+1);
//        byte[] b2r=a.getBaseRange(taxon2, startBase, endBase+1);
        for (int s = startBase; s <= endBase; s++) {
            byte b1=a.getBase(taxon1, s);
            byte b2=a.getBase(taxon2, s);
            if((b1!=Alignment.UNKNOWN_DIPLOID_ALLELE)&&(b2!=Alignment.UNKNOWN_DIPLOID_ALLELE)) {                
                testCnt[s]++;
                if(b1==b2) {
                    correctCnt[s]++;
                    if(b1==minorGenotype[s]) minorCorrectCnt[s]++;
                } else {errorCnt[s]++;}   
            }
        }
    }
    
    public static double computeHetBitDistances(long[] iMj, long[] iMn, long[] jMj, long[] jMn,
            int startBlock, int endBlock, int minSitesCompared, int maxDiffCnt) {
        int sameCnt = 0, diffCnt = 0, hetCnt = 0;
        for (int x = startBlock; x <= endBlock; x++) {
            long same = (iMj[x] & jMj[x]) | (iMn[x] & jMn[x]);
            long diff = (iMj[x] & jMn[x]) | (iMn[x] & jMj[x]);
            long hets = same & diff;
            sameCnt += BitUtil.pop(same);
            diffCnt += BitUtil.pop(diff);
            hetCnt += BitUtil.pop(hets);
            if(diffCnt>maxDiffCnt) return 1.0;
        }
        double identity = (double) (sameCnt + (hetCnt / 2)) / (double) (sameCnt + diffCnt + hetCnt);
        double dist = 1 - identity;
        int sites = sameCnt + diffCnt - hetCnt;
        if (sites > minSitesCompared) {
            return dist;
        } else {
            return Double.NaN;
          //  return sites;
        }
    }
    
    private Alignment removeHighErrorSites(Alignment a, MutableNucleotideAlignment msa, double maxErrorRate, boolean removeUntested) {
        String REMOVED_STATUS = "Removed";
        String LOG_TEST = "Remove High Error Sites";
        String LOG_TEST_UNTESTED = "Remove High Error Sites Untested";

       // MutableNucleotideAlignment msa = MutableNucleotideAlignment.getInstance(a);
        int sitesWithHighError = 0, untestedSites = 0;
        for (int s = 0; s < a.getSiteCount(); s++) {
            if (Double.isNaN(errorRate[s])) {
                if (removeUntested) {
                    msa.clearSiteForRemoval(s);
                    if (snpLogging != null) {
                        snpLogging.writeEntry(a, s, null, null, LOG_TEST_UNTESTED, REMOVED_STATUS, "NaN", "");
                    }
                    untestedSites++;
                }
            } else if (errorRate[s] > maxErrorRate) {
                msa.clearSiteForRemoval(s);
                if (snpLogging != null) {
                    snpLogging.writeEntry(a, s, null, null, LOG_TEST, REMOVED_STATUS, String.valueOf(errorRate[s]), String.valueOf(maxErrorRate));
                }
                sitesWithHighError++;
            }
        }
        System.out.printf("Initial Sites %d  ErrorRate %g  HighErrorSites %d UntestedSites %d %n", a.getSiteCount(),
                maxErrorRate, sitesWithHighError, untestedSites);
        msa.clean();
        System.out.printf("Final Sites %d  ErrorRate %g  HighErrorSites %d %n", msa.getSiteCount(), maxErrorRate, sitesWithHighError);
        return msa;
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
            String outfile=root+"anno.hmp.h5";
            annotateMAF(a,  outfile);
            annotateIBSError(a,  outfile);
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
    
    private void initDiploids(Alignment a) {
        minorGenotype=new byte[a.getSiteCount()];
        for (int i = 0; i < minorGenotype.length; i++) {
            minorGenotype[i]=AlignmentUtils.getDiploidValue(a.getMinorAllele(i), a.getMinorAllele(i));
        }
    }

    private void annotateIBSError(Alignment oA, String hdfFile) {
        oA.optimizeForTaxa(null);
        initDiploids(oA);
        double realDist = AlignmentFilterByGBSUtils.getErrorRateForDuplicatedTaxa(oA, true, false, false);
        double randomDist = AlignmentFilterByGBSUtils.getErrorRateForDuplicatedTaxa(oA, true, true, false);
        System.out.println("Ratio of RandomToReal:" + randomDist / realDist);
        MutableNucleotideAlignment mna = MutableNucleotideAlignment.getInstance(oA);
        calculateErrorByIBS(oA,mna,hdfFile);
            
        reportPercentilesOfErrorRates();
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
