package net.maizegenetics.analysis.imputation;
/*
 * MinorWindowViterbiImputationPlugin
 */


import net.maizegenetics.dna.snp.*;
import net.maizegenetics.dna.snp.io.ProjectionAlignmentIO;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.analysis.popgen.DonorHypoth;
import net.maizegenetics.analysis.distance.IBSDistanceMatrix;
import net.maizegenetics.prefs.TasselPrefs;
import net.maizegenetics.util.*;
import net.maizegenetics.util.BitSet;

import org.apache.commons.lang.ArrayUtils;
import org.apache.log4j.Logger;

import javax.swing.*;

import java.awt.*;
import java.io.BufferedOutputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.text.DecimalFormat;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import net.maizegenetics.dna.map.PositionList;

import static net.maizegenetics.dna.snp.GenotypeTable.UNKNOWN_DIPLOID_ALLELE;
import static net.maizegenetics.dna.snp.GenotypeTable.WHICH_ALLELE.Major;
import static net.maizegenetics.dna.snp.GenotypeTable.WHICH_ALLELE.Minor;
import static net.maizegenetics.dna.snp.GenotypeTableUtils.isHeterozygous;
import static net.maizegenetics.dna.snp.NucleotideAlignmentConstants.GAP_DIPLOID_ALLELE;
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTableBuilder;



/**
 * Imputation approach that relies on nearest neighbor searches of defined haplotypes,
 * followed by HMM Viterbi resolution or block-based resolution.
 * In June 2013, it is the best approach for substantially unrelated taxa.
 * <p>
 * The algorithm relies on having a series of donor haplotypes.  If phased haplotypes are
 * already known, they can be used directly.  If they need to be reconstructed, the
 * {@link FindMergeHaplotypesPlugin} can be used to create haplotypes for windows
 * across the genome.
 * <p>
 * Imputation is done one taxon at a time using the donor haplotypes.  The strategy is as
 *  follows:
 * <li> Every 64 sites is considered a block (or a word in the bitset terminology - length of a long).  There is
 * a separate initial search for nearest neighbors centered around each block. The flanks
 * of the scan window vary but always contain the number of minor alleles specified (a key parameter).
 * Minor alleles approximate the information content of the search.
 * <li> Calculate distance between the target taxon and the donor haplotypes for every window, and rank
 * the top 10 donor haplotypes for each 64 site focus block.
 * <li> Evaluate by Viterbi whether 1 or 2 haplotypes will explain all the sites across the
 * donor haplotype window.  If successful, move to the next region for donor haplotypes.
 * <li> Resolve each focus block by nearest neighbors.  If inbred NN are not close enough,
 * then do a hybrid search seeded with one parent being the initial 10 haplotypes.
 * Set bases based on what is resolved first inbred or hybrid.
 *
 *<p>
 * Error rates are bounded away from zero, but adding 0.5 error to all error
 * rates that that were observed to be zero.
 * <p>
 * TODO:
 * <li>Move accuracy to one method outside of setAlignmentWithDonors
 *
 * @author Edward Buckler
 * @author Kelly Swarts
 * @cite
 */
//@Citation("Two papers: Viterbi from Bradbury, et al (in prep) Recombination patterns in maize\n"+
//        "NearestNeighborSearch Swarts,...,Buckler (in prep) Imputation with large genotyping by sequencing data\n")
public class FILLINImputationPlugin extends AbstractPlugin {
    private String hmpFile;
    private String donorFile;
    private String outFileBase;
    private String errFile=null;
    private boolean hybridNN=true;//if true, uses combination mode in focus block, else set don't impute (default is true)
    private int minMinorCnt=20;
    private int minMajorRatioToMinorCnt=10;  //refinement of minMinorCnt to account for regions with all major
    private int maxDonorHypotheses=20;  //number of hypotheses of record from an inbred or hybrid search of a focus block
    private boolean isOutputProjection=false;

    private double maximumInbredError=0.01;  //inbreds are tested first, if too much error hybrids are tested.
    private double maxHybridErrorRate=0.003;
    private int minTestSites=100;  //minimum number of compared sites to find a hit

    //kelly options
    private boolean twoWayViterbi= true;//if true, the viterbi runs in both directions (the longest path length wins, if inconsistencies)
    private double maxHybridErrFocusHomo= .001;//max error rate for discrepacy between two haplotypes for the focus block. it's default is higher because calculating for fewer sites
    private double maxInbredErrFocusHomo= .003;
    private double maxSmashErrFocusHomo= .01;
    private double maxInbredErrFocusHet= .001;//the error rate for imputing one haplotype in focus block for a het taxon
    private double maxSmashErrFocusHet= .01;
    private double hetThresh= 0.02;//threshold for whether a taxon is considered heterozygous
    
    //options for masking and calculating accuracy
    private String maskKeyFile= null;
    private double propSitesMask= .001;
    private GenotypeTable maskKey= null;
    private double[] MAFClass= null;//new double[]{0,.02,.05,.10,.20,.3,.4,.5,1};
    private int[] MAF= null;
    private static String MAFFile;
    private double[][][] mafAll= null;//sam as all, but first array holds MAF category
    private double[][] all= new double[3][5]; //arrays held ("columns"): 0-maskedMinor, 1-maskedHet, 2-maskedMajor; each array ("rows"):0-to minor, 1-to het, 2-to major, 3-unimp, 4-total for known type
        
    
    private GenotypeTable unimpAlign;  //the unimputed alignment to be imputed, unphased
    private int testing=0;  //level of reporting to stdout
    //major and minor alleles can be differ between the donor and unimp alignment
    private boolean isSwapMajorMinor=true;  //if swapped try to fix it

    private boolean resolveHetIfUndercalled=false;//if true, sets genotypes called to a het to a het, even if already a homozygote
    //variables for tracking accuracy
    private int totalRight=0, totalWrong=0, totalHets=0; //global variables tracking errors on the fly
    private int[] siteErrors, siteCorrectCnt, taxonErrors, taxonCorrectCnt;  //error recorded by sites

    private boolean verboseOutput=true;

    //initialize the transition matrix (T1)
    double[][] transition = new double[][] {
            {.999,.0001,.0003,.0001,.0005},
            {.0002,.999,.00005,.00005,.0002},
            {.0002,.00005,.999,.00005,.0002},
            {.0002,.00005,.00005,.999,.0002},
            {.0005,.0001,.0003,.0001,.999}
    };
    //i tried to optimize the transition matrix for the focus blocks, but peter's worked the best
    double[][] transitionFocus = transition;
    //            new double[][] {
//                    {.98,.04,.08,.04,.003},
//                    {.02,.98,.05,.0005,.0002},
//                    {.0002,.0005,.98,.0005,.0002},
//                    {.0002,.0005,.005,.98,.0002},
//                    {.003,.04,.08,.04,.98}
//        };
    //initialize the emission matrix, states (5) in rows, observations (3) in columns
    double[][] emission = new double[][] {
                    {.998,.001,.001},
                    {.6,.2,.2},
                    {.4,.2,.4},
                    {.2,.2,.6},
                    {.001,.001,.998}
    };


    //    int[] donorIndices;
    private static ArgsEngine engine = new ArgsEngine();
    private static final Logger myLogger = Logger.getLogger(FILLINImputationPlugin.class);

    public FILLINImputationPlugin() {
        super(null, false);
    }

    public FILLINImputationPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    /**
     *
     * @param donorFile should be phased haplotypes
     * @param unImpTargetFile sites must match exactly with donor file
     * @param exportFile output file of imputed sites
     * @param minMinorCnt determines the size of the search window, low recombination 20-30, high recombination 10-15
     * @param minTestSites
     * @param minSitesPresent
     * @param maxHybridErrorRate
     * @param isOutputProjection
     * @param imputeDonorFile
     */
    public void runMinorWindowViterbiImputation(String donorFile, String unImpTargetFile,
                                                String exportFile, int minMinorCnt, int minTestSites, int minSitesPresent,
                                                double maxHybridErrorRate, boolean isOutputProjection, boolean imputeDonorFile) {
        long time=System.currentTimeMillis();
        TasselPrefs.putAlignmentRetainRareAlleles(false);
        System.out.println("Retain Rare alleles is:"+TasselPrefs.getAlignmentRetainRareAlleles());
        this.minTestSites=minTestSites;
        this.isOutputProjection=isOutputProjection;
        unimpAlign=ImportUtils.readGuessFormat(unImpTargetFile);
        GenotypeTable[] donorAlign;
        if(donorFile.contains(".gX")) {donorAlign=loadDonors(donorFile, unimpAlign, minTestSites, verboseOutput);}
        else {donorAlign=loadDonors(donorFile, unimpAlign, minTestSites, verboseOutput);}
        if (maskKeyFile!=null) {
            System.out.println("File already masked. Use input key file for calculating accuracy");
            GenotypeTable inMaskKey= ImportUtils.readGuessFormat(maskKeyFile);
            if (inMaskKey.positions()!=unimpAlign.positions()) maskKey= filterKey(inMaskKey, unimpAlign);
            else maskKey= inMaskKey;
            }
        else if (unimpAlign.depth()!=null) unimpAlign= maskFileByDepth(unimpAlign,7, 7);
        else unimpAlign= maskPropSites(unimpAlign,propSitesMask);

        OpenBitSet[][] conflictMasks=createMaskForAlignmentConflicts(unimpAlign, donorAlign, verboseOutput);

        siteErrors=new int[unimpAlign.numberOfSites()];
        siteCorrectCnt=new int[unimpAlign.numberOfSites()];
        taxonErrors=new int[unimpAlign.numberOfTaxa()];
        taxonCorrectCnt=new int[unimpAlign.numberOfTaxa()];
        
        System.out.printf("Unimputed taxa:%d sites:%d %n",unimpAlign.numberOfTaxa(),unimpAlign.numberOfSites());
        System.out.println("Creating mutable alignment");
        Object mna=null;
        if(isOutputProjection) {
            mna=new ProjectionBuilder(ImportUtils.readGuessFormat(donorFile));
        } else {
            if(exportFile.contains("hmp.h5")) {
                mna= GenotypeTableBuilder.getTaxaIncremental(this.unimpAlign.positions(),exportFile);
            }else {
                mna= GenotypeTableBuilder.getTaxaIncremental(this.unimpAlign.positions());
            }

        }
        
        int numThreads = Runtime.getRuntime().availableProcessors();
        double notMissing= 0;
        double het= 0;
        for (int taxon = 0; taxon < unimpAlign.numberOfTaxa(); taxon++) {het+= (double)unimpAlign.heterozygousCountForTaxon(taxon); notMissing+= (double)unimpAlign.totalNonMissingForTaxon(taxon);}
        double avgHet= het/notMissing;
        double sqDif= 0;
        for (int taxon = 0; taxon < unimpAlign.numberOfTaxa(); taxon++) {double x= avgHet-((double)unimpAlign.heterozygousCountForTaxon(taxon)/(double)unimpAlign.totalNonMissingForTaxon(taxon)); sqDif+= (x*x);}
        System.out.println("Average Heterozygosity: "+avgHet+" plus/minus std error: "+Math.sqrt(sqDif/(unimpAlign.numberOfTaxa()-1))/Math.sqrt(unimpAlign.numberOfTaxa()));
        System.out.println("Time to read in files and generate masks: "+((System.currentTimeMillis()-time)/1000)+" sec");
        ExecutorService pool = Executors.newFixedThreadPool(numThreads);
        if (MAFFile!=null && MAFClass!=null) {MAF= readInMAFFile(MAFFile,unimpAlign, MAFClass); mafAll= new double[MAFClass.length][3][5]; System.out.println("Calculating accuracy within supplied MAF categories.");}
        for (int taxon = 0; taxon < unimpAlign.numberOfTaxa(); taxon+=1) {
            byte[] taxonKey= maskKey.genotypeAllSites(maskKey.taxa().indexOf(unimpAlign.taxaName(taxon)));//this is sloppy, but need to figure out how to filter taxon and sites. Key file does have to be in the same ordre as unimputed but need to contain all taxa.
            int[] trackBlockNN= new int[5];//global variable to track number of focus blocks solved in NN search for system out; index 0 is inbred, 1 is viterbi, 2 is smash, 3 is not solved, 4 is total for all modes
            ImputeOneTaxon theTaxon= (((double)unimpAlign.heterozygousCountForTaxon(taxon)/(double)unimpAlign.totalNonMissingForTaxon(taxon))<hetThresh)?
                new ImputeOneTaxon(taxon, donorAlign, taxonKey, minSitesPresent, conflictMasks,imputeDonorFile, mna, trackBlockNN, maxInbredErrFocusHomo, maxHybridErrFocusHomo, maxSmashErrFocusHomo, true):
                    new ImputeOneTaxon(taxon, donorAlign, taxonKey, minSitesPresent, conflictMasks,imputeDonorFile, mna, trackBlockNN, maxInbredErrFocusHet, 0, maxSmashErrFocusHet, false);
            //        theTaxon.run();
            pool.execute(theTaxon);
        }
        pool.shutdown();
        try{
            if (!pool.awaitTermination(48, TimeUnit.HOURS)) {
                System.out.println("processing threads timed out.");
            }
        }catch(Exception e) {
            System.out.println("Error processing threads");
        }
        System.out.println("");
        StringBuilder s=new StringBuilder();
        s.append(String.format("%s %s MinMinor:%d ", donorFile, unImpTargetFile, minMinorCnt));
        System.out.println(s.toString());
        //double errRate=(double)totalWrong/(double)(totalRight+totalWrong);
        //System.out.printf("TotalRight %d  TotalWrong %d TotalHets: %d ErrRateExcHet:%g %n",totalRight, totalWrong, totalHets, errRate);
        
        if(isOutputProjection) {
            ProjectionAlignmentIO.writeToFile(exportFile, ((ProjectionBuilder)mna).build());
        } else {
            GenotypeTableBuilder ab=(GenotypeTableBuilder)mna;
            if(ab.isHDF5()) {
                ab.build();
            } else {
                ExportUtils.writeToHapmap(ab.build(), false, exportFile, '\t', null);
            }
        }
        System.out.printf("%d %g %d %n",minMinorCnt, maximumInbredError, maxDonorHypotheses);
        double runtime= (double)(System.currentTimeMillis()-time)/(double)1000;
        accuracyOut(all, runtime);
        if (MAFClass!=null) accuracyMAFOut(mafAll);
        System.out.println("Time to read in files, impute target genotypes, and calculate accuracy: "+runtime+" seconds");
    }

    private class ImputeOneTaxon implements Runnable{
        int taxon;
        GenotypeTable[] donorAlign;
        byte[] taxonKey;
        int minSitesPresent;
        OpenBitSet[][] conflictMasks;
        boolean imputeDonorFile;
        GenotypeTableBuilder alignBuilder=null;
        ProjectionBuilder projBuilder=null;

        int[] trackBlockNN;//global variable to track number of focus blocks solved in NN search for system out; index 0 is inbred, 1 is viterbi, 2 is smash, 3 is not solved, 4 is total for all modes
        double focusInbredErr; //threshold for one haplotype imputation in focus block mode
        double focusHybridErr; //threshold for Viterbi imputation in focus block mode
        double focusSmashErr; //threshold for haplotype combination in focus block mode
        boolean hetsMiss; //for inbred lines in two haplotype combination, set hets to missing because likely error. for heterozygous, impute estimated hets in focus block mode
        double[][] allOneTaxon= new double[3][5]; //arrays held ("columns"): 0-maskedMinor, 1-maskedHet, 2-maskedMajor; each array ("rows"):0-to minor, 1-to het, 2-to major, 3-unimp, 4-total for known type
        double[][][] mafTaxon;
        
        public ImputeOneTaxon(int taxon, GenotypeTable[] donorAlign, byte[] taxonKey, int minSitesPresent, OpenBitSet[][] conflictMasks,
            boolean imputeDonorFile, Object mna, int[] trackBlockNN, double focusInbErr, double focusHybridErr, double focusSmashErr, boolean hetsToMissing) {
            this.taxon=taxon;
            this.donorAlign=donorAlign;
            this.taxonKey=taxonKey;
            this.minSitesPresent=minSitesPresent;
            this.conflictMasks=conflictMasks;
            this.imputeDonorFile=imputeDonorFile;
            this.trackBlockNN=trackBlockNN;
            this.focusInbredErr=focusInbErr;
            this.focusHybridErr=focusHybridErr;
            this.focusSmashErr=focusSmashErr;
            this.hetsMiss=hetsToMissing;
            if(mna instanceof GenotypeTableBuilder) {alignBuilder=(GenotypeTableBuilder)mna;}
            else if(mna instanceof ProjectionBuilder) {projBuilder=(ProjectionBuilder)mna;}
            else {throw new IllegalArgumentException("Only Aligmnent or Projection Builders may be used.");}
            
        }



        @Override
        public void run() {
            StringBuilder sb=new StringBuilder();
            String name=unimpAlign.taxaName(taxon);
            ImputedTaxon impTaxon=new ImputedTaxon(taxon, unimpAlign.genotypeAllSites(taxon),isOutputProjection);
            boolean het= (focusHybridErr==0)?true:false;
            int[] unkHets=countUnknownAndHets(impTaxon.getOrigGeno());
            sb.append(String.format("Imputing %d:%s AsHet:%b Mj:%d, Mn:%d Unk:%d Hets:%d... ", taxon,name,het,
                    unimpAlign.allelePresenceForAllSites(taxon, Major).cardinality(),
                    unimpAlign.allelePresenceForAllSites(taxon, Minor).cardinality(), unkHets[0], unkHets[1]));
            boolean enoughData=(unimpAlign.totalNonMissingForTaxon(taxon)>minSitesPresent);
//                System.out.println("Too much missing data");
//                continue;
//            }
            int countFullLength=0;
            int countByFocus= 0;
            for (int da = 0; (da < donorAlign.length)&&enoughData ; da++) {
                int donorOffset=unimpAlign.siteOfPhysicalPosition(donorAlign[da].chromosomalPosition(0), donorAlign[da].chromosome(0));
                int blocks=donorAlign[da].allelePresenceForAllSites(0, Major).getNumWords();
                BitSet[] maskedTargetBits=arrangeMajorMinorBtwAlignments(unimpAlign, taxon, donorOffset,
                        donorAlign[da].numberOfSites(),conflictMasks[da][0],conflictMasks[da][1]);
                int[] donorIndices;
                if(imputeDonorFile){
                    donorIndices=new int[donorAlign[da].numberOfTaxa()-1];
                    for (int i = 0; i < donorIndices.length; i++) {donorIndices[i]=i; if(i>=taxon) donorIndices[i]++;}
                } else {
                    donorIndices=new int[donorAlign[da].numberOfTaxa()];
                    for (int i = 0; i < donorIndices.length; i++) {donorIndices[i]=i;}
                }
                DonorHypoth[][] regionHypthInbred=new DonorHypoth[blocks][maxDonorHypotheses];
                calcInbredDist(impTaxon,maskedTargetBits, donorAlign[da]);
                for (int focusBlock = 0; focusBlock < blocks; focusBlock++) {
                    int[] resultRange=getBlockWithMinMinorCount(maskedTargetBits[0].getBits(),maskedTargetBits[1].getBits(), focusBlock, minMinorCnt);
                    if(resultRange==null) continue; //no data in the focus Block
                    //search for the best inbred donors for a segment
                    regionHypthInbred[focusBlock]=getBestInbredDonors(taxon, impTaxon, resultRange[0],resultRange[2], focusBlock, donorAlign[da], donorIndices);
                }
                impTaxon.setSegmentSolved(false);
                //tries to solve the entire donorAlign region with 1 or 2 donor haplotypes
                impTaxon=apply1or2Haplotypes(taxon, donorAlign[da], donorOffset, regionHypthInbred,  impTaxon, maskedTargetBits, maxHybridErrorRate);
                if(impTaxon.isSegmentSolved()) {
//                    System.out.printf("VertSolved da:%d L:%s%n",da, donorAlign[da].getLocus(0));
                    countFullLength++; continue;}
                //resorts to solving block by block, first by inbred, then by viterbi, and then by hybrid
                impTaxon=applyBlockNN(impTaxon, taxon, donorAlign[da], donorOffset, regionHypthInbred, hybridNN, maskedTargetBits,
                            minMinorCnt, focusInbredErr, focusHybridErr, focusSmashErr, donorIndices, trackBlockNN, hetsMiss);
                if(impTaxon.isSegmentSolved()) {
//                    System.out.printf("VertSolved da:%d L:%s%n",da, donorAlign[da].getLocus(0));
                    countByFocus++; continue;
                } else continue;
            }
            double[][] accForTaxon= taxonAccuracy(taxonKey,impTaxon.getImpGeno(),unimpAlign, MAF, MAFClass);
            synchronized (all) {//add per taxon accuracy to total count. This prevents problems with multithreading, I think.
                for (int i = 0; i < accForTaxon.length; i++) { for (int j = 0; j < accForTaxon[i].length; j++) { all[i][j]+= accForTaxon[i][j];}}
            }
            if (MAFClass!=null) synchronized (mafAll) {//same as all, but sorted by MAF class
                for (int i = 0; i < mafTaxon.length; i++) { for (int j = 0; j < mafTaxon[0].length; j++) { for (int k = 0; k < mafTaxon[0][0].length; k++) {
                        mafAll[i][j][k]+= mafTaxon[i][j][k];}}}
            }
            double totalFocus= (double)trackBlockNN[3]+(double)trackBlockNN[4];
            sb.append(String.format("InbredOrViterbi:%d FocusBlock:%d PropFocusInbred:%f PropFocusViterbi:%f PropFocusSmash:%f PropFocusMissing:%f BlocksSolved:%d ",
                    countFullLength, countByFocus, (double)trackBlockNN[0]/totalFocus, (double)trackBlockNN[1]/totalFocus, (double)trackBlockNN[2]/totalFocus, (double)trackBlockNN[3]/totalFocus, impTaxon.getBlocksSolved()));
            //sb.append(String.format("InbredOrViterbi:%d FocusBlock:%d BlocksSolved:%d ", countFullLength, countByFocus, impTaxon.getBlocksSolved()));
            int[] unk=countUnknownAndHets(impTaxon.resolveGeno);
            sb.append(String.format("Unk:%d PropMissing:%g ", unk[0], (double) unk[0] / (double) impTaxon.getOrigGeno().length));
            sb.append(String.format("Het:%d PropHet:%g ", unk[1], (double)unk[1]/(double)impTaxon.getOrigGeno().length));
            if(!isOutputProjection) {
                alignBuilder.addTaxon(unimpAlign.taxa().get(taxon), impTaxon.resolveGeno);
            } else {
                projBuilder.addTaxon(unimpAlign.taxa().get(taxon),impTaxon.getBreakPoints());
            }
//            double rate=(double)taxon/(double)(System.currentTimeMillis()-time);
//            double remaining=(unimpAlign.getSequenceCount()-taxon)/(rate*1000);
//            System.out.printf("TimeLeft:%.1fs %n", remaining);
            if(verboseOutput) System.out.println(sb.toString());
//            if(isOutputProjection) {
//                // System.out.println(breakPoints.toString());
//                ((ProjectionAlignment)mna).setCompositionOfTaxon(taxon, impTaxon.breakPoints);
//            }
        }
        
        public double[][] taxonAccuracy(byte[] key, byte[] imputed, GenotypeTable unimpAlign, int[] MAF, double[] MAFClass) {
            byte diploidN= GenotypeTable.UNKNOWN_DIPLOID_ALLELE;
            double[][] all= new double[3][5];
            boolean use= false; boolean mafOn= false; int maf= -1;
            if (MAF!=null && MAFClass!=null) {mafTaxon= new double[MAFClass.length][3][5]; use= true; mafOn= true;}
            for (int site = 0; site < imputed.length; site++) {
                use= (mafOn && MAF[site] > -1)?true:false;
                if (use) maf= MAF[site];
                byte known = key[site];
                if (known == diploidN) continue;
                byte imp = imputed[site];
                if (GenotypeTableUtils.isHeterozygous(known) == true) {
                    all[1][4]++; if (use) mafTaxon[maf][1][4]++;
                    if (imp == diploidN) {all[1][3]++; if (use) mafTaxon[maf][1][3]++;}
                    else if (GenotypeTableUtils.isEqual(imp, known) == true) {all[1][1]++; if (use) mafTaxon[maf][1][1]++;}
                    else if (GenotypeTableUtils.isHeterozygous(imp) == false && GenotypeTableUtils.isPartiallyEqual(imp, unimpAlign.minorAllele(site)) == true) {all[1][0]++; if (use) mafTaxon[maf][1][0]++;}//to minor 
                    else if (GenotypeTableUtils.isHeterozygous(imp) == false && GenotypeTableUtils.isPartiallyEqual(imp, unimpAlign.majorAllele(site)) == true) {all[1][2]++; if (use) mafTaxon[maf][1][2]++;}
                    else {all[1][4]--;  if (use) mafTaxon[maf][1][4]--;}//implies >2 allele states at given genotype
                } else if (known == GenotypeTableUtils.getDiploidValue(unimpAlign.minorAllele(site),unimpAlign.minorAllele(site))) {
                    all[0][4]++; if (use) mafTaxon[maf][0][4]++;
                    if (imp == diploidN) {all[0][3]++;  if (use) mafTaxon[maf][0][3]++;}
                    else if (GenotypeTableUtils.isEqual(imp, known) == true) {all[0][0]++;  if (use) mafTaxon[maf][0][0]++;}
                    else if (GenotypeTableUtils.isHeterozygous(imp) == true && GenotypeTableUtils.isPartiallyEqual(imp, known) == true) {all[0][1]++;  if (use) mafTaxon[maf][0][1]++;}
                    else {all[0][2]++;  if (use) mafTaxon[maf][0][3]++;}
                } else if (known == GenotypeTableUtils.getDiploidValue(unimpAlign.majorAllele(site),unimpAlign.majorAllele(site))) {
                    all[2][4]++;  if (use) mafTaxon[maf][2][4]++;
                    if (imp == diploidN) {all[2][3]++;  if (use) mafTaxon[maf][2][3]++;}
                    else if (GenotypeTableUtils.isEqual(imp, known) == true) {all[2][2]++; if (use) mafTaxon[maf][2][2]++;}
                    else if (GenotypeTableUtils.isHeterozygous(imp) == true && GenotypeTableUtils.isPartiallyEqual(imp, known) == true) {all[2][1]++;  if (use) mafTaxon[maf][2][1]++;}
                    else {all[2][0]++; if (use) mafTaxon[maf][2][0]++;}
                } else continue;
            }
            return all;
        }

    }

     public static GenotypeTable[] loadDonors(String donorFileRoot, GenotypeTable unimpAlign, int minTestSites, boolean verboseOutput){
        File theDF=new File(donorFileRoot);
        String prefilter=theDF.getName().split(".gX.")[0]+".gc"; //grabs the left side of the file
        String prefilterOld=theDF.getName().split("s\\+")[0]+"s"; //grabs the left side of the file
        ArrayList<File> d=new ArrayList<File>();
        for (File file : theDF.getParentFile().listFiles()) {
            if(file.getName().equals(theDF.getName())) {d.add(file);}
             if(file.getName().startsWith(prefilter)) {d.add(file);}
             if(file.getName().startsWith(prefilterOld)) {d.add(file);}
         }
        ArrayList<Integer> removeDonors= new ArrayList<>();
        PositionList donU= unimpAlign.positions();
        GenotypeTable[] donorAlign=new GenotypeTable[d.size()];
        for (int i = 0; i < donorAlign.length; i++) {
            if(verboseOutput) System.out.println("Starting Read");
            donorAlign[i]=ImportUtils.readFromHapmap(d.get(i).getPath(), (ProgressListener)null);
            ArrayList<Integer> subSites= new ArrayList<>();
            PositionList donP= donorAlign[i].positions();
            for (int j = 0; j < donorAlign[i].numberOfSites(); j++) {if (donU.indexOf(donP.get(j)) > -1) subSites.add(j);} //if unimputed contains donorAlign position keep in donor align
            if (subSites.size()!=donorAlign[i].numberOfSites()){
                if (subSites.size()<2) {
                    if(verboseOutput) System.out.printf("Donor file contains <2 matching sites and will not be used:%s",d.get(i).getPath());
                    removeDonors.add(i);
                }
                else {
                    donorAlign[i]= FilterGenotypeTable.getInstance(donorAlign[i], ArrayUtils.toPrimitive(subSites.toArray(new Integer[subSites.size()])));
                    if(verboseOutput) System.out.printf("Donor file sites filtered to match target:%s taxa:%d sites:%d %n",d.get(i).getPath(), donorAlign[i].numberOfTaxa(),donorAlign[i].numberOfSites());
                }
                if (subSites.size() < minTestSites*2 && verboseOutput) System.out.println("This donor alignment contains marginally sufficient matching snp positions. Region unlikely to impute well.");
            }
            else if(verboseOutput) System.out.printf("Donor file shares all sites with target:%s taxa:%d sites:%d %n",d.get(i).getPath(), donorAlign[i].numberOfTaxa(),donorAlign[i].numberOfSites());
            //createMaskForAlignmentConflicts(unimpAlign,donorAlign[i],true);
        }
        if (removeDonors.isEmpty()==false) for(GenotypeTable gt:donorAlign) {donorAlign= (GenotypeTable[])ArrayUtils.removeElement(donorAlign, gt);}
        MAFFile= donorFileRoot.substring(0, donorFileRoot.indexOf(".gX"))+"MAF.txt";
        if (new File(MAFFile).isFile()==false) generateMAF(MAFFile, donorAlign);
        return donorAlign;
    }
    //outputs a text file with minor allele frequencies by locus and position in the donor files. for accuracy
    public static void generateMAF(String MAFFileName, GenotypeTable[] donorAlign) {
        System.out.println("Generating MAF file for input donor files");
        File outputFile = new File(MAFFile);
        try {
            
            DataOutputStream outStream = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outputFile)));
            DecimalFormat df = new DecimalFormat("0.########");
            for (GenotypeTable don:donorAlign) {
                for (int site = 0; site < don.numberOfSites(); site++) {
                    outStream.writeBytes(don.chromosomeName(site)+"\t"+don.chromosomalPosition(site)+"\t"+df.format(don.minorAlleleFrequency(site))+"\n");
                }
            }
        }
        catch (Exception e) {
            System.out.println("Problem writing MAF file");
        }
    }
    
    //returns an array the length of key.numberOfSites() that contains the MAF class index of the donor file or -1 if not included
    public static int[] readInMAFFile (String donorMAFFile, GenotypeTable unimpAlign, double[] MAFClass) {
        int[] MAF= new int[unimpAlign.numberOfSites()];
        for (int i = 0; i < MAF.length; i++) {MAF[i]= -1;} //initialize to -1
        try {
            FileInputStream fis= new FileInputStream(donorMAFFile);
            Scanner scanner= new Scanner(fis);
            do {
                String next= scanner.nextLine();
                if (next.isEmpty()) break;
                String[] vals= next.split("\t");
                int site= unimpAlign.siteOfPhysicalPosition(Integer.parseInt(vals[1]),unimpAlign.chromosome(vals[0]));
                if (site<0) continue;
                int currClass= Arrays.binarySearch(MAFClass, Double.parseDouble(vals[2]));
                MAF[site]= currClass<0?Math.abs(currClass)-1:currClass;
            }
            while (scanner.hasNextLine());
            scanner.close();
            fis.close();
        }
        catch (Exception e) {
            System.out.println("Problem reading in taxa names");
        }
        return MAF;
    }

    private ImputedTaxon apply1or2Haplotypes(int taxon, GenotypeTable donorAlign, int donorOffset,
                                             DonorHypoth[][] regionHypth, ImputedTaxon impT,
                                             BitSet[] maskedTargetBits, double maxHybridErrorRate) {

        int blocks=maskedTargetBits[0].getNumWords();
        //do flanking search
        if(testing==1) System.out.println("Starting complete hybrid search");
        int[] d=getAllBestDonorsAcrossChromosome(regionHypth,blocks/20);  //TODO
        DonorHypoth[] best2donors=getBestHybridDonors(taxon, maskedTargetBits[0].getBits(),
                maskedTargetBits[1].getBits(), 0, blocks-1, blocks/2, donorAlign, d, d, true);
        if(testing==1) System.out.println(Arrays.toString(best2donors));
        ArrayList<DonorHypoth> goodDH=new ArrayList<DonorHypoth>();
        for (DonorHypoth dh : best2donors) {
            if((dh!=null)&&(dh.getErrorRate()<maxHybridErrorRate)) {
                if(dh.isInbred()==false){
                    dh=getStateBasedOnViterbi(dh, donorOffset, donorAlign, twoWayViterbi, transition);
                }
                if(dh!=null) goodDH.add(dh);
            }
        }
        if(goodDH.isEmpty()) return impT;
        DonorHypoth[] vdh=new DonorHypoth[goodDH.size()];
        for (int i = 0; i < vdh.length; i++) {vdh[i]=goodDH.get(i);}
        impT.setSegmentSolved(true);
        return setAlignmentWithDonors(donorAlign,vdh, donorOffset,false,impT, false, false);
    }

    /**
     * Create mask for all sites where major & minor are swapped in alignments
     * returns in this order [goodMask, swapMjMnMask, errorMask, invariantMask]
     */
    private OpenBitSet[][] createMaskForAlignmentConflicts(GenotypeTable unimpAlign, GenotypeTable[] donorAlign, boolean print) {
        OpenBitSet[][] result=new OpenBitSet[donorAlign.length][4];
        for (int da = 0; da < result.length; da++) {
//            int donorOffset=unimpAlign.siteOfPhysicalPosition(donorAlign[da].chromosomalPosition(0), donorAlign[da].chromosome(0), donorAlign[da].siteName(0));
            int donorOffset=unimpAlign.positions().indexOf(donorAlign[da].positions().get(0));
            OpenBitSet goodMask=new OpenBitSet(donorAlign[da].numberOfSites());
            OpenBitSet swapMjMnMask=new OpenBitSet(donorAlign[da].numberOfSites());
            OpenBitSet errorMask=new OpenBitSet(donorAlign[da].numberOfSites());
            OpenBitSet invariantMask=new OpenBitSet(donorAlign[da].numberOfSites());
            int siteConflicts=0, swaps=0, invariant=0, good=0;
            for (int i = 0; i < donorAlign[da].numberOfSites(); i++) {
                /*we have three classes of data:  invariant in one alignment, conflicts about minor and minor,
                *swaps of major and minor.  Adding the invariant reduces imputation accuracy.
                *the major/minor swaps should be flipped in the comparisons
                */
                byte tMj=unimpAlign.majorAllele(i + donorOffset);
                byte tMn=unimpAlign.minorAllele(i + donorOffset);
                byte daMj=donorAlign[da].majorAllele(i);
                byte daMn=donorAlign[da].minorAllele(i);
                if(daMn==GenotypeTable.UNKNOWN_ALLELE) {
                    invariant++;
                    invariantMask.set(i);
                    goodMask.set(i);
                } else
                if((daMj==tMn)&&(daMn==tMj)) {
                    swaps++;
                    swapMjMnMask.set(i);
                    goodMask.set(i);
                } else
                if((daMj!=tMj)) {
                    siteConflicts++;
                    errorMask.set(i);
                    goodMask.set(i);
                }

            }
            goodMask.not();
            if(print) System.out.println("Donor:"+da+"invariant in donor:"+invariant+" swapConflicts:"+swaps+" errors:"+siteConflicts);
            result[da]=new OpenBitSet[] {goodMask, swapMjMnMask, errorMask, invariantMask};
        }
        return result;
    }

    //This assumes the input key contains all of the sites in unimpAlign, plus potentially more of either
    private GenotypeTable filterKey(GenotypeTable maskKey, GenotypeTable unimpAlign) {
        System.out.println("Filtering user input key file...\nSites in original Key file: "+maskKey.numberOfSites());
        String[] unimpNames= new String[unimpAlign.numberOfSites()];
        for (int site = 0; site < unimpAlign.numberOfSites(); site++) {unimpNames[site]= unimpAlign.siteName(site);}
        ArrayList<String> keepSites= new ArrayList<>();for (int site = 0; site < maskKey.numberOfSites(); site++) {
            if (unimpAlign.positions().indexOf(maskKey.positions().get(site))>-1) keepSites.add(maskKey.siteName(site));
        }
        FilterGenotypeTable filter= FilterGenotypeTable.getInstance(maskKey, keepSites.toArray(new String[keepSites.size()]));
        GenotypeTable newMask= GenotypeTableBuilder.getGenotypeCopyInstance(filter);
        System.out.println("Sites in new mask: "+newMask.numberOfSites());
        return newMask;
    }
//for depths 4 or more, requires hets to be called by more than one for less depth allele. returns 2-array whre index 0 is mask and 1 is key
    private GenotypeTable maskFileByDepth(GenotypeTable a, int depthToMask, int maskDenom) {
        System.out.println("Masking file using depth\nSite depth to mask: "+depthToMask+"Divisor for physical positions to be masked: "+maskDenom);
        GenotypeCallTableBuilder mask= GenotypeCallTableBuilder.getInstance(a.numberOfTaxa(), a.numberOfSites());
        GenotypeCallTableBuilder key= GenotypeCallTableBuilder.getInstance(a.numberOfTaxa(), a.numberOfSites());
        
        int cnt= 0;
        for (int taxon = 0; taxon < a.numberOfTaxa(); taxon++) {
            int taxaCnt= 0;
            mask.setBaseRangeForTaxon(taxon, 0, a.genotypeAllSites(taxon));
            for (int site = 0; site < a.numberOfSites(); site++) {
                if (GenotypeTableUtils.isEqual(NucleotideGenotypeTable.UNKNOWN_DIPLOID_ALLELE, a.genotype(taxon, site))) continue;
                if (a.physicalPositions()[site]%maskDenom!=0) continue;
                int[] currD= a.depthForAlleles(taxon, site);
                if (currD[0]+currD[1]!=depthToMask) continue;
                else if ((a.isHeterozygous(taxon, site)==false) ||
                            (depthToMask > 3 && currD[0] > 1 && currD[1] > 1)|| 
                            (depthToMask < 4)) {
                    mask.setBase(taxon, site, NucleotideGenotypeTable.UNKNOWN_DIPLOID_ALLELE); key.setBase(taxon, site, a.genotype(taxon, site)); taxaCnt++;
                }
            }
            System.out.println(taxaCnt+" sites masked for "+a.taxaName(taxon)); cnt+= taxaCnt;
        }
        System.out.println(cnt+" sites masked at a depth of "+depthToMask+" (site numbers that can be divided by "+maskDenom+")");
        maskKey= GenotypeTableBuilder.getInstance(key.build(), a.positions(), a.taxa());
        return GenotypeTableBuilder.getInstance(mask.build(), a.positions(), a.taxa());
    }
    
    private GenotypeTable maskPropSites(GenotypeTable a, double propSitesMask) {
        System.out.println("Masking file without depth\nMasking "+propSitesMask+" proportion of sites");
        GenotypeCallTableBuilder mask= GenotypeCallTableBuilder.getInstance(a.numberOfTaxa(), a.numberOfSites());
        GenotypeCallTableBuilder key= GenotypeCallTableBuilder.getInstance(a.numberOfTaxa(), a.numberOfSites());
        
        int presGenos= 0;
        for (int taxon = 0; taxon < a.numberOfTaxa(); taxon++) {presGenos+= a.totalNonMissingForTaxon(taxon);}
        int expected= (int)(propSitesMask*(double)presGenos);
        int cnt= 0;
        for (int taxon = 0; taxon < a.numberOfTaxa(); taxon++) {
            int taxaCnt= 0;
            mask.setBaseRangeForTaxon(taxon, 0, a.genotypeAllSites(taxon));
            for (int site = 0; site < a.numberOfSites(); site++) {
                if (Math.random()<propSitesMask && GenotypeTableUtils.isEqual(NucleotideGenotypeTable.UNKNOWN_DIPLOID_ALLELE, a.genotype(taxon, site))==false) {
                    mask.setBase(taxon, site, NucleotideGenotypeTable.UNKNOWN_DIPLOID_ALLELE); key.setBase(taxon, site, a.genotype(taxon, site)); taxaCnt++;
                }
            }
            cnt+= taxaCnt;
        }
        System.out.println(cnt+" sites masked randomly not based on depth ("+expected+" expected at "+propSitesMask+")");
        maskKey= GenotypeTableBuilder.getInstance(key.build(), a.positions(), a.taxa());
        return GenotypeTableBuilder.getInstance(mask.build(), a.positions(), a.taxa());
    }
    
    //this is the sample multiple r2.
    private static double pearsonR2(double[][] all) {
        int size= 0;
        for (int x = 0; x < 3; x++) {size+= (all[x][4]-all[x][3]);}
        double[][] xy= new double[2][size]; //0 is x, 1 is y
        int last= 0;//the last index filled
        for (double x = 0; x < 3; x++) { for (double y = 0; y < 3; y++) {
                for (int fill = last; fill < last+all[(int)x][(int)y]; fill++) {
                    xy[0][fill]= x;
                    xy[1][fill]= y;
                }
                last= last+(int)all[(int)x][(int)y];
            }}
        double meanX= 0; double meanY= 0; double varX= 0; double varY= 0; double covXY= 0; double r2= 0.0;
        for (int i = 0; i < xy[0].length; i++) {meanX+=xy[0][i]; meanY+= xy[1][i];}
        meanX= meanX/(xy[0].length-1); meanY= meanY/(xy[1].length-1);
        double currX, currY;
        for (int i = 0; i < xy[0].length; i++) {
            currX= xy[0][i]-meanX; currY= xy[1][i]-meanY;
            varX+= currX*currX; varY+= currY*currY;
            covXY+= currX*currY;
        }
        r2= (covXY/(Math.sqrt(varX)*Math.sqrt(varY)))*(covXY/(Math.sqrt(varX)*Math.sqrt(varY)));
        System.out.println("Unadjusted R2 value for "+size+" comparisons: "+r2);
        return r2;
    }
    
    private void accuracyOut(double[][] all, double time) {
        DecimalFormat df = new DecimalFormat("0.########");
        double r2= pearsonR2(all);
        try {
            File outputFile = new File(outFileBase.substring(0, outFileBase.indexOf(".hmp")) + "DepthAccuracy.txt");
            DataOutputStream outStream = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outputFile)));
            outStream.writeBytes("##Taxon\tTotalSitesMasked\tTotalSitesCompared\tTotalPropUnimputed\tNumMinor\tCorrectMinor\tMinorToHet\tMinorToMajor\tUnimpMinor"
                    + "\tNumHets\tHetToMinor\tCorrectHet\tHetToMajor\tUnimpHet\tNumMajor\tMajorToMinor\tMajorToHet\tCorrectMajor\tUnimpMajor\tR2\n");
            outStream.writeBytes("##TotalByImputed\t"+(all[0][4]+all[1][4]+all[2][4])+"\t"+(all[0][4]+all[1][4]+all[2][4]-all[0][3]-all[1][3]-all[2][3])+"\t"+
                    ((all[0][3]+all[1][3]+all[2][3])/(all[0][4]+all[1][4]+all[2][4]))+"\t"+all[0][4]+"\t"+all[0][0]+"\t"+all[0][1]+"\t"+all[0][2]+"\t"+all[0][3]+
                    "\t"+all[1][4]+"\t"+all[1][0]+"\t"+all[1][1]+"\t"+all[1][2]+"\t"+all[1][3]+"\t"+all[2][4]+"\t"+all[2][0]+"\t"+all[2][1]+"\t"+all[2][2]+
                    "\t"+all[2][3]+"\t"+r2+"\n");
            outStream.writeBytes("#Minor=0,Het=1,Major=2;x is masked(known), y is predicted\nx\ty\tN\tprop\n"
                    +0+"\t"+0+"\t"+all[0][0]+"\t"+df.format((all[0][0])/(all[0][0]+all[0][1]+all[0][2]))+"\n"
                    +0+"\t"+.5+"\t"+all[0][1]+"\t"+df.format((all[0][1])/(all[0][0]+all[0][1]+all[0][2]))+"\n"
                    +0+"\t"+1+"\t"+all[0][2]+"\t"+df.format((all[0][2])/(all[0][0]+all[0][1]+all[0][2]))+"\n"
                    +.5+"\t"+0+"\t"+all[1][0]+"\t"+df.format((all[1][0])/(all[1][0]+all[1][1]+all[1][2]))+"\n"
                    +.5+"\t"+.5+"\t"+all[1][1]+"\t"+df.format((all[1][1])/(all[1][0]+all[1][1]+all[1][2]))+"\n"
                    +.5+"\t"+1+"\t"+all[1][2]+"\t"+df.format((all[1][2])/(all[1][0]+all[1][1]+all[1][2]))+"\n"
                    +1+"\t"+0+"\t"+all[2][0]+"\t"+df.format((all[2][0])/(all[2][0]+all[2][1]+all[2][2]))+"\n"
                    +1+"\t"+.5+"\t"+all[2][1]+"\t"+df.format((all[2][1])/(all[2][0]+all[2][1]+all[2][2]))+"\n"
                    +1+"\t"+1+"\t"+all[2][2]+"\t"+df.format((all[2][2])/(all[2][0]+all[2][1]+all[2][2]))+"\n");
            outStream.writeBytes("#Proportion unimputed:\n#minor <- "+all[0][3]/all[0][4]+"\n#het<- "+all[1][3]/all[1][4]+"\n#major<- "+all[2][3]/all[2][4]+"\n");
            outStream.writeBytes("#Time to impute and calculate accuracy: "+time+" seconds");
            System.out.println("##Taxon\tTotalSitesMasked\tTotalSitesCompared\tTotalPropUnimputed\tNumMinor\tCorrectMinor\tMinorToHet\tMinorToMajor\tUnimpMinor"
                    + "\tNumHets\tHetToMinor\tCorrectHet\tHetToMajor\tUnimpHet\tNumMajor\tMajorToMinor\tMajorToHet\tCorrectMajor\tUnimpMajor\tR2");
            System.out.println("TotalByImputed\t"+(all[0][4]+all[1][4]+all[2][4])+"\t"+(all[0][4]+all[1][4]+all[2][4]-all[0][3]-all[1][3]-all[2][3])+"\t"+
                    ((all[0][3]+all[1][3]+all[2][3])/(all[0][4]+all[1][4]+all[2][4]))+"\t"+all[0][4]+"\t"+all[0][0]+"\t"+all[0][1]+"\t"+all[0][2]+"\t"+all[0][3]+
                    "\t"+all[1][4]+"\t"+all[1][0]+"\t"+all[1][1]+"\t"+all[1][2]+"\t"+all[1][3]+"\t"+all[2][4]+"\t"+all[2][0]+"\t"+all[2][1]+"\t"+all[2][2]+
                    "\t"+all[2][3]+"\t"+r2);
            System.out.println("Proportion unimputed:\nminor: "+all[0][3]/all[0][4]+"\nhet: "+all[1][3]/all[1][4]+"\nmajor: "+all[2][3]/all[2][4]);
            System.out.println("#Minor=0,Het=1,Major=2;x is masked(known), y is predicted\nx\ty\tN\tprop\n"
                    +0+"\t"+0+"\t"+all[0][0]+"\t"+(all[0][0])/(all[0][0]+all[0][1]+all[0][2])+"\n"
                    +0+"\t"+.5+"\t"+all[0][1]+"\t"+(all[0][1])/(all[0][0]+all[0][1]+all[0][2])+"\n"
                    +0+"\t"+1+"\t"+all[0][2]+"\t"+(all[0][2])/(all[0][0]+all[0][1]+all[0][2])+"\n"
                    +.5+"\t"+0+"\t"+all[1][0]+"\t"+(all[1][0])/(all[1][0]+all[1][1]+all[1][2])+"\n"
                    +.5+"\t"+.5+"\t"+all[1][1]+"\t"+(all[1][1])/(all[1][0]+all[1][1]+all[1][2])+"\n"
                    +.5+"\t"+1+"\t"+all[1][2]+"\t"+(all[1][2])/(all[1][0]+all[1][1]+all[1][2])+"\n"
                    +1+"\t"+0+"\t"+all[2][0]+"\t"+(all[2][0])/(all[2][0]+all[2][1]+all[2][2])+"\n"
                    +1+"\t"+.5+"\t"+all[2][1]+"\t"+(all[2][1])/(all[2][0]+all[2][1]+all[2][2])+"\n"
                    +1+"\t"+1+"\t"+all[2][2]+"\t"+(all[2][2])/(all[2][0]+all[2][1]+all[2][2])+"\n");
            outStream.close();
        } catch (Exception e) {
            System.out.println(e);
        }
    }
    
    private void accuracyMAFOut(double[][][] mafAll) {
        DecimalFormat df = new DecimalFormat("0.########");
        if (MAF!=null && MAFClass!=null) try {
            File outputFile = new File(outFileBase.substring(0, outFileBase.indexOf(".hmp")) + "DepthAccuracyMAF.txt");
            DataOutputStream outStream = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outputFile)));
            outStream.writeBytes("##\tMAFClass\tTotalSitesMasked\tTotalSitesCompared\tTotalPropUnimputed\tNumHets\tHetToMinor\tHetToMajor\tCorrectHet\tUnimpHet\tNumMinor\tMinorToMajor\tMinorToHet\tCorrectMinor\t"
                    + "UnimpMinor\tNumMajor\tMajorToMinor\tMajorToHet\tCorrectMajor\tUnimputedMajor\tr2\n");
            for (int i= 0; i<MAFClass.length;i++) {
                outStream.writeBytes("##TotalByImputed\t"+MAFClass[i]+"\t"+(mafAll[i][0][4]+mafAll[i][1][4]+mafAll[i][2][4])+"\t"+(mafAll[i][0][4]+mafAll[i][1][4]+mafAll[i][2][4]-mafAll[i][0][3]-mafAll[i][1][3]-mafAll[i][2][3])+"\t"+
                    ((mafAll[i][0][3]+mafAll[i][1][3]+mafAll[i][2][3])/(mafAll[i][0][4]+mafAll[i][1][4]+mafAll[i][2][4]))+"\t"+mafAll[i][0][4]+"\t"+mafAll[i][0][0]+"\t"+mafAll[i][0][1]+"\t"+mafAll[i][0][2]+"\t"+mafAll[i][0][3]+
                    "\t"+mafAll[i][1][4]+"\t"+mafAll[i][1][0]+"\t"+mafAll[i][1][1]+"\t"+mafAll[i][1][2]+"\t"+mafAll[i][1][3]+"\t"+mafAll[i][2][4]+"\t"+mafAll[i][2][0]+"\t"+mafAll[i][2][1]+"\t"+mafAll[i][2][2]+
                    "\t"+mafAll[i][2][3]+"\t"+pearsonR2(mafAll[i])+"\n");
            }
            outStream.writeBytes("#MAFClass,Minor=0,Het=1,Major=2;x is masked(known), y is predicted\nMAF\tx\ty\tN\tprop\n");
            for (int i= 0; i<MAFClass.length;i++) { outStream.writeBytes(
                    MAFClass[i]+"\t"+0+"\t"+0+"\t"+mafAll[i][0][0]+"\t"+df.format((mafAll[i][0][0])/(mafAll[i][0][0]+mafAll[i][0][1]+mafAll[i][0][2]))+"\n"
                    +MAFClass[i]+"\t"+0+"\t"+.5+"\t"+mafAll[i][0][1]+"\t"+df.format((mafAll[i][0][1])/(mafAll[i][0][0]+mafAll[i][0][1]+mafAll[i][0][2]))+"\n"
                    +MAFClass[i]+"\t"+0+"\t"+1+"\t"+mafAll[i][0][2]+"\t"+df.format((mafAll[i][0][2])/(mafAll[i][0][0]+mafAll[i][0][1]+mafAll[i][0][2]))+"\n"
                    +MAFClass[i]+"\t"+.5+"\t"+0+"\t"+mafAll[i][1][0]+"\t"+df.format((mafAll[i][1][0])/(mafAll[i][1][0]+mafAll[i][1][1]+mafAll[i][1][2]))+"\n"
                    +MAFClass[i]+"\t"+.5+"\t"+.5+"\t"+mafAll[i][1][1]+"\t"+df.format((mafAll[i][1][1])/(mafAll[i][1][0]+mafAll[i][1][1]+mafAll[i][1][2]))+"\n"
                    +MAFClass[i]+"\t"+.5+"\t"+1+"\t"+mafAll[i][1][2]+"\t"+df.format((mafAll[i][1][2])/(mafAll[i][1][0]+mafAll[i][1][1]+mafAll[i][1][2]))+"\n"
                    +MAFClass[i]+"\t"+1+"\t"+0+"\t"+mafAll[i][2][0]+"\t"+df.format((mafAll[i][2][0])/(mafAll[i][2][0]+mafAll[i][2][1]+mafAll[i][2][2]))+"\n"
                    +MAFClass[i]+"\t"+1+"\t"+.5+"\t"+mafAll[i][2][1]+"\t"+df.format((mafAll[i][2][1])/(mafAll[i][2][0]+mafAll[i][2][1]+mafAll[i][2][2]))+"\n"
                    +MAFClass[i]+"\t"+1+"\t"+1+"\t"+mafAll[i][2][2]+"\t"+df.format((mafAll[i][2][2])/(mafAll[i][2][0]+mafAll[i][2][1]+mafAll[i][2][2]))+"\n");
            }
            outStream.writeBytes("#Proportion unimputed:\n#MAF\tminor\thet\tmajor\n");
            for (int i= 0; i<MAFClass.length;i++) { 
                outStream.writeBytes("#"+MAFClass[i]+"\t"+mafAll[i][0][3]/mafAll[i][0][4]+"\t"+mafAll[i][1][3]/mafAll[i][1][4]+"\t"+mafAll[i][2][3]/mafAll[i][2][4]+"\n");
            }
            outStream.flush();
            outStream.close();
        } catch (Exception e) {
            System.out.println(e);
        }
    }

    private int[] countUnknownAndHets(byte[] a) {
        int cnt=0, cntHets=0;
        for (int i = 0; i < a.length; i++) {
            if(a[i]==UNKNOWN_DIPLOID_ALLELE) {cnt++;}
            else if(isHeterozygous(a[i])) {cntHets++;}
        }
        return new int[]{cnt,cntHets};
    }

    private BitSet[] arrangeMajorMinorBtwAlignments(GenotypeTable unimpAlign, int bt, int donorOffset, int donorLength,
            OpenBitSet goodMask, OpenBitSet swapMjMnMask) {
        int unimpAlignStartBlock=donorOffset/64;
        int shift=(donorOffset-(unimpAlignStartBlock*64));
        int unimpAlignEndBlock=unimpAlignStartBlock+((donorLength+shift-1)/64);
        OpenBitSet mjUnImp=new OpenBitSet(unimpAlign.allelePresenceForSitesBlock(bt, Major, unimpAlignStartBlock, unimpAlignEndBlock + 1));
        OpenBitSet mnUnImp=new OpenBitSet(unimpAlign.allelePresenceForSitesBlock(bt, Minor, unimpAlignStartBlock, unimpAlignEndBlock + 1));
        OpenBitSet mjTbs=new OpenBitSet(donorLength);
        OpenBitSet mnTbs=new OpenBitSet(donorLength);
        for (int i = 0; i < donorLength; i++) {
            if(mjUnImp.fastGet(i+shift)) mjTbs.set(i);
            if(mnUnImp.fastGet(i+shift)) mnTbs.set(i);
        }
        OpenBitSet newmj=new OpenBitSet(mjTbs);
        OpenBitSet newmn=new OpenBitSet(mnTbs);
        mjTbs.and(goodMask);
        mnTbs.and(goodMask);
        //       System.out.printf("mjTbs:%d Goodmask:%d swapMjMnMask:%d",mjTbs.getNumWords(),goodMask.getNumWords(), swapMjMnMask.getNumWords());
        if(isSwapMajorMinor) {
            newmj.and(swapMjMnMask);
            newmn.and(swapMjMnMask);
            mjTbs.or(newmn);
            mnTbs.or(newmj);
        }
//        System.out.printf("Arrange shift:%d finalEnd:%d totalBits:%d DesiredLength:%d %n", shift, finalEnd, totalBits, donorLength);
        BitSet[] result={mjTbs,mnTbs};
        return result;
    }

    /**
     * If the target regions has Mendelian errors that it looks for overlapping regional
     * solutions that are better.
     * TODO:  This can be done much more robustly.
     * @param impT
     * @param targetTaxon
     * @param regionHypth
     */
    private ImputedTaxon applyBlockNN(ImputedTaxon impT, int targetTaxon,
            GenotypeTable donorAlign, int donorOffset, DonorHypoth[][] regionHypth, boolean hybridMode, BitSet[] maskedTargetBits, 
            int minMinorCnt, double focusInbredErr, double focusHybridErr, double focusSmashErr, int[] donorIndices, int[] blockNN, boolean hetsToMiss) {
        int[] currBlocksSolved= new int[5];//track number of focus blocks solved in NN search for system out; index 0 is inbred, 1 is viterbi, 2 is smash, 3 is not solved, 4 is total for all modes
        int blocks=maskedTargetBits[0].getNumWords();
        for (int focusBlock = 0; focusBlock < blocks; focusBlock++) {
            boolean vit= false;//just if viterbi works for this block
            int[] d= getUniqueDonorsForBlock(regionHypth,focusBlock);//get the donor indices for that block from regionHypoth
            int[] resultRange=getBlockWithMinMinorCount(maskedTargetBits[0].getBits(),maskedTargetBits[1].getBits(), focusBlock, minMinorCnt);
            if(resultRange==null) {currBlocksSolved[3]++; continue; } //no data in the focus Block
            if(d==null) {currBlocksSolved[3]++; continue; } //no donors for the focus block
            //search for the best hybrid donors for a segment
            DonorHypoth[] best2donors=getBestHybridDonors(targetTaxon, maskedTargetBits[0].getBits(resultRange[0], resultRange[2]),
                        maskedTargetBits[1].getBits(resultRange[0], resultRange[2]), resultRange[0], resultRange[2], focusBlock, donorAlign, d, d, true);
            if(best2donors[0]==null) {currBlocksSolved[3]++; continue; } //no good hybrid donors for the focus block

            //check for Viterbi and inbred
            ArrayList<DonorHypoth> goodDH= new ArrayList<DonorHypoth>(); //KLS0201
            if (best2donors[0].getErrorRate()<focusInbredErr) {
                for (DonorHypoth dh : best2donors) {
                    if((dh!=null)) {
                        if((dh.isInbred()==false)&&(dh.getErrorRate()<focusHybridErr)){
                            dh=getStateBasedOnViterbi(dh, donorOffset, donorAlign, twoWayViterbi, transition);
                            if(dh!=null) vit= true;
                        }
                    if(dh!=null&&(dh.getErrorRate()<focusInbredErr)) goodDH.add(dh);
                    }
                }
                if (goodDH.size()!=0) {
                    DonorHypoth[] vdh=new DonorHypoth[goodDH.size()];
                    for (int i = 0; i < vdh.length; i++) {vdh[i]=goodDH.get(i);}
                    regionHypth[focusBlock]= vdh;
                    impT= setAlignmentWithDonors(donorAlign, regionHypth[focusBlock], donorOffset, true, impT, false, false);
                        impT.incBlocksSolved();
                    if(vit==true) {currBlocksSolved[1]++;} else {currBlocksSolved[0]++;}
                    currBlocksSolved[4]++; continue;
                    }
                }
            //if fails, try smash mode
            else if (hybridMode&&(best2donors[0].getErrorRate()<focusSmashErr)) {//smash mode goes into hybrid block NN
                for (DonorHypoth dh:best2donors) {
                    if(dh!=null&&dh.getErrorRate()<focusSmashErr) {
                        goodDH.add(dh);
                }
                    if (goodDH.size()!=0) {
                        DonorHypoth[] vdh=new DonorHypoth[goodDH.size()];
                        for (int i = 0; i < vdh.length; i++) {vdh[i]=goodDH.get(i);}
                        regionHypth[focusBlock]= vdh;
                        impT= setAlignmentWithDonors(donorAlign, regionHypth[focusBlock], donorOffset, true,impT, true, hetsToMiss);//only set donors for focus block //KLS0201
                        impT.incBlocksSolved(); currBlocksSolved[2]++; currBlocksSolved[4]++; continue;
            }
        }
            //if fails, do not impute this focus block    
            } else currBlocksSolved[3]++;
        }
        
        if (currBlocksSolved[4]!=0) impT.setSegmentSolved(true);
        else impT.setSegmentSolved(false);
        int leftNullCnt=currBlocksSolved[3];
        if(testing==1)
            System.out.printf("targetTaxon:%d hybridError:%g block:%d proportionBlocksImputed:%d null:%d inbredDone:%d viterbiDone:%d hybridDone:%d noData:%d %n",
                targetTaxon, focusHybridErr, blocks,currBlocksSolved[4]/blocks, leftNullCnt, currBlocksSolved[0], currBlocksSolved[1], currBlocksSolved[2], currBlocksSolved[3]);
        for (int i = 0; i < currBlocksSolved.length; i++) {blockNN[i]+= currBlocksSolved[i];}
        return impT;
    }

    private DonorHypoth getStateBasedOnViterbi(DonorHypoth dh, int donorOffset, GenotypeTable donorAlign, boolean forwardReverse, double[][] trans) {
        TransitionProbability tpF = new TransitionProbability();
        EmissionProbability ep = new EmissionProbability();
        tpF.setTransitionProbability(trans);
        ep.setEmissionProbability(emission);
        int startSite=dh.startBlock*64;
        int endSite=(dh.endBlock*64)+63;
        if(endSite>=donorAlign.numberOfSites()) endSite=donorAlign.numberOfSites()-1;
        int sites=endSite-startSite+1;
        byte[] callsF=new byte[sites];
        byte[] callsR=new byte[sites];
        //System.out.printf("%d %d %d %n",dh.donor1Taxon, startSite, endSite+1);
        byte[] d1b=donorAlign.genotypeRange(dh.donor1Taxon, startSite, endSite + 1);
        byte[] d2b=donorAlign.genotypeRange(dh.donor2Taxon, startSite, endSite + 1);
        byte[] t1b=unimpAlign.genotypeRange(dh.targetTaxon, startSite + donorOffset, endSite + 1 + donorOffset);
        int informSites=0, nonMendel=0;
        ArrayList<Byte> nonMissingObs = new ArrayList<Byte>();
        ArrayList<Integer> snpPositions = new ArrayList<Integer>();
        for(int cs=0; cs<sites; cs++) {
            if(t1b[cs]==UNKNOWN_DIPLOID_ALLELE) continue;
            if(d1b[cs]==UNKNOWN_DIPLOID_ALLELE) continue;
            if(d2b[cs]==UNKNOWN_DIPLOID_ALLELE) continue;
            if(t1b[cs]==GAP_DIPLOID_ALLELE) continue;
            if(d1b[cs]==GAP_DIPLOID_ALLELE) continue;
            if(d2b[cs]==GAP_DIPLOID_ALLELE) continue;
            if(d1b[cs]==d2b[cs]) {
                if(t1b[cs]!=d1b[cs]) nonMendel++;
                continue;
            }
            informSites++;
            byte state=1;
            if(t1b[cs]==d1b[cs]) {state=0;}
            else if(t1b[cs]==d2b[cs]) {state=2;}
            nonMissingObs.add(state);
            snpPositions.add(cs+startSite);
        }
        if(informSites<10) return null;
        double nonMendRate=(double)nonMendel/(double)informSites;
        if(testing==1) System.out.printf("NonMendel:%d InformSites:%d ErrorRate:%g %n",nonMendel, informSites, nonMendRate);
        if(nonMendRate>this.maximumInbredError*5) return null;
        byte[] informStatesF=new byte[informSites];
        for (int i = 0; i < informStatesF.length; i++) informStatesF[i]=nonMissingObs.get(i);
        int[] pos=new int[informSites];
        for (int i = 0; i < pos.length; i++) pos[i]=snpPositions.get(i);
        int chrlength = donorAlign.chromosomalPosition(endSite) - donorAlign.chromosomalPosition(startSite);
        tpF.setAverageSegmentLength( chrlength / sites );
        tpF.setPositions(pos);

        double probHeterozygous=0.5;
        double phom = (1 - probHeterozygous) / 2;
        double[] pTrue = new double[]{phom, .25*probHeterozygous ,.5 * probHeterozygous, .25*probHeterozygous, phom};
        ViterbiAlgorithm vaF = new ViterbiAlgorithm(informStatesF, tpF, ep, pTrue);
        vaF.calculate();
        if(testing==1) System.out.println("Input:"+Arrays.toString(informStatesF));
        byte[] resultStatesF=vaF.getMostProbableStateSequence();
        if(testing==1) System.out.println("Resul:"+Arrays.toString(resultStatesF));
        DonorHypoth dh2=new DonorHypoth(dh.targetTaxon,dh.donor1Taxon,
                dh.donor2Taxon, dh.startBlock, dh.focusBlock, dh.endBlock);
        int currPos=0;
        for(int cs=0; cs<sites; cs++) {
            callsF[cs]=(resultStatesF[currPos]==1)?(byte)1:(byte)(resultStatesF[currPos]/2); //converts the scale back to 0,1,2 from 0..4
            if((pos[currPos]<cs+startSite)&&(currPos<resultStatesF.length-1)) currPos++;
        }

        if (forwardReverse==true) {
            TransitionProbability tpR = new TransitionProbability();
            tpR.setTransitionProbability(transition);
            byte[] informStatesR= informStatesF;
            ArrayUtils.reverse(informStatesR);
            int[] posR= pos;
            ArrayUtils.reverse(posR);
            tpR.setAverageSegmentLength( chrlength / sites );
            tpR.setPositions(posR);
            ViterbiAlgorithm vaR = new ViterbiAlgorithm(informStatesR, tpR, ep, pTrue);
            vaR.calculate();
            byte[] resultStatesR=vaR.getMostProbableStateSequence();//this sequence is backwards/from the reverse viterbi
            ArrayUtils.reverse(resultStatesR); //flip the reverse viterbi calls to the same orientation as forward
            int currPosR=0;
            for(int cs=0; cs<sites; cs++) { //converts the scale back to 0,1,2 from 0..4
                callsR[cs]=(resultStatesR[currPosR]==1)?(byte)1:(byte)(resultStatesR[currPosR]/2);
                if((pos[currPosR]<cs+startSite)&&(currPosR<resultStatesF.length-1)) currPosR++;
            }
            //compare the forward and reverse viterbi, use the one with the longest path length if they contradict
            byte[] callsC=callsF;
            for(int i= 0;i<pos.length;i++) {
                int cs= pos[i]-startSite;
                if (callsF[cs]!=callsR[cs]&&i<pos.length/2) callsC[cs]= callsR[cs];
            }
            if (testing==1) {
                if (resultStatesF[0]!=resultStatesR[0]||resultStatesF[resultStatesF.length-1]!=resultStatesR[resultStatesR.length-1]) System.out.println("FR:\n"+Arrays.toString(informStatesF)+"\n"+Arrays.toString(informStatesR)+"\n"+Arrays.toString(resultStatesF)+"\n"+Arrays.toString(resultStatesR)+"\n"+Arrays.toString(callsF)+"\n"+Arrays.toString(callsR)+"\n"+Arrays.toString(callsC));
            }
            dh2.phasedResults= callsC;
            return dh2;
        }

        dh2.phasedResults= callsF;
        return dh2;
    }



    /**
     * Given a start 64 site block, it expands to the left and right until it hits
     * the minimum Minor Site count in the target taxon
     * @param mnT - minor allele bit presence in a series of longs
     * @param focusBlock
     * @param minMinorCnt
     * @return arrays of blocks {startBlock, focusBlock, endBlock}
     */
    private int[] getBlockWithMinMinorCount(long[] mjT, long[] mnT, int focusBlock, int minMinorCnt) {
        int blocks=mjT.length;
        int majorCnt=Long.bitCount(mjT[focusBlock]);
        int minorCnt=Long.bitCount(mnT[focusBlock]);
        int endBlock=focusBlock, startBlock=focusBlock;
        int minMajorCnt=minMinorCnt*minMajorRatioToMinorCnt;
        while((minorCnt<minMinorCnt)&&(majorCnt<minMajorCnt)) {
            boolean preferMoveStart=(focusBlock-startBlock<endBlock-focusBlock)?true:false;
            if(startBlock==0) {preferMoveStart=false;}
            if(endBlock==blocks-1) {preferMoveStart=true;}
            if((startBlock==0)&&(endBlock==blocks-1)) break;
            if(preferMoveStart) {//expand start
                startBlock--;
                minorCnt+=Long.bitCount(mnT[startBlock]);
                majorCnt+=Long.bitCount(mjT[startBlock]);
            } else { //expand end
                endBlock++;
                minorCnt+=Long.bitCount(mnT[endBlock]);
                majorCnt+=Long.bitCount(mjT[startBlock]);
            }
        }
        int[] result={startBlock, focusBlock, endBlock};
        return result;
    }


    /**
     * Calculates an inbred distance
     * TODO:  Should be converted to the regular same, diff, hets, & sum
     * @param modBitsOfTarget
     * @return array with sites, same count, diff count, het count
     */
    private void calcInbredDist(ImputedTaxon impT, BitSet[] modBitsOfTarget, GenotypeTable donorAlign) {
        int blocks=modBitsOfTarget[0].getNumWords();
        impT.allDist=new byte[donorAlign.numberOfTaxa()][4][blocks];
        long[] iMj=modBitsOfTarget[0].getBits();
        long[] iMn=modBitsOfTarget[1].getBits();
        for (int donor1 = 0; donor1 < impT.allDist.length; donor1++) {
            long[] jMj=donorAlign.allelePresenceForAllSites(donor1, Major).getBits();
            long[] jMn=donorAlign.allelePresenceForAllSites(donor1, Minor).getBits();
            for (int i = 0; i <blocks; i++) {
                long same = (iMj[i] & jMj[i]) | (iMn[i] & jMn[i]);
                long diff = (iMj[i] & jMn[i]) | (iMn[i] & jMj[i]);
                long hets = same & diff;
                int sameCnt = BitUtil.pop(same);
                int diffCnt = BitUtil.pop(diff);
                int hetCnt = BitUtil.pop(hets);
                int sites = sameCnt + diffCnt - hetCnt;
                impT.allDist[donor1][2][i]=(byte)diffCnt;
                impT.allDist[donor1][3][i]=(byte)hetCnt;
                impT.allDist[donor1][0][i]=(byte)sites;
                impT.allDist[donor1][1][i]=(byte)sameCnt;
            }
        }
    }

    /**
     * Simple algorithm that tests every possible two donor combination to minimize
     * the number of unmatched informative alleles.  Currently, there is litte tie
     * breaking, longer matches are favored.
     * @param targetTaxon
     * @param startBlock
     * @param endBlock
     * @param focusBlock
     * @return int[] array of {donor1, donor2, testSites}
     */
    private DonorHypoth[] getBestInbredDonors(int targetTaxon, ImputedTaxon impT, int startBlock, int endBlock,
            int focusBlock, GenotypeTable donorAlign, int[] donor1indices) {
        TreeMap<Double,DonorHypoth> bestDonors=new TreeMap<Double,DonorHypoth>();
        double lastKeytestPropUnmatched=1.0;
        double inc=1e-9;
        int donorTaxaCnt=donorAlign.numberOfTaxa();
        for (int d1 : donor1indices) {
            int testSites=0;
            int sameCnt = 0, diffCnt = 0, hetCnt = 0;
            for (int i = startBlock; i <=endBlock; i++) {
                sameCnt+=impT.allDist[d1][1][i];
                diffCnt+=impT.allDist[d1][2][i];
                hetCnt+=impT.allDist[d1][3][i];
            }
            testSites= sameCnt + diffCnt - hetCnt;
            if(testSites<minTestSites) continue;
            double testPropUnmatched = 1.0-(((double) (sameCnt) - (double)(0.5*hetCnt)) / (double) (testSites));
            int totalMendelianErrors=(int)((double)testSites*testPropUnmatched);
            //prob of this given poisson with lambda=1, total Mendelian errors=k (WE SHOULD DO THIS)
            inc+=1e-9;  //this looks strange, but makes the keys unique and ordered
            testPropUnmatched+=inc;  //this looks strange, but makes the keys unique and ordered
            if(testPropUnmatched<lastKeytestPropUnmatched) {
                DonorHypoth theDH=new DonorHypoth(targetTaxon, d1, d1, startBlock,
                        focusBlock, endBlock, testSites, totalMendelianErrors);
                DonorHypoth prev=bestDonors.put(new Double(testPropUnmatched), theDH);
                if(prev!=null) System.out.println("Inbred TreeMap index crash:"+testPropUnmatched);
                if(bestDonors.size()>maxDonorHypotheses) {
                    bestDonors.remove(bestDonors.lastKey());
                    lastKeytestPropUnmatched=bestDonors.lastKey();
                }
            }
        }
        //TODO consider reranking by Poisson.  The calculating Poisson on the fly is too slow.
        DonorHypoth[] result=new DonorHypoth[maxDonorHypotheses];
        int count=0;
        for (DonorHypoth dh : bestDonors.values()) {
            result[count]=dh;
            count++;
        }
        return result;
    }

    /**
     * Produces a sort list of most prevalent donorHypotheses across the donorAlign.
     * Currently it is only looking at the very best for each focus block.
     * @param allDH
     * @param minHypotheses
     * @return
     */
    private int[] getAllBestDonorsAcrossChromosome(DonorHypoth[][] allDH, int minHypotheses) {
        TreeMap<Integer,Integer> bd=new TreeMap<Integer,Integer>();
        for (int i = 0; i < allDH.length; i++) {
            DonorHypoth dh=allDH[i][0];
            if(dh==null) continue;
            if(bd.containsKey(dh.donor1Taxon)) {
                bd.put(dh.donor1Taxon, bd.get(dh.donor1Taxon)+1);
            } else {
                bd.put(dh.donor1Taxon, 1);
            }
        }
        if(testing==1) System.out.println(bd.size()+":"+bd.toString());
        if(testing==1) System.out.println("");
        int highDonors=0;
        for (int i : bd.values()) {if(i>minHypotheses) highDonors++;}
        int[] result=new int[highDonors];
        highDonors=0;
        for (Map.Entry<Integer,Integer> e: bd.entrySet()) {
            if(e.getValue()>minHypotheses) result[highDonors++]=e.getKey();}
        if(testing==1) System.out.println(Arrays.toString(result));
        return result;
    }

    private int[] getUniqueDonorsForBlock(DonorHypoth[][] regionHypth, int block) {//change so only adding those donors that are not identical for target blocks
        TreeSet<Integer> donors= new TreeSet<>();
        for (int h = 0; h < regionHypth[block].length; h++) {
            if (regionHypth[block][h]!=null) {
                donors.add(regionHypth[block][h].donor1Taxon);
                donors.add(regionHypth[block][h].donor2Taxon);
            }
        }
        if (donors.isEmpty()) return null;
        int[] result= ArrayUtils.toPrimitive(donors.toArray(new Integer[donors.size()]));
        return result;
    }

    /**
     * Simple algorithm that tests every possible two donor combination to minimize
     * the number of unmatched informative alleles.  Currently, there is litte tie
     * breaking, longer matches are favored.
     * @param targetTaxon
     * @param startBlock
     * @param endBlock
     * @param focusBlock
     * @return int[] array of {donor1, donor2, testSites}
     */
    private DonorHypoth[] getBestHybridDonors(int targetTaxon, long[] mjT, long[] mnT,
                                              int startBlock, int endBlock, int focusBlock, GenotypeTable donorAlign, int[] donor1Indices, int[] donor2Indices,
                                              boolean viterbiSearch) {
        TreeMap<Double,DonorHypoth> bestDonors=new TreeMap<Double,DonorHypoth>();
        double lastKeytestPropUnmatched=1.0;
        double inc=1e-9;
        double[] donorDist;
        for (int d1: donor1Indices) {
            long[] mj1=donorAlign.allelePresenceForSitesBlock(d1, Major, startBlock, endBlock + 1);
            long[] mn1=donorAlign.allelePresenceForSitesBlock(d1, Minor, startBlock, endBlock + 1);
            for (int d2 : donor2Indices) {
                if((!viterbiSearch)&&(d1==d2)) continue;
                long[] mj2=donorAlign.allelePresenceForSitesBlock(d2, Major, startBlock, endBlock + 1);
                long[] mn2=donorAlign.allelePresenceForSitesBlock(d2, Minor, startBlock, endBlock + 1);
                if(viterbiSearch) {
                    donorDist=IBSDistanceMatrix.computeHetBitDistances(mj1, mn1, mj2, mn2, minTestSites);
                    if((d1!=d2)&&(donorDist[0]<this.maximumInbredError)) continue;
                }
                int[] mendErr=mendelErrorComparison(mjT, mnT, mj1, mn1, mj2, mn2);
                if(mendErr[1]<minTestSites) continue;
                double testPropUnmatched=(double)(mendErr[0])/(double)mendErr[1];
                inc+=1e-9;
                testPropUnmatched+=inc;
                if(testPropUnmatched<lastKeytestPropUnmatched) {
                    DonorHypoth theDH=new DonorHypoth(targetTaxon, d1, d2, startBlock,
                            focusBlock, endBlock, mendErr[1], mendErr[0]);
                    DonorHypoth prev=bestDonors.put(new Double(testPropUnmatched), theDH);
                    if(prev!=null) {
                        System.out.println("Hybrid TreeMap index crash:"+testPropUnmatched);
                    }
                    if(bestDonors.size()>maxDonorHypotheses) {
                        bestDonors.remove(bestDonors.lastKey());
                        lastKeytestPropUnmatched=bestDonors.lastKey();
                    }
                }
            }
        }
        DonorHypoth[] result=new DonorHypoth[maxDonorHypotheses];
        int count=0;
        for (DonorHypoth dh : bestDonors.values()) {
            result[count]=dh;
            count++;
        }
        return result;
    }

    private int[] mendelErrorComparison(long[] mjT, long[] mnT, long[] mj1, long[] mn1,
                                        long[] mj2, long[] mn2) {
        int mjUnmatched=0;
        int mnUnmatched=0;
        int testSites=0;
        for (int i = 0; i < mjT.length; i++) {
            long siteMask=(mjT[i]|mnT[i])&(mj1[i]|mn1[i])&(mj2[i]|mn2[i]);
            mjUnmatched+=Long.bitCount(siteMask&mjT[i]&(mjT[i]^mj1[i])&(mjT[i]^mj2[i]));
            mnUnmatched+=Long.bitCount(siteMask&mnT[i]&(mnT[i]^mn1[i])&(mnT[i]^mn2[i]));
            testSites+=Long.bitCount(siteMask);
        }
        int totalMendelianErrors=mjUnmatched+mnUnmatched;
        // double testPropUnmatched=(double)(totalMendelianErrors)/(double)testSites;
        return (new int[] {totalMendelianErrors, testSites});
    }


    private ImputedTaxon setAlignmentWithDonors(GenotypeTable donorAlign, DonorHypoth[] theDH, int donorOffset,
                                                boolean setJustFocus, ImputedTaxon impT, boolean smashOn, boolean hetsMiss) {
        if(theDH[0].targetTaxon<0) return impT;
        boolean print=false;
        int startSite=(setJustFocus)?theDH[0].getFocusStartSite():theDH[0].startSite;
        int endSite=(setJustFocus)?theDH[0].getFocusEndSite():theDH[0].endSite;
        if(endSite>=donorAlign.numberOfSites()) endSite=donorAlign.numberOfSites()-1;
//        if (print) System.out.println("B:"+mna.getBaseAsStringRange(theDH[0].targetTaxon, startSite, endSite));
//        int[] prevDonors;
        int[] prevDonors=new int[]{-1, -1};
        if(theDH[0].getPhaseForSite(startSite)==0) {prevDonors=new int[]{theDH[0].donor1Taxon, theDH[0].donor1Taxon};}
        else if(theDH[0].getPhaseForSite(startSite)==2) {prevDonors=new int[]{theDH[0].donor2Taxon, theDH[0].donor2Taxon};}
        else if(theDH[0].getPhaseForSite(startSite)==1) {prevDonors=new int[]{theDH[0].donor1Taxon, theDH[0].donor2Taxon};}
//        if(impT.areBreakPointsEmpty()) {
//            prevDonors=new int[]{-1, -1};
//        } else {
//            prevDonors=impT.breakPoints.lastEntry().getValue();
//        }
        int[] currDonors=prevDonors;

        //this generates an alignment of the closest donors for the taxon block from the results of the hybrid donor search (to use for frequencies using 2PQ)
//        TreeSet<Identifier> closeDonors;
//        closeDonors= new TreeSet<>();
//        for (DonorHypoth dh:theDH) {
//            closeDonors.add(donorAlign.getIdGroup().getIdentifier(dh.donor1Taxon)); closeDonors.add(donorAlign.getIdGroup().getIdentifier(dh.donor2Taxon));
//        }
//        IdGroup ids= new SimpleIdGroup(closeDonors.toArray(new Identifier[closeDonors.size()]));
//            Alignment bestDonors= FilterAlignment.getInstance(donorAlign, ids);
////
        for(int cs=startSite; cs<=endSite; cs++) {
            byte donorEst=UNKNOWN_DIPLOID_ALLELE;
            byte neighbor=0;
            for (int i = 0; (i < theDH.length) && (donorEst==UNKNOWN_DIPLOID_ALLELE); i++) {
                neighbor++;
                if((theDH[i]==null)||(theDH[i].donor1Taxon<0)) continue;
                if(theDH[i].getErrorRate()>this.maximumInbredError) continue;
                byte bD1=donorAlign.genotype(theDH[i].donor1Taxon, cs);
                if(theDH[i].getPhaseForSite(cs)==0) {
                    donorEst=bD1;
                    if(i==0) currDonors=new int[]{theDH[0].donor1Taxon, theDH[0].donor1Taxon};
                }
                else {
                    byte bD2=donorAlign.genotype(theDH[i].donor2Taxon, cs);
                    if(theDH[i].getPhaseForSite(cs)==2) {
                        donorEst=bD2;
                        if(i==0) currDonors=new int[]{theDH[0].donor2Taxon, theDH[0].donor2Taxon};
                    }
                    else {
                        donorEst=GenotypeTableUtils.getUnphasedDiploidValueNoHets(bD1, bD2);
                        if(i==0) currDonors=new int[]{theDH[0].donor1Taxon, theDH[0].donor2Taxon};
                    }
                }
            }
            byte knownBase=impT.getOrigGeno(cs+donorOffset);
            String knownBaseString= NucleotideAlignmentConstants.getNucleotideIUPAC(knownBase);
            if(!Arrays.equals(prevDonors, currDonors)) {
      //          impT.breakPoints.put(donorAlign.chromosomalPosition(cs), currDonors);    //TODO fix this for projection
                prevDonors=currDonors;
            }
            if(theDH[0].phasedResults==null) {impT.chgHis[cs+donorOffset]=(byte)-neighbor;}
            else {impT.chgHis[cs+donorOffset]=(byte)neighbor;}

            impT.impGeno[cs+donorOffset]= donorEst;  //predicted based on neighbor
            String donorEstString= NucleotideAlignmentConstants.getNucleotideIUPAC(donorEst);
            if(knownBase==UNKNOWN_DIPLOID_ALLELE) {
                if (isHeterozygous(donorEst)) {
                    if (smashOn && hetsMiss) {//if imputing a heterozygote, just set to missing
                        impT.resolveGeno[cs+donorOffset]= knownBase;
                        //System.out.println("SetToMissing"+":"+theDH[0].targetTaxon+":"+cs+donorOffset+":"+knownBaseString+":"+donorEstString);
                    }
                    else impT.resolveGeno[cs+donorOffset]= donorEst; //if not in hybrid, set to het
                }
                else {//if not imputed to a het
                    impT.resolveGeno[cs+donorOffset]= donorEst;}
            } else if (isHeterozygous(donorEst)){
                if(resolveHetIfUndercalled&&GenotypeTableUtils.isPartiallyEqual(knownBase, donorEst)&&smashOn==false){//if smash off, set homozygotes imputed to het to het
                    System.out.println("ResolveHet:"+theDH[0].targetTaxon+":"+cs+donorOffset+":"+knownBaseString+":"+donorEstString);
                    impT.resolveGeno[cs+donorOffset]= donorEst;
                }
            }
        } //end of cs loop
        //enter a stop of the DH at the beginning of the next block
        int lastDApos=donorAlign.chromosomalPosition(endSite);
        int nextSite=unimpAlign.siteOfPhysicalPosition(lastDApos, donorAlign.chromosome(0))+1;
       // if(nextSite<unimpAlign.numberOfSites()) impT.breakPoints.put(unimpAlign.chromosomalPosition(nextSite), new int[]{-1,-1});       //TODO fix this for projection
        //    if (print) System.out.println("E:"+mna.getBaseAsStringRange(theDH[0].targetTaxon, startSite, endSite));
        return impT;

        /**
         * Code that records results to projectionGenotype.  Some of this may be needed above.
         if(!Arrays.equals(prevDonors, currDonors)) {
         DonorHaplotypes dhaps=new DonorHaplotypes(donorAlign.chromosome(prevDonorStart), donorAlign.chromosomalPosition(prevDonorStart),
         donorAlign.chromosomalPosition(cs),prevDonors[0],prevDonors[1]);
         impT.addBreakPoint(dhaps);
         prevDonors=currDonors;
         prevDonorStart=cs;
         }
         if(theDH[0].phasedResults==null) {impT.chgHis[cs+donorOffset]=(byte)-neighbor;}
         else {impT.chgHis[cs+donorOffset]=(byte)neighbor;}

         impT.impGeno[cs+donorOffset]= donorEst;  //predicted based on neighbor
         if(knownBase==UNKNOWN_DIPLOID_ALLELE) {impT.resolveGeno[cs+donorOffset]= donorEst;}
         else {if(isHeterozygous(donorEst)) {
         if(resolveHetIfUndercalled&&GenotypeTableUtils.isPartiallyEqual(knownBase,donorEst))
         {//System.out.println(theDH[0].targetTaxon+":"+knownBase+":"+donorEst);
         impT.resolveGeno[cs+donorOffset]= donorEst;}
         }}
         } //end of cs loop
         DonorHaplotypes dhaps=new DonorHaplotypes(donorAlign.chromosome(prevDonorStart), donorAlign.chromosomalPosition(prevDonorStart),
         donorAlign.chromosomalPosition(endSite),prevDonors[0],prevDonors[1]);
         impT.addBreakPoint(dhaps);
         */
    }

    private double calcErrorForTaxonAndSite(ImputedTaxon impT) {
        for (int cs = 0; cs < impT.getOrigGeno().length; cs++) {
            byte donorEst=impT.getImpGeno(cs);
            byte knownBase=impT.getOrigGeno(cs);
            if((knownBase!=UNKNOWN_DIPLOID_ALLELE)&&(donorEst!=UNKNOWN_DIPLOID_ALLELE)) {
                if(isHeterozygous(donorEst)||isHeterozygous(knownBase)) {
                    totalHets++;
                } else if(knownBase==donorEst) {
                    totalRight++;
                    siteCorrectCnt[cs]++;
                    taxonCorrectCnt[impT.taxon()]++;
                } else {
                    totalWrong++;
                    siteErrors[cs]++;
                    taxonErrors[impT.taxon()]++;
                }
            }

        }
        return (double)taxonErrors[impT.taxon()]/(double)(taxonCorrectCnt[impT.taxon()]+taxonErrors[impT.taxon()]);
    }

    public static int[] compareAlignment(String origFile, String maskFile, String impFile, boolean noMask) {
        boolean taxaOut=false;
        GenotypeTable oA=ImportUtils.readGuessFormat(origFile);
        System.out.printf("Orig taxa:%d sites:%d %n",oA.numberOfTaxa(),oA.numberOfSites());
        GenotypeTable mA=null;
        if(noMask==false) {mA=ImportUtils.readGuessFormat(maskFile);
            System.out.printf("Mask taxa:%d sites:%d %n",mA.numberOfTaxa(),mA.numberOfSites());
        }
        GenotypeTable iA=ImportUtils.readGuessFormat(impFile);
        System.out.printf("Imp taxa:%d sites:%d %n",iA.numberOfTaxa(),iA.numberOfSites());
        int correct=0;
        int errors=0;
        int unimp=0;
        int hets=0;
        int gaps=0;
        for (int t = 0; t < iA.numberOfTaxa(); t++) {
            int e=0,c=0,u=0,h=0;
            int oATaxa=oA.taxa().indexOf(iA.taxaName(t));
            for (int s = 0; s < iA.numberOfSites(); s++) {
                if(noMask||(oA.genotype(oATaxa, s)!=mA.genotype(t, s))) {
                    byte ib=iA.genotype(t, s);
                    byte ob=oA.genotype(oATaxa, s);
                    if((ib==UNKNOWN_DIPLOID_ALLELE)||(ob==UNKNOWN_DIPLOID_ALLELE)) {unimp++; u++;}
                    else if(ib==GAP_DIPLOID_ALLELE) {gaps++;}
                    else if(ib==ob) {
                        correct++;
                        c++;
                    } else {
                        if(isHeterozygous(ob)||isHeterozygous(ib)) {hets++; h++;}
                        else {errors++;
                            e++;
//                            if(t==0) System.out.printf("%d %d %s %s %n",t,s,oA.getBaseAsString(oATaxa, s), iA.getBaseAsString(t, s));
                        }
                    }
                }
            }
            if(taxaOut) System.out.printf("%s %d %d %d %d %n",iA.taxaName(t),u,h,c,e);
        }
        System.out.println("MFile\tIFile\tGap\tUnimp\tUnimpHets\tCorrect\tErrors");
        System.out.printf("%s\t%s\t%d\t%d\t%d\t%d\t%d%n",maskFile, impFile, gaps, unimp,hets,correct,errors);
        return new int[]{gaps, unimp,hets,correct,errors};
    }

    @Override
    public void setParameters(String[] args) {
        if (args.length == 0) {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }

        engine.add("-hmp", "-hmpFile", true);
        engine.add("-o", "--outFile", true);
        engine.add("-d", "--donorH", true);
        engine.add("-maskKeyFile", "--maskKeyFile", true);
        engine.add("-propSitesMask", "--propSitesMask", true);
        engine.add("-mxHet", "--hetThresh", true);
        engine.add("-sC", "--startChrom", false); //TODO why set to false
        engine.add("-eC", "--endChrom", false); //TODO why set to false
        engine.add("-minMnCnt", "--minMnCnt", true);
        engine.add("-mxInbErr", "--mxInbErr", true);
        engine.add("-mxHybErr", "--mxHybErr", true);
        engine.add("-mxVitFocusErr", "--mxVitFocusErr", true);
        engine.add("-mxInbFocusErr", "--mbIndFocusErr", true);
        engine.add("-mxComFocusErr", "--mxComFocusErr", true);
        engine.add("-mxInbFocusErrHet", "--mxInbFocusErrHet", true);
        engine.add("-mxComFocusErrHet", "--mxComFocusErrHet", true);
        engine.add("-hybNNOff", "--hybNNOff", true);
        engine.add("-mxDonH", "--mxDonH", true);
        engine.add("-mnTestSite", "--mnTestSite", true);
        engine.add("-projA", "--projAlign", false);
        engine.add("-runChrMode", "--runChrMode", false);
        engine.add("-nV", "--nonVerbose",false);
        engine.parse(args);
        hmpFile = engine.getString("-hmp");
        outFileBase = engine.getString("-o");
        donorFile = engine.getString("-d");
        maskKeyFile = engine.getString("-maskKeyFile");
        if(engine.getBoolean("-mxHet")) {
            hetThresh = Double.parseDouble(engine.getString("-mxHet"));
        }
        if(engine.getBoolean("-mxHet")) {
            propSitesMask = Double.parseDouble(engine.getString("-propSitesMask"));
        }
        if (engine.getBoolean("-mxInbErr")) {
            maximumInbredError = Double.parseDouble(engine.getString("-mxInbErr"));
        }
        if (engine.getBoolean("-mxHybErr")) {
            maxHybridErrorRate = Double.parseDouble(engine.getString("-mxHybErr"));
        }
        if (engine.getBoolean("-mxVitFocusErr")) {
            maxHybridErrFocusHomo = Double.parseDouble(engine.getString("-mxVitFocusErr"));
        }
        if (engine.getBoolean("-mxInbFocusErr")) {
            maxInbredErrFocusHomo = Double.parseDouble(engine.getString("-mxInbFocusErr"));
        }
        if (engine.getBoolean("-mxComFocusErr")) {
            maxSmashErrFocusHomo = Double.parseDouble(engine.getString("-mxComFocusErr"));
        }
        if (engine.getBoolean("-mxInbFocusErrHet")) {
            maxInbredErrFocusHet = Double.parseDouble(engine.getString("-mxInbFocusErrHet"));
        }
        if (engine.getBoolean("-mxComFocusErrHet")) {
            maxSmashErrFocusHet = Double.parseDouble(engine.getString("-mxComFocusErrHet"));
        }
        if (engine.getBoolean("-minMnCnt")) {
            minMinorCnt = Integer.parseInt(engine.getString("-minMnCnt"));
        }
        if (engine.getBoolean("-hybNNOff")) hybridNN=false;
        if (engine.getBoolean("-mxDonH")) {
            maxDonorHypotheses = Integer.parseInt(engine.getString("-mxDonH"));
        }
        if (engine.getBoolean("-mnTestSite")) {
            minTestSites = Integer.parseInt(engine.getString("-mnTestSite"));
        }
        if (engine.getBoolean("-projA")) isOutputProjection=true;
        if (engine.getBoolean("-nV")) verboseOutput=false;
    }



    private void printUsage() {
        myLogger.info(
                "\n\n\nAvailable options for the FILLINImputationPlugin are as follows:\n"
                        + "-hmp   Input HapMap file of target genotypes to impute. Accepts all file types supported by TASSEL5\n"
                        + "-d    Donor haplotype files from output of FILLINFindHaplotypesPlugin. Use “.gX” in the input filename to denote the substring “.gc#s#” found in donor files\n"
                        + "-o     Output file; hmp.txt.gz and .hmp.h5 accepted. Required\n"
                        + "-maskKeyFile An optional key file to indicate that file is already masked for accuracy calculation. Non-missing genotypes indicate masked sites. Else, will generate own mask\n"
                        + "-propSitesMask   The proportion of non missing sites to mask for accuracy calculation if depth is not available (default:"+propSitesMask+"\n"
                        + "-mxHet   Threshold per taxon heterozygosity for treating taxon as heterozygous (no Viterbi, het thresholds). (default:"+hetThresh+"\n"
                        + "-minMnCnt    Minimum number of informative minor alleles in the search window (or "+minMajorRatioToMinorCnt+"X major)\n"
                        + "-mxInbErr    Maximum error rate for applying one haplotype to entire site window (default:"+maximumInbredError+"\n"
                        + "-mxHybErr    Maximum error rate for applying Viterbi with to haplotypes to entire site window (default:"+maxHybridErrorRate+"\n"
                        + "-mxVitFocusErr    Maximum error rate for applying Viterbi with to haplotypes to entire site window  (default:"+maxHybridErrFocusHomo+")\n"
                        + "-mxInbFocusErr    Maximum error rate to apply one haplotype for inbred (heterozygosity below mxHet) taxon  for 64-site focus blocks (default:"+maxInbredErrFocusHomo+")\n"
                        + "-mxComFocusErr    Maximum error rate to apply two haplotypes modeled as a heterozygote for inbred (heterozygosity below mxHet) taxon  for 64-site focus blocks (default:"+maxSmashErrFocusHomo+")\n"
                        + "-mxInbFocusErrHet    Maximum error rate to apply one haplotype for outbred (heterozygosity above mxHet) taxon  for 64-site focus blocks (default:"+maxInbredErrFocusHet+")\n"
                        + "-mxComFocusErrHet    Maximum error rate to apply two haplotypes modeled as a heterozygote for outbred (heterozygosity above mxHet) taxon  for 64-site focus blocks (default:"+maxSmashErrFocusHet+")\n"
                        + "-hybNNOff    Whether to model two haplotypes as heterozygotic for focus blocks (default:"+hybridNN+")\n"
                        + "-mxDonH   Maximum number of donor hypotheses to be explored (default: "+maxDonorHypotheses+")\n"
                        + "-mnTestSite   Minimum number of sites to test for IBS between haplotype and target in focus block  (default:"+minTestSites+")\n"
                        + "-projA   Create a projection alignment for high density markers (default off)\n"
        );
    }

    @Override
    public DataSet performFunction(DataSet input) {
        runMinorWindowViterbiImputation(donorFile, hmpFile, outFileBase,
                minMinorCnt, minTestSites, 100, maxHybridErrorRate, isOutputProjection, false);
        return null;
    }



    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "ImputeByFILLIN";
    }

    @Override
    public String getToolTipText() {
        return "Imputation that relies on a combination of HMM and Nearest Neighbor";
    }
}

