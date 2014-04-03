package net.maizegenetics.analysis.imputation;

import net.maizegenetics.analysis.distance.IBSDistanceMatrix;
import net.maizegenetics.analysis.popgen.DonorHypoth;
import net.maizegenetics.dna.map.DonorHaplotypes;
import net.maizegenetics.dna.map.PositionList;
import net.maizegenetics.dna.snp.*;
import net.maizegenetics.dna.snp.io.ProjectionGenotypeIO;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.prefs.TasselPrefs;
import net.maizegenetics.util.ArgsEngine;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;
import net.maizegenetics.util.ProgressListener;
import org.apache.commons.lang.ArrayUtils;
import org.apache.log4j.Logger;

import javax.swing.*;
import java.awt.*;
import java.io.File;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import static net.maizegenetics.dna.snp.GenotypeTable.UNKNOWN_DIPLOID_ALLELE;
import static net.maizegenetics.dna.snp.GenotypeTable.WHICH_ALLELE.Major;
import static net.maizegenetics.dna.snp.GenotypeTable.WHICH_ALLELE.Minor;
import static net.maizegenetics.dna.snp.GenotypeTableUtils.isHeterozygous;
import static net.maizegenetics.dna.snp.NucleotideAlignmentConstants.GAP_DIPLOID_ALLELE;



/**
 * FILLIN imputation relies on a libary of haplotypes and uses nearest neighbor searches
 * followed by HMM Viterbi resolution or block-based resolution.
 * It is the best approach for substantially unrelated taxa in TASSEL.  BEAGLE4 is a better
 * approach currently for landraces, while FILLIN outperforms if there is a good reference
 * set of haplotypes.
 * <p></p>
 * The algorithm relies on having a series of donor haplotypes.  If phased haplotypes are
 * already known, they can be used directly.  If they need to be reconstructed, the
 * {@link FindMergeHaplotypesPlugin} can be used to create haplotypes for windows
 * across the genome.
 * <p></p>
 * Imputation is done one taxon at a time using the donor haplotypes.  The strategy is as
 *  follows:
 * <li></li>Every 64 sites is considered a block (or a word in the bitset terminology - length of a long).  There is
 * a separate initial search for nearest neighbors centered around each block. The flanks
 * of the scan window vary but always contain the number of minor alleles specified (a key parameter).
 * Minor alleles approximate the information content of the search.
 * <li></li> Calculate distance between the target taxon and the donor haplotypes for every window, and rank
 * the top 10 donor haplotypes for each 64 site focus block.
 * <li></li> Evaluate by Viterbi whether 1 or 2 haplotypes will explain all the sites across the
 * donor haplotype window.  If successful, move to the next region for donor haplotypes.
 * <li></li> Resolve each focus block by nearest neighbors.  If inbred NN are not close enough,
 * then do a hybrid search seeded with one parent being the initial 10 haplotypes.
 * Set bases based on what is resolved first inbred or hybrid.
 *
 *<p></p>
 * Error rates are bounded away from zero, but adding 0.5 error to all error
 * rates that that were observed to be zero.
 * <p></p>

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
    private double minimumHybridDonorDistance=maximumInbredError*5;

    //options for focus blocks
    private double maxHybridErrFocusHomo= .3333*maxHybridErrorRate;////max error rate for discrepacy between two haplotypes for the focus block. it's default is higher because calculating for fewer sites
    private double maxInbredErrFocusHomo= .3*maximumInbredError;//.003;
    private double maxSmashErrFocusHomo= maximumInbredError;//.01;
    private double maxInbredErrFocusHet= .1*maximumInbredError;//.001;//the error rate for imputing one haplotype in focus block for a het taxon
    private double maxSmashErrFocusHet= maximumInbredError;//.01;
    private double hetThresh= 0.02;//threshold for whether a taxon is considered heterozygous
    
    //options for masking and calculating accuracy
    private boolean accuracyOn= true;
    private String maskKeyFile= null;
    private double propSitesMask= .001;
    private GenotypeTable maskKey= null;
    private double[] MAFClass= null;//new double[]{0,.02,.05,.10,.20,.3,.4,.5,1};
        
    
    public static GenotypeTable unimpAlign;  //the unimputed alignment to be imputed, unphased
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

    //initialize the emission matrix, states (5) in rows, observations (3) in columns
    double[][] emission = new double[][] {
                    {.998,.001,.001},
                    {.6,.2,.2},
                    {.4,.2,.4},
                    {.2,.2,.6},
                    {.001,.001,.998}
    };


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
        GenotypeTable[] donorAlign= null;
        try {
            if (donorFile.contains(".gX")) {donorAlign=loadDonors(donorFile, unimpAlign, minTestSites, verboseOutput);}
            else if(true) {donorAlign=loadDonors(donorFile, unimpAlign, minTestSites, 4000,true);}  //todo add option for size
            //todo break into lots of little alignments
            else throw new IllegalArgumentException();
        }
        catch (IllegalArgumentException e){
            throw new IllegalArgumentException("Incorrect donor file supplied. Must contain '.gX' within the file name,\nand match a set of files output from FILLINFindHaplotypesPlugin()");
        }

        OpenBitSet[][] conflictMasks=createMaskForAlignmentConflicts(unimpAlign, donorAlign, verboseOutput);

        siteErrors=new int[unimpAlign.numberOfSites()];
        siteCorrectCnt=new int[unimpAlign.numberOfSites()];
        taxonErrors=new int[unimpAlign.numberOfTaxa()];
        taxonCorrectCnt=new int[unimpAlign.numberOfTaxa()];
        
        System.out.printf("Unimputed taxa:%d sites:%d %n",unimpAlign.numberOfTaxa(),unimpAlign.numberOfSites());
        System.out.println("Creating Export GenotypeTable:"+exportFile);
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
        characterizeHeterozygosity();
        int numThreads = Runtime.getRuntime().availableProcessors();
        System.out.println("Time to read in files and generate masks: "+((System.currentTimeMillis()-time)/1000)+" sec");
        ExecutorService pool = Executors.newFixedThreadPool(numThreads);
        
        for (int taxon = 0; taxon < unimpAlign.numberOfTaxa(); taxon+=1) {
            int[] trackBlockNN= new int[5];//global variable to track number of focus blocks solved in NN search for system out; index 0 is inbred, 1 is viterbi, 2 is smash, 3 is not solved, 4 is total for all modes
            ImputeOneTaxon theTaxon= (((double)unimpAlign.heterozygousCountForTaxon(taxon)/(double)unimpAlign.totalNonMissingForTaxon(taxon))<hetThresh)?
                new ImputeOneTaxon(taxon, donorAlign, minSitesPresent, conflictMasks,imputeDonorFile, mna, trackBlockNN, maxInbredErrFocusHomo, maxHybridErrFocusHomo, maxSmashErrFocusHomo, true):
                    new ImputeOneTaxon(taxon, donorAlign, minSitesPresent, conflictMasks,imputeDonorFile, mna, trackBlockNN, maxInbredErrFocusHet, 0, maxSmashErrFocusHet, false);
            //        theTaxon.run(); //retained to provide a quick way to debug.  uncomment to help in debugging.
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
        
        double runtime= (double)(System.currentTimeMillis()-time)/(double)1000;
        if(isOutputProjection) {
            ProjectionGenotypeIO.writeToFile(exportFile, ((ProjectionBuilder) mna).build());
        } else {
            GenotypeTableBuilder ab=(GenotypeTableBuilder)mna;
            if(ab.isHDF5()) {
                ab.build();
            } else {
                ExportUtils.writeToHapmap(ab.build(), false, exportFile, '\t', null);
            }
        }
        System.out.printf("%d %g %d %n",minMinorCnt, maximumInbredError, maxDonorHypotheses);
        System.out.println("Time to read in files, impute target genotypes, and calculate accuracy: "+runtime+" seconds");
    }

    private class ImputeOneTaxon implements Runnable{
        int taxon;
        GenotypeTable[] donorAlign;
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
        
        public ImputeOneTaxon(int taxon, GenotypeTable[] donorAlign, int minSitesPresent, OpenBitSet[][] conflictMasks,
            boolean imputeDonorFile, Object mna, int[] trackBlockNN, double focusInbErr, double focusHybridErr, double focusSmashErr, boolean hetsToMissing) {
            this.taxon=taxon;
            this.donorAlign=donorAlign;
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
            int countFullLength=0;
            int countByFocus= 0;
            for (int da = 0; (da < donorAlign.length)&&enoughData ; da++) {
                int donorOffset=unimpAlign.siteOfPhysicalPosition(donorAlign[da].chromosomalPosition(0), donorAlign[da].chromosome(0));
                int blocks=donorAlign[da].allelePresenceForAllSites(0, Major).getNumWords();
                BitSet[] maskedTargetBits=arrangeMajorMinorBtwAlignments(unimpAlign, taxon, donorOffset,
                        donorAlign[da].numberOfSites(),conflictMasks[da][0],conflictMasks[da][1]);

                //if imputing the donor file, these donor indices prevent self imputation
                int[] donorIndices;
                if(imputeDonorFile){
                    donorIndices=new int[donorAlign[da].numberOfTaxa()-1];
                    for (int i = 0; i < donorIndices.length; i++) {donorIndices[i]=i; if(i>=taxon) donorIndices[i]++;}
                } else {
                    donorIndices=new int[donorAlign[da].numberOfTaxa()];
                    for (int i = 0; i < donorIndices.length; i++) {donorIndices[i]=i;}
                }

                //Finds the best haplotype donors for each focus block within a donorGenotypeTable
                DonorHypoth[][] regionHypthInbred=new DonorHypoth[blocks][maxDonorHypotheses];
                byte[][][] targetToDonorDistances=FILLINImputationUtils.calcAllelePresenceCountsBtwTargetAndDonors(maskedTargetBits,
                        donorAlign[da]);
                for (int focusBlock = 0; focusBlock < blocks; focusBlock++) {
                    int[] resultRange=getBlockWithMinMinorCount(maskedTargetBits[0].getBits(),maskedTargetBits[1].getBits(),
                            focusBlock, minMinorCnt);
                    if(resultRange==null) continue; //no data in the focus Block
                    //search for the best inbred donors for a segment
                    regionHypthInbred[focusBlock]=FILLINImputationUtils.findHomozygousDonorHypoth(taxon, resultRange[0], resultRange[2],
                            focusBlock, donorIndices, targetToDonorDistances, minTestSites, maxDonorHypotheses);
                }
                impTaxon.setSegmentSolved(false);
                //tries to solve the entire donorAlign region with 1 or 2 donor haplotypes by Virterbi
                impTaxon=apply1or2Haplotypes(taxon, donorAlign[da], donorOffset, regionHypthInbred,  impTaxon, maskedTargetBits,
                        maxHybridErrorRate, targetToDonorDistances);
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
            double totalFocus= (double)trackBlockNN[3]+(double)trackBlockNN[4];
            sb.append(String.format("InbredOrViterbi:%d FocusBlock:%d PropFocusInbred:%f PropFocusViterbi:%f PropFocusSmash:%f PropFocusMissing:%f BlocksSolved:%d ",
                    countFullLength, countByFocus, (double)trackBlockNN[0]/totalFocus, (double)trackBlockNN[1]/totalFocus, (double)trackBlockNN[2]/totalFocus, (double)trackBlockNN[3]/totalFocus, impTaxon.getBlocksSolved()));
            //sb.append(String.format("InbredOrViterbi:%d FocusBlock:%d BlocksSolved:%d ", countFullLength, countByFocus, impTaxon.getBlocksSolved()));
            int[] unk=countUnknownAndHets(impTaxon.resolveGeno);
            sb.append(String.format("Unk:%d PropMissing:%g ", unk[0], (double) unk[0] / (double) impTaxon.getOrigGeno().length));
            sb.append(String.format("Het:%d PropHet:%g ", unk[1], (double)unk[1]/(double)impTaxon.getOrigGeno().length));
            sb.append(" BreakPoints:"+impTaxon.getBreakPoints().size());
            if(!isOutputProjection) {
                alignBuilder.addTaxon(unimpAlign.taxa().get(taxon), impTaxon.resolveGeno);
            } else {
                projBuilder.addTaxon(unimpAlign.taxa().get(taxon),impTaxon.getBreakPoints());
            }
            if(verboseOutput) System.out.println(sb.toString());
        }
    }

    private void characterizeHeterozygosity() {
        double notMissing= 0;
        double het= 0;
        for (int taxon = 0; taxon < unimpAlign.numberOfTaxa(); taxon++) {
            het+= (double)unimpAlign.heterozygousCountForTaxon(taxon);
            notMissing+= (double)unimpAlign.totalNonMissingForTaxon(taxon);}
        double avgHet= het/notMissing;
        double sqDif= 0;
        for (int taxon = 0; taxon < unimpAlign.numberOfTaxa(); taxon++) {
            double x= avgHet-((double)unimpAlign.heterozygousCountForTaxon(taxon)/(double)unimpAlign.totalNonMissingForTaxon(taxon));
            sqDif+= (x*x);}
        System.out.println("Average Heterozygosity: "+avgHet+" plus/minus std error: "+Math.sqrt(sqDif/(unimpAlign.numberOfTaxa()-1))/Math.sqrt(unimpAlign.numberOfTaxa()));

    }

    private static GenotypeTable[] loadDonors(String donorFileRoot, GenotypeTable unimpAlign, int minTestSites, boolean verboseOutput){
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
            //for (int j = 0; j < donorAlign[i].numberOfSites(); j++) {if (donU.indexOf(donP.get(j)) > -1) subSites.add(j);} //if unimputed contains donorAlign position keep in donor align
            for (int j = 0; j < donorAlign[i].numberOfSites(); j++) {if (donU.siteOfPhysicalPosition(donP.physicalPositions()[j], donP.chromosome(j)) > -1) subSites.add(j);} //if unimputed contains donorAlign position keep in donor align
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
        return donorAlign;
    }

    private static GenotypeTable[] loadDonors(String donorFile, GenotypeTable unimpAlign, int minTestSites, int appoxSitesPerHaplotype, boolean verboseOutput){
        GenotypeTable donorMasterGT=ImportUtils.readGuessFormat(donorFile);
        donorMasterGT=GenotypeTableBuilder.getHomozygousInstance(donorMasterGT);
        int[][] donorFirstLastSites=FILLINFindHaplotypesPlugin.divideChromosome(donorMasterGT,appoxSitesPerHaplotype,true);
        PositionList donU= unimpAlign.positions();
        GenotypeTable[] donorAlign=new GenotypeTable[donorFirstLastSites.length];
        for (int i = 0; i < donorAlign.length; i++) {
            if(verboseOutput) System.out.println("Starting Read");
            donorAlign[i]=FilterGenotypeTable.getInstance(donorMasterGT,donorFirstLastSites[i][0],donorFirstLastSites[i][1]);
            }
        return donorAlign;
    }

    /**
     * Solve the entire donor alignment window with either one or two haplotypes.  Generally, it is solved with the
     * HMM Viterbi algorithm.  The Viterbi algorithm is relatively slow the big choice is how to choose the potential
     * donors.
     *
     *
     * @param taxon
     * @param donorAlign
     * @param donorOffset
     * @param regionHypoth
     * @param impT
     * @param maskedTargetBits
     * @param maxHybridErrorRate
     * @return
     */
    private ImputedTaxon apply1or2Haplotypes(int taxon, GenotypeTable donorAlign, int donorOffset,
                                             DonorHypoth[][] regionHypoth, ImputedTaxon impT,
                                             BitSet[] maskedTargetBits, double maxHybridErrorRate, byte[][][] targetToDonorDistances) {

        int blocks=maskedTargetBits[0].getNumWords();
        if(testing==1) System.out.println("Starting complete hybrid search");
        //create a list of the best donors based on showing up frequently high in many focus blocks
        //in the test data set the best results achieved with on the best hypothesis recovered.
//        int[] d=FILLINImputationUtils.mostFrequentDonorsAcrossFocusBlocks(regionHypoth, maxDonorHypotheses);
        //Alternative test is find best donors
       int[] d=FILLINImputationUtils.bestDonorsAcrossEntireRegion(targetToDonorDistances, minTestSites,maxDonorHypotheses);
        int[] testList=FILLINImputationUtils.fillInc(0,donorAlign.numberOfTaxa()-1);
        int[] bestDonorList=Arrays.copyOfRange(d,0,Math.min(d.length,5));
        DonorHypoth[] bestDBasedOnBest=FILLINImputationUtils.findHeterozygousDonorHypoth(taxon, maskedTargetBits[0].getBits(),
                maskedTargetBits[1].getBits(), 0, blocks-1, blocks/2, donorAlign, bestDonorList, testList, maxDonorHypotheses, minTestSites);

        //make all combination of best donor and find the the pairs that minimize errors
        //with the true switch also will make inbreds
        DonorHypoth[] best2Dsearchdonors=FILLINImputationUtils.findHeterozygousDonorHypoth(taxon, maskedTargetBits[0].getBits(), maskedTargetBits[1].getBits(), 0, blocks-1, blocks/2, donorAlign, d, d, maxDonorHypotheses, minTestSites);
        DonorHypoth[] best2donors=FILLINImputationUtils.combineDonorHypothArrays(maxDonorHypotheses,bestDBasedOnBest,best2Dsearchdonors);
        if(testing==1) System.out.println(Arrays.toString(best2donors));
        ArrayList<DonorHypoth> goodDH=new ArrayList<DonorHypoth>();
        for (DonorHypoth dh : best2donors) {
            if((dh!=null)&&(dh.getErrorRate()<maxHybridErrorRate)) {
                if(dh.isInbred()==false){  //if not inbred then Virterbi is needed to determine segments
                    dh=getStateBasedOnViterbi(dh, donorOffset, donorAlign, twoWayViterbi, transition);
                }
                if(dh!=null) goodDH.add(dh);
            }
        }
        if(goodDH.isEmpty()) {
//            System.out.printf("%s %s SS:%d LS:%d %s %s %n",unimpAlign.taxaName(taxon),donorAlign.chromosome(0).getName(),
//                    donorAlign.chromosomalPosition(0),donorAlign.chromosomalPosition(donorAlign.numberOfSites()-1),
//                    Arrays.toString(d), Arrays.toString(d1));
//            for (DonorHypoth donorHypoth : bestDBasedOnBest) {
//                System.out.println(donorHypoth.toString());
//            }
            return impT;
        }
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
            int donorOffset=unimpAlign.positions().siteOfPhysicalPosition(donorAlign[da].positions().physicalPositions()[0], donorAlign[da].positions().chromosome(0));
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

    private int[] countUnknownAndHets(byte[] a) {
        int cnt=0, cntHets=0;
        for (int i = 0; i < a.length; i++) {
            if(a[i]==UNKNOWN_DIPLOID_ALLELE) {cnt++;}
            else if(isHeterozygous(a[i])) {cntHets++;}
        }
        return new int[]{cnt,cntHets};
    }

    /*
    Major and minor allele get swapped between alignment.  This method flips these, and sets the bad sites to missing
     */
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
        //Selects only the informatives sites for the Viterbi algorithm
        //Perhaps hets should be removed from the donor
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
            if(isHeterozygous(d1b[cs]) || isHeterozygous(d2b[cs])) continue;
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
        if(testing>0) System.out.println("Input:"+Arrays.toString(informStatesF)+" Swaps"+countSwaps(informStatesF));
        byte[] resultStatesF=vaF.getMostProbableStateSequence();
        if(testing>0) System.out.println("Resul:"+Arrays.toString(resultStatesF)+" Swaps"+countSwaps(resultStatesF));
        DonorHypoth dh2=new DonorHypoth(dh.targetTaxon,dh.donor1Taxon,
                dh.donor2Taxon, dh.startBlock, dh.focusBlock, dh.endBlock);
        int currPos=0;
        for(int cs=0; cs<sites; cs++) {  //converts to informative states back to all states
            callsF[cs]=(resultStatesF[currPos]==1)?(byte)1:(byte)(resultStatesF[currPos]/2); //converts the scale back to 0,1,2 from 0..4
            if((pos[currPos]<cs+startSite)&&(currPos<resultStatesF.length-1)) currPos++;
        }
        if(testing>1) System.out.println("callsF:"+Arrays.toString(callsF)+" Swaps"+countSwaps(callsF));

        if (forwardReverse==true) {
            TransitionProbability tpR = new TransitionProbability();
            tpR.setTransitionProbability(transition);
            byte[] informStatesR=Arrays.copyOf(informStatesF,informStatesF.length);
            ArrayUtils.reverse(informStatesR);
            int[] posR= Arrays.copyOf(pos,pos.length);
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
            if(testing>1) System.out.println("callsR:"+Arrays.toString(callsR)+" Swaps"+countSwaps(callsR));
            //compare the forward and reverse viterbi, use the one with the longest path length if they contradict
            byte[] callsC=Arrays.copyOf(callsF,callsF.length);  //this does not copy it is just a pointer assignment
            for(int i= 0;i<pos.length;i++) {
                int cs= pos[i]-startSite;
                if (callsF[cs]!=callsR[cs]&&i<pos.length/2) callsC[cs]= callsR[cs];
            }
            if (testing>0) {
                if (resultStatesF[0]!=resultStatesR[0]||resultStatesF[resultStatesF.length-1]!=resultStatesR[resultStatesR.length-1]) {
                    System.out.println("FR:\n"+Arrays.toString(informStatesF)+"\n"+Arrays.toString(informStatesR)+"\n"+
                            Arrays.toString(resultStatesF)+"\n"+Arrays.toString(resultStatesR)+"\n"+Arrays.toString(callsF)+"\n"+
                            Arrays.toString(callsR)+"\n"+Arrays.toString(callsC));
                }
            }
            if(testing>1) System.out.println("callsC:"+Arrays.toString(callsC)+" Swaps"+countSwaps(callsC));
            if(testing>1) System.out.println("DonorOffset:"+donorOffset+" Reverse:"+countSwaps(callsC)+" dh:"+dh.toString());
            dh2.phasedResults= callsC;
            return dh2;
        }
        System.out.println("Forward:"+countSwaps(callsF)+" dh:"+dh.toString());
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
     * the number of unmatched informative alleles.  Currently, there is little tie
     * breaking, longer matches are favored.
     * @param targetTaxon
     * @param startBlock
     * @param endBlock
     * @param focusBlock
     * @return int[] array of {donor1, donor2, testSites}  sorted by  testPropUnmatched
     */
    private DonorHypoth[] getBestHybridDonors(int targetTaxon, long[] mjT, long[] mnT,
                                              int startBlock, int endBlock, int focusBlock, GenotypeTable donorAlign, int[] donor1Indices, int[] donor2Indices,
                                              boolean viterbiSearch) {
        TreeMap<Double,DonorHypoth> bestDonors=new TreeMap<Double,DonorHypoth>();
        Set<Long> testedDonorPairs=new HashSet<>();
        double lastKeytestPropUnmatched=1.0;
        double inc=1e-9;
        double[] donorDist;
        for (int d1: donor1Indices) {
            long[] mj1=donorAlign.allelePresenceForSitesBlock(d1, Major, startBlock, endBlock + 1);
            long[] mn1=donorAlign.allelePresenceForSitesBlock(d1, Minor, startBlock, endBlock + 1);
            for (int d2 : donor2Indices) {
                if((!viterbiSearch)&&(d1==d2)) continue;
                if(testedDonorPairs.contains(((long)d1<<32)+(long)d2)) continue;
                if(testedDonorPairs.contains(((long)d2<<32)+(long)d1)) continue;
                testedDonorPairs.add(((long)d1<<32)+(long)d2);
                long[] mj2=donorAlign.allelePresenceForSitesBlock(d2, Major, startBlock, endBlock + 1);
                long[] mn2=donorAlign.allelePresenceForSitesBlock(d2, Minor, startBlock, endBlock + 1);
                if(viterbiSearch) {
                    donorDist=IBSDistanceMatrix.computeHetBitDistances(mj1, mn1, mj2, mn2, minTestSites);
                    if((d1!=d2)&&(donorDist[0]<this.minimumHybridDonorDistance)) continue;
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
        //prevDonors is used to track recombination breakpoint start and stop
        int[] prevDonors=new int[]{-1, -1};
        if(theDH[0].getPhaseForSite(startSite)==0) {prevDonors=new int[]{theDH[0].donor1Taxon, theDH[0].donor1Taxon};}
        else if(theDH[0].getPhaseForSite(startSite)==2) {prevDonors=new int[]{theDH[0].donor2Taxon, theDH[0].donor2Taxon};}
        else if(theDH[0].getPhaseForSite(startSite)==1) {prevDonors=new int[]{theDH[0].donor1Taxon, theDH[0].donor2Taxon};}
        int prevDonorStart=startSite;

        int[] currDonors=prevDonors;

        for(int cs=startSite; cs<=endSite; cs++) {
            byte donorEst=UNKNOWN_DIPLOID_ALLELE;
            byte neighbor=0;
            //site imputation will try to use all donor hypotheses; breakpoints only use the best (e.g. i==0 tests)
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
            //record recombination breakpoint if it changes
            if(!Arrays.equals(prevDonors, currDonors)) {
                DonorHaplotypes dhaps=new DonorHaplotypes(donorAlign.chromosome(prevDonorStart), donorAlign.chromosomalPosition(prevDonorStart),
                        donorAlign.chromosomalPosition(cs),prevDonors[0],prevDonors[1]);
                impT.addBreakPoint(dhaps);
                prevDonors=currDonors;
                prevDonorStart=cs;
            }
            //chgHis records the section of the imputation pipeline that made the change, Viterbi negative, blockNN positive
            if(theDH[0].phasedResults==null) {impT.chgHis[cs+donorOffset]=(byte)-neighbor;}
            else {impT.chgHis[cs+donorOffset]=(byte)neighbor;}

            impT.impGeno[cs+donorOffset]= donorEst;  //predicted based on neighbor
            //if genotype is unknown or het undercall then resolves
            //todo there is currently not an option if the predicted disagrees with a homozygous call
            if(knownBase==UNKNOWN_DIPLOID_ALLELE) {
                if (isHeterozygous(donorEst)) {
                    if (smashOn && hetsMiss) {//if imputing a heterozygote, just set to missing
                        impT.resolveGeno[cs+donorOffset]= knownBase;
                    }
                    else impT.resolveGeno[cs+donorOffset]= donorEst; //if not in hybrid, set to het
                }
                else {//if not imputed to a het
                    impT.resolveGeno[cs+donorOffset]= donorEst;}
            } else if (isHeterozygous(donorEst)){
                if(resolveHetIfUndercalled&&GenotypeTableUtils.isPartiallyEqual(knownBase, donorEst)&&smashOn==false){//if smash off, set homozygotes imputed to het to het
                    //System.out.println("ResolveHet:"+theDH[0].targetTaxon+":"+cs+donorOffset+":"+knownBaseString+":"+donorEstString);
                    impT.resolveGeno[cs+donorOffset]= donorEst;
                }
            }
        } //end of cs loop
        //enter a stop of the DH at the beginning of the next block
        DonorHaplotypes dhaps=new DonorHaplotypes(donorAlign.chromosome(prevDonorStart), donorAlign.chromosomalPosition(prevDonorStart),
                donorAlign.chromosomalPosition(endSite),prevDonors[0],prevDonors[1]);
        impT.addBreakPoint(dhaps);
        return impT;

    }

    private int countSwaps(byte[] swaps) {
        int nSwap=0;
        byte lastValue=swaps[0];
        for (byte swap : swaps) {
            if(lastValue!=swap) {
                nSwap++;
                lastValue=swap;
            }
        }
        return nSwap;
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

    @Override
    public void setParameters(String[] args) {
        if (args.length == 0) {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }

        engine.add("-hmp", "-hmpFile", true);
        engine.add("-o", "--outFile", true);
        engine.add("-d", "--donorH", true);
        engine.add("-accuracyOff", "--accuracyOff", true);
        engine.add("-maskKeyFile", "--maskKeyFile", true);
        engine.add("-propSitesMask", "--propSitesMask", true);
        engine.add("-mxHet", "--hetThresh", true);
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
        if(engine.getBoolean("-propSitesMask")) {
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
        if (engine.getBoolean("-accuracyOff")) accuracyOn=false;
        if (engine.getBoolean("-projA")) isOutputProjection=true;
        if (engine.getBoolean("-nV")) verboseOutput=false;
    }



    private void printUsage() {
        myLogger.info(
                "\n\n\nAvailable options for the FILLINImputationPlugin are as follows:\n"
                        + "-hmp   Input HapMap file of target genotypes to impute. Accepts all file types supported by TASSEL5\n"
                        + "-d    Donor haplotype files from output of FILLINFindHaplotypesPlugin. Use .gX in the input filename to denote the substring .gc#s# found in donor files\n"
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
                        + "-accuracyOff   Do not calculate accuracy for imputation (default on)\n"
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

