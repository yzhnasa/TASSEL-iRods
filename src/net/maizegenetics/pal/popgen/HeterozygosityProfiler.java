/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.pal.popgen;


import java.util.Arrays;
import java.util.Random;
import net.maizegenetics.gwas.imputation.EmissionProbability;
import net.maizegenetics.gwas.imputation.TransitionProbability;
import net.maizegenetics.gwas.imputation.ViterbiAlgorithm;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.ImportUtils;
import net.maizegenetics.util.ProgressListener;

/**
 *
 * @author edbuckler
 */
public class HeterozygosityProfiler {
    private Alignment tA;
    private float[][] blockHetProb=null;

    public HeterozygosityProfiler(String unImpTargetFile, String exportFile) {
   
        tA=ImportUtils.readFromHapmap(unImpTargetFile, false, (ProgressListener)null);
        tA.optimizeForTaxa(null);
        tA.optimizeForSites(null);
        int blocks=tA.getAllelePresenceForAllSites(0, 0).getNumWords();
        blockHetProb=new float[tA.getSequenceCount()][blocks];
        
        int totalBlocks=0, totalhetBlocks=0;
        TransitionProbability tp=setTransitionProb(tA, blocks);
        double[][] hetExp=simulateHetInformation(-1, true);
        double[][] homoWithErrorExp=simulateHetInformation(0.002, false);
        for (int bt = 0; bt < tA.getSequenceCount(); bt+=1) {
            int taxaImpCnt=0;
            //DonorHypoth[][] regionHypth=new DonorHypoth[blocks][10];
            String name=tA.getIdGroup().getIdentifier(bt).getFullName();
  //          System.out.printf("Imputing %d:%s ... %n", bt,name);
            long time=System.currentTimeMillis();
            long[] mjT=tA.getAllelePresenceForAllSites(bt, 0).getBits();
            long[] mnT=tA.getAllelePresenceForAllSites(bt, 1).getBits();
            int[][] hetCnt=profileTaxaForHets(mjT, mnT);
            int sum=0, hsum=0;
            for (int i = 0; i < hetCnt[0].length; i++) {
                sum+=hetCnt[0][i];
                totalBlocks++;
                hsum+=hetCnt[1][i];
                if(hetCnt[1][i]>0) totalhetBlocks++;
            }
            int avgCovPerBlock=sum/hetCnt[0].length;
            byte[] obs = new byte[blocks];
            byte[] obsRev = new byte[blocks];
            for (int b = 0; b < blocks; b++) {
                obs[b]=(hetCnt[1][b]>0)?(byte)1:(byte)0;
                obs[b]=(byte)hetCnt[1][b];
                obsRev[blocks-b-1]=obs[b];
            }
            System.out.println("avgCovPerBlock:"+avgCovPerBlock);
            EmissionProbability ep=setEmmisionProb(hetExp[avgCovPerBlock],homoWithErrorExp[avgCovPerBlock]);  //in the future needs to include coverage poisson
            double expectedPhet=0.8;
            double[] initialProb = new double[]{1 - expectedPhet, expectedPhet};
        //    double[] initialProb = new double[]{1,0};
            ViterbiAlgorithm va = new ViterbiAlgorithm(obs, tp, ep, initialProb);
            va.calculate();
            System.out.printf("Hets\t %d:%s \t %d \t %d %n", bt,name, sum, hsum);
            byte[] states = va.getMostProbableStateSequence();
            System.out.println(Arrays.toString(obs));
            System.out.println(Arrays.toString(states));
            va = new ViterbiAlgorithm(obsRev, tp, ep, initialProb);
            va.calculate();
            byte[] statesRev = va.getMostProbableStateSequence();
            for (int b = 0; b < blocks; b++) {
                states[b]=statesRev[blocks-b-1];
            }
            System.out.println(Arrays.toString(states));
            for (int focusBlock = 0; focusBlock < blocks; focusBlock++) {
                blockHetProb[bt][focusBlock]=Float.NaN;
//                int[] resultRange=getBlockWithMinMinorCount(mnT, focusBlock, minMinorCnt);
//                total++;
//                regionHypth[focusBlock]=getBestInbredDonors(bt, resultRange[0],resultRange[2], focusBlock);

            }
            
        }
        System.out.printf("Totals Blocks\t %d  HetBlocks\t %d \t %g %n", totalBlocks, 
                totalhetBlocks, (double)totalhetBlocks/(double)totalBlocks);
        
    }
    
    private double[][] simulateHetInformation(double maf, boolean realMAF) {
        int maxDepth=65;
        double[][] hetPerCoverage=new double[65][65];
        int[][] hetCoverageCnt=new int[65][65];
        int[] coverageCnt=new int[65];
        Random r=new Random();
        int mafcnt=0;
        for(int x=0; x<100000; x++) {
            byte[] mjA=new byte[64];
            byte[] mnA=new byte[64];
            for (int c = 0; c < x%256; c++) {
                int site=r.nextInt(64);
                if(realMAF) maf=tA.getMinorAlleleFrequency((mafcnt++)%tA.getSiteCount());
                if(r.nextDouble()<maf) {mnA[site]++;}
                else {mjA[site]++;}
            }
            int cnt=0, hets=0;
            for (int i = 0; i < mnA.length; i++) {
                int total=mjA[i]+mnA[i];
                if(total>0) {
                    cnt++;
                    if((mjA[i]>0)&&(mnA[i]>0)) hets++;
                }
            }
            coverageCnt[cnt]++;
            hetCoverageCnt[cnt][hets]++;
           // if(hets>1) hetCoverageCnt[1][cnt]++;
        }
        for (int i = 0; i < hetPerCoverage.length; i++) {
            for (int j = 0; j < maxDepth; j++) {
                hetPerCoverage[i][j]=(double)hetCoverageCnt[i][j]/(double)coverageCnt[i]+0.000001;     
            }
//           hetPerCoverage[0][i]=(double)hetCoverageCnt[0][i]/(double)coverageCnt[i]; 
//           hetPerCoverage[1][i]=(double)hetCoverageCnt[1][i]/(double)coverageCnt[i]; 
        }
//        for (int j = 0; j < maxDepth; j++) {
//                System.out.println("Cvg:"+j+" simulateHetInformation:"+Arrays.toString(hetPerCoverage[j]));   
//            }
//        System.out.println("simulateHetInformation:"+Arrays.toString(hetPerCoverage[0]));
//        System.out.println("simulateHetInformation:"+Arrays.toString(hetPerCoverage[1]));
        return hetPerCoverage;
    }
    
    private EmissionProbability setEmmisionProb(double[] hetExp, double[] homoExp) {
        EmissionProbability ep = new EmissionProbability();
        double hetGivenHet=0.20;
        double coverageCorrector=1;
        //states rows, observations columns
        double[][] probMatrix = new double[2][2];

//        // homozygous state
//        probMatrix[0][0] = .95 * (coverageCorrector); //observe hom
//        probMatrix[0][1] = .05 * (coverageCorrector);//observe het
//        //heterozygous state
//        probMatrix[1][0] = (1 - hetGivenHet) * (coverageCorrector); //observe hom
//        probMatrix[1][1] = hetGivenHet * (coverageCorrector);//observe het

        
        // homozygous state
//        probMatrix[0][0] = 1-homoExp; //observe hom
//        probMatrix[0][1] = homoExp;//observe het
        probMatrix[0]=homoExp;
        //heterozygous state
//        probMatrix[1][0] = 1-hetExp; //observe hom
//        probMatrix[1][1] = hetExp;//observe het
        probMatrix[1]=hetExp;
        //System.out.printf("hetExp:%g homoExp:%g %n",hetExp, homoExp);
//        System.out.println(Arrays.deepToString(probMatrix));
        ep.setEmissionProbability(probMatrix);
        return ep;
    }
    
    private TransitionProbability setTransitionProb(Alignment a, int blocks) {
        TransitionProbability tp = new TransitionProbability();
    	int nsites = blocks;
    	int ntaxa = a.getSequenceCount();
    	int chrlen = blocks;
    	
    	double phet = 0.50;
    	int totalTransitions = (nsites - 1) * ntaxa /1;//this is the key parameter that needs to be thought through
    	int hetHet = (int) Math.floor(phet*totalTransitions);
    	int hetHom = 2 * ntaxa;
    	//int[][] transCount = new int[][]{{totalTransitions - hetHet - 2*hetHom, hetHom},{hetHom, hetHet}};
        int[][] transCount = new int[][]{{30000, 500},{500, 30000}};
    	
    	tp.setTransitionCounts(transCount, chrlen, ntaxa);
        int[] pos=new int[nsites];
        for (int i = 0; i < pos.length; i++) {
            pos[i]=i;
        }
    	tp.setPositions(pos);
        
        return tp;
    }
    
    private int[][] profileTaxaForHets(long[] mjT, long[] mnT) {
        int[][] result=new int[2][mjT.length];
        for (int i = 0; i < mjT.length; i++) {
            result[0][i]=Long.bitCount(mjT[i]|mnT[i]);
            result[1][i]=Long.bitCount(mjT[i]&mnT[i]);
        }
        return result;
    }
    
    
    
    /**
     *
     * @param args
     */
    public static void main(String[] args) {
      
      String root="/Users/edbuckler/SolexaAnal/GBS/build20120110/imp/";
//        String root="/Volumes/LaCie/build20120110/imp/";

        String donorFile=root+"NAMfounder20120110.imp.hmp.txt";
 //       String donorFile=root+"DTMAfounder20120110.imp.hmp.txt";
        String unImpTargetFile=root+"ZeaSyn20120110.hmp.txt";
        String impTargetFile=root+"ZeaSyn20120110.imp.hmp.txt";

        boolean buildInput=true;
        boolean filterTrue=true;
        if(buildInput) {KnownParentMinorWindowImputation.createSynthetic(donorFile, unImpTargetFile, 2000, 0.4, 0.5, 1000);}

        HeterozygosityProfiler hp=new HeterozygosityProfiler(unImpTargetFile, impTargetFile);
        hp.simulateHetInformation(-1, true);
        hp.simulateHetInformation(0.01, false);
        hp.simulateHetInformation(0.002, false);
    //    HeterozygosityProfiler hp=new HeterozygosityProfiler(donorFile, impTargetFile);

    }
}
