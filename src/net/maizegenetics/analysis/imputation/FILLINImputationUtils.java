package net.maizegenetics.analysis.imputation;

import net.maizegenetics.analysis.popgen.DonorHypoth;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.BitUtil;

import java.util.TreeMap;

import static net.maizegenetics.dna.snp.GenotypeTable.WHICH_ALLELE.Major;
import static net.maizegenetics.dna.snp.GenotypeTable.WHICH_ALLELE.Minor;

/**
 * Defines xxxx
 *
 * @author Ed Buckler
 */
public class FILLINImputationUtils {

    /**
     * Counts union and intersection of major and minor alleles between the target genotype and potential
     * donor genotypes.  The counts are done by 64 sites block.  These counts can be quickly used to estimate
     * distance for the set of blocks.
     * @param modBitsOfTarget major and minor presence bits for target genotype (must be aligned same as donor)
     * @param donorAlign genotypeTable with potential donor genotypes
     * @return array with [donor index][sites, same count, diff count, het count index][block index]
     */
    public static byte[][][] calcAllelePresenceCountsBtwTargetAndDonors(BitSet[] modBitsOfTarget, GenotypeTable donorAlign) {
        int blocks=modBitsOfTarget[0].getNumWords();
        byte[][][] allDist=new byte[donorAlign.numberOfTaxa()][4][blocks];
        long[] iMj=modBitsOfTarget[0].getBits();
        long[] iMn=modBitsOfTarget[1].getBits();
        for (int donor1 = 0; donor1 < allDist.length; donor1++) {
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
                allDist[donor1][2][i]=(byte)diffCnt;
                allDist[donor1][3][i]=(byte)hetCnt;
                allDist[donor1][0][i]=(byte)sites;
                allDist[donor1][1][i]=(byte)sameCnt;
            }
        }
        return allDist;
    }


    /**
     *Simple algorithm that tests every possible two donor combination to minimize
     * the number of unmatched informative alleles.  Currently, there is litte tie
     * breaking, longer matches are favored.
     * @param targetTaxon
     * @param startBlock
     * @param endBlock
     * @param focusBlock
     * @param donor1indices
     * @param targetToDonorDistances
     * @param minTestSites
     * @param maxDonorHypotheses
     * @return
     */
    public static DonorHypoth[] getBestInbredDonors(int targetTaxon, int startBlock, int endBlock, int focusBlock,
                                              int[] donor1indices, byte[][][] targetToDonorDistances,
                                              int minTestSites, int maxDonorHypotheses) {
        TreeMap<Double,DonorHypoth> bestDonors=new TreeMap<>();   //todo change this to a heap
        double lastKeytestPropUnmatched=1.0;
        double inc=1e-9;
        for (int d1 : donor1indices) {
            int testSites=0;
            int sameCnt = 0, diffCnt = 0, hetCnt = 0;
            for (int i = startBlock; i <=endBlock; i++) {
                sameCnt+=targetToDonorDistances[d1][1][i];
                diffCnt+=targetToDonorDistances[d1][2][i];
                hetCnt+=targetToDonorDistances[d1][3][i];
            }
            testSites= sameCnt + diffCnt - hetCnt;
            if(testSites<minTestSites) continue;
            double testPropUnmatched = 1.0-(((double) (sameCnt) - (double)(0.5*hetCnt)) / (double) (testSites));
            int totalMendelianErrors=(int)((double)testSites*testPropUnmatched);
            double rankingErrors=((double)totalMendelianErrors+1.0)/(double)testSites;  //todo resolve best ranking

            //prob of this given poisson with lambda=1, total Mendelian errors=k (WE SHOULD DO THIS)
            inc+=1e-9;  //this looks strange, but makes the keys unique and ordered
            testPropUnmatched+=inc;  //this looks strange, but makes the keys unique and ordered
            if(testPropUnmatched<lastKeytestPropUnmatched) {
                DonorHypoth theDH=new DonorHypoth(targetTaxon, d1, d1, startBlock,
                        focusBlock, endBlock, testSites, totalMendelianErrors);
                DonorHypoth prev=bestDonors.put(new Double(testPropUnmatched), theDH);
   //             DonorHypoth prev=bestDonors.put(rankingErrors, theDH);
                if(prev!=null) System.out.println("Inbred TreeMap index crash:"+testPropUnmatched);
                if(bestDonors.size()>maxDonorHypotheses) {
                    bestDonors.remove(bestDonors.lastKey());
                    lastKeytestPropUnmatched=bestDonors.lastKey();
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
}
