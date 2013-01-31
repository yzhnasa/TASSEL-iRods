/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.pal.popgen;

import java.util.ArrayList;
import java.util.Random;
import java.util.TreeMap;
import net.maizegenetics.pal.alignment.*;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.BitUtil;
import net.maizegenetics.util.OpenBitSet;
import net.maizegenetics.util.ProgressListener;


/**
 *  Imputation methods that relies on a list of possible parents, and optimized for
 * highly heterozygous samples.  It works at the scale of 64 sites to accelerate searching
 * through all possible parental combinations for a window.  
 * 
 * It only begins each 64 site block and expand outward to find set number of minor
 * alleles in the target sequence, and then it looks for all possible parents.
 * 
 * @author edbuckler, kswarts, aromero
 */
public class MinorWindowImputation {
    private Alignment unimpAlign;
    private Alignment donorAlign;
    private int testing=0;  //level of reporting to stdout
    private OpenBitSet swapMjMnMask, invariantMask, errorMask, goodMask;  //mask of sites that major and minor are swapped
    private int parentsRight=0, parentsWrong=0;
    private int blocks=-1;
    int[] hSite, hTaxon;
    byte[] hState;
    int maskSitCnt=0;
    boolean maskAndTest=true;
    boolean isSwapMajorMinor=false;
    private int maxDonorHypotheses=20;
    private double maximumInbredError=0.02;  //inbreds are tested first, if too much error hybrids are tested.
    private int minTestSites=20;  //minimum number of compared sites to find a hit
    private byte[][][] allDist;
    
    private static int highMask=0xF0;
    private static int lowMask=0x0F;
    
    int totalRight=0, totalWrong=0;

    /**
     * 
     * @param donorFile should be imputed inbreds that were founders of population
     * @param unImpTargetFile sites must match exactly with donor file
     * @param exportFile 
     * @param minMinorCnt determines the size of the search window, low recombination 20-30, high recombination 10-15
     * @param resolveMethod 0=use top ten; 1=uses only the best
     */
    public MinorWindowImputation(String donorFile, String unImpTargetFile, 
            String exportFile, int minMinorCnt, int resolveMethod) {
        donorAlign=ImportUtils.readFromHapmap(donorFile, false, (ProgressListener)null);
        donorAlign.optimizeForTaxa(null);
        System.out.printf("Donor taxa:%d sites:%d %n",donorAlign.getSequenceCount(),donorAlign.getSiteCount());        
        unimpAlign=ImportUtils.readFromHapmap(unImpTargetFile, false, (ProgressListener)null);
        unimpAlign.optimizeForTaxa(null);
        System.out.printf("Donor taxa:%d sites:%d %n",unimpAlign.getSequenceCount(),unimpAlign.getSiteCount());
        
        if(maskAndTest) maskSites(200);  //mask sites to test imputation accuracy
        blocks=unimpAlign.getAllelePresenceForAllSites(0, 0).getNumWords();
        
        //Create mask for all sites where major & minor are swapped in alignments
        goodMask=new OpenBitSet(unimpAlign.getSiteCount());
        errorMask=new OpenBitSet(unimpAlign.getSiteCount());
        swapMjMnMask=new OpenBitSet(unimpAlign.getSiteCount());
        invariantMask=new OpenBitSet(unimpAlign.getSiteCount());
        int siteConflicts=0, swaps=0, invariant=0, good=0;
        for (int i = 0; i < unimpAlign.getSiteCount(); i++) {
            /*we have three classes of data:  invariant in one alignment, conflicts about minor and minor,
            *swaps of major and minor.  Adding the invariant reduces imputation accuracy.
            *the major/minor swaps should be flipped in the comparisons
            *
            */
            if(donorAlign.getMinorAllele(i)==Alignment.UNKNOWN_ALLELE) {
                invariant++;
                invariantMask.set(i);
                goodMask.set(i);
            } else 
            if((donorAlign.getMajorAllele(i)==unimpAlign.getMinorAllele(i))&&(donorAlign.getMinorAllele(i)==unimpAlign.getMajorAllele(i))) {
                swaps++;
                swapMjMnMask.set(i);
                goodMask.set(i);
            } else
            if((donorAlign.getMajorAllele(i)!=unimpAlign.getMajorAllele(i))) {
                siteConflicts++;
                errorMask.set(i);
                goodMask.set(i);
            } 
            
        }
        goodMask.not();
        System.out.println("invariant in donor:"+invariant+" swapConflicts"+swaps+" errors:"+siteConflicts);
        MutableNucleotideAlignment mna=MutableNucleotideAlignment.getInstance(this.unimpAlign);
        Random r=new Random(0);
        int total=0, hybrid=0;
        long time=System.currentTimeMillis();
        for (int bt = 0; bt < unimpAlign.getSequenceCount(); bt+=1) {
            System.out.println("");
            DonorHypoth[][] regionHypth=new DonorHypoth[blocks][maxDonorHypotheses];
            String name=unimpAlign.getIdGroup().getIdentifier(bt).getFullName();
            BitSet[] modBits=arrangeMajorMinorBtwAlignments(unimpAlign,bt);
            System.out.printf("Imputing %d:%s Mj:%d, Mn:%d Unk:%d ... ", bt,name,modBits[0].cardinality(),
                    modBits[1].cardinality(), countUnknown(mna,bt));
 //           System.out.println(mna.getBaseAsStringRow(bt));
            if(modBits[0].cardinality()<1000) continue;
            long[] mnT=modBits[1].getBits();
            calcInbredDist(modBits);
            for (int focusBlock = 0; focusBlock < blocks; focusBlock++) {
                int[] resultRange=getBlockWithMinMinorCount(mnT, focusBlock, minMinorCnt);
                total++;
                regionHypth[focusBlock]=getBestInbredDonors(bt, resultRange[0],resultRange[2], focusBlock, modBits);
//                for (int i = 0; i < regionHypth[focusBlock].length; i++) {
//                    if(regionHypth[focusBlock][i]==null) continue;
//                    regionHypth[focusBlock][i].donor1Taxon=r.nextInt(donorAlign.getSequenceCount());
//                    regionHypth[focusBlock][i].donor2Taxon=regionHypth[focusBlock][i].donor1Taxon;
//                }
            }
            solveRegionally(mna, bt, regionHypth,resolveMethod);
            int unk=countUnknown(mna,bt);
            System.out.printf("Done Unk: %d Prop:%g", unk, (double)unk/(double)mna.getSiteCount());
          //  System.out.println(bt+mna.getBaseAsStringRow(bt));
        }
        System.out.println("");
        System.out.println("Time:"+(time-System.currentTimeMillis()));
        StringBuilder s=new StringBuilder();
        s.append(String.format("%s %s MinMinor:%d ", donorFile, unImpTargetFile, minMinorCnt));
        if(maskAndTest) s.append(compareSites(mna));
        if(testing>0) s.append(String.format("ParentsRight:%d Wrong %d", parentsRight, parentsWrong));
        System.out.println(s.toString());
        System.out.printf("TotalRight %d  TotalWrong %d Rate%n",totalRight, totalWrong);
        
        ExportUtils.writeToHapmap(mna, false, exportFile, '\t', null);
    }
    
    private int countUnknown(Alignment a, int taxon) {
        int cnt=0;
        for (int i = 0; i < a.getSiteCount(); i++) {
            if(a.getBase(taxon, i)==Alignment.UNKNOWN_DIPLOID_ALLELE) cnt++;
        }
        return cnt;
    }
    
    private BitSet[] arrangeMajorMinorBtwAlignments(Alignment unimpAlign, int bt) { 
        OpenBitSet mjTbs=new OpenBitSet(unimpAlign.getAllelePresenceForAllSites(bt, 0));
        OpenBitSet mnTbs=new OpenBitSet(unimpAlign.getAllelePresenceForAllSites(bt, 1));
        mjTbs.and(goodMask);
        mnTbs.and(goodMask);
        if(testing>3) System.out.printf("PostGood %d mjTbs %d mnTbs %d %n",mjTbs.size(),mjTbs.cardinality(),mnTbs.cardinality());
        OpenBitSet newmj=new OpenBitSet(unimpAlign.getAllelePresenceForAllSites(bt, 1));
        OpenBitSet newmn=new OpenBitSet(unimpAlign.getAllelePresenceForAllSites(bt, 0));
        if(testing>3) System.out.printf("%d newmj %d newmn %d %n",newmj.size(),newmj.cardinality(),newmn.cardinality());
        newmj.and(swapMjMnMask);
        newmn.and(swapMjMnMask);
        if(testing>3) System.out.printf("newmj %d newmn %d %n",newmj.cardinality(),newmn.cardinality());
        if(isSwapMajorMinor) {
            mjTbs.or(newmj);
            mnTbs.or(newmn);
        }
        if(testing>3) System.out.printf("%d mjTbs %d mnTbs %d %n",newmj.size(),mjTbs.cardinality(),mnTbs.cardinality());
        BitSet[] result={mjTbs,mnTbs};
        return result;
    }
    
    /**
     * If the target regions has Mendelian errors that it looks for overlapping regional
     * solutions that are better.
     * @param mna
     * @param targetTaxon
     * @param regionHypth 
     */
    private void solveRegionally(MutableNucleotideAlignment mna, int targetTaxon, 
            DonorHypoth[][] regionHypth, int resolveMethod) {
 //       System.out.printf("Start Tx:%d F:Before Unk:%d %n",targetTaxon, countUnknown(mna, targetTaxon));
        for (int focusBlock = 0; focusBlock < blocks; focusBlock++) {
            DonorHypoth cbh=regionHypth[focusBlock][0];
//            System.out.printf("%d %d %d %g %n", regionHypth[focusBlock][0].targetTaxon, 
//                    regionHypth[focusBlock][0].focusBlock, regionHypth[focusBlock][0].mendelianErrors, 
//                    regionHypth[focusBlock][0].getErrorRate());
            if(regionHypth[focusBlock][0].mendelianErrors==0) {
                //setAlignmentWithDonors(cbh,focusBlock,mna);
                if(resolveMethod==0) {setAlignmentWithDonors(regionHypth[focusBlock],focusBlock,mna);}
                if(resolveMethod==1) {setAlignmentWithDonors(regionHypth[focusBlock][0],focusBlock,mna);}
            }
            else {
                int minMendelErrors=cbh.mendelianErrors;
                int currBestBlock=cbh.focusBlock;
                int bestBlockDistance=Integer.MAX_VALUE;
                for (int i = cbh.startBlock; i <= cbh.endBlock; i++) {
                    if((regionHypth[i][0].startBlock<=focusBlock)&&(regionHypth[i][0].endBlock>=focusBlock)) {
                        if(regionHypth[i][0].mendelianErrors<minMendelErrors) {
                            currBestBlock=i;
                            minMendelErrors=regionHypth[i][0].mendelianErrors;
                            bestBlockDistance=Math.abs(i-focusBlock);
                        } 
                        else if((regionHypth[i][0].mendelianErrors==minMendelErrors)&&
                                (Math.abs(i-focusBlock)<bestBlockDistance)) {
                            currBestBlock=i;
                            minMendelErrors=regionHypth[i][0].mendelianErrors;
                            bestBlockDistance=Math.abs(i-focusBlock);
                        }
                    }
                    
                }
//                System.out.printf("c: %d %d %d %g %n", regionHypth[currBestBlock][0].targetTaxon, 
//                    regionHypth[currBestBlock][0].focusBlock, regionHypth[currBestBlock][0].mendelianErrors,
//                    regionHypth[focusBlock][0].getErrorRate());
                if(regionHypth[currBestBlock][0].getErrorRate()<this.maximumInbredError) {
                    if(resolveMethod==0) {setAlignmentWithDonors(regionHypth[currBestBlock],focusBlock,mna);}
                    if(resolveMethod==1) {setAlignmentWithDonors(regionHypth[currBestBlock][0],focusBlock,mna);}
                }
            }
    //        System.out.printf("Tx:%d F:%d Unk:%d %n",targetTaxon,focusBlock, countUnknown(mna, targetTaxon));
        }
    }
    
    /**
     * Given a start 64 site block, it expands to the left and right until it hits
     * the minimum Minor Site count in the target taxon
     * @param mnT - minor allele bit presence in a series of longs
     * @param focusBlock
     * @param minMinorCnt
     * @return arrays of blocks {startBlock, focusBlock, endBlock}
     */
    private int[] getBlockWithMinMinorCount(long[] mnT, int focusBlock, int minMinorCnt) {
        int minorCnt=Long.bitCount(mnT[focusBlock]);
        int endBlock=focusBlock, startBlock=focusBlock;
        while(minorCnt<minMinorCnt) {
            boolean preferMoveStart=(focusBlock-startBlock<endBlock-focusBlock)?true:false;
          //  System.out.println("bool startBlock"+startBlock+" focus"+focusBlock);
            if(startBlock==0) {preferMoveStart=false;}
            if(endBlock==blocks-1) {preferMoveStart=true;}
            if((startBlock==0)&&(endBlock==blocks-1)) break;
            if(preferMoveStart) {//expand start
                startBlock--;
              //  System.out.println("startBlock"+startBlock+" minorCnt"+minorCnt);
                minorCnt+=Long.bitCount(mnT[startBlock]);
            } else { //expand end
                endBlock++;
               // System.out.println("endBlock"+endBlock+" minorCnt"+minorCnt);
                minorCnt+=Long.bitCount(mnT[endBlock]);  
            } 
        }
        int[] result={startBlock, focusBlock, endBlock};
        return result;
    }
    
        /**
     * 
     * @param donor1
     * @param modBitsOfTarget
     * @return array with testSites, minorSites, unmatched (errors)
     */
    private void calcInbredDist(BitSet[] modBitsOfTarget) {
        allDist=new byte[donorAlign.getSequenceCount()][3][blocks];
    //    byte[][] result=new byte[3][blocks];
        long[] mjT=modBitsOfTarget[0].getBits();
        long[] mnT=modBitsOfTarget[1].getBits();
        for (int donor1 = 0; donor1 < allDist.length; donor1++) {
            long[] mj1=donorAlign.getAllelePresenceForAllSites(donor1, 0).getBits();
            long[] mn1=donorAlign.getAllelePresenceForAllSites(donor1, 1).getBits();
            for (int i = 0; i <blocks; i++) {
                long siteMask=(mjT[i]|mnT[i])&(mj1[i]|mn1[i]);
                allDist[donor1][2][i]+=(byte)Long.bitCount(siteMask&mjT[i]&(mjT[i]^mj1[i]));
                allDist[donor1][2][i]+=(byte)Long.bitCount(siteMask&mnT[i]&(mnT[i]^mn1[i]));
                allDist[donor1][0][i]=(byte)Long.bitCount(siteMask);
                allDist[donor1][1][i]=(byte)Long.bitCount(siteMask&mnT[i]);  
            }
        }
      //  return result;
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
    private DonorHypoth[] getBestInbredDonors(int targetTaxon, int startBlock, int endBlock, int focusBlock,
            BitSet[] modBits) {
//       long[] mjT=modBits[0].getBits(startBlock, endBlock);
//       long[] mnT=modBits[1].getBits(startBlock, endBlock);
        int[] donors={-1,-1,0};
        double minPropUnmatched=1.0;
        int maxTestSites=0;
        int donorTieCnt=0;
        //int[] rDonors=parseDonorsNamesForRegion(unimpAlign.getTaxaName(targetTaxon),startBlock*64);
//        if(testing>3) System.out.printf("StartSite %d EndSite %d RealD1 %d RealD2 %d %n",startBlock*64, 
//                (endBlock*64+63),rDonors[0],rDonors[1]);

        TreeMap<Double,DonorHypoth> bestDonors=new TreeMap<Double,DonorHypoth>();
        bestDonors.put(1.0, new DonorHypoth());
        double lastKeytestPropUnmatched=1.0;
        int donorTaxaCnt=donorAlign.getSequenceCount();
        for (int d1 = 0; d1 < donorTaxaCnt; d1++) {
            int mjUnmatched=0;
            int mnUnmatched=0;
            int testSites=0;
            for (int i = startBlock; i <=endBlock; i++) {
                mjUnmatched+=allDist[d1][2][i];
                testSites+=allDist[d1][0][i];
                if(mjUnmatched>5) break;
            }
            if(testSites<minTestSites) continue;
            int totalMendelianErrors=mjUnmatched+mnUnmatched;
            double testPropUnmatched=(double)(totalMendelianErrors)/(double)testSites;
            if((bestDonors.size()<maxDonorHypotheses)||(testPropUnmatched<lastKeytestPropUnmatched)) {
                DonorHypoth theDH=new DonorHypoth(targetTaxon, d1, d1, startBlock,
                    focusBlock, endBlock, testSites, totalMendelianErrors);
                bestDonors.put(new Double(testPropUnmatched), theDH);
                if(bestDonors.size()>maxDonorHypotheses) {
                    bestDonors.remove(bestDonors.lastKey());
                    lastKeytestPropUnmatched=bestDonors.lastKey();
                }
            }
//            if((testing>1)&&(rDonors[0]==d1)&&(rDonors[1]==d1))
//                System.out.printf("Donor %d %d %d %d %d block %d-%d-%d %n", d1, d1, mjUnmatched, mnUnmatched,
//                        testSites, startBlock, focusBlock, endBlock);

            if((testPropUnmatched<minPropUnmatched)||
                    ((testPropUnmatched==minPropUnmatched)&&(testSites>maxTestSites))) {
                donors[0]=d1;
                donors[1]=d1;
                minPropUnmatched=testPropUnmatched;
                donors[2]=maxTestSites=testSites;
                donorTieCnt=0;
            } else if(testPropUnmatched==minPropUnmatched) donorTieCnt++;
            
        }
        DonorHypoth[] result=new DonorHypoth[maxDonorHypotheses];
        int count=0;
        for (DonorHypoth dh : bestDonors.values()) {
            result[count]=dh; 
            count++;
        }
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
     * @return int[] array of {donor1, donor2, testSites}
     */
    private DonorHypoth[] getBestDonors(int targetTaxon, int startBlock, int endBlock, int focusBlock) {
        long[] mjT=unimpAlign.getAllelePresenceForAllSites(targetTaxon, 0).getBits();
        long[] mnT=unimpAlign.getAllelePresenceForAllSites(targetTaxon, 1).getBits();
        int[] donors={-1,-1,0};
        double minPropUnmatched=1.0;
        int maxTestSites=0;
        int donorTieCnt=0;
        int[] rDonors=parseDonorsNamesForRegion(unimpAlign.getTaxaName(targetTaxon),startBlock*64);
        if(testing>3) System.out.printf("StartSite %d EndSite %d RealD1 %d RealD2 %d %n",startBlock*64, 
                (endBlock*64+63),rDonors[0],rDonors[1]);
        long[] swapMask=this.swapMjMnMask.getBits();
        TreeMap<Double,DonorHypoth> bestDonors=new TreeMap<Double,DonorHypoth>();
        bestDonors.put(1.0, new DonorHypoth());
        
        for (int d1 = 0; d1 < donorAlign.getSequenceCount(); d1++) {
            long[] mj1=donorAlign.getAllelePresenceForAllSites(d1, 0).getBits();
            long[] mn1=donorAlign.getAllelePresenceForAllSites(d1, 1).getBits();
            for (int d2 = d1; d2 < donorAlign.getSequenceCount(); d2++) {
                long[] mj2=donorAlign.getAllelePresenceForAllSites(d2, 0).getBits();
                long[] mn2=donorAlign.getAllelePresenceForAllSites(d2, 1).getBits();
                int mjUnmatched=0;
                int mnUnmatched=0;
                int testSites=0;
                int testTargetMajor=0;
                int testTargetMinor=0;
                for (int i = startBlock; i <= endBlock; i++) {
                    long siteMask=swapMask[i]&(mjT[i]|mnT[i])&(mj1[i]|mn1[i])&(mj2[i]|mn2[i]);
                    mjUnmatched+=Long.bitCount(siteMask&mjT[i]&(mjT[i]^mj1[i])&(mjT[i]^mj2[i]));
                    mnUnmatched+=Long.bitCount(siteMask&mnT[i]&(mnT[i]^mn1[i])&(mnT[i]^mn2[i]));
                    testSites+=Long.bitCount(siteMask);
                    testTargetMajor+=Long.bitCount(siteMask&mjT[i]);
                    testTargetMinor+=Long.bitCount(siteMask&mnT[i]);
                }
                int totalMendelianErrors=mjUnmatched+mnUnmatched;
                double testPropUnmatched=(double)(totalMendelianErrors)/(double)testSites;
                if(testPropUnmatched<bestDonors.lastKey()) {
                    DonorHypoth theDH=new DonorHypoth(targetTaxon, d1, d2, startBlock,
                        focusBlock, endBlock, testSites, totalMendelianErrors);
                    bestDonors.put(new Double(testPropUnmatched), theDH);
                    if(bestDonors.size()>maxDonorHypotheses) bestDonors.remove(bestDonors.lastKey());
                }
                if((testing>1)&&(rDonors[0]==d1)&&(rDonors[1]==d2))
                    System.out.printf("Donor %d %d %d %d %d %d %d block %d-%d-%d %n", d1, d2, mjUnmatched, mnUnmatched,
                            testSites, testTargetMajor, testTargetMinor, startBlock, focusBlock, endBlock);
                
                if((testPropUnmatched<minPropUnmatched)||
                        ((testPropUnmatched==minPropUnmatched)&&(testSites>maxTestSites))) {
                    donors[0]=d1;
                    donors[1]=d2;
                    minPropUnmatched=testPropUnmatched;
                    donors[2]=maxTestSites=testSites;
                    donorTieCnt=0;
                } else if(testPropUnmatched==minPropUnmatched) donorTieCnt++;
                
            }
            
        }
        if(testing>1){
            if((rDonors[0]==donors[0])&&(rDonors[1]==donors[1])) {
                System.out.printf("Correct ties:%d bestMatch:%g %n",donorTieCnt,minPropUnmatched);}
            else {System.out.printf("WRONG D1:%d D2:%d ties:%d bestMatch:%g %n",donors[0], donors[1],donorTieCnt,minPropUnmatched);}
        }
        if((rDonors[0]==donors[0])&&(rDonors[1]==donors[1])) {parentsRight++;}
            else {parentsWrong++;}
        DonorHypoth[] result=new DonorHypoth[maxDonorHypotheses];
        int count=0;
        for (DonorHypoth dh : bestDonors.values()) {
            result[count]=dh; 
            count++;
        }
        return result;
    }
  
    
    /**
     * Takes a donor hypothesis and applies it to the output alignment 
     * @param theDH
     * @param mna 
     */
    private void setAlignmentWithDonors(DonorHypoth theDH, int focusBlock, MutableNucleotideAlignment mna) {
        if(theDH.donor1Taxon<0) return;
        int startSite=focusBlock*64;
        int endSite=startSite+63;
        if(endSite>=unimpAlign.getSiteCount()) endSite=unimpAlign.getSiteCount()-1;
        for(int cs=startSite; cs<=endSite; cs++) {
            byte bD1=donorAlign.getBase(theDH.donor1Taxon, cs);
            byte bD2=donorAlign.getBase(theDH.donor2Taxon, cs);
            byte donorEst=Alignment.UNKNOWN_DIPLOID_ALLELE;
            if(bD1==Alignment.UNKNOWN_DIPLOID_ALLELE) {donorEst=bD2;}
            else if(bD2==Alignment.UNKNOWN_DIPLOID_ALLELE) {donorEst=bD1;}
            else {donorEst=(byte)((bD1&highMask)|(bD2&lowMask));
            }
            //need to check whether the heterozygote is put together properly
            //need to change to 
            byte knownBase=mna.getBase(theDH.targetTaxon, cs);
            if((knownBase!=Alignment.UNKNOWN_DIPLOID_ALLELE)&&(knownBase!=donorEst)&&(donorEst!=Alignment.UNKNOWN_DIPLOID_ALLELE)) {
//                System.out.printf("Error %d %s %s %n",theDH.targetTaxon, mna.getBaseAsString(theDH.targetTaxon, cs),
//                        NucleotideAlignmentConstants.getNucleotideIUPAC(donorEst));      
                if(!AlignmentUtils.isHeterozygous(donorEst)) totalWrong++;
            }
            if((knownBase!=Alignment.UNKNOWN_DIPLOID_ALLELE)&&(knownBase==donorEst)&&(donorEst!=Alignment.UNKNOWN_DIPLOID_ALLELE)) {
//                System.out.printf("Error %d %s %s %n",theDH.targetTaxon, mna.getBaseAsString(theDH.targetTaxon, cs),
//                        NucleotideAlignmentConstants.getNucleotideIUPAC(donorEst));
                totalRight++;
            }
            mna.setBase(theDH.targetTaxon, cs, donorEst);
//            if(mna.getBase(theDH.targetTaxon, cs)==Alignment.UNKNOWN_DIPLOID_ALLELE) {
//                    mna.setBase(theDH.targetTaxon, cs, donorEst);
////                    impSiteCnt++;
////                    taxaImpCnt++;
//                }
        } //end of cs loop

    }
    
        /**
     * Takes a donor hypothesis and applies it to the output alignment 
     * @param theDH
     * @param mna 
     */
    private void setAlignmentWithDonors(DonorHypoth[] theDH, int focusBlock, MutableNucleotideAlignment mna) {
        if(theDH[0].targetTaxon<0) return;
        int startSite=focusBlock*64;
        int endSite=startSite+63;
        if(endSite>=unimpAlign.getSiteCount()) endSite=unimpAlign.getSiteCount()-1;
//        System.out.println("B:"+mna.getBaseAsStringRange(theDH[0].targetTaxon, startSite, endSite));
        for(int cs=startSite; cs<=endSite; cs++) {
            byte donorEst=Alignment.UNKNOWN_DIPLOID_ALLELE; 
            for (int i = 0; i < theDH.length; i++) {
                if((theDH[i]==null)||(theDH[i].donor1Taxon<0)) continue;
                byte bD1=donorAlign.getBase(theDH[i].donor1Taxon, cs);
                byte bD2=donorAlign.getBase(theDH[i].donor2Taxon, cs);
                donorEst=Alignment.UNKNOWN_DIPLOID_ALLELE;
                if(bD1==Alignment.UNKNOWN_DIPLOID_ALLELE) {donorEst=bD2;}
                else if(bD2==Alignment.UNKNOWN_DIPLOID_ALLELE) {donorEst=bD1;}
                else {donorEst=(byte)((bD1&highMask)|(bD2&lowMask));}
                if(donorEst!=Alignment.UNKNOWN_DIPLOID_ALLELE) break;
            }
            byte knownBase=mna.getBase(theDH[0].targetTaxon, cs);
            if((knownBase!=Alignment.UNKNOWN_DIPLOID_ALLELE)&&(knownBase!=donorEst)&&(donorEst!=Alignment.UNKNOWN_DIPLOID_ALLELE)) {
//                System.out.printf("Error %d %s %s %n",theDH.targetTaxon, mna.getBaseAsString(theDH.targetTaxon, cs),
//                        NucleotideAlignmentConstants.getNucleotideIUPAC(donorEst));      
                if(!AlignmentUtils.isHeterozygous(donorEst)) totalWrong++;
            }
            if((knownBase!=Alignment.UNKNOWN_DIPLOID_ALLELE)&&(knownBase==donorEst)&&(donorEst!=Alignment.UNKNOWN_DIPLOID_ALLELE)) {
//                System.out.printf("Error %d %s %s %n",theDH.targetTaxon, mna.getBaseAsString(theDH.targetTaxon, cs),
//                        NucleotideAlignmentConstants.getNucleotideIUPAC(donorEst));
                totalRight++;
            }
            mna.setBase(theDH[0].targetTaxon, cs, donorEst);
        } //end of cs loop
//        System.out.println("E:"+mna.getBaseAsStringRange(theDH[0].targetTaxon, startSite, endSite));
    }
    
    private int[] parseDonorsNamesForRegion(String taxonName, int site) {
        String[] sb=taxonName.split("\\|");
        int[] donors={-1,-1};
        for (int i = 1; i < sb.length; i++) {
            String[] ss=sb[i].split("s");
            int leftSite=Integer.parseInt(ss[1]);
            if(leftSite<=site) {
                String[] sd=ss[0].split("_");
                donors[0]=Integer.parseInt(sd[0]);
                donors[1]=Integer.parseInt(sd[1]);
              //  break;
            }
        }
        return donors;
    }
    
    /**
     * Used for testing imputation accuracy.  Sets every sampIntensity sites to missing
     * so they can be compared at the end of imputation.
     * @param sampIntensity 
     */
    private void maskSites(int sampIntensity) {
        System.out.println("Beginning to mask sites");
        MutableNucleotideAlignment mna=MutableNucleotideAlignment.getInstance(unimpAlign);
        int maxSites=mna.getSequenceCount()*((mna.getSiteCount()/sampIntensity)+1);
        hSite=new int[maxSites];
        hTaxon=new int[maxSites];
        hState=new byte[maxSites];
        int cnt=0;
        for (int t = 0; t < mna.getSequenceCount(); t++) {
            for (int s = t%sampIntensity; s < mna.getSiteCount(); s+=sampIntensity) {
                hSite[cnt]=s;
                hTaxon[cnt]=t;
                hState[cnt]=mna.getBase(t, s);
                mna.setBase(t, s, Alignment.UNKNOWN_DIPLOID_ALLELE);
                cnt++;
            }
 //           System.out.println(t+":"+cnt);
        }
        maskSitCnt=cnt;
        mna.clean();
        compareSites(unimpAlign);
        unimpAlign=BitAlignment.getInstance(mna, false);
        compareSites(unimpAlign);
        System.out.println("Sites masked");
    }
    
    private String compareSites(Alignment a) {
        int missingCnt=0, correctCnt=0, errorCnt=0, notImp=0, hetCnt=0;
        for (int i = 0; i < maskSitCnt; i++) {
            if(hState[i]==Alignment.UNKNOWN_DIPLOID_ALLELE) {
                missingCnt++;
                continue;
            }
            byte impb=a.getBase(hTaxon[i], hSite[i]);
            if(AlignmentUtils.isHeterozygous(impb)) {
                hetCnt++;
                byte[] orig=AlignmentUtils.getDiploidValues(hState[i]);
                byte[] imphet=AlignmentUtils.getDiploidValues(impb);
                if((orig[0]==imphet[0])||(orig[0]==imphet[1])||(impb==hState[i])) {correctCnt++;}
                else {errorCnt++;}
            } else if(impb==Alignment.UNKNOWN_DIPLOID_ALLELE) {
                notImp++;
            } else if(impb==hState[i]) {
                correctCnt++;
            } else {errorCnt++;}
        }
        double errRate=(double)errorCnt/(double)(errorCnt+correctCnt);
        return String.format("Missing: %d Het: %d NotImp: %d Error: %d Correct: %d ErrorRate: %g ", 
                missingCnt, hetCnt, notImp, errorCnt, correctCnt, errRate);
        
    }

   
    
      
    
    /**
     *
     * @param args
     */
    public static void main(String[] args) {
//      String root="/Users/edbuckler/SolexaAnal/GBS/build20120110/imp/";
//        String root="/Volumes/LaCie/build20120110/imp/";
        String root="/Users/edbuckler/SolexaAnal/GBS/build20120701/06_HapMap/";

        String donorFile=root+"TaxaWMRG8FINAL_chr10.hmp.txt.gz";
 //       String donorFile=root+"AllTaxa_BPEC_AllZea_GBS_Build_July_2012_FINAL_chr10.hmp.txt.gz";
 //       String donorFile=root+"DTMAfounder20120110.imp.hmp.txt";
        String unImpTargetFile=root+"TaxaW142FINAL_chr10.hmp.txt.gz";
        String impTargetFile=root+"Test20120110.imp.hmp.txt";

       // boolean buildInput=false;
       // boolean filterTrue=true;
       // if(buildInput) {KnownParentMinorWindowImputation.createSynthetic(donorFile, unImpTargetFile, 2000, 0.4, 0.5, 1000);}

        System.out.println("Resolve Method 0");
        MinorWindowImputation e64NNI;
        e64NNI=new MinorWindowImputation(donorFile,
                unImpTargetFile, impTargetFile,15,0);
        System.out.println("Resolve Method 1");
         e64NNI=new MinorWindowImputation(donorFile,
                unImpTargetFile, impTargetFile,15,1);


        
//        for (int recSize = 512; recSize < 4000; recSize+=(recSize/2)) {
//            for (int mm = 5; mm < 40; mm+=5) {
//                System.out.println("Rec size"+recSize);
//                unImpTargetFile=root+recSize+"ZeaSyn20120110.hmp.txt";
//                if(buildInput) {createSynthetic(donorFile, unImpTargetFile, recSize, 0.4, -1, 1000);}
//                e64NNI=new KnownParentMinorWindowImputation(donorFile, unImpTargetFile, impTargetFile,mm,1);
//                }
//        }
    }
    
}

