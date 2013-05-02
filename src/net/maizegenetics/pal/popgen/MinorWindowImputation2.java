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
 * Imputation methods that relies on an alignment of possible gametes (donorAlign), and optimized for
 * highly heterozygous samples.  It works at the scale of 64 sites to accelerate searching
 * through all possible parental combinations for a window.  
 * 
 * It only begins each 64 site block and expand outward to find set number of minor
 * alleles in the target sequence, and then it looks for all possible parents.
 * 
 * TODO: change taxa by position, remove bad sites, test centromere, add hets back in, fuse overlapping 
 * haplotypes into chromosome long ones
 * 
 * @author edbuckler, kswarts, aromero
 */
public class MinorWindowImputation2 {
    private Alignment unimpAlign;
    private Alignment donorAlign;
    private int testing=0;  //level of reporting to stdout
    private OpenBitSet swapMjMnMask, invariantMask, errorMask, goodMask;  //mask of sites that major and minor are swapped
    private int parentsRight=0, parentsWrong=0;
    private int blocks=-1;
    private int[] siteErrors, siteCallCnt;
    int minMinorCnt=15;
    boolean isSwapMajorMinor=false;
    private int maxDonorHypotheses=10;
    private double maximumInbredError=0.02;  //inbreds are tested first, if too much error hybrids are tested.
    private int minTestSites=20;  //minimum number of compared sites to find a hit
    private byte[][][] allDist;
    private boolean overLapHapFile=true;
    
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
    public MinorWindowImputation2(String donorFile, String unImpTargetFile, 
            String exportFile, int minMinorCnt) {
        this.minMinorCnt=minMinorCnt;
        donorAlign=ImportUtils.readFromHapmap(donorFile, false, (ProgressListener)null);
        donorAlign.optimizeForTaxa(null);
        System.out.printf("Donor taxa:%d sites:%d %n",donorAlign.getSequenceCount(),donorAlign.getSiteCount());        
        unimpAlign=ImportUtils.readFromHapmap(unImpTargetFile, false, (ProgressListener)null);
        unimpAlign.optimizeForTaxa(null);
        siteErrors=new int[unimpAlign.getSiteCount()];
        siteCallCnt=new int[unimpAlign.getSiteCount()];
        System.out.printf("Unimputed taxa:%d sites:%d %n",unimpAlign.getSequenceCount(),unimpAlign.getSiteCount());
        
        System.out.println("Creating mutable alignment");
        blocks=unimpAlign.getAllelePresenceForAllSites(0, 0).getNumWords();
        createMaskForAlignmentConflicts();
        MutableNucleotideAlignment mna=MutableNucleotideAlignment.getInstance(this.unimpAlign);
        int total=0, hybrid=0;
        long time=System.currentTimeMillis();
        for (int bt = 0; bt < unimpAlign.getSequenceCount(); bt+=1) {
            hybrid=0;
            System.out.println("");
            DonorHypoth[][] regionHypth=new DonorHypoth[blocks][maxDonorHypotheses];
            String name=unimpAlign.getIdGroup().getIdentifier(bt).getFullName();
            BitSet[] modBits=arrangeMajorMinorBtwAlignments(unimpAlign,bt);
            System.out.printf("Imputing %d:%s Mj:%d, Mn:%d Unk:%d ... ", bt,name,modBits[0].cardinality(),
                    modBits[1].cardinality(), countUnknown(mna,bt));
            if(modBits[0].cardinality()<100) continue;
            long[] mnT=modBits[1].getBits();
            calcInbredDist(modBits);
            for (int focusBlock = 0; focusBlock < blocks; focusBlock++) {
                int[] resultRange=getBlockWithMinMinorCount(mnT, focusBlock, minMinorCnt);
                total++;
                regionHypth[focusBlock]=getBestInbredDonors(bt, resultRange[0],resultRange[2], focusBlock);
                //if inbreds fail then try hybrids
//                if(regionHypth[focusBlock][0].getErrorRate()>maximumInbredError) {
//                    long[] mjTRange=modBits[0].getBits(resultRange[0],resultRange[2]);
//                    long[] mnTRange=modBits[1].getBits(resultRange[0],resultRange[2]);
//                    regionHypth[focusBlock]=getBestHybridDonors(bt, mjTRange, mnTRange, resultRange[0],resultRange[2], 
//                            focusBlock, regionHypth[focusBlock]);
//                    hybrid++;
//                }
            }
            solveRegionally(mna, bt, regionHypth);
            int unk=countUnknown(mna,bt);
            System.out.printf("Hybrid:%d Done Unk: %d Prop:%g", hybrid, unk, (double)unk/(double)mna.getSiteCount());
          //  System.out.println(bt+mna.getBaseAsStringRow(bt));
        }
        System.out.println("");
        System.out.println("Time:"+(time-System.currentTimeMillis()));
        StringBuilder s=new StringBuilder();
        s.append(String.format("%s %s MinMinor:%d ", donorFile, unImpTargetFile, minMinorCnt));
        if(testing>0) s.append(String.format("ParentsRight:%d Wrong %d", parentsRight, parentsWrong));
        System.out.println(s.toString());
        System.out.printf("TotalRight %d  TotalWrong %d Rate%n",totalRight, totalWrong);
//        for (int i = 0; i < siteErrors.length; i++) {
//            System.out.printf("%d %d %d %g %g %n",i,siteCallCnt[i],siteErrors[i], 
//                    (double)siteErrors[i]/(double)siteCallCnt[i], unimpAlign.getMinorAlleleFrequency(i));
//        }
        ExportUtils.writeToHapmap(mna, false, exportFile, '\t', null);
        System.out.printf("%d %g %d ",minMinorCnt, maximumInbredError, maxDonorHypotheses);
    }
    
    /**
     * Create mask for all sites where major & minor are swapped in alignments
     */
    private void createMaskForAlignmentConflicts() {
        goodMask=new OpenBitSet(unimpAlign.getSiteCount());
        errorMask=new OpenBitSet(unimpAlign.getSiteCount());
        swapMjMnMask=new OpenBitSet(unimpAlign.getSiteCount());
        invariantMask=new OpenBitSet(unimpAlign.getSiteCount());
        int siteConflicts=0, swaps=0, invariant=0, good=0;
        for (int i = 0; i < unimpAlign.getSiteCount(); i++) {
            /*we have three classes of data:  invariant in one alignment, conflicts about minor and minor,
            *swaps of major and minor.  Adding the invariant reduces imputation accuracy.
            *the major/minor swaps should be flipped in the comparisons
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
        System.out.println("invariant in donor:"+invariant+" swapConflicts:"+swaps+" errors:"+siteConflicts);
    }
    
    private int countUnknown(Alignment a, int taxon) {
        int cnt=0;
        for (int i = 0; i < a.getSiteCount(); i++) {
            if(a.getBase(taxon, i)==Alignment.UNKNOWN_DIPLOID_ALLELE) cnt++;
        }
        return cnt;
    }
    
    private int[] findTaxaOffsetForOverLap(int totalTaxa, int site) {
        int tOff=totalTaxa/3;
        int mod=((site-2048)/2048)%3;
        int start=tOff*mod;
        return (new int[] {start, start+tOff});
    }
    
    private BitSet[] arrangeMajorMinorBtwAlignments(Alignment unimpAlign, int bt) { 
        OpenBitSet mjTbs=new OpenBitSet(unimpAlign.getAllelePresenceForAllSites(bt, 0));
        OpenBitSet mnTbs=new OpenBitSet(unimpAlign.getAllelePresenceForAllSites(bt, 1));
        mjTbs.and(goodMask);
        mnTbs.and(goodMask);
        if(isSwapMajorMinor) {
            OpenBitSet newmj=new OpenBitSet(unimpAlign.getAllelePresenceForAllSites(bt, 1));
            OpenBitSet newmn=new OpenBitSet(unimpAlign.getAllelePresenceForAllSites(bt, 0));
            newmj.and(swapMjMnMask);
            newmn.and(swapMjMnMask);
            mjTbs.or(newmj);
            mnTbs.or(newmn);
        }
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
            DonorHypoth[][] regionHypth) {
        for (int focusBlock = 0; focusBlock < blocks; focusBlock++) {
            DonorHypoth cbh=regionHypth[focusBlock][0];
            if(regionHypth[focusBlock][0].mendelianErrors==0) {
                setAlignmentWithDonors(regionHypth[focusBlock],focusBlock,mna);
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
                if(regionHypth[currBestBlock][0].getErrorRate()<this.maximumInbredError) {
                    setAlignmentWithDonors(regionHypth[currBestBlock],focusBlock,mna);
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
            if(startBlock==0) {preferMoveStart=false;}
            if(endBlock==blocks-1) {preferMoveStart=true;}
            if((startBlock==0)&&(endBlock==blocks-1)) break;
            if(preferMoveStart) {//expand start
                startBlock--;
                minorCnt+=Long.bitCount(mnT[startBlock]);
            } else { //expand end
                endBlock++;
                minorCnt+=Long.bitCount(mnT[endBlock]);  
            } 
        }
        int[] result={startBlock, focusBlock, endBlock};
        return result;
    }
    
    /**
     * Calculates an inbred distance
     * TODO:  Should be converted to the regular same, diff, hets, & sum
     * @param modBitsOfTarget
     * @return array with testSites, minorSites, unmatched (errors)
     */
    private void calcInbredDist(BitSet[] modBitsOfTarget) {
        allDist=new byte[donorAlign.getSequenceCount()][3][blocks];
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
    private DonorHypoth[] getBestInbredDonors(int targetTaxon, int startBlock, int endBlock, int focusBlock) {
        TreeMap<Double,DonorHypoth> bestDonors=new TreeMap<Double,DonorHypoth>();
        bestDonors.put(1.0, new DonorHypoth());
        double lastKeytestPropUnmatched=1.0;
        int donorTaxaCnt=donorAlign.getSequenceCount();
        int startDonor=0, endDonor=donorTaxaCnt;
        if(overLapHapFile) {  //gets the taxa stagger for donor file
            int[] tRange=findTaxaOffsetForOverLap(donorTaxaCnt, (focusBlock*64)+32);
            startDonor=tRange[0];
            endDonor=tRange[1];
        }
        for (int d1 = startDonor; d1 < endDonor; d1++) {
            int misMatch=0;
            int testSites=0;
            for (int i = startBlock; i <=endBlock; i++) {
                misMatch+=allDist[d1][2][i];
                testSites+=allDist[d1][0][i];
                if(misMatch>5) break;  //TODO consider whether this is good
            }
            if(testSites<minTestSites) continue;
            int totalMendelianErrors=misMatch;
            double testPropUnmatched=(double)(totalMendelianErrors)/(double)testSites;
 //this is the key line
            if(testPropUnmatched>this.maximumInbredError) continue; //toss all matches above the inbred threshold 
            testPropUnmatched+=(double)d1/(double)(donorTaxaCnt*1000);  //this looks strange, but makes the keys unique and ordered
            if((bestDonors.size()<maxDonorHypotheses)) {
                DonorHypoth theDH=new DonorHypoth(targetTaxon, d1, d1, startBlock,
                    focusBlock, endBlock, testSites, totalMendelianErrors);
                DonorHypoth prev=bestDonors.put(new Double(testPropUnmatched), theDH);
                if(prev!=null) System.out.println("TreeMap index crash:"+testPropUnmatched);
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
            int startBlock, int endBlock, 
            int focusBlock, DonorHypoth[] inbredHypoth) {
        TreeMap<Double,DonorHypoth> bestDonors=new TreeMap<Double,DonorHypoth>();
        bestDonors.put(1.0, new DonorHypoth());
        double lastKeytestPropUnmatched=1.0;
        int donorTaxaCnt=donorAlign.getSequenceCount();
        int startDonor=0, endDonor=donorTaxaCnt;
        if(overLapHapFile) {  //gets the taxa stagger for donor file
            int[] tRange=findTaxaOffsetForOverLap(donorTaxaCnt, (focusBlock*64)+32);
            startDonor=tRange[0];
            endDonor=tRange[1];
        }
        for (int iH = 0; iH < inbredHypoth.length; iH++) {
            int d1 = inbredHypoth[iH].donor1Taxon; 
            if(d1<0)    { 
       //         System.out.printf("d1: %d startBlock:%d  endBlock:%d mjT.len:%d %n",d1, startBlock,  endBlock, mjT.length);
                continue;
            }
            long[] mj1=donorAlign.getAllelePresenceForSitesBlock(d1, 0, startBlock, endBlock+1);
            long[] mn1=donorAlign.getAllelePresenceForSitesBlock(d1, 1, startBlock, endBlock+1);
            for (int d2 = startDonor; d2 < endDonor; d2++) {
                long[] mj2=donorAlign.getAllelePresenceForSitesBlock(d2, 0, startBlock, endBlock+1);
                long[] mn2=donorAlign.getAllelePresenceForSitesBlock(d2, 1, startBlock, endBlock+1);
                
                //if(testSites<minTestSites) continue;
         //       int totalMendelianErrors=misMatch;
                int[] mendErr=mendelErrorComparison(mjT, mnT, mj1, mn1, mj2, mn2);
                double testPropUnmatched=(double)(mendErr[0])/(double)mendErr[1];
                if(testPropUnmatched>this.maximumInbredError) continue; //toss all matches above the inbred threshold 
                testPropUnmatched+=(double)d2/(double)(donorTaxaCnt*1000);  //this looks strange, but makes the keys unique and ordered
                testPropUnmatched+=(double)d1/(double)(donorTaxaCnt*1000);
                if((bestDonors.size()<maxDonorHypotheses)) {
                    DonorHypoth theDH=new DonorHypoth(targetTaxon, d1, d2, startBlock,
                        focusBlock, endBlock, mendErr[1], mendErr[0]);
                    DonorHypoth prev=bestDonors.put(new Double(testPropUnmatched), theDH);
                    if(prev!=null) System.out.println("TreeMap index crash:"+testPropUnmatched);
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
//        long[] mj1=donorAlign.getAllelePresenceForAllSites(d1, 0).getBits();
//        long[] mn1=donorAlign.getAllelePresenceForAllSites(d1, 1).getBits();
//        long[] mj2=donorAlign.getAllelePresenceForAllSites(d2, 0).getBits();
//        long[] mn2=donorAlign.getAllelePresenceForAllSites(d2, 1).getBits();
        int mjUnmatched=0;
        int mnUnmatched=0;
        int testSites=0;
        int testTargetMajor=0;
        int testTargetMinor=0;
        for (int i = 0; i < mjT.length; i++) {
            long siteMask=(mjT[i]|mnT[i])&(mj1[i]|mn1[i])&(mj2[i]|mn2[i]);
            mjUnmatched+=Long.bitCount(siteMask&mjT[i]&(mjT[i]^mj1[i])&(mjT[i]^mj2[i]));
            mnUnmatched+=Long.bitCount(siteMask&mnT[i]&(mnT[i]^mn1[i])&(mnT[i]^mn2[i]));
            testSites+=Long.bitCount(siteMask);
            testTargetMajor+=Long.bitCount(siteMask&mjT[i]);
            testTargetMinor+=Long.bitCount(siteMask&mnT[i]);
        }
        int totalMendelianErrors=mjUnmatched+mnUnmatched;
       // double testPropUnmatched=(double)(totalMendelianErrors)/(double)testSites;
        return (new int[] {totalMendelianErrors, testSites});       
    }
    
 
    
        /**
     * Takes a donor hypothesis and applies it to the output alignment 
     * @param theDH
     * @param mna 
     */
    private void setAlignmentWithDonors(DonorHypoth[] theDH, int focusBlock, MutableNucleotideAlignment mna) {
        if(theDH[0].targetTaxon<0) return;
        boolean print=false;
        int startSite=focusBlock*64;
        int endSite=startSite+63;
        if(endSite>=unimpAlign.getSiteCount()) endSite=unimpAlign.getSiteCount()-1;
        if (print) System.out.println("B:"+mna.getBaseAsStringRange(theDH[0].targetTaxon, startSite, endSite));
        for(int cs=startSite; cs<=endSite; cs++) {
            byte donorEst=Alignment.UNKNOWN_DIPLOID_ALLELE; 
            for (int i = 0; (i < theDH.length) && (donorEst==Alignment.UNKNOWN_DIPLOID_ALLELE); i++) {
                if((theDH[i]==null)||(theDH[i].donor1Taxon<0)) continue;
                byte bD1=donorAlign.getBase(theDH[i].donor1Taxon, cs);
                byte bD2=donorAlign.getBase(theDH[i].donor2Taxon, cs);
                donorEst=Alignment.UNKNOWN_DIPLOID_ALLELE;
                if(bD1==Alignment.UNKNOWN_DIPLOID_ALLELE) {donorEst=bD2;}
                else if(bD2==Alignment.UNKNOWN_DIPLOID_ALLELE) {donorEst=bD1;}
                else {donorEst=(byte)((bD1&highMask)|(bD2&lowMask));}
            //    if(donorEst!=Alignment.UNKNOWN_DIPLOID_ALLELE) break;
            }
            byte knownBase=mna.getBase(theDH[0].targetTaxon, cs);
            if((knownBase!=Alignment.UNKNOWN_DIPLOID_ALLELE)&&(knownBase!=donorEst)&&(donorEst!=Alignment.UNKNOWN_DIPLOID_ALLELE)) {
//                System.out.printf("Error %d %s %s %n",theDH.targetTaxon, mna.getBaseAsString(theDH.targetTaxon, cs),
//                        NucleotideAlignmentConstants.getNucleotideIUPAC(donorEst));      
                if(!AlignmentUtils.isHeterozygous(donorEst)) {
                    totalWrong++;
                    siteErrors[cs]++;
                }
            }
            if((knownBase!=Alignment.UNKNOWN_DIPLOID_ALLELE)&&(knownBase==donorEst)&&(donorEst!=Alignment.UNKNOWN_DIPLOID_ALLELE)) {
//                System.out.printf("Error %d %s %s %n",theDH.targetTaxon, mna.getBaseAsString(theDH.targetTaxon, cs),
//                        NucleotideAlignmentConstants.getNucleotideIUPAC(donorEst));
                totalRight++;
                siteCallCnt[cs]++;
            }
            mna.setBase(theDH[0].targetTaxon, cs, donorEst);
        } //end of cs loop
        if (print) System.out.println("E:"+mna.getBaseAsStringRange(theDH[0].targetTaxon, startSite, endSite));
    }
    
    private static void compareAlignment(String origFile, String maskFile, String impFile) {
        Alignment oA=ImportUtils.readFromHapmap(origFile, false, (ProgressListener)null);
        System.out.printf("Orig taxa:%d sites:%d %n",oA.getSequenceCount(),oA.getSiteCount());        
        Alignment mA=ImportUtils.readFromHapmap(maskFile, false, (ProgressListener)null);
        System.out.printf("Mask taxa:%d sites:%d %n",mA.getSequenceCount(),mA.getSiteCount());
        Alignment iA=ImportUtils.readFromHapmap(impFile, false, (ProgressListener)null);
        System.out.printf("Imp taxa:%d sites:%d %n",iA.getSequenceCount(),iA.getSiteCount());
        int correct=0;
        int errors=0;
        int unimp=0;
        int hets=0;
        for (int t = 0; t < oA.getSequenceCount(); t++) {
            for (int s = 0; s < oA.getSiteCount(); s++) {
                if(oA.getBase(t, s)!=mA.getBase(t, s)) {
                    byte ib=iA.getBase(t, s);
                    byte ob=oA.getBase(t, s);
                    if(ib==Alignment.UNKNOWN_DIPLOID_ALLELE) {unimp++;}
                    else if(ib==ob) {
                        correct++;
                    } else {
                        if(AlignmentUtils.isHeterozygous(ob)) {hets++;}
                        else {errors++;}
                    }
                }       
            }
        }
        System.out.printf("Unimp:%d Hets:%d Correct:%d Errors:%d %n",unimp,hets,correct,errors);
    }
      
    
    /**
     *
     * @param args
     */
    public static void main(String[] args) {
      String root="/Users/edbuckler/SolexaAnal/GBS/build20120701/impResults/";

 //       String donorFile=root+"DTMAfounder20120110.imp.hmp.txt";
//        String origFile=root+"04_PivotMergedTaxaTBT.c10_s0_s4095.hmp.txt.gz";
//        String donorFile=root+"masked4096Merge20130429.hmp.txt.gz";
//        String unImpTargetFile=root+"04_PivotMergedTaxaTBT.c10_s0_s4095_masked.hmp.txt.gz";
//        String impTargetFile=root+"4E20Minor2MxDon.c10_s0_s4095_masked.imp.hmp.txt.gz";
        
        String origFile=root+"all25.c10_s0_s8191.hmp.txt.gz";
        String donorFile=root+"masked8191Merge20130429.hmp.txt.gz";
        String unImpTargetFile=root+"all25.c10_s0_s8191_masked.hmp.txt.gz";
        String impTargetFile=root+"fixswap2Err20Minor10MxDon.c10_s0_s8191_masked.imp.hmp.txt.gz";

        MinorWindowImputation2 e64NNI;
        System.out.println("Resolve Method 0: Minor 20");
        e64NNI=new MinorWindowImputation2(donorFile, unImpTargetFile, impTargetFile,20);
        compareAlignment(origFile,unImpTargetFile,impTargetFile);
//        e64NNI=new MinorWindowImputation(donorFile,
//                unImpTargetFile, impTargetFile,15,1);
//         e64NNI=new MinorWindowImputation(donorFile2,
//                unImpTargetFile, impTargetFile,15,1);


    }
    
}

