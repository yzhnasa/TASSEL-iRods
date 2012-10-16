/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.pal.popgen;

import java.util.ArrayList;
import java.util.Random;
import net.maizegenetics.pal.alignment.*;
import net.maizegenetics.pal.ids.Identifier;
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
public class KnownParentMinorWindowImputation {
    private Alignment unimpAlign;
    private Alignment donorAlign;
    private int testing=2;  //level of reporting to stdout
    private OpenBitSet swapMjMnMask;
    private int parentsRight=0, parentsWrong=0;
     int blocks=-1;
    int[] hSite, hTaxon;
    byte[] hState;
    int maskSitCnt=0;
    int maxNN=50;
    boolean maskAndTest=true;
    
    private static int highMask=0xF0;
    private static int lowMask=0x0F;

    public KnownParentMinorWindowImputation(String donorFile, String unImpTargetFile, String exportFile, int minMinorCnt) {
        donorAlign=ImportUtils.readFromHapmap(donorFile, false, (ProgressListener)null);
        donorAlign.optimizeForTaxa(null);
        unimpAlign=ImportUtils.readFromHapmap(unImpTargetFile, false, (ProgressListener)null);
        unimpAlign.optimizeForTaxa(null);
        
        if(maskAndTest) maskSites(300);  //mask sites to test imputation accuracy
        blocks=unimpAlign.getAllelePresenceForAllSites(0, 0).getNumWords();
        
        //Create mask for all sites where major & minor are swapped in alignments
        swapMjMnMask=new OpenBitSet(unimpAlign.getSiteCount());
        int swapConflicts=0;
        for (int i = 0; i < unimpAlign.getSiteCount(); i++) {
            if((donorAlign.getMajorAllele(i)!=unimpAlign.getMajorAllele(i))||(donorAlign.getMinorAllele(i)!=unimpAlign.getMinorAllele(i))) {
                swapConflicts++;
                swapMjMnMask.set(i);
            }  
        }
        swapMjMnMask.not();
        System.out.println("swapConflicts"+swapConflicts+" same:"+swapMjMnMask.cardinality());
        MutableNucleotideAlignment mna=MutableNucleotideAlignment.getInstance(this.unimpAlign);
        Random r=new Random(0);
        for (int bt = 0; bt < unimpAlign.getSequenceCount(); bt+=1) {
            int taxaImpCnt=0;
            String name=unimpAlign.getIdGroup().getIdentifier(bt).getFullName();
            System.out.printf("Imputing %d:%s ... %n", bt,name);
            long time=System.currentTimeMillis();
            long[] mjT=unimpAlign.getAllelePresenceForAllSites(bt, 0).getBits();
            long[] mnT=unimpAlign.getAllelePresenceForAllSites(bt, 1).getBits();
            for (int focusBlock = 0; focusBlock < blocks; focusBlock++) {
                int[] resultRange=getBlockWithMinMinorCount(mnT, focusBlock, minMinorCnt);
                int[] donors=getBestDonors(bt, resultRange[0],resultRange[2], focusBlock);
//                donors[0]=r.nextInt(donorAlign.getSequenceCount());
//                donors[1]=r.nextInt(donorAlign.getSequenceCount());
                DonorHypoth theDH=new DonorHypoth(bt, donors[0], donors[1], resultRange[0],
                        resultRange[2], focusBlock);
                setAlignmentWithDonors(theDH,mna);


            }

//            System.out.printf("Finished %d Imp %d %d %n", System.currentTimeMillis()-time, impSiteCnt, taxaImpCnt);
//            if(bt%10==0) compareSites(mna);
        }
        StringBuilder s=new StringBuilder();
        s.append(String.format("%s %s MinMinor:%d ", donorFile, unImpTargetFile, minMinorCnt));
        if(maskAndTest) s.append(compareSites(mna));
        if(testing>0) s.append(String.format("ParentsRight:%d Wrong %d", parentsRight, parentsWrong));
        System.out.println(s.toString());
        
        ExportUtils.writeToHapmap(mna, false, exportFile, '\t', null);
      
       //if we put the share size in the tree map, only remove those at a transiti0n boundary 
       // System.out.printf("L%d R%d p:%g %n",ss.startBlock, ss.endBlock, ss.p);
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
            if(startBlock==0) preferMoveStart=false;
            if(endBlock==blocks-1) preferMoveStart=true;
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
     * Simple algorithm that tests every possible two donor combination to minimize
     * the number of unmatched informative alleles.  Currently, there is litte tie
     * breaking, longer matches are favored.
     * @param targetTaxon
     * @param startBlock
     * @param endBlock
     * @param focusBlock
     * @return int[] array of {donor1, donor2, testSites}
     */
    private int[] getBestDonors(int targetTaxon, int startBlock, int endBlock, int focusBlock) {
        long[] mjT=unimpAlign.getAllelePresenceForAllSites(targetTaxon, 0).getBits();
        long[] mnT=unimpAlign.getAllelePresenceForAllSites(targetTaxon, 1).getBits();
        int[] donors={-1,-1,0};
        double minPropUnmatched=1.0;
        int maxTestSites=0;
        int[] rDonors=parseDonorsNamesForRegion(unimpAlign.getTaxaName(targetTaxon),startBlock*64);
        if(testing>3) System.out.printf("StartSite %d EndSite %d RealD1 %d RealD2 %d %n",startBlock*64, 
                (endBlock*64+63),rDonors[0],rDonors[1]);
        long[] swapMask=this.swapMjMnMask.getBits();
        
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
                double testPropUnmatched=(double)(mjUnmatched+mnUnmatched)/(double)testSites;
                if((testing>1)&&(rDonors[0]==d1)&&(rDonors[1]==d2))
                    System.out.printf("Donor %d %d %d %d %d %d %d block %d-%d-%d %n", d1, d2, mjUnmatched, mnUnmatched,
                            testSites, testTargetMajor, testTargetMinor, startBlock, focusBlock, endBlock);
                if((testPropUnmatched<minPropUnmatched)||
                        ((testPropUnmatched==minPropUnmatched)&&(testSites>maxTestSites))) {
                    donors[0]=d1;
                    donors[1]=d2;
                    minPropUnmatched=testPropUnmatched;
                    donors[2]=maxTestSites=testSites;
                }
                
            }
            
        }
        if(testing>1){
            if((rDonors[0]==donors[0])&&(rDonors[1]==donors[1])) {System.out.println("Correct");}
            else {System.out.println("WRONG");}
        }
        if((rDonors[0]==donors[0])&&(rDonors[1]==donors[1])) {parentsRight++;}
            else {parentsWrong++;}
        return donors;
    }
    
    /**
     * Takes a donor hypothesis and applies it to the output alignment 
     * @param theDH
     * @param mna 
     */
    private void setAlignmentWithDonors(DonorHypoth theDH, MutableNucleotideAlignment mna) {
        int startSite=theDH.focusBlock*64;
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
            mna.setBase(theDH.targetTaxon, cs, donorEst);
//            if(mna.getBase(theDH.targetTaxon, cs)==Alignment.UNKNOWN_DIPLOID_ALLELE) {
//                    mna.setBase(theDH.targetTaxon, cs, donorEst);
////                    impSiteCnt++;
////                    taxaImpCnt++;
//                }
        } //end of cs loop
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
    
    private String reportTaxaMakeUp(ArrayList<Integer>[] data) {
        StringBuilder s=new StringBuilder();
        for (int i = 0; i < data[0].size(); i++) {
            s.append(data[0].get(i));
            s.append(":");
            s.append(data[1].get(i));
            s.append("\t");
        }
        return s.toString();
    }
   
    
    /**
     * 
     * @param iMj
     * @param iMn
     * @param jMj
     * @param jMn
     * @return [0]= number of sites in comparison,
     * [1]=number of minor alleles in comparison from taxon i  
     * [2]=sites that agree
     */
    private int[] getIdentity(long iMj, long iMn, long jMj, long jMn) {
        int[] results=new int[3];
        long iMnjMn=iMn&jMn;
        long iMnjMj=iMn&jMj;
        long iMnComps=iMnjMn|iMnjMj;
        long sameL=(iMj&jMj)|(iMnjMn);
        long diffL=(iMj&jMn)|(iMnjMj);
        long hetsL=sameL&diffL;
        int same=(int)BitUtil.pop(sameL);
        int diff=(int)BitUtil.pop(diffL);
        int hets=(int)BitUtil.pop(hetsL);
        results[1]=(int)BitUtil.pop(iMnComps);  
        int sum=same+diff+hets;
        results[0]=sum-(2*hets);
        results[2]=same-(hets/2); //check the minus sign
        return results;
    }
    
    private static void createSynthetic(String donorFile, String unImpTargetFile, int blockSize,
            double propPresent, double homoProp, int taxaNumber) {
        Alignment a=ImportUtils.readFromHapmap(donorFile, (ProgressListener)null);
        System.out.printf("Read %s Sites %d Taxa %d %n", donorFile, a.getSiteCount(), a.getSequenceCount());
        MutableNucleotideAlignment mna= MutableNucleotideAlignment.getInstance(a, taxaNumber, a.getSiteCount());
        Random r=new Random();
        for (int t = 0; t < taxaNumber; t++) {
            StringBuilder tName=new StringBuilder("ZM"+t);
            for (int b = 0; b < a.getSiteCount(); b+=blockSize) {  //change to bp?
                int p1=r.nextInt(a.getSequenceCount());
                int p2=r.nextInt(a.getSequenceCount());
//                p2=t%a.getSequenceCount();  //only put one crossover in
                if(p2<p1) {int temp=p1; p1=p2; p2=temp;}
                tName.append("|"+p1+"_"+p2+"s"+b);
                for (int s = b; (s < b+blockSize) && (s<a.getSiteCount()); s++) {
                    if(r.nextDouble()<propPresent) {
                        if(r.nextDouble()<0.5) {
                            mna.setBase(t, s, a.getBase(p1, s));
                        } else {
                            mna.setBase(t, s, a.getBase(p2, s));
                        }
                    } else {
                        mna.setBase(t, s, Alignment.UNKNOWN_DIPLOID_ALLELE);
                    }
                    if(r.nextDouble()<0.000) {  //method to add error to simulation 0.002 is reasonable
                        mna.setBase(t, s, AlignmentUtils.getDiploidValue(a.getMinorAllele(s),a.getMinorAllele(s)));
                    }
                }//end of site        
            } //end of blocks
            System.out.println(tName.toString());
            mna.setTaxonName(t, new Identifier(tName.toString()));
        }
        mna.clean();
        ExportUtils.writeToHapmap(mna, false, unImpTargetFile, '\t', null);
    }
    
    
    /**
     *
     * @param args
     */
    public static void main(String[] args) {
//      String root="/Users/edbuckler/SolexaAnal/GBS/build20120110/imp/";
        String root="/Volumes/LaCie/build20120110/imp/";

 //       String donorFile=root+"NAMfounder20120110.imp.hmp.txt";
        String donorFile=root+"DTMAfounder20120110.imp.hmp.txt";
        String unImpTargetFile=root+"ZeaSyn20120110.hmp.txt";
        String impTargetFile=root+"ZeaSyn20120110.imp.hmp.txt";

        boolean buildInput=true;
        boolean filterTrue=true;
        if(buildInput) {createSynthetic(donorFile, unImpTargetFile, 2000, 0.4, -1, 1000);}

        KnownParentMinorWindowImputation e64NNI=new KnownParentMinorWindowImputation(donorFile,
                unImpTargetFile, impTargetFile,30);
        
        for (int recSize = 128; recSize < 10000; recSize+=(recSize/2)) {
            for (int mm = 5; mm < 60; mm+=5) {
                System.out.println("Rec size"+recSize);
                unImpTargetFile=root+recSize+"ZeaSyn20120110.hmp.txt";
                if(buildInput) {createSynthetic(donorFile, unImpTargetFile, recSize, 0.4, -1, 1000);}
                e64NNI=new KnownParentMinorWindowImputation(donorFile, unImpTargetFile, impTargetFile,mm);
                }
        }
    }
    
}
 class DonorHypoth {
    int targetTaxon=-1;
    int donor1Taxon=-1;
    int donor2Taxon=-1;
    int startBlock=-1;
    int focusBlock=-1;
    int endBlock=-1;
    double pError=1;
    double pHeterozygous=-1, pHomoD1=-1, pHomoD2=-11;
    int totalSites=0;
    int mendelianSites=0;

    public DonorHypoth(int targetTaxon, int donor1Taxon, int donor2Taxon, int startBlock, 
            int focusBlock, int endBlock, int totalSites, int mendelianSites) {
        this(targetTaxon, donor1Taxon, donor2Taxon, startBlock, focusBlock, endBlock);
        this.totalSites=totalSites;
        this.mendelianSites=mendelianSites;
    }
    
    public DonorHypoth(int targetTaxon, int donor1Taxon, int donor2Taxon, int startBlock, 
            int focusBlock, int endBlock) {
        this.targetTaxon=targetTaxon;
        if(donor1Taxon<donor2Taxon) {
            this.donor1Taxon=donor1Taxon;
            this.donor2Taxon=donor2Taxon;
        } else {
            this.donor1Taxon=donor2Taxon;
            this.donor2Taxon=donor1Taxon;
        }
        this.startBlock=startBlock;
        this.focusBlock=focusBlock;
        this.endBlock=endBlock;
        
    }
    
    public String toString() {
        return String.format("FTx:%d D1Tx:%d D1Tx:%d SBk:%d FBk:%d EBk:%d TS:%d MS:%s ", targetTaxon,
                donor1Taxon, startBlock, focusBlock, endBlock, totalSites, mendelianSites);
    }
}
