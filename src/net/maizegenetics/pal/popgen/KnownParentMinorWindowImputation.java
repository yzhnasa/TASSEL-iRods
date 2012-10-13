/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.pal.popgen;

import java.util.ArrayList;
import java.util.Random;
import net.maizegenetics.pal.alignment.*;
import net.maizegenetics.pal.ids.Identifier;
import net.maizegenetics.pal.statistics.ApproxFastChiSquareDistribution;
import net.maizegenetics.util.BitUtil;
import net.maizegenetics.util.OpenBitSet;
import net.maizegenetics.util.ProgressListener;


/**
 * Finds the nearest neighbor for every 64 site window.  In case of ties, it 
 * extends to neighboring 64bp windows.  
 * 
 * @author edbuckler
 */
public class KnownParentMinorWindowImputation {
    private Alignment unimpAlign;
    private Alignment donorAlign;
    private int testing=2;
    private OpenBitSet swapMjMnMask;
//    int minSites=256;
//    int maxWindow=2048/64;
    double minIdentityDiff=0.01;
//    int[][][] null64Share; //region, site cnt, siteIdentity
//    float[][][] null64ShareProb; //region, site cnt, siteIdentity
    int blocks=-1;
    int[] hSite, hTaxon;
    byte[] hState;
    int maskSitCnt=0;
    int maxNN=50;
    double minProb=0.0001;
    boolean maskAndTest=true;
    ApproxFastChiSquareDistribution fcs=new ApproxFastChiSquareDistribution(1000,200);

    public KnownParentMinorWindowImputation(String donorFile, String unImpTargetFile, String exportFile, int minMinorCnt) {
        donorAlign=ImportUtils.readFromHapmap(donorFile, false, (ProgressListener)null);
        donorAlign.optimizeForTaxa(null);
        unimpAlign=ImportUtils.readFromHapmap(unImpTargetFile, false, (ProgressListener)null);
        unimpAlign.optimizeForTaxa(null);
        
        if(maskAndTest) maskSites(300);
        blocks=unimpAlign.getAllelePresenceForAllSites(0, 0).getNumWords();
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
        int impSiteCnt=0;
        Random r=new Random(0);
        for (int bt = 0; bt < unimpAlign.getSequenceCount(); bt++) {
            int taxaImpCnt=0;
            String name=unimpAlign.getIdGroup().getIdentifier(bt).getFullName();
            System.out.printf("Imputing %d:%s ... %n", bt,name);
            long time=System.currentTimeMillis();
            long[] mjT=unimpAlign.getAllelePresenceForAllSites(bt, 0).getBits();
            long[] mnT=unimpAlign.getAllelePresenceForAllSites(bt, 1).getBits();
            for (int startBlock = 0; startBlock < blocks; startBlock++) {
                int minorCnt=Long.bitCount(mnT[startBlock]);
                int majorCnt=Long.bitCount(mjT[startBlock]);
                int endBlock=startBlock;
                while((minorCnt<minMinorCnt)&&(endBlock<(blocks-1))) {
                    minorCnt+=Long.bitCount(mnT[++endBlock]);
                    majorCnt+=Long.bitCount(mjT[endBlock]);
                }
//                System.out.printf("Taxa %d Name: %s blocks: %d end: %d MjCnt %d MnCnt %d %n", bt, name,
//                        (1+endBlock-startBlock), endBlock, minorCnt, majorCnt);
                int[] donors=getBestDonors(bt, startBlock, endBlock);
//                donors[0]=r.nextInt(donorAlign.getSequenceCount());
//                donors[0]=r.nextInt(donorAlign.getSequenceCount());
                int startSite=startBlock*64;//TODO move window to middle
                int endSite=startSite+63;
                if(endSite>=unimpAlign.getSiteCount()) endSite=unimpAlign.getSiteCount()-1;
                int highMask=15<<4;
                int lowMask=15;
                for(int cs=startSite; cs<=endSite; cs++) {
                    byte bD1=donorAlign.getBase(donors[0], cs);
                    byte bD2=donorAlign.getBase(donors[1], cs);
                    byte donorEst=Alignment.UNKNOWN_DIPLOID_ALLELE;
                    if(bD1==Alignment.UNKNOWN_DIPLOID_ALLELE) {donorEst=bD2;}
                    else if(bD2==Alignment.UNKNOWN_DIPLOID_ALLELE) {donorEst=bD1;}
                    else {donorEst=(byte)((bD1&highMask)|(bD2&lowMask));
//                        byte donorEst2=(bD1<bD2)?AlignmentUtils.getDiploidValue(bD1, bD2):AlignmentUtils.getDiploidValue(bD2, bD1);
//                        if(donorEst!=donorEst2) 
//                            System.out.printf("bd1:%d bd2:%d donorEst: %d donorEst2 %d %n",bD1,bD2,donorEst,donorEst2);
//                           System.out.printf("bd1:%s bd2:%s donorEst: %s donorEst2 %s %n",
//                                   NucleotideAlignmentConstants.getNucleotideIUPAC(bD1),NucleotideAlignmentConstants.getNucleotideIUPAC(bD2),
//                                   NucleotideAlignmentConstants.getNucleotideIUPAC(donorEst),NucleotideAlignmentConstants.getNucleotideIUPAC(donorEst2));
                    }

                    
                    if(mna.getBase(bt, cs)==Alignment.UNKNOWN_DIPLOID_ALLELE) {
                            mna.setBase(bt, cs, donorEst);
                            impSiteCnt++;
                            taxaImpCnt++;
                        }
                }
 //               System.out.printf("Best Donors D1: %d D2: %d maxTestSites: %d %n", donors[0], donors[1], donors[2]);
                if(endBlock==(blocks-1)) {
                    break;
                }
            }

//            System.out.printf("Finished %d Imp %d %d %n", System.currentTimeMillis()-time, impSiteCnt, taxaImpCnt);
//            if(bt%10==0) compareSites(mna);
        }
        if(maskAndTest) compareSites(mna);
        ExportUtils.writeToHapmap(mna, false, exportFile, '\t', null);
      
       //if we put the share size in the tree map, only remove those at a transiti0n boundary 
       // System.out.printf("L%d R%d p:%g %n",ss.left, ss.right, ss.p);
    }
    
    private int[] getBestDonors(int targetTaxon, int startBlock, int endBlock) {
        long[] mjT=unimpAlign.getAllelePresenceForAllSites(targetTaxon, 0).getBits();
        long[] mnT=unimpAlign.getAllelePresenceForAllSites(targetTaxon, 1).getBits();
        int[] donors={-1,-1,0};
        double minPropUnmatched=1.0;
        int maxTestSites=0;
        int[] rDonors=getDonorsForRegion(unimpAlign.getTaxaName(targetTaxon),startBlock*64);
        if(testing>3) System.out.printf("StartSite %d EndSite %d RealD1 %d RealD2 %d %n",startBlock*64, 
                (endBlock*64+63),rDonors[0],rDonors[1]);
//        System.out.println("T :"+unimpAlign.getBaseAsStringRow(targetTaxon));
//        System.out.println("D1:"+donorAlign.getBaseAsStringRow(rDonors[0]));
//        System.out.println("D2:"+donorAlign.getBaseAsStringRow(rDonors[1]));
//        BitUtil.printBitLong(mjT[0]);
//        BitUtil.printBitLong(mnT[0]);
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
                    if((testing>5)&&(rDonors[0]==d1)&&(rDonors[1]==d2)) {
                        System.out.printf("mjUn %d mnUn %d %n",mjUnmatched, mnUnmatched);
                        System.out.print(i+"J:"); BitUtil.printBitLong(mjT[i]);
                        System.out.print(i+"J:"); BitUtil.printBitLong(mj1[i]);
                        System.out.print(i+"J:"); BitUtil.printBitLong(mj2[i]);
                        System.out.print(i+"N:"); BitUtil.printBitLong(mnT[i]);
                        System.out.print(i+"N:"); BitUtil.printBitLong(mn1[i]);
                        System.out.print(i+"N:"); BitUtil.printBitLong(mn2[i]);
                        System.out.println("T1:"+unimpAlign.getBaseAsStringRange(targetTaxon, i*64, (i*64+64)));
                        System.out.println("D1:"+donorAlign.getBaseAsStringRange(d1, i*64, (i*64+64)));
                        System.out.println("D2:"+donorAlign.getBaseAsStringRange(d2, i*64, (i*64+64)));
                    }
                    testSites+=Long.bitCount(siteMask);
                    testTargetMajor+=Long.bitCount(siteMask&mjT[i]);
                    testTargetMinor+=Long.bitCount(siteMask&mnT[i]);
                }
                double testPropUnmatched=(double)(mjUnmatched+mnUnmatched)/(double)testSites;
                if((testing>1)&&(rDonors[0]==d1)&&(rDonors[1]==d2))
                    System.out.printf("Donor %d %d %d %d %d %d %d block %d %n", d1, d2, mjUnmatched, mnUnmatched,
                            testSites, testTargetMajor, testTargetMinor, startBlock);
                if(testPropUnmatched<minPropUnmatched) {
                    donors[0]=d1;
                    donors[1]=d2;
                    minPropUnmatched=testPropUnmatched;
                    donors[2]=maxTestSites=testSites;
//                    System.out.printf("Better %d %d %d %d %d %d %d %n", d1, d2, mjUnmatched, mnUnmatched,
//                            testSites, testTargetMajor, testTargetMinor);
                }
                
            }
            
        }
        return donors;
    }
    
    private int[] getDonorsForRegion(String taxonName, int site) {
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
    
    private void compareSites(Alignment a) {
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
        System.out.printf("Missing: %d Het: %d NotImp: %d Error: %d Correct: %d ErrorRate: %g %n", 
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
            for (int b = 0; b < a.getSiteCount(); b+=blockSize) {
                int p1=r.nextInt(a.getSequenceCount());
                int p2=r.nextInt(a.getSequenceCount());
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
                    
                }//end of site        
            } //end of blocks
            System.out.println(tName.toString());
            mna.setTaxonName(t, new Identifier(tName.toString()));
        }
        mna.clean();
        ExportUtils.writeToHapmap(mna, false, unImpTargetFile, '\t', null);
    }
    
    
    public static void main(String[] args) {
      String root="/Users/edbuckler/SolexaAnal/GBS/build20120110/imp/";
//        String root="/Volumes/LaCie/build20120110/imp/";

        String donorFile=root+"NAMfounder20120110seg.imp.hmp.txt";
        String unImpTargetFile=root+"ZeaSyn20120110.hmp.txt";
        String impTargetFile=root+"ZeaSyn20120110.imp.hmp.txt";

        boolean buildInput=true;
        boolean filterTrue=true;
        if(buildInput) {createSynthetic(donorFile, unImpTargetFile, 1900, 0.4, -1, 1000);}

       // System.out.println(AlignmentUtils.getDiploidValue()
        KnownParentMinorWindowImputation e64NNI=new KnownParentMinorWindowImputation(donorFile,
                unImpTargetFile, impTargetFile,20);
    }
    
}
 class ShareSizeX {
    int baseTaxon=-1;
    int compTaxon=-1;
    int left=-1;
    int right=-1;
    double p=1;
    double fsum=0;
    int df=0;

    public ShareSizeX(int baseTaxon, int compTaxon, int left, int right, double p, double fsum) {
        this(baseTaxon, compTaxon, left, right);
        this.p=p;
        this.fsum=fsum;
    }
    
    public ShareSizeX(int baseTaxon, int compTaxon, int left, int right) {
        this.baseTaxon=baseTaxon;
        this.compTaxon=compTaxon;
        this.left=left;
        this.right=right;
        df=(right-left+1)*2;
    }
    
    public void moveLeft(double p, double fsum) {
        left--;
        this.p=p;
        this.fsum=fsum;
        df=df+2;
    }
    
    public void moveRight(double p, double fsum) {
        right++;
        this.p=p;
        this.fsum=fsum;
        df=df+2;
    }
    
    public String toString() {
        return String.format("BTx:%d CTx:%d L:%d R:%d FSum:%g P:%g ", baseTaxon,
                compTaxon, left, right, fsum, p);
    }
}
