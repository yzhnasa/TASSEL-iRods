/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.pal.popgen;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import java.util.TreeMap;
import net.maizegenetics.pal.alignment.*;
import net.maizegenetics.pal.statistics.ChiSquareDistribution;
import net.maizegenetics.pipeline.EdTests;
import net.maizegenetics.util.BitUtil;
import net.maizegenetics.util.ProgressListener;

/**
 * Finds the nearest neighbor for every 64 site window.  In case of ties, it 
 * extends to neighboring 64bp windows.  It can be restricted to only find 
 * neighbors for a taxa within a high density genotype set.  This is useful for 
 * projection.  This class is expected to be used with Projection alignment.
 * 
 * @author edbuckler
 */
public class Expanding64NNImputation {
    private byte[][] same, diff, hets;
    private TBitAlignment ldAlign;
    int minSites=256;
    int maxWindow=2048/64;
    double minIdentityDiff=0.01;
    int[][][] null64Share; //region, site cnt, siteIdentity
    float[][][] null64ShareProb; //site cnt, minor cnt, siteIdentity
    int blocks=-1;
    int[] hSite, hTaxon;
    byte[] hState;
    int maskSitCnt=0;
    int maxNN=10;
    double minProb=0.01;
    boolean maskAndTest=true;

    public Expanding64NNImputation(TBitAlignment ldAlign, String exportFile) {
        this.ldAlign=ldAlign;
        if(maskAndTest) maskSites(300);
        blocks=ldAlign.getAllelePresenceForAllSites(0, 0).getNumWords();
        this.createNull64Share(ldAlign, 40000000);
        MutableNucleotideAlignment mna=MutableNucleotideAlignment.getInstance(this.ldAlign);
        for (int bt = 0; bt < ldAlign.getSequenceCount(); bt++) {
            float[][] idp=getTaxaIdentityProbMatrix(bt);
            for (int x = 0; x < idp[0].length; x++) {
                int startSite=x*64;
                int endSite=startSite+63;
                if(endSite>=ldAlign.getSiteCount()) endSite=ldAlign.getSiteCount()-1;
                TreeMap<Double, ShareSize> bestTaxa=new TreeMap<Double, ShareSize>();
                for (int t = 0; t < idp.length; t++) {
                    ShareSize xss=new ShareSize(bt,t, x, x);
                    ShareSize fss=getMaxShare(idp,xss);
                    if(fss.p>minProb) continue;
                    if((bestTaxa.size()<maxNN)||(fss.p<bestTaxa.lastEntry().getKey())) {
                        bestTaxa.put(fss.p, fss);
                        if(bestTaxa.size()>maxNN) bestTaxa.remove(bestTaxa.lastEntry().getKey());
                    }
                    //System.out.printf("%g\t",idp[t][x]);   
                }
                for(ShareSize c: bestTaxa.values()) {
                    int ct=c.compTaxon;
                    for(int cs=startSite; cs<=endSite; cs++) {
                        if(mna.getBase(bt, cs)==Alignment.UNKNOWN_DIPLOID_ALLELE) {
                            mna.setBase(bt, cs, ldAlign.getBase(ct, cs));
                        }
                    }
                }
            }
        }
        if(maskAndTest) compareSites(mna);
        ExportUtils.writeToHapmap(mna, false, exportFile, '\t', null);
      
       //if we put the share size in the tree map, only remove those at a transiti0n boundary 
       // System.out.printf("L%d R%d p:%g %n",ss.left, ss.right, ss.p);
    }
    
    private void maskSites(int sampIntensity) {
        System.out.println("Beginning to mask sites");
        MutableNucleotideAlignment mna=MutableNucleotideAlignment.getInstance(ldAlign);
        int maxSites=mna.getSequenceCount()*((mna.getSiteCount()/sampIntensity)+1);
        hSite=new int[maxSites];
        hTaxon=new int[maxSites];
        hState=new byte[maxSites];
        int cnt=0;
        for (int t = 0; t < mna.getSequenceCount(); t++) {
            for (int s = t; s < mna.getSiteCount(); s+=sampIntensity) {
                hSite[cnt]=s;
                hTaxon[cnt]=t;
                hState[cnt]=mna.getBase(t, s);
                mna.setBase(t, s, Alignment.UNKNOWN_DIPLOID_ALLELE);
                cnt++;
            }  
        }
        maskSitCnt=cnt;
        mna.clean();
        compareSites(ldAlign);
        ldAlign=TBitAlignment.getInstance(mna);
        compareSites(ldAlign);
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
    
    private ShareSize getMaxShare(float[][] idp, ShareSize currShare) {
        if(currShare.left==currShare.right) {
            currShare.fsum=idp[currShare.compTaxon][currShare.left];
            currShare.p=1.0-ChiSquareDistribution.cdf(currShare.fsum, 2);
        }
        double tL=-1, tR=-1;
        if(currShare.left>0) {
            tL=idp[currShare.compTaxon][currShare.left-1];
        }
        if(currShare.right<idp[0].length-1) {
            tR=idp[currShare.compTaxon][currShare.right+1];
        }
        if(tL>tR) {
            double testFsum=currShare.fsum+tL;
            double testp=1.0-ChiSquareDistribution.cdf(testFsum, currShare.df+2);
            if(testp<currShare.p) {
                currShare.moveLeft(testp, testFsum);
                return getMaxShare(idp, currShare);
            }
        } else {
            double testFsum=currShare.fsum+tR;
            double testp=1.0-ChiSquareDistribution.cdf(testFsum, currShare.df+2);
            if(testp<currShare.p) {
                currShare.moveRight(testp, testFsum); 
                return getMaxShare(idp, currShare);
            }
        }
        return currShare;
    }

    private void createNull64Share(TBitAlignment a, int maxSampling) {
        null64Share=new int[blocks][65][65];
        for (int i = 0; i < blocks; i++) {
            for (int j = 0; j < null64Share[0].length; j++) {
                for (int k = 0; k <=j; k++) {
                    null64Share[i][j][k]=1;
                }                
            } 
        }
        Random r=new Random(0);
        int samplingPerTaxon=maxSampling/(a.getSequenceCount()*a.getSequenceCount()/2);
        System.out.println("samplingPerTaxonContrast:"+samplingPerTaxon);
        for (int t1 = 0; t1 < a.getSequenceCount(); t1++) {
            long[] iMj=a.getAllelePresenceForAllSites(t1, 0).getBits();
            long[] iMn=a.getAllelePresenceForAllSites(t1, 1).getBits();
            for (int t2 = 0; t2 < a.getSequenceCount(); t2++) {
                if(t1==t2) continue;
                long[] jMj=a.getAllelePresenceForAllSites(t2, 0).getBits();
                long[] jMn=a.getAllelePresenceForAllSites(t2, 1).getBits();
                for (int sN = 0; sN < blocks; sN++) {
               //     int br=r.nextInt(numBins);
                    int b=sN;
                    int[] results=this.getIdentity(iMj[b], iMn[b], jMj[b], jMn[b]);
              //      if((sN==0)&&((t1==39)||(t2==39))) System.out.printf("%d %d %d %s %n",sN, t1, t2, Arrays.toString(results));
                    null64Share[sN][results[0]][results[2]]++;
                }
            }
        }
        null64ShareProb=new float[blocks][65][65];
        for (int i = 0; i < blocks; i++) {
            for (int j = 0; j < null64Share[0].length; j++) {
                int sum=0, bsum=0;
                for (int k = 0; k <=j; k++) {
                    sum+=null64Share[i][j][k];
                }
                for (int k = j; k >=0; k--) {
                    bsum+=null64Share[i][j][k];
                  //  null64ShareProb[i][j][k]=(float)bsum/(float)sum;
                    null64ShareProb[i][j][k]=(float)(-2*Math.log((double)bsum/(double)sum));
                }                
            } 
        }
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
    
    private float[][] getTaxaIdentityProbMatrix(int taxa) {
        long[] iMj=ldAlign.getAllelePresenceForAllSites(taxa, 0).getBits();
        long[] iMn=ldAlign.getAllelePresenceForAllSites(taxa, 1).getBits();
        int sections=iMj.length;
        float[][] result=new float[ldAlign.getSequenceCount()][sections];
        for (int t = 0; t < ldAlign.getSequenceCount(); t++) {
            long[] jMj=ldAlign.getAllelePresenceForAllSites(t, 0).getBits();
            long[] jMn=ldAlign.getAllelePresenceForAllSites(t, 1).getBits();
            for(int x=0; x<sections; x++) {
                int[] results=this.getIdentity(iMj[x], iMn[x], jMj[x], jMn[x]);
                result[t][x]=null64ShareProb[x][results[0]][results[2]];
            }
        }
        return result;
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
    
    
    public static void main(String[] args) {
        String root="/Users/edbuckler/SolexaAnal/bigprojection/";
        String gFile=root+"282_Zea20120110_scv10mF8maf002_mgs_E1pLD5kpUn.c10.hmp.txt";
        String hFile=root+"maizeHapMapV2_B73RefGenV2_201203028_chr10.hmp.txt";
        String mFile=root+"merge.chr10.hmp.txt";
        String exFile=root+"merge_01.chr10.imp.hmp.txt";
        boolean buildInput=false;
        if(buildInput) {
            TBitAlignment gbsMap=TBitAlignment.getInstance(ImportUtils.readFromHapmap(gFile, (ProgressListener)null));
            System.out.println("GBS Map Read");
    //        SBitAlignment hapMap=(SBitAlignment)readGZOfSBit(hapFileAGP1, true);
            SBitAlignment hapMap=(SBitAlignment)ImportUtils.readFromHapmap(hFile, (ProgressListener)null);
            System.out.println("HapMap Read");
            hapMap=(SBitAlignment)EdTests.fixHapMapNames(hapMap);  //adds tags so that HapMapNames are recognizable
            System.out.println("HapMap Names Fixed");
            MutableNucleotideAlignment mna=EdTests.combineAlignments(hapMap, gbsMap);
            System.out.println("HapMap and GBS combined");
            mna.clean();
            ExportUtils.writeToHapmap(mna, false, mFile, '\t', null);
        }
        TBitAlignment mergeMap=TBitAlignment.getInstance(ImportUtils.readFromHapmap(mFile, (ProgressListener)null));
        Expanding64NNImputation e64NNI=new Expanding64NNImputation(mergeMap, exFile);
    //    TBitAlignment mergeMap=TBitAlignment.getInstance(mna);
    }
    
}

class ShareSize {
    int baseTaxon=-1;
    int compTaxon=-1;
    int left=-1;
    int right=-1;
    double p=1;
    double fsum=0;
    int df=0;

    public ShareSize(int baseTaxon, int compTaxon, int left, int right, double p, double fsum) {
        this(baseTaxon, compTaxon, left, right);
        this.p=p;
        this.fsum=fsum;
    }
    
    public ShareSize(int baseTaxon, int compTaxon, int left, int right) {
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
