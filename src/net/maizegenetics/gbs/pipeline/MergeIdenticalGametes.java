/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.gbs.pipeline;

import java.io.BufferedWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map.Entry;
import java.util.TreeMap;
import java.util.TreeSet;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.AlignmentUtils;
import net.maizegenetics.pal.alignment.ExportUtils;
import net.maizegenetics.pal.alignment.GeneticMap;
import net.maizegenetics.pal.alignment.ImportUtils;
import net.maizegenetics.pal.alignment.MutableNucleotideAlignment;
import net.maizegenetics.pal.alignment.NucleotideAlignmentConstants;
import net.maizegenetics.pal.alignment.ReadPolymorphismUtils;
import net.maizegenetics.pal.distance.IBSDistanceMatrix;
import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.Identifier;
import net.maizegenetics.pal.ids.SimpleIdGroup;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.BitUtil;
import net.maizegenetics.util.ExceptionUtils;
import net.maizegenetics.util.OpenBitSet;
import net.maizegenetics.util.ProgressListener;
import net.maizegenetics.util.Utils;

/**
 *
 * 
 * TODO IDEAS:  
 * Chunk chromosomes into 50Mb segments
 * Make Base HapMapV2 - use real names
 * Try to make full length haplotypes
 * plus add short inbred segments not present full ones
 * 
 * 
 * @author edbuckler
 */
public class MergeIdenticalGametes {
    private Alignment inAlign;
    private double[] propMissing;
    private int[] siteErrors, siteCallCnt;
    private double minJointGapProb=0.01;
    private boolean callGaps=false;
    private GeneticMap gm=null;
    private BitSet badMask=null;
    private double minimumMissing=0.4;
    private int maxHaplotypes=2000;
    private int siteWindow=1024*4*6;
    private int minSitesForSectionComp=50;
    private double maxHetFreq=0.02;
    private int minTaxaInGroup=2;

    public MergeIdenticalGametes(String inFile, String exportFile, String errorInFile,
            String errorExportFile, double maxDistance, int minSites) {        
        inAlign=ImportUtils.readFromHapmap(inFile, false, (ProgressListener)null);
        inAlign.optimizeForTaxa(null);
        try{
            gm=ReadPolymorphismUtils.readGeneticMapFile(errorInFile);
            badMask=maskBadSites(gm,inAlign);
            badMask=null;
        }
        catch(Exception e) {
            System.out.println("No map file with bad sites");
            badMask=null;
        }
        
        System.out.printf("In taxa:%d sites:%d %n",inAlign.getSequenceCount(),inAlign.getSiteCount());
        int sites=inAlign.getSiteCount();
        siteErrors=new int[inAlign.getSiteCount()];
        siteCallCnt=new int[inAlign.getSiteCount()];
        propMissing=new double[inAlign.getSequenceCount()];
        //Develop lists of the modest coverage, homozygous taxa
//        TreeMap<Integer,Integer> presentRanking=new TreeMap<Integer,Integer>(Collections.reverseOrder());
//        for (int i = 0; i < inAlign.getSequenceCount(); i++) {
//            int totalSitesNotMissing = inAlign.getTotalNotMissingForTaxon(i);
//            double hetFreq=(double)inAlign.getHeterozygousCountForTaxon(i)/(double)totalSitesNotMissing;
//            propMissing[i]=(double)(1+sites-totalSitesNotMissing)/(double)sites; //1 prevents error in joint missing calculations
//            double propPresent=1.0-propMissing[i];
////            System.out.printf("%s %d %g %g %n",inAlign.getTaxaName(i),totalSitesNotMissing,hetFreq, propMissing[i]);
//            if((hetFreq>maxHetFreq)||(totalSitesNotMissing<minSites)) continue;
//            int index=(1000000*((int)(propPresent*100)))+i;
////            System.out.println(index);
//            presentRanking.put(index, i);
//        }
        
        MutableNucleotideAlignment mna=createEmptyHaplotypeAlignment(maxHaplotypes);
        int siteWindow=this.siteWindow;
        int startBlock=0;
        int blockWindow=siteWindow/64;
        int blockStep=blockWindow/2;
        int sections=blockWindow;// sometimes /4
        int lastBlock= inAlign.getAllelePresenceForAllSites(0, 0).getNumWords()-1;
        for(startBlock=0; startBlock<=(lastBlock-blockWindow+1); startBlock+=blockStep) {
            int startSite=startBlock*64;
            int endBlock=startBlock+blockWindow-1;
            int mod=(startBlock/blockStep)%3;
            if((lastBlock-endBlock)<blockStep) endBlock=lastBlock;
            TreeMap<Integer,Integer> presentRanking=createPresentRankingForWindow(startBlock, endBlock, minSites, maxHetFreq);
            System.out.printf("Block %d Inbred and modest coverage:%d %n",startBlock,presentRanking.size());
            System.out.printf("Current Site %d Current block %d EndBlock: %d LastBlock:%d %n",startSite, startBlock, endBlock, lastBlock);
            TreeMap<Integer,byte[]> results=mergeWithinWindow( presentRanking, startBlock, endBlock, sections, maxDistance);
            int index=mod*maxHaplotypes;
            for (byte[] calls : results.values()) {
                mna.setBaseRange(index, startSite, calls);
                index++;
            }
        }
        mna.clean();
        ExportUtils.writeToHapmap(mna, false, exportFile, '\t', null);
        if(errorExportFile!=null) exportBadSites(errorExportFile, 0.01);
//        for (int i = 0; i < siteErrors.length; i++) {
//            System.out.printf("%d %d %d %g %g %n",i,siteCallCnt[i],siteErrors[i], 
//                    (double)siteErrors[i]/(double)siteCallCnt[i], inAlign.getMinorAlleleFrequency(i));
//        }
    }
    
    private TreeMap<Integer,Integer> createPresentRankingForWindow(int startBlock, int endBlock, 
            int minSites, double maxHetFreq) {
        int sites=64*(endBlock-startBlock+1);
        TreeMap<Integer,Integer> presentRanking=new TreeMap<Integer,Integer>(Collections.reverseOrder());
        for (int i = 0; i < inAlign.getSequenceCount(); i++) {
            long[] mj=inAlign.getAllelePresenceForSitesBlock(i, 0, startBlock, endBlock);
            long[] mn=inAlign.getAllelePresenceForSitesBlock(i, 1, startBlock, endBlock);
            int totalSitesNotMissing = 0;
            int hetCnt=0;
            for (int j = 0; j < mj.length; j++) {
                totalSitesNotMissing+=BitUtil.pop(mj[j]|mn[j]);
                hetCnt+=BitUtil.pop(mj[j]&mn[j]);
            }
            double hetFreq=(double)hetCnt/(double)totalSitesNotMissing;
            propMissing[i]=(double)(1+sites-totalSitesNotMissing)/(double)sites; //1 prevents error in joint missing calculations
            double propPresent=1.0-propMissing[i];
//            System.out.printf("%s %d %g %g %n",inAlign.getTaxaName(i),totalSitesNotMissing,hetFreq, propMissing[i]);
            if((hetFreq>maxHetFreq)||(totalSitesNotMissing<minSites)) continue;
            int index=(1000000*((int)(propPresent*100)))+i;
//            System.out.println(index);
            presentRanking.put(index, i);
        }
        return presentRanking;
    }
    
    private MutableNucleotideAlignment createEmptyHaplotypeAlignment(int maxHaplotypes) {
        IdGroup outIDG=new SimpleIdGroup(maxHaplotypes);
        for (int i = 0; i < maxHaplotypes; i++) {
            outIDG.setIdentifier(i, new Identifier("Hap"+i+"mod0"));
//            outIDG.setIdentifier(i+maxHaplotypes, new Identifier("Hap"+i+"mod1"));
//            outIDG.setIdentifier(i+(2*maxHaplotypes), new Identifier("Hap"+i+"mod2"));
        }
        MutableNucleotideAlignment mna=MutableNucleotideAlignment.getInstance(outIDG, inAlign.getSiteCount());
        for (int i = 0; i < inAlign.getSiteCount(); i++) {
            mna.addSite(i);
            mna.setLocusOfSite(i, inAlign.getLocus(i));
            mna.setPositionOfSite(i, inAlign.getPositionInLocus(i));
        }
        return mna;
    }
    
    private BitSet maskBadSites(GeneticMap gm, Alignment a) {
        OpenBitSet obs=new OpenBitSet(a.getSiteCount());
        int count=0;
        for (int i = 0; i < gm.getNumberOfMarkers(); i++) {
            int site=a.getSiteOfPhysicalPosition(gm.getPhysicalPosition(i), null);
            if(site>0) {obs.set(site);}
            
        }
        System.out.println("Bad Sites matched:"+obs.cardinality());
        obs.not();  //change all bad sites to 0, good to 1
        
        return obs;
    }
    
    private void exportBadSites(String exportMap, double errorThreshold) {
        BufferedWriter bw = null;
        try {
            String fullFileName = Utils.addSuffixIfNeeded(exportMap, ".txt", new String[]{".txt"});
            bw = Utils.getBufferedWriter(fullFileName);
            bw.write("<Map>\n");
            for (int i = 0; i < inAlign.getSiteCount(); i++) {
                double errorsRate=(double)siteErrors[i]/(double)siteCallCnt[i];
                if(errorsRate<errorThreshold) continue;
                bw.write(inAlign.getSNPID(i)+"\t");
                bw.write(inAlign.getLocusName(i) +"\t");
                bw.write(i+"\t"); //dummy for genetic position
                bw.write(inAlign.getPositionInLocus(i) +"\n"); //dummy for genetic position
            } 
            bw.close();
            
        } catch (Exception e) {
            e.printStackTrace();
            throw new IllegalArgumentException("Error writing GeneticMap file: " + exportMap + ": " + ExceptionUtils.getExceptionCauses(e));
        }
    }
    
    private TreeMap<Integer,byte[]> mergeWithinWindow(TreeMap<Integer,Integer> presentRanking,
            int startBlock, int endBlock, int sections, double maxDistance){
        int startSite=startBlock*64;
        int endSite=63+(endBlock*64);
        if(endSite>=inAlign.getSiteCount()) endSite=inAlign.getSiteCount()-1;
        TreeMap<Integer,ArrayList> mergeSets=new TreeMap<Integer,ArrayList>();
        TreeMap<Integer,byte[]> results=new TreeMap<Integer,byte[]>(Collections.reverseOrder());
        TreeSet<Integer> unmatched=new TreeSet<Integer>(presentRanking.values());
        IdGroup inIDG=inAlign.getIdGroup();
        for (Entry<Integer,Integer> e : presentRanking.entrySet()) {
            int taxon1=e.getValue();
            if(unmatched.contains(taxon1)==false) continue;//already included in another group
            ArrayList<Integer> hits=new ArrayList<Integer>();
            unmatched.remove(taxon1);
            for(int taxon2 : unmatched) {
               double[][] dist=IBSDistanceMatrix.computeHetBitDistances(inAlign, 
                       taxon1, taxon2, minSitesForSectionComp, startBlock, endBlock, sections, badMask);
               double maxDist=0;
               for (double[] ds : dist) {
                   if(Double.isNaN(ds[0])) {maxDist=1.0;}
                   else if(ds[0]>maxDist) maxDist=ds[0];  
               }
               if(maxDist<maxDistance) {
                   hits.add(taxon2);
               }
            }
            byte[] calls=null;
            //System.out.println("Unk"+countUnknown(calls));
            if((hits.size()+1)<this.minTaxaInGroup) continue;
            if(hits.size()>0) {
                ArrayList<String> mergeNames=new ArrayList<String>();
                mergeNames.add(inIDG.getIdentifier(taxon1).getFullName());
                mergeSets.put(taxon1, hits);         
//                System.out.printf("Taxa1 %d %s %n",taxon1,hits.toString());
                System.out.print(inAlign.getFullTaxaName(taxon1)+"=");
                for (Integer taxon2 : hits) {
                    unmatched.remove(taxon2);
                    System.out.print(inAlign.getFullTaxaName(taxon2)+"=");
                    mergeNames.add(inIDG.getIdentifier(taxon2).getFullName());
                }
                System.out.println("");              
                calls=consensusGameteCalls(inAlign, mergeNames, startSite, endSite, 0.05);
            } else {
                calls=inAlign.getBaseRange(taxon1, startSite, endSite+1);
            }
            int[] unkCnt=countUnknown(calls);
            double missingFreq=(double)unkCnt[0]/(double)inAlign.getSiteCount();
            double hetFreq=(double)unkCnt[1]/(double)(inAlign.getSiteCount()-unkCnt[0]);
            if((missingFreq<minimumMissing)&&(hetFreq<0.01)) {
                int index=(hits.size()*200000)+taxon1;
                System.out.printf("Output %s plus %d missingF:%g hetF:%g index: %d %n",inIDG.getIdentifier(taxon1).getFullName(),
                        hits.size(), missingFreq, hetFreq, index);   
                results.put(index, calls);
            }
            if(results.size()>=maxHaplotypes) break;
        }
        return results;
    }
    
    private byte[] consensusGameteCalls(Alignment a, List<String> taxa, int startSite, 
            int endSite, double maxError) {
        int[] taxaIndex = new int[taxa.size()];
        for (int t = 0; t < taxaIndex.length; t++) {  //why are we working with names rather than numbers
            taxaIndex[t] = a.getIdGroup().whichIdNumber(taxa.get(t));
        }
        byte[] calls = new byte[endSite-startSite+1];
        Arrays.fill(calls, Alignment.UNKNOWN_DIPLOID_ALLELE);
        for (int s = startSite; s <= endSite; s++) {
            byte mjAllele = a.getMajorAllele(s);
            byte mnAllele = a.getMinorAllele(s);
            byte mj=AlignmentUtils.getDiploidValue(mjAllele,mjAllele);
            byte mn=AlignmentUtils.getDiploidValue(mnAllele,mnAllele);
            byte het = AlignmentUtils.getDiploidValue(mjAllele, mnAllele);
            int mjCnt=0, mnCnt=0;
            for (int t = 0; t < taxaIndex.length; t++) {
                byte ob = a.getBase(taxaIndex[t], s);
                if (ob == Alignment.UNKNOWN_DIPLOID_ALLELE) {
                    continue;
                }
                if (ob == mj) {
                    mjCnt++;
                } else if (ob == mn) {
                    mnCnt++;
                } else if (ob == het) {
                    mjCnt++;
                    mnCnt++;
                }
            }
            int totalCnt = mjCnt + mnCnt;
            
            if (totalCnt == 0) {
                double missingProp=1.0;
                for (int t : taxaIndex) {missingProp*=propMissing[t];}
                if(callGaps&(missingProp<minJointGapProb)) calls[s-startSite]=NucleotideAlignmentConstants.GAP_DIPLOID_ALLELE;
                continue;
            }
            double errorRate;
            if(totalCnt>1) siteCallCnt[s]+=totalCnt;
            if(mjCnt < mnCnt) {
                errorRate=(double)mjCnt/(double)totalCnt;
                if(errorRate<maxError) {calls[s-startSite] = mn;}
                else {siteErrors[s]+=mjCnt;}
            } else {
                errorRate=(double)mnCnt/(double)totalCnt;
                if(errorRate<maxError) {calls[s-startSite] = mj;}
                else {siteErrors[s]+=mnCnt;}
            }
        }
        return calls;
    }
    
    public static ArrayList<Integer> maxMajorAllelesTaxa(Alignment a, int numMaxTaxa, int alleleNumber) {
        ArrayList<Integer> maxTaxa=new ArrayList<Integer>();
        OpenBitSet curMj=new OpenBitSet(a.getSiteCount());
        long maxMjCnt=curMj.cardinality();
        for (int i = 0; i < numMaxTaxa; i++) {
            long bestCnt=0;
            int bestAddTaxa=-1;
            for (int t = 0; t < a.getSequenceCount(); t++) {
                BitSet testMj=new OpenBitSet(a.getAllelePresenceForAllSites(t, alleleNumber));
                testMj.union(curMj);
                long cnt=testMj.cardinality();
                if(cnt>bestCnt) {
                    bestCnt=cnt;
                    bestAddTaxa=t;
                }
            }
            if(maxMjCnt==bestCnt) continue;
            curMj.union(a.getAllelePresenceForAllSites(bestAddTaxa, alleleNumber));
            maxMjCnt=curMj.cardinality();
            maxTaxa.add(bestAddTaxa);
            System.out.printf("Allele:%d Taxa: %s %d %n",alleleNumber,a.getTaxaName(bestAddTaxa),maxMjCnt);
        }
        return maxTaxa;
    }
    
    private int[] countUnknown(byte[] b) {
        int cnt=0, cntHet=0;
        for (int i = 0; i < b.length; i++) {
            if(b[i]==Alignment.UNKNOWN_DIPLOID_ALLELE) {cnt++;}
            else if(AlignmentUtils.isHeterozygous(b[i])) {cntHet++;}
        }
        int[] result={cnt,cntHet};
        return result;
    }
    
    
    
    
    
    /**
     *
     * @param args
     */
    public static void main(String[] args) {
//      String root="/Users/edbuckler/SolexaAnal/GBS/build20120110/imp/";
        String root="/Users/edbuckler/SolexaAnal/GBS/build20120701/impOrig/";
//        String root="/Users/edbuckler/SolexaAnal/GBS/build20120701/06_HapMap/";

        String infile=root+"Z0NE00N_chr10.hmp.txt.gz";
        String infile12228=root+"all25f12288.c10.hmp.txt.gz";
        String infile24K=root+"04_PivotMergedTaxaTBT.c10_s0_s24575_masked.hmp.txt.gz";
        String infile4097=root+"all25.c10_s4096_s8191_masked.hmp.txt.gz";
        String infileCN=root+"CNNAM_Filt_2_6.c10.hmp.txt.gz";
       // String infileS=root+"Z0NE00N_chr10S.hmp.txt.gz";
        String infileS=root+"All_chr10S.hmp.txt.gz";
        String infileL=root+"All_chr10L.hmp.txt.gz";
 //       String donorFile=root+"AllTaxa_BPEC_AllZea_GBS_Build_July_2012_FINAL_chr10.hmp.txt.gz";
 //       String donorFile=root+"DTMAfounder20120110.imp.hmp.txt";
 //       String mergeFile=root+"maskedMerge20130429.hmp.txt.gz";
        String mergeFile=root+"w24575Of24KMerge20130513.hmp.txt.gz";
        String errorFile=root+"mcErrorXMerge20130425.txt";
        String errorFile2=root+"mcErrorXMerge20130425.txt";
        errorFile=errorFile2=null;
 //       String mergeFile=root+"AllHFMerge_chr10L.hmp.txt.gz";

        System.out.println("Resolve Method 0");
        MergeIdenticalGametes e64NNI;
//        Alignment a=ImportUtils.readFromHapmap(infileCN, false, (ProgressListener)null);
//        a.optimizeForTaxa(null);
//        ArrayList<Identifier> al=new ArrayList<Identifier>();
//        int sites=a.getSiteCount();
//        for (int i = 0; i < a.getSequenceCount(); i++) {
//            int totalSitesNotMissing = a.getTotalNotMissingForTaxon(i);
//            double propMissing=(double)(1+sites-totalSitesNotMissing)/(double)sites; //1 prevents error in joint missing calculations
//            if(propMissing<0.8) al.add(a.getIdGroup().getIdentifier(i));
//        }
//        SimpleIdGroup ids=new SimpleIdGroup(al.size());
//        for (int i = 0; i < al.size(); i++) {
//            ids.setIdentifier(i, al.get(i));
//        }
//        Alignment fa=FilterAlignment.getInstance(a, ids,false);
//    //    fa.optimizeForSites(null);
//        ExportUtils.writeToHapmap(fa, false, root+"hcovCN10Merge20130505.hmp.txt.gz", '\t', null);
//        System.exit(0);
        
      //  e64NNI=new MergeIdenticalGametes(infile, mergeFile,0.02,1000);
//        e64NNI=new MergeIdenticalGametes(infile12228, mergeFile,0.02,1000);
//        e64NNI=new MergeIdenticalGametes(infile4097, mergeFile, errorFile, errorFile2, 0.01,500);
        e64NNI=new MergeIdenticalGametes(infile24K, mergeFile, errorFile, errorFile2, 0.01,500);
//        e64NNI=new MergeIdenticalGametes(infile4097, mergeFile, errorFile, errorFile2, 0.01,100);
//        e64NNI=new MergeIdenticalGametes(infileCN, mergeFile, errorFile, errorFile2, 0.01,500);
//        e64NNI=new MergeIdenticalGametes(infileS, mergeFile,0.02,500);
//        e64NNI=new MergeIdenticalGametes(infileL, mergeFile,0.02,1000);

    }
    
}
