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
import net.maizegenetics.pal.alignment.BitAlignment;
import net.maizegenetics.pal.alignment.BitAlignmentHDF5;
import net.maizegenetics.pal.alignment.ExportUtils;
import net.maizegenetics.pal.alignment.FilterAlignment;
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
 * TODO:
 * 1.  Add names back in
 * 2.  Cluster and choose haplotypes by cluster and number of new minor alleles (or information)
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
    private int minSitesForSectionComp=50;
    private double maxHetFreq=0.02;
    private double maxErrorInCreatingConsensus=0.05;
    private int minTaxaInGroup=2;
    private String hapMapNames="CrappyMap";

    public MergeIdenticalGametes(String inFile, String exportFile, String errorInFile,
            String errorExportFile, double maxDistance, int minSites, int appoxSitesPerHaplotype, String hapMapNames) {
        this.hapMapNames=hapMapNames;
        Alignment baseAlign=null;
        System.out.println("Reading: "+inFile);
        if(inFile.contains(".h5")) {
            baseAlign=BitAlignmentHDF5.getInstance(inFile);
        } else {
            baseAlign=ImportUtils.readFromHapmap(inFile, false, (ProgressListener)null);
        }
        //inAlign.optimizeForTaxa(null);
        int[][] divisions=divideChromosome(baseAlign, appoxSitesPerHaplotype);
       // System.exit(0);
        inAlign=baseAlign;
        
        
        System.out.printf("In taxa:%d sites:%d %n",inAlign.getSequenceCount(),inAlign.getSiteCount());
        
        for (int i = 0; i < divisions.length; i++) {
            MutableNucleotideAlignment mna=createHaplotypeAlignment(divisions[i][0], divisions[i][1], baseAlign,
             minSites,  maxDistance);
            String newExport=exportFile.replace("s+.hmp", "s"+i+".hmp");
            ExportUtils.writeToHapmap(mna, false, newExport, '\t', null);
  //          if(errorExportFile!=null) exportBadSites(errorExportFile, 0.01);   
        }
        
//        for (int i = 0; i < siteErrors.length; i++) {
//            System.out.printf("%d %d %d %g %g %n",i,siteCallCnt[i],siteErrors[i], 
//                    (double)siteErrors[i]/(double)siteCallCnt[i], inAlign.getMinorAlleleFrequency(i));
//        }
    }
    
    private MutableNucleotideAlignment createHaplotypeAlignment(int startSite, int endSite, Alignment baseAlign,
            int minSites, double maxDistance) {
        Alignment fa=FilterAlignment.getInstance(baseAlign, startSite, endSite);
        inAlign=BitAlignment.getInstance(fa, true);
        inAlign.optimizeForTaxa(null);
        int sites=inAlign.getSiteCount();
        System.out.printf("SubInAlign taxa:%d sites:%d %n",inAlign.getSequenceCount(),inAlign.getSiteCount());
        siteErrors=new int[inAlign.getSiteCount()];
        siteCallCnt=new int[inAlign.getSiteCount()];
        propMissing=new double[inAlign.getSequenceCount()];
        int startBlock=0;
        int endBlock=inAlign.getAllelePresenceForAllSites(0, 0).getNumWords();
        int sections=endBlock;  //need to remove sections
        TreeMap<Integer,Integer> presentRanking=createPresentRankingForWindow(startBlock, endBlock, minSites, maxHetFreq);
        System.out.printf("Block %d Inbred and modest coverage:%d %n",startBlock,presentRanking.size());
        System.out.printf("Current Site %d Current block %d EndBlock: %d %n",startSite, startBlock, endBlock);
        TreeMap<Integer,byte[]> results=mergeWithinWindow(presentRanking, startBlock, endBlock, sections, maxDistance);
        MutableNucleotideAlignment mna=createEmptyHaplotypeAlignment(results.size());
        int index=0;
        for (byte[] calls : results.values()) {
            mna.setBaseRange(index, 0, calls);
            index++;
        }
        mna.clean();
        return mna;
    }
    
    private int[][] divideChromosome(Alignment a, int appoxSitesPerHaplotype) {
        int subAlignCnt=(int)Math.round((double)a.getSiteCount()/(double)appoxSitesPerHaplotype);
        int prefBlocks=(a.getSiteCount()/(subAlignCnt*64));
        int[][] divisions=new int[subAlignCnt][2];
        int cnt=0;
        for (int i = 0; i < subAlignCnt; i++) {
            divisions[i][0]=i*prefBlocks*64;
            divisions[i][1]=divisions[i][0]+(prefBlocks*64)-1;
        }
        divisions[subAlignCnt-1][1]=a.getSiteCount()-1;
        System.out.printf("Alignment Sites:%d ApproxSites:%d RealSites:%d %n",a.getSiteCount(),appoxSitesPerHaplotype, prefBlocks*64);
        System.out.println("Chromosome Divisions:"+Arrays.deepToString(divisions));
        return divisions;
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
            if(inAlign.getFullTaxaName(i).contains(hapMapNames)) propPresent+=1.0;  //Puts HapMap at the top of the sort
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
            boolean isHapMap=inAlign.getFullTaxaName(taxon1).contains(hapMapNames);
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
 //              if(inAlign.getFullTaxaName(taxon1).contains("BKN009")) System.out.printf("BKN009 %s %g %n",inAlign.getFullTaxaName(taxon2), maxDist);
            }
            byte[] calls=null;
            //System.out.println("Unk"+countUnknown(calls));
            if(((hits.size()+1)<this.minTaxaInGroup)&&(!isHapMap)) continue;
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
                calls=consensusGameteCalls(inAlign, mergeNames, startSite, endSite, maxErrorInCreatingConsensus);
            } else {
                calls=inAlign.getBaseRange(taxon1, startSite, endSite+1);
            }
            int[] unkCnt=countUnknown(calls);
            double missingFreq=(double)unkCnt[0]/(double)inAlign.getSiteCount();
            double hetFreq=(double)unkCnt[1]/(double)(inAlign.getSiteCount()-unkCnt[0]);
            if(isHapMap||((missingFreq<minimumMissing)&&(hetFreq<0.01))) {
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
       // String root="/Users/edbuckler/SolexaAnal/GBS/build20120701/impOrig/";
        String root="/Volumes/LaCie/build20120701/IMP26/orig/";

  //      String root="/Volumes/LaCie/build20120701/impOrig/";

   //     String infile24K=root+"04_PivotMergedTaxaTBT.c10_s0_s24575_masked.hmp.txt.gz";
  //      String infile24K=root+"all26HM2.c10_s0_24575.hmp.txt.gz";
  //      String infileH5=root+"all26HM2.c10_s0_24575.hmp.h5";
        String infileH5=root+"all26HM2.hmp.h5";
        String infile24K=root+"JustHMw24575OfHM224KMerge20130517b.hmp.txt.gz";
        String infile26v=root+"AllZeaGBS_v2.6_MERGEDUPSNPS_20130513_chr+.hmp.txt.gz";
        String secFile26v=root+"all26_8k.c+s+.hmp.txt.gz";

        String mergeFile=root+"all26HM2_8k.c10.hmp.txt.gz";
        String errorFile=root+"mcErrorXMerge20130425.txt";
        String errorFile2=root+"mcErrorXMerge20130425.txt";
        errorFile=errorFile2=null;

        System.out.println("Resolve Method 0");
        MergeIdenticalGametes e64NNI;


//        e64NNI=new MergeIdenticalGametes(infile24K, mergeFile, errorFile, errorFile2, 0.01,500,"CappyMap");
        for (int chr = 9; chr <= 9; chr++) {
            String infile26vMod=infile26v.replace("chr+.hmp", "chr"+chr+".hmp");
            String outfiles=secFile26v.replace(".c+s+.", ".c"+chr+"s+.");
            e64NNI=new MergeIdenticalGametes(infile26vMod, outfiles, errorFile, errorFile2, 0.01,500,8192, "CappyMap");
            
        }
       // e64NNI=new MergeIdenticalGametes(infileH5, mergeFile, errorFile, errorFile2, 0.01,500,8192, "CappyMap");


    }
    
}
