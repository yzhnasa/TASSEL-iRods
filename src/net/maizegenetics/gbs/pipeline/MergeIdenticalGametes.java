/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.gbs.pipeline;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Map.Entry;
import java.util.TreeMap;
import java.util.TreeSet;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.AlignmentUtils;
import net.maizegenetics.pal.alignment.BitAlignment;
import net.maizegenetics.pal.alignment.ExportUtils;
import net.maizegenetics.pal.alignment.ImportUtils;
import net.maizegenetics.pal.alignment.MutableNucleotideAlignment;
import net.maizegenetics.pal.distance.IBSDistanceMatrix;
import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.Identifier;
import net.maizegenetics.pal.ids.SimpleIdGroup;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;
import net.maizegenetics.util.ProgressListener;

/**
 *
 * @author edbuckler
 */
public class MergeIdenticalGametes {
    private Alignment inAlign;
    private Alignment donorAlign;

    public MergeIdenticalGametes(String inFile, String exportFile, double maxDistance, int minSites) {        
        inAlign=ImportUtils.readFromHapmap(inFile, false, (ProgressListener)null);
        inAlign.optimizeForTaxa(null);
        IdGroup inIDG=inAlign.getIdGroup();
        System.out.printf("In taxa:%d sites:%d %n",inAlign.getSequenceCount(),inAlign.getSiteCount());
        //Develop lists of the modest coverage, homozygous taxa
        TreeMap<Integer,Integer> presentRanking=new TreeMap<Integer,Integer>(Collections.reverseOrder());
        TreeSet<Integer> unmatched=new TreeSet<Integer>();
        int weirdB73=-1;
        for (int i = 0; i < inAlign.getSequenceCount(); i++) {
            int totalSitesNotMissing = inAlign.getTotalNotMissingForTaxon(i);
            double hetFreq=(double)inAlign.getHeterozygousCountForTaxon(i)/(double)totalSitesNotMissing;
            //System.out.printf("%s %d %g %n",inAlign.getTaxaName(i),totalSitesNotMissing,hetFreq);
            if((hetFreq>0.006)||(totalSitesNotMissing<2000)) continue;
            presentRanking.put(totalSitesNotMissing, i);
            unmatched.add(i);
            if(inAlign.getFullTaxaName(i).contains("B73:81FE7ABXX:3:250031244")) {
                System.out.println("B73:81FE7ABXX:3:250031244 i:"+i);
                weirdB73=i;
            }
        }
        System.out.println("Inbred and modest coverage:"+unmatched.size());
        TreeMap<Integer,ArrayList> mergeSets=new TreeMap<Integer,ArrayList>();
        TreeMap<Identifier,byte[]> results=new TreeMap<Identifier,byte[]>();
        for (Entry<Integer,Integer> e : presentRanking.entrySet()) {
            //System.out.println(e.toString());
           // System.out.println("B73:81FE7ABXX:3:250031244 in list:"+unmatched.contains(weirdB73));
            int taxon1=e.getValue();
            if(unmatched.contains(taxon1)==false) continue;//already included in another group
            ArrayList<Integer> hits=new ArrayList<Integer>();
            unmatched.remove(taxon1);
            for(int taxon2 : unmatched) {
               double dist=IBSDistanceMatrix.computeHetBitDistances(inAlign, taxon1, taxon2, minSites, false)[0];
               if(inAlign.getTaxaName(taxon1).contains("REF")) {
                   if(inAlign.getTaxaName(taxon2).contains("B73")) {
                       System.out.printf("Dist %s %s %g %n",inAlign.getFullTaxaName(taxon1), inAlign.getFullTaxaName(taxon2),dist);
                   }
               }
               if(dist<maxDistance) {
                   hits.add(taxon2);
               }
            }
            byte[] calls=inAlign.getBaseRow(taxon1);
            //System.out.println("Unk"+countUnknown(calls));
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
                calls=MergeIdenticalTaxaPlugin.consensusCalls(inAlign, mergeNames, false, 0.8);
               // System.out.println("UnkACalls:"+countUnknown(calls));
            }
            int[] unkCnt=countUnknown(calls);
            double missingFreq=(double)unkCnt[0]/(double)inAlign.getSiteCount();
            double hetFreq=(double)unkCnt[1]/(double)(inAlign.getSiteCount()-unkCnt[0]);
            if((missingFreq<0.5)&&(hetFreq<0.01)) {
                System.out.printf("Output %s plus %d missingF:%g hetF:%g %n",inIDG.getIdentifier(taxon1).getFullName(),
                        hits.size(), missingFreq, hetFreq);
                results.put(new Identifier(inIDG.getIdentifier(taxon1).getFullName()+":P"+hits.size()), calls);
            }
            //System.out.println("Unk"+countUnknown(calls));
        }
     //   IdGroup outIDG=new SimpleIdGroup(results.size());
        
        IdGroup outIDG=new SimpleIdGroup((Identifier[])results.keySet().toArray(new Identifier[0]));
        MutableNucleotideAlignment mna=MutableNucleotideAlignment.getInstance(outIDG, inAlign.getSiteCount());
//        mna.getLoci()
        for (int i = 0; i < inAlign.getSiteCount(); i++) {
            mna.addSite(i);
            mna.setLocusOfSite(i, inAlign.getLocus(i));
            mna.setPositionOfSite(i, inAlign.getPositionInLocus(i));
        }
        for (int i = 0; i < outIDG.getIdCount(); i++) {
            mna.setBaseRange(i, 0, results.get(outIDG.getIdentifier(i)));
        }
        mna.clean();
        ExportUtils.writeToHapmap(mna, false, exportFile, '\t', null);
        Alignment ba=ImportUtils.readFromHapmap(inFile, false, (ProgressListener)null);
        System.out.printf("Mrg taxa:%d sites:%d %n",ba.getSequenceCount(),ba.getSiteCount());
       // Alignment ba=BitAlignment.getInstance(mna, 2, false, false);
        maxMajorAllelesTaxa(ba,50,0);
        System.out.println();
        maxMajorAllelesTaxa(ba,500,1);
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
//        String root="/Volumes/LaCie/build20120110/imp/";
        String root="/Users/edbuckler/SolexaAnal/GBS/build20120701/06_HapMap/";

        String infile=root+"Z0NE00N_chr10.hmp.txt.gz";
        String infile12228=root+"all25f12288.c10.hmp.txt.gz";
       // String infileS=root+"Z0NE00N_chr10S.hmp.txt.gz";
        String infileS=root+"All_chr10S.hmp.txt.gz";
        String infileL=root+"All_chr10L.hmp.txt.gz";
 //       String donorFile=root+"AllTaxa_BPEC_AllZea_GBS_Build_July_2012_FINAL_chr10.hmp.txt.gz";
 //       String donorFile=root+"DTMAfounder20120110.imp.hmp.txt";
        String mergeFile=root+"Merge20130421.hmp.txt";
 //       String mergeFile=root+"AllHFMerge_chr10L.hmp.txt.gz";

        System.out.println("Resolve Method 0");
        MergeIdenticalGametes e64NNI;
      //  e64NNI=new MergeIdenticalGametes(infile, mergeFile,0.02,1000);
        e64NNI=new MergeIdenticalGametes(infile12228, mergeFile,0.02,1000);
//        e64NNI=new MergeIdenticalGametes(infileS, mergeFile,0.02,500);
//        e64NNI=new MergeIdenticalGametes(infileL, mergeFile,0.02,1000);

    }
    
}
