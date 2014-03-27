/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.analysis.imputation;

import java.awt.Frame;
import net.maizegenetics.dna.snp.*;
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTableBuilder;

import java.io.BufferedOutputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import javax.swing.ImageIcon;
import net.maizegenetics.dna.map.Chromosome;

import static net.maizegenetics.dna.snp.GenotypeTable.UNKNOWN_DIPLOID_ALLELE;
import static net.maizegenetics.dna.snp.GenotypeTableUtils.isHeterozygous;
import static net.maizegenetics.dna.snp.NucleotideAlignmentConstants.GAP_DIPLOID_ALLELE;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.util.ArgsEngine;
import org.apache.commons.lang.ArrayUtils;
import org.apache.log4j.Logger;

/**
 *
 * @author kls283
 */
public class FILLINImputationAccuracyPlugin extends AbstractPlugin {
    public int[] MAF= null;
    private double[][] all= new double[3][5]; //arrays held ("columns"): 0-maskedMinor, 1-maskedHet, 2-maskedMajor; each array ("rows"):0-to minor, 1-to het, 2-to major, 3-unimp, 4-total for known type
    private double[][][] mafAll= null;//sam as all, but first array holds MAF category
    private String maskKeyFile;
    private String unimpFile;
    private String donorFile;
    private String outFile;
    private double propSitesMask;
    private GenotypeTable maskKey= null;
    private double[] MAFClass= null;//new double[]{0,.02,.05,.10,.20,.3,.4,.5,1};
    private boolean verboseOutput= true;
    String[] FILLINargs;
    
    private static ArgsEngine engine = new ArgsEngine();
    private static final Logger myLogger = Logger.getLogger(FILLINImputationPlugin.class);
    
    public FILLINImputationAccuracyPlugin(Frame parentFrame) {
       super(parentFrame, false);
   }
    
    public void initiateAccuracy(String unimpFile, String maskKeyFile, double[] mafClass, double propSitesMask) {
        GenotypeTable unimp= ImportUtils.readGuessFormat(unimpFile);
//        if (mafClass!=null) {//if mafClass not null, assign MAF categories to sites in unimputed
//            generateMAF(donors, unimp, mafClass);
//            mafAll= new double[mafClass.length][3][5];
//            if (verboseOutput) System.out.println("Calculating accuracy within supplied MAF categories.");
//        }
        //mask file if not already masked (depth or proportion) and generate key
        if (maskKeyFile!=null) {
            if (verboseOutput) System.out.println("File already masked. Use input key file for calculating accuracy");
            maskKey= ImportUtils.readGuessFormat(maskKeyFile);
            if (Arrays.equals(maskKey.physicalPositions(), unimp.physicalPositions())==false) {
                maskKey= filterKey(maskKey, unimp);
            }
            if (maskKey==null) {
                FILLINImputationPlugin.unimpAlign= maskPropSites(unimp,propSitesMask);
            }
        }
        //else if (unimp.depth()!=null) FILLINImputationPlugin.unimpAlign= maskFileByDepth(unimp,7, 7);
        else FILLINImputationPlugin.unimpAlign= maskPropSites(unimp,propSitesMask);
    }
    
        public void generateMAF(GenotypeTable[] donorAlign, GenotypeTable unimpAlign, double[] mafClass) {
        MAF= new int[unimpAlign.numberOfSites()];
        for (GenotypeTable don:donorAlign) {
            for (int site = 0; site < don.numberOfSites(); site++) {
                int unimpSite= unimpAlign.positions().indexOf(don.positions().get(site));
                if (unimpSite < 0) {MAF[unimpSite]= -1; continue;}
                int search= Arrays.binarySearch(mafClass, don.minorAlleleFrequency(site));
                MAF[unimpSite]= search<0?Math.abs(search)-1:search;
            }
        }
    }

    //for depths 4 or more, requires hets to be called by more than one for less depth allele. returns 2-array whre index 0 is mask and 1 is key
    private GenotypeTable maskFileByDepth(GenotypeTable a, int depthToMask, int maskDenom) {
        if (verboseOutput) System.out.println("Masking file using depth\nSite depth to mask: "+depthToMask+"Divisor for physical positions to be masked: "+maskDenom);
        GenotypeCallTableBuilder mask= GenotypeCallTableBuilder.getInstance(a.numberOfTaxa(), a.numberOfSites());
        GenotypeCallTableBuilder key= GenotypeCallTableBuilder.getInstance(a.numberOfTaxa(), a.numberOfSites());
        
        int cnt= 0;
        for (int taxon = 0; taxon < a.numberOfTaxa(); taxon++) {
            int taxaCnt= 0;
            mask.setBaseRangeForTaxon(taxon, 0, a.genotypeAllSites(taxon));
            for (int site = 0; site < a.numberOfSites(); site++) {
                if (GenotypeTableUtils.isEqual(NucleotideGenotypeTable.UNKNOWN_DIPLOID_ALLELE, a.genotype(taxon, site))) continue;
                if (a.physicalPositions()[site]%maskDenom!=0) continue;
                int[] currD= a.depthForAlleles(taxon, site);
                if (currD[0]+currD[1]!=depthToMask) continue;
                else if ((a.isHeterozygous(taxon, site)==false) ||
                            (depthToMask > 3 && currD[0] > 1 && currD[1] > 1)|| 
                            (depthToMask < 4)) {
                    mask.setBase(taxon, site, NucleotideGenotypeTable.UNKNOWN_DIPLOID_ALLELE); key.setBase(taxon, site, a.genotype(taxon, site)); taxaCnt++;
                }
            }
            if (verboseOutput) System.out.println(taxaCnt+" sites masked for "+a.taxaName(taxon)); cnt+= taxaCnt;
        }
        if (verboseOutput) System.out.println(cnt+" sites masked at a depth of "+depthToMask+" (site numbers that can be divided by "+maskDenom+")");
        maskKey= GenotypeTableBuilder.getInstance(key.build(), a.positions(), a.taxa());
        return GenotypeTableBuilder.getInstance(mask.build(), a.positions(), a.taxa());
    }
    
    private GenotypeTable maskPropSites(GenotypeTable a, double propSitesMask) {
        if (verboseOutput) System.out.println("Masking file without depth\nMasking "+propSitesMask+" proportion of sites");
        GenotypeCallTableBuilder mask= GenotypeCallTableBuilder.getInstance(a.numberOfTaxa(), a.numberOfSites());
        GenotypeCallTableBuilder key= GenotypeCallTableBuilder.getInstance(a.numberOfTaxa(), a.numberOfSites());
        
        int presGenos= 0;
        for (int taxon = 0; taxon < a.numberOfTaxa(); taxon++) {presGenos+= a.totalNonMissingForTaxon(taxon);}
        int expected= (int)(propSitesMask*(double)presGenos);
        int cnt= 0;
        for (int taxon = 0; taxon < a.numberOfTaxa(); taxon++) {
            int taxaCnt= 0;
            mask.setBaseRangeForTaxon(taxon, 0, a.genotypeAllSites(taxon));
            for (int site = 0; site < a.numberOfSites(); site++) {
                if (Math.random()<propSitesMask && GenotypeTableUtils.isEqual(NucleotideGenotypeTable.UNKNOWN_DIPLOID_ALLELE, a.genotype(taxon, site))==false) {
                    mask.setBase(taxon, site, NucleotideGenotypeTable.UNKNOWN_DIPLOID_ALLELE); key.setBase(taxon, site, a.genotype(taxon, site)); taxaCnt++;
                }
            }
            cnt+= taxaCnt;
        }
        if (verboseOutput) System.out.println(cnt+" sites masked randomly not based on depth ("+expected+" expected at "+propSitesMask+")");
        maskKey= GenotypeTableBuilder.getInstance(key.build(), a.positions(), a.taxa());
        return GenotypeTableBuilder.getInstance(mask.build(), a.positions(), a.taxa());
    }
    
    //filters for site position and chromosome. returns null if mask has fewer chromosomes or positions in chromosomes than unimp
    private GenotypeTable filterKey(GenotypeTable maskKey, GenotypeTable unimp) {
        if (verboseOutput) System.out.println("Filtering user input key file...\nSites in original Key file: "+maskKey.numberOfSites());
        String[] unimpNames= new String[unimp.numberOfSites()];
        for (int site = 0; site < unimp.numberOfSites(); site++) {unimpNames[site]= unimp.siteName(site);}
        int[] unimpPos;
        int[] keyPos;
        ArrayList<String> keepSites= new ArrayList<>();
        if (Arrays.equals(unimp.chromosomes(), maskKey.chromosomes())==false) maskKey= matchChromosomes(unimp, maskKey);
        if (maskKey==null) return null;
        for (Chromosome chr:unimp.chromosomes()) {
            int[] startEndUnimp= unimp.startAndEndOfChromosome(chr); int[] startEndKey= maskKey.startAndEndOfChromosome(chr);
            unimpPos= Arrays.copyOfRange(unimp.physicalPositions(), startEndUnimp[0], startEndUnimp[1]+1);
            keyPos= Arrays.copyOfRange(maskKey.physicalPositions(), startEndKey[0], startEndKey[1]+1);
            for (int posOnChr = 0; posOnChr < unimpPos.length; posOnChr++) {//if input hapmap sites not in key, return null
                if (Arrays.binarySearch(keyPos, unimpPos[posOnChr])<0) return null;
            }
            for (int posOnChr = 0; posOnChr < keyPos.length; posOnChr++) {//if key site in input hapmap, retain
                if (Arrays.binarySearch(unimpPos, keyPos[posOnChr])>-1) {
                    keepSites.add(maskKey.siteName(startEndKey[0]+posOnChr));
                }
            }
        }
        FilterGenotypeTable filter= FilterGenotypeTable.getInstance(maskKey, keepSites.toArray(new String[keepSites.size()]));
        GenotypeTable newMask= GenotypeTableBuilder.getGenotypeCopyInstance(filter);
        if (verboseOutput) System.out.println("Sites in new mask: "+newMask.numberOfSites());
        return newMask;
    }
    
    private GenotypeTable matchChromosomes(GenotypeTable unimp, GenotypeTable maskKey) {
        Chromosome[] unimpChr= unimp.chromosomes();
        Chromosome[] keyChr= maskKey.chromosomes();
        ArrayList<Integer> keepSites= new ArrayList<>();
        for (Chromosome chr:unimpChr) { //if any of the chromosomes in input do not exist in key, return null (which then masks proportion)
            if (Arrays.binarySearch(keyChr, chr)<0) return null;
        }
        for (Chromosome chr:keyChr) { //keep sites on key that are on matching chromosomes
            if (Arrays.binarySearch(unimpChr, chr)>-1) {
                for (int site = maskKey.startAndEndOfChromosome(chr)[0]; site < maskKey.startAndEndOfChromosome(chr)[1]; site++) {
                    keepSites.add(site);
                }
            }
        }
        FilterGenotypeTable filter= FilterGenotypeTable.getInstance(maskKey, ArrayUtils.toPrimitive(keepSites.toArray(new Integer[keepSites.size()])));
        GenotypeTable matchChr= GenotypeTableBuilder.getGenotypeCopyInstance(filter);
        return matchChr;
    }
    
    //this is the sample multiple r2.
    private double pearsonR2(double[][] all) {
        int size= 0;
        for (int x = 0; x < 3; x++) {size+= (all[x][4]-all[x][3]);}
        double[][] xy= new double[2][size]; //0 is x, 1 is y
        int last= 0;//the last index filled
        for (double x = 0; x < 3; x++) { for (double y = 0; y < 3; y++) {
                for (int fill = last; fill < last+all[(int)x][(int)y]; fill++) {
                    xy[0][fill]= x;
                    xy[1][fill]= y;
                }
                last= last+(int)all[(int)x][(int)y];
            }}
        double meanX= 0; double meanY= 0; double varX= 0; double varY= 0; double covXY= 0; double r2= 0.0;
        for (int i = 0; i < xy[0].length; i++) {meanX+=xy[0][i]; meanY+= xy[1][i];}
        meanX= meanX/(xy[0].length-1); meanY= meanY/(xy[1].length-1);
        double currX, currY;
        for (int i = 0; i < xy[0].length; i++) {
            currX= xy[0][i]-meanX; currY= xy[1][i]-meanY;
            varX+= currX*currX; varY+= currY*currY;
            covXY+= currX*currY;
        }
        r2= (covXY/(Math.sqrt(varX)*Math.sqrt(varY)))*(covXY/(Math.sqrt(varX)*Math.sqrt(varY)));
        if (verboseOutput) System.out.println("Unadjusted R2 value for "+size+" comparisons: "+r2);
        return r2;
    }
    
    private void accuracyOut(double[][] all, double time) {
        DecimalFormat df = new DecimalFormat("0.########");
        double r2= pearsonR2(all);
        try {
            File outputFile = new File(outFile.substring(0, outFile.indexOf(".hmp")) + "DepthAccuracy.txt");
            DataOutputStream outStream = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outputFile)));
            outStream.writeBytes("##Taxon\tTotalSitesMasked\tTotalSitesCompared\tTotalPropUnimputed\tNumMinor\tCorrectMinor\tMinorToHet\tMinorToMajor\tUnimpMinor"
                    + "\tNumHets\tHetToMinor\tCorrectHet\tHetToMajor\tUnimpHet\tNumMajor\tMajorToMinor\tMajorToHet\tCorrectMajor\tUnimpMajor\tR2\n");
            outStream.writeBytes("##TotalByImputed\t"+(all[0][4]+all[1][4]+all[2][4])+"\t"+(all[0][4]+all[1][4]+all[2][4]-all[0][3]-all[1][3]-all[2][3])+"\t"+
                    ((all[0][3]+all[1][3]+all[2][3])/(all[0][4]+all[1][4]+all[2][4]))+"\t"+all[0][4]+"\t"+all[0][0]+"\t"+all[0][1]+"\t"+all[0][2]+"\t"+all[0][3]+
                    "\t"+all[1][4]+"\t"+all[1][0]+"\t"+all[1][1]+"\t"+all[1][2]+"\t"+all[1][3]+"\t"+all[2][4]+"\t"+all[2][0]+"\t"+all[2][1]+"\t"+all[2][2]+
                    "\t"+all[2][3]+"\t"+r2+"\n");
            outStream.writeBytes("#Minor=0,Het=1,Major=2;x is masked(known), y is predicted\nx\ty\tN\tprop\n"
                    +0+"\t"+0+"\t"+all[0][0]+"\t"+df.format((all[0][0])/(all[0][0]+all[0][1]+all[0][2]))+"\n"
                    +0+"\t"+.5+"\t"+all[0][1]+"\t"+df.format((all[0][1])/(all[0][0]+all[0][1]+all[0][2]))+"\n"
                    +0+"\t"+1+"\t"+all[0][2]+"\t"+df.format((all[0][2])/(all[0][0]+all[0][1]+all[0][2]))+"\n"
                    +.5+"\t"+0+"\t"+all[1][0]+"\t"+df.format((all[1][0])/(all[1][0]+all[1][1]+all[1][2]))+"\n"
                    +.5+"\t"+.5+"\t"+all[1][1]+"\t"+df.format((all[1][1])/(all[1][0]+all[1][1]+all[1][2]))+"\n"
                    +.5+"\t"+1+"\t"+all[1][2]+"\t"+df.format((all[1][2])/(all[1][0]+all[1][1]+all[1][2]))+"\n"
                    +1+"\t"+0+"\t"+all[2][0]+"\t"+df.format((all[2][0])/(all[2][0]+all[2][1]+all[2][2]))+"\n"
                    +1+"\t"+.5+"\t"+all[2][1]+"\t"+df.format((all[2][1])/(all[2][0]+all[2][1]+all[2][2]))+"\n"
                    +1+"\t"+1+"\t"+all[2][2]+"\t"+df.format((all[2][2])/(all[2][0]+all[2][1]+all[2][2]))+"\n");
            outStream.writeBytes("#Proportion unimputed:\n#minor <- "+all[0][3]/all[0][4]+"\n#het<- "+all[1][3]/all[1][4]+"\n#major<- "+all[2][3]/all[2][4]+"\n");
            outStream.writeBytes("#Time to impute and calculate accuracy: "+time+" seconds");
            if (verboseOutput) System.out.println("##Taxon\tTotalSitesMasked\tTotalSitesCompared\tTotalPropUnimputed\tNumMinor\tCorrectMinor\tMinorToHet\tMinorToMajor\tUnimpMinor"
                    + "\tNumHets\tHetToMinor\tCorrectHet\tHetToMajor\tUnimpHet\tNumMajor\tMajorToMinor\tMajorToHet\tCorrectMajor\tUnimpMajor\tR2");
            if (verboseOutput) System.out.println("TotalByImputed\t"+(all[0][4]+all[1][4]+all[2][4])+"\t"+(all[0][4]+all[1][4]+all[2][4]-all[0][3]-all[1][3]-all[2][3])+"\t"+
                    ((all[0][3]+all[1][3]+all[2][3])/(all[0][4]+all[1][4]+all[2][4]))+"\t"+all[0][4]+"\t"+all[0][0]+"\t"+all[0][1]+"\t"+all[0][2]+"\t"+all[0][3]+
                    "\t"+all[1][4]+"\t"+all[1][0]+"\t"+all[1][1]+"\t"+all[1][2]+"\t"+all[1][3]+"\t"+all[2][4]+"\t"+all[2][0]+"\t"+all[2][1]+"\t"+all[2][2]+
                    "\t"+all[2][3]+"\t"+r2);
            if (verboseOutput) System.out.println("Proportion unimputed:\nminor: "+all[0][3]/all[0][4]+"\nhet: "+all[1][3]/all[1][4]+"\nmajor: "+all[2][3]/all[2][4]);
            if (verboseOutput) System.out.println("#Minor=0,Het=1,Major=2;x is masked(known), y is predicted\nx\ty\tN\tprop\n"
                    +0+"\t"+0+"\t"+all[0][0]+"\t"+(all[0][0])/(all[0][0]+all[0][1]+all[0][2])+"\n"
                    +0+"\t"+.5+"\t"+all[0][1]+"\t"+(all[0][1])/(all[0][0]+all[0][1]+all[0][2])+"\n"
                    +0+"\t"+1+"\t"+all[0][2]+"\t"+(all[0][2])/(all[0][0]+all[0][1]+all[0][2])+"\n"
                    +.5+"\t"+0+"\t"+all[1][0]+"\t"+(all[1][0])/(all[1][0]+all[1][1]+all[1][2])+"\n"
                    +.5+"\t"+.5+"\t"+all[1][1]+"\t"+(all[1][1])/(all[1][0]+all[1][1]+all[1][2])+"\n"
                    +.5+"\t"+1+"\t"+all[1][2]+"\t"+(all[1][2])/(all[1][0]+all[1][1]+all[1][2])+"\n"
                    +1+"\t"+0+"\t"+all[2][0]+"\t"+(all[2][0])/(all[2][0]+all[2][1]+all[2][2])+"\n"
                    +1+"\t"+.5+"\t"+all[2][1]+"\t"+(all[2][1])/(all[2][0]+all[2][1]+all[2][2])+"\n"
                    +1+"\t"+1+"\t"+all[2][2]+"\t"+(all[2][2])/(all[2][0]+all[2][1]+all[2][2])+"\n");
            outStream.close();
        } catch (Exception e) {
            if (verboseOutput) System.out.println(e);
        }
    }
    
    private void accuracyMAFOut(double[][][] mafAll) {
        DecimalFormat df = new DecimalFormat("0.########");
        if (MAF!=null && MAFClass!=null) try {
            File outputFile = new File(outFile.substring(0, outFile.indexOf(".hmp")) + "DepthAccuracyMAF.txt");
            DataOutputStream outStream = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outputFile)));
            outStream.writeBytes("##\tMAFClass\tTotalSitesMasked\tTotalSitesCompared\tTotalPropUnimputed\tNumHets\tHetToMinor\tHetToMajor\tCorrectHet\tUnimpHet\tNumMinor\tMinorToMajor\tMinorToHet\tCorrectMinor\t"
                    + "UnimpMinor\tNumMajor\tMajorToMinor\tMajorToHet\tCorrectMajor\tUnimputedMajor\tr2\n");
            for (int i= 0; i<MAFClass.length;i++) {
                outStream.writeBytes("##TotalByImputed\t"+MAFClass[i]+"\t"+(mafAll[i][0][4]+mafAll[i][1][4]+mafAll[i][2][4])+"\t"+(mafAll[i][0][4]+mafAll[i][1][4]+mafAll[i][2][4]-mafAll[i][0][3]-mafAll[i][1][3]-mafAll[i][2][3])+"\t"+
                    ((mafAll[i][0][3]+mafAll[i][1][3]+mafAll[i][2][3])/(mafAll[i][0][4]+mafAll[i][1][4]+mafAll[i][2][4]))+"\t"+mafAll[i][0][4]+"\t"+mafAll[i][0][0]+"\t"+mafAll[i][0][1]+"\t"+mafAll[i][0][2]+"\t"+mafAll[i][0][3]+
                    "\t"+mafAll[i][1][4]+"\t"+mafAll[i][1][0]+"\t"+mafAll[i][1][1]+"\t"+mafAll[i][1][2]+"\t"+mafAll[i][1][3]+"\t"+mafAll[i][2][4]+"\t"+mafAll[i][2][0]+"\t"+mafAll[i][2][1]+"\t"+mafAll[i][2][2]+
                    "\t"+mafAll[i][2][3]+"\t"+pearsonR2(mafAll[i])+"\n");
            }
            outStream.writeBytes("#MAFClass,Minor=0,Het=1,Major=2;x is masked(known), y is predicted\nMAF\tx\ty\tN\tprop\n");
            for (int i= 0; i<MAFClass.length;i++) { outStream.writeBytes(
                    MAFClass[i]+"\t"+0+"\t"+0+"\t"+mafAll[i][0][0]+"\t"+df.format((mafAll[i][0][0])/(mafAll[i][0][0]+mafAll[i][0][1]+mafAll[i][0][2]))+"\n"
                    +MAFClass[i]+"\t"+0+"\t"+.5+"\t"+mafAll[i][0][1]+"\t"+df.format((mafAll[i][0][1])/(mafAll[i][0][0]+mafAll[i][0][1]+mafAll[i][0][2]))+"\n"
                    +MAFClass[i]+"\t"+0+"\t"+1+"\t"+mafAll[i][0][2]+"\t"+df.format((mafAll[i][0][2])/(mafAll[i][0][0]+mafAll[i][0][1]+mafAll[i][0][2]))+"\n"
                    +MAFClass[i]+"\t"+.5+"\t"+0+"\t"+mafAll[i][1][0]+"\t"+df.format((mafAll[i][1][0])/(mafAll[i][1][0]+mafAll[i][1][1]+mafAll[i][1][2]))+"\n"
                    +MAFClass[i]+"\t"+.5+"\t"+.5+"\t"+mafAll[i][1][1]+"\t"+df.format((mafAll[i][1][1])/(mafAll[i][1][0]+mafAll[i][1][1]+mafAll[i][1][2]))+"\n"
                    +MAFClass[i]+"\t"+.5+"\t"+1+"\t"+mafAll[i][1][2]+"\t"+df.format((mafAll[i][1][2])/(mafAll[i][1][0]+mafAll[i][1][1]+mafAll[i][1][2]))+"\n"
                    +MAFClass[i]+"\t"+1+"\t"+0+"\t"+mafAll[i][2][0]+"\t"+df.format((mafAll[i][2][0])/(mafAll[i][2][0]+mafAll[i][2][1]+mafAll[i][2][2]))+"\n"
                    +MAFClass[i]+"\t"+1+"\t"+.5+"\t"+mafAll[i][2][1]+"\t"+df.format((mafAll[i][2][1])/(mafAll[i][2][0]+mafAll[i][2][1]+mafAll[i][2][2]))+"\n"
                    +MAFClass[i]+"\t"+1+"\t"+1+"\t"+mafAll[i][2][2]+"\t"+df.format((mafAll[i][2][2])/(mafAll[i][2][0]+mafAll[i][2][1]+mafAll[i][2][2]))+"\n");
            }
            outStream.writeBytes("#Proportion unimputed:\n#MAF\tminor\thet\tmajor\n");
            for (int i= 0; i<MAFClass.length;i++) { 
                outStream.writeBytes("#"+MAFClass[i]+"\t"+mafAll[i][0][3]/mafAll[i][0][4]+"\t"+mafAll[i][1][3]/mafAll[i][1][4]+"\t"+mafAll[i][2][3]/mafAll[i][2][4]+"\n");
            }
            outStream.flush();
            outStream.close();
        } catch (Exception e) {
            if (verboseOutput) System.out.println(e);
        }
    }
    
    public void calcAccuracy(GenotypeTable imputed, GenotypeTable unimpAlign, double runtime) {
        byte diploidN= GenotypeTable.UNKNOWN_DIPLOID_ALLELE;
        boolean use= false; boolean mafOn= false; int maf= -1;
        if (mafAll!=null) {use= true; mafOn= true;}
        for (int taxon = 0; taxon < imputed.numberOfTaxa(); taxon++) {
            int keyTaxon= maskKey.taxa().indexOf(imputed.taxaName(taxon));//if key file contains fewer taxa, or different numbers, or different order
            if (keyTaxon<0) continue;//if doesn't exist, skip
            for (int site = 0; site < imputed.numberOfSites(); site++) {
                use= (mafOn && MAF[site] > -1)?true:false;
                if (use) maf= MAF[site];
                byte known = maskKey.genotype(keyTaxon, site);
                if (known == diploidN) continue;
                byte imp = imputed.genotype(taxon, site);
                if (GenotypeTableUtils.isHeterozygous(known) == true) {
                    all[1][4]++; if (use) mafAll[maf][1][4]++;
                    if (imp == diploidN) {all[1][3]++; if (use) mafAll[maf][1][3]++;}
                    else if (GenotypeTableUtils.isEqual(imp, known) == true) {all[1][1]++; if (use) mafAll[maf][1][1]++;}
                    else if (GenotypeTableUtils.isHeterozygous(imp) == false && GenotypeTableUtils.isPartiallyEqual(imp, unimpAlign.minorAllele(site)) == true) {all[1][0]++; if (use) mafAll[maf][1][0]++;}//to minor 
                    else if (GenotypeTableUtils.isHeterozygous(imp) == false && GenotypeTableUtils.isPartiallyEqual(imp, unimpAlign.majorAllele(site)) == true) {all[1][2]++; if (use) mafAll[maf][1][2]++;}
                    else {all[1][4]--;  if (use) mafAll[maf][1][4]--;}//implies >2 allele states at given genotype
                } else if (known == GenotypeTableUtils.getDiploidValue(unimpAlign.minorAllele(site),unimpAlign.minorAllele(site))) {
                    all[0][4]++; if (use) mafAll[maf][0][4]++;
                    if (imp == diploidN) {all[0][3]++;  if (use) mafAll[maf][0][3]++;}
                    else if (GenotypeTableUtils.isEqual(imp, known) == true) {all[0][0]++;  if (use) mafAll[maf][0][0]++;}
                    else if (GenotypeTableUtils.isHeterozygous(imp) == true && GenotypeTableUtils.isPartiallyEqual(imp, known) == true) {all[0][1]++;  if (use) mafAll[maf][0][1]++;}
                    else {all[0][2]++;  if (use) mafAll[maf][0][3]++;}
                } else if (known == GenotypeTableUtils.getDiploidValue(unimpAlign.majorAllele(site),unimpAlign.majorAllele(site))) {
                    all[2][4]++;  if (use) mafAll[maf][2][4]++;
                    if (imp == diploidN) {all[2][3]++;  if (use) mafAll[maf][2][3]++;}
                    else if (GenotypeTableUtils.isEqual(imp, known) == true) {all[2][2]++; if (use) mafAll[maf][2][2]++;}
                    else if (GenotypeTableUtils.isHeterozygous(imp) == true && GenotypeTableUtils.isPartiallyEqual(imp, known) == true) {all[2][1]++;  if (use) mafAll[maf][2][1]++;}
                    else {all[2][0]++; if (use) mafAll[maf][2][0]++;}
                } else continue;
            }
        }
        accuracyOut(all, runtime);
        if (MAFClass!=null) accuracyMAFOut(mafAll);
    }
    

    @Override
    public void setParameters(String[] args) {
        FILLINargs= Arrays.copyOf(args, args.length);
        if (args.length == 0) {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }

        engine.add("-hmp", "-hmpFile", true);
        engine.add("-o", "--outFile", true);
        engine.add("-d", "--donorH", true);
        engine.add("-maskKeyFile", "--maskKeyFile", true);
        engine.add("-propSitesMask", "--propSitesMask", true);
        engine.add("-mxHet", "--hetThresh", true);
        engine.add("-minMnCnt", "--minMnCnt", true);
        engine.add("-mxInbErr", "--mxInbErr", true);
        engine.add("-mxHybErr", "--mxHybErr", true);
        engine.add("-hybNNOff", "--hybNNOff", true);
        engine.add("-mxDonH", "--mxDonH", true);
        engine.add("-mnTestSite", "--mnTestSite", true);
        engine.add("-projA", "--projAlign", false);
        engine.add("-runChrMode", "--runChrMode", false);
        engine.add("-nV", "--nonVerbose",false);
        engine.parse(args);
        unimpFile = engine.getString("-hmp");
        outFile = engine.getString("-o");
        maskKeyFile = engine.getString("-maskKeyFile");
        if(engine.getBoolean("-propSitesMask")) {
            propSitesMask = Double.parseDouble(engine.getString("-propSitesMask"));
        }
        if (engine.getBoolean("-nV")) verboseOutput=false;
    }



    private void printUsage() {
        myLogger.info(
                "\n\n\nThis plugin masks files, calls FILLINImputationPlugin and calculates accuracy by unadjusted R2 at masked sites.\n"
                + "If not provided with a masked file and associated key file, random sites masked. \n"
                + "If masked data provided, please set masked sites to missing in the file for imputation and indicate masked sites bynon-missing genotypes in the associated key file.\n"
                + "Available options for the FILLINImputationAccuracyPlugin are as follows:\n"
                        + "-hmp   Input HapMap file of target genotypes to impute. Accepts all file types supported by TASSEL5\n"
                        + "-d    Donor haplotype files from output of FILLINFindHaplotypesPlugin. Use .gX in the input filename to denote the substring .gc#s# found in donor files\n"
                        + "-o     Output file; hmp.txt.gz and .hmp.h5 accepted. Required\n"
                        + "-maskKeyFile An optional key file to indicate that file is already masked for accuracy calculation. Non-missing genotypes indicate masked sites. Else, will generate own mask\n"
                        + "-propSitesMask   The proportion of non missing sites to mask for accuracy calculation if depth is not available (default:"+propSitesMask+"\n"
                        + "-mxHet   Threshold per taxon heterozygosity for treating taxon as heterozygous (no Viterbi, het thresholds)\n"
                        + "-minMnCnt    Minimum number of informative minor alleles in the search window\n"
                        + "-mxInbErr    Maximum error rate for applying one haplotype to entire site window\n"
                        + "-mxHybErr    Maximum error rate for applying Viterbi with to haplotypes to entire site window\n"
                        + "-hybNNOff    Whether to model two haplotypes as heterozygotic for focus blocks\n"
                        + "-mxDonH   Maximum number of donor hypotheses to be explored\n"
                        + "-mnTestSite   Minimum number of sites to test for IBS between haplotype and target in focus block\n"
                        + "-projA   Create a projection alignment for high density markers (default off)\n"
                        + "-nV   Supress system out if flagged\n"
        );
    }

    @Override
    public DataSet performFunction(DataSet input) {
        try {
            long time=System.currentTimeMillis();
            FILLINImputationPlugin FILLIN = new FILLINImputationPlugin();
            FILLIN.setParameters(FILLINargs);
            initiateAccuracy(unimpFile, maskKeyFile,null,propSitesMask);
            FILLIN.performFunction(null);
            double runtime= (double)(System.currentTimeMillis()-time)/(double)1000;
            calcAccuracy(ImportUtils.readGuessFormat(outFile),FILLINImputationPlugin.unimpAlign,runtime);
        }
        finally {
            fireProgress(100);
        }
        return null;
    }
    
    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "ImputeByFILLIN";
    }

    @Override
    public String getToolTipText() {
        return "Imputation that relies on a combination of HMM and Nearest Neighbor";
    }

    /**
     * Legacy approach for measuring accuracy, but need to maintain some tests.
     * @deprecated
     */
    @Deprecated
    public static int[] compareAlignment(String origFile, String maskFile, String impFile, boolean noMask) {
        boolean taxaOut=false;
        GenotypeTable oA=ImportUtils.readGuessFormat(origFile);
        System.out.printf("Orig taxa:%d sites:%d %n",oA.numberOfTaxa(),oA.numberOfSites());
        GenotypeTable mA=null;
        if(noMask==false) {mA=ImportUtils.readGuessFormat(maskFile);
            System.out.printf("Mask taxa:%d sites:%d %n",mA.numberOfTaxa(),mA.numberOfSites());
        }
        GenotypeTable iA=ImportUtils.readGuessFormat(impFile);
        System.out.printf("Imp taxa:%d sites:%d %n",iA.numberOfTaxa(),iA.numberOfSites());
        int correct=0;
        int errors=0;
        int unimp=0;
        int hets=0;
        int gaps=0;
        for (int t = 0; t < iA.numberOfTaxa(); t++) {
            int e=0,c=0,u=0,h=0;
            int oATaxa=oA.taxa().indexOf(iA.taxaName(t));
            for (int s = 0; s < iA.numberOfSites(); s++) {
                if(noMask||(oA.genotype(oATaxa, s)!=mA.genotype(t, s))) {
                    byte ib=iA.genotype(t, s);
                    byte ob=oA.genotype(oATaxa, s);
                    if((ib==UNKNOWN_DIPLOID_ALLELE)||(ob==UNKNOWN_DIPLOID_ALLELE)) {unimp++; u++;}
                    else if(ib==GAP_DIPLOID_ALLELE) {gaps++;}
                    else if(ib==ob) {
                        correct++;
                        c++;
                    } else {
                        if(isHeterozygous(ob)||isHeterozygous(ib)) {hets++; h++;}
                        else {errors++;
                            e++;
//                            if(t==0) System.out.printf("%d %d %s %s %n",t,s,oA.getBaseAsString(oATaxa, s), iA.getBaseAsString(t, s));
                        }
                    }
                }
            }
            if(taxaOut) System.out.printf("%s %d %d %d %d %n",iA.taxaName(t),u,h,c,e);
        }
        System.out.println("MFile\tIFile\tGap\tUnimp\tUnimpHets\tCorrect\tErrors");
        System.out.printf("%s\t%s\t%d\t%d\t%d\t%d\t%d%n",maskFile, impFile, gaps, unimp,hets,correct,errors);
        return new int[]{gaps, unimp,hets,correct,errors};
    }
}
