/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.analysis.imputation;

import java.io.BufferedOutputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import net.maizegenetics.dna.map.Chromosome;
import net.maizegenetics.dna.snp.FilterGenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.GenotypeTableUtils;
import net.maizegenetics.dna.snp.ImportUtils;
import net.maizegenetics.dna.snp.NucleotideGenotypeTable;
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTable;
import net.maizegenetics.dna.snp.genotypecall.GenotypeCallTableBuilder;

/**
 *
 * @author kls283
 */
public class ImputationAccuracy {
    public static int[] MAF= null;
    private static double[][] all= new double[3][5]; //arrays held ("columns"): 0-maskedMinor, 1-maskedHet, 2-maskedMajor; each array ("rows"):0-to minor, 1-to het, 2-to major, 3-unimp, 4-total for known type
    private static double[][][] mafAll= null;//sam as all, but first array holds MAF category
    private static GenotypeTable maskKey= null;
    private static double[] MAFClass= null;//new double[]{0,.02,.05,.10,.20,.3,.4,.5,1};
    private static String outFileBase;
    
    
    public static void initiateAccuracy(GenotypeTable unimp, GenotypeTable[] donors, String maskKeyFile, double[] mafClass, double propSitesMask, String outFile) {
        outFileBase= outFile;
        if (mafClass!=null) {//if mafClass not null, assign MAF categories to sites in unimputed
            generateMAF(donors, unimp, mafClass);
            mafAll= new double[mafClass.length][3][5];
            System.out.println("Calculating accuracy within supplied MAF categories.");
        }
        //mask file if not already masked (depth or proportion) and generate key
        if (maskKeyFile!=null) {
            System.out.println("File already masked. Use input key file for calculating accuracy");
            GenotypeTable inMaskKey= ImportUtils.readGuessFormat(maskKeyFile);
            if (inMaskKey.positions()!=unimp.positions()) maskKey= filterKey(inMaskKey, unimp);
            else maskKey= inMaskKey;
        }
        //else if (unimp.depth()!=null) FILLINImputationPlugin.unimpAlign= maskFileByDepth(unimp,7, 7);
        else FILLINImputationPlugin.unimpAlign= maskPropSites(unimp,propSitesMask);
    }
    
        public static void generateMAF(GenotypeTable[] donorAlign, GenotypeTable unimpAlign, double[] mafClass) {
        ImputationAccuracy.MAF= new int[unimpAlign.numberOfSites()];
        for (GenotypeTable don:donorAlign) {
            for (int site = 0; site < don.numberOfSites(); site++) {
                int unimpSite= unimpAlign.positions().indexOf(don.positions().get(site));
                if (unimpSite < 0) {ImputationAccuracy.MAF[unimpSite]= -1; continue;}
                int search= Arrays.binarySearch(mafClass, don.minorAlleleFrequency(site));
                ImputationAccuracy.MAF[unimpSite]= search<0?Math.abs(search)-1:search;
            }
        }
    }

    //for depths 4 or more, requires hets to be called by more than one for less depth allele. returns 2-array whre index 0 is mask and 1 is key
    private static GenotypeTable maskFileByDepth(GenotypeTable a, int depthToMask, int maskDenom) {
        System.out.println("Masking file using depth\nSite depth to mask: "+depthToMask+"Divisor for physical positions to be masked: "+maskDenom);
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
            System.out.println(taxaCnt+" sites masked for "+a.taxaName(taxon)); cnt+= taxaCnt;
        }
        System.out.println(cnt+" sites masked at a depth of "+depthToMask+" (site numbers that can be divided by "+maskDenom+")");
        maskKey= GenotypeTableBuilder.getInstance(key.build(), a.positions(), a.taxa());
        return GenotypeTableBuilder.getInstance(mask.build(), a.positions(), a.taxa());
    }
    
    private static GenotypeTable maskPropSites(GenotypeTable a, double propSitesMask) {
        System.out.println("Masking file without depth\nMasking "+propSitesMask+" proportion of sites");
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
        System.out.println(cnt+" sites masked randomly not based on depth ("+expected+" expected at "+propSitesMask+")");
        maskKey= GenotypeTableBuilder.getInstance(key.build(), a.positions(), a.taxa());
        return GenotypeTableBuilder.getInstance(mask.build(), a.positions(), a.taxa());
    }
    
    //This assumes the input key contains all of the sites in unimpAlign, plus potentially more of either
    private static GenotypeTable filterKey(GenotypeTable maskKey, GenotypeTable unimpAlign) {
        System.out.println("Filtering user input key file...\nSites in original Key file: "+maskKey.numberOfSites());
        String[] unimpNames= new String[unimpAlign.numberOfSites()];
        for (int site = 0; site < unimpAlign.numberOfSites(); site++) {unimpNames[site]= unimpAlign.siteName(site);}
        ArrayList<String> keepSites= new ArrayList<>();for (int site = 0; site < maskKey.numberOfSites(); site++) {
            if (unimpAlign.positions().siteOfPhysicalPosition(maskKey.positions().physicalPositions()[site], maskKey.positions().chromosome(site))>-1) keepSites.add(maskKey.siteName(site));
        }
        FilterGenotypeTable filter= FilterGenotypeTable.getInstance(maskKey, keepSites.toArray(new String[keepSites.size()]));
        GenotypeTable newMask= GenotypeTableBuilder.getGenotypeCopyInstance(filter);
        System.out.println("Sites in new mask: "+newMask.numberOfSites());
        return newMask;
    }
    
    //this is the sample multiple r2.
    private static double pearsonR2(double[][] all) {
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
        System.out.println("Unadjusted R2 value for "+size+" comparisons: "+r2);
        return r2;
    }
    
    private static void accuracyOut(double[][] all, double time) {
        DecimalFormat df = new DecimalFormat("0.########");
        double r2= pearsonR2(all);
        try {
            File outputFile = new File(outFileBase.substring(0, outFileBase.indexOf(".hmp")) + "DepthAccuracy.txt");
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
            System.out.println("##Taxon\tTotalSitesMasked\tTotalSitesCompared\tTotalPropUnimputed\tNumMinor\tCorrectMinor\tMinorToHet\tMinorToMajor\tUnimpMinor"
                    + "\tNumHets\tHetToMinor\tCorrectHet\tHetToMajor\tUnimpHet\tNumMajor\tMajorToMinor\tMajorToHet\tCorrectMajor\tUnimpMajor\tR2");
            System.out.println("TotalByImputed\t"+(all[0][4]+all[1][4]+all[2][4])+"\t"+(all[0][4]+all[1][4]+all[2][4]-all[0][3]-all[1][3]-all[2][3])+"\t"+
                    ((all[0][3]+all[1][3]+all[2][3])/(all[0][4]+all[1][4]+all[2][4]))+"\t"+all[0][4]+"\t"+all[0][0]+"\t"+all[0][1]+"\t"+all[0][2]+"\t"+all[0][3]+
                    "\t"+all[1][4]+"\t"+all[1][0]+"\t"+all[1][1]+"\t"+all[1][2]+"\t"+all[1][3]+"\t"+all[2][4]+"\t"+all[2][0]+"\t"+all[2][1]+"\t"+all[2][2]+
                    "\t"+all[2][3]+"\t"+r2);
            System.out.println("Proportion unimputed:\nminor: "+all[0][3]/all[0][4]+"\nhet: "+all[1][3]/all[1][4]+"\nmajor: "+all[2][3]/all[2][4]);
            System.out.println("#Minor=0,Het=1,Major=2;x is masked(known), y is predicted\nx\ty\tN\tprop\n"
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
            System.out.println(e);
        }
    }
    
    private static void accuracyMAFOut(double[][][] mafAll) {
        DecimalFormat df = new DecimalFormat("0.########");
        if (MAF!=null && MAFClass!=null) try {
            File outputFile = new File(outFileBase.substring(0, outFileBase.indexOf(".hmp")) + "DepthAccuracyMAF.txt");
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
            System.out.println(e);
        }
    }
    
    public static void calcAccuracy(GenotypeTable imputed, GenotypeTable unimpAlign, double runtime) {
        byte diploidN= GenotypeTable.UNKNOWN_DIPLOID_ALLELE;
        boolean use= false; boolean mafOn= false; int maf= -1;
        if (mafAll!=null) {use= true; mafOn= true;}
        for (int taxon = 0; taxon < imputed.numberOfTaxa(); taxon++) {
            for (int site = 0; site < imputed.numberOfSites(); site++) {
                use= (mafOn && MAF[site] > -1)?true:false;
                if (use) maf= MAF[site];
                byte known = maskKey.genotype(taxon, site);
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
}
