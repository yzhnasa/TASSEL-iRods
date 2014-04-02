/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.dna.map;

/**
 * Stores variables from tag genetic mapping from GWAS. This class is used for I/O of HDF5. 
 * @author Fei Lu
 */
public class TagGWASMapInfo {
    /**Tag count in master tagCount file, unknown = Integer.MIN_VALUE*/
    public int readCount = Integer.MIN_VALUE;
    /**Physical chromosome of tag, unknown = Integer.MIN_VALUE*/
    public int pChr = Integer.MIN_VALUE;
    /**Physical position of tag, unknown = Integer.MIN_VALUE*/
    public int pPos = Integer.MIN_VALUE;
    /**If the tag is mapped by aligner*/
    public boolean ifMap = false;
    /**If the tag is reference tag*/
    public boolean ifRef = false;
    /**If the tag is unique to one position in genome*/
    public boolean ifUnique = false;
    /**Genetic mapping chromosome of a tag, unknown = Integer.MIN_VALUE*/
    public int gChr = Integer.MIN_VALUE;
    /**Genetic mapping position of a tag, unknown  = Integer.MIN_VALUE*/
    public int gPos = Integer.MIN_VALUE;
    /**P-value of genetic mapping of a tag*/
    public double gwasPValue = 1;
    /**Total number of significant site, unknown  = Integer.MIN_VALUE*/
    public int numSigSite = Integer.MIN_VALUE;
    /**Number of taxa where tag exist, unknown  = Integer.MIN_VALUE*/
    public int tagTaxaCount = Integer.MIN_VALUE;
    /**Total number of significant chromosome, unknown  = Integer.MIN_VALUE*/
    public int numSigChr = Integer.MIN_VALUE;
    /**Likelihood ratio of the most significant chromosome Vs the second most significant chromosome, log10 value, unknown  = 0*/
    public double lRatioSB = 0;
    /**Likelihood ratio of the most significant chromosome Vs the median most significant chromosome, log10 value, unknown  = 0*/
    public double lRatioMB = 0;
    /**Number of site on best chromosome having more significant p value than the most significant p value on the second best chromosome, unknown = Integer.MIN_VALUE*/
    public int numSiteOnBestChrThanSecondBest = Integer.MIN_VALUE;
    /**Starting significant site on the best chromosome, unknown  = Integer.MIN_VALUE*/
    public int sigSiteStart = Integer.MIN_VALUE;
    /**Ending significant site on the best chromosome, inclusive, unknown  = Integer.MIN_VALUE*/
    public int sigSiteEnd = Integer.MIN_VALUE;
    
    public TagGWASMapInfo () {}
    
    public TagGWASMapInfo (int readCount, int pChr, int pPos, boolean ifMap, boolean ifRef, boolean ifUnique, int gChr, int gPos, double gwasPValue, int numSigSite, int tagTaxaCount, int numSigChr,
                             double lRatioSB, double lRatioMB, int numSiteOnBestChrThanSecondBest, int sigSiteStart, int sigSiteEnd) {
        this.readCount = readCount; this.pChr = pChr; this.pPos = pPos; this.ifMap = ifMap; this.ifRef = ifRef; this.ifUnique = ifUnique; this.gChr = gChr; this.gPos = gPos;
        this.gwasPValue = gwasPValue; this.numSigSite = numSigSite; this.tagTaxaCount = tagTaxaCount; this.numSigChr = numSigChr; this.lRatioSB = lRatioSB;
        this.lRatioMB = lRatioMB; this.numSiteOnBestChrThanSecondBest = numSiteOnBestChrThanSecondBest; this.sigSiteStart = sigSiteStart; this.sigSiteEnd = sigSiteEnd;
    }
    
}
