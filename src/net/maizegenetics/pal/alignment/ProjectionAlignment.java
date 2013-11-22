/*
 * ProjectionAlignment
 */
package net.maizegenetics.pal.alignment;

import net.maizegenetics.pal.alignment.bit.BitStorage;
import net.maizegenetics.pal.position.PositionList;
import net.maizegenetics.pal.taxa.TaxaList;
import net.maizegenetics.pal.taxa.TaxaListBuilder;
import net.maizegenetics.pal.taxa.Taxon;
import net.maizegenetics.util.DonorHaplotypes;

import java.util.Map;
import java.util.SortedSet;

/**
 * This class projects high Density markers on large group of taxa through a
 * look up table system. The lookup table generally needs be built through some
 * imputation approach.
 *
 * @author ed, terry
 * @deprecated  Use the ProjectionBuilder that returns a CoreAlignment with ProjectionGenotype
 */
@Deprecated
public class ProjectionAlignment {
    private BitStorage myFreqBitStorage;
    private BitStorage myReferenceBitStorage;
    private final TaxaList myTaxaList;
    private PositionList myPositionList;
    private int[][] mySiteBreaks;  //temporary - saving not needed
    private int[][][] myHDTaxa;  //taxa ids should be saved
    private int[][] myPosBreaks;  //positions should be saved
    private final Alignment myBaseAlignment;  //high density marker alignment that is being projected.
    private Map<Taxon,SortedSet<DonorHaplotypes>> allBreakPoints;
    private int myNumSites;
    
    private int[][] cacheTaxonSiteBound, cacheTaxonDonors; 
    public long cacheUseCnt=0, lookupCnt=0;

    
    ProjectionAlignment(Alignment hdAlign, Map<Taxon,SortedSet<DonorHaplotypes>> allBreakPoints) {
        myTaxaList=new TaxaListBuilder().addAll(allBreakPoints.keySet()).build();
        myBaseAlignment = hdAlign;
        this.allBreakPoints=allBreakPoints;

//
//        mySiteBreaks = new int[getSequenceCount()][];
//        for (int taxon = 0; taxon < myHDTaxa.length; taxon++) {
//            if(myHDTaxa[taxon]==null) continue;
//            mySiteBreaks[taxon] = new int[myPosBreaks[taxon].length];
//            for (int i = 0; i < myPosBreaks[taxon].length; i++) {
//                int site = myBaseAlignment.getSiteOfPhysicalPosition(myPosBreaks[taxon][i], null);
//                if (site < 0) {
//                    site = -(site + 1);
//                }
//                this.mySiteBreaks[taxon][i] = site;
//            }
//        }
//        init();
    }
    

    
//
//
//
//    private void init() {
//        myNumSites=myBaseAlignment.getSiteCount();
//        cacheTaxonSiteBound=new int[getSequenceCount()][2];
//        cacheTaxonDonors=new int[getSequenceCount()][2];
//        for (int i = 0; i < getSequenceCount(); i++) {
//         //   translateTaxon(i,0);
//            cacheNewTaxonSiteRange(i,0);
//        }
////        System.out.println(Arrays.toString(mySiteBreaks[0]));
////        for (int i = 0; i < 16000; i++) {
////            System.out.println(i);
////            System.out.println("TOrig"+Arrays.toString(translateTaxon(0,i)));
////            System.out.println("TNew"+Arrays.toString(translateTaxonX(0,i)));
////            System.out.println("cacheB:"+Arrays.toString(cacheTaxonSiteBound[0]));
////        }
//    }


//    private int[] translateTaxon(int taxon, int site) {
//        if (mySiteBreaks[taxon] == null) {
//            return null;
//        }
//        if((cacheTaxonSiteBound[taxon][0]<=site)&&(site<=cacheTaxonSiteBound[taxon][1])) {
//            cacheUseCnt++;
//            return cacheTaxonDonors[taxon];
//        } else {
//            lookupCnt++;
//            return cacheNewTaxonSiteRange(taxon, site);
//        }
////        int b = Arrays.binarySearch(mySiteBreaks[taxon], site);
////        if (b < 0) {
////            b = -(b + 2);  //this will not work if it does not start with zero.
////        }
////        return myHDTaxa[taxon][b];
//    }
//
//    private int[] translateTaxonX(int taxon, int site) {
//        if (mySiteBreaks[taxon] == null) {
//            return null;
//        }
//        int b = Arrays.binarySearch(mySiteBreaks[taxon], site);
//        if (b < 0) {
//            b = -(b + 2);  //this will not work if it does not start with zero.
//        }
//        return myHDTaxa[taxon][b];
//    }
//
//    private int[] cacheNewTaxonSiteRange(int taxon, int site){
//        if (mySiteBreaks[taxon] == null) return null;
//        int b = Arrays.binarySearch(mySiteBreaks[taxon], site);
//        if (b < 0) {
//            b = -(b + 2);  //this will not work if it does not start with zero.
//        }
//        cacheTaxonSiteBound[taxon][0]=mySiteBreaks[taxon][b];
//        if((b+1)<mySiteBreaks[taxon].length) {
//            cacheTaxonSiteBound[taxon][1]=mySiteBreaks[taxon][b+1];
//        } else {
//            cacheTaxonSiteBound[taxon][1]=myNumSites;
//        }
//        cacheTaxonDonors[taxon]=myHDTaxa[taxon][b];
//        return myHDTaxa[taxon][b];
//    }


}
