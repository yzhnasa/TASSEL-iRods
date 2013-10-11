package net.maizegenetics.pal.alignment.genotype;

import com.google.common.cache.CacheBuilder;
import com.google.common.cache.CacheLoader;
import com.google.common.cache.LoadingCache;
import com.google.common.collect.ImmutableList;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.AlignmentUtils;
import net.maizegenetics.pal.alignment.NucleotideAlignmentConstants;
import net.maizegenetics.pal.position.Chromosome;
import net.maizegenetics.pal.util.DonorHaplotypes;

import java.util.Arrays;
import java.util.Iterator;
import java.util.NavigableSet;
import java.util.concurrent.ExecutionException;

/**
 * Defines xxxx
 *
 * @author Ed Buckler
 */
public class ProjectionGenotype extends AbstractGenotype {
    private final Alignment myBaseAlignment;  //high density marker alignment that is being projected.
    private ImmutableList<NavigableSet<DonorHaplotypes>> allBreakPoints;

    private int myNumSites;
    private int myTaxaCount;

    private static final int SHIFT_AMOUNT = 16;
    private static final int HDF5_GENOTYPE_BLOCK_SIZE = 1 << SHIFT_AMOUNT;
    public static final int SITE_BLOCK_MASK = ~(HDF5_GENOTYPE_BLOCK_SIZE - 1);
    private long lastKey=Long.MIN_VALUE;
    private byte[] lastData=null;

    private final LoadingCache<Long, byte[]> myGenoCache;
    private final CacheLoader<Long, byte[]> myGenoLoader = new CacheLoader<Long, byte[]>() {
        public byte[] load(Long key) {
    //        System.out.println("loading "+getTaxonFromKey(key));
            int startSite = getSiteStartFromKey(key) << SHIFT_AMOUNT;
            Chromosome startChr=myBaseAlignment.getChromosome(startSite);
            int startPosition=myBaseAlignment.getPositionInChromosome(startSite);
            NavigableSet<DonorHaplotypes> theBP=allBreakPoints.get(getTaxonFromKey(key));
            int length = Math.min(HDF5_GENOTYPE_BLOCK_SIZE, getSiteCount() - startSite);
            byte[] data=new byte[length];
            Arrays.fill(data,Alignment.UNKNOWN_DIPLOID_ALLELE);
            DonorHaplotypes searchDH=new DonorHaplotypes(startChr,startPosition,-1,-1,-1);
            Iterator<DonorHaplotypes> bpIter=theBP.tailSet(theBP.lower(searchDH),true).iterator();
            DonorHaplotypes currDH=bpIter.next();
            int[] currDHSiteRange=siteRangeForDonor(currDH);
            for (int i=startSite; i<startSite+length; i++) {
                if(i<currDHSiteRange[0]) continue;
                if(i>currDHSiteRange[1]) {
                   // System.out.println(i);
                    if(bpIter.hasNext()) {
                        currDH=bpIter.next();
                        currDHSiteRange=siteRangeForDonor(currDH);
                        i=currDHSiteRange[0];
                        if(i>=startSite+length) break;
                    } else {
                      break;
                    }
                }
                byte p1=myBaseAlignment.getBase(currDH.getParent1index(),i);
                byte p2=myBaseAlignment.getBase(currDH.getParent2index(),i);
                data[i-startSite]=AlignmentUtils.getUnphasedDiploidValueNoHets(p1, p2);
//                data[i-startSite]=(byte)(i%64);
            }
            return data;
        }
    };

    public ProjectionGenotype(Alignment hdAlign, ImmutableList<NavigableSet<DonorHaplotypes>> allBreakPoints) {
        super(allBreakPoints.size(), hdAlign.getSiteCount(), false, NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES);
        myBaseAlignment = hdAlign;
        this.allBreakPoints=allBreakPoints;
        myTaxaCount=allBreakPoints.size();

        myGenoCache = CacheBuilder.newBuilder()
                .maximumSize((3 * getTaxaCount()) / 2)
                .build(myGenoLoader);
    }

    private static long getCacheKey(int taxon, int site) {
        return ((long) taxon << 33) + (site / HDF5_GENOTYPE_BLOCK_SIZE);
    }

    private static int getTaxonFromKey(long key) {
        return (int) (key >>> 33);
    }

    private static int getSiteStartFromKey(long key) {
        return (int) ((key << 33) >>> 33);
    }

    private static int getStartSite(int site) {
        return site & SITE_BLOCK_MASK;
    }

    private int[] siteRangeForDonor(DonorHaplotypes dh) {
        int start=myBaseAlignment.getSiteOfPhysicalPosition(dh.getStartPosition(),dh.getChromosome());
        if(start<0) start=-(start+1);
        int end=myBaseAlignment.getSiteOfPhysicalPosition(dh.getEndPosition(),dh.getChromosome());
        if(end<0) end=-(end+1);
        return new int[]{start,end};
    }

    @Override
    public byte getBase(int taxon, int site) {
        long key = getCacheKey(taxon, site);
        if(key==lastKey) {
            return lastData[site % HDF5_GENOTYPE_BLOCK_SIZE];
        }
        try {
            byte[] data = myGenoCache.get(key);
            lastKey=key;
            lastData=data;
            return data[site % HDF5_GENOTYPE_BLOCK_SIZE];
        } catch (ExecutionException e) {
            e.printStackTrace();
            throw new IllegalStateException("HDF5ByteGenotype: getBase: Error getting base from cache.");
        }
    }

    @Override
    public String getBaseAsString(int taxon, int site) {
        return NucleotideAlignmentConstants.getNucleotideIUPAC(getBase(taxon, site));
    }

    @Override
    public String getDiploidAsString(int site, byte value) {
        return NucleotideAlignmentConstants.getNucleotideIUPAC(value);
    }


//    @Override
//    public byte[] getGenotypeForAllSites(int taxon) {
//        return myGenotype.getAllColumns(taxon);
//    }
//
//    @Override
//    public byte[] getGenotypeForSiteRange(int taxon, int start, int end) {
//        return myGenotype.getColumnRange(taxon, start, end);
//    }
//
//    @Override
//    public byte[] getGenotypeForAllTaxa(int site) {
//        return myGenotype.getAllRows(site);
//    }

    @Override
    public void transposeData(boolean siteInnerLoop) {

//        if (siteInnerLoop) {
//            if (mySiteInnerLoop == null) {
//                mySiteInnerLoop = SuperByteMatrixBuilder.getInstanceTranspose(myTaxonInnerLoop);
//            }
//            myGenotype = mySiteInnerLoop;
//        } else {
//            if (myTaxonInnerLoop == null) {
//                myTaxonInnerLoop = SuperByteMatrixBuilder.getInstanceTranspose(mySiteInnerLoop);
//            }
//            myGenotype = myTaxonInnerLoop;
//        }

    }

//    private void init() {
//        myNumSites=myBaseAlignment.getSiteCount();
//        cacheTaxonSiteBound=new int[myTaxaCount][2];
//        cacheTaxonDonors=new int[myTaxaCount][2];
//        for (int i = 0; i < myTaxaCount; i++) {
//            //   translateTaxon(i,0);
//            cacheNewTaxonSiteRange(i,0);
//        }
//    }
//
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
//
//    @Override
//    public String getBaseAsString(int taxon, int site) {
//        return NucleotideAlignmentConstants.getNucleotideIUPAC(getBase(taxon, site));
//    }
//
//    @Override
//    public String getDiploidAsString(int site, byte value) {
//        return NucleotideAlignmentConstants.getNucleotideIUPAC(value);
//    }



//    @Override
//    public byte[] getBaseRow(int taxon) {
//
//        int numBreaks = mySiteBreaks[taxon].length;
//        byte[] result = new byte[myNumSites];
//        for (int i = 0; i < numBreaks - 1; i++) {
//            int[] hdTaxon = myHDTaxa[taxon][i];
//            for (int j = mySiteBreaks[taxon][i]; j < mySiteBreaks[taxon][i + 1]; j++) {
//               // result[j] = myBaseAlignment.getBase(hdTaxon, j);
//                result[j]=AlignmentUtils.getDiploidValue(myBaseAlignment.getBase(hdTaxon[0], j), myBaseAlignment.getBase(hdTaxon[1], j));
//            }
//        }
//
//        int[] hdTaxon = myHDTaxa[taxon][numBreaks - 1];
//        for (int j = mySiteBreaks[taxon][numBreaks - 1], n = getSiteCount(); j < n; j++) {
// //           result[j] = myBaseAlignment.getBase(hdTaxon, j);
//            result[j]=AlignmentUtils.getDiploidValue(myBaseAlignment.getBase(hdTaxon[0], j), myBaseAlignment.getBase(hdTaxon[1], j));
//        }
//
//        return result;
//    }






}
