package net.maizegenetics.pal.alignment.genotype;

import com.google.common.cache.CacheBuilder;
import com.google.common.cache.CacheLoader;
import com.google.common.cache.LoadingCache;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Range;
import com.google.common.collect.RangeMap;
import com.google.common.collect.TreeRangeMap;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.AlignmentUtils;
import net.maizegenetics.pal.alignment.NucleotideAlignmentConstants;
import net.maizegenetics.pal.position.Chromosome;
import net.maizegenetics.pal.util.DonorHaplotypes;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.NavigableSet;
import java.util.concurrent.ExecutionException;

/**
 * Projection genotype use defined haplotypes and breakpoints that point to high density
 * genotypes (baseAlignment).  These are used to efficiently store and connect low density maps with imputed high density genotypes.
 * <p></p>
 * The alignment built by this builder is a CoreAlignment with a ProjectionGenotype.  The taxa indice come from the
 * projection alignment file, while the site indices are the same as the base alignment.
 * TODO this implement a projection interface with the getDonorHaplotypes method
 *
 * @author Ed Buckler
 */
public class ProjectionGenotype extends AbstractGenotype {
    private final Alignment myBaseAlignment;  //high density marker alignment that is being projected.
    private ImmutableList<NavigableSet<DonorHaplotypes>> allBreakPoints;

    private static final int SHIFT_AMOUNT = 16;
    private static final int HDF5_GENOTYPE_BLOCK_SIZE = 1 << SHIFT_AMOUNT;
    public static final int SITE_BLOCK_MASK = ~(HDF5_GENOTYPE_BLOCK_SIZE - 1);
    private long lastKey=Long.MIN_VALUE;
    private byte[] lastData=null;

    private final LoadingCache<Long, byte[]> myGenoCache;
    private final CacheLoader<Long, byte[]> myGenoLoader = new CacheLoader<Long, byte[]>() {
        public byte[] load(Long key) {
            if(getTaxonFromKey(key)<4) System.out.println("loading "+getTaxonFromKey(key));
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
            if(currDHSiteRange[1]>startSite+length) currDHSiteRange[1]=startSite+length;
//            byte[] p1a=myBaseAlignment.getBaseRange(currDH.getParent1index(),currDHSiteRange[0], currDHSiteRange[1]);
//            byte[] p2a=myBaseAlignment.getBaseRange(currDH.getParent2index(),currDHSiteRange[0], currDHSiteRange[1]);
            for (int i=startSite; i<startSite+length; i++) {
                if(i<currDHSiteRange[0]) continue;
                if(i>currDHSiteRange[1]) {
                   // System.out.println(i);
                    if(bpIter.hasNext()) {
                        currDH=bpIter.next();
                        currDHSiteRange=siteRangeForDonor(currDH);
                        i=currDHSiteRange[0];
                        if(i>=startSite+length) break;
                        if(currDHSiteRange[1]>startSite+length) currDHSiteRange[1]=startSite+length;
//                        p1a=myBaseAlignment.getBaseRange(currDH.getParent1index(),currDHSiteRange[0], currDHSiteRange[1]);
//                        p2a=myBaseAlignment.getBaseRange(currDH.getParent2index(),currDHSiteRange[0], currDHSiteRange[1]);
                    } else {
                      break;
                    }
                }
                byte p1=myBaseAlignment.getBase(currDH.getParent1index(),i);
                byte p2=myBaseAlignment.getBase(currDH.getParent2index(),i);
//                byte p1=p1a[i-currDHSiteRange[0]];
//                byte p2=p2a[i-currDHSiteRange[0]];
                data[i-startSite]=AlignmentUtils.getUnphasedDiploidValueNoHets(p1, p2);
//                data[i-startSite]=(byte)(i%64);
            }
            return data;
        }
    };

    private ArrayList<RangeMap<Integer,DonorSiteHaps>> breakMaps;
    private DonorSiteHaps[] currentDSH;

    private final LoadingCache<Integer, byte[]> mySiteCache;
    private final CacheLoader<Integer, byte[]> mySiteLoader = new CacheLoader<Integer, byte[]>() {
        @Override
        public byte[] load(Integer site) throws Exception {
            byte[] baseGeno=new byte[myBaseAlignment.getTaxaCount()];
            for (int t=0; t<baseGeno.length; t++) baseGeno[t]=myBaseAlignment.getBase(t,site);
            byte[] data=new byte[getTaxaCount()];
            for (int t=0; t<data.length; t++) {
                if((currentDSH[t]==null)||(!currentDSH[t].containsSite(site))) {
                    currentDSH[t]=breakMaps.get(t).get(site);
                    //TODO consider null
                }
                byte p1=baseGeno[currentDSH[t].getParent1index()];
                byte p2=baseGeno[currentDSH[t].getParent2index()];
                data[t]=AlignmentUtils.getUnphasedDiploidValueNoHets(p1, p2);
            }
            return data;
        }
    };

    public ProjectionGenotype(Alignment hdAlign, ImmutableList<NavigableSet<DonorHaplotypes>> allBreakPoints) {
        super(allBreakPoints.size(), hdAlign.getSiteCount(), false, NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES);
        myBaseAlignment = hdAlign;
        this.allBreakPoints=allBreakPoints;
        myGenoCache = CacheBuilder.newBuilder()
                .maximumSize((10 * getTaxaCount()) / 2)
                .build(myGenoLoader);
        breakMaps=new ArrayList<>(getTaxaCount());
        for (NavigableSet<DonorHaplotypes> allBreakPoint : allBreakPoints) {
            RangeMap<Integer,DonorSiteHaps> tRM=TreeRangeMap.create();
            for (DonorHaplotypes dh : allBreakPoint) {
                int[] siteRange=siteRangeForDonor(dh);
                DonorSiteHaps dsh=new DonorSiteHaps(siteRange[0],siteRange[1],dh.getParent1index(),dh.getParent2index());
                tRM.put(Range.closed(siteRange[0],siteRange[1]),dsh);
                //TODO consider putting in blank range maps
            }
            breakMaps.add(tRM);
        }
        currentDSH=new DonorSiteHaps[getTaxaCount()];
        mySiteCache = CacheBuilder.newBuilder()
                .maximumSize(100_000)
                .build(mySiteLoader);
    }

    public NavigableSet<DonorHaplotypes> getDonorHaplotypes(int taxon) {
        return allBreakPoints.get(taxon);
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

//    @Override
//    public byte getBase(int taxon, int site) {
//        try {
//            byte[] data = mySiteCache.get(site);
//            return data[taxon];
//        } catch (ExecutionException e) {
//            e.printStackTrace();
//            throw new IllegalStateException("HDF5ByteGenotype: getBase: Error getting base from cache.");
//        }
//    }

    public byte[] getAllBaseForSite(int site) {
        try {
            byte[] data = mySiteCache.get(site);
            return data;
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




   private class DonorSiteHaps {
        private final int startSite;
        private final int endSite;
        private final int parent1index;
        private final int parent2index;

       private DonorSiteHaps(int startSite, int endSite, int parent1index, int parent2index) {
           this.startSite=startSite;
           this.endSite=endSite;
           this.parent1index=parent1index;
           this.parent2index=parent2index;
       }

       private int getStartSite() {
           return startSite;
       }

       private int getEndSite() {
           return endSite;
       }

       private int getParent1index() {
           return parent1index;
       }

       private int getParent2index() {
           return parent2index;
       }

       private boolean containsSite(int site) {
           if((site<startSite)||(site>endSite)) return false;
           return true;
       }
   }

}
