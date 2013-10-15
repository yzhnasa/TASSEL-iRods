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
 * The alignment built by this builder is a CoreAlignment with a ProjectionGenotype.  The taxa indices come from the
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
    private byte[] donorForCachedSite;
    private byte[] projForCachedTaxon;
    private int cachedSite=-1;
    int[] primDSH; //startSite,endSite,parent1,parent2 array for the


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
        primDSH=new int[myTaxaCount*4];
        Arrays.fill(primDSH,Integer.MIN_VALUE);
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

//        @Override
//    public byte getBase(int taxon, int site) {
//        if((currentDSH[taxon]==null)||(!currentDSH[taxon].containsSite(site))) {
//            currentDSH[taxon]=breakMaps.get(taxon).get(site);
//            //TODO consider null
//        }
//        byte p1=myBaseAlignment.getBase(currentDSH[taxon].getParent1index(),site);
//        byte p2=myBaseAlignment.getBase(currentDSH[taxon].getParent2index(),site);
//        return AlignmentUtils.getUnphasedDiploidValueNoHets(p1, p2);
//    }



    @Override
    public byte getBase(int taxon, int site) {
        //test transpose problems
        if(site!=cachedSite) {
            donorForCachedSite=myBaseAlignment.getGenotypeMatrix().getGenotypeForAllTaxa(site);
            cachedSite=site;
        }
        int primPos=taxon<<2;
        if((site<primDSH[primPos++])||(site>primDSH[primPos++])) {
            DonorSiteHaps currentDSH=breakMaps.get(taxon).get(site);
            primPos=taxon<<2;
            primDSH[primPos++]=currentDSH.getStartSite();
            primDSH[primPos++]=currentDSH.getEndSite();
            primDSH[primPos++]=currentDSH.getParent1index();
            primDSH[primPos]=currentDSH.getParent2index();
            primPos=(taxon<<2)+2;
            //TODO consider null
        }
 //       if(primDSH[primPos]==primDSH[primPos+1]) return donorForCachedSite[primDSH[primPos]];
        return AlignmentUtils.getUnphasedDiploidValueNoHets(donorForCachedSite[primDSH[primPos]], donorForCachedSite[primDSH[primPos+1]]);
    }



    @Override
    public String getBaseAsString(int taxon, int site) {
        return NucleotideAlignmentConstants.getNucleotideIUPAC(getBase(taxon, site));
    }

    @Override
    public String getDiploidAsString(int site, byte value) {
        return NucleotideAlignmentConstants.getNucleotideIUPAC(value);
    }


    @Override
    public void transposeData(boolean siteInnerLoop) {
        myBaseAlignment.getGenotypeMatrix().transposeData(siteInnerLoop);
    }


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
