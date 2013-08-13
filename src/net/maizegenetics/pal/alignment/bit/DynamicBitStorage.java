package net.maizegenetics.pal.alignment.bit;

import com.google.common.cache.CacheBuilder;
import com.google.common.cache.CacheLoader;
import com.google.common.cache.LoadingCache;
import net.maizegenetics.pal.alignment.Alignment.ALLELE_SCOPE_TYPE;
import net.maizegenetics.pal.alignment.AlignmentUtils;
import net.maizegenetics.util.BitSet;

import java.util.concurrent.ExecutionException;

/**
 * Defines xxxx
 *
 * @author Ed Buckler
 */
public class DynamicBitStorage implements BitStorage {
    private byte[][] myGenotype; //[taxa][sites]  TODO convert this to external object.
    //currently I am testing storing the alleles as homozgyous genotypes, there may be some performance benefits
    //however this is not standard.
    private byte[] myMajorHomo; //[taxa][sites]  TODO convert this to external object.
    private byte[] myMinorHomo; //[taxa][sites]  TODO convert this to external object.
    private byte[] myHet; //[taxa][sites]  TODO convert this to external object.
    private enum SB {
        TAXA (0), SITE (1);
        public final int index;
        SB(int index) {
            this.index=index;
        }
    };

    private LoadingCache<Long,BitSet[]> bitCache; //taxon id number, BitSet[2] = {MajorBitSet, MinorBitSet}
    private CacheLoader<Long,BitSet[]> bitLoader= new CacheLoader<Long,BitSet[]>() {
        public BitSet[] load(Long key) {
            BitSet[] bs;
            if(getDirectionFromKey(key)==SB.TAXA) {
                byte[] a1=myMajorHomo;
                byte[] a2=myMinorHomo;
                int taxon=getSiteOrTaxonFromKey(key);
                long time=System.nanoTime();
                bs=AlignmentUtils.calcBitPresenceFromGenotype(myGenotype[taxon], a1, a2); //allele comp
               // bs=AlignmentUtils.calcBitPresenceFromGenotype(myGenotype[taxon], a1, a2,myHet);  //homozygyous genotype testing
                double timePerObj=(double)(System.nanoTime()-time)/(double)(myMajorHomo.length);
                double timeT=(double)(System.nanoTime()-time)/(double)(1e9);
                if(key%100==0) System.out.printf("Key:%d  Taxon:%d Time:%gs  AvgPerObj:%gns %n",key, taxon,timeT,timePerObj);
            } else {
                int site=getSiteOrTaxonFromKey(key);
                byte a1=myMajorHomo[site];
                byte a2=myMinorHomo[site];
                byte[] geno=new byte[myGenotype.length];
                for (int i=0; i<geno.length; i++) {geno[i]=myGenotype[i][site];}
                bs=AlignmentUtils.calcBitPresenceFromGenotype(geno, a1, a2);
//                  bs=null;
            }

            return bs;
        } };

    private long getKey(SB direction, ALLELE_SCOPE_TYPE aT, int siteOrTaxon) {
        return ((long)direction.index<<48) | ((long)aT.ordinal()<<32) | (long)siteOrTaxon;
    }

    private int getSiteOrTaxonFromKey(long key) {
        return (int)((key<<32)>>>32);
    }

    private SB getDirectionFromKey(long key) {
        if(key>>>48==SB.TAXA.index) return SB.TAXA;
        return SB.SITE;
    }

    @Override
    public BitSet getAllelePresenceForAllSites(int taxon, int alleleNumber) {
        try{return bitCache.get(getKey(SB.TAXA, ALLELE_SCOPE_TYPE.Global_Frequency,taxon))[alleleNumber];}
        catch (ExecutionException e) {
            e.printStackTrace();
            return null;
        }
    }

    @Override
    public BitSet getAllelePresenceForAllTaxa(int site, int alleleNumber) {
        try{return bitCache.get(getKey(SB.SITE, ALLELE_SCOPE_TYPE.Global_Frequency,site))[alleleNumber];}
        catch (ExecutionException e) {
            e.printStackTrace();
            return null;
        }
    }

    @Override
    public long[] getAllelePresenceForSitesBlock(int taxon, int alleleNumber, int startBlock, int endBlock) {
        return new long[0];  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public BitSet getPhasedAllelePresenceForAllSites(int taxon, boolean firstParent, int alleleNumber) {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public BitSet getPhasedAllelePresenceForAllTaxa(int site, boolean firstParent, int alleleNumber) {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public long[] getPhasedAllelePresenceForSitesBlock(int taxon, boolean firstParent, int alleleNumber, int startBlock, int endBlock) {
        return new long[0];  //To change body of implemented methods use File | Settings | File Templates.
    }


    public DynamicBitStorage(byte[][] myGenotype, byte[] myMajor, byte[] myMinor) {
        this.myGenotype=myGenotype;
        this.myMajorHomo=new byte[myMajor.length];
        this.myMinorHomo=new byte[myMajor.length];
        this.myHet=new byte[myMajor.length];
        for (int i=0; i<myMajor.length; i++) {
            myMajorHomo[i]=AlignmentUtils.getDiploidValue(myMajor[i],myMajor[i]);
            myMinorHomo[i]=AlignmentUtils.getDiploidValue(myMinor[i],myMinor[i]);
            myHet[i]=AlignmentUtils.getDiploidValue(myMajor[i],myMinor[i]);
        }
        bitCache= CacheBuilder.newBuilder()
                .maximumSize(1000000)
                .build(bitLoader);
    }
}
