package net.maizegenetics.pal.alignment.bit;

import com.google.common.cache.CacheBuilder;
import com.google.common.cache.CacheLoader;
import com.google.common.cache.LoadingCache;
import net.maizegenetics.pal.alignment.Alignment.ALLELE_SCOPE_TYPE;
import net.maizegenetics.pal.alignment.AlignmentUtils;
import net.maizegenetics.pal.alignment.genotype.Genotype;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.ProgressListener;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.ExecutionException;

/**
 * Defines xxxx
 *
 * @author Ed Buckler
 */
public class DynamicBitStorage implements BitStorage {
    private Genotype myGenotype; //[taxa][sites]  TODO convert this to external object.
    //currently I am testing storing the alleles as homozgyous genotypes, there may be some performance benefits
    //however this is not standard.
    private ALLELE_SCOPE_TYPE myPreferredScope=ALLELE_SCOPE_TYPE.Frequency;
    private byte[] myPrefMajor; //[taxa][sites]  TODO convert this to external object.
    private byte[] myPrefMinor; //[taxa][sites]  TODO convert this to external object.

    private byte[] myHet; //[taxa][sites]  TODO convert this to external object.
    private final int myTaxaCount;
    private final int mySiteCount;
    public static final int SBoff=58;
    public static final int AllScopeoff=52;


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
                byte[] a1=myPrefMajor;
                byte[] a2=myPrefMinor;
                int taxon=getSiteOrTaxonFromKey(key);
                bs=AlignmentUtils.calcBitPresenceFromGenotype(myGenotype.getBaseRow(taxon), a1, a2); //allele comp
                return bs;
            } else {
                ArrayList toFill=new ArrayList<Integer>();
                toFill.add(key);
                try {
                    bitCache.putAll(loadAll(toFill));
                    return bitCache.get(key);
                } catch (Exception e) {
                    e.printStackTrace();
                    return null;
                }
//                  bs=null;
            }
        }

        @Override
        public Map<Long, BitSet[]> loadAll(Iterable<? extends Long> keys) throws Exception {
            long key=keys.iterator().next();
            //This pivoting code is needed if myGenotype is store in taxa direction
            //It runs about 7 times faster than getting base sequentially across taxa.
            HashMap<Long, BitSet[]> result=new HashMap<Long, BitSet[]>(64);
            int site=getSiteOrTaxonFromKey(key);
            int length=(mySiteCount-site<64)?mySiteCount-site:64;
            byte[][] myGenotypeTBlock=new byte[length][myTaxaCount];
            for (int t=0; t<myTaxaCount; t++) {
                for (int s=0; s<myGenotypeTBlock.length; s++) {
                    myGenotypeTBlock[s][t]=myGenotype.getBase(t,site+s);
                }
            }
            for (int i=0; i<length; i++) {
                byte a1=myPrefMajor[site+i];
                byte a2=myPrefMinor[site+i];
                BitSet[] bs=AlignmentUtils.calcBitPresenceFromGenotype(myGenotypeTBlock[i], a1, a2);
                result.put(getKey(SB.SITE,myPreferredScope,site+i),bs);
            }
            return result;
        }
    };

    private long getKey(SB direction, ALLELE_SCOPE_TYPE aT, int siteOrTaxon) {
        return ((long)direction.index<<SBoff) | ((long)aT.ordinal()<<AllScopeoff) | (long)siteOrTaxon;
    }

    private int getSiteOrTaxonFromKey(long key) {
        return (int)((key<<32)>>>32);
    }

    private SB getDirectionFromKey(long key) {
        if(key>>>SBoff==SB.TAXA.index) return SB.TAXA;
        return SB.SITE;
    }

    @Override
    public BitSet getAllelePresenceForAllSites(int taxon, int alleleNumber) {
        try{return bitCache.get(getKey(SB.TAXA, myPreferredScope,taxon))[alleleNumber];}
        catch (ExecutionException e) {
            e.printStackTrace();
            return null;
        }
    }

    @Override
    public BitSet getAllelePresenceForAllTaxa(int site, int alleleNumber) {
        try{return bitCache.get(getKey(SB.SITE, myPreferredScope,site))[alleleNumber];}
        catch (ExecutionException e) {
            e.printStackTrace();
            return null;
        }
    }

    @Override
    public long[] getAllelePresenceForSitesBlock(int taxon, int alleleNumber, int startBlock, int endBlock) {
        BitSet result=getAllelePresenceForAllSites(taxon, alleleNumber);
        if(result==null) return new long[0];
        return result.getBits(startBlock, endBlock-1);   //BitSet is inclusive, while this method is exclusive.
    }

    @Override
    public BitSet getPhasedAllelePresenceForAllSites(int taxon, boolean firstParent, int alleleNumber) {
        throw new UnsupportedOperationException("Not implemented yet");
    }

    @Override
    public BitSet getPhasedAllelePresenceForAllTaxa(int site, boolean firstParent, int alleleNumber) {
        throw new UnsupportedOperationException("Not implemented yet");
    }

    @Override
    public long[] getPhasedAllelePresenceForSitesBlock(int taxon, boolean firstParent, int alleleNumber, int startBlock, int endBlock) {
        throw new UnsupportedOperationException("Not implemented yet");
    }
    
    @Override
    public boolean isSBitFriendly() {
        return true;
    }

    @Override
    public boolean isTBitFriendly() {
        return true;
    }

    @Override
    public void optimizeForTaxa(ProgressListener listener) {
        for (int t=0; t<myTaxaCount; t++) {
            try{bitLoader.load(getKey(SB.TAXA,myPreferredScope,t)); }
            catch (Exception e) {
                e.printStackTrace();
            }
            if (listener != null) {
                listener.progress((int) (((double) (t + 1) / (double)myTaxaCount) * 100.0), null);
            }
        }
        System.out.println("optimizeForSites:"+bitCache.size());
    }

    @Override
    public void optimizeForSites(ProgressListener listener) {
        for (int s=0; s<mySiteCount; s+=64) {
            try{bitCache.get(getKey(SB.SITE, myPreferredScope, s)); }
            catch (Exception e) {
                e.printStackTrace();
            }
            if ((listener!=null)&&(s%1000==0)) {
                listener.progress((int) (((double) (s + 1) / (double)mySiteCount) * 100.0), null);
            }
        }
        System.out.println("optimizeForSites:"+bitCache.size());
    }


    public DynamicBitStorage(Genotype myGenotype, byte[] myMajor, byte[] myMinor) {
        this.myGenotype=myGenotype;
        mySiteCount=myGenotype.getSiteCount();
        myTaxaCount=myGenotype.getTaxaCount();
        this.myPrefMajor=Arrays.copyOf(myMajor,myMajor.length);
        this.myPrefMinor=Arrays.copyOf(myMinor,myMinor.length);
        bitCache= CacheBuilder.newBuilder()
                .maximumSize(5_000_000)
                .build(bitLoader);
    }
}
