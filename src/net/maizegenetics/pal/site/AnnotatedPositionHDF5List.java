package net.maizegenetics.pal.site;

import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Reader;
import com.google.common.cache.CacheBuilder;
import com.google.common.cache.CacheLoader;
import com.google.common.cache.LoadingCache;
import com.google.common.collect.Range;
import com.google.common.collect.RangeMap;
import com.google.common.collect.TreeRangeMap;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.HapMapHDF5Constants;

import java.lang.reflect.Array;
import java.nio.IntBuffer;
import java.util.*;
import java.util.concurrent.ExecutionException;

/**
 * HDF5 immutable instance of {@link AnnotatedPositionList}.  Use the {@link AnnotatedPositionArrayList.Builder}
 * to create the list.
 *
 * @author Ed Buckler
 */
public final class AnnotatedPositionHDF5List implements AnnotatedPositionList {
    private final IHDF5Reader reader;
    private final int numPositions;
//    private final byte[] refAlleles;
//    private final byte[] majorAlleles;
//    private final byte[] ancAlleles;
//    private final byte[] hiDepAlleles;
    private final Map<Chromosome,ChrOffPos> myChrOffPosTree;
    private final Map<String,Chromosome> myChrNameHash;
    private final RangeMap<Integer, Chromosome> rangeMap; //map of the start and end indices of chromosomes

    /*Byte representations of DNA sequences are stored in blocks of 65536 sites*/
    public static final int BLOCKSIZE=1<<16;
    public static final int blockMask=BLOCKSIZE-1;
    public static final int siteMask=~(BLOCKSIZE-1);

    private LoadingCache<Integer,AnnotatedPosition> mySiteList; //key site > AnnoPos
    private CacheLoader<Integer,AnnotatedPosition> annoPosLoader = new CacheLoader<Integer,AnnotatedPosition>()  {
        @Override
        public AnnotatedPosition load(Integer key) {
            ArrayList toFill=new ArrayList<Integer>();
            toFill.add(key);
            try {
                mySiteList.putAll(loadAll(toFill));
                return get(key);
            } catch (Exception e) {
                return null;
            }
        }

        @Override
        public Map<Integer, AnnotatedPosition> loadAll(Iterable<? extends Integer> keys) throws Exception {
            int key=keys.iterator().next();
            HashMap<Integer, AnnotatedPosition> result=new HashMap<Integer, AnnotatedPosition>(BLOCKSIZE);
            byte[][] afOrder;
            float[] maf;
            float[] paf;
            String[] snpIDs;
            int startSite=key&siteMask;
            int length=((numPositions-startSite)<BLOCKSIZE)?numPositions-startSite:BLOCKSIZE;
 //           System.out.println("Reading from HDF5 site anno:"+startSite);
            synchronized(reader) {
                afOrder = reader.readByteMatrixBlockWithOffset(HapMapHDF5Constants.ALLELE_FREQ_ORD, 2, length, 0l, startSite);
                maf= reader.readFloatArrayBlockWithOffset(HapMapHDF5Constants.MAF,length, startSite);
                paf= reader.readFloatArrayBlockWithOffset(HapMapHDF5Constants.SITECOV,length, startSite);
                snpIDs=reader.readStringArrayBlockWithOffset(HapMapHDF5Constants.SNP_IDS,length, startSite);
            }
            for (int i=0; i<length; i++) {
                int site=i+startSite;
                Chromosome chr=getChromosome(site);
                ChrOffPos cop=myChrOffPosTree.get(chr);
                int pos=cop.position[site-cop.startSiteOff];
                Position p=new CorePosition.Builder(chr,pos).snpName(snpIDs[i]).build();
                AnnotatedPosition ap=new CoreAnnotatedPosition.Builder(p)
                        .majAllele(afOrder[0][i])
                        .maf(maf[i]).siteCoverage(paf[i]).build();
                result.put(site,ap);
            }
            return result;
        }
    };

    private class ChrOffPos {
        final int startSiteOff;
        final int endSiteOff;
        final int[] position;
        private ChrOffPos(int startSiteOff, int endSiteOff, int[] position) {
            this.startSiteOff=startSiteOff;
            this.endSiteOff=endSiteOff;
            this.position=position;
        }
    }

    private AnnotatedPositionHDF5List(IHDF5Reader reader) {
        this.reader=reader;
        int[] variableSites = reader.readIntArray(HapMapHDF5Constants.POSITIONS);
        this.numPositions=variableSites.length;
        String[] lociStrings = reader.readStringArray(HapMapHDF5Constants.LOCI);
        ArrayList<Chromosome> chrs=new ArrayList<Chromosome>();
        for (String ls : lociStrings) {
            chrs.add(new Chromosome(ls));
        }
        int[] locusIndices = reader.readIntArray(HapMapHDF5Constants.LOCUS_INDICES);
        myChrOffPosTree=new TreeMap<Chromosome,ChrOffPos>();
        myChrNameHash=new HashMap<String,Chromosome>();
        rangeMap = TreeRangeMap.create();
        int currStart=0;
        int currLocusIndex=locusIndices[0];
        for (int i=0; i<locusIndices.length; i++) {
            if((i==(locusIndices.length-1))||currLocusIndex!=locusIndices[i]) {
                int end=(i==locusIndices.length-1)?i:i-1;
                int[] cPos=Arrays.copyOfRange(variableSites,currStart,end+1);
                Chromosome currChr=chrs.get(currLocusIndex);
                myChrOffPosTree.put(currChr, new ChrOffPos(currStart, end, cPos));
                myChrNameHash.put(currChr.getName(),currChr);
                rangeMap.put(Range.closed(currStart,end),currChr);
                currStart=i;
                currLocusIndex=locusIndices[i];
            }
        }
        mySiteList= CacheBuilder.newBuilder()
                .maximumSize(1000000)
                .build(annoPosLoader);
    }

    @Override
    public byte getReferenceAllele(int site) {
        try {
            return mySiteList.get(site).getReferenceAllele();
        } catch (ExecutionException e) {
            e.printStackTrace();
            return Alignment.UNKNOWN_ALLELE;
        }
    }

    @Override
    public byte[] getReference(int startSite, int endSite) {
        byte[] result = new byte[endSite - startSite];
        //System.arraycopy(refAlleles,startSite,result,0, result.length);
        return result;
    }

    @Override
    public byte[] getReference() {
        return null;
       // return Arrays.copyOf(refAlleles,refAlleles.length);
    }

    @Override
    public boolean hasReference() {
        return true;
    }

    @Override
    public String[] getSNPIDs() {
        try {
            String[] theIDs=new String[numPositions];
            for (int i = 0; i < theIDs.length; i++) {
                theIDs[i]=mySiteList.get(i).getSNPID();
            }
            return theIDs;
        } catch (ExecutionException e) {
            e.printStackTrace();
            return null;
        }
    }

    @Override
    public String getSNPID(int site) {
        try {
            return mySiteList.get(site).getSNPID();
        } catch (ExecutionException e) {
            e.printStackTrace();
            return null;
        }
    }

    @Override
    public int getSiteCount() {
        return numPositions;
    }

    @Override
    public int getChromosomeSiteCount(Chromosome chromosome) {
        return myChrOffPosTree.get(chromosome).position.length;
    }

    @Override
    public int[] getStartAndEndOfChromosome(Chromosome chromosome) {
        ChrOffPos cop=myChrOffPosTree.get(chromosome);
        if(cop==null) return null;
        return new int[]{cop.startSiteOff,cop.endSiteOff};
    }

    @Override
    public int getPositionInChromosome(int site) {
        Chromosome chr=rangeMap.get(site);
        ChrOffPos cop=myChrOffPosTree.get(chr);
        return cop.position[site-cop.startSiteOff];
    }

    @Override
    public int getSiteOfPhysicalPosition(int physicalPosition, Chromosome chromosome) {
        ChrOffPos cop=myChrOffPosTree.get(chromosome);
        if(cop==null) return Integer.MIN_VALUE;
        int i=Arrays.binarySearch(cop.position, physicalPosition); //AvgPerObj:227.5715ns  for 2million positions
        i+=(i<0)?-cop.startSiteOff:cop.startSiteOff;
        while((i>0)&&(physicalPosition==get(i-1).getPosition())) {i--;} //backup to the first position if there are duplicates
        return i;
    }

    @Override
    public int getSiteOfPhysicalPosition(int physicalPosition, Chromosome chromosome, String snpID) {
        int result=getSiteOfPhysicalPosition(physicalPosition, chromosome);
        if (result < 0) {return result;}
        else {
            if (snpID.equals(getSNPID(result))) {return result;
            } else {
                int index=result;
                while ((index < numPositions) && (getPositionInChromosome(index) == physicalPosition)) {
                    if (snpID.equals(getSNPID(index))) {return index;}
                    result++;
                }
                return -result - 1;
            }
        }
    }

    @Override
    public int[] getPhysicalPositions() {
        int[] result=new int[numPositions];
        IntBuffer ib=IntBuffer.wrap(result);
        for (ChrOffPos cop: myChrOffPosTree.values()) {
            ib.put(cop.position);
        }
        return result;
    }

    @Override
    public String getChromosomeName(int site) {
        return rangeMap.get(site).getName();
    }

    @Override
    public Chromosome getChromosome(int site) {
        return rangeMap.get(site);
    }

    @Override
    public Chromosome getChromosome(String name) {
        return myChrNameHash.get(name);
    }

    @Override
    public Chromosome[] getChromosomes() {
        return myChrOffPosTree.keySet().toArray(new Chromosome[0]);
    }

    @Override
    public int getNumChromosomes() {
        return myChrOffPosTree.size();
    }

    @Override
    public int[] getChromosomesOffsets() {
        int[] result=new int[myChrOffPosTree.size()];
        int index=0;
        for (ChrOffPos cop: myChrOffPosTree.values()) {
            result[index++]=cop.startSiteOff;
        }
        return result;
    }

    @Override
    public int getIndelSize(int site) {
        try{return mySiteList.get(site).getKnownVariants()[1].length();}
    catch (ExecutionException e) {
        e.printStackTrace();
        return -1;
    }
    }

    @Override
    public boolean isIndel(int site) {
        try{return mySiteList.get(site).isIndel();}
        catch (ExecutionException e) {
            e.printStackTrace();
            return false;
        }
    }

    @Override
    public String getGenomeAssembly() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public boolean isPositiveStrand(int site) {
        return true;
    }

    // List methods

    @Override
    public int size() {
        return numPositions;
    }

    @Override
    public boolean isEmpty() {
        return (numPositions==0);
    }

    @Override
    public boolean contains(Object o) {
        if(o instanceof Position) {
            Position p=(Position)o;
            int site=getSiteOfPhysicalPosition(p.getPosition(),p.getChromosome());
            //test for SNP ID also?
            if(site>=0) return true;
        }
        return false;
    }

    @Override
    public Iterator<AnnotatedPosition> iterator() {
        //TODO don't know the best way to do this
        return null;
       // return mySiteList.iterator();
    }

    @Override
    public Object[] toArray() {
        AnnotatedPosition[] aps=new AnnotatedPosition[numPositions];
        for (int i=0; i<numPositions; i++) {
            aps[i]=get(i);
        }
        return aps;
    }

    @Override
    public <AnnotatedPosition> AnnotatedPosition[] toArray(AnnotatedPosition[] a) {
        if (a.length < numPositions) {
            // If array is too small, allocate the new one with the same component type
            a = (AnnotatedPosition[])Array.newInstance(a.getClass().getComponentType(), numPositions);
        } else if (a.length > numPositions) {
            // If array is to large, set the first unassigned element to null
            a[numPositions] = null;
        }
        for (int i=0; i<numPositions; i++) {
            a[i]=(AnnotatedPosition)get(i);
        }
        return a;
    }

    /**Not supported immutable class*/
    @Override@Deprecated
    public boolean add(AnnotatedPosition e) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    /**Not supported immutable class*/
    @Override@Deprecated
    public boolean remove(Object o) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    @Override
    public boolean containsAll(Collection<?> c) {
        throw new UnsupportedOperationException("Not implemented yet.");
    }

    /**Not supported immutable class*/
    @Override@Deprecated
    public boolean addAll(Collection<? extends AnnotatedPosition> c) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    /**Not supported immutable class*/
    @Override@Deprecated
    public boolean addAll(int index, Collection<? extends AnnotatedPosition> c) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    /**Not supported immutable class*/
    @Override@Deprecated
    public boolean removeAll(Collection<?> c) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    /**Not supported immutable class*/
    @Override@Deprecated
    public boolean retainAll(Collection<?> c) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    /**Not supported immutable class*/
    @Override@Deprecated
    public void clear() {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    @Override
    public AnnotatedPosition get(int index) {
        try {
            return mySiteList.get(index);
        } catch (ExecutionException e) {
            return null;
        }
    }

    /**Not supported immutable class*/
    @Override@Deprecated
    public AnnotatedPosition set(int index, AnnotatedPosition element) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    /**Not supported immutable class*/
    @Override@Deprecated
    public void add(int index, AnnotatedPosition element) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    /**Not supported immutable class*/
    @Override@Deprecated
    public AnnotatedPosition remove(int index) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    @Override
    public int indexOf(Object o) {
        throw new UnsupportedOperationException("Not implemented yet.");
        //return mySiteList.indexOf(o);
    }

    @Override
    public int lastIndexOf(Object o) {
        throw new UnsupportedOperationException("Not implemented yet.");
       // return mySiteList.lastIndexOf(o);
    }

    @Override
    public ListIterator<AnnotatedPosition> listIterator() {
        throw new UnsupportedOperationException("Not implemented yet.");
     //   return mySiteList.listIterator();
    }

    @Override
    public ListIterator<AnnotatedPosition> listIterator(int index) {
        throw new UnsupportedOperationException("Not implemented yet.");
       // return mySiteList.listIterator(index);
    }

    @Override
    public List<AnnotatedPosition> subList(int fromIndex, int toIndex) {
        throw new UnsupportedOperationException("Not implemented yet.");
        //return mySiteList.subList(fromIndex, toIndex);
    }

    /**
     * A builder for creating immutable AnnotatedPositionArrayList
     *
     * <p>Example:
     * <pre>   {@code
     *   AnnotatedPositionArrayList.Builder b=new AnnotatedPositionArrayList.Builder();
     *   for (int i = 0; i <size; i++) {
     *       AnnotatedPosition ap=new CoreAnnotatedPosition.Builder(chr[chrIndex[i]],pos[i]).refAllele(refSeq[i]).build();
     *       b.add(ap);
     *       }
     *   instance=b.build();}
     *
     * <p>Builder instances can be reused - it is safe to call {@link #build}
     * multiple times to build multiple lists in series. Each new list
     * contains the one created before it.
     */
    public static class Builder {
        private final String hdf5FileName;
        private final IHDF5Reader reader;
        /**
         * Creates a new builder. The returned builder is equivalent to the builder
         * generated by {@link }.
         */
        public Builder(String hdf5FileName) {
            this.hdf5FileName=hdf5FileName;
            this.reader = HDF5Factory.openForReading(hdf5FileName);
        }

        /**
         * Returns a newly-created {@code ImmutableList} based on the contents of
         * the {@code Builder}.
         */
        public AnnotatedPositionList build() {
            return new AnnotatedPositionHDF5List(reader);
        }
    }
}


