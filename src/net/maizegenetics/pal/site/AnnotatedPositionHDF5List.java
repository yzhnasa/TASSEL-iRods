package net.maizegenetics.pal.site;

import com.google.common.base.Preconditions;
import com.google.common.cache.CacheLoader;
import com.google.common.cache.LoadingCache;

import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: edbuckler
 * Date: 7/30/13
 * Time: 5:22 PM
 */
public final class AnnotatedPositionHDF5List implements AnnotatedPositionList {

    private final List<AnnotatedPosition> mySiteList;

    /*Byte representations of DNA sequences are stored in blocks of 65536 sites*/
    public static final int BLOCKSIZE=1<<16;
    public static final int blockMask=BLOCKSIZE-1;
    public static final int siteMask=~(BLOCKSIZE-1);

    private enum AlleleT {REF(3), MAJ(7), ANC(13), HIDEP(17);
        private int value;
        private AlleleT(int value) {this.value = value;}
    };
    private LoadingCache<Integer,byte[]> myScopeAlleleCache; //key (alleleType << 17)|startBlock
    private CacheLoader<Integer,byte[]> genoLoader = new CacheLoader<Integer,byte[]>() {
        public byte[] load(Integer key) {
            int block=key&blockMask;
            int startSite=block<<16;
            int alleleTint=key>>>17;
            int length=(mySiteList.size()-startSite<BLOCKSIZE)?mySiteList.size()-startSite:BLOCKSIZE;
            byte[] data=new byte[length];
            if(alleleTint==AlleleT.REF.value) {
                for (int i = startSite; i < startSite+length; i++) data[i]=mySiteList.get(i).getReferenceAllele();
            } else  if(alleleTint==AlleleT.MAJ.value) {
                for (int i = startSite; i < startSite+length; i++) data[i]=mySiteList.get(i).getGlobalMajorAllele();
            } if(alleleTint==AlleleT.ANC.value) {
                for (int i = startSite; i < startSite+length; i++) data[i]=mySiteList.get(i).getAncestralAllele();
            } if(alleleTint==AlleleT.HIDEP.value) {
                for (int i = startSite; i < startSite+length; i++) data[i]=mySiteList.get(i).getHighDepthAllele();
            }
            return data;
        } };

    private AnnotatedPositionHDF5List(ArrayList<AnnotatedPosition> builderList) {
        mySiteList=new ArrayList<AnnotatedPosition>(builderList);
    }

    private static final int getCacheKey(AlleleT aT, int site) {
        return (aT.value<<17)|(site>>16);
    }

    @Override
    public byte getReferenceAllele(int site) {
        return mySiteList.get(site).getReferenceAllele();
    }
    
    @Override
    public byte[] getReference(int startSite, int endSite) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public byte[] getReference() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public boolean hasReference() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String[] getSNPIDs() {
        String[] theIDs=new String[mySiteList.size()];
        for (int i = 0; i < theIDs.length; i++) {
            theIDs[i]=mySiteList.get(i).getSNPID();
        }
        return theIDs;
    }

    @Override
    public String getSNPID(int site) {
        return mySiteList.get(site).getSNPID();
    }

    @Override
    public int getSiteCount() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int getChromosomeSiteCount(Chromosome chromosome) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int[] getStartAndEndOfChromosome(Chromosome chromosome) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int getPositionInChromosome(int site) {
        return mySiteList.get(site).getPosition();
    }

    @Override
    public int getSiteOfPhysicalPosition(int physicalPosition, Chromosome chromosome) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int getSiteOfPhysicalPosition(int physicalPosition, Chromosome chromosome, String snpID) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int[] getPhysicalPositions() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String getChromosomeName(int site) {
        return mySiteList.get(site).getChromosome().getName();
    }

    @Override
    public Chromosome getChromosome(int site) {
        return mySiteList.get(site).getChromosome();
    }

    @Override
    public Chromosome getChromosome(String name) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public Chromosome[] getChromosomes() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int getNumChromosomes() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int[] getChromosomesOffsets() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int getIndelSize(int site) {
        return mySiteList.get(site).getKnownVariants()[1].length();
    }

    @Override
    public boolean isIndel(int site) {
        return mySiteList.get(site).isIndel();
    }

    @Override
    public String getGenomeAssembly() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public boolean isPositiveStrand(int site) {
        return (1==mySiteList.get(site).getStrand());
    }
    
    // List methods

    @Override
    public int size() {
        return mySiteList.size();
    }

    @Override
    public boolean isEmpty() {
        return mySiteList.isEmpty();
    }

    @Override
    public boolean contains(Object o) {
        return mySiteList.contains(o);
    }

    @Override
    public Iterator<AnnotatedPosition> iterator() {
        return mySiteList.iterator();
    }

    @Override
    public Object[] toArray() {
        return mySiteList.toArray();
    }

    @Override
    public <AnnotatedSite> AnnotatedSite[] toArray(AnnotatedSite[] a) {
        return mySiteList.toArray(a);
    }

    @Override
    public boolean add(AnnotatedPosition e) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    @Override
    public boolean remove(Object o) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    @Override
    public boolean containsAll(Collection<?> c) {
        return mySiteList.containsAll(c);
    }

    @Override
    public boolean addAll(Collection<? extends AnnotatedPosition> c) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    @Override
    public boolean addAll(int index, Collection<? extends AnnotatedPosition> c) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    @Override
    public boolean removeAll(Collection<?> c) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    @Override
    public boolean retainAll(Collection<?> c) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    @Override
    public void clear() {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    @Override
    public AnnotatedPosition get(int index) {
        return mySiteList.get(index);
    }

    @Override
    public AnnotatedPosition set(int index, AnnotatedPosition element) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    @Override
    public void add(int index, AnnotatedPosition element) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    @Override
    public AnnotatedPosition remove(int index) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    @Override
    public int indexOf(Object o) {
        return mySiteList.indexOf(o);
    }

    @Override
    public int lastIndexOf(Object o) {
        return mySiteList.lastIndexOf(o);
    }

    @Override
    public ListIterator<AnnotatedPosition> listIterator() {
        return mySiteList.listIterator();
    }

    @Override
    public ListIterator<AnnotatedPosition> listIterator(int index) {
        return mySiteList.listIterator(index);
    }

    @Override
    public List<AnnotatedPosition> subList(int fromIndex, int toIndex) {
        return mySiteList.subList(fromIndex, toIndex);
    }

    /**
     * A builder for creating immutable list instances, especially
     * {@code public static final} lists ("constant lists").
     *
     * <p>Example:
     * <pre>   {@code
     *   public static final ImmutableList<Color> GOOGLE_COLORS
     *       = new ImmutableList.Builder<Color>()
     *           .addAll(WEBSAFE_COLORS)
     *           .add(new Color(0, 191, 255))
     *           .build();}</pre>
     *
     * <p>Builder instances can be reused - it is safe to call {@link #build}
     * multiple times to build multiple lists in series. Each new list
     * contains the one created before it.
     */
    public static class Builder {
        private final ArrayList<AnnotatedPosition> contents = new ArrayList<AnnotatedPosition>();

        /**
         * Creates a new builder. The returned builder is equivalent to the builder
         * generated by {@link }.
         */
        public Builder() {}

        /**
         * Adds {@code element} to the {@code AnnotatedPositionList}.
         *
         * @param element the element to add
         * @return this {@code Builder} object
         * @throws NullPointerException if {@code element} is null
         */
        public Builder add(AnnotatedPosition element) {
            Preconditions.checkNotNull(element, "element cannot be null");
            contents.add(element);
            return this;
        }

        /**
         * Adds each element of {@code elements} to the {@code AnnotatedPositionList}.
         *
         * @param elements the {@code Iterable} to add to the {@code AnnotatedPositionList}
         * @return this {@code Builder} object
         * @throws NullPointerException if {@code elements} is or contains null
         */
        public Builder addAll(Iterable<? extends AnnotatedPosition> elements) {
            if (elements instanceof Collection) {
                @SuppressWarnings("unchecked")
                Collection<? extends AnnotatedPosition> collection = (Collection<? extends AnnotatedPosition>) elements;
                contents.ensureCapacity(contents.size() + collection.size());
            }
            for (AnnotatedPosition elem : elements) {
                Preconditions.checkNotNull(elem, "elements contains a null");
                contents.add(elem);
            }
            return this;
        }

        /**
         * Returns a newly-created {@code ImmutableList} based on the contents of
         * the {@code Builder}.
         */
        public AnnotatedPositionList build() {
            return new AnnotatedPositionHDF5List(contents);
        }
    }

}
