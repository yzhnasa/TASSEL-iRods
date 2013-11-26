package net.maizegenetics.dna.map;

import com.google.common.collect.ArrayListMultimap;
import net.maizegenetics.dna.map.Position.Allele;

import java.nio.IntBuffer;
import java.util.*;

/**
 * In memory immutable instance of {@link PositionList}.  Use the {@link PositionListBuilder}
 * to create the list.  This list is sorted by position.
 *
 * @author Ed Buckler
 */
final class PositionArrayList implements PositionList {
    private final List<Position> mySiteList;
    private final int numPositions;
    private final byte[][] alleles;
    private final Map<Chromosome,ChrOffPos> myChrOffPosTree;
    private final Map<String,Chromosome> myChrNameHash;

    private static class ChrOffPos {
        final int startSiteOff;
        final int endSiteOff;
        final int[] position;
        private ChrOffPos(int startSiteOff, int endSiteOff, int[] position) {
            this.startSiteOff=startSiteOff;
            this.endSiteOff=endSiteOff;
            this.position=position;
        }
    }

    PositionArrayList(ArrayList<Position> builderList) {
        this.numPositions=builderList.size();
        alleles=new byte[Allele.COUNT][numPositions];
//        refAlleles=new byte[numPositions];
//        majorAlleles=new byte[numPositions];
//        ancAlleles=new byte[numPositions];
//        hiDepAlleles=new byte[numPositions];
        ArrayListMultimap<Chromosome,Integer> pTS=ArrayListMultimap.create();
        mySiteList=new ArrayList<Position>(builderList.size());
        myChrOffPosTree=new TreeMap<Chromosome,ChrOffPos>();
        myChrNameHash=new HashMap();
        int currStart=0;
        Chromosome currChr=builderList.get(0).getChromosome();
        for (int i=0; i<builderList.size(); i++) {
            Position ap=builderList.get(i);
            for (Allele allele : Allele.values()) {
              alleles[allele.index()][i]=ap.getAllele(allele);
            }
            mySiteList.add(ap);
            if((i==(builderList.size()-1))||!ap.getChromosome().equals(currChr)) {
                int end=(i==builderList.size()-1)?i:i-1;
                myChrOffPosTree.put(currChr, new ChrOffPos(currStart, i-1, null));
                currChr=ap.getChromosome();
                currStart=i;
            }
            pTS.put(ap.getChromosome(),ap.getPosition());
        }
        for (Chromosome chr: pTS.keySet()) {
            List<Integer> p=pTS.get(chr);
            int[] intP=new int[p.size()];
            for (int i=0; i<intP.length; i++) {intP[i]=p.get(i);}
            ChrOffPos currOff=myChrOffPosTree.get(chr);
            myChrOffPosTree.put(chr, new ChrOffPos(currOff.startSiteOff, currOff.endSiteOff, intP));
            myChrNameHash.put(chr.getName(),chr);
            }
        pTS=null;
    }

    @Override
    public byte getReferenceAllele(int site) {
        return mySiteList.get(site).getAllele(Allele.REF);
    }
    
    @Override
    public byte[] getReference(int startSite, int endSite) {
        byte[] result = new byte[endSite - startSite];
        System.arraycopy(alleles[Allele.REF.index()],startSite,result,0, result.length);
        return result;
    }

    @Override
    public byte[] getReference() {
        return Arrays.copyOf(alleles[Allele.REF.index()],alleles[Allele.REF.index()].length);
    }

    @Override
    public boolean hasReference() {
        return true;
    }

    @Override
    public String siteName(int site) {
        return mySiteList.get(site).getSNPID();
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
        return mySiteList.get(site).getPosition();
    }

    @Override
    public int siteOfPhysicalPosition(int physicalPosition, Chromosome chromosome) {
        ChrOffPos cop=myChrOffPosTree.get(chromosome);
        if(cop==null) return Integer.MIN_VALUE;
        int i=Arrays.binarySearch(cop.position, physicalPosition); //AvgPerObj:227.5715ns  for 2million positions
        while((i>0)&&(physicalPosition==cop.position[i-1])) {i--;} //backup to the first position if there are duplicates
        i+=(i<0)?-cop.startSiteOff:cop.startSiteOff;

        return i;
    }

    @Override
    public int siteOfPhysicalPosition(int physicalPosition, Chromosome chromosome, String snpName) {
        int result=siteOfPhysicalPosition(physicalPosition, chromosome);
        if (result < 0) {return result;}
        else {
            if (snpName.equals(siteName(result))) {return result;
            } else {
                int index=result;
                while ((index < numPositions) && (getPositionInChromosome(index) == physicalPosition)) {
                    if (snpName.equals(siteName(index))) {return index;}
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
        return mySiteList.get(site).getChromosome().getName();
    }

    @Override
    public Chromosome getChromosome(int site) {
        return mySiteList.get(site).getChromosome();
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
        return mySiteList.get(site).getKnownVariants()[1].length();
    }

    @Override
    public boolean isIndel(int site) {
        return mySiteList.get(site).isIndel();
    }

    @Override
    public String genomeVersion() {
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
    public Iterator<Position> iterator() {
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

    /**Not supported immutable class*/
    @Override@Deprecated
    public boolean add(Position e) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    /**Not supported immutable class*/
    @Override@Deprecated
    public boolean remove(Object o) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    @Override
    public boolean containsAll(Collection<?> c) {
        return mySiteList.containsAll(c);
    }

    /**Not supported immutable class*/
    @Override@Deprecated
    public boolean addAll(Collection<? extends Position> c) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    /**Not supported immutable class*/
    @Override@Deprecated
    public boolean addAll(int index, Collection<? extends Position> c) {
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
    public Position get(int index) {
        return mySiteList.get(index);
    }

    /**Not supported immutable class*/
    @Override@Deprecated
    public Position set(int index, Position element) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    /**Not supported immutable class*/
    @Override@Deprecated
    public void add(int index, Position element) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    /**Not supported immutable class*/
    @Override@Deprecated
    public Position remove(int index) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    @Override
    public int indexOf(Object o) {
        return Collections.binarySearch(mySiteList,(Position)o);
 //       return mySiteList.indexOf(o);
    }

    @Override
    public int lastIndexOf(Object o) {
        return mySiteList.lastIndexOf(o);
    }

    @Override
    public ListIterator<Position> listIterator() {
        return listIterator(0);
    }

    @Override
    public ListIterator<Position> listIterator(final int index) {
        return new ListIterator<Position>() {
            private final ListIterator<Position> i= mySiteList.listIterator(index);
            public boolean hasNext()     {return i.hasNext();}
            public Position next()              {return i.next();}
            public boolean hasPrevious() {return i.hasPrevious();}
            public Position previous()          {return i.previous();}
            public int nextIndex()       {return i.nextIndex();}
            public int previousIndex()   {return i.previousIndex();}
            public void remove() {throw new UnsupportedOperationException();}
            public void set(Position e) {throw new UnsupportedOperationException();}
            public void add(Position e) {throw new UnsupportedOperationException();}
        };
    }

    @Override
    public List<Position> subList(int fromIndex, int toIndex) {
        return Collections.unmodifiableList(mySiteList.subList(fromIndex, toIndex));
    }

}
