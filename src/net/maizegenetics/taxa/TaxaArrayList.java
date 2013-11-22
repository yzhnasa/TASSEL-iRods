package net.maizegenetics.taxa;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;
import org.apache.log4j.Logger;

import java.util.*;

/**
 * In memory immutable instance of {@link TaxaList}. Basic list of taxa
 * (samples) that are used in Alignments and other purposes.
 *
 * Use {@link TaxaListBuilder} to instantiate.
 *
 * @author Ed Buckler
 *
 */
class TaxaArrayList implements TaxaList {

    private static final Logger myLogger = Logger.getLogger(TaxaArrayList.class);
    private final List<Taxon> myTaxaList;
    private final int myNumTaxa;
    private final Multimap<String, Integer> myNameToIndex;

    TaxaArrayList(TaxaListBuilder builder) {

        List<Taxon> srcList = builder.getImmutableList();
        myTaxaList = new ArrayList<Taxon>(srcList.size());
        myNumTaxa = srcList.size();
        myNameToIndex = HashMultimap.create(srcList.size() * 2, 1);
        int index = 0;
        for (Taxon Taxon : srcList) {
            myTaxaList.add(Taxon);
            if (myNameToIndex.containsKey(Taxon.getName())) {
                myLogger.warn("init: Taxa name is duplicated :" + Taxon.getName());
            }
            myNameToIndex.put(Taxon.getName(), index);

            //TODO Ed, we need to talk about this. -Terry
            //if (!Taxon.getFullName().equals(Taxon.getName())) {
            //    myNameToIndex.put(Taxon.getName(), index);
            //}

            index++;
        }
    }

    @Override
    public int getTaxaCount() {
        return myNumTaxa;
    }

    @Override
    public String getTaxaName(int index) {
        return myTaxaList.get(index).getName();
    }

    @Override
    public int size() {
        return myNumTaxa;
    }

    @Override
    public List<Integer> getIndicesMatchingTaxon(String name) {
        return new ArrayList<Integer>(myNameToIndex.get(name));
    }

    @Override
    public List<Integer> getIndicesMatchingTaxon(Taxon taxon) {
        return new ArrayList<Integer>(myNameToIndex.get(taxon.getName()));
    }

    @Override
    public boolean isEmpty() {
        return myTaxaList.isEmpty();
    }

    @Override
    public boolean contains(Object o) {
        return myTaxaList.contains(o);
    }

    @Override
    public Iterator<Taxon> iterator() {
        return myTaxaList.iterator();
    }

    @Override
    public Object[] toArray() {
        return myTaxaList.toArray();
    }

    @Override
    public <T> T[] toArray(T[] a) {
        return myTaxaList.toArray(a);
    }

    @Override
    public boolean add(Taxon Taxon) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    @Override
    public boolean remove(Object o) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    @Override
    public boolean containsAll(Collection<?> c) {
        return myTaxaList.containsAll(c);
    }

    @Override
    public boolean addAll(Collection<? extends Taxon> c) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    @Override
    public boolean addAll(int index, Collection<? extends Taxon> c) {
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
    public Taxon get(int index) {
        return myTaxaList.get(index);
    }

    @Override
    public Taxon set(int index, Taxon element) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    @Override
    public void add(int index, Taxon element) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    @Override
    public Taxon remove(int index) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    @Override
    public int indexOf(Object o) {
        Taxon at=(Taxon)o;   //uses the hashMap to do this quickly
        Collection<Integer> result=myNameToIndex.get(at.getName());
        for (Integer i : result) {
            if(myTaxaList.get(i).equals(at)) return i;
        }
        return -1;
    }

    @Override
    public int lastIndexOf(Object o) {
        return myTaxaList.lastIndexOf(o);
    }

    @Override
    public ListIterator<Taxon> listIterator() {
        return listIterator(0);
    }

    @Override
    public ListIterator<Taxon> listIterator(final int index) {
        return new ListIterator<Taxon>() {
            private final ListIterator<Taxon> i = myTaxaList.listIterator(index);

            public boolean hasNext() {
                return i.hasNext();
            }

            public Taxon next() {
                return i.next();
            }

            public boolean hasPrevious() {
                return i.hasPrevious();
            }

            public Taxon previous() {
                return i.previous();
            }

            public int nextIndex() {
                return i.nextIndex();
            }

            public int previousIndex() {
                return i.previousIndex();
            }

            public void remove() {
                throw new UnsupportedOperationException();
            }

            public void set(Taxon e) {
                throw new UnsupportedOperationException();
            }

            public void add(Taxon e) {
                throw new UnsupportedOperationException();
            }
        };
    }

    @Override
    public List<Taxon> subList(int fromIndex, int toIndex) {
        return Collections.unmodifiableList(myTaxaList.subList(fromIndex, toIndex));
    }
}
