package net.maizegenetics.pal.site;

import net.maizegenetics.pal.alignment.Locus;

import java.util.*;

/**
 * Created with IntelliJ IDEA.
 * User: edbuckler
 * Date: 7/30/13
 * Time: 5:22 PM
 */
public class SimplePositionList implements PositionList {

    private final List<Position> mySiteList = new ArrayList<Position>();

    @Override
    public byte getReferenceAllele(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
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
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String getSNPID(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int getSiteCount() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int getLocusSiteCount(Locus locus) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int[] getStartAndEndOfLocus(Locus locus) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int getPositionInLocus(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int getSiteOfPhysicalPosition(int physicalPosition, Locus locus) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int getSiteOfPhysicalPosition(int physicalPosition, Locus locus, String snpID) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int[] getPhysicalPositions() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public byte getPositionType(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public byte[] getPositionTypes() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String getLocusName(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public Locus getLocus(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public Locus getLocus(String name) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public Locus[] getLoci() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int getNumLoci() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int[] getLociOffsets() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int getIndelSize(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public boolean isIndel(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public byte getMajorAllele(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String getMajorAlleleAsString(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public byte getMinorAllele(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String getMinorAlleleAsString(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public byte[] getMinorAlleles(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public byte[] getAlleles(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double getMinorAlleleFrequency(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public double getMajorAlleleFrequency(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String getGenomeAssembly() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public boolean isPositiveStrand(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    @Override
    public int[][] getAllelesSortedByFrequency(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public Object[][] getDiploidssSortedByFrequency(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
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

    @Override
    public boolean add(Position e) {
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
    public boolean addAll(Collection<? extends Position> c) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    @Override
    public boolean addAll(int index, Collection<? extends Position> c) {
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
    public Position get(int index) {
        return mySiteList.get(index);
    }

    @Override
    public Position set(int index, Position element) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    @Override
    public void add(int index, Position element) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    @Override
    public Position remove(int index) {
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
    public ListIterator<Position> listIterator() {
        return mySiteList.listIterator();
    }

    @Override
    public ListIterator<Position> listIterator(int index) {
        return mySiteList.listIterator(index);
    }

    @Override
    public List<Position> subList(int fromIndex, int toIndex) {
        return mySiteList.subList(fromIndex, toIndex);
    }

}
