/*
 * AbstractTags
 */
package net.maizegenetics.gbs.tagdist;

import cern.colt.GenericSorting;
import java.util.Arrays;

/**
 * Basic methods for working with PE Tags, including sorting and search.
 * @author Fei Lu
 */
public abstract class AbstractPETags implements PETags {

    protected int tagLengthInLong;  //TODO fully implement on reading
    //tagLengthInLong should be much larger than normal e.g. 8
    protected long[][] tagsF;  // for memory efficiency the rows first and second half of the read
                               // columns are the index of the reads.
    protected long[][] tagsB;
    protected short[] tagFLength;  // length of tag (number of bases)  // 1 byte
    protected short[] tagBLength;
    protected long[][] contig;
    protected short[] f2blength;  //entire length of PE tags (cutsite to cutsite), negative could be use to indicate minimum when not overlapping
    
    
    protected String[] taxaName; //TODO REMOVE

   
    @Override  //TODO REMOVE
    public int getTaxaCount () {
        return taxaName.length;
    }
    
    @Override  //TODO REMOVE
    public String getTaxaName (int index) {
        return taxaName[index];
    }
    
    @Override
    public long[] getTagF(int index) {
        long[] theTag = new long[tagLengthInLong];
        for (int i = 0; i < tagLengthInLong; i++) {
            theTag[i] = tagsF[i][index];
        }
        return theTag;
    }
    
    @Override
    public long[] getTagB(int index) {
        long[] theTag = new long[tagLengthInLong];
        for (int i = 0; i < tagLengthInLong; i++) {
            theTag[i] = tagsB[i][index];
        }
        return theTag;
    }

    @Override
    public int getTagCount() {
        return tagsF[0].length;
    }
    
    @Override
    public int getTagIndex(long[] tagF, long[] tagB) {
        //code inspired by COLT lower bound function
        int first = 0;
        int len = tagsF[0].length - first;
        int comp = 0;
        while (len > 0) {
            int half = len / 2;
            int middle = first + half;
            if ((comp = compareTags(middle, tagF, tagB)) < 0) {
                first = middle + 1;
                len -= half + 1;
            } else {
                len = half;
            }
        }
        if ((first < tagsF[0].length) && (compareTags(first, tagF, tagB) == 0)) {
            return first;
        }
        return -(first + 1);
    }

    public void orderTagsFTagsB () {
        for (int i = 0; i < this.getTagCount(); i++) {
            this.orderTagFTagB(i);
        }
    }
    
    private boolean orderTagFTagB (int index) {
        if (this.compareTagFTagB(index) == 1) {
            this.switchTagFTagB(index);
            return true;
        }
        return false;
    }
    
    private void switchTagFTagB (int index) {
        for (int i = 0; i < tagLengthInLong; i++) {
            long temp = tagsF[i][index];
            tagsF[i][index] = tagsB[i][index];
            tagsB[i][index] = temp;
        }
    }
    
    private int compareTagFTagB (int index) {
        for (int i = 0; i < tagLengthInLong; i++) {
            if (tagsF[i][index] < tagsB[i][index]) return -1;
            if (tagsF[i][index] > tagsB[i][index]) return 1;
        }
        return 0;
    }
    
    private int compareTags(int index1, long[] tagF, long[] tagB) {
        for (int i = 0; i < tagLengthInLong; i++) {
            if (tagsF[i][index1] < tagF[i]) {
                return -1;
            }
            if (tagsF[i][index1] > tagF[i]) {
                return 1;
            }
        }
        for (int i = 0; i < tagLengthInLong; i++) {
            if (tagsB[i][index1] < tagB[i]) {
                return -1;
            }
            if (tagsB[i][index1] > tagB[i]) {
                return 1;
            }
        }
        return 0;
    }

    
    @Override
    public short getTagFLength(int index) {
        return tagFLength[index];
    }
    
    @Override
    public short getTagBLength(int index) {
        return tagBLength[index];
    }

    @Override
    public int getTagSizeInLong() {
        return tagLengthInLong;
    }

    protected void iniMatrix (int tagLengthInLong, int tagNum, int taxaNum) {
        this.tagLengthInLong = tagLengthInLong;
        tagsF = new long[tagLengthInLong][tagNum];
        tagsB = new long[tagLengthInLong][tagNum];
        tagFLength = new short[tagNum];
        tagBLength = new short[tagNum];
        taxaName = new String[taxaNum];
    }
    
    @Override
    public void swap(int index1, int index2) {
        long temp;
        for (int i = 0; i < tagLengthInLong; i++) {
            temp = tagsF[i][index1];
            tagsF[i][index1] = tagsF[i][index2];
            tagsF[i][index2] = temp;
            temp = tagsB[i][index1];
            tagsB[i][index1] = tagsB[i][index2];
            tagsB[i][index2] = temp;
        }
        short tl;
        tl = tagFLength[index1];
        tagFLength[index1] = tagFLength[index2];
        tagFLength[index2] = tl;
        tl = tagBLength[index1];
        tagBLength[index1] = tagBLength[index2];
        tagBLength[index2] = tl;
    }

    @Override
    public int compare(int index1, int index2) {
        for (int i = 0; i < tagLengthInLong; i++) {
            if (tagsF[i][index1] < tagsF[i][index2]) {
                return -1;
            }
            if (tagsF[i][index1] > tagsF[i][index2]) {
                return 1;
            }
        }
        for (int i = 0; i < tagLengthInLong; i++) {
            if (tagsB[i][index1] < tagsB[i][index2]) {
                return -1;
            }
            if (tagsB[i][index1] > tagsB[i][index2]) {
                return 1;
            }
        }
        return 0;
    }

    public void sort() {
        System.out.println("Position index sort begin.");
        GenericSorting.quickSort(0, this.getTagCount(), this, this);
        System.out.println("Position index sort end.");
    }
    
}
