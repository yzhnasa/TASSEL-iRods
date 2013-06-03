package net.maizegenetics.gbs.tagdist;

import cern.colt.Swapper;
import cern.colt.function.IntComparator;
import net.maizegenetics.gbs.tagdist.TagsByTaxa.FilePacking;

/**
 * Basic interface for holding sets of sequence tag (these are compressed into 2-bit codings
 * of tags).  The length of tag in good sequence is also tracked.
 * @author Fei Lu
 */
public interface PETags extends Swapper, IntComparator {

    /**
     * Returns the number of long use to represent the sequence
     * @return
     */
    public int getTagSizeInLong();

    /**
     * Returns the length in bp of a particular forward tag
     * @param index
     * @return length
     */
    public short getTagFLength(int index);
    
    /**
     * Returns the length in bp of a particular backward tag
     * @param index
     * @return length
     */
    public short getTagBLength(int index);

    /**
     * Get the compressed forward read sequence in a long array for a given index
     * @param index
     * @return compressed read sequence in long array
     */
    public long[] getTagF(int index);
    
    /**
     * Get the compressed backward read sequence in a long array for a given index
     * @param index
     * @return compressed read sequence in long array
     */
    public long[] getTagB(int index);

    /**
     * Gets the first index of a read (the only one if a unique list).
     * If the read is not found then it return
     * a negative value indicating its insertion point.
     * @param read as a compressed long array
     * @return index of the read in the array
     */
    public int getTagIndex(long[] tagF, long[] tagB);


    /**
     * This is the number of different tags in the list (NOT THE SUM OF THE COUNTS)
     * The index will vary from 0 to (ReadTotal-1)
     * This is the number of distinct tags if readUnique is true
     * @return total number of tags
     */
    public int getTagCount();
    
    /*
     * @return contig of a PE tag, if there is no contig, return null
     */
    public long[] getContig (int index);
    
    /*
     * @return contig length
     */
    public short getContigLength (int index);
    
    /*
     * @return contig lengthInLong
     */
    public byte getContigLengthInLong (int index);
    
    /**
     * Read PE tag distribution
     * @param infileS is the file name
     * @param format is the file format in FilePacking
     */
    
    public void readDistFile (String infileS, FilePacking format);
    
    /**
     * Write PE tag distribution
     * @param outfileS is the file name
     * @param format is the file format in FilePacking
     * @param minCount is the minimum count of a PE tag
     */
    public void writeDistFile (String outfileS, FilePacking format, int minCount);
}
