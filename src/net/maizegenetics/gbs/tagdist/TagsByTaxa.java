/*
 * TagsByTaxa
 */
package net.maizegenetics.gbs.tagdist;

import net.maizegenetics.util.OpenBitSet;

import java.io.File;

/**
 * Tags by Taxa is the interface for tracking the distribution of tags across
 * taxa.  It can be thought of as a matrix of nTaxa by nTags.  Taxa names are
 * also maintained by this interface.
 *
 *
 * @author edbuckler
 */
public interface TagsByTaxa extends Tags {

    /**
     * Types of file packing used by GBS data structures.  Bit compresses all information
     * to presence/absence, while Byte/Short/Int tradeoff memory size with depth of coverage.
     * Text is only used for file storage.
     */
    public static enum FilePacking {

        Bit, Byte, Short, Int, Text
    };

    void addReadsToTagTaxon(int tagIndex, int taxaIndex, int addValue);

    int getIndexOfTaxaName(String taxon);

    int getReadCountForTagTaxon(int tagIndex, int taxaIndex);

    byte[] getTaxaReadCountsForTag(int readIndex);

    public void truncateTaxonNames();

    public void setMethodByRows(boolean rowSetMethod);

    /**
     * The presence/absence of the taxa in a bitSet format
     * @param readIndex
     * @return distribution of taxa
     */
    OpenBitSet getTaxaReadBitsForTag(int readIndex);

    int getTaxaCount();

    /**
     * Get the number of taxa with a given tag 
     * It is count of the taxa with readCount>0
     * @param readIndex
     * @return number of taxa with a read
     */
    int getNumberOfTaxaWithTag(int readIndex);

    String getTaxaName(int taxaIndex);

    String[] getTaxaNames();

    void setReadCountForTagTaxon(int tagIndex, int taxaIndex, int value);

    void writeDistFile(File outFile, FilePacking fileType, int minCount);

    void initMatrices(int taxaNum, int tagNum);

    /**
     * Add taxa to the TagsByTaxa matrix, they will all be set to distribution value of
     * zero.
     * @param addTaxaNames
     */
    void addTaxa(String[] addTaxaNames);

    /**In implementations that use a RandomAccessFile for storage, this clears the RAM buffer of any
    remaining data, writes it to the file on disk, and closes the file.*/
    public void getFileReadyForClosing();
}
