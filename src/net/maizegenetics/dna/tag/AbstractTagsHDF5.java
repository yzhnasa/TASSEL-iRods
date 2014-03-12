/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.dna.tag;

import ch.systemsx.cisd.hdf5.HDF5FloatStorageFeatures;
import ch.systemsx.cisd.hdf5.HDF5GenericStorageFeatures;
import ch.systemsx.cisd.hdf5.HDF5IntStorageFeatures;
import ch.systemsx.cisd.hdf5.IHDF5Writer;
import net.maizegenetics.util.Tassel5HDF5;
import net.maizegenetics.util.Tassel5HDF5Constants;

/**
 * Basic implementations of HDF5 tags. This is designed to annotate tags with a bunch of attributes, to solve the memory issues
 * @author Fei Lu
 */
public abstract class AbstractTagsHDF5 extends AbstractTags implements Tassel5HDF5 {
    protected IHDF5Writer h5 = null;
    protected int currentBlockIndex = -1;
    protected int currentIndex = -1;
    
    @Override
    public int getBlockSize() {
        return Tassel5HDF5Constants.BLOCK_SIZE;
    }

    @Override
    public int getBlockNum() {
        int num = super.getTagCount()/this.getBlockSize();
        if (super.getTagCount()%this.getBlockSize() == 0) return num;
        else return num+1;
    }

    @Override
    public int getBlockIndex(int currentIndex) {
        return currentIndex/this.getBlockSize();
    }

    @Override
    public int getCurrentBlockIndex() {
        return currentBlockIndex;
    }

    @Override
    public int getCurrentIndex() {
        return currentIndex;
    } 

    @Override
    public int getCurrentIndexWithinBlock () {
        return this.getCurrentIndex()%this.getBlockSize();
    }
    
    /**
     * Initialize tag matrix
     * @param tagCount
     * @param tagLengthInLong 
     */
    protected void initializeMatrix (int tagCount, int tagLengthInLong) {
        tagLengthInLong = tagLengthInLong;
        tags = new long[tagLengthInLong][tagCount];
        tagLength = new byte[tagCount];
    }
    
    @Override
    public boolean isInCurrentBlock(int queryIndex) {
        int queryBlockIndex = queryIndex/this.getBlockSize();
        if (queryBlockIndex == currentBlockIndex) return true;
        return false;
    }
    
    @Override
    public HDF5IntStorageFeatures getIntStorageFeatures() {
        return Tassel5HDF5Constants.intDeflation;
    }

    @Override
    public HDF5GenericStorageFeatures getGenericStorageFeatures() {
        return Tassel5HDF5Constants.genDeflation;
    }

    @Override
    public HDF5FloatStorageFeatures getFloatStorageFeatures() {
        return Tassel5HDF5Constants.floatDeflation;
    }
    
    @Override
    public void sort () {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    @Override
    public void swap(int index1, int index2) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
}
