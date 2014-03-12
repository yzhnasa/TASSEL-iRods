/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.util;

import ch.systemsx.cisd.hdf5.HDF5FloatStorageFeatures;
import ch.systemsx.cisd.hdf5.HDF5GenericStorageFeatures;
import ch.systemsx.cisd.hdf5.HDF5IntStorageFeatures;

/**
 * Interface of HDF5 operations
 * @author Fei Lu
 */
public interface Tassel5HDF5 {
    
    /**
     * Returns the size of HDF5 block
     * @return 
     */
    public int getBlockSize();
    
    /**
     * Returns the number of block
     * @return 
     */
    public int getBlockNum ();
    
    /**
     * Returns the block index based on the current record index in a full list
     * @param currentIndex as the record index in the full list
     * @return 
     */
    public int getBlockIndex (int currentIndex);
    
    /**
     * Returns current block index
     * @return 
     */
    public int getCurrentBlockIndex ();
    
    /**
     * Returns current index
     * @return 
     */
    public int getCurrentIndex ();
    
    /**
     * Returns current index relative to current block
     * @return 
     */
    public int getCurrentIndexWithinBlock ();
    
    /**
     * Return boolean value if an index belongs to current block
     * @param queryIndex
     * @return 
     */
    public boolean isInCurrentBlock(int queryIndex);
    
    /**
     * Returns default deflation level of boolean, byte, short, int, long, Enum
     * See {@link net.maizegenetics.util.Tassel5HDF5Constants}
     * @return 
     */
    public HDF5IntStorageFeatures getIntStorageFeatures ();
    
    /**
     * Returns default deflation level of String and Class
     * See {@link net.maizegenetics.util.Tassel5HDF5Constants}
     * @return 
     */
    public HDF5GenericStorageFeatures getGenericStorageFeatures();
    
    /**
     * Returns default deflation level of float and double
     * See {@link net.maizegenetics.util.Tassel5HDF5Constants}
     * @return 
     */
    public HDF5FloatStorageFeatures getFloatStorageFeatures();
    /**
     * Populate a block in memory with default values, update current block index at the same time
     * @param blockIndex
     */
    public void populateBlock (int blockIndex);
    
    /**
     * Initialize a HDF5 file
     * @param hdf5FileS 
     */
    public void initializeHDF5 (String hdf5FileS);
    
    /**
     * Read in HDF5 file
     * @param hdf5FileS 
     */
    public void readHDF5 (String hdf5FileS);
    
    /**
     * Read in a block from HDF5 file
     * @param blockIndex 
     */
    public void readBlock (int blockIndex);
    
    /**
     * Write current block to HDF5 file
     * @param blockIndex 
     */
    public void writeBlock (int blockIndex);
}
