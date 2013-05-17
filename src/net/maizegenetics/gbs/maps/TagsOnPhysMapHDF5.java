/*
 * TagsOnPhysMapHDF5
 */
package net.maizegenetics.gbs.maps;

import ch.systemsx.cisd.hdf5.*;
import java.io.File;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;
import net.maizegenetics.gbs.tagdist.AbstractTags;
import net.maizegenetics.pal.alignment.Locus;
import org.apache.log4j.Logger;

/**
 *
 * @author Ed, Terry
 */
public class TagsOnPhysMapHDF5 extends AbstractTags implements TOPMInterface {

    private static final Logger myLogger = Logger.getLogger(TagsOnPhysMapHDF5.class);
    private static final int NUM_UNITS_TO_CACHE_ON_GET = 64;
    private static final int MAX_CACHE_SIZE = 20;
    private static final int BITS_TO_SHIFT_FOR_CHUNK = 16;
    private static final int CHUNK_SIZE = 1 << BITS_TO_SHIFT_FOR_CHUNK;
    private final static byte BYTE_MISSING = Byte.MIN_VALUE;
    private final static int intMissing = Integer.MIN_VALUE;
    private int myMaxVariants = 8;
    private final Map<Integer, byte[][]> myCachedVariantDefs = new LinkedHashMap<Integer, byte[][]>((MAX_CACHE_SIZE * 3) / 2) {
        @Override
        protected boolean removeEldestEntry(Map.Entry eldest) {
            boolean result = size() > MAX_CACHE_SIZE;
            if (result) {
                writeCachedVariantDefs();
            }
            return result;
        }
    };
    private final Map<Integer, byte[][]> myCachedVariantOffsets = new LinkedHashMap<Integer, byte[][]>((MAX_CACHE_SIZE * 3) / 2) {
        @Override
        protected boolean removeEldestEntry(Map.Entry eldest) {
            boolean result = size() > MAX_CACHE_SIZE;
            if (result) {
                writeCachedVariantOffsets();
            }
            return result;
        }
    };
    private Set<Integer> myDirtyVariantDefs = new TreeSet<Integer>();
    private Set<Integer> myDirtyVariantOffsets = new TreeSet<Integer>();
    private int maxMapping = 4;
    private byte[] multimaps;  // number of locations this tagSet maps to; unknown = Byte.MIN_VALUE; multiple, but unknown number = 99
    private IHDF5Writer myHDF5 = null;
    private int cachedMappingIndex = -1;
    private TagMappingInfo cachedTMI = null;
    private int cachedMappingBlock = -1;
    private TagMappingInfo[][] cachedTMIBlock = null;
    private boolean cleanMap = true;
    private boolean cacheAllMappingBlocks = false;
    private HDF5CompoundType<TagMappingInfo> tmiType = null;
    private int[] indicesOfSortByPosition;
    private int myNumTags = 0;
    private boolean redundantTags = true;  // this field has not been utilized yet

    public static void createFile(TagsOnPhysicalMap inTags, String newHDF5file, int maxMapping, int maxVariants) {

        int tagLengthInLong = inTags.getTagSizeInLong();
        int tagCount = inTags.getTagCount();
        long[][] tags = new long[tagLengthInLong][tagCount];
        byte[] tagLength = new byte[tagCount];
        for (int i = 0; i < tagCount; i++) {
            long[] ct = inTags.getTag(i);
            for (int j = 0; j < tagLengthInLong; j++) {
                tags[j][i] = ct[j];
            }
            tagLength[i] = (byte) inTags.getTagLength(i);
        }

        IHDF5Writer h5 = null;
        try {
            myLogger.info("Creating HDF5 File: " + newHDF5file);
            IHDF5WriterConfigurator config = HDF5Factory.configure(new File(newHDF5file));
            config.overwrite();
            // config.dontUseExtendableDataTypes();
            config.useUTF8CharacterEncoding();
            h5 = config.writer();

            h5.setIntAttribute("/", "tagCount", tagCount);
            h5.setIntAttribute("/", "maxVariants", maxVariants);
            h5.setIntAttribute("/", "maxMapping", maxMapping);
            h5.setIntAttribute("/", "tagLengthInLong", tagLengthInLong);
            //create tag matrix
            h5.createLongMatrix("tags", inTags.getTagSizeInLong(), tagCount, inTags.getTagSizeInLong(), tagCount);
            h5.writeLongMatrix("tags", tags);
            System.out.println("...Tags written");
            h5.createByteArray("tagLength", tagCount);
            h5.writeByteArray("tagLength", tagLength);
            System.out.println("...Tags lengths written");
            h5.createByteArray("multimaps", tagCount);
            h5.writeByteArray("multimaps", inTags.multimaps);
            System.out.println("...multimapping  written");
            HDF5CompoundType<TagMappingInfo> tmiType = h5.compounds().getInferredType(TagMappingInfo.class);
            myLogger.info("Chunk Size for Tags: " + CHUNK_SIZE);
            int numOfChunks = tagsToChunks(tagCount);
            myLogger.info("Number of Chunks: " + numOfChunks);
            int numTagsPadded = numOfChunks * CHUNK_SIZE;
            for (int mi = 0; mi < maxMapping; mi++) {
                h5.compounds().createArray("map" + mi, tmiType, numTagsPadded, CHUNK_SIZE);
                TagMappingInfo[] thTMI = new TagMappingInfo[CHUNK_SIZE];
                int block = 0;
                for (int i = 0; i < numTagsPadded; i++) {
                    if ((mi == 0) && (i < tagCount)) {
                        thTMI[i % CHUNK_SIZE] = (new TagMappingInfo(inTags.getChromosome(i),
                                inTags.getStrand(i), inTags.getStartPosition(i),
                                inTags.getEndPosition(i), inTags.getDivergence(i)));
                    } else {
                        thTMI[i % CHUNK_SIZE] = new TagMappingInfo();
                    }
                    if ((i + 1) % CHUNK_SIZE == 0) {
                        h5.compounds().writeArrayBlock("map" + mi, tmiType, thTMI, block);
                        thTMI = new TagMappingInfo[CHUNK_SIZE];
                        block++;
                        myLogger.info("Tag locations written: " + (i + 1));
                    }
                }
                myLogger.info("...map" + mi + " positions written");
            }

            h5.createByteMatrix("variantDef", NUM_UNITS_TO_CACHE_ON_GET, maxVariants);
            byte[][] variantDef = inTags.getVariantDef();
            int numVariants = variantDef[0].length;
            if (numVariants > maxVariants) {
                throw new IllegalArgumentException("TagsOnPhysMapHDF5: createFile: max variants can't be less than original TOPM Variant Defs: " + numVariants);
            }
            if (variantDef.length != tagCount) {
                throw new IllegalStateException("TagsOnPhysMapHDF5: createFile: variant def dimension: " + variantDef.length + " doesn't match tag count: " + tagCount);
            }
            byte[][] temp = new byte[tagCount][maxVariants];
            for (int i = 0; i < tagCount; i++) {
                for (int j = 0; j < maxVariants; j++) {
                    if (j < numVariants) {
                        temp[i][j] = variantDef[i][j];
                    } else {
                        temp[i][j] = BYTE_MISSING;
                    }
                }
            }
            h5.writeByteMatrix("variantDef", temp);

            h5.createByteMatrix("variantPosOff", NUM_UNITS_TO_CACHE_ON_GET, maxVariants);
            byte[][] variantPosOff = inTags.getVariantOff();
            numVariants = variantPosOff[0].length;
            if (numVariants > maxVariants) {
                throw new IllegalArgumentException("TagsOnPhysMapHDF5: createFile: max variants can't be less than original TOPM Variant Offsets: " + numVariants);
            }
            if (variantPosOff.length != tagCount) {
                throw new IllegalStateException("TagsOnPhysMapHDF5: createFile: variant def offset dimension: " + variantPosOff.length + " doesn't match tag count: " + tagCount);
            }
            temp = new byte[tagCount][maxVariants];
            for (int i = 0; i < tagCount; i++) {
                for (int j = 0; j < maxVariants; j++) {
                    if (j < numVariants) {
                        temp[i][j] = variantPosOff[i][j];
                    } else {
                        temp[i][j] = BYTE_MISSING;
                    }
                }
            }
            h5.writeByteMatrix("variantPosOff", temp);

        } finally {
            try {
                h5.close();
            } catch (Exception e) {
                // do nothing
            }
        }
    }

    private static int tagsToChunks(long tags) {
        return (int) (((tags - 1) >>> BITS_TO_SHIFT_FOR_CHUNK) + 1);
    }

    public TagsOnPhysMapHDF5(String filename) {
        myHDF5 = HDF5Factory.open(filename);
        myMaxVariants = myHDF5.getIntAttribute("/", "maxVariants");
        myNumTags = myHDF5.getIntAttribute("/", "tagCount");
    }

    public TagsOnPhysMapHDF5(String theHDF5file, boolean cacheAllMappingBlocks) {

        this.cacheAllMappingBlocks = cacheAllMappingBlocks;
        System.out.println("Opening :" + theHDF5file);
        myHDF5 = HDF5Factory.open(theHDF5file);
        myNumTags = myHDF5.getIntAttribute("/", "tagCount");
        this.tagLengthInLong = myHDF5.getIntAttribute("/", "tagLengthInLong");
        this.tags = myHDF5.readLongMatrix("tags");
        this.tagLength = myHDF5.readByteArray("tagLength");
        this.multimaps = myHDF5.readByteArray("multimaps");
        this.maxMapping = myHDF5.getIntAttribute("/", "maxMapping");
        tmiType = myHDF5.compounds().getInferredType(TagMappingInfo.class);
        cachedTMIBlock = new TagMappingInfo[maxMapping][];
        cacheMappingInfo(0);
        System.out.println(theHDF5file + " read with tags:" + myNumTags);

    }

    private byte[] getCachedVariantDefs(int tag) {

        int startTag = (tag / NUM_UNITS_TO_CACHE_ON_GET) * NUM_UNITS_TO_CACHE_ON_GET;

        byte[][] result = myCachedVariantDefs.get(startTag);
        if (result != null) {
            return result[tag - startTag];
        }

        int numToRetrieve = Math.min(NUM_UNITS_TO_CACHE_ON_GET, myNumTags - startTag);

        synchronized (myHDF5) {
            result = myHDF5.readByteMatrixBlockWithOffset("variantDef", numToRetrieve, myMaxVariants, startTag, 0);
        }

        myCachedVariantDefs.put(startTag, result);

        return result[tag - startTag];

    }

    private byte[] getCachedVariantOffsets(int tag) {

        int startTag = (tag / NUM_UNITS_TO_CACHE_ON_GET) * NUM_UNITS_TO_CACHE_ON_GET;

        byte[][] result = myCachedVariantOffsets.get(startTag);
        if (result != null) {
            return result[tag - startTag];
        }

        int numToRetrieve = Math.min(NUM_UNITS_TO_CACHE_ON_GET, myNumTags - startTag);

        synchronized (myHDF5) {
            result = myHDF5.readByteMatrixBlockWithOffset("variantPosOff", numToRetrieve, myMaxVariants, startTag, 0);
        }

        myCachedVariantOffsets.put(startTag, result);

        return result[tag - startTag];

    }

    private void writeCachedVariantDefs() {
        synchronized (myDirtyVariantDefs) {
            if (myDirtyVariantDefs.size() == 0) {
                return;
            }

            Iterator<Integer> itr = myDirtyVariantDefs.iterator();
            while (itr.hasNext()) {
                int startTag = itr.next();
                byte[][] current = myCachedVariantDefs.get(startTag);
                synchronized (myHDF5) {
                    myHDF5.writeByteMatrixBlockWithOffset("variantDef", current, startTag, 0);
                }
            }

            myDirtyVariantDefs.clear();
        }
    }

    private void writeCachedVariantOffsets() {
        synchronized (myDirtyVariantOffsets) {
            if (myDirtyVariantOffsets.size() == 0) {
                return;
            }

            Iterator<Integer> itr = myDirtyVariantOffsets.iterator();
            while (itr.hasNext()) {
                int startTag = itr.next();
                byte[][] current = myCachedVariantOffsets.get(startTag);
                synchronized (myHDF5) {
                    myHDF5.writeByteMatrixBlockWithOffset("variantPosOff", current, startTag, 0);
                }
            }

            myDirtyVariantOffsets.clear();
        }
    }

    private void cacheMappingInfo(int index) {
        if (index == cachedMappingIndex) {
            return;
        }
        int block = index >> BITS_TO_SHIFT_FOR_CHUNK;
        if (cachedMappingBlock != block) {
            if (cleanMap == false) {
                saveCacheBackToFile();
            }
            for (int mi = 0; mi < maxMapping; mi++) {
                cachedTMIBlock[mi] = myHDF5.compounds().readArrayBlock("map0", tmiType, CHUNK_SIZE, block);
                cachedMappingBlock = block;
                if (cacheAllMappingBlocks == false) {
                    break;
                }
            }
        }
        this.cachedTMI = cachedTMIBlock[0][index % CHUNK_SIZE];
        this.cachedMappingIndex = index;
    }

    private void saveCacheBackToFile() {
        int block = cachedMappingIndex >> BITS_TO_SHIFT_FOR_CHUNK;
        if (cachedMappingBlock != block) {
            for (int mi = 0; mi < maxMapping; mi++) {
                if (cleanMap == false) {
                    myHDF5.compounds().writeArrayBlock("map0", tmiType, cachedTMIBlock[mi], block);
                }
                if (cacheAllMappingBlocks == false) {
                    break;
                }
            }
            if (cleanMap == false) {
                myHDF5.writeByteArray("multimaps", multimaps);
            }  //this could be made more efficient by just writing the block
            cleanMap = true;
        }
    }

    public void getFileReadyForClosing() {
        writeCachedVariantDefs();
        writeCachedVariantOffsets();
        saveCacheBackToFile();
    }

    //    private void cacheMappingInfo(int index) {
    //        TagMappingInfo[] aTMI=h5.compounds().readArrayBlockWithOffset("map0", tmiType, 1, index);
    //        this.cachedTMI=aTMI[0];
    //        this.cachedMappingIndex=index;
    //    }
    public TagMappingInfo getAlternateTagMappingInfo(int index, int mapIndex) {
        if (cacheAllMappingBlocks == false) {
            cacheMappingInfo(index);
        }
        if (cachedMappingIndex != index) {
            cacheMappingInfo(index);
        }
        return cachedTMIBlock[mapIndex][index % CHUNK_SIZE];
    }

    public void setAlternateTagMappingInfo(int index, int mapIndex, TagMappingInfo theTMI) {
        if (cacheAllMappingBlocks == false) {
            cacheMappingInfo(index);
        }
        if (cachedMappingIndex != index) {
            cacheMappingInfo(index);
        }
        cachedTMIBlock[mapIndex][index % CHUNK_SIZE] = theTMI;
        if (multimaps[index] >= mapIndex) {
            multimaps[index] = (byte) (mapIndex + 1);
        }
        cleanMap = false;
    }

    public void swapTagMappingInfo(int index, int mapIndex, int mapIndex2) {
        if (cacheAllMappingBlocks == false) {
            cacheMappingInfo(index);
        }
        if (cachedMappingIndex != index) {
            cacheMappingInfo(index);
        }
        TagMappingInfo tempTMI = cachedTMIBlock[mapIndex][index % CHUNK_SIZE];
        cachedTMIBlock[mapIndex][index % CHUNK_SIZE] = cachedTMIBlock[mapIndex2][index % CHUNK_SIZE];
        cachedTMIBlock[mapIndex2][index % CHUNK_SIZE] = tempTMI;
        cleanMap = false;
    }

    @Override
    public int addVariant(int tagIndex, byte offset, byte base) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    @Override
    public byte getByteMissing() {
        return BYTE_MISSING;
    }

    @Override
    public int getChromosome(int index) {
        if (cachedMappingIndex != index) {
            cacheMappingInfo(index);
        }
        return cachedTMI.chromosome;
    }

    @Override
    public int[] getChromosomes() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public byte getDcoP(int index) {
        if (cachedMappingIndex != index) {
            cacheMappingInfo(index);
        }
        return cachedTMI.dcoP;
    }

    @Override
    public byte getDivergence(int index) {
        if (cachedMappingIndex != index) {
            cacheMappingInfo(index);
        }
        return cachedTMI.divergence;
    }

    @Override
    public int getEndPosition(int index) {
        if (cachedMappingIndex != index) {
            cacheMappingInfo(index);
        }
        return cachedTMI.divergence;
    }

    @Override
    public Locus[] getLoci() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public Locus getLocus(int tagIndex) {
        if (cachedMappingIndex != tagIndex) {
            cacheMappingInfo(tagIndex);
        }
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public byte getMapP(int index) {
        if (cachedMappingIndex != index) {
            cacheMappingInfo(index);
        }
        return cachedTMI.mapP;
    }

    @Override
    public byte getMultiMaps(int index) {
        return this.multimaps[index];
    }

    @Override
    public int[] getPositionArray(int index) {
        if (cachedMappingIndex != index) {
            cacheMappingInfo(index);
        }
        int[] r = {cachedTMI.chromosome, cachedTMI.strand, cachedTMI.startPosition};
        return r;
    }

    @Override
    public int getReadIndexForPositionIndex(int posIndex) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int getSize() {
        return myNumTags;
    }

    @Override
    public int getStartPosition(int index) {
        if (cachedMappingIndex != index) {
            cacheMappingInfo(index);
        }
        return cachedTMI.startPosition;
    }

    @Override
    public byte getStrand(int tagIndex) {
        if (cachedMappingIndex != tagIndex) {
            cacheMappingInfo(tagIndex);
        }
        return cachedTMI.strand;
    }

    @Override
    public int[] getUniquePositions(int chromosome) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
    
    @Override
    public byte getVariantDef(int tagIndex, int variantIndex) {
        return getCachedVariantDefs(tagIndex)[variantIndex];
    }

    @Override
    public byte[] getVariantDefArray(int tagIndex) {
        return getCachedVariantDefs(tagIndex);
    }

    @Override
    public byte getVariantPosOff(int tagIndex, int variantIndex) {
        return getCachedVariantOffsets(tagIndex)[variantIndex];
    }

    @Override
    public byte[] getVariantPosOffArray(int tagIndex) {
        return getCachedVariantOffsets(tagIndex);
    }

    public void setMultimaps(int index, byte multimaps) {
        this.multimaps[index] = multimaps;
    }

    @Override
    public void setChromoPosition(int index, int chromosome, byte strand, int positionMin, int positionMax) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public void setDivergence(int index, byte divergence) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public void setMapP(int index, byte mapP) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public void setMapP(int index, double mapP) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public void setVariantDef(int tagIndex, int variantIndex, byte def) {
        synchronized (myDirtyVariantDefs) {
            int startTag = (tagIndex / NUM_UNITS_TO_CACHE_ON_GET) * NUM_UNITS_TO_CACHE_ON_GET;
            getCachedVariantDefs(tagIndex)[variantIndex] = def;
            myDirtyVariantDefs.add(startTag);
        }
    }

    @Override
    public void setVariantPosOff(int tagIndex, int variantIndex, byte offset) {
        synchronized (myDirtyVariantOffsets) {
            int startTag = (tagIndex / NUM_UNITS_TO_CACHE_ON_GET) * NUM_UNITS_TO_CACHE_ON_GET;
            getCachedVariantOffsets(tagIndex)[variantIndex] = offset;
            myDirtyVariantOffsets.add(startTag);
        }
    }

    @Override
    public int getMaxNumVariants() {
        return myMaxVariants;
    }

    public static void main(String[] args) {
        String filename = "/Volumes/ThunderboltSSD/TOPMs_Merging20130114/topm_100_lines.topm.h5";
        TagsOnPhysMapHDF5 topm = new TagsOnPhysMapHDF5(filename);
        int tagCount = topm.getSize();
        System.out.println("tag count: " + tagCount);
        int maxVariants = topm.getMaxNumVariants();
        System.out.println("max variants: " + maxVariants);
        for (int i = 0; i < tagCount; i++) {
            for (int j = 0; j < maxVariants; j++) {
                byte offset = topm.getVariantPosOff(i, j);
                byte def = topm.getVariantDef(i, j);
                System.out.println("tag: " + i + "  variant: " + j + "  offset: " + offset + "  def: " + def);
                if (offset != BYTE_MISSING) {
                    topm.setVariantPosOff(i, j, (byte) (offset + 1));
                    //topm.setVariantPosOff(i, j, (byte) i);
                    topm.setVariantDef(i, j, (byte) i);
                }
            }
        }
        topm.getFileReadyForClosing();
    }

    @Override
    public byte[][] getVariantOff() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public byte[][] getVariantDef() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public void clearVariants() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
}
