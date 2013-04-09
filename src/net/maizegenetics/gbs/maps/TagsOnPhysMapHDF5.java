/*
 * TagsOnPhysMapHDF5
 */
package net.maizegenetics.gbs.maps;

import ch.systemsx.cisd.hdf5.*;
import java.io.File;
import net.maizegenetics.gbs.tagdist.AbstractTags;
import net.maizegenetics.pal.alignment.Locus;
import org.apache.log4j.Logger;

/**
 *
 * @author edbuckler
 */
public class TagsOnPhysMapHDF5 extends AbstractTags implements TOPMInterface {

    private static final Logger myLogger = Logger.getLogger(TagsOnPhysMapHDF5.class);
    private static final int bitShiftChunk = 16;
    private static final int chunkSize = 1 << bitShiftChunk;
    private final static byte byteMissing = Byte.MIN_VALUE;
    private final static int intMissing = Integer.MIN_VALUE;
    private int maxVariants = 8;
    private int maxMapping = 4;
    private byte[] multimaps;  // number of locations this tagSet maps to; unknown = Byte.MIN_VALUE; multiple, but unknown number = 99
    private IHDF5Writer h5 = null;
    private int cachedMappingIndex = -1;
    private TagMappingInfo cachedTMI = null;
    private int cachedMappingBlock = -1;
    private TagMappingInfo[][] cachedTMIBlock = null;
    private boolean cleanMap = true;
    private boolean cacheAllMappingBlocks = false;
    private HDF5CompoundType<TagMappingInfo> tmiType = null;
    private int[] indicesOfSortByPosition;
    private int tagNum = 0;
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

        myLogger.info("Creating HDF5 File: " + newHDF5file);
        IHDF5WriterConfigurator config = HDF5Factory.configure(new File(newHDF5file));
        config.overwrite();
        // config.dontUseExtendableDataTypes();
        config.useUTF8CharacterEncoding();
        IHDF5Writer h5 = config.writer();

        h5.setIntAttribute("/", "tagCount", inTags.getTagCount());
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
        myLogger.info("Chunk Size for Tags: " + chunkSize);
        int numOfChunks = tagsToChunks(tagCount);
        myLogger.info("Number of Chunks: " + numOfChunks);
        int numTagsPadded = numOfChunks * chunkSize;
        for (int mi = 0; mi < maxMapping; mi++) {
            h5.compounds().createArray("map" + mi, tmiType, numTagsPadded, chunkSize);
            TagMappingInfo[] thTMI = new TagMappingInfo[chunkSize];
            int block = 0;
            for (int i = 0; i < numTagsPadded; i++) {
                if ((mi == 0) && (i < tagCount)) {
                    thTMI[i % chunkSize] = (new TagMappingInfo(inTags.getChromosome(i),
                            inTags.getStrand(i), inTags.getStartPosition(i),
                            inTags.getEndPosition(i), inTags.getDivergence(i)));
                } else {
                    thTMI[i % chunkSize] = new TagMappingInfo();
                }
                if ((i + 1) % chunkSize == 0) {
                    h5.compounds().writeArrayBlock("map" + mi, tmiType, thTMI, block);
                    thTMI = new TagMappingInfo[chunkSize];
                    block++;
                    myLogger.info("Tag locations written: " + (i + 1));
                }
            }
            myLogger.info("...map" + mi + " positions written");
        }

        h5.createByteMatrix("variantDef", maxVariants, 1);
        h5.createByteMatrix("variantPosOff", maxVariants, 1);
        h5.writeByteMatrix("variantDef", inTags.variantDef);
        h5.writeByteMatrix("variantPosOff", inTags.variantPosOff);
        System.out.println("...variant Positions written");
        h5.close();
    }

    private static int tagsToChunks(long tags) {
        return (int) (((tags - 1) >>> bitShiftChunk) + 1);
    }

    public TagsOnPhysMapHDF5(String theHDF5file, boolean cacheAllMappingBlocks) {
        this.cacheAllMappingBlocks = cacheAllMappingBlocks;
        System.out.println("Opening :" + theHDF5file);
        h5 = HDF5Factory.open(theHDF5file);
        this.tagNum = h5.getIntAttribute("/", "tagCount");
        this.tagLengthInLong = h5.getIntAttribute("/", "tagLengthInLong");
        this.tags = h5.readLongMatrix("tags");
        this.tagLength = h5.readByteArray("tagLength");
        this.multimaps = h5.readByteArray("multimaps");
        this.maxMapping = h5.getIntAttribute("/", "maxMapping");
        tmiType = h5.compounds().getInferredType(TagMappingInfo.class);
        cachedTMIBlock = new TagMappingInfo[maxMapping][];
        cacheMappingInfo(0);
        System.out.println(theHDF5file + " read with tags:" + tagNum);
    }

    private void cacheMappingInfo(int index) {
        if (index == cachedMappingIndex) {
            return;
        }
        int block = index >> bitShiftChunk;
        if (cachedMappingBlock != block) {
            if (cleanMap == false) {
                saveCacheBackToFile();
            }
            for (int mi = 0; mi < maxMapping; mi++) {
                cachedTMIBlock[mi] = h5.compounds().readArrayBlock("map0", tmiType, chunkSize, block);
                cachedMappingBlock = block;
                if (cacheAllMappingBlocks == false) {
                    break;
                }
            }
        }
        this.cachedTMI = cachedTMIBlock[0][index % chunkSize];
        this.cachedMappingIndex = index;
    }

    private void saveCacheBackToFile() {
        int block = cachedMappingIndex >> bitShiftChunk;
        if (cachedMappingBlock != block) {
            for (int mi = 0; mi < maxMapping; mi++) {
                if (cleanMap == false) {
                    h5.compounds().writeArrayBlock("map0", tmiType, cachedTMIBlock[mi], block);
                }
                if (cacheAllMappingBlocks == false) {
                    break;
                }
            }
            if (cleanMap == false) {
                h5.writeByteArray("multimaps", multimaps);
            }  //this could be made more efficient by just writing the block
            cleanMap = true;
        }
    }

    public void getFileReadyForClosing() {
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
        return cachedTMIBlock[mapIndex][index % chunkSize];
    }

    public void setAlternateTagMappingInfo(int index, int mapIndex, TagMappingInfo theTMI) {
        if (cacheAllMappingBlocks == false) {
            cacheMappingInfo(index);
        }
        if (cachedMappingIndex != index) {
            cacheMappingInfo(index);
        }
        cachedTMIBlock[mapIndex][index % chunkSize] = theTMI;
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
        TagMappingInfo tempTMI = cachedTMIBlock[mapIndex][index % chunkSize];
        cachedTMIBlock[mapIndex][index % chunkSize] = cachedTMIBlock[mapIndex2][index % chunkSize];
        cachedTMIBlock[mapIndex2][index % chunkSize] = tempTMI;
        cleanMap = false;
    }

    @Override
    public int addVariant(int tagIndex, byte offset, byte base) {
        throw new UnsupportedOperationException("Not supported yet.");
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
        return this.tagNum;
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
    public byte getVariantDef(int tagIndex, int variantIndex) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public byte[] getVariantDefArray(int tagIndex) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public byte getVariantPosOff(int tagIndex, int variantIndex) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public byte[] getVariantPosOffArray(int tagIndex) {
        throw new UnsupportedOperationException("Not supported yet.");
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
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public void setVariantPosOff(int tagIndex, int variantIndex, byte offset) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public static void main(String[] args) {
        String inTOPMFile = "/Users/terry/GBSWorkshopOutput/topm/rice.topm.bin";
        String theTOPMH5 = "/Users/terry/GBSWorkshopOutput/topm/rice.topm.h5";
        //String inTOPMFile = "/Users/edbuckler/SolexaAnal/GBS/build20120701/test/topm100.txt";
        //String theTOPMH5 = "/Users/edbuckler/SolexaAnal/GBS/build20120701/test/topm100.h5";
        //String inTOPMFile = "/Users/edbuckler/SolexaAnal/GBS/build20120701/04_TOPM/Bowtie2/AllZeaMasterTags_c10_20120703.topm";
        //String theTOPMH5 = "/Users/edbuckler/SolexaAnal/GBS/build20120701/test/AllZeaMasterTags_c10_20120703.topm.h5";
        boolean binary = !inTOPMFile.endsWith(".txt");
        TagsOnPhysicalMap inTOPM = new TagsOnPhysicalMap(inTOPMFile, binary);
        TagsOnPhysMapHDF5.createFile(inTOPM, theTOPMH5, 4, 8);
        TagsOnPhysMapHDF5 h5TOPM = new TagsOnPhysMapHDF5(theTOPMH5, false);
        int count = 0, error = 0;
        for (int i = 0; i < inTOPM.getSize(); i += 1000) {
            if (inTOPM.getChromosome(i) > -1) {
                count++;
            }
            if (inTOPM.getChromosome(i) != h5TOPM.getChromosome(i)) {
                error++;
            }
            // if(i!=h5TOPM.getChromosome(i)) error++;
            System.out.printf("%d %d %d %n", i, inTOPM.getChromosome(i), h5TOPM.getChromosome(i));

        }
        System.out.println(count + " error=" + error);
    }
}
