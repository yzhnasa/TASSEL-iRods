/*
 * TagsOnPhysMapHDF5
 */
package net.maizegenetics.gbs.maps;

import ch.systemsx.cisd.hdf5.*;
import java.io.File;
import java.util.Arrays;
import java.util.Random;
import java.util.TreeMap;
import java.util.TreeSet;
import net.maizegenetics.gbs.tagdist.AbstractTags;
import net.maizegenetics.pal.alignment.Locus;
import org.apache.log4j.Logger;

/**
 *
 * TODO:
 * <li> createFile - needs to deal with inputing just TAGS or TOPM with compression
 * <li> create indicesOfSortByPosition, myChromosomes, myUniquePositions to work
 * 
 * 
 * @author Ed, Terry
 */
public class TagsOnPhysMapHDF5 extends AbstractTags implements TOPMInterface {

    private static final Logger myLogger = Logger.getLogger(TagsOnPhysMapHDF5.class);
    private static final int NUM_UNITS_TO_CACHE_ON_GET = 64;
    private static final int MAX_CACHE_SIZE = 20;
    private static final int BITS_TO_SHIFT_FOR_CHUNK = 16;
    private static final int CHUNK_SIZE = 1 << BITS_TO_SHIFT_FOR_CHUNK;
    private static HDF5GenericStorageFeatures genoFeatures = HDF5GenericStorageFeatures.createDeflation(HDF5GenericStorageFeatures.MAX_DEFLATION_LEVEL); //used by mapping object
    private static HDF5IntStorageFeatures vectorFeatures = HDF5IntStorageFeatures.createDeflation(HDF5IntStorageFeatures.DEFAULT_DEFLATION_LEVEL); //used by vectors
    private int myMaxVariants = 16;
    private int maxMapping = 4;
    private byte[] multimaps;  // number of locations this tagSet maps to; unknown = Byte.MIN_VALUE; multiple, but unknown number = 99
    private byte[] bestStrand;
    private int[] bestChr;
    private int[] bestStartPos;
    private byte[][] variantDefs;
    private byte[][] variantOffsets;
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
    private boolean hasDetailedMapping=false;
    private boolean redundantTags = true;  // this field has not been utilized yet
    private int[] myChromosomes=null;  //sort ascending
    private int[][] myUniquePositions=null;  //dimensions [chromosome][positions] note ragged, and sort ascending

    public static void createFile(TOPMInterface inTags, String newHDF5file, int maxMapping, int maxVariants) {

        int tagLengthInLong = inTags.getTagSizeInLong();
        int tagCount = inTags.getTagCount();
        tagCount=tagCount/200;
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
            h5.createLongMatrix("tags", inTags.getTagSizeInLong(), tagCount, inTags.getTagSizeInLong(), tagCount, vectorFeatures);
            h5.writeLongMatrix("tags", tags, vectorFeatures);
            System.out.println("...Tags written");
            h5.createByteArray("tagLength", tagCount, vectorFeatures);
            h5.writeByteArray("tagLength", tagLength, vectorFeatures);
            System.out.println("...Tags lengths written");
            
            //need to create branch between instance of Tags and TOPM
            
            
            //creating fast bestChr and bestPosition
            byte[] mmOut=new byte[tagCount];
            byte[] strandOut=new byte[tagCount];
            int[] chrOut=new int[tagCount];
            int[] posOut=new int[tagCount];
            
            for (int i = 0; i < tagCount; i++) {
                mmOut[i]=inTags.getMultiMaps(i);
                strandOut[i]=inTags.getStrand(i);
                chrOut[i]=inTags.getChromosome(i);
                posOut[i]=inTags.getStartPosition(i);
            }
            h5.createByteArray("multimaps", tagCount);
            h5.writeByteArray("multimaps", mmOut, vectorFeatures);
            h5.createByteArray("bestStrand", tagCount);
            h5.writeByteArray("bestStrand", strandOut, vectorFeatures);
            h5.createIntArray("bestChr", tagCount, vectorFeatures);
            h5.writeIntArray("bestChr", chrOut, vectorFeatures);
            h5.createIntArray("bestStartPos", tagCount, vectorFeatures);
            h5.writeIntArray("bestStartPos", posOut, vectorFeatures);
            System.out.println("...multimapping, strand, chr, position  written");
                        
            HDF5CompoundType<TagMappingInfo> tmiType = h5.compounds().getInferredType(TagMappingInfo.class);
            myLogger.info("Chunk Size for Tags: " + CHUNK_SIZE);
            int numOfChunks = tagsToChunks(tagCount);
            myLogger.info("Number of Chunks: " + numOfChunks);
            int numTagsPadded = numOfChunks * CHUNK_SIZE;
            for (int mi = 0; mi < maxMapping; mi++) {
 //               h5.compounds().createArray("map" + mi, tmiType, numTagsPadded, CHUNK_SIZE);
                h5.compounds().createArray("map" + mi, tmiType, numTagsPadded, CHUNK_SIZE,genoFeatures);
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
//                        h5.compounds().writeArrayBlock("map" + mi, tmiType, thTMI, block);
                        h5.compounds().writeArrayBlock("map" + mi, tmiType, thTMI, block);
                        thTMI = new TagMappingInfo[CHUNK_SIZE];
                        block++;
                        myLogger.info("Tag locations written: " + (i + 1));
                    }
                }
                myLogger.info("...map" + mi + " positions written");
            }
            
            h5.createByteMatrix("variantDef", NUM_UNITS_TO_CACHE_ON_GET, maxVariants);
            int numVariants = inTags.getMaxNumVariants();
            if (numVariants > maxVariants) {
                throw new IllegalArgumentException("TagsOnPhysMapHDF5: createFile: max variants can't be less than original TOPM Variant Defs: " + numVariants);
            }
            byte[][] temp = new byte[tagCount][maxVariants];
            for (int i = 0; i < tagCount; i++) {
                for (int j = 0; j < maxVariants; j++) {
                    if (j < numVariants) {
                        temp[i][j] = inTags.getVariantDef(i, j);
                    } else {
                        temp[i][j] = BYTE_MISSING;
                    }
                }
            }
            h5.writeByteMatrix("variantDef", temp);

            h5.createByteMatrix("variantPosOff", NUM_UNITS_TO_CACHE_ON_GET, maxVariants);
            temp = new byte[tagCount][maxVariants];
            for (int i = 0; i < tagCount; i++) {
                for (int j = 0; j < maxVariants; j++) {
                    if (j < numVariants) {
                        temp[i][j] = inTags.getVariantPosOff(i, j);
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
        this(filename, true);
    }

    public TagsOnPhysMapHDF5(String theHDF5file, boolean cacheAllMappingBlocks) {
        this.cacheAllMappingBlocks = cacheAllMappingBlocks;
        System.out.println("Opening :" + theHDF5file);
        myHDF5 = HDF5Factory.open(theHDF5file);
        hasDetailedMapping=myHDF5.exists("map0");
        myNumTags = myHDF5.getIntAttribute("/", "tagCount");
        this.tagLengthInLong = myHDF5.getIntAttribute("/", "tagLengthInLong");
        this.tags = myHDF5.readLongMatrix("tags");
        this.tagLength = myHDF5.readByteArray("tagLength");
        this.multimaps = myHDF5.readByteArray("multimaps");
        this.maxMapping = myHDF5.getIntAttribute("/", "maxMapping");
        tmiType = myHDF5.compounds().getInferredType(TagMappingInfo.class);
        cachedTMIBlock = new TagMappingInfo[maxMapping][];
        if(hasDetailedMapping) cacheMappingInfo(0);
        if(!myHDF5.exists("bestStrand")) {populateBestMappings();}
        else {
            bestStrand=myHDF5.readByteArray("bestStrand");
            bestChr= myHDF5.readIntArray("bestChr");
            bestStartPos= myHDF5.readIntArray("bestStartPos");
        }
        loadVariantsIntoMemory();
        System.out.println(theHDF5file + " read with tags:" + myNumTags);
        populateChrAndPositions();
        System.gc();
        Random r=new Random();
//        for (int i = 0; i < 10; i++) {
//            int tag=r.nextInt(myNumTags);
//            System.out.println(Arrays.toString(getVariantDefArray(tag)));
//        }
    }
    
    private boolean populateBestMappings() {
        bestStrand=new byte[myNumTags];
        bestChr=new int[myNumTags];
        bestStartPos=new int[myNumTags];
        if(myHDF5.exists("map0")==false) {
            Arrays.fill(bestStrand, TOPMInterface.BYTE_MISSING);
            Arrays.fill(bestChr, TOPMInterface.INT_MISSING);
            Arrays.fill(bestStartPos, TOPMInterface.INT_MISSING);
            return false;
        }
        for (int i = 0; i < myNumTags; i++) {
            bestStrand[i]=getStrand(i);
            bestChr[i]=getChromosome(i);
            bestStartPos[i]=getStartPosition(i);
        }
        return true;
    }
    
    private boolean loadVariantsIntoMemory() {
        int howManyDef=0;
        int readBlock=4096*16;
        variantDefs=new byte[myNumTags][];
        variantOffsets=new byte[myNumTags][];
        if(!myHDF5.exists("variantDef")) return false;
        for (int blockStep = 0; blockStep < myNumTags; blockStep+=readBlock) {
            int blockSize=(myNumTags-blockStep<readBlock)?myNumTags-blockStep:readBlock;
            byte[][] vd=myHDF5.readByteMatrixBlockWithOffset("variantDef",blockSize,myMaxVariants,blockStep,0);
//            System.out.println(Arrays.toString(vd[0]));
            byte[][] vo=myHDF5.readByteMatrixBlockWithOffset("variantPosOff",blockSize,myMaxVariants,blockStep,0);
            for (int j = 0; j < blockSize; j++) {
                int cnt=0;
                for (byte bs : vd[j]) {if (bs!=TOPMInterface.BYTE_MISSING) cnt++;}
                byte[] vdReDim=new byte[cnt];
                byte[] voReDim=new byte[cnt];
                for (int i = 0; i < cnt; i++) {
                    vdReDim[i]=vd[j][i];
                    voReDim[i]=vo[j][i];
                    howManyDef++;
                }
                variantDefs[blockStep+j]=vdReDim;
                variantOffsets[blockStep+j]=voReDim;
            }
            
            //byte[] vd=myHDF5.readByteArrayBlockWithOffset(null, i, i)
        }
        System.out.println("Real Variant Defs:"+howManyDef);
        return true;
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
//        writeCachedVariantDefs();
//        writeCachedVariantOffsets();
//        saveCacheBackToFile();
    }

    //    private void cacheMappingInfo(int index) {
    //        TagMappingInfo[] aTMI=h5.compounds().readArrayBlockWithOffset("map0", tmiType, 1, index);
    //        this.cachedTMI=aTMI[0];
    //        this.cachedMappingIndex=index;
    //    }
    public TagMappingInfo getAlternateTagMappingInfo(int index, int mapIndex) {
        if(hasDetailedMapping) return null;
        if (cacheAllMappingBlocks == false) {
            cacheMappingInfo(index);
        }
        if (cachedMappingIndex != index) {
            cacheMappingInfo(index);
        }
        return cachedTMIBlock[mapIndex][index % CHUNK_SIZE];
    }

    public void setAlternateTagMappingInfo(int index, int mapIndex, TagMappingInfo theTMI) {
        if(!hasDetailedMapping) throw new IllegalStateException("Detailed mapping not present");
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
        if(!hasDetailedMapping) throw new IllegalStateException("Detailed mapping not present");
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
    public int getChromosome(int index) {
        if(bestChr!=null) return index;
        if (cachedMappingIndex != index) {
            cacheMappingInfo(index);
        }
        return cachedTMI.chromosome;
    }

    @Override
    public int[] getChromosomes() {
        if(myChromosomes==null) populateChrAndPositionsX();
        return myChromosomes;
    }
    
    private void populateChrAndPositions() {
        long chrSum=0;
        System.out.println("chrSum"+chrSum);
        TreeMap<Integer,TreeSet<Integer>> theChrs=new TreeMap<Integer,TreeSet<Integer>>();
        for (int i = 0; i < myNumTags; i++) {
            int chr=getChromosome(i);
            if(chr!=TOPMInterface.INT_MISSING) {
                if(!theChrs.containsKey(chr)) {
                    TreeSet<Integer> thePos=new TreeSet<Integer>();
                    thePos.add(getStartPosition(i));
                    theChrs.put(chr, thePos);
                } else {
                    theChrs.get(chr).add(getStartPosition(i));
                }
            }
        }
        myChromosomes=new int[theChrs.size()];
        int cnt=0;
      //  for (int aChr : theChrs) {myChromosomes[cnt++]=aChr;}
        
    }
    
    private void populateChrAndPositionsX() {
        TreeSet<Integer> theChrs=new TreeSet<Integer>();
        for (int i = 0; i < myNumTags; i++) {
            int chr=getChromosome(i);
            if(chr!=TOPMInterface.INT_MISSING) {
                if(!theChrs.contains(i)) {
                    theChrs.add(chr);
                }
            }
        }
        myChromosomes=new int[theChrs.size()];
        int cnt=0;
        for (int aChr : theChrs) {myChromosomes[cnt++]=aChr;}
        
    }

    @Override
    public byte getDcoP(int index) {
        if(!hasDetailedMapping) throw new IllegalStateException("Detailed mapping not present");
        if (cachedMappingIndex != index) {
            cacheMappingInfo(index);
        }
        return cachedTMI.dcoP;
    }

    @Override
    public byte getDivergence(int index) {
        if(!hasDetailedMapping) throw new IllegalStateException("Detailed mapping not present");
        if (cachedMappingIndex != index) {
            cacheMappingInfo(index);
        }
        return cachedTMI.divergence;
    }

    @Override
    public int getEndPosition(int index) {
        if(!hasDetailedMapping) throw new IllegalStateException("Detailed mapping not present");
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
        if(!hasDetailedMapping) throw new IllegalStateException("Detailed mapping not present");
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
        if(bestStartPos!=null)  return bestStartPos[index];
        if (cachedMappingIndex != index) {
            cacheMappingInfo(index);
        }
        return cachedTMI.startPosition;
    }

    @Override
    public byte getStrand(int tagIndex) {
        if(bestStrand!=null)  return bestStrand[tagIndex];
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
        if(variantDefs[tagIndex]==null) return TOPMInterface.BYTE_MISSING;
        return variantDefs[tagIndex][variantIndex];
    }

    @Override
    public byte[] getVariantDefArray(int tagIndex) {
        return variantDefs[tagIndex];
    }

    @Override
    public byte getVariantPosOff(int tagIndex, int variantIndex) {
        if(variantDefs[tagIndex]==null) return TOPMInterface.BYTE_MISSING;
        return variantOffsets[tagIndex][variantIndex];
    }

    @Override
    public byte[] getVariantPosOffArray(int tagIndex) {
        return variantOffsets[tagIndex];
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
    public synchronized void setVariantDef(int tagIndex, int variantIndex, byte def) {
        myHDF5.writeByteMatrixBlockWithOffset("variantDef", new byte[][]{{def},}, tagIndex, variantIndex);
        variantDefs[tagIndex][variantIndex]=def;
    }

    @Override
    public synchronized void setVariantPosOff(int tagIndex, int variantIndex, byte offset) {
        myHDF5.writeByteMatrixBlockWithOffset("variantPosOff", new byte[][]{{offset},}, tagIndex, variantIndex);
        variantOffsets[tagIndex][variantIndex]=offset;
    }
    
    /**
     * Preferred method for setting variant information
     * @param tagIndex
     * @param defAndOffset Two dimension [0=def, 1=offset][upto 16 bytes for each SNP]
     */
    public synchronized void setAllVariantInfo(int tagIndex, byte[][] defAndOffset) {
        myHDF5.writeByteMatrixBlockWithOffset("variantDef", new byte[][]{defAndOffset[0]}, tagIndex, 0);
        myHDF5.writeByteMatrixBlockWithOffset("variantPosOff", new byte[][]{defAndOffset[1]}, tagIndex, 0);
        variantDefs[tagIndex]=defAndOffset[1];
        variantOffsets[tagIndex]=defAndOffset[1];
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
