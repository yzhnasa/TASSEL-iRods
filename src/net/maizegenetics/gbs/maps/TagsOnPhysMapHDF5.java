/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.gbs.maps;

import cern.colt.GenericSorting;
import cern.colt.Swapper;
import cern.colt.function.IntComparator;
import ch.systemsx.cisd.hdf5.*;
import java.io.File;
import java.util.Arrays;
import net.maizegenetics.gbs.tagdist.AbstractTags;
import net.maizegenetics.pal.alignment.Locus;
import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.SimpleIdGroup;

/**
 *
 * @author edbuckler
 */
public class TagsOnPhysMapHDF5 extends AbstractTags implements TOPMInterface {

    static int chunkSize = 1;// << 16;
    public final static byte byteMissing = Byte.MIN_VALUE;
    public final static int intMissing = Integer.MIN_VALUE;
    public int maxVariants = 8;
    byte[] multimaps;  // number of locations this tagSet maps to; unknown = Byte.MIN_VALUE; multiple, but unknown number = 99
    private IHDF5Reader h5=null;
    private int cachedMappingIndex=-1;
    private TagMappingInfo cachedTMI=null;
    //    int[] chromosome;  // 4 bytes
//    byte[] strand; // 1 = same sense as reference FASTA file.  -1 = opposite sense.  unknown = Byte.MIN_VALUE  // 1 byte
//    int[] startPosition;  // chromosomal position of the barcoded end of the tag  // 4 bytes
//    int[] endPosition;  // chromosomal position of the common adapter end of the tag (smaller than startPosition if tag matches minus strand)  // 4 bytes
//    byte[] divergence;  // number of diverging bp (edit distance) from reference, unknown = Byte.MIN_VALUE
//    byte[][] variantPosOff;  // offset from position minimum, maximum number of variants is defined above  // maxVariants bytes
//    byte[][] variantDef; // allele state - A, C, G, T or some indel definition  // maxVariants bytes
//    byte[] dcoP, mapP;  //Round(Log2(P)), unknown Byte.MIN_VALUE
    //if these disagree with the location, then set the p to negative
    // 1+4+1+4+4+1+8+8+1+1 = 33 bytes per position + 16 bytes for a two long tag + 1 byte for tagLength in bases = 50 bytes
    // ~50 bytes per position.  If we have 10 million tags then this will be at 500M byte data structure.
    HDF5CompoundType<TagMappingInfo> tmiType=null;
    int[] indicesOfSortByPosition;
    public int tagNum = 0;
    boolean redundantTags = true;  // this field has not been utilized yet

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
        IHDF5WriterConfigurator config = HDF5Factory.configure(new File(newHDF5file));
        System.out.println("Creating HDF5 file: " + newHDF5file);
        config.overwrite();
//        config.dontUseExtendableDataTypes();
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
        for (int i = 0; i < maxMapping; i++) {
            h5.compounds().createArray("map"+i, tmiType, tagCount);  
        }
        TagMappingInfo[] thTMI=new TagMappingInfo[chunkSize];
        for(int i=0; i<inTags.getTagCount(); i++) {
            thTMI[i%chunkSize] = (new TagMappingInfo(inTags.getChromosome(i),
                    inTags.getStrand(i), inTags.getStartPosition(i), 
                    inTags.getEndPosition(i), inTags.getDivergence(i)));
            if((i+1)%chunkSize==0) {
                h5.compounds().writeArrayBlockWithOffset("map0", tmiType, thTMI, i);
                thTMI=new TagMappingInfo[chunkSize];
            }
        }
        System.out.println("...map0 position written");
        h5.createByteMatrix("variantDef", maxVariants, 1);
        h5.createByteMatrix("variantPosOff", maxVariants, 1);
        h5.writeByteMatrix("variantDef", inTags.variantDef);
        h5.writeByteMatrix("variantPosOff", inTags.variantPosOff);
        System.out.println("...variant Positions written");
        h5.close();
    }
    
    public TagsOnPhysMapHDF5(String theHDF5file) {
        System.out.println("Opening :"+theHDF5file);
        h5=HDF5Factory.openForReading(theHDF5file);
        this.tagNum=h5.getIntAttribute("/", "tagCount");
        this.tagLengthInLong=h5.getIntAttribute("/", "tagLengthInLong");
        this.tags=h5.readLongMatrix("tags");
        this.tagLength=h5.readByteArray("tagLength");
        this.multimaps=h5.readByteArray("multimaps");
        tmiType = h5.compounds().getInferredType(TagMappingInfo.class);
        cacheMappingInfo(0);
        System.out.println(theHDF5file+" read with tags:"+tagNum);
    }
    
    private void cacheMappingInfo(int index) {
        TagMappingInfo[] aTMI=h5.compounds().readArrayBlockWithOffset("map0", tmiType, 1, index);
        this.cachedTMI=aTMI[0];
        this.cachedMappingIndex=index;
    }

    @Override
    public int addVariant(int tagIndex, byte offset, byte base) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int getChromosome(int index) {
        if(cachedMappingIndex!=index) cacheMappingInfo(index);
        return cachedTMI.chromosome;
    }

    @Override
    public int[] getChromosomes() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public byte getDcoP(int index) {
        if(cachedMappingIndex!=index) cacheMappingInfo(index);
        return cachedTMI.dcoP;
    }

    @Override
    public byte getDivergence(int index) {
        if(cachedMappingIndex!=index) cacheMappingInfo(index);
        return cachedTMI.divergence;
    }

    @Override
    public int getEndPosition(int index) {
        if(cachedMappingIndex!=index) cacheMappingInfo(index);
        return cachedTMI.divergence;
    }

    @Override
    public Locus[] getLoci() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public Locus getLocus(int tagIndex) {
        if(cachedMappingIndex!=tagIndex) cacheMappingInfo(tagIndex);
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public byte getMapP(int index) {
        if(cachedMappingIndex!=index) cacheMappingInfo(index);
        return cachedTMI.mapP;
    }

    @Override
    public byte getMultiMaps(int index) {
        return this.multimaps[index];
    }

    @Override
    public int[] getPositionArray(int index) {
        if(cachedMappingIndex!=index) cacheMappingInfo(index);
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
        if(cachedMappingIndex!=index) cacheMappingInfo(index);
        return cachedTMI.startPosition;
    }

    @Override
    public byte getStrand(int tagIndex) {
        if(cachedMappingIndex!=tagIndex) cacheMappingInfo(tagIndex);
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
 /*   
     void initPhysicalSort() {
        System.out.println("initPhysicalSort");
        indicesOfSortByPosition = new int[tagNum];
        for (int i = 0; i < indicesOfSortByPosition.length; i++) {
            indicesOfSortByPosition[i] = i;
        }
        Swapper swapperPos = new Swapper() {
            public void swap(int a, int b) {
                int t1;
                t1 = indicesOfSortByPosition[a];
                indicesOfSortByPosition[a] = indicesOfSortByPosition[b];
                indicesOfSortByPosition[b] = t1;
            }
        };
        IntComparator compPos = new IntComparator() {

            public int compare(int a, int b) {
                int index1 = indicesOfSortByPosition[a];
                int index2 = indicesOfSortByPosition[b];
                if (chromosome[index1] < chromosome[index2]) {
                    return -1;
                }
                if (chromosome[index1] > chromosome[index2]) {
                    return 1;
                }
                if (startPosition[index1] < startPosition[index2]) {
                    return -1;
                }
                if (startPosition[index1] > startPosition[index2]) {
                    return 1;
                }
                if (strand[index1] < strand[index2]) {
                    return -1;
                }
                if (strand[index1] > strand[index2]) {
                    return 1;
                }
                for (int i = 0; i < tagLengthInLong; i++) {
                    if (tags[i][index1] < tags[i][index2]) {
                        return -1;
                    }
                    if (tags[i][index1] > tags[i][index2]) {
                        return 1;
                    }
                }
                return 0;
            }
        };
        System.out.println("Position index sort begin.");
        GenericSorting.quickSort(0, indicesOfSortByPosition.length, compPos, swapperPos);
        System.out.println("Position index sort end.");
    }
*/
    public static void main(String[] args) {
        String inTOPMFile = "/Users/edbuckler/SolexaAnal/GBS/build20120701/test/topm100.txt";
        String theTOPMH5 = "/Users/edbuckler/SolexaAnal/GBS/build20120701/test/topm100.h5";
//        String inTOPMFile = "/Users/edbuckler/SolexaAnal/GBS/build20120701/04_TOPM/Bowtie2/AllZeaMasterTags_c10_20120703.topm";
//        String theTOPMH5 = "/Users/edbuckler/SolexaAnal/GBS/build20120701/test/AllZeaMasterTags_c10_20120703.topm.h5";
        TagsOnPhysicalMap inTOPM = new TagsOnPhysicalMap(inTOPMFile, false);
        TagsOnPhysMapHDF5.createFile(inTOPM, theTOPMH5,4,8);
        TagsOnPhysMapHDF5 h5TOPM=new TagsOnPhysMapHDF5(theTOPMH5);
        int count=0, error=0;
        for (int i = 0; i < 99000; i++) {
            if(inTOPM.getChromosome(i)>-1)
                count++;
            if(inTOPM.getChromosome(i)!=h5TOPM.getChromosome(i)) error++;
            System.out.printf("%d %d %d %n",i,inTOPM.getChromosome(i),h5TOPM.getChromosome(i));
            
        }
        System.out.println(count+" error="+error);
    }
}
