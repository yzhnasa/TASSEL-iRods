/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.gbs.pipeline;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.TreeSet;
import net.maizegenetics.gbs.maps.TagMappingInfoV3;
import net.maizegenetics.gbs.maps.TagsOnPhysicalMapV3;
import net.maizegenetics.gbs.util.BaseEncoder;
import net.maizegenetics.gbs.util.SAMUtils;
import net.maizegenetics.util.MultiMemberGZIPInputStream;

/**
 * Methods to annotate TOPM file, including adding mapping info from aligners, adding PE tag position and genetic position, model prediction for the best position
 * Code of mappingSource. 0: Bowtie2; 1: BWA; 2: BLAST; 3: PE; 4: Genetic Mapping
 * @author Fei Lu
 */
public class AnnotateTOPM {
    /**TOPM file that will be annotated*/
    TagsOnPhysicalMapV3 topm;
    /**Record Sam record. When tmiBuffers[0] is full, output to TOPM block. Substitute tmiBuffers[i] with tmiBuffers[i+1]. Size = numBuffers×maxMappingNum(-K option)×TOPM CHUNK_SIZE*/
    TagMappingInfoV3[][][] tmiBuffers = null;
    /**check if tmiBuffer[0] is full for output*/
    boolean[][] bufferLights = null;
    /**Number of buffers*/
    int bufferNum = 2;
    /**record how mapping tags are imported in each buffer*/
    int[] lightCounts;
    /**Actual tag index of each tmiBuffers[i]*/
    int[] bufferStartTagIndex = null;
    /**Tag index range of the tmiBuffers*/
    int[] bufferTagIndexRange = null;
    /**When the second buffer has filled at least with this cutoff, update buffer*/
    int updateBufferCountCutoff;
    
    public AnnotateTOPM (TagsOnPhysicalMapV3 topm) {
        this.topm = topm;   
    }
    
    public void annotateWithBowtie2 (String samFileS, int maxMappingNum) {
        String[] dataSetNames = topm.creatTagMappingInfoDatasets(topm.getMaxMappingCount(), maxMappingNum);
        this.iniTMIBuffers(bufferNum, maxMappingNum);
        System.out.println("Reading SAM format tag alignment (Bowtie2) from: " + samFileS);
        System.out.println("Coverting SAM to TOPMHDF5...");
        byte mappingSource = 0;
        try {
            BufferedReader br;
            if (samFileS.endsWith(".gz")) {
                br = new BufferedReader(new InputStreamReader(new MultiMemberGZIPInputStream(new FileInputStream(new File(samFileS)))));
            } else {
                br = new BufferedReader(new FileReader(new File(samFileS)), 65536);
            }
            while(br.readLine().startsWith("@PG")==false) {};
            String inputStr = null;
            while((inputStr = br.readLine())!=null) {
                String[] temp =inputStr.split("\\s");
                int orientiation=Integer.parseInt(temp[1]);
                int chr = Integer.MIN_VALUE;
                byte strand = Byte.MIN_VALUE;
                int startPos = Integer.MIN_VALUE;
                int endPos = Integer.MIN_VALUE;
                short mappingScore = Short.MIN_VALUE;
                byte divergence = Byte.MIN_VALUE;
                String seqS = temp[9];
                if (orientiation == 4) {
                    
                }
                else if (orientiation == 16 || orientiation == 272) {
                    seqS = BaseEncoder.getReverseComplement(seqS);
                    chr = Integer.parseInt(temp[2]);
                    strand = -1;
                    int[] alignSpan = SAMUtils.adjustCoordinates(temp[5], Integer.parseInt(temp[3]));
                    startPos = alignSpan[1];
                    endPos = alignSpan[0];
                    mappingScore = Short.parseShort(temp[11].split(":")[2]);
                    if (temp[17].startsWith("NM")) {
                        divergence = Byte.parseByte(temp[17].split(":")[2]);
                    }
                    else {
                        divergence = Byte.parseByte(temp[16].split(":")[2]);
                    }
                }
                else {
                    chr = Integer.parseInt(temp[2]);
                    strand = 1;
                    int[] alignSpan = SAMUtils.adjustCoordinates(temp[5], Integer.parseInt(temp[3]));
                    startPos = alignSpan[0];
                    endPos = alignSpan[1];
                    mappingScore = Short.parseShort(temp[11].split(":")[2]);
                    if (temp[17].startsWith("NM")) {
                        divergence = Byte.parseByte(temp[17].split(":")[2]);
                    }
                    else {
                        divergence = Byte.parseByte(temp[16].split(":")[2]);
                    }
                }
                TagMappingInfoV3 theTMI = new TagMappingInfoV3(chr, strand, startPos, endPos, divergence, mappingSource, mappingScore);
                long[] seq = BaseEncoder.getLongArrayFromSeq(seqS,topm.getTagSizeInLong()*32);
                int tagIndex = topm.getTagIndex(seq);
                if (tagIndex < this.bufferTagIndexRange[0] || tagIndex >= this.bufferTagIndexRange[bufferNum-1]) {
                    System.out.println("The index of the tag from sam file is out of buffer range. Program quits.");
                    System.out.println("Please increase the buffer number");
                    System.exit(1);
                }
                int bufferIndex = Arrays.binarySearch(bufferStartTagIndex, tagIndex);
                if (bufferIndex < 0) bufferIndex = -bufferIndex-2;
                int bufferTagIndex = tagIndex % topm.getChunkSize();
                int mappingBlockIndex = this.getMappingBlockIndex(bufferIndex, bufferTagIndex);
                if (mappingBlockIndex == Integer.MIN_VALUE) continue;
                tmiBuffers[bufferIndex][mappingBlockIndex][bufferTagIndex] = theTMI;
                if (bufferLights[bufferIndex][bufferTagIndex] == false) {
                    lightCounts[bufferIndex]++;
                }
                bufferLights[bufferIndex][bufferTagIndex] = true;
                if (lightCounts[0] == topm.getChunkSize() && lightCounts[1] > this.updateBufferCountCutoff) {
                    this.saveTMIBufferToTOPM(tmiBuffers[0], dataSetNames, this.bufferStartTagIndex[0]/topm.getChunkSize(), mappingSource);
                    this.updateTMIBuffer();
                }   
            }
            this.saveTMIBufferToTOPM(tmiBuffers[0], dataSetNames, this.bufferStartTagIndex[0]/topm.getChunkSize(), mappingSource);
            topm.setMaxMapping(topm.getMaxMappingCount()+maxMappingNum);
            br.close();
        } catch (Exception e) {
            
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    public void annotateWithBWA (String samFileS, int maxMappingNum) {
        String[] dataSetNames = topm.creatTagMappingInfoDatasets(topm.getMaxMappingCount(), maxMappingNum);
        this.iniTMIBuffers(bufferNum, maxMappingNum);
        System.out.println("Reading SAM format tag alignment (BWA) from: " + samFileS);
        System.out.println("Coverting SAM to TOPMHDF5...");
        byte mappingSource = 1;
        try {
            BufferedReader br;
            if (samFileS.endsWith(".gz")) {
                br = new BufferedReader(new InputStreamReader(new MultiMemberGZIPInputStream(new FileInputStream(new File(samFileS)))));
            } else {
                br = new BufferedReader(new FileReader(new File(samFileS)), 65536);
            }
            String inputStr = null;
            while ((inputStr=br.readLine()).startsWith("@")) {};
            while(inputStr!=null) {
                String[] temp =inputStr.split("\\s");
                int orientiation=Integer.parseInt(temp[1]);
                int chr = Integer.MIN_VALUE;
                byte strand = Byte.MIN_VALUE;
                int startPos = Integer.MIN_VALUE;
                int endPos = Integer.MIN_VALUE;
                short mappingScore = Short.MIN_VALUE;
                byte divergence = Byte.MIN_VALUE;
                String seqS = temp[9];
                String XAString = null;
                if (temp.length == 20) {
                    XAString = temp[19].replaceFirst("XA:Z:", "");
                }
                if (orientiation == 4) {
                    
                }
                else if (orientiation == 16) {
                    seqS = BaseEncoder.getReverseComplement(seqS);
                    chr = Integer.parseInt(temp[2]);
                    strand = -1;
                    int[] alignSpan = SAMUtils.adjustCoordinates(temp[5], Integer.parseInt(temp[3]));
                    startPos = alignSpan[1];
                    endPos = alignSpan[0];
                    divergence = Byte.parseByte(temp[12].split(":")[2]);
                }
                else {
                    chr = Integer.parseInt(temp[2]);
                    strand = 1;
                    int[] alignSpan = SAMUtils.adjustCoordinates(temp[5], Integer.parseInt(temp[3]));
                    startPos = alignSpan[0];
                    endPos = alignSpan[1];
                    divergence = Byte.parseByte(temp[12].split(":")[2]);
                }
                TagMappingInfoV3 theTMI = new TagMappingInfoV3(chr, strand, startPos, endPos, divergence, mappingSource, mappingScore);
                long[] seq = BaseEncoder.getLongArrayFromSeq(seqS,topm.getTagSizeInLong()*32);
                int tagIndex = topm.getTagIndex(seq);
                if (tagIndex < this.bufferTagIndexRange[0] || tagIndex >= this.bufferTagIndexRange[bufferNum-1]) {
                    System.out.println("The index of the tag from sam file is out of buffer range. Program quits.");
                    System.out.println("Please increase the buffer number");
                    System.exit(1);
                }
                int bufferIndex = Arrays.binarySearch(bufferStartTagIndex, tagIndex);
                if (bufferIndex < 0) bufferIndex = -bufferIndex-2;
                int bufferTagIndex = tagIndex % topm.getChunkSize();
                int mappingBlockIndex = this.getMappingBlockIndex(bufferIndex, bufferTagIndex);
                if (mappingBlockIndex == Integer.MIN_VALUE) continue;
                tmiBuffers[bufferIndex][mappingBlockIndex][bufferTagIndex] = theTMI;
                if (bufferLights[bufferIndex][bufferTagIndex] == false) {
                    lightCounts[bufferIndex]++;
                }
                if (XAString != null) {
                    temp = XAString.split(";");
                    for (int i = 0; i < temp.length; i++) {
                        mappingBlockIndex = this.getMappingBlockIndex(bufferIndex, bufferTagIndex);
                        if (mappingBlockIndex == Integer.MIN_VALUE) break;
                        String[] tem = temp[i].split(",");
                        chr = Integer.parseInt(tem[0]);
                        if (tem[1].startsWith("+")) {
                            strand = 1;
                            int[] alignSpan = SAMUtils.adjustCoordinates(tem[2], Integer.parseInt(tem[1].substring(1)));
                            startPos = alignSpan[1];
                            endPos = alignSpan[0];
                        }
                        else {
                            strand = -1;
                            int[] alignSpan = SAMUtils.adjustCoordinates(tem[2], Integer.parseInt(tem[1].substring(1)));
                            startPos = alignSpan[0];
                            endPos = alignSpan[1];
                        }
                        divergence = Byte.parseByte(tem[3]);
                        theTMI = new TagMappingInfoV3(chr, strand, startPos, endPos, divergence, mappingSource, mappingScore);
                        tmiBuffers[bufferIndex][mappingBlockIndex][bufferTagIndex] = theTMI;
                    }
                }
                bufferLights[bufferIndex][bufferTagIndex] = true;
                if (lightCounts[0] == topm.getChunkSize() && lightCounts[1] > this.updateBufferCountCutoff) {
                    this.saveBWATMIBufferToTOPM(tmiBuffers[0], dataSetNames, this.bufferStartTagIndex[0]/topm.getChunkSize(), mappingSource);
                    this.updateTMIBuffer();
                }
                inputStr=br.readLine();
            }
            this.saveBWATMIBufferToTOPM(tmiBuffers[0], dataSetNames, this.bufferStartTagIndex[0]/topm.getChunkSize(), mappingSource);
            topm.setMaxMapping(topm.getMaxMappingCount()+maxMappingNum);
            br.close();
        } catch (Exception e) {
            
            e.printStackTrace();
            System.exit(1);
        }       
    }
            
    public void annotateWithBLAST (String blastM8FileS, int maxMappingNum) {
        String[] dataSetNames = topm.creatTagMappingInfoDatasets(topm.getMaxMappingCount(), maxMappingNum);
        this.iniTMIBuffers(2, maxMappingNum);
        System.out.println("Reading SAM format tag alignment (BLAST) from: " + blastM8FileS);
        System.out.println("Coverting SAM to TOPMHDF5...");
        byte mappingSource = 2;
        try {
            BufferedReader br;
            if (blastM8FileS.endsWith(".gz")) {
                br = new BufferedReader(new InputStreamReader(new MultiMemberGZIPInputStream(new FileInputStream(new File(blastM8FileS)))));
            } else {
                br = new BufferedReader(new FileReader(new File(blastM8FileS)), 65536);
            }
            String inputStr = null;
            while((inputStr = br.readLine())!=null) {
                String[] temp =inputStr.split("\\s+");
                int chr = Integer.parseInt(temp[1]);
                byte strand = Byte.MIN_VALUE;
                int startPos = Integer.parseInt(temp[8]);
                int endPos = Integer.parseInt(temp[9]);
                if (startPos < endPos) {
                    strand = 1;
                }
                else {
                    strand = -1;
                }
                short mappingScore = Short.parseShort(temp[11].replaceAll("\\..+", ""));
                byte divergence = Byte.MIN_VALUE;
                TagMappingInfoV3 theTMI = new TagMappingInfoV3(chr, strand, startPos, endPos, divergence, mappingSource, mappingScore);
                int tagIndex = Integer.parseInt(temp[0]);
                if (tagIndex < this.bufferTagIndexRange[0] || tagIndex >= this.bufferTagIndexRange[bufferNum-1]) {
                    System.out.println("The index of the tag from sam file is out of buffer range. Program quits.");
                    System.out.println("Please increase the buffer number");
                    System.exit(1);
                }
                int bufferIndex = Arrays.binarySearch(bufferStartTagIndex, tagIndex);
                if (bufferIndex < 0) bufferIndex = -bufferIndex-2;
                int bufferTagIndex = tagIndex % topm.getChunkSize();
                int mappingBlockIndex = this.getMappingBlockIndex(bufferIndex, bufferTagIndex);
                if (mappingBlockIndex == Integer.MIN_VALUE) continue;
                tmiBuffers[bufferIndex][mappingBlockIndex][bufferTagIndex] = theTMI;
                if (bufferLights[bufferIndex][bufferTagIndex] == false) {
                    lightCounts[bufferIndex]++;
                }
                bufferLights[bufferIndex][bufferTagIndex] = true;
                if (lightCounts[0] == topm.getChunkSize() && lightCounts[1] > this.updateBufferCountCutoff) {
                    this.saveTMIBufferToTOPM(tmiBuffers[0], dataSetNames, this.bufferStartTagIndex[0]/topm.getChunkSize(), mappingSource);
                    this.updateTMIBuffer();
                }   
            }
            this.saveTMIBufferToTOPM(tmiBuffers[0], dataSetNames, this.bufferStartTagIndex[0]/topm.getChunkSize(), mappingSource);
            topm.setMaxMapping(topm.getMaxMappingCount()+maxMappingNum);
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    public void annotateWithPE (String PETOPMFileS) {
        
    }
    
    public void annotateWithGM (String TOGMFileS) {
        
    }
    
    /**
     * BWA doesn't provide mapping score, so the rank is just based on the order provided for multiple hits
     * @param tmiBuffer
     * @param dataSetNames
     * @param chunkIndex
     * @param mappingSource 
     */
    private void saveBWATMIBufferToTOPM (TagMappingInfoV3[][] tmiBuffer, String[] dataSetNames, int chunkIndex, byte mappingSource) {
        for (int i = 0; i < tmiBuffer.length; i++) {
            for (int j = 0; j < tmiBuffer[i].length; j++) {
                if (tmiBuffer[i][j] != null) continue;
                tmiBuffer[i][j] = new TagMappingInfoV3();
                tmiBuffer[i][j].setMappingSource(mappingSource);
            }
        }
        for (int i = 0; i < tmiBuffer.length; i++) {
            for (int j = 0; j < tmiBuffer[i].length; j++) {
                tmiBuffer[i][j].setMappingRank((byte)i);
            }
        }
        topm.writeTagMappingInfoDataSets(dataSetNames, tmiBuffer, chunkIndex);
    }
    
    private void saveTMIBufferToTOPM (TagMappingInfoV3[][] tmiBuffer, String[] dataSetNames, int chunkIndex, byte mappingSource) {
        for (int i = 0; i < tmiBuffer.length; i++) {
            for (int j = 0; j < tmiBuffer[i].length; j++) {
                if (tmiBuffer[i][j] != null) continue;
                tmiBuffer[i][j] = new TagMappingInfoV3();
                tmiBuffer[i][j].setMappingSource(mappingSource);
            }
        }
        for (int i = 0; i < tmiBuffer[0].length; i++) {
            TreeSet<Short> set = new TreeSet();
            for (int j = 0; j < tmiBuffer.length; j++) {
                set.add(tmiBuffer[j][i].mappingScore);
            }
            Short[] sA = set.toArray(new Short[set.size()]);
            byte[] rank = new byte[tmiBuffer.length];
            for (int j = 0; j < rank.length; j++) {
                rank[j] = (byte)(sA.length - Arrays.binarySearch(sA, tmiBuffer[j][i].mappingScore) -1);
                tmiBuffer[j][i].setMappingRank(rank[j]);
            }
        }
        topm.writeTagMappingInfoDataSets(dataSetNames, tmiBuffer, chunkIndex);
    }
    
    private void updateTMIBuffer () {
        for (int i = 0; i < tmiBuffers.length-1; i++) {
            tmiBuffers[i] = tmiBuffers[i+1];
            bufferStartTagIndex[i] = bufferStartTagIndex[i+1];
            bufferLights[i] = bufferLights[i+1];
            lightCounts[i] = lightCounts[i+1];
        }
        tmiBuffers[tmiBuffers.length-1] = new TagMappingInfoV3[tmiBuffers[0].length][topm.getChunkSize()];
        bufferStartTagIndex[tmiBuffers.length-1]+=topm.getChunkSize();
        bufferLights[tmiBuffers.length-1] = new boolean[topm.getChunkSize()];
        lightCounts[tmiBuffers.length-1] = 0;
        this.calBufferTagIndexRange();
    }
    
    private void iniTMIBuffers (int bufferNum, int maxMappingNum) {
        tmiBuffers = new TagMappingInfoV3[bufferNum][maxMappingNum][topm.getChunkSize()];
        bufferStartTagIndex = new int[bufferNum];
        for (int i = 0; i < bufferNum; i++) {
            bufferStartTagIndex[i] = i*topm.getChunkSize();
        }
        this.calBufferTagIndexRange();
        this.bufferLights = new boolean[bufferNum][topm.getChunkSize()];
        this.lightCounts = new int[bufferNum];
        updateBufferCountCutoff = (int)((double)topm.getChunkSize() * 0.2);
    }
    
    private void calBufferTagIndexRange () {
        this.bufferTagIndexRange = new int[2];
        bufferTagIndexRange[0] = this.bufferStartTagIndex[0];
        bufferTagIndexRange[1] = this.bufferStartTagIndex[bufferStartTagIndex.length-1]+topm.getChunkSize();
    }
    
    private int getMappingBlockIndex (int bufferIndex, int bufferTagIndex) {
        for (int i = 0; i < tmiBuffers[0].length; i++) {
            if (tmiBuffers[bufferIndex][i][bufferTagIndex] == null) {
                return i;
            }
        }
        return Integer.MIN_VALUE;
    }
}
