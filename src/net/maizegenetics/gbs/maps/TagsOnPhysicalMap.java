package net.maizegenetics.gbs.maps;

import cern.colt.GenericSorting;
import cern.colt.Swapper;
import cern.colt.function.IntComparator;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.InputStreamReader;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

import net.maizegenetics.gbs.tagdist.AbstractTags;
import net.maizegenetics.gbs.tagdist.TagCountMutable;
import net.maizegenetics.gbs.tagdist.Tags;

import net.maizegenetics.gbs.util.SAMUtils;
import net.maizegenetics.gbs.util.BaseEncoder;

import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.ImportUtils;
import net.maizegenetics.pal.alignment.Locus;
import net.maizegenetics.pal.alignment.NucleotideAlignmentConstants;
import net.maizegenetics.util.MultiMemberGZIPInputStream;

/**
 * Holds tag data compressed in long and a physical position. This can have two
 * variations either include redundant positions or only unique positions. If
 * redundant than tags that map to multiple regions should be placed adjacently
 *
 * Default 40 bytes per position. If we have 10,000,000 million positions then
 * this will be at 400M byte data structure.
 *
 * User: ed
 */
public class TagsOnPhysicalMap extends AbstractTags implements TOPMInterface {

    public int maxVariants = 8;
    byte[] multimaps;  // number of locations this tagSet maps to; unknown = Byte.MIN_VALUE; multiple, but unknown number = 99
    int[] chromosome;  // 4 bytes
    byte[] strand; // 1 = same sense as reference FASTA file.  -1 = opposite sense.  unknown = Byte.MIN_VALUE  // 1 byte
    int[] startPosition;  // chromosomal position of the barcoded end of the tag  // 4 bytes
    int[] endPosition;  // chromosomal position of the common adapter end of the tag (smaller than startPosition if tag matches minus strand)  // 4 bytes
    byte[] divergence;  // number of diverging bp (edit distance) from reference, unknown = Byte.MIN_VALUE
    byte[][] variantPosOff;  // offset from position minimum, maximum number of variants is defined above  // maxVariants bytes [tag][variant]
    byte[][] variantDef; // allele state - A, C, G, T or some indel definition  // maxVariants bytes [tag][variant]  
    byte[] dcoP, mapP;  //Round(Log2(P)), unknown Byte.MIN_VALUE
    //if these disagree with the location, then set the p to negative
    // 1+4+1+4+4+1+8+8+1+1 = 33 bytes per position + 16 bytes for a two long tag + 1 byte for tagLength in bases = 50 bytes
    // ~50 bytes per position.  If we have 10 million tags then this will be at 500M byte data structure.
    int[] indicesOfSortByPosition;
    public int tagNum = 0;
    boolean redundantTags = true;  // this field has not been utilized yet

    public static enum SAMFormat {

        BWA, BOWTIE2
    };  // Supported SAM formats (each program defines some custom features)
    private SAMFormat mySAMFormat = SAMFormat.BWA;  // BWA by default

    public TagsOnPhysicalMap() {
    }

    public TagsOnPhysicalMap(String inFile, boolean binary) {
        if (binary) {
            readBinaryFile(new File(inFile));
        } else {
            readTextFile(new File(inFile));
        }
        initPhysicalSort();
    }

    public TagsOnPhysicalMap(int rows) {
        initMatrices(rows);
    }

    public TagsOnPhysicalMap(int rows, int tagLengthInLong, int maxVariants) {
        this.tagLengthInLong = tagLengthInLong;
        this.maxVariants = maxVariants;
        initMatrices(rows);
    }

    public TagsOnPhysicalMap(Tags readList) {
        tagLengthInLong = readList.getTagSizeInLong();
        initMatrices(readList.getTagCount());
        for (int i = 0; i < readList.getTagCount(); i++) {
            for (int j = 0; j < tagLengthInLong; j++) {
                tags[j][i] = readList.getTag(i)[j];
            }
        }
    }

    public TagsOnPhysicalMap(TagsOnPhysicalMap oldTOPM, boolean filterDuplicates) {
        this.tagLengthInLong = oldTOPM.tagLengthInLong;
        this.maxVariants = oldTOPM.maxVariants;
        oldTOPM.sortTable(true);
        int uniqueCnt = 1;
        for (int i = 1; i < oldTOPM.getSize(); i++) {
            if (!Arrays.equals(oldTOPM.getTag(i - 1), oldTOPM.getTag(i))) {
                uniqueCnt++;
            }
        }
        System.out.println("The Physical Map File has UniqueTags:" + uniqueCnt + " TotalLocations:" + oldTOPM.getSize());
        initMatrices(uniqueCnt);
        this.tagNum = uniqueCnt;
        uniqueCnt = 0;
        this.copyTagMapRow(oldTOPM, 0, 0, filterDuplicates);
        for (int i = 1; i < oldTOPM.getTagCount(); i++) {
            //  System.out.println(this.getTag(uniqueCnt)[1]+":"+ oldTOPM.getTag(i)[1]);
            //  System.out.println(":"+ oldTOPM.getTag(i)[1]);
            if (!Arrays.equals(this.getTag(uniqueCnt), oldTOPM.getTag(i))) {
                uniqueCnt++;
                //          copyTagMapRow(oldTOPM, i, uniqueCnt, filterDuplicates);
            } else {
                // System.out.printf("i=%d uniqueCnt=%d %n",i, uniqueCnt);
            }
            copyTagMapRow(oldTOPM, i, uniqueCnt, filterDuplicates);
        }
        initPhysicalSort();
    }

    void initMatrices(int rows) {
        tags = new long[tagLengthInLong][rows];  // 16 bytes
        tagLength = new byte[rows];   // length of tag (number of bases)  // 1 byte
        multimaps = new byte[rows];   // number of locations this tagSet maps to, unknown = Byte.MIN_VALUE
        chromosome = new int[rows];   // 4 bytes
        strand = new byte[rows];      // 1 = same sense as reference FASTA file.  -1 = opposite sense.  unknown = Byte.MIN_VALUE  // 1 byte
        startPosition = new int[rows];  // chromosomal position of the barcoded end of the tag  // 4 bytes
        endPosition = new int[rows];  // chromosomal position of the common adapter end of the tag (smaller than startPosition if tag matches minus strand)  // 4 bytes
        divergence = new byte[rows];  // number of diverging bp from reference, unknown = Byte.MIN_VALUE
        variantPosOff = new byte[rows][maxVariants];  // offset from position minimum, maximum number of variants is defined above  // maxVariants bytes
        variantDef = new byte[rows][maxVariants];     // allele state - A, C, G, T or some indel definition  // maxVariants bytes
        dcoP = new byte[rows];
        mapP = new byte[rows];  // Round(Log2(P)), unknown = Byte.MIN_VALUE;  if these disagree with the location, then set the p to negative
        tagNum = rows;
    }

    public void expandMaxVariants(int newMaxVariants) {
        if (newMaxVariants <= this.maxVariants) {
            System.out.println("TagsOnPhysicalMap.expandMaxVariants(" + newMaxVariants + ") not performed because newMaxVariants (" + newMaxVariants
                    + ") <= current maxVariants (" + this.maxVariants + ")");
            return;
        }
        int oldMaxVariants = this.maxVariants;
        byte[][] newVariantPosOff = new byte[tagNum][newMaxVariants];
        byte[][] newVariantDef = new byte[tagNum][newMaxVariants];
        for (int t = 0; t < tagNum; ++t) {
            for (int v = 0; v < this.maxVariants; ++v) {
                newVariantPosOff[t][v] = this.variantPosOff[t][v];
                newVariantDef[t][v] = this.variantDef[t][v];
            }
            for (int v = this.maxVariants; v < newMaxVariants; ++v) {
                newVariantPosOff[t][v] = Byte.MIN_VALUE;
                newVariantDef[t][v] = Byte.MIN_VALUE;
            }
        }
        this.maxVariants = newMaxVariants;
        this.variantPosOff = newVariantPosOff;
        this.variantDef = newVariantDef;
        System.out.println("TagsOnPhysicalMap maxVariants expanded from " + oldMaxVariants + " to " + newMaxVariants);
    }

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

    /**
     * This helps collapse identical reads from different regions of the genome
     * together.
     *
     * @param sourceTOPM
     * @param sourceRow
     * @param destRow
     * @param merge
     */
    public void copyTagMapRow(TagsOnPhysicalMap sourceTOPM, int sourceRow, int destRow, boolean merge) {
        boolean overwrite = true;
        long[] ctag = sourceTOPM.getTag(sourceRow);
        if (Arrays.equals(ctag, this.getTag(destRow)) && merge) {
            overwrite = false;
        }
        for (int i = 0; i < tagLengthInLong; i++) {
            tags[i][destRow] = ctag[i];
        }
        if (overwrite) {
            tagLength[destRow] = sourceTOPM.tagLength[sourceRow];
            multimaps[destRow] = sourceTOPM.multimaps[sourceRow];
            chromosome[destRow] = sourceTOPM.chromosome[sourceRow];
            strand[destRow] = sourceTOPM.strand[sourceRow];
            startPosition[destRow] = sourceTOPM.startPosition[sourceRow];
            endPosition[destRow] = sourceTOPM.endPosition[sourceRow];
            divergence[destRow] = sourceTOPM.divergence[sourceRow];
            for (int j = 0; j < maxVariants; j++) {
                variantPosOff[destRow][j] = sourceTOPM.variantPosOff[sourceRow][j];
                variantDef[destRow][j] = sourceTOPM.variantPosOff[sourceRow][j];
            }
            dcoP[destRow] = sourceTOPM.dcoP[sourceRow];
            mapP[destRow] = sourceTOPM.mapP[sourceRow];
        } else {
            //tagLength[destRow]=tagLength[sourceRow];
            if ((chromosome[destRow] != sourceTOPM.chromosome[sourceRow])
                    || (strand[destRow] != sourceTOPM.strand[sourceRow])
                    || (startPosition[destRow] != sourceTOPM.startPosition[sourceRow])
                    || (endPosition[destRow] != sourceTOPM.endPosition[sourceRow])) {
                multimaps[destRow] += sourceTOPM.multimaps[sourceRow];
                chromosome[destRow] = strand[destRow] = Byte.MIN_VALUE;
                startPosition[destRow] = endPosition[destRow] = Integer.MIN_VALUE;
            }
            //dcoP[destRow]=dcoP[sourceRow];
            //mapP[destRow]=mapP[sourceRow];
        }
    }

    public String printRow(int row) {
        StringBuilder sb = new StringBuilder();
        sb.append(sb);
        //long
        sb.append(BaseEncoder.getSequenceFromLong(this.getTag(row)) + "\t");
        sb.append(printWithMissing(tagLength[row]) + "\t");
        sb.append(printWithMissing(multimaps[row]) + "\t");
        sb.append(printWithMissing(chromosome[row]) + "\t");
        sb.append(printWithMissing(strand[row]) + "\t");
        sb.append(printWithMissing(startPosition[row]) + "\t");
        sb.append(printWithMissing(endPosition[row]) + "\t");
        sb.append(printWithMissing(divergence[row]) + "\t");
        for (int j = 0; j < maxVariants; j++) {
            sb.append(printWithMissing(variantPosOff[row][j]) + "\t");
            sb.append(printWithMissing(variantDef[row][j]) + "\t");
        }
        sb.append(printWithMissing(dcoP[row]) + "\t");
        sb.append(printWithMissing(mapP[row]) + "\t");
        return sb.toString();
    }

    public String printRow(int row, boolean byPosition) {
        if (byPosition) {
            return printRow(indicesOfSortByPosition[row]);
        }
        return printRow(row);
    }

    public void printRows(int numRows) {
        for (int i = 0; i < numRows; i++) {
            System.out.println(printRow(i));
        }
    }

    public void printRows(int numRows, boolean requirePhysPosition, boolean byPosition) {
        int outCount = 0;
        for (int i = 0; outCount < numRows; i++) {
            int r = (byPosition) ? indicesOfSortByPosition[i] : i;
            if ((requirePhysPosition == true) && (chromosome[r] < 1)) {
                continue;
            }
            System.out.println(printRow(r));
            outCount++;
        }
    }

    public void printRows(int numRows, boolean requirePhysPosition, int printChr) {
        int outCount = 0;
        boolean byPosition = true;
        for (int i = 0; outCount < numRows; i++) {
            int r = (byPosition) ? indicesOfSortByPosition[i] : i;
            if ((requirePhysPosition == true) && (chromosome[r] != printChr)) {
                continue;
            }
            System.out.println(printRow(r));
            outCount++;
        }
    }

    public long sortTable(boolean byHaplotype) {
        System.out.print("Starting Read Table Sort ...");
        if (byHaplotype == false) {
            //TODO change the signature at some time
            System.out.print("ERROR:  Position sorting has been eliminated ...");
            return -1;
        }
        long time = System.currentTimeMillis();
        GenericSorting.quickSort(0, tags[0].length, this, this);
        long totalTime = System.currentTimeMillis() - time;
        System.out.println("Done in " + totalTime + "ms");
        initPhysicalSort();
        return totalTime;
    }

    protected void readBinaryFile(File currentFile) {
        int tagsInput = 0;
        try {
            DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(currentFile), 65536));
            System.out.println("File = " + currentFile);
            tagNum = dis.readInt();
            tagLengthInLong = dis.readInt();
            maxVariants = dis.readInt();
            initMatrices(tagNum);
            for (int row = 0; row < tagNum; row++) {
                for (int j = 0; j < tagLengthInLong; j++) {
                    tags[j][row] = dis.readLong();
                }
                tagLength[row] = dis.readByte();
                multimaps[row] = dis.readByte();
                chromosome[row] = dis.readInt();
                strand[row] = dis.readByte();
                startPosition[row] = dis.readInt();
                endPosition[row] = dis.readInt();
                divergence[row] = dis.readByte();
                for (int j = 0; j < maxVariants; j++) {
                    variantPosOff[row][j] = dis.readByte();
                    variantDef[row][j] = dis.readByte();
                }
                dcoP[row] = dis.readByte();
                mapP[row] = dis.readByte();
                tagsInput++;
                if (row % 1000000 == 0) {
                    System.out.println("TagMapFile Row Read:" + row);
                }
            }
            dis.close();
        } catch (Exception e) {
            System.out.println("Error tagsInput=" + tagsInput + " e=" + e);
        }
        System.out.println("Count of Tags=" + tagsInput);
    }

    public boolean variantsDefined(int tagIndex) {
        for (int i = 0; i < maxVariants; i++) {
            if ((variantPosOff[tagIndex][i] != Byte.MIN_VALUE) && (variantDef[tagIndex][i] != Byte.MIN_VALUE)) {
                return true;
            }
        }
        return false;
    }

    public void writeBinaryWVariantsFile(File outFile) {

        int hapsOutput = 0;
        try {
            int[] numTagsWithDefinedVariantsPerChr = new int[20];
            int numTagsWithDefinedVariants = 0;
            for (int row = 0; row < tagNum; row++) {
                if (variantsDefined(row)) {
                    numTagsWithDefinedVariantsPerChr[getChromosome(row)]++;
                    numTagsWithDefinedVariants++;
                }
            }

            System.out.println("writeBinaryWVariantsFile: number tags with defined variants: " + numTagsWithDefinedVariants);
            for (int i = 0; i < numTagsWithDefinedVariantsPerChr.length; i++) {
                if (numTagsWithDefinedVariantsPerChr[i] != 0) {
                    System.out.println("writeBinaryWVariantsFile: Chromosome: " + i + " Has Number Tags with Variants: " + numTagsWithDefinedVariantsPerChr[i]);
                }
            }

            DataOutputStream fw = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outFile), 4000000));

            fw.writeInt(numTagsWithDefinedVariants);
            fw.writeInt(tagLengthInLong);
            fw.writeInt(maxVariants);
            for (int row = 0; row < tagNum; row++) {
                if (variantsDefined(row)) {
                    for (int j = 0; j < tagLengthInLong; j++) {
                        fw.writeLong(tags[j][row]);
                    }
                    fw.writeByte(tagLength[row]);
                    fw.writeByte(multimaps[row]);
                    fw.writeInt(chromosome[row]);
                    fw.writeByte(strand[row]);
                    fw.writeInt(startPosition[row]);
                    fw.writeInt(endPosition[row]);
                    fw.writeByte(divergence[row]);
                    for (int j = 0; j < maxVariants; j++) {
                        fw.writeByte(variantPosOff[row][j]);
                        fw.writeByte(variantDef[row][j]);
                    }
                    fw.writeByte(dcoP[row]);
                    fw.writeByte(mapP[row]);
                    hapsOutput++;
                }
            }
            fw.flush();
            fw.close();
            System.out.println("writeBinaryWVariantsFile: Tag positions written to: " + outFile.toString());
            System.out.println("writeBinaryWVariantsFile: Number of tags in file: " + hapsOutput);
        } catch (Exception e) {
            e.printStackTrace();
            System.out.println("Catch in writing output file e=" + e);
        }
    }

    public void writeBinaryFile(File outFile, int minResolution, boolean requirePhysPosition,
            boolean requireDCOMap, float minDCOP, boolean binary) {
        int hapsOutput = 0;
        try {
            DataOutputStream fw = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outFile), 4000000));
            if (requirePhysPosition) {
                fw.writeInt(mappedTags()[0]);
            } // the index 0 provides the number of tags with unique positions
            else {
                fw.writeInt(tagNum);
            }
            fw.writeInt(tagLengthInLong);
            fw.writeInt(maxVariants);
            for (int row = 0; row < tagNum; row++) {
                if ((requirePhysPosition == true) && (chromosome[row] == Integer.MIN_VALUE)) {
                    continue;
                }
                for (int j = 0; j < tagLengthInLong; j++) {
                    fw.writeLong(tags[j][row]);
                }
                fw.writeByte(tagLength[row]);
                fw.writeByte(multimaps[row]);
                fw.writeInt(chromosome[row]);
                fw.writeByte(strand[row]);
                fw.writeInt(startPosition[row]);
                fw.writeInt(endPosition[row]);
                fw.writeByte(divergence[row]);
                for (int j = 0; j < maxVariants; j++) {
                    fw.writeByte(variantPosOff[row][j]);
                    fw.writeByte(variantDef[row][j]);
                }
                fw.writeByte(dcoP[row]);
                fw.writeByte(mapP[row]);
                hapsOutput++;
            }
            fw.flush();
            fw.close();
            System.out.println("Tag positions written to:" + outFile.toString());
            System.out.println("Number of tags in file:" + hapsOutput);
        } catch (Exception e) {
            System.out.println("Catch in writing output file e=" + e);
        }
    }

    /**
     * @return An int[] result where : result[0] = The number of tags with a
     * unique physical positions in this file (i.e. , tags for which the
     * chromosome number is known). result[1] = The number of tags which align
     * to multiple positions (i.e., where multimaps[tagIndex] > 0)
     *
     */
    public int[] mappedTags() {
        int[] result = {0, 0};
        int unique = 0, multi = 1;  // the indices of result
        for (int row = 0; row < tagNum; row++) {
            if (chromosome[row] == Integer.MIN_VALUE) {
                if (multimaps[row] > 0) {
                    result[multi]++;
                }
            } else {
                result[unique]++;
            }
        }
        return result;
    }

    public void writeBinaryFile(File outFile) {
        writeBinaryFile(outFile, Integer.MAX_VALUE, false, false, Float.NaN, true);
    }

    protected void writeBinaryFile(File outFile, boolean binary) {
        writeBinaryFile(outFile, Integer.MAX_VALUE, false, false, Float.NaN, binary);
    }

    public void readTextFile(File inFile) {
        System.out.println("Reading tag alignment from:" + inFile.toString());
        String[] inputLine = {"NotRead"};
        try {
            BufferedReader br = new BufferedReader(new FileReader(inFile), 65536);
            inputLine = br.readLine().split("\t");

            //Parse header line (Number of tags, number of Long ints per tag, maximum variant bases per tag).
            this.tagNum = Integer.parseInt(inputLine[0]);
            this.tagLengthInLong = Integer.parseInt(inputLine[1]);
            this.maxVariants = Integer.parseInt(inputLine[2]);
            // initMatrices(9000000);
            initMatrices(tagNum);
            //Loop through remaining lines, store contents in a series of arrays indexed by row number
            for (int row = 0; row < tagNum; row++) {
                inputLine = br.readLine().split("\t");
                int c = 0;
                long[] tt = BaseEncoder.getLongArrayFromSeq(inputLine[c++]);
                for (int j = 0; j < tt.length; j++) {
                    tags[j][row] = tt[j];
                }
                tagLength[row] = parseByteWMissing(inputLine[c++]);
                multimaps[row] = parseByteWMissing(inputLine[c++]);
                chromosome[row] = parseIntWMissing(inputLine[c++]);
                strand[row] = parseByteWMissing(inputLine[c++]);
                startPosition[row] = parseIntWMissing(inputLine[c++]);
                endPosition[row] = parseIntWMissing(inputLine[c++]);
                divergence[row] = parseByteWMissing(inputLine[c++]);
                for (int j = 0; j < maxVariants; j++) {
                    variantPosOff[row][j] = parseByteWMissing(inputLine[c++]);
                    variantDef[row][j] = parseByteWMissing(inputLine[c++]);
                }
                dcoP[row] = parseByteWMissing(inputLine[c++]);
                mapP[row] = parseByteWMissing(inputLine[c++]);
                if (row % 1000000 == 0) {
                    System.out.println("Row Read:" + row);
                }
            }
        } catch (Exception e) {
            System.out.println("Catch in reading TagOnPhysicalMap file e=" + e);
            e.printStackTrace();
            System.out.println(Arrays.toString(inputLine));
        }
        System.out.println("Number of tags in file:" + tagNum);
    }

    public void writeTextFile(File outfile) {
        try {
            DataOutputStream fw = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outfile), 65536));
            fw.writeBytes(tagNum + "\t" + tagLengthInLong + "\t" + maxVariants + "\n");
            for (int row = 0; row < tagNum; row++) {
                fw.writeBytes(printRow(row) + "\n");
            }
            fw.flush();
            fw.close();
        } catch (Exception e) {
            System.out.println("Catch in writeTextFile file e=" + e);
        }
        System.out.println("Number of tags in file:" + tagNum);

    }

    @Override
    public int getReadIndexForPositionIndex(int posIndex) {
        return indicesOfSortByPosition[posIndex];
    }

    @Override
    public int[] getPositionArray(int index) {
        int[] r = {chromosome[index], strand[index], startPosition[index]};
        return r;
    }

//    public int getReadIndex(byte chr, byte strand, int posMin) {
//   //     dep_ReadWithPositionInfo querySHP=new dep_ReadWithPositionInfo(chromosome, strand, posMin);
//   //     PositionComparator posCompare = new PositionComparator();
//    //    int hit1=Arrays.binarySearch(shp, querySHP, posCompare);
//        int hit1=-1;  //redo
//        return hit1;
//    }
    public static String printWithMissing(byte b) {
        if (b == Byte.MIN_VALUE) {
            return "*";
        }
        return Byte.toString(b);
    }

    public static String printWithMissing(int i) {
        if (i == Integer.MIN_VALUE) {
            return "*";
        }
        return Integer.toString(i);
    }

    public static byte parseByteWMissing(String s) {
        if (s.equals("*")) {
            return Byte.MIN_VALUE;
        }
        int i;
        try {
            i = Integer.parseInt(s);
            if (i > 127) {
                return 127;
            }
            return (byte) i;
        } catch (NumberFormatException nf) {
            return Byte.MIN_VALUE;
        }
    }

    public static int parseIntWMissing(String s) {
        if (s.equals("*")) {
            return Integer.MIN_VALUE;
        }
        int i;
        try {
            i = Integer.parseInt(s);
            return i;
        } catch (NumberFormatException nf) {
            return Integer.MIN_VALUE;
        }
    }

    @Override
    public int getSize() {
        return tagNum;
    }

    @Override
    public byte getMultiMaps(int index) {
        return multimaps[index];
    }
    
    @Override
    public int getChromosome(int index) {
        return chromosome[index];
    }

    @Override
    public int getStartPosition(int index) {
        return startPosition[index];
    }

    @Override
    public int getEndPosition(int index) {
        return endPosition[index];
    }

    @Override
    public byte getDivergence(int index) {
        return divergence[index];
    }

    @Override
    public byte getMapP(int index) {
        return mapP[index];
    }

    @Override
    public byte getDcoP(int index) {
        return dcoP[index];
    }

    public byte[][] getVariantOff() {
        byte[][] result = new byte[getTagCount()][maxVariants];
        for (int i = 0; i < getTagCount(); i++) {
            System.arraycopy(variantPosOff[i], 0, result[i], 0, maxVariants);
        }
        return result;
    }

    public byte[][] getVariantDef() {
        byte[][] result = new byte[getTagCount()][maxVariants];
        for (int i = 0; i < getTagCount(); i++) {
            for (int j = 0; j < maxVariants; j++) {
                result[i][j] = getVariantDef(i, j);
            }
        }
        return result;
    }

    @Override
    public void setChromoPosition(int index, int chromosome, byte strand, int positionMin,
            int positionMax) {
        this.chromosome[index] = chromosome;
        this.strand[index] = strand;
        this.startPosition[index] = positionMin;
        this.endPosition[index] = positionMax;
    }

    @Override
    public void setDivergence(int index, byte divergence) {
        this.divergence[index] = divergence;
    }

    @Override
    public void setMapP(int index, byte mapP) {
        this.mapP[index] = mapP;
    }

    @Override
    public void setMapP(int index, double mapP) {
        if (Double.isInfinite(mapP)) {
            this.mapP[index] = Byte.MAX_VALUE;
            return;
        }
        if (Double.isNaN(mapP) || (mapP < 0) || (mapP > 1)) {
            this.mapP[index] = Byte.MIN_VALUE;
            return;
        }
        if (mapP < 1e-126) {
            this.mapP[index] = Byte.MAX_VALUE;
            return;
        }
        this.mapP[index] = (byte) (-Math.round(Math.log10(mapP)));
    }

    @Override
    public int addVariant(int tagIndex, byte offset, byte base) {
        for (int i = 0; i < maxVariants; i++) {
            if ((variantPosOff[tagIndex][i] == Byte.MIN_VALUE) && (variantDef[tagIndex][i] == Byte.MIN_VALUE)) {
                variantPosOff[tagIndex][i] = offset;
                variantDef[tagIndex][i] = base;
                return i;
            }
        }
        return -1; //no free space
    }

    /**
     * Decodes bitwise flags from the code in SAM field 3. The code is a 32-bit
     * integer, so I boolean AND the flag value with the code, and compare the
     * result with the flag value. If they are equal, that means all bits set in
     * the flag number were set in the code (i.e., that flag is turned on).
     */
    private boolean flagSet(int code, int flag) {
        int flagValue = 1 << flag; //1<<flag is equivalent to 2^flag
        return ((code & flagValue) == flagValue);
    }

    /**
     * Reads SAM files output from BWA or bowtie2
     */
    public void readSAMFile(String inputFileName, int tagLengthInLong) {
        System.out.println("Reading SAM format tag alignment from: " + inputFileName);
        this.tagLengthInLong = tagLengthInLong;
        String inputStr = "Nothing has been read from the file yet";
        int nHeaderLines = countTagsInSAMfile(inputFileName); // detects if the file is bowtie2, initializes topm matrices
        int tagIndex = Integer.MIN_VALUE;
        try {
            BufferedReader br;
            if (inputFileName.endsWith(".gz")) {
                br = new BufferedReader(new InputStreamReader(new MultiMemberGZIPInputStream(new FileInputStream(new File(inputFileName)))));
            } else {
                br = new BufferedReader(new FileReader(new File(inputFileName)), 65536);
            }
            for (int i = 0; i < nHeaderLines; i++) {
                br.readLine();
            } // Skip over the header
            for (tagIndex = 0; tagIndex < tagNum; tagIndex++) {
                inputStr = br.readLine();
                parseSAMAlignment(inputStr, tagIndex);
                if (tagIndex % 1000000 == 0) {
                    System.out.println("Read " + tagIndex + " tags.");
                }
            }
            br.close();
        } catch (Exception e) {
            System.out.println("\n\nCatch in reading SAM alignment file at tag " + tagIndex + ":\n\t" + inputStr + "\nError: " + e + "\n\n");
            e.printStackTrace();
            System.exit(1);
        }
    }

    private int countTagsInSAMfile(String inputFileName) {
        mySAMFormat = SAMFormat.BWA;  // format is BWA by default
        tagNum = 0;
        int nHeaderLines = 0;
        String currLine = null;
        try {
            String[] inputLine;
            ArrayList<String> chrNames = new ArrayList<String>();
            BufferedReader br;
            if (inputFileName.endsWith(".gz")) {
                br = new BufferedReader(new InputStreamReader(new MultiMemberGZIPInputStream(new FileInputStream(new File(inputFileName)))));
            } else {
                br = new BufferedReader(new FileReader(new File(inputFileName)), 65536);
            }
            while ((currLine = br.readLine()) != null) {
                inputLine = currLine.split("\\s");
                if (inputLine[0].contains("@")) {
                    //SAM files produced by Bowtie2 contain the string "@PG     ID:bowtie2      PN:bowtie2 "
                    if (inputLine[1].contains("bowtie2")) {
                        mySAMFormat = SAMFormat.BOWTIE2;
                    }
                    nHeaderLines++;
                } else {
                    String chr = inputLine[2];
                    if (!chrNames.contains(chr)) {
                        chrNames.add(chr);
                    }
                    tagNum++;
                    if (tagNum % 1000000 == 0) {
                        System.out.println("Counted " + tagNum + " tags.");
                    }
                }
            }
            br.close();
            System.out.println("Found " + tagNum + " tags in SAM file.  Assuming " + mySAMFormat + " file format.");
        } catch (Exception e) {
            System.out.println("Catch in counting lines of alignment file at line " + currLine + ": " + e);
            e.printStackTrace();
            System.exit(1);
        }
        initMatrices(tagNum);
        return nHeaderLines;
    }

    private void parseSAMAlignment(String inputStr, int tagIndex) {
        String[] inputLine = inputStr.split("\t");
        int name = 0, flag = 1, chr = 2, pos = 3, cigar = 5, tagS = 9; // column indices in inputLine
        String nullS = this.getNullTag();
        if ((Integer.parseInt(inputLine[flag]) & 4) == 4) {  // bit 0x4 (= 2^2 = 4) is set: NO ALIGNMENT
            recordLackOfSAMAlign(tagIndex, inputLine[tagS], inputLine[name], nullS);
        } else {  // aligns to one or more positions
            HashMap<String, Integer> SAMFields = parseOptionalFieldsFromSAMAlignment(inputLine);
            byte bestHits = (byte) Math.min(SAMFields.get("nBestHits"), Byte.MAX_VALUE);
            byte editDist = (byte) Math.min(SAMFields.get("editDist"), Byte.MAX_VALUE);
            byte currStrand = (Integer.parseInt(inputLine[flag]) == 16) ? (byte) -1 : (byte) 1;
            recordSAMAlign(
                    tagIndex,
                    inputLine[tagS],
                    inputLine[name],
                    nullS,
                    bestHits,
                    inputLine[chr],
                    currStrand,
                    Integer.parseInt(inputLine[pos]),
                    inputLine[cigar],
                    editDist);
        }
    }

    private HashMap<String, Integer> parseOptionalFieldsFromSAMAlignment(String[] inputLine) {
        HashMap<String, Integer> SAMFields = new HashMap<String, Integer>();
        if (mySAMFormat == SAMFormat.BWA) {
            for (int field = 11; field < inputLine.length; field++) { // Loop through all the optional field of the SAM alignment
                if (inputLine[field].regionMatches(0, "X0", 0, 2)) {        // X0 = SAM format for # of "high-quality" alignments of this query.  Specific to BWA.
                    SAMFields.put("nBestHits", Integer.parseInt(inputLine[field].split(":")[2]));
                } else if (inputLine[field].regionMatches(0, "NM", 0, 2)) { // NM = SAM format for edit distance to the reference.  Common to BWA and Bowtie2.
                    SAMFields.put("editDist", Integer.parseInt(inputLine[field].split(":")[2]));
                }
            }
        } else {  // bowtie2 -M format
            for (int field = 11; field < inputLine.length; field++) { // Loop through all the optional field of the SAM alignment
                if (inputLine[field].regionMatches(0, "AS", 0, 2)) {        // AS = SAM format for alignment score of the best alignment.  Specific to bowtie2.
                    SAMFields.put("bestScore", Integer.parseInt(inputLine[field].split(":")[2]));
                } else if (inputLine[field].regionMatches(0, "XS", 0, 2)) { // XS = SAM format for alignment score of 2nd best alignment.  Specific to bowtie2.
                    SAMFields.put("nextScore", Integer.parseInt(inputLine[field].split(":")[2]));
                } else if (inputLine[field].regionMatches(0, "NM", 0, 2)) { // NM = SAM format for edit distance to the reference.  Common to BWA and Bowtie2.
                    SAMFields.put("editDist", Integer.parseInt(inputLine[field].split(":")[2]));
                }
            }
            if (SAMFields.containsKey("bestScore")) {
                if (SAMFields.containsKey("nextScore")) {
                    if (SAMFields.get("bestScore") > SAMFields.get("nextScore")) {
                        SAMFields.put("nBestHits", 1);
                    } else {
                        SAMFields.put("nBestHits", 99);  // 99 will stand for an unknown # of multiple hits
                    }
                } else {
                    SAMFields.put("nBestHits", 1);
                }
            }
        }
        return SAMFields;
    }

    private void recordLackOfSAMAlign(int tagIndex, String tagS, String tagName, String nullS) {
        recordTagFromSAMAlignment(tagIndex, tagS, tagName, nullS, (byte) 1);
        multimaps[tagIndex] = 0;   // or should this be unknown = Byte.MIN_VALUE???
        chromosome[tagIndex] = Integer.MIN_VALUE;
        strand[tagIndex] = Byte.MIN_VALUE;
        startPosition[tagIndex] = Integer.MIN_VALUE;
        endPosition[tagIndex] = Integer.MIN_VALUE;
        divergence[tagIndex] = Byte.MIN_VALUE;
        for (int var = 0; var < maxVariants; var++) {
            variantPosOff[tagIndex][var] = Byte.MIN_VALUE;
            variantDef[tagIndex][var] = Byte.MIN_VALUE;
        }
        dcoP[tagIndex] = Byte.MIN_VALUE;
        mapP[tagIndex] = Byte.MIN_VALUE;
    }

    private void recordSAMAlign(int tagIndex, String tagS, String tagName, String nullS, byte nBestHits, String chrS, byte strand, int pos, String cigar, byte editDist) {
        recordTagFromSAMAlignment(tagIndex, tagS, tagName, nullS, strand);
        multimaps[tagIndex] = nBestHits;
        if (nBestHits == 1) {
            chromosome[tagIndex] = parseChrString(chrS);
            this.strand[tagIndex] = strand;
            recordStartEndPostionFromSAMAlign(tagIndex, strand, pos, cigar);
        } else {
            chromosome[tagIndex] = Integer.MIN_VALUE;
            this.strand[tagIndex] = Byte.MIN_VALUE;
            startPosition[tagIndex] = Integer.MIN_VALUE;
            endPosition[tagIndex] = Integer.MIN_VALUE;
        }
        divergence[tagIndex] = editDist;
        for (int var = 0; var < maxVariants; var++) {
            variantPosOff[tagIndex][var] = Byte.MIN_VALUE;
            variantDef[tagIndex][var] = Byte.MIN_VALUE;
        }
        dcoP[tagIndex] = Byte.MIN_VALUE;
        mapP[tagIndex] = Byte.MIN_VALUE;
    }

    private void recordTagFromSAMAlignment(int tagIndex, String tagS, String tagName, String nullS, byte strand) {
        if (strand == -1) {
            tagS = BaseEncoder.getReverseComplement(tagS);
        }
        if (tagS.length() < tagLengthInLong * 32) {  // pad with polyA
            tagS = tagS + nullS;
            tagS = tagS.substring(0, (tagLengthInLong * 32));
        }
        long[] tagSequence = BaseEncoder.getLongArrayFromSeq(tagS);
        for (int chunk = 0; chunk < tagLengthInLong; chunk++) {
            tags[chunk][tagIndex] = tagSequence[chunk];
        }
        tagName = tagName.replaceFirst("count=[0-9]+", "");
        tagName = tagName.replaceFirst("length=", "");
        tagLength[tagIndex] = Byte.parseByte(tagName);
    }

    private int parseChrString(String chrS) {
        int chr = Integer.MIN_VALUE;
        chrS = chrS.replace("chr", "");
        try {
            chr = Integer.parseInt(chrS);
        } catch (NumberFormatException e) {
            System.out.println("\n\nSAMConverterPlugin detected a non-numeric chromosome name: " + chrS
                    + "\n\nPlease change the FASTA headers in your reference genome sequence to integers "
                    + "(>1, >2, >3, etc.) OR to 'chr' followed by an integer (>chr1, >chr2, >chr3, etc.)\n\n");
            System.exit(1);
        }
        return chr;
    }

    private void recordStartEndPostionFromSAMAlign(int tagIndex, byte strand, int pos, String cigar) {
        int[] alignSpan = SAMUtils.adjustCoordinates(cigar, pos);
        try {
            if (strand == 1) {
                startPosition[tagIndex] = alignSpan[0];
                endPosition[tagIndex] = alignSpan[1];
            } else if (strand == -1) {
                startPosition[tagIndex] = alignSpan[1];
                endPosition[tagIndex] = alignSpan[0];
            } else {
                throw new Exception("Unexpected value for strand: " + strand + "(expect 1 or -1)");
            }
        } catch (Exception e) {
            System.out.println("Error in recordStartEndPostionFromSAMAlign: " + e);
            e.printStackTrace();
            System.exit(1);
        }
    }

    @Override
    public void swap(int index1, int index2) {
        long tl;
        for (int i = 0; i < tagLengthInLong; i++) {
            tl = tags[i][index1];
            tags[i][index1] = tags[i][index2];
            tags[i][index2] = tl;
        }
        int tb;
        tb = tagLength[index1];
        tagLength[index1] = tagLength[index2];
        tagLength[index2] = (byte) tb;
        tb = multimaps[index1];
        multimaps[index1] = multimaps[index2];
        multimaps[index2] = (byte) tb;
        tb = chromosome[index1];
        chromosome[index1] = chromosome[index2];
        chromosome[index2] = tb;
        tb = strand[index1];
        strand[index1] = strand[index2];
        strand[index2] = (byte) tb;
        int ti;
        ti = startPosition[index1];
        startPosition[index1] = startPosition[index2];
        startPosition[index2] = ti;
        ti = endPosition[index1];
        endPosition[index1] = endPosition[index2];
        endPosition[index2] = ti;
        tb = divergence[index1];
        divergence[index1] = divergence[index2];
        divergence[index2] = (byte) tb;
        for (int j = 0; j < maxVariants; j++) {
            tb = variantPosOff[index1][j];
            variantPosOff[index1][j] = variantPosOff[index2][j];
            variantPosOff[index2][j] = (byte) tb;
            tb = variantDef[index1][j];
            variantDef[index1][j] = variantDef[index2][j];
            variantDef[index2][j] = (byte) tb;
        }
        tb = dcoP[index1];
        dcoP[index1] = dcoP[index2];
        dcoP[index2] = (byte) tb;
        tb = mapP[index1];
        mapP[index1] = mapP[index2];
        mapP[index2] = (byte) tb;
    }

    @Override
    public int compare(int index1, int index2) {
        for (int i = 0; i < tagLengthInLong; i++) {
            if (tags[i][index1] < tags[i][index2]) {
                return -1;
            }
            if (tags[i][index1] > tags[i][index2]) {
                return 1;
            }
        }
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
        return 0;
    }

    /**
     * Returns an array whose <i>values</i> are the distinct chromosomes in this
     * file, as stored in the chromosome[] array. The indices are arbitrary.
     */
    @Override
    public int[] getChromosomes() {
        ArrayList<Integer> chrs = new ArrayList<Integer>();
        for (int i = 0; i < chromosome.length; i++) {
            if (!chrs.contains(chromosome[i]) && chromosome[i] != Integer.MIN_VALUE) {
                chrs.add(chromosome[i]);
                int x = 1;
            }
        }
        int[] result = new int[chrs.size()];
        for (int i = 0; i < chrs.size(); i++) {
            result[i] = chrs.get(i);
        }
        return result;
    }

    @Override
    public Locus[] getLoci() {
        int[] chrs = getChromosomes();
        Locus[] result = new Locus[chrs.length];

        for (int i = 0; i < result.length; i++) {
            result[i] = new Locus(chrs[i] + "", chrs[i] + "", -1, -1, null, null);
        }
        return result;
    }

    @Override
    public Locus getLocus(int tagIndex) {
        if (chromosome[tagIndex] == TOPMInterface.INT_MISSING) {
            return null;
        } //Return null for unmapped tags
        return new Locus(chromosome[tagIndex] + "", chromosome[tagIndex] + "", -1, -1, null, null);
    }

    /**
     * Returns an array containing all variant position offsets for the tag at
     * the supplied index.
     */
    @Override
    public byte[] getVariantPosOffArray(int tagIndex) {
        return variantPosOff[tagIndex];
    }

    /**
     * Returns an array containing all variant definitions for the tag at the
     * supplied index.
     */
    @Override
    public byte[] getVariantDefArray(int tagIndex) {
        byte[] result = new byte[maxVariants];
        for (int i = 0; i < maxVariants; i++) {
            result[i] = getVariantDef(tagIndex, i);
        }
        return result;
    }

    @Override
    public byte getVariantDef(int tagIndex, int variantIndex) {
        return asciiToNucleotideByte(variantDef[tagIndex][variantIndex]);
    }

    @Override
    public void setVariantDef(int tagIndex, int variantIndex, byte def) {
        variantDef[tagIndex][variantIndex] = (byte) NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES[0][def].charAt(0);
    }

    @Override
    public byte getVariantPosOff(int tagIndex, int variantIndex) {
        return variantPosOff[tagIndex][variantIndex];
    }

    @Override
    public void setVariantPosOff(int tagIndex, int variantIndex, byte offset) {
        variantPosOff[tagIndex][variantIndex] = offset;
    }

    @Override
    public int getMaxNumVariants() {
        return maxVariants;
    }

    @Override
    public byte getStrand(int tagIndex) {
        return strand[tagIndex];
    }

    @Override
    public int[] getUniquePositions(int chromosome) {
        HashSet<Integer> positions = new HashSet<Integer>();
//        int[] positions = new int[getTagCount()];
        for (int i = 0; i < getTagCount(); i++) {
            int tagIndex = getReadIndexForPositionIndex(i);
            if (getChromosome(tagIndex) != chromosome) {
                continue;
            }
            int startPos = getStartPosition(tagIndex);
            for (int varIndex = 0; varIndex < maxVariants; varIndex++) {
                int offset = getVariantPosOff(tagIndex, varIndex);
                if (offset == Byte.MIN_VALUE) {
                    continue;
                }
                int position = startPos + offset;
                positions.add(position);
            }
        }

        int[] result = new int[positions.size()];//new int[positions.size()];
        int i = 0;
        for (Integer position : positions) {
            result[i] = position.intValue();
            i++;
        }
        return result;
    }

    /**
     * First, find unique tags by concatenating them into a TagCountMutable
     * object and collapsing duplicate rows.
     */
    private static TagsOnPhysicalMap uniqueTags(String[] filenames) {

        //Find dimensions of concatenated file
        int tagLengthInLong = 0, maxVariants = 0, tagNum = 0;
        for (String name : filenames) {
            TagsOnPhysicalMap file = new TagsOnPhysicalMap(name, true);
            tagNum += file.tagNum;
            if (file.maxVariants > maxVariants) {
                maxVariants = file.maxVariants;
            }
            if (file.getTagSizeInLong() > tagLengthInLong) {
                tagLengthInLong = file.getTagSizeInLong();
            }
        }

        //Create new concatenated file
        TagCountMutable tc = new TagCountMutable(tagLengthInLong, tagNum);
        for (String name : filenames) {
            TagsOnPhysicalMap file = new TagsOnPhysicalMap(name, true);
            for (int i = 0; i < file.tagNum; i++) {
                tc.addReadCount(file.getTag(i), file.getTagLength(i), 1);
            }
        }

        //Reduce concatenated file to unique tags
        tc.collapseCounts();
        tc.shrinkToCurrentRows();

        TagsOnPhysicalMap result = new TagsOnPhysicalMap(tc);
        result.expandMaxVariants(maxVariants);
        result.clearVariants();

        return result;
    }

    /**
     * Merges several TOPM files into one, removing duplicate rows. Variants
     * from matching tags are combined in the output file. If there are more
     * variants in the input files than can fit in the output, extra variants
     * are discarded with no particular order.
     */
    public static TagsOnPhysicalMap merge(String[] filenames) {

        TagsOnPhysicalMap output = uniqueTags(filenames);

        System.out.println(
                "Output file will contain "
                + output.tagNum + " unique tags, "
                + output.getTagSizeInLong() + " longs/tag, "
                + output.maxVariants + " variants. ");

        for (String name : filenames) {
            TagsOnPhysicalMap file = new TagsOnPhysicalMap(name, true);

            int varsAdded = 0, varsSkipped = 0;
            for (int inTag = 0; inTag < file.tagNum; inTag++) {  //Loop over tags in input

                //Find index of corresponding tag in output file, and copy attributes from input to output
                int outTag = output.getTagIndex(file.getTag(inTag));
                copyTagAttributes(file, inTag, output, outTag);

                //For corresponding tags, compare variants
                for (int outVar = 0; outVar < output.maxVariants; outVar++) {

                    byte outOff = output.getVariantPosOff(outTag, outVar);
                    if (outOff != BYTE_MISSING) {//Skip filled output variants or re-initialize them
                        varsSkipped++;
                        continue;
                    }

                    for (int inVar = 0; inVar < file.maxVariants; inVar++) {
                        byte offset = file.getVariantPosOff(inTag, outVar);
                        byte def = file.getVariantDef(inTag, outVar);
                        if (offset == BYTE_MISSING) {
                            continue;                            //Skip blank input variants
                        }
                        //If we get here, output variant is blank and input variant is non-blank at the same tag & variant indices
                        varsAdded++;
                        output.setVariantPosOff(outTag, outVar, offset);
                        output.setVariantDef(outTag, outVar, def);
                        file.setVariantPosOff(inTag, inVar, Byte.MIN_VALUE);  //Erase this variant so it isn't encountered again
                        break; //Go to next output variant after copying
                    }
                }
            }
            System.out.println(varsAdded + " variants added.");
            System.out.println(varsSkipped + " variants skipped.");
        }
        return output;
    }

    /**
     * Copies values of everything BUT tag sequence and variant data (i.e. tag
     * attributes) from one TOPM file to another.
     */
    private static void copyTagAttributes(TagsOnPhysicalMap input, int inTag, TagsOnPhysicalMap output, int outTag) {
        output.tagLength[outTag] = input.tagLength[inTag];
        output.multimaps[outTag] = input.multimaps[inTag];
        output.chromosome[outTag] = input.chromosome[inTag];
        output.strand[outTag] = input.strand[inTag];
        output.startPosition[outTag] = input.startPosition[inTag];
        output.endPosition[outTag] = input.endPosition[inTag];
        output.divergence[outTag] = input.divergence[inTag];
        output.dcoP[outTag] = input.dcoP[inTag];
        output.mapP[outTag] = input.mapP[inTag];
    }

    /**
     * Fills the variant definition & offset arrays with the value of
     * "byteMissing".
     */
    public void clearVariants() {
        for (int i = 0; i < getTagCount(); i++) {
            Arrays.fill(variantDef[i], BYTE_MISSING);
            Arrays.fill(variantPosOff[i], BYTE_MISSING);
        }
    }

    /**
     * Fills the variant definition & offset at the given indices with the value
     * of "byteMissing".
     */
    private void clearVariant(int tag, int variant) {
        setVariantDef(tag, variant, BYTE_MISSING);
        setVariantPosOff(tag, variant, BYTE_MISSING);
    }

    /**
     * Clears variant sites that are not found in the supplied alignments.
     */
    public void filter(String[] filenames) {
        HashMap<String, Integer> hapmapSites = new HashMap<String, Integer>();

        //Map all site positions from all alignments to their index.
        for (String filename : filenames) {
            System.out.println("Filtering out sites from " + filename + ".");
            hapmapSites.putAll(hapmapSites(ImportUtils.readFromHapmap(filename, null)));
        }
        System.out.println("There are " + hapmapSites.size() + " sites in the hapmap files.");

        //Map all tag variant positions to their index.
        HashMap<String, Integer> topmSites = uniqueSites(); //Map tag variant positions to bases

        //Sites should be a subset of tag variants, so check that there are fewer of them.
        System.out.println("Found " + topmSites.size() + " unique sites in " + tagNum + " tags in TOPM.");
        if (topmSites.size() < hapmapSites.size()) {
            System.out.println("Warning: more unique sites exist in hapmap file.");
        }

        //Next, check that each alignment site is present in   If not, either there was no more room or something is wrong.
        HashSet<Integer> fullSites = fullTagPositions(); //Tags which already have the maximum number of variants
        ArrayList<String> insertedSites = new ArrayList<String>(); //Sites that don't correspond to a real tag
        ArrayList<String> skippedSites = new ArrayList<String>(); //Real sites that were skipped due to "full tags"

        int basePerTag = tagLengthInLong * 32;
        for (String hapmapSNP : hapmapSites.keySet().toArray(new String[hapmapSites.size()])) {

            int chr = Integer.parseInt(hapmapSNP.split("\\t")[0]);
            int pos = Integer.parseInt(hapmapSNP.split("\\t")[1]);
            if (topmSites.get(hapmapSNP) == null) {

//                System.out.print("Warning: SNP "+chr+":"+pos+" is not in TOPM.  ");

                boolean inRange = false;
                for (int i = -basePerTag; i < basePerTag; i++) {
                    if (fullSites.contains(i + pos)) {
                        inRange = true;
                        break;
                    }
                }

                if (inRange) {
                    skippedSites.add(hapmapSNP);
//                    System.out.println("However, it is within range of a tag with the max. number of variants.");
                } else {
                    insertedSites.add(hapmapSNP);
                    System.out.println();
                }
                hapmapSites.remove(hapmapSNP);
                continue;
            }
        }

        System.out.println("The following SNPs were not in the TOPM, but are within range of a tag with the max. number of variants:");
        for (String site : skippedSites) {
            System.out.println(site);
        }

        System.out.println("The following SNPs were not in the TOPM, and do not correspond to any known tag:");
        for (String site : insertedSites) {
            System.out.println(site);
        }

        //Remove any sites from the TOPM that are absent in the alignment
        int removedSites = 0;
        for (int tag = 0; tag < tagNum; tag++) {
            int chr = getChromosome(tag);
            int pos = getStartPosition(tag);
            for (int variant = 0; variant < maxVariants; variant++) {

                byte off = getVariantPosOff(tag, variant);
                String site = chr + "\t" + (pos + off);
                if (!hapmapSites.containsKey(site)) {
                    clearVariant(tag, variant);
                    removedSites++;
                }
            }
        }

        topmSites = uniqueSites();
        System.out.println("Removed " + removedSites + " TOPM sites not present in alignment and ignored " + insertedSites.size() + " alignment sites not present in TOPM.");
        System.out.println("There are " + topmSites.size() + " sites in the TOPM now, as compared to " + hapmapSites.size() + " sites in the alignment.");
        if (topmSites.size() != hapmapSites.size()) {
            System.out.println("Warning: number of filtered sites does not match number of alignment sites.");
        }
    }

    /**
     * Returns the start positions of tags whose variant arrays are full.
     */
    public HashSet<Integer> fullTagPositions() {
        HashSet<Integer> result = new HashSet<Integer>();
        for (int i = 0; i < tagNum; i++) {
            boolean variantsFull = true;
            for (int j = 0; j < maxVariants; j++) {
                byte off = getVariantPosOff(i, j);
                if (off == Byte.MIN_VALUE) {
                    variantsFull = false;
                    break;
                }
            }

            if (variantsFull) {
                result.add(getStartPosition(i));
            }
        }
        return result;
    }

    /**
     * Maps unique sites (i.e. chromosome and position) to the indices of the
     * tags in which they are found.
     */
    public HashMap<String, Integer> uniqueSites() {
        HashMap<String, Integer> snps = new HashMap<String, Integer>();
        for (int tag = 0; tag < tagNum; tag++) {        //Visit each tag in TOPM

            for (int variant = 0; variant < maxVariants; variant++) {                //Visit each variant in TOPM
                byte off = getVariantPosOff(tag, variant);
                if (off == BYTE_MISSING) {
                    continue;
                }

                int a = getStartPosition(tag);
                int b = getVariantPosOff(tag, variant);
                int c = a + b;
                String pos = getChromosome(tag) + "\t" + c;

                snps.put(pos, tag);
            }
        }
        return snps;
    }

    public static HashMap<String, Integer> hapmapSites(Alignment file) {
        HashMap<String, Integer> snps = new HashMap<String, Integer>();
        for (int site = 0; site < file.getSiteCount(); site++) {        //Visit each site
            String pos = (file.getLocusName(site) + "\t"
                    + file.getPositionInLocus(site));
            if (file.getPositionInLocus(site) > 2000000000) {
                System.out.println(pos);
            }
            snps.put(pos, site);
        }
        return snps;
    }

    /**
     * Returns an array whose indices are the number of mappings and whose
     * elements are the number of tags with that mapping.
     */
    public int[] mappingDistribution() {
        int[] result = new int[128]; //Only up to 127 multiple mappings are stored.
        for (int i = 0; i < tagNum; i++) {
            if (multimaps[i] > (result.length - 1)) {
                result[127]++;
            }
            if (multimaps[i] == BYTE_MISSING) {
                result[0]++;
                continue;
            } else {
                result[multimaps[i]]++;
            }
        }
        return result;
    }

    private byte asciiToNucleotideByte(byte ascii) {
        switch (ascii) {
            case 'A':
                return NucleotideAlignmentConstants.A_ALLELE;
            case 'C':
                return NucleotideAlignmentConstants.C_ALLELE;
            case 'G':
                return NucleotideAlignmentConstants.G_ALLELE;
            case 'T':
                return NucleotideAlignmentConstants.T_ALLELE;
            case 'N':
                return Alignment.UNKNOWN_ALLELE;
            case '-':
                return NucleotideAlignmentConstants.GAP_ALLELE;
            case '+':
                return NucleotideAlignmentConstants.INSERT_ALLELE;
            case BYTE_MISSING:
                return BYTE_MISSING;
            default:
                throw new IllegalArgumentException("asciiToNucleotideByte: unknown ascii: " + Byte.toString(ascii));
        }
    }
}
