/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.gbs.maps;

import cern.colt.GenericSorting;
import cern.colt.Swapper;
import cern.colt.function.IntComparator;
import java.io.BufferedOutputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.util.Map.Entry;
import java.util.TreeMap;
import java.util.TreeSet;
import net.maizegenetics.gbs.tagdist.AbstractTags;
import net.maizegenetics.gbs.util.BaseEncoder;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.Locus;
import net.maizegenetics.pal.alignment.NucleotideAlignmentConstants;

/**
 *
 * @author edbuckler
 */
public abstract class AbstractTagsOnPhysicalMap extends AbstractTags implements TOPMInterface {
    protected int[] bestChr; // 4 bytes
    // 4 bytes
    //if these disagree with the location, then set the p to negative
    // 1+4+1+4+4+1+8+8+1+1 = 33 bytes per position + 16 bytes for a two long tag + 1 byte for tagLength in bases = 50 bytes
    // ~50 bytes per position.  If we have 10 million tags then this will be at 500M byte data structure.
    protected int[] indicesOfSortByPosition;
    protected int myMaxVariants = 8;
    protected byte[] multimaps; // number of locations this tagSet maps to; unknown = Byte.MIN_VALUE; multiple, but unknown number = 99
    // number of locations this tagSet maps to; unknown = Byte.MIN_VALUE; multiple, but unknown number = 99
    protected int[] bestStartPos; // chromosomal position of the barcoded end of the tag  // 4 bytes
    // chromosomal position of the barcoded end of the tag  // 4 bytes
    protected byte[] bestStrand; // 1 = same sense as reference FASTA file.  -1 = opposite sense.  unknown = Byte.MIN_VALUE  // 1 byte
    // 1 = same sense as reference FASTA file.  -1 = opposite sense.  unknown = Byte.MIN_VALUE  // 1 byte
    protected int myNumTags = 0;
    //TODO
    public int tagNum;  //remove this and set to above
    public int maxVariants;
    
    protected byte[][] variantDefs; // allele state - A, C, G, T or some indel definition  // myMaxVariants bytes [tag][variant]
    // allele state - A, C, G, T or some indel definition  // myMaxVariants bytes [tag][variant]
    protected byte[][] variantOffsets; // offset from position minimum, maximum number of variants is defined above  // myMaxVariants bytes [tag][variant]
    // offset from position minimum, maximum number of variants is defined above  // myMaxVariants bytes [tag][variant]
    protected int[] myChromosomes = null; //sort ascending
    protected int[][] myUniquePositions = null; //dimensions [chromosome][positions] note ragged, and sort ascending

    public AbstractTagsOnPhysicalMap() {
    }

    
    @Override
    public int getMaxNumVariants() {
        return myMaxVariants;
    }
    
    @Override
    public byte getMultiMaps(int index) {
        return multimaps[index];
    }
    
    @Override
    public int getChromosome(int index) {
        return bestChr[index];
    }

    @Override
    public int getSize() {
        return myNumTags;
    }

    @Override
    public int getStartPosition(int index) {
        return bestStartPos[index];
    }

    @Override
    public byte getStrand(int tagIndex) {
        return bestStrand[tagIndex];
    }

    public byte[][] getVariantDef() {
        byte[][] result = new byte[getTagCount()][myMaxVariants];
        for (int i = 0; i < getTagCount(); i++) {
            for (int j = 0; j < myMaxVariants; j++) {
                result[i][j] = getVariantDef(i, j);
            }
        }
        return result;
    }
    

    @Override
    public byte getVariantDef(int tagIndex, int variantIndex) {
        if((variantDefs[tagIndex]==null)||(variantDefs[tagIndex].length<=variantIndex)) return TOPMInterface.BYTE_MISSING;
        return variantDefs[tagIndex][variantIndex];
    }

    /**
     * Returns an array containing all variant definitions for the tag at the
     * supplied index.
     */
    @Override
    public byte[] getVariantDefArray(int tagIndex) {
        if(variantDefs[tagIndex]==null) return null;
        byte[] result = new byte[variantDefs[tagIndex].length];
        for (int i = 0; i < variantDefs[tagIndex].length; i++) {
            result[i] = getVariantDef(tagIndex, i);
        }
        return result;
    }

    @Override
    public byte[][] getVariantOff() {
        byte[][] result = new byte[getTagCount()][myMaxVariants];
        for (int i = 0; i < getTagCount(); i++) {
            System.arraycopy(variantOffsets[i], 0, result[i], 0, myMaxVariants);
        }
        return result;
    }

    @Override
    public byte getVariantPosOff(int tagIndex, int variantIndex) {
        if((variantOffsets[tagIndex]==null)||(variantOffsets[tagIndex].length<=variantIndex)) return TOPMInterface.BYTE_MISSING;
        return variantOffsets[tagIndex][variantIndex];
    }

    /**
     * Returns an array containing all variant position offsets for the tag at
     * the supplied index.
     */
    @Override
    public byte[] getVariantPosOffArray(int tagIndex) {
        return variantOffsets[tagIndex];
    }
    
    public String printRow(int row) {
        StringBuilder sb = new StringBuilder();
        sb.append(sb);
        //long
        sb.append(BaseEncoder.getSequenceFromLong(this.getTag(row)) + "\t");
        sb.append(printWithMissing(tagLength[row]) + "\t");
        sb.append(printWithMissing(multimaps[row]) + "\t");
        sb.append(printWithMissing(bestChr[row]) + "\t");
        sb.append(printWithMissing(bestStrand[row]) + "\t");
        sb.append(printWithMissing(bestStartPos[row]) + "\t");
        sb.append(printWithMissing(getEndPosition(row)) + "\t");
        sb.append(printWithMissing(getDivergence(row)) + "\t");
        for (int j = 0; j < myMaxVariants; j++) {
            sb.append(printWithMissing(getVariantPosOff(row, j)) + "\t");
            byte vd=getVariantDef(row, j);
            if(vd==TOPMInterface.BYTE_MISSING) {sb.append(printWithMissing(vd) + "\t");}
            else {sb.append(NucleotideAlignmentConstants.getNucleotideIUPAC(vd) + "\t");}
        }
        sb.append(printWithMissing(getDcoP(row)) + "\t");
        sb.append(printWithMissing(getMapP(row)) + "\t");
 //       System.out.println("Line:"+row+":"+sb.toString());
        return sb.toString();
    }
    
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

    public void writeTextFile(File outfile) {
        try {
            DataOutputStream fw = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outfile), 65536));
            fw.writeBytes(myNumTags + "\t" + tagLengthInLong + "\t" + myMaxVariants + "\n");
            for (int row = 0; row < myNumTags; row++) {
                fw.writeBytes(printRow(row) + "\n");
            }
            fw.flush();
            fw.close();
        } catch (Exception e) {
            System.out.println("Catch in writeTextFile file e=" + e);
            e.printStackTrace();
        }
        System.out.println("Number of tags in file:" + myNumTags);
    }

    @Override
    public int[] getChromosomes() {
        if(myChromosomes==null) populateChrAndVarPositions();
        return myChromosomes;
    }
        
    @Override
    public int[] getUniquePositions(int chromosome) {
        if(myUniquePositions==null) populateChrAndVarPositions();
        return myUniquePositions[chromosome];
    }
    
    protected void populateChrAndVarPositions() {
        long chrSum = 0;
        System.out.println("chrSum" + chrSum);
        TreeMap<Integer, TreeSet<Integer>> theChrs = new TreeMap<Integer, TreeSet<Integer>>();
        for (int i = 0; i < myNumTags; i++) {
            int chr = getChromosome(i);
            if (chr != TOPMInterface.INT_MISSING) {
                if (!theChrs.containsKey(chr)) {
                    theChrs.put(chr, new TreeSet<Integer>());
                }
                TreeSet<Integer> thePos = theChrs.get(chr);
                int startPos = getStartPosition(i);
                byte[] varOffs = getVariantPosOffArray(i);
                for (byte b : varOffs) {
                    thePos.add((int) (startPos + b));
                }
            }
        }
        myChromosomes = new int[theChrs.size()];
        myUniquePositions = new int[theChrs.size()][];
        int cnt = 0;
        for (Entry<Integer, TreeSet<Integer>> aChr : theChrs.entrySet()) {
            myUniquePositions[cnt] = new int[aChr.getValue().size()];
            int p = 0;
            for (int ls : aChr.getValue()) {
                myUniquePositions[cnt][p++] = ls;
            }
            myChromosomes[cnt++] = aChr.getKey();
            System.out.printf("Chr:%d TagStart:%d %n", myChromosomes[cnt - 1], myUniquePositions[cnt - 1].length);
        }
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
        if (bestChr[tagIndex] == TOPMInterface.INT_MISSING) {
            return null;
        } //Return null for unmapped tags
        return new Locus(bestChr[tagIndex] + "", bestChr[tagIndex] + "", -1, -1, null, null);
    }
    
   void initPhysicalSort() {
        System.out.println("initPhysicalSort");
        indicesOfSortByPosition = new int[myNumTags];
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
                if (bestChr[index1] < bestChr[index2]) {
                    return -1;
                }
                if (bestChr[index1] > bestChr[index2]) {
                    return 1;
                }
                if (bestStartPos[index1] < bestStartPos[index2]) {
                    return -1;
                }
                if (bestStartPos[index1] > bestStartPos[index2]) {
                    return 1;
                }
                if (bestStrand[index1] < bestStrand[index2]) {
                    return -1;
                }
                if (bestStrand[index1] > bestStrand[index2]) {
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
    
}
