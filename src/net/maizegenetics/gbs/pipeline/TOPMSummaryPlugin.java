/*
 * TOPMSummaryPlugin
 */
package net.maizegenetics.gbs.pipeline;

import java.awt.Frame;
import java.io.BufferedWriter;
import java.io.File;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import javax.swing.ImageIcon;
import net.maizegenetics.gbs.maps.TagsOnPhysicalMap;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.util.ArgsEngine;
import net.maizegenetics.util.Utils;
import org.apache.log4j.Logger;

/**
 *
 * @author terry
 */
public class TOPMSummaryPlugin extends AbstractPlugin {

    private final Logger myLogger = Logger.getLogger(TOPMSummaryPlugin.class);
    private ArgsEngine myArgsEngine = null;
    private String myInputFilename = null;
    private TagsOnPhysicalMap myInputTOPM = null;
    private int myTagCount = 0;
    private int[] myChromosomes;
    private Map<Integer, Integer>[] myVariantsPerPosition;
    private Map<Integer, Set<Byte>>[] myVariantDefsPerPosition;
    private int myNumUndefinedStrandedTags = 0;
    private Set<Byte> myUndefinedStrandValues = new HashSet<Byte>();
    private String myOutputFilename = null;
    private int[] myNumSNPsPerChromosome;
    private int[] myNumTagsPerVariantsDefined;

    public TOPMSummaryPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    @Override
    public DataSet performFunction(DataSet input) {

        myInputTOPM = new TagsOnPhysicalMap(myInputFilename, true);
        myTagCount = myInputTOPM.getTagCount();
        myLogger.info("performFunction: Number of Tags: " + myTagCount);

        myChromosomes = myInputTOPM.getChromosomes();
        Arrays.sort(myChromosomes);

        myNumSNPsPerChromosome = new int[myChromosomes.length];

        myNumTagsPerVariantsDefined = new int[myInputTOPM.maxVariants + 1];

        myVariantsPerPosition = new TreeMap[myChromosomes.length];
        for (int m = 0; m < myChromosomes.length; m++) {
            myVariantsPerPosition[m] = new TreeMap<Integer, Integer>();
        }
        myVariantDefsPerPosition = new TreeMap[myChromosomes.length];
        for (int m = 0; m < myChromosomes.length; m++) {
            myVariantDefsPerPosition[m] = new TreeMap<Integer, Set<Byte>>();
        }

        for (int i = 0; i < myTagCount; i++) {
            int startPos = myInputTOPM.getStartPosition(i);
            int endPos = myInputTOPM.getEndPosition(i);
            byte strand = myInputTOPM.getStrand(i);
            int chrom = myInputTOPM.getChromosome(i);
            int index = Arrays.binarySearch(myChromosomes, chrom);

            int tagLength = myInputTOPM.getTagLength(i);
            //String tag = BaseEncoder.getSequenceFromLong(myInputTOPM.getTag(i));
            if (strand == 1) {
                if (index < 0) {
                    myLogger.error("performFunction: tag: " + i + " chromosome: " + chrom + " not reported by getChromosomes()");
                    continue;
                }
                if (startPos > endPos) {
                    myLogger.error("performFunction: tag: " + i + " invalid state: strand: " + strand + "  start position: " + startPos + "  end position: " + endPos);
                    continue;
                }
                //int startEndLength = endPos - startPos + 1;
                //if (startEndLength != tagLength) {
                //    myLogger.warn("performFunction: tag: " + i + " tag length: " + tagLength + " doesn't equal (end pos - start pos + 1): " + startEndLength);
                //}
                int numDefinedVariants = 0;
                for (int j = 0; j < myInputTOPM.maxVariants; j++) {
                    int offset = myInputTOPM.getVariantPosOff(i, j);
                    byte def = myInputTOPM.getVariantDef(i, j);
                    //if (offset != Byte.MIN_VALUE) {
                    //    myLogger.info("performFunction: Tag: " + i + " Positive Strand: Defined Variant: Offset: " + offset + "  def: " + (char) def);
                    //}
                    if ((offset >= 0) && (def >= 0)) {
                        numDefinedVariants++;
                        myNumSNPsPerChromosome[index]++;
                        //if ((char) def != tag.charAt(offset)) {
                        //    myLogger.error("performFunction: Mismatch: Sequence From Long: " + tag + "  offset: " + offset + " def: " + (char) def + " from tag: " + tag.charAt(offset));
                        //}
                        int position = startPos + offset;
                        if (position > endPos) {
                            myLogger.warn("performFunction: tag: " + i + " tag length: " + tagLength + " on chromosome: " + chrom + " has invalid offset: " + offset + " puts physical postion: " + position + " outside range: " + startPos + " to " + endPos);
                        }
                        Integer count = myVariantsPerPosition[index].get(position);
                        if (count == null) {
                            myVariantsPerPosition[index].put(position, 1);
                            Set<Byte> temp = new HashSet<Byte>();
                            temp.add(def);
                            myVariantDefsPerPosition[index].put(position, temp);
                        } else {
                            myVariantsPerPosition[index].put(position, count + 1);
                            Set<Byte> temp = myVariantDefsPerPosition[index].get(position);
                            temp.add(def);
                        }
                    }
                }
                myNumTagsPerVariantsDefined[numDefinedVariants]++;
            } else if (strand == -1) {
                if (index < 0) {
                    myLogger.error("performFunction: tag: " + i + " chromosome: " + chrom + " not reported by getChromosomes()");
                    continue;
                }
                if (startPos < endPos) {
                    myLogger.error("performFunction: tag: " + i + " invalid state: strand: " + strand + "  start position: " + startPos + "  end position: " + endPos);
                    continue;
                }
                //int startEndLength = startPos - endPos + 1;
                //if (startEndLength != tagLength) {
                //    myLogger.warn("performFunction: tag: " + i + " tag length: " + tagLength + " doesn't equal (start pos - end pos + 1): " + startEndLength);
                //}
                int numDefinedVariants = 0;
                for (int j = 0; j < myInputTOPM.maxVariants; j++) {
                    int offset = myInputTOPM.getVariantPosOff(i, j);
                    byte def = myInputTOPM.getVariantDef(i, j);
                    //if (offset != Byte.MIN_VALUE) {
                    //    myLogger.info("performFunction: Tag: " + i + " Negative Strand: Defined Variant: Offset: " + offset + "  def: " + (char) def);
                    //}
                    if ((offset >= 0) && (def >= 0)) {
                        numDefinedVariants++;
                        myNumSNPsPerChromosome[index]++;
                        //if ((char) def != tag.charAt(offset)) {
                        //    myLogger.error("performFunction: Mismatch: Sequence From Long: " + tag + "  offset: " + offset + " def: " + (char) def + " from tag: " + tag.charAt(offset));
                        //}
                        int position = startPos + offset;
                        if (position < endPos) {
                            myLogger.warn("performFunction: tag: " + i + " tag length: " + tagLength + " on chromosome: " + chrom + " has invalid offset: " + offset + " puts physical postion: " + position + " outside range: " + startPos + " to " + endPos);
                        }
                        Integer count = myVariantsPerPosition[index].get(position);
                        if (count == null) {
                            myVariantsPerPosition[index].put(position, 1);
                            Set<Byte> temp = new HashSet<Byte>();
                            temp.add(def);
                            myVariantDefsPerPosition[index].put(position, temp);
                        } else {
                            myVariantsPerPosition[index].put(position, count + 1);
                            Set<Byte> temp = myVariantDefsPerPosition[index].get(position);
                            temp.add(def);
                        }
                    }
                }
                myNumTagsPerVariantsDefined[numDefinedVariants]++;
            } else {
                myNumUndefinedStrandedTags++;
                myUndefinedStrandValues.add(strand);
            }
        }

        myLogger.info("performFunction: Number of Tags with Undefined Strands: " + myNumUndefinedStrandedTags);
        Iterator itr = myUndefinedStrandValues.iterator();
        while (itr.hasNext()) {
            myLogger.info("performFunction: Undefined Strand Value: " + itr.next());
        }

        for (int i = 0; i < myChromosomes.length; i++) {
            myLogger.info("performFunction: Chromosome: " + myChromosomes[i] + " Number of SNPs: " + myNumSNPsPerChromosome[i]);
        }

        for (int i = 0; i <= myInputTOPM.maxVariants; i++) {
            myLogger.info("performFunction: Number of Tags: " + myNumTagsPerVariantsDefined[i] + " Has: " + i + " Variants Defined");
        }

        printSummary();
        return null;
    }

    private void printSummary() {
        BufferedWriter writer = null;

        try {
            writer = Utils.getBufferedWriter(myOutputFilename);
            writer.append("Chromosome\tPosition\tNum Variants\tVariant Defs\n");
            for (int c = 0; c < myChromosomes.length; c++) {
                Iterator itr = myVariantsPerPosition[c].entrySet().iterator();
                while (itr.hasNext()) {
                    Map.Entry entry = (Map.Entry) itr.next();
                    writer.append(myChromosomes[c] + "\t" + entry.getKey() + "\t" + entry.getValue() + "\t");
                    Set<Byte> defSet = myVariantDefsPerPosition[c].get(entry.getKey());
                    Iterator itr2 = defSet.iterator();
                    boolean notFirst = false;
                    while (itr2.hasNext()) {
                        if (notFirst) {
                            writer.append(",");
                        } else {
                            notFirst = true;
                        }
                        writer.append((char) ((Byte) itr2.next()).byteValue());
                    }
                    writer.append("\n");
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
        } finally {
            try {
                writer.close();
            } catch (Exception ex) {
                // do nothing
            }
        }
    }

    private void printUsage() {
        myLogger.info(
                "\nThe options for the TOPMSummaryPlugin:\n"
                + "-input Input TOPM\n"
                + "-output Output Filename\n");
    }

    @Override
    public void setParameters(String[] args) {

        if (args.length == 0) {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }

        if (myArgsEngine == null) {
            myArgsEngine = new ArgsEngine();
            myArgsEngine.add("-input", "-input", true);
            myArgsEngine.add("-output", "-output", true);
        }
        myArgsEngine.parse(args);

        myInputFilename = myArgsEngine.getString("-input");
        if ((myInputFilename == null) || (myInputFilename.length() == 0)) {
            printUsage();
            throw new IllegalArgumentException("TOPMSummaryPlugin: setParameters: Must define input file");
        }
        File inputFile = new File(myInputFilename);
        if (!inputFile.exists()) {
            printUsage();
            throw new IllegalArgumentException("TOPMSummaryPlugin: setParameters: The input file doesn't exist: " + myInputFilename);
        }

        myOutputFilename = myArgsEngine.getString("-output");
        if ((myOutputFilename == null) || (myOutputFilename.length() == 0)) {
            printUsage();
            throw new IllegalArgumentException("TOPMSummaryPlugin: setParameters: Must define output file");
        }
        File outputFile = new File(myOutputFilename);
        if (outputFile.exists()) {
            printUsage();
            throw new IllegalArgumentException("TOPMSummaryPlugin: setParameters: The output file already exist: " + myInputFilename);
        }

    }

    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String getToolTipText() {
        throw new UnsupportedOperationException("Not supported yet.");
    }
}
