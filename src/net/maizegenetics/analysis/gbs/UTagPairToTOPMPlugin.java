/*
 * Plugin for the UNEAK pipeline to convert the Tag Pair file into a (fake) TOPM file for reintegration with the normal GBS pipeline
 */
package net.maizegenetics.analysis.gbs;

import net.maizegenetics.dna.tag.UTagPairs;
import java.awt.Frame;
import java.io.File;
import javax.swing.ImageIcon;
import net.maizegenetics.dna.map.TagsOnPhysicalMap;
import net.maizegenetics.dna.tag.AbstractTags;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.util.ArgsEngine;
import org.apache.log4j.Logger;

/**
 *
 * @author Jason Wallace
 *
 * This plugin takes the TagPairs file produced by earlier steps in the UNEAK
 * pipeline and converts it to a TOPM alignment file that can be used by
 * DiscoverySNPCallerPlugin to call SNPs as if from a reference genome.
 *
 */
public class UTagPairToTOPMPlugin extends AbstractPlugin {

    private ArgsEngine engine = null;
    private Logger logger = Logger.getLogger(UTagPairToTOPMPlugin.class);

    //Input and output options
    String infile, out_textfile, out_binfile;
    private boolean outputAsText = false, outputAsBinary = false;

    //Internal data for converting tag pairs
    private UTagPairs tp;   //Structure to read in tag pairs
    int startChrom = 1;	//Default dummy chromosome
    int padding = 1000;	//How far to set the dummy coordinates apart
    byte myStrand = 1, myMultimaps = 1, myDcoP = Byte.MIN_VALUE, myMapP = Byte.MIN_VALUE;   //Dummy values to feed to the TOPM structure 
    //int myMaxMapping = 1, myMaxVariants = 8;    //HDF5-specific dummy values for output; not currently in use since DiscoverySNPCallerPlugin doesn't take HDF5-TOPM as input

    public UTagPairToTOPMPlugin() {
        super(null, false);
    }

    public UTagPairToTOPMPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    @Override
    public DataSet performFunction(DataSet input) {
        tp = new UTagPairs(infile);

        TagsOnPhysicalMap myTopm = makeTopmFromTagPairs(tp);

        if (outputAsText) {
            logger.info("Outputting TOPM in text format to " + out_textfile);
            myTopm.writeTextFile(new File(out_textfile));
        }
        if (outputAsBinary) {
            logger.info("Outputting TOPM in binary format to " + out_binfile);
            myTopm.writeBinaryFile(new File(out_binfile));
        }

        return null;
    }

    /**
     * Convert the TagPair object to a TOPM object in memory, filling in things
     * with dummy data where needed
     */
    private TagsOnPhysicalMap makeTopmFromTagPairs(UTagPairs tp) {
        HelperTags tempTags = new HelperTags(tp.getTagNum(), tp.getTag(0).length);
        //Load in all the data for each tag pair; having trouble finding the function to load actual sequences
        for (int i = 0; i < tp.getTagNum(); i++) {
            tempTags.setTag(i, tp.getTag(i), tp.getTagLength(i));
        }
        TagsOnPhysicalMap myTopm = new TagsOnPhysicalMap(tempTags);
        int currPos = 1;
        int currChrom = startChrom;
        for (int i = 0; i < tp.getTagNum(); i++) {

            myTopm.setChromoPosition(i, currChrom, myStrand, currPos, currPos + tp.getTagLength(i) - 1);
            myTopm.setDivergence(i, (byte) 0);

            //These may not be necessary; don't know
            myTopm.setDcoP(i, Byte.MIN_VALUE);  //May have to alter TOPM class to do this
            myTopm.setMapP(i, Byte.MIN_VALUE);
            myTopm.setMultimaps(i, (byte) 1);

            //Increment position after odd-numbered tags (so pairs are at the same position)
            if (i % 2 == 1) {
                currPos += padding;
            }
            //If over max interger value, increment chromosome and start over
            if (currPos >= Integer.MAX_VALUE) {
                currChrom++;
                currPos = 1;
            }
        }
        return myTopm;
    }

    private void printUsage() {
        logger.info(
                "\n\nUsage is as follows:\n"
                + " -c (or --chrom)     Default chromosome number to start at (default: " + startChrom + " )\n"
                + " -d (or --distance)  Distance to pad tag pairs by (default: " + padding + " )\n"
                + " -i (or --input)     TagPair file to convert (required) \n"
                + " -b (or --binary)    Output file name for binary format \n"
                + " -t (or --text)      Output file name for text format\n"
                + "\nNOTE: Must supply at least one of -b or -t\n");
    }

    @Override
    public void setParameters(String[] args) {
        if (args.length == 0) {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease use the above arguments/options.\n\n");
        }
        if (engine == null) {
            engine = new ArgsEngine();
            engine.add("-c", "--chrom", true);
            engine.add("-d", "--distance", true);
            engine.add("-i", "--input", true);
            engine.add("-t", "--text", true);
            engine.add("-b", "--binary", true);
            engine.parse(args);
        }
        //Chromosome
        if (engine.getBoolean("-c")) {
            startChrom = Integer.parseInt(engine.getString("-c"));
        }

        //Padding distance
        if (engine.getBoolean("-d")) {
            padding = Integer.parseInt(engine.getString("-d"));
        }

        //Input file
        if (engine.getBoolean("-i")) {
            infile = engine.getString("-i");
        } else {
            printUsage();
            throw new IllegalArgumentException("\n\nPlease specify an input file.\n\n");
        }

        //Various output files
        if (engine.getBoolean("-b")) {
            out_binfile = engine.getString("-b");
            outputAsBinary = true;
        }
        if (engine.getBoolean("-t")) {
            out_textfile = engine.getString("-t");
            outputAsText = true;
        }

        //Check that at least one output format selected
        if (!(outputAsBinary || outputAsText)) {
            throw new IllegalArgumentException("\n\nPlease specify at least one output format .\n\n");
        }
    }

    @Override
    public ImageIcon getIcon() {
        throw new UnsupportedOperationException("Not supported yet.");
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

/*
 * A small helper class that exists solely to pass tag info to TOPM
 */
class HelperTags extends AbstractTags {

    /* Inherited
     protected int tagLengthInLong;  //
     protected long[][] tags;  //Index in rows, actual tag components in columns
     protected byte[] tagLength;  // length of tag (number of bases)  // 1 byte
     */
    public HelperTags(int numtags, int myTagLengthInLong) {
        tagLengthInLong = myTagLengthInLong;
        tags = new long[tagLengthInLong][numtags];
        tagLength = new byte[numtags];
    }

    public void setTag(int index, long[] tagValue, byte myTagLength) {
        tagLength[index] = myTagLength;
        for (int i = 0; i < tagValue.length; i++) {
            tags[i][index] = tagValue[i];
        }
    }
}
