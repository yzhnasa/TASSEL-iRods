/*
 * BinaryToText
 */
package net.maizegenetics.analysis.gbs;

import com.google.common.collect.Range;

import net.maizegenetics.dna.map.TagsOnPhysicalMap;
import net.maizegenetics.dna.tag.TagCounts;
import net.maizegenetics.dna.tag.TagsByTaxa.FilePacking;
import net.maizegenetics.dna.tag.TagsByTaxaByte;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.PluginParameterTerry;

import org.apache.log4j.Logger;

import javax.swing.*;
import java.awt.*;
import java.io.File;
import java.util.Arrays;

/**
 *
 * @author Terry Casstevens
 */
public class BinaryToTextPluginNew extends ParameterConceptPlugin {

    private Logger myLogger = Logger.getLogger(BinaryToTextPluginNew.class);

    public static enum FILE_TYPES {

        TOPM, TagCounts, TBTByte
    };

    public enum PARAMETERS {

        inputFile, outputFile, fileType
    };

    protected PluginParameterTerry<String> myInputFile = new PluginParameterTerry.Builder<String>(PARAMETERS.inputFile, null, String.class).required(true).inFile().build();
    protected PluginParameterTerry<String> myOutputFile = new PluginParameterTerry.Builder<String>(PARAMETERS.outputFile, null, String.class).required(true).outFile().build();
    protected PluginParameterTerry<FILE_TYPES> myFileType = new PluginParameterTerry.Builder<FILE_TYPES>(PARAMETERS.fileType, FILE_TYPES.TOPM, FILE_TYPES.class).range(Range.encloseAll(Arrays.asList(FILE_TYPES.values()))).build();

    public BinaryToTextPluginNew(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    @Override
    public DataSet processData(DataSet input) {

        switch (getType()) {
            case TOPM:
                TagsOnPhysicalMap topm = new TagsOnPhysicalMap(getInput(), true);
                topm.writeTextFile(new File(getOutput()));
                break;
            case TagCounts:
                TagCounts tc = new TagCounts(getInput(), FilePacking.Byte);
                tc.writeTagCountFile(getOutput(), FilePacking.Text, 0);
                break;
            case TBTByte:
                TagsByTaxaByte tbtbyte = new TagsByTaxaByte(getInput(), FilePacking.Byte);
                tbtbyte.writeDistFile(new File(getOutput()), FilePacking.Text, 0);
                break;
        }

        return null;

    }

    public void setInput(String filename) {
        setParameter(PARAMETERS.inputFile, filename);
    }

    public String getInput() {
        return myInputFile.value();
    }

    public void setOutput(String filename) {
        setParameter(PARAMETERS.outputFile, filename);
    }

    public String getOutput() {
        return myOutputFile.value();
    }

    public void setType(FILE_TYPES type) {
        setParameter(PARAMETERS.fileType, type.toString());
    }

    public FILE_TYPES getType() {
        return myFileType.value();
    }

    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "Binary to Text";
    }

    @Override
    public String getToolTipText() {
        return "Binary to Text";
    }
}
