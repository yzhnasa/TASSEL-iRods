/*
 * SeparatePlugin.java
 *
 * Created on July 31, 2010
 *
 */
package net.maizegenetics.baseplugins;

import net.maizegenetics.pal.alignment.*;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginEvent;

import java.awt.Frame;

import javax.swing.*;

import java.net.URL;

import java.util.ArrayList;
import java.util.List;

import org.apache.log4j.Logger;

/**
 *
 * @author terry
 */
public class SeparatePlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(SeparatePlugin.class);

    /** Creates a new instance of SeparatePlugin */
    public SeparatePlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    public DataSet performFunction(DataSet input) {

        try {
            List<Datum> inputs = input.getDataSet();
            List<DataSet> result = new ArrayList();

            for (Datum current : inputs) {
                Object currentValue = current.getData();
                if (currentValue instanceof MarkerPhenotype) {

                    MarkerPhenotype mp = (MarkerPhenotype) currentValue;
                    Phenotype pheno = mp.getPhenotype();
                    String phenoName = current.getName() + "_pheno";
                    Datum phenoDatum = new Datum(phenoName, pheno, null);

                    Alignment align = mp.getAlignment();
                    String alignName = current.getName() + "_align";
                    Datum alignDatum = new Datum(alignName, align, null);

                    DataSet tds = new DataSet(new Datum[]{phenoDatum, alignDatum}, this);
                    result.add(tds);

                } else if (currentValue instanceof Alignment) {

                    Alignment align = (Alignment) currentValue;
                    Alignment[] alignments = align.getAlignments();
                    for (int i = 0; i < alignments.length; i++) {
                        int[] offsets = alignments[i].getLociOffsets();
                        if (offsets.length > 1) {
                            Locus[] loci = alignments[i].getLoci();
                            for (int j = 0; j < offsets.length; j++) {
                                String name = current.getName() + "_chrom" + loci[j];
                                int endSite;
                                try {
                                    endSite = offsets[j + 1] - 1;
                                } catch (Exception e) {
                                    endSite = alignments[i].getSiteCount() - 1;
                                }
                                Datum td = new Datum(name, FilterAlignment.getInstance(alignments[i], offsets[j], endSite), null);
                                DataSet tds = new DataSet(td, null);
                                result.add(tds);
                            }
                        } else {
                            if (alignments.length > 1) {
                                String name = current.getName() + "_chrom" + alignments[i].getLocus(0);
                                Datum td = new Datum(name, alignments[i], null);
                                DataSet tds = new DataSet(td, null);
                                result.add(tds);
                            }
                        }
                    }
                }
            }

            if (result.isEmpty()) {
                if (isInteractive()) {
                    JOptionPane.showMessageDialog(getParentFrame(), "Nothing to Separate");
                } else {
                    myLogger.warn("performFunction: Nothing to Separate.");
                }
                return null;
            } else {
                DataSet resultDataSet = DataSet.getDataSet(result, this);
                fireDataSetReturned(new PluginEvent(resultDataSet, SeparatePlugin.class));
                return resultDataSet;
            }
        } finally {
            fireProgress(100);
        }

    }

    /**
     * Icon for this plugin to be used in buttons, etc.
     *
     * @return ImageIcon
     */
    public ImageIcon getIcon() {
        URL imageURL = SeparatePlugin.class.getResource("images/Separate.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    /**
     * Button name for this plugin to be used in buttons, etc.
     *
     * @return String
     */
    public String getButtonName() {
        return "Separate";
    }

    /**
     * Tool Tip Text for this plugin
     *
     * @return String
     */
    public String getToolTipText() {
        return "Separate Data into Components";
    }
}
