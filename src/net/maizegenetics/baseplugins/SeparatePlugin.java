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
                fireDataSetReturned(new PluginEvent(tds, SeparatePlugin.class));

            } else if (currentValue instanceof Alignment) {
                Alignment align = (Alignment) current.getData();
                Alignment[] alignments = align.getAlignments();
                if (alignments.length > 1) {
                    for (int i = 0; i < alignments.length; i++) {
                        if (alignments[i] != null) {
                            String name = current.getName() + "_align" + i;
                            Datum td = new Datum(name, alignments[i], null);
                            DataSet tds = new DataSet(td, null);
                            result.add(tds);
                            fireDataSetReturned(new PluginEvent(tds, SeparatePlugin.class));
                        }
                    }
                }
            }
        }

        if (result.size() == 0) {
            if (isInteractive()) {
                JOptionPane.showMessageDialog(getParentFrame(), "Nothing to Separate");
            } else {
                myLogger.warn("performFunction: Nothing to Separate.");
            }
            return null;
        } else {
            return DataSet.getDataSet(result, this);
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
