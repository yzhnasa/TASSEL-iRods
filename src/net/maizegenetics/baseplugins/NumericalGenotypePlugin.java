/*
 * NumericalGenotypePlugin.java
 *
 * Created on May 8, 2008
 *
 */
package net.maizegenetics.baseplugins;

import java.util.ArrayList;
import java.util.List;

import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.SimplePhenotype;
import net.maizegenetics.pal.alignment.Trait;

import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;

import javax.swing.*;

/**
 * Numerical Genotype Plugin.
 *
 * @author terryc
 */
public class NumericalGenotypePlugin extends AbstractPlugin {

    public static enum TRANSFORM_TYPE {

        colapse, separated
    };
    private TRANSFORM_TYPE myTransformType = TRANSFORM_TYPE.colapse;

    /**
     * Principle Components Transform Plugin.
     */
    public NumericalGenotypePlugin() {
        super(null, false);
    }

    /**
     *
     */
    public DataSet performFunction(DataSet dataSet) {

        List data = dataSet.getDataOfType(Alignment.class);
        if ((data != null) && (data.size() == 1)) {
            Datum datum = (Datum) data.get(0);
            if (myTransformType == TRANSFORM_TYPE.colapse) {
                Alignment input = (Alignment) datum.getData();
                SimplePhenotype result = collapseTransform(input);
                DataSet tds = new DataSet(new Datum(datum.getName() + "_Collapse", result, null), this);
                fireDataSetReturned(tds);
                return tds;
            } else if (myTransformType == TRANSFORM_TYPE.separated) {
                Alignment input = (Alignment) datum.getData();
                SimplePhenotype result = separatedTransform(input);
                DataSet tds = new DataSet(new Datum(datum.getName() + "_Separated", result, null), this);
                fireDataSetReturned(tds);
                return tds;
            } else {
                throw new IllegalArgumentException("NumericalGenotypePlugin: performFunction: unknown transform type.");
            }
        }

        return null;

    }

    public static SimplePhenotype separatedTransform(Alignment input) {

        int seqCount = input.getSequenceCount();
        int siteCount = input.getSiteCount();

        byte[][] alleleCounts = new byte[siteCount][];
        for (int i = 0; i < siteCount; i++) {
            alleleCounts[i] = input.getAlleles(i);
        }

        int columnCount = 0;
        for (int i = 0; i < siteCount; i++) {
            columnCount = columnCount + alleleCounts[i].length;
        }

        double[][] pcValues = new double[seqCount][columnCount];
        List<Trait> traitNames = new ArrayList<Trait>();

        int offset = 0;
        for (int i = 0; i < siteCount; i++) {

            int currentSize = alleleCounts[i].length;
            for (int k = 0; k < currentSize; k++) {
                traitNames.add(new Trait("S" + i, false, Trait.TYPE_DATA));
            }

            for (int j = 0; j < seqCount; j++) {
                byte[] current = input.getBaseArray(j, i);
                if ((current[0] == Alignment.UNKNOWN_ALLELE) && (current[1] == Alignment.UNKNOWN_ALLELE)) {
                    for (int k = 0; k < currentSize; k++) {
                        pcValues[j][k + offset] = Double.NaN;
                    }
                } else {
                    for (int k = 0; k < currentSize; k++) {
                        if (current[0] == alleleCounts[i][k]) {
                            pcValues[j][currentSize - k - 1 + offset] = 1.0;
                        }
                        if (current[1] == alleleCounts[i][k]) {
                            pcValues[j][currentSize - k - 1 + offset] = 1.0;
                        }
                    }
                }

            }

            offset = offset + currentSize;

        }

        return new SimplePhenotype(input.getIdGroup(), traitNames, pcValues);

    }

    public static SimplePhenotype collapseTransform(Alignment input) {

        int seqCount = input.getSequenceCount();
        int siteCount = input.getSiteCount();

        double[][] pcValues = new double[seqCount][siteCount];
        List<Trait> traitNames = new ArrayList<Trait>();

        for (int i = 0; i < siteCount; i++) {

            traitNames.add(new Trait("S" + i, false, Trait.TYPE_DATA));

            byte[] allelesByFreq = input.getAlleles(i);

            for (int j = 0; j < seqCount; j++) {
                byte[] current = input.getBaseArray(j, i);
                if ((current[0] == allelesByFreq[0]) || (current[1] == allelesByFreq[0])) {
                    pcValues[j][i] = 0.0;
                } else if ((current[0] == Alignment.UNKNOWN_ALLELE) && (current[1] == Alignment.UNKNOWN_ALLELE)) {
                    pcValues[j][i] = Double.NaN;
                } else {
                    pcValues[j][i] = 1.0;
                }
            }

        }

        return new SimplePhenotype(input.getIdGroup(), traitNames, pcValues);

    }

    public void setTransformType(TRANSFORM_TYPE type) {
        myTransformType = type;
    }

    public TRANSFORM_TYPE getTransformType() {
        return myTransformType;
    }

    public String getToolTipText() {
        return "";
    }

    public ImageIcon getIcon() {
        return null;
    }

    public String getButtonName() {
        return "";
    }
}
