/*
 * NumericalGenotypePlugin.java
 *
 * Created on May 8, 2008
 *
 */
package net.maizegenetics.baseplugins;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

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

        Map[] alleleCounts = new Map[siteCount];

        for (int i = 0; i < siteCount; i++) {

            Map charCount = new HashMap();
            for (int j = 0; j < seqCount; j++) {

                char current = input.getBaseChar(j, i);
                if ((current != 'N') && (current != '?')) {
                    Integer count = (Integer) charCount.get(current);
                    if (count != null) {
                        charCount.put(current, new Integer(count.intValue() + 1));
                    } else {
                        charCount.put(current, 1);
                    }

                }

            }

            alleleCounts[i] = charCount;

        }

        int columnCount = 0;
        for (int i = 0; i < siteCount; i++) {
            columnCount = columnCount + alleleCounts[i].size();
        }

        double[][] pcValues = new double[seqCount][columnCount];
        //String[][] traitNames = new String[2][columnCount];
        List<Trait> traitNames = new ArrayList<Trait>();

        int offset = 0;
        for (int i = 0; i < siteCount; i++) {

            Map currentHash = alleleCounts[i];
            int currentSize = currentHash.size();
            Character[] sortChars = new Character[currentSize];
            Iterator itr = currentHash.keySet().iterator();
            int count = 0;
            while (itr.hasNext()) {
                sortChars[count++] = (Character) itr.next();
            }

            boolean change = true;
            while (change) {

                change = false;

                for (int k = 0; k < currentSize - 1; k++) {
                    Integer first = (Integer) currentHash.get(sortChars[k]);
                    Integer second = (Integer) currentHash.get(sortChars[k + 1]);
                    if (first.compareTo(second) > 0) {
                        Character temp = sortChars[k];
                        sortChars[k] = sortChars[k + 1];
                        sortChars[k + 1] = temp;
                        change = true;
                    }
                }

            }

            for (int k = 0; k < currentSize; k++) {
                //traitNames[0][k + offset] = "S" + i;
                //traitNames[1][k + offset] = sortChars[k].toString();
                traitNames.add(new Trait("S" + i, false, Trait.TYPE_DATA));
                currentHash.put(sortChars[k], k);
            }



            for (int j = 0; j < seqCount; j++) {
                char current = input.getBaseChar(j, i);
                if ((current == 'N') || (current == '?')) {
                    for (int k = 0; k < currentSize; k++) {
                        pcValues[j][k + offset] = Double.NaN;
                    }
                } else {
                    int position = ((Integer) currentHash.get(current)).intValue();
                    pcValues[j][position + offset] = 1.0;
                }
            }

            offset = offset + currentSize;

        }

        //SimplePhenotype result = new SimplePhenotype(input, pcValues, new String[]{"Trait", "Env"}, traitNames);
        SimplePhenotype result = new SimplePhenotype(input.getIdGroup(), traitNames, pcValues);
        return result;

    }

    public static SimplePhenotype collapseTransform(Alignment input) {

        int seqCount = input.getSequenceCount();
        int siteCount = input.getSiteCount();

        double[][] pcValues = new double[seqCount][siteCount];
        //String[] traitNames = new String[siteCount];
        List<Trait> traitNames = new ArrayList<Trait>();

        for (int i = 0; i < siteCount; i++) {

            //traitNames[i] = "S" + i;
            traitNames.add(new Trait("S" + i, false, Trait.TYPE_DATA));

            Map charCount = new HashMap();
            char highestChar = 'z';
            for (int j = 0; j < seqCount; j++) {

                char current = input.getBaseChar(j, i);
                if ((current != 'N') && (current != '?')) {
                    Integer count = (Integer) charCount.get(current);
                    if (count != null) {
                        charCount.put(current, new Integer(count.intValue() + 1));
                    } else {
                        charCount.put(current, 1);
                    }

                }

            }

            Iterator itr = charCount.keySet().iterator();
            Integer highestCount = -1;
            while (itr.hasNext()) {
                char currentChar = ((Character) itr.next()).charValue();
                Integer currentCount = (Integer) charCount.get(currentChar);
                if (currentCount.compareTo(highestCount) > 0) {
                    highestCount = currentCount;
                    highestChar = currentChar;
                }

            }


            for (int j = 0; j < seqCount; j++) {
                char current = input.getBaseChar(j, i);
                if ((current == 'N') || (current == '?')) {
                    pcValues[j][i] = Double.NaN;
                } else if (current == highestChar) {
                    pcValues[j][i] = 0.0;
                } else {
                    pcValues[j][i] = 1.0;
                }
            }

        }

        //return new SimplePhenotype(input, pcValues, traitNames);
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
