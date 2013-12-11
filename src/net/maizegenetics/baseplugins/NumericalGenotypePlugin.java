/*
 * NumericalGenotypePlugin.java
 *
 * Created on May 8, 2008
 *
 */
package net.maizegenetics.baseplugins;

import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableUtils;
import net.maizegenetics.trait.SimplePhenotype;
import net.maizegenetics.trait.Trait;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;

import javax.swing.*;
import java.util.*;

/**
 * Numerical Genotype Plugin.
 *
 * @author terryc
 */
public class NumericalGenotypePlugin extends AbstractPlugin {

    public static enum TRANSFORM_TYPE {

        collapse, separated
    };
    private TRANSFORM_TYPE myTransformType = TRANSFORM_TYPE.collapse;

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

        List data = dataSet.getDataOfType(GenotypeTable.class);
        if ((data != null) && (data.size() == 1)) {
            Datum datum = (Datum) data.get(0);
            if (myTransformType == TRANSFORM_TYPE.collapse) {
                GenotypeTable input = (GenotypeTable) datum.getData();
                SimplePhenotype result = collapseTransform(input);
                DataSet tds = new DataSet(new Datum(datum.getName() + "_Collapse", result, null), this);
                fireDataSetReturned(tds);
                return tds;
            } else if (myTransformType == TRANSFORM_TYPE.separated) {
                GenotypeTable input = (GenotypeTable) datum.getData();
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

    public static SimplePhenotype separatedTransform(GenotypeTable input) {

        int seqCount = input.numberOfTaxa();
        int siteCount = input.numberOfSites();

        Map[] alleleCounts = new Map[siteCount];

        for (int i = 0; i < siteCount; i++) {

            Map<Byte, Integer> charCount = new HashMap<Byte, Integer>();
            for (int j = 0; j < seqCount; j++) {

                byte[] current = input.genotypeArray(j, i);
                byte orderedValue = GenotypeTableUtils.getUnphasedDiploidValue(current[0], current[1]);
                if (orderedValue != GenotypeTable.UNKNOWN_DIPLOID_ALLELE) {
                    Integer count = (Integer) charCount.get(orderedValue);
                    if (count != null) {
                        charCount.put(orderedValue, new Integer(count.intValue() + 1));
                    } else {
                        charCount.put(orderedValue, 1);
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
        List<Trait> traitNames = new ArrayList<Trait>();

        int offset = 0;
        for (int i = 0; i < siteCount; i++) {

            Map<Byte, Integer> currentHash = alleleCounts[i];
            int currentSize = currentHash.size();
            byte[] sortChars = new byte[currentSize];
            Iterator itr = currentHash.keySet().iterator();
            int count = 0;
            while (itr.hasNext()) {
                sortChars[count++] = (Byte) itr.next();
            }

            boolean change = true;
            while (change) {

                change = false;

                for (int k = 0; k < currentSize - 1; k++) {
                    Integer first = (Integer) currentHash.get(sortChars[k]);
                    Integer second = (Integer) currentHash.get(sortChars[k + 1]);
                    if (first.compareTo(second) > 0) {
                        byte temp = sortChars[k];
                        sortChars[k] = sortChars[k + 1];
                        sortChars[k + 1] = temp;
                        change = true;
                    }
                }

            }

            for (int k = 0; k < currentSize; k++) {
                traitNames.add(new Trait("S" + i, false, Trait.TYPE_DATA));
                currentHash.put(sortChars[k], k);
            }

            for (int j = 0; j < seqCount; j++) {
                byte[] current = input.genotypeArray(j, i);
                byte orderedValue = GenotypeTableUtils.getUnphasedDiploidValue(current[0], current[1]);
                if (orderedValue == GenotypeTable.UNKNOWN_DIPLOID_ALLELE) {
                    for (int k = 0; k < currentSize; k++) {
                        pcValues[j][k + offset] = Double.NaN;
                    }
                } else {
                    int position = ((Integer) currentHash.get(orderedValue)).intValue();
                    pcValues[j][position + offset] = 1.0;
                }
            }

            offset = offset + currentSize;

        }

        return new SimplePhenotype(input.taxa(), traitNames, pcValues);
    }

    public static SimplePhenotype collapseTransform(GenotypeTable input) {

        int seqCount = input.numberOfTaxa();
        int siteCount = input.numberOfSites();

        double[][] pcValues = new double[seqCount][siteCount];
        List<Trait> traitNames = new ArrayList<Trait>();

        for (int i = 0; i < siteCount; i++) {

            traitNames.add(new Trait("S" + i, false, Trait.TYPE_DATA));

            byte[] allelesByFreq = input.alleles(i);

            for (int j = 0; j < seqCount; j++) {
                byte[] current = input.genotypeArray(j, i);
                if ((current[0] == GenotypeTable.UNKNOWN_ALLELE) && (current[1] == GenotypeTable.UNKNOWN_ALLELE)) {
                    pcValues[j][i] = Double.NaN;
                } else {
                    pcValues[j][i] = 1.0;
                    if (current[0] == allelesByFreq[0]) {
                        pcValues[j][i] -= 0.5;
                    }
                    if (current[1] == allelesByFreq[0]) {
                        pcValues[j][i] -= 0.5;
                    }
                }
            }

        }

        return new SimplePhenotype(input.taxa(), traitNames, pcValues);

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
