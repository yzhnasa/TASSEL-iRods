/*
 * GenotypeSummaryPlugin
 */
package net.maizegenetics.baseplugins;

import java.util.Arrays;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

import java.awt.Frame;

import java.net.URL;

import javax.swing.ImageIcon;
import javax.swing.JOptionPane;

import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.report.SimpleTableReport;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginEvent;

import org.apache.log4j.Logger;

/**
 *
 * @author terry
 */
public class GenotypeSummaryPlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(GenotypeSummaryPlugin.class);
    private static final String NA = "NA";
    private static final Double ZERO_DOUBLE = 0.0;
    private static final int ZERO_INT = 0;
    private static final Integer ONE_INTEGER = 1;
    private int myNumGametesMissing = 0;
    private int myNumHeterozygous = 0;

    public GenotypeSummaryPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    public DataSet performFunction(DataSet input) {

        try {

            myNumGametesMissing = 0;
            myNumHeterozygous = 0;

            List<Datum> alignInList = input.getDataOfType(Alignment.class);

            if (alignInList.size() != 1) {
                String gpMessage = "Invalid selection.  Please select one genotype alignment.";
                if (isInteractive()) {
                    JOptionPane.showMessageDialog(getParentFrame(), gpMessage);
                } else {
                    myLogger.error(gpMessage);
                }
                return null;
            }

            Datum current = alignInList.get(0);
            Alignment alignment = (Alignment) current.getData();
            String name = current.getName();

            SimpleTableReport siteSummary = getSiteSummary(alignment);
            //SimpleTableReport taxaSummary = getTaxaSummary(alignment);
            SimpleTableReport overallSummary = getOverallSummary(alignment);

            List<Datum> summaryTables = new ArrayList<Datum>();
            summaryTables.add(new Datum(name + "_OverallSummary", overallSummary, "Overall Summary of " + name));
            summaryTables.add(new Datum(name + "_SiteSummary", siteSummary, "Site Summary of " + name));
            //summaryTables.add(new Datum(name + "_TaxaSummary", taxaSummary, "Taxa Summary of " + name));

            if (summaryTables.isEmpty()) {
                return null;
            }

            DataSet output = new DataSet(summaryTables, this);
            fireDataSetReturned(new PluginEvent(output, GenotypeSummaryPlugin.class));

            return output;

        } finally {
            fireProgress(100);
        }

    }

    private SimpleTableReport getOverallSummary(Alignment alignment) {

        Object[] firstColumnNames = new String[]{"Stat Type", "Value", "Diploid Value", "Number", "Proportion", "Frequency"};

        int numSites = alignment.getSiteCount();
        int numTaxa = alignment.getSequenceCount();

        HashMap<String, Integer> diploidValueCounts = new HashMap<String, Integer>();
        for (int r = 0; r < numTaxa; r++) {
            for (int c = 0; c < numSites; c++) {
                String current = alignment.getBaseAsString(r, c);
                Integer num = diploidValueCounts.get(current);
                if (num == null) {
                    diploidValueCounts.put(current, ONE_INTEGER);
                } else {
                    diploidValueCounts.put(current, ++num);
                }
            }
        }

        int totalGametes = numSites * numTaxa * 2;
        int totalGametesNotMissing = totalGametes - myNumGametesMissing;

        int numDiploidsMissing = 0;
        Integer numMissingInt = diploidValueCounts.get(Alignment.UNKNOWN_ALLELE_STR);
        if (numMissingInt == null) {
            numMissingInt = diploidValueCounts.get(Alignment.UNKNOWN_DIPLOID_ALLELE_STR);
            if (numMissingInt != null) {
                numDiploidsMissing = numMissingInt.intValue();
            }
        } else {
            numDiploidsMissing = numMissingInt.intValue();
        }

        int totalDiploids = numSites * numTaxa;
        int totalDiploidsNotMissing = totalDiploids - numDiploidsMissing;

        int count = 0;
        int numRows = Math.max(diploidValueCounts.size(), 13);
        Object[][] data = new Object[numRows][firstColumnNames.length];

        data[count][0] = "Number of Taxa";
        data[count++][1] = (double) numTaxa;

        data[count][0] = "Number of Sites";
        data[count++][1] = (double) numSites;

        data[count][0] = "Sites x Taxa";
        data[count++][1] = (double) totalDiploids;

        data[count][0] = "Number Not Missing";
        data[count++][1] = (double) totalDiploidsNotMissing;

        data[count][0] = "Proportion Not Missing";
        data[count++][1] = (double) totalDiploidsNotMissing / (double) totalDiploids;

        data[count][0] = "Number Missing";
        data[count++][1] = (double) numDiploidsMissing;

        data[count][0] = "Proportion Missing";
        data[count++][1] = (double) numDiploidsMissing / (double) totalDiploids;

        data[count][0] = "Number Gametes";
        data[count++][1] = (double) totalGametes;

        data[count][0] = "Gametes Not Missing";
        data[count++][1] = (double) totalGametesNotMissing;

        data[count][0] = "Proportion Gametes Not Missing";
        data[count++][1] = (double) totalGametesNotMissing / (double) totalGametes;

        data[count][0] = "Gametes Missing";
        data[count++][1] = (double) myNumGametesMissing;

        data[count][0] = "Proportion Gametes Missing";
        data[count++][1] = (double) myNumGametesMissing / (double) totalGametes;

        data[count][0] = "Proportion Heterozygous";
        data[count++][1] = (double) myNumHeterozygous / (double) totalGametes;

        count = 0;
        Iterator itr = diploidValueCounts.keySet().iterator();
        while (itr.hasNext()) {
            String value = (String) itr.next();
            Integer numValue = (Integer) diploidValueCounts.get(value);
            data[count][2] = value;
            data[count][3] = numValue.intValue();
            data[count][4] = numValue.doubleValue() / (double) totalDiploids;
            data[count++][5] = numValue.doubleValue() / (double) totalDiploidsNotMissing;
        }

        return new SimpleTableReport("Overall Summary", firstColumnNames, data);
    }

    private SimpleTableReport getSiteSummary(Alignment alignment) {

        String[] firstColumnNames = new String[]{"Site Number", "Site Name", "Physical Position", "Number of Taxa", "Major Allele", "Major Allele Gametes", "Major Allele Proportion", "Major Allele Frequency",
            "Minor Allele", "Minor Allele Gametes", "Minor Allele Proportion", "Minor Allele Frequency"};
        String[] lastColumnNames = new String[]{"Gametes Missing", "Proportion Missing", "Proportion Heterozygous",
            "Inbreeding Coefficient", "Inbreeding Coefficient Scaled by Missing"};

        List columnNames = new ArrayList(Arrays.asList(firstColumnNames));

        int maxAlleles = alignment.getMaxNumAlleles();
        if (alignment.retainsRareAlleles()) {
            maxAlleles++;
        }
        for (int i = 2; i < maxAlleles; i++) {
            String alleleHeading = "Allele " + (i + 1);
            columnNames.add(alleleHeading);
            columnNames.add(alleleHeading + " Gametes");
            columnNames.add(alleleHeading + " Proportion");
            columnNames.add(alleleHeading + " Frequency");
        }

        columnNames.addAll(Arrays.asList(lastColumnNames));

        int numSites = alignment.getSiteCount();
        Object[][] data = new Object[numSites][columnNames.size()];

        int totalGametes = alignment.getSequenceCount() * 2;

        for (int i = 0; i < numSites; i++) {

            int totalGametesNotMissing = alignment.getTotalGametesNotMissing(i);
            int count = 0;

            data[i][count++] = i;
            data[i][count++] = alignment.getSNPID(i);
            data[i][count++] = alignment.getPositionInLocus(i);
            data[i][count++] = alignment.getSequenceCount();

            int[][] alleles = alignment.getAllelesSortedByFrequency(i);
            int numAlleles = alleles[0].length;

            for (int a = 0; a < numAlleles; a++) {
                data[i][count++] = alignment.getBaseAsString(i, (byte) alleles[0][a]);
                data[i][count++] = alleles[1][a];
                data[i][count++] = (double) alleles[1][a] / (double) totalGametes;
                data[i][count++] = (double) alleles[1][a] / (double) totalGametesNotMissing;
            }

            for (int b = 0; b < (maxAlleles - numAlleles); b++) {
                data[i][count++] = NA;
                data[i][count++] = ZERO_INT;
                data[i][count++] = ZERO_DOUBLE;
                data[i][count++] = ZERO_DOUBLE;
            }

            int totalGametesMissing = totalGametes - totalGametesNotMissing;
            myNumGametesMissing = myNumGametesMissing + totalGametesMissing;
            data[i][count++] = totalGametesMissing;
            data[i][count++] = (double) totalGametesMissing / (double) totalGametes;

            int numHeterozygous = alignment.getHeterozygousCount(i);
            myNumHeterozygous = myNumHeterozygous + numHeterozygous;
            data[i][count++] = (double) numHeterozygous / (double) totalGametes;

            data[i][count++] = "TBD";
            data[i][count++] = "TBD";

        }

        String[] columnNameStrings = new String[columnNames.size()];
        columnNames.toArray(columnNameStrings);
        return new SimpleTableReport("Site Summary", columnNameStrings, data);

    }

    private SimpleTableReport getTaxaSummary(Alignment alignment) {

        Object[] columnNames = new String[]{"Taxa", "Number of Sites", "Proportion Missing",
            "Proportion Heterozygous", "Inbreeding Coefficient",
            "Inbreeding Coefficient Scaled by Missing"};
        int numSites = alignment.getSiteCount();
        int numTaxa = alignment.getSequenceCount();
        Object[][] data = new Object[numTaxa][columnNames.length];

        int totalGametes = numSites * 2;
        for (int i = 0; i < numTaxa; i++) {
            int count = 0;
            data[i][count++] = i;
            data[i][count++] = numSites;
            data[i][count++] = "Proportion Missing";
            data[i][count++] = "Proportion Heterozygous";
            data[i][count++] = "Inbreeding Coefficient";
            data[i][count++] = "ICSBM";
        }

        return new SimpleTableReport("Taxa Summary", columnNames, data);

    }

    public ImageIcon getIcon() {
        URL imageURL = GenotypeSummaryPlugin.class.getResource("images/summary.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    public String getButtonName() {
        return "Geno Summary";
    }

    public String getToolTipText() {
        return "Genotype Summary";
    }
}
