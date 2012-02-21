/*
 * GenotypeSummaryPlugin
 */
package net.maizegenetics.baseplugins;

import java.util.ArrayList;
import java.util.List;

import java.awt.Frame;

import java.util.Arrays;
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
    private static final Double ZERO_DOUBLE = new Double(0.0);
    private int myNumGametesMissing = 0;

    public GenotypeSummaryPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    public DataSet performFunction(DataSet input) {

        try {

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
            SimpleTableReport taxaSummary = getTaxaSummary(alignment);
            SimpleTableReport overallSummary = getOverallSummary(alignment);

            List<Datum> summaryTables = new ArrayList<Datum>();
            summaryTables.add(new Datum(name + "_OverallSummary", overallSummary, "Overall Summary of " + name));
            summaryTables.add(new Datum(name + "_SiteSummary", siteSummary, "Site Summary of " + name));
            summaryTables.add(new Datum(name + "_TaxaSummary", taxaSummary, "Taxa Summary of " + name));

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

        Object[] columnNames = new String[]{"Number of Taxa", "Number of Sites", "Proportion Missing",
            "Proportion Heterozygous"};

        int numSites = alignment.getSiteCount();
        int totalGametes = numSites * alignment.getSequenceCount() * 2;
        int count = 0;
        Object[][] data = new Object[1][columnNames.length];
        data[0][count++] = alignment.getSequenceCount();
        data[0][count++] = numSites;
        data[0][count++] = (double) myNumGametesMissing / (double) totalGametes;
        data[0][count++] = "TBD";

        return new SimpleTableReport("Overall Summary", columnNames, data);
    }

    private SimpleTableReport getSiteSummary(Alignment alignment) {

        String[] firstColumnNames = new String[]{"Site Number", "Number of Taxa", "Major Allele", "Major Proportion", "Major Allele Frequency",
            "Minor Allele", "Minor Proportion", "Minor Allele Frequency"};
        String[] lastColumnNames = new String[]{"Proportion Missing", "Proportion Heterozygous",
            "Inbreeding Coefficient", "Inbreeding Coefficient Scaled by Missing"};

        List columnNames = new ArrayList(Arrays.asList(firstColumnNames));

        int maxAlleles = alignment.getMaxNumAlleles();
        for (int i = 2; i < maxAlleles; i++) {
            String alleleHeading = "Allele " + (i + 1);
            columnNames.add(alleleHeading);
            columnNames.add(alleleHeading + "Proportion");
            columnNames.add(alleleHeading + "Frequency");
        }
        
        if (alignment.retainsRareAlleles()) {
            maxAlleles++;
            columnNames.add("Rare Alleles");
            columnNames.add("Rare Proportion");
            columnNames.add("Rare Frequency");
        }

        columnNames.addAll(Arrays.asList(lastColumnNames));

        int numSites = alignment.getSiteCount();
        Object[][] data = new Object[numSites][columnNames.size()];

        int totalGametes = alignment.getSequenceCount() * 2;

        for (int i = 0; i < numSites; i++) {

            int totalGametesNotMissing = alignment.getTotalGametesNotMissing(i);
            int count = 0;

            data[i][count++] = i;
            data[i][count++] = alignment.getSequenceCount();

            int[][] alleles = alignment.getAllelesSortedByFrequency(i);
            int numAlleles = alleles[0].length;
            
            System.out.println("site" + i + "  max alleles: " + maxAlleles + "  num alleles: " + numAlleles);

            for (int a = 0; a < numAlleles; a++) {
                data[i][count++] = alignment.getBaseAsString(i, (byte) alleles[0][a]);
                data[i][count++] = (double) alleles[1][a] / (double) totalGametes;
                data[i][count++] = (double) alleles[1][a] / (double) totalGametesNotMissing;
            }

            for (int b = 0; b < (maxAlleles - numAlleles); b++) {
                data[i][count++] = NA;
                data[i][count++] = ZERO_DOUBLE;
                data[i][count++] = ZERO_DOUBLE;
            }

            data[i][count++] = (double) (totalGametes - totalGametesNotMissing) / (double) totalGametes;
            data[i][count++] = "TBD";
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
        return null;
    }

    public String getButtonName() {
        return "Genotype Summary";
    }

    public String getToolTipText() {
        return "Genotype Summary";
    }
}
