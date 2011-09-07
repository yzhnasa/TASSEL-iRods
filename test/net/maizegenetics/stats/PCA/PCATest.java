package net.maizegenetics.stats.PCA;


import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;
import java.util.regex.Pattern;

import javax.swing.JTable;
import javax.swing.ListSelectionModel;

import junit.framework.Assert;

import net.maizegenetics.baseplugins.NumericalGenotypePlugin;
import net.maizegenetics.baseplugins.numericaltransform.Conversion;
import net.maizegenetics.baseplugins.numericaltransform.ImputePanel;
import net.maizegenetics.baseplugins.numericaltransform.NumericalTransformPlugin;
import net.maizegenetics.baseplugins.numericaltransform.PCAPanel;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.AnnotatedAlignmentUtils;
import net.maizegenetics.pal.alignment.ImportUtils;
import net.maizegenetics.pal.alignment.Phenotype;
import net.maizegenetics.pal.alignment.SimplePhenotype;
import net.maizegenetics.pal.alignment.Trait;
import net.maizegenetics.pal.report.TableReport;
import net.maizegenetics.plugindef.Datum;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;

public class PCATest {
	private Alignment geno;
	private String genotypeFile = "test/datafiles/mdp_genotype.hmp.txt";
	Datum imputedDatum;
	
	@Before
	public void setUp() throws Exception {
		geno = ImportUtils.readFromHapmap(genotypeFile);
		geno = AnnotatedAlignmentUtils.removeSitesBasedOnFreqIgnoreGapsMissing(geno, 0.05, 210);
		SimplePhenotype temp = NumericalGenotypePlugin.collapseTransform(geno);
		JTable myTblTraits = createTraitTable(temp);
		int nrows = myTblTraits.getRowCount();
		myTblTraits.setRowSelectionInterval(0, nrows-1);
		Datum phenoDatum = new Datum("name", temp, "comment");
		ImputePanel imputePanel = new ImputePanel(phenoDatum);
		imputedDatum = imputePanel.createImputedData(myTblTraits);
	}

	@After
	public void tearDown() throws Exception {
	}
	
	@Test
	public void testPCA() throws Exception {
		PCAPanel pcaPanel = new PCAPanel(imputedDatum, null);
		JTable myTblTraits = createTraitTable((Phenotype) imputedDatum.getData());
		int nrows = myTblTraits.getRowCount();
		myTblTraits.setRowSelectionInterval(0, nrows-1);
		
		List<Datum> theResults = pcaPanel.createPCAData(myTblTraits);
		
		//There should be three tables
		Assert.assertEquals(3, theResults.size());
		
		//the phenotype -- pca axes
		TableReport pcaTable = (TableReport) theResults.get(0).getData();
		String filename = "test/datafiles/pcatable.txt";
		
		//write 
		writeTable(pcaTable, filename);
		
		//compare
		String message = matchesTable(pcaTable, filename);
		Assert.assertTrue("pca table does not match: " + message, message.length() == 0);
		
		//eigenvectors
		TableReport eigenvectorTable = (TableReport) theResults.get(1).getData();
		filename = "test/datafiles/eigenvector_table.txt";
		
		//write 
		writeTable(eigenvectorTable, filename);
		
		//compare
		message = matchesTable(eigenvectorTable, filename);
		Assert.assertTrue("eigenvector table does not match: " + message, message.length() == 0);
		
		//compare cumulative variance
		TableReport varianceTable = (TableReport) theResults.get(2).getData();
		filename = "test/datafiles/cumulative_variance_table.txt";
		
		//write 
		writeTable(varianceTable, filename);
		
		//compare
		message = matchesTable(varianceTable, filename);
		Assert.assertTrue("cumulative variance table does not match: " + message, message.length() == 0);
	}
	
	private String matchesTable(TableReport report, String filename) {
		Pattern tab = Pattern.compile("\t");
		try {
			BufferedReader br = new BufferedReader(new FileReader(filename));
			int nrecords = report.getRowCount();
			int ncols = report.getColumnCount();
			for (int r = 0; r < nrecords; r++) {
				String[] data = tab.split(br.readLine());
				for (int c = 0; c < ncols; c++) {
					if(!data[c].equals(report.getValueAt(r, c).toString())) {
						return "Data did not match at row " + r + ", column " + c;
					}
				}
				
			}
			return "";
		} catch (IOException e) {
			e.printStackTrace();
		}
		return "IO error";
	}
	
	private void writeTable(TableReport report, String filename) {
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(filename));
			int nrecords = report.getRowCount();
			int ncols = report.getColumnCount();
			for (int r = 0; r < nrecords; r++) {
				StringBuilder sb = new StringBuilder(report.getValueAt(r, 0).toString());
				for (int c = 1; c < ncols; c++) {
					sb.append("\t").append(report.getValueAt(r, c).toString());
				}
				bw.write(sb.toString());
				bw.newLine();
			}
			bw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
    private JTable createTraitTable(Phenotype ca) {
        String[] tableColumnNames = {"Column", "Percent Missing Data"};
        Object[] missingData = Conversion.getPercentMissingData(ca, 2);
        Object[][] tableData = new Object[ca.getNumberOfTraits()][2];
        for (int i = 0; i < ca.getNumberOfTraits(); i++) {
            tableData[i][0] = ca.getTrait(i).getName() + "." + ca.getTrait(i).getProperty(Trait.FACTOR_ENV);
            tableData[i][1] = missingData[i];
        }
        JTable tblAvailableColumns = new JTable(tableData, tableColumnNames);
        tblAvailableColumns.setEnabled(true);
        tblAvailableColumns.setRowSelectionAllowed(true);
        tblAvailableColumns.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
        return tblAvailableColumns;
    }

}
