package net.maizegenetics.stats.GLM;


import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.regex.Pattern;

import junit.framework.Assert;

import net.maizegenetics.baseplugins.FixedEffectLMPlugin;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.AnnotatedAlignmentUtils;
import net.maizegenetics.pal.alignment.CombinePhenotype;
import net.maizegenetics.pal.alignment.FilterPhenotype;
import net.maizegenetics.pal.alignment.MarkerPhenotype;
import net.maizegenetics.pal.alignment.Phenotype;
import net.maizegenetics.pal.alignment.ReadPhenotypeUtils;
import net.maizegenetics.pal.alignment.ReadSequenceAlignmentUtils;
import net.maizegenetics.pal.report.TableReport;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;

public class GLMTest {
	Alignment d8;
	Phenotype pop;
	Phenotype traits;
	
	@Before
	public void setUp() throws Exception {
		d8 = ReadSequenceAlignmentUtils.readBasicAlignments("test/datafiles/d8_sequence.phy", 40);
		d8 = AnnotatedAlignmentUtils.removeSitesOutsideRange(d8, 1000, 2465);
		d8 = AnnotatedAlignmentUtils.removeSitesBasedOnFreqIgnoreGapsMissing(d8, 0.05, 50);
		
		pop = ReadPhenotypeUtils.readGenericFile("test/datafiles/mdp_population_structure.txt");
		pop = FilterPhenotype.getInstance(pop, null, new int[]{0,1});
		traits = ReadPhenotypeUtils.readNumericalAlignment("test/datafiles/mdp_traits.txt");
	}

	/**
	 * Tests GLM with traits and population covariates
	 */
	@Test
	public void testGLMwithPopulation() {
		CombinePhenotype ph = CombinePhenotype.getInstance(traits, pop, false);
		MarkerPhenotype mp = MarkerPhenotype.getInstance(d8, ph, false);
		
		FixedEffectLMPlugin glm = new FixedEffectLMPlugin(null, false);
		DataSet input = new DataSet(new Datum("glm test", mp, "glm test"), null);
		DataSet ds = glm.performFunction(input);
		
		Datum datum = ds.getData(0);
		TableReport stats = (TableReport) datum.getData();
		String fname = "test/datafiles/glm_stats_output.txt";
		
		//used to create a new validation result as necessary, normally commented out
//		writeTable(stats, fname);
		
		String message = matchesTable(stats, fname);
		Assert.assertTrue("stats table: " + message, message.length() == 0);
		
		datum = ds.getData(1);
		TableReport estimates = (TableReport) datum.getData();
		fname = "test/datafiles/glm_alleles_output.txt";
		
		//used to create a new validation result as necessary, normally commented out
//		writeTable(estimates, fname);
		
		message = matchesTable(estimates, fname);
		Assert.assertTrue("alleles table: " + message, message.length() == 0);
	}
	
	/**
	 * Tests GLM with traits and population covariates
	 */
	@Test
	public void testGLMwithPopulationPermutations() {
		CombinePhenotype ph = CombinePhenotype.getInstance(traits, pop, false);
		MarkerPhenotype mp = MarkerPhenotype.getInstance(d8, ph, false);
		
		FixedEffectLMPlugin glm = new FixedEffectLMPlugin(null, false);
		DataSet input = new DataSet(new Datum("glm test", mp, "glm test"), null);
		glm.setPermute(true);
		glm.setNumberOfPermutations(100);
		glm.setRandomizer(100);
		DataSet ds = glm.performFunction(input);
		
		Datum datum = ds.getData(0);
		TableReport stats = (TableReport) datum.getData();
		String fname = "test/datafiles/glm_permutation_stats_output.txt";
		
		//used to create a new validation result as necessary, normally commented out
//		writeTable(stats, fname);
		
		String message = matchesTable(stats, fname);
		Assert.assertTrue("stats table: " + message, message.length() == 0);
		
		datum = ds.getData(1);
		TableReport estimates = (TableReport) datum.getData();
		fname = "test/datafiles/glm_alleles_output.txt";
		
		//used to create a new validation result as necessary, normally commented out
//		writeTable(estimates, fname);
		
		message = matchesTable(estimates, fname);
		Assert.assertTrue("alleles table: " + message, message.length() == 0);
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
	
	@After
	public void tearDown() throws Exception {
		
	}

}
