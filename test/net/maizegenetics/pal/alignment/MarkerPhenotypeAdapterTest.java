package net.maizegenetics.pal.alignment;

import static org.junit.Assert.*;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;

import junit.framework.Assert;

import net.maizegenetics.pal.alignment.CombinePhenotype;
import net.maizegenetics.pal.alignment.Phenotype;
import net.maizegenetics.pal.alignment.ReadPhenotypeUtils;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

public class MarkerPhenotypeAdapterTest {
	MarkerPhenotypeAdapter theAdapter;
	ArrayList<double[]> expectedPhenotypes;
	ArrayList<double[]> expectedCovariates;
	ArrayList<String[]> expectedMarkers;
	
	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
	}

	@AfterClass
	public static void tearDownAfterClass() throws Exception {
	}

	@Before
	public void setUp() throws Exception {
		String traitFile = "test/datafiles/test_traits.txt";
		String markerFile = "test/datafiles/diploid_SSR.txt";
		String popFile = "test/datafiles/test_qmatrix.txt";
		expectedPhenotypes = new ArrayList<double[]>();
		expectedCovariates = new ArrayList<double[]>();
		expectedMarkers = new ArrayList<String[]>();
		
		BufferedReader br = new BufferedReader(new FileReader(traitFile));
		br.readLine();
		br.readLine();
		br.readLine();
		String inputline = br.readLine();
		while (inputline != null) {
			String[] parsedline = inputline.split("\\s+");
			double[] data = new double[5];
			for (int i = 0; i < 5; i++) data[i] = Double.parseDouble(parsedline[i+1]);
			expectedPhenotypes.add(data);
			inputline = br.readLine();
		}
		br.close();
		
		br = new BufferedReader(new FileReader(popFile));
		br.readLine();
		br.readLine();
		inputline = br.readLine();
		while (inputline != null) {
			String[] parsedline = inputline.split("\\s+");
			double[] data = new double[3];
			for (int i = 0; i < 3; i++) data[i] = Double.parseDouble(parsedline[i+1]);
			expectedCovariates.add(data);
			inputline = br.readLine();
		}
		br.close();

		br = new BufferedReader(new FileReader(markerFile));
		br.readLine();
		br.readLine();
		for  (int i = 0; i < 30; i++) {
			inputline = br.readLine();
			String[] parsedline = inputline.split("\\s+");
			expectedMarkers.add(parsedline);
		}
		br.close();

		Phenotype testpheno = ReadPhenotypeUtils.readNumericalAlignment(traitFile);
		Phenotype testmarkers = ReadPhenotypeUtils.readPolymorphismAlignment(markerFile);
		Phenotype testcov = ReadPhenotypeUtils.readNumericalAlignment(popFile);
		for (int i = 0; i < 3; i++) {
			testcov.getTrait(i).setType(Trait.TYPE_COVARIATE);
		}
				
		Phenotype combinedData = CombinePhenotype.getInstance(new Phenotype[]{testpheno, testcov, testmarkers}, false);
		theAdapter = new MarkerPhenotypeAdapter(combinedData);
	}

	@After
	public void tearDown() throws Exception {
	}

	@Test
	public void testGetNumberOfRows() {
		;
	}

	@Test
	public void testGetNumberOfPhenotypes() {
		assertEquals(90, theAdapter.getNumberOfRows(0));
		assertEquals(60, theAdapter.getNumberOfRows(1));
	}

	@Test
	public void testGetPhenotypeValues() {
		double[] values = theAdapter.getPhenotypeValues(0);
		double[] expectedValues = new double[90];
		int count = 0;
		for (double[] dbldata : expectedPhenotypes) {
			expectedValues[count] = dbldata[0];
			expectedValues[count + 30] = dbldata[1];
			expectedValues[count++ + 60] = dbldata[2];
		}
		for (int i = 0; i < 90; i++) assertEquals(expectedValues[i], values[i], 1e-10);
		
		count = 0;
		values = theAdapter.getPhenotypeValues(1);
		for (double[] dbldata : expectedPhenotypes) {
			expectedValues[count] = dbldata[3];
			expectedValues[count + 30] = dbldata[4];
			expectedValues[count++ + 60] = Double.NaN;
		}
		for (int i = 0; i < 60; i++) assertEquals(expectedValues[i], values[i], 1e-10);
	}

	@Test
	public void testGetPhenotypeName() {
		Assert.assertEquals("trt1",theAdapter.getPhenotypeName(0));
		Assert.assertEquals("trt2", theAdapter.getPhenotypeName(1));
	}

	@Test
	public void testGetNumberOfFactors() {
		Assert.assertEquals(1, theAdapter.getNumberOfFactors());
	}

	@Test
	public void testGetFactorValues() {
		String[] factorval = theAdapter.getFactorValues(0,0);
		for (int i = 0; i < 30; i++) assertEquals("loc1", factorval[i]);
		for (int i = 30; i < 60; i++) assertEquals("loc2", factorval[i]);
		for (int i = 60; i < 90; i++) assertEquals("loc3", factorval[i]);
		factorval = theAdapter.getFactorValues(1, 0);
		for (int i = 0; i < 30; i++) assertEquals("loc1", factorval[i]);
		for (int i = 30; i < 60; i++) assertEquals("loc2", factorval[i]);
	}

	@Test
	public void testGetFactorName() {
		assertEquals("env", theAdapter.getFactorName(0));
	}

	@Test
	public void testGetNumberOfCovariates() {
		assertEquals(3,theAdapter.getNumberOfCovariates());
	}

	@Test
	public void testGetCovariateValues() {
		double[] values = theAdapter.getCovariateValues(0,0);
		double[] expectedValues = new double[90];
		int count = 0;
		for (double[] dbldata : expectedCovariates) {
			expectedValues[count] = dbldata[0];
			expectedValues[count + 30] = dbldata[0];
			expectedValues[count++ + 60] = dbldata[0];
		}
		for (int i = 0; i < 90; i++) assertEquals(expectedValues[i], values[i], 1e-10);
		
		values = theAdapter.getCovariateValues(0,1);
		count = 0;
		for (double[] dbldata : expectedCovariates) {
			expectedValues[count] = dbldata[1];
			expectedValues[count + 30] = dbldata[1];
			expectedValues[count++ + 60] = dbldata[1];
		}
		for (int i = 0; i < 90; i++) assertEquals(expectedValues[i], values[i], 1e-10);
		
		values = theAdapter.getCovariateValues(0,2);
		count = 0;
		for (double[] dbldata : expectedCovariates) {
			expectedValues[count] = dbldata[2];
			expectedValues[count + 30] = dbldata[2];
			expectedValues[count++ + 60] = dbldata[2];
		}
		for (int i = 0; i < 90; i++) assertEquals(expectedValues[i], values[i], 1e-10);
		
	}

	@Test
	public void testGetCovariateName() {
		assertEquals("Q1", theAdapter.getCovariateName(0));
		assertEquals("Q2", theAdapter.getCovariateName(1));
		assertEquals("Q3", theAdapter.getCovariateName(2));
	}

	@Test
	public void testGetNumberOfMarkers() {
		int nmarkers = theAdapter.getNumberOfMarkers();
		assertEquals(6, nmarkers);
	}

	@Test
	public void testGetMarkerName() {
		assertEquals("S01", theAdapter.getMarkerName(0));
		assertEquals("S03", theAdapter.getMarkerName(2));
		assertEquals("S06", theAdapter.getMarkerName(5));
	}

	@Test
	public void testGetMarkerValues() {
		Object[] values = theAdapter.getMarkerValue(0, 0);
		String[] expectedValues = new String[90];
		int count = 0;
		for (String[] strdata : expectedMarkers) {
			String val = strdata[1];
			if (val.contains("?")) val = "?";
			expectedValues[count] = val;
			expectedValues[count + 30] = val;
			expectedValues[count++ + 60] = val;
		}
		for (int i = 0; i < 90; i++) {
			String msg = "failed at i = " + i;
			assertEquals(msg, expectedValues[i], values[i]);
		}

	}

}
