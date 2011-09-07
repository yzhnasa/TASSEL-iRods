package net.maizegenetics.gwas.LeafTraits;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.concurrent.Callable;
import java.util.regex.Pattern;

import cern.colt.matrix.DoubleMatrix1D;

import net.maizegenetics.jGLiM.AbstractLinearModel;
import net.maizegenetics.jGLiM.CovariateModelEffect;
import net.maizegenetics.jGLiM.LinearModelWithSweep;
import net.maizegenetics.jGLiM.ModelEffect;

public class ModelFitter extends Thread {
	HashMap<String, Integer> genotypeMap;
	double[][] genotypes;
	int firstMarker;
	static AGPMap theAGPMap = null;
	double[][] residuals;
	String[] residualSamples;
	int chromosome;
	int traitnumber;
	String thisTraitName;
	int totalSnps;
	int snpcount;
	FileNames files = null;
	
	
	public static final int[] chromosomeLength = new int[]{300239041,
		234752839,
		230558137,
		247095508,
		216915529,
		169254300,
		170974187,
		174515299,
		152350485,
		149686046};

//	public static final String namGenoFile = "C:/Projects/NAM/leaf traits/markers/ImputedMarkerGenotypes_leaf_traits_082809.txt";
	
	public ModelFitter(int chromosome, int trait) {
		traitnumber = trait;
		this.chromosome = chromosome;
		getAGPMap(null);
	}
	
	public ModelFitter(int chromosome, int trait, FileNames files) {
		traitnumber = trait;
		this.chromosome = chromosome;
		getAGPMap(files.agpmap);
		this.files = files;
	}

	public ModelFitter(int chromosome, FileNames files) {
		this.chromosome = chromosome;
		getAGPMap(files.agpmap);
		this.files = files;
	}
	
	public void run() {
		init();
		long start = System.currentTimeMillis();

		testSnps();
		
		long elapsed = System.currentTimeMillis() - start;
		System.out.println("Elapsed time for leaf traits = " + elapsed);
	}
	
	public void init() {
		if (files == null) importResiduals();
		else if (files.residuals != null) importResiduals(files.residuals);
		else importResiduals();
		importNamMarkersforMap_alt();
	}
	
	public synchronized void getAGPMap(File mapfile) {
		if (theAGPMap == null) {
			theAGPMap = new AGPMap(mapfile);
		}
	}
	
	protected void importNamMarkersforMap() {
		String dir = "C:/Projects/NAM/leaf traits/markers/";
		String rils = "imputedMarkersGWAS.chr";
		String ibm = "imputedIBMMarkersGWAS.chr";
		String ending = ".082809.txt";
		
		String rilfilename = dir + rils + chromosome + ending;
		String ibmfilename = dir + ibm + chromosome + ending;
		
		genotypeMap = new HashMap<String, Integer>();
		try {
			BufferedReader br;
			if (files == null) br = new BufferedReader(new FileReader(rilfilename));
			else br = new BufferedReader(new FileReader(files.namMarkersByChr[chromosome - 1]));
			int nRows = 0;
			String[] parsedLine = br.readLine().split("\t");
			int nCols = parsedLine.length - 1;
			firstMarker = Integer.parseInt(parsedLine[2].substring(1));
			while (br.readLine() != null) nRows++;
			br.close();
						
			if (files == null) br = new BufferedReader(new FileReader(ibmfilename));
			else br = new BufferedReader(new FileReader(files.ibmMarkersByChr[chromosome - 1]));
			while (br.readLine() != null) nRows++;
			br.close();
			
			genotypes = new double[nRows][nCols];
			if (files == null) br = new BufferedReader(new FileReader(rilfilename));
			else br = new BufferedReader(new FileReader(files.namMarkersByChr[chromosome - 1]));
			br.readLine();
			String inputLine;
			int count = 0;
			while ((inputLine = br.readLine()) != null) {
				parsedLine = inputLine.split("\t");
				genotypeMap.put(parsedLine[0], count);
				for (int c = 0; c < nCols; c++) genotypes[count][c] = Double.parseDouble(parsedLine[c + 1]);
				
				count++;
			}
			br.close();
			
			if (files == null) br = new BufferedReader(new FileReader(ibmfilename));
			else br = new BufferedReader(new FileReader(files.ibmMarkersByChr[chromosome - 1]));
			br.readLine();
			while ((inputLine = br.readLine()) != null) {
				parsedLine = inputLine.split("\t");
				genotypeMap.put(parsedLine[0], count);
				for (int i = 0; i < nCols; i++) genotypes[count][i] = Double.parseDouble(parsedLine[i + 1]);
				count++;
			}
			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
			System.exit(-1);
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		
	}
	
	protected void importNamMarkersforMap_alt() {
		HashMap<String, Integer> sampleNameMap = new HashMap<String, Integer>();
		int nRows = residualSamples.length;
		for (int i = 0; i < nRows; i++) sampleNameMap.put(residualSamples[i], new Integer(i));

		genotypeMap = new HashMap<String, Integer>();
		try {
			BufferedReader br;
//			br = new BufferedReader(new FileReader(files.namMarkersByChr[chromosome - 1]));
//			int nRows = 0;
//			String[] parsedLine = br.readLine().split("\t");
//			int nCols = parsedLine.length - 1;
//			firstMarker = Integer.parseInt(parsedLine[2].substring(1));
//			while (br.readLine() != null) nRows++;
//			br.close();
						
			br = new BufferedReader(new FileReader(files.namMarkersByChr[chromosome - 1]));
			String inputLine = br.readLine();
			String[] parsedLine = inputLine.split("\t");
			int nMarkers = parsedLine.length - 1;
			firstMarker = Integer.parseInt(parsedLine[2].substring(1));
			genotypes = new double[nRows][nMarkers];
			for (int i = 0; i < nRows; i++) {
				for (int j = 0; j < nMarkers; j++) {
					genotypes[i][j] = Double.NaN;
				}
			}
			
			int count = 0;
			while ((inputLine = br.readLine()) != null) {
				parsedLine = inputLine.split("\t");
				Integer sampleNumber = sampleNameMap.get(parsedLine[0]);
				if (sampleNumber != null) {
					int n = sampleNumber.intValue();
					for (int j = 0; j < nMarkers; j++) {
						genotypes[n][j] = Double.parseDouble(parsedLine[j + 1]);
					}
				}

//				for (int c = 0; c < nCols; c++) genotypes[count][c] = Double.parseDouble(parsedLine[c + 1]);
				
//				count++;
			}
			br.close();
			
			int totalcount = 0;
			for (int i = 0; i < nRows; i++) {
				count = 0;
				for (int j = 0; j < nMarkers; j++) {
					if (Double.isNaN(genotypes[i][j])) {
						count++;
						totalcount++;
					}
				}
				if (count > 0) System.out.println("missing values for genotype " + i + ", name = " + residualSamples[i]);
			}
			
			System.out.println("count of nans = " + totalcount);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
			System.exit(-1);
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		
	}

	protected void importResiduals() {
//		String[] traitname = new String[]{"leaf_length","leaf_width","leaf_angle","transformed_angle"};
//		String beginning = "C:/Projects/NAM/leaf traits/residuals/residuals_distancebased_for_";
//		String ending = ".txt";
//		String filename = beginning + traitname[traitnumber] + ending;
//		thisTraitName = traitname[traitnumber];
		thisTraitName = "angle";
		String filename = "C:/Projects/NAM/leaf traits/residuals/residuals_distancebased_for_transformed_angle.txt";
		File residuals = new File(filename);
		importResiduals(residuals);
	}
	
	protected void importResiduals(File residualfile) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(residualfile));
			int nRows = 0;
			br.readLine();
			while (br.readLine() != null) nRows++;
			br.close();
			
			residuals = new double[nRows][10];
			residualSamples = new String[nRows];
			br = new BufferedReader(new FileReader(residualfile));
			br.readLine();
			for (int i = 0; i < nRows; i++) {
				String[] parsedLine = br.readLine().split("\t");
				residualSamples[i] = parsedLine[0];
				for (int c = 0; c < 10; c++) {
					residuals[i][c] = Double.parseDouble(parsedLine[c + 1]);
				}
			}
			br.close();
			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
			System.exit(-1);
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}

	}
	
	public void testSnps() {
		SnpData snpdata = new SnpData(chromosome);
		totalSnps = snpdata.getNumberOfSnps();
		int maxmarker = genotypes[0].length - 1;
		int nSamples = residualSamples.length;
		double[] snpvalues = new double[nSamples];
		BufferedWriter bw = openOutputFile("C:/Projects/NAM/leaf traits/results/testSnps_" + thisTraitName + "_chr" + chromosome + ".txt");
		StringBuilder sb = new StringBuilder("Position\tF\tp\tlog(1/p)\teffect");
		for (int i = 1; i <= 26; i++) sb.append("\tpop").append(i);
		writeToOutput(sb.toString(), bw);
		
		snpcount = 0;
		while (snpdata.next()) {
			snpcount++;
			
			double[] parents = snpdata.getGenotype();
			int pos = snpdata.getPosition();
			Object[][] markers = theAGPMap.getInterval(chromosome, pos);
			int left = AGPMap.getMarkerPosition(markers[0]);
			if (left == -1) left = 0;
			int right = AGPMap.getMarkerPosition(markers[1]);
			if (right == -1) right = chromosomeLength[chromosome - 1];
			
			//proportion of distance of snp between left and right markers
			double pd = ((double)(pos - left)) / ((double)(right - left));
			
			int leftmarker = AGPMap.getMarkerNumber(markers[0]);
			if (leftmarker == -1) leftmarker = 0; //telomere
			else leftmarker = leftmarker - firstMarker + 1;
			int rightmarker = AGPMap.getMarkerNumber(markers[1]);
			if (rightmarker == -1) rightmarker = maxmarker; //telomere
			else rightmarker = rightmarker - firstMarker + 1;
			
			ArrayList<Integer> nonMissingEntries = new ArrayList<Integer>();
			ArrayList<Integer> nonMissingPop = new ArrayList<Integer>();
			
			for (int i = 0; i < nSamples; i++) {
				String sample = residualSamples[i];
				int g = genotypeMap.get(sample);
				int thispop = getPopulation(sample);
				
				snpvalues[i] = genotypes[g][leftmarker] * (1 - pd) + genotypes[g][rightmarker] * pd;
				snpvalues[i] = snpvalues[i] * parents[thispop - 1];
				if (!Double.isNaN(snpvalues[i])) {
					nonMissingEntries.add(i);
					nonMissingPop.add(thispop);
				}
			}
			
			//build and solve the model
			int N = nonMissingEntries.size();
			ArrayList<ModelEffect> effects = new ArrayList<ModelEffect>();
			int[] mean = new int[N];
			double[] snp = new double[N];
			double[] y = new double[N];
			
			int count = 0;
			for (Integer i : nonMissingEntries) {
				snp[count] = snpvalues[i];
				y[count] = residuals[i][chromosome - 1];
				count++;
			}
			effects.add(new ModelEffect(mean));
			int[] poplevels = ModelEffect.getIntegerLevels(nonMissingPop);
			ModelEffect mepop = new ModelEffect(poplevels);
			effects.add(mepop);
			effects.add(new CovariateModelEffect(snp));
			
			LinearModelWithSweep lmws = new LinearModelWithSweep(effects, y);
			double errSS = lmws.getErrorSS();
			double errdf = lmws.getErrordf();
			ArrayList effectSS = lmws.incrementalEffectSS();
			DoubleMatrix1D beta = lmws.getBeta();
			int n = beta.size();
			double effectsize = beta.getQuick(n - 1);
			
			double F = ((Double)effectSS.get(2)).doubleValue() / errSS * errdf;
			double p;
			try{p = AbstractLinearModel.Ftest(F, 1, errdf);} catch (Exception e) {p = Double.NaN;}
			
			sb = new StringBuilder();
			sb.append(pos);
			sb.append("\t").append(F);
			sb.append("\t").append(p);
			sb.append("\t").append(Math.log10(1/p));
			sb.append("\t").append(effectsize);
			for (int i = 0; i < 26; i++) sb.append("\t").append(parents[i]);
			
			writeToOutput(sb.toString(), bw);
		}
		
		closeOutput(bw);
		System.out.println("Finished chromosome " + chromosome);
		
	}
	
	protected BufferedWriter openOutputFile(String filename) {
		try {
			return new BufferedWriter(new FileWriter(filename));
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		return null;
	}

	protected BufferedWriter openOutputFile(File file) {
		try {
			return new BufferedWriter(new FileWriter(file));
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		return null;
	}
	
	protected BufferedWriter openOutputFile(File file, boolean append) {
		try {
			return new BufferedWriter(new FileWriter(file, append));
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		return null;
	}
	
	protected void writeToOutput(String out, BufferedWriter bw) {
		try {
			bw.write(out);
			bw.newLine();
			bw.flush();
		} catch (IOException e) {
			e.printStackTrace();
			try{bw.close();} catch (Exception e2) {}
			System.exit(-1);
		}
	}

	protected void closeOutput(BufferedWriter bw) {
		try {
			bw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public static int getPopulation(String sample) {
		if (sample.charAt(0) == 'M') return 17;
		return Integer.parseInt(sample.substring(1, 4));
	}
	
	public double getProgress() {
		return ((double) snpcount) / ((double) totalSnps);
	}

	public String getThisTraitName() {
		return thisTraitName;
	}

	public void setThisTraitName(String thisTraitName) {
		this.thisTraitName = thisTraitName;
	}
}
