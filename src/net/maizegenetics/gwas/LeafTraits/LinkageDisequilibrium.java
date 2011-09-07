package net.maizegenetics.gwas.LeafTraits;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.regex.Pattern;

import net.maizegenetics.jGLiM.ModelEffect;

import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.linalg.QRDecomposition;
import cern.jet.math.Functions;

public class LinkageDisequilibrium {
	ArrayList<int[]> founderGenotypes;
	double[][] genotypes;
	final int numberOfTaxa = 26;
	final int numberOfColumns = 38;
	final int poscolumn = 3;
	AGPMap theAGPMap = null;
	int chromosome;
	int firstMarker;
	int maxMarker;
	int[] pops;
	QRDecomposition theQR;
	ModelEffect mepop;
	
	//input files
	public final static String[] inputFiles = new String[]{
		"C:/Projects/NAM/association/Fastphase/fastphase_chr1.txt",
		"C:/Projects/NAM/association/Fastphase/fastphase_chr2.txt",
		"C:/Projects/NAM/association/Fastphase/fastphase_chr3.txt",
		"C:/Projects/NAM/association/Fastphase/fastphase_chr4.txt",
		"C:/Projects/NAM/association/Fastphase/fastphase_chr5.txt",
		"C:/Projects/NAM/association/Fastphase/fastphase_chr6.txt",
		"C:/Projects/NAM/association/Fastphase/fastphase_chr7.txt",
		"C:/Projects/NAM/association/Fastphase/fastphase_chr8.txt",
		"C:/Projects/NAM/association/Fastphase/fastphase_chr9.txt",
		"C:/Projects/NAM/association/Fastphase/fastphase_chr10.txt"};
	
	public final static String namMapFile = "C:/Projects/NAM/leaf traits/data/NAM_Map_20090730.txt";
//	public final static String outfilename = "C:/Projects/NAM/leaf traits/ld/pairsAdjustedForPop.txt";
	public final static String outfilename = "C:/Projects/NAM/leaf traits/ld/pairsAdjustedForPop_jl.txt";
//	public final static String outfilename = "C:/Projects/NAM/leaf traits/ld/pairsAdjustedForPop.txt";
	
	public LinkageDisequilibrium(int chromosome) {
		this.chromosome = chromosome;
		getAGPMap();
		importSNPs();
		importNamMarkersforMap();
//		calculatePairwiseLDandDistance();
		calculateCorrectedLD_restricted();
	}
	
	public void calculatePairwiseLDandDistance() {
		//shuffle genotypes
		Collections.shuffle(founderGenotypes);
		
		//calculate LD and distance for pairs between sets of size n;
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(outfilename));
			
			int n = 500;
			int[] set = new int[n * 2];
			int count = 0;
			int genoNumber = 0;
			
			while (count < n*2) {
//				int[] parents = int2Genotype(founderGenotypes.get(genoNumber)[1]);
				set[count++] = genoNumber;
				genoNumber++;
			}
			for (int i = 0; i < n; i++) {
				double[] geno1 = projectSNP(founderGenotypes.get(set[i]));
				int dist1 = founderGenotypes.get(set[i])[0];
				double dist1cm = theAGPMap.getCmFromPosition(chromosome, dist1);
//				double[] geno1 = projectSNP(new int[]{dist1,134217727});
				for (int j = n; j < 2*n; j++) {
					double[] geno2 = projectSNP(founderGenotypes.get(set[j]));
					int dist2 = founderGenotypes.get(set[j])[0];
					double dist2cm = theAGPMap.getCmFromPosition(chromosome, dist2);
//					double[] geno2 = projectSNP(new int[]{dist2,134217727});
					double r2 = rsquare(geno1, geno2);
					bw.write(Integer.toString(dist1));
					bw.write(",");
					bw.write(Integer.toString(dist2));
					bw.write(",");
					bw.write(Integer.toString(Math.abs(dist1 - dist2)));
					bw.write(",");
					bw.write(Double.toString(r2));
					bw.write(",");
					bw.write(Double.toString(dist1cm));
					bw.write(",");
					bw.write(Double.toString(dist2cm));
					bw.write(",");
					bw.write(Double.toString(Math.abs(dist1cm - dist2cm)));
					bw.newLine();
				}
			}
			bw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	public void calculateLDandDistance() {
		//shuffle genotypes
		Collections.shuffle(founderGenotypes);
		
		//calculate LD and distance for pairs between sets of size n;
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(outfilename));
			
			int n = 10000;
			int[] set = new int[n];
			int count = 0;
			int genoNumber = 0;
			int incr = 10;
			while (count < n) {
				int[] parents = int2Genotype(founderGenotypes.get(genoNumber)[1]);
				set[count] = genoNumber;
				genoNumber += incr;
				count++;
			}
			int dist1 = founderGenotypes.get(set[3000])[0];
			double dist1cm = theAGPMap.getCmFromPosition(10, dist1);
			double[] geno1 = projectSNP(new int[]{dist1,134217727});
			
			for (int i = 0; i < n; i++) {
				int dist2 = founderGenotypes.get(set[i])[0];
				double dist2cm = theAGPMap.getCmFromPosition(10, dist2);
				double[] geno2 = projectSNP(new int[]{dist2,134217727});
				double r2 = rsquare(geno1, geno2);
				bw.write(Integer.toString(dist1));
				bw.write(",");
				bw.write(Integer.toString(dist2));
				bw.write(",");
				bw.write(Integer.toString(dist2 - dist1));
				bw.write(",");
				bw.write(Double.toString(r2));
				bw.write(",");
				bw.write(Double.toString(dist1cm));
				bw.write(",");
				bw.write(Double.toString(dist2cm));
				bw.write(",");
				bw.write(Double.toString(dist2cm - dist1cm));
				bw.newLine();

			}
			bw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	public void calculateCorrectedLD() {
		//shuffle genotypes
		Collections.shuffle(founderGenotypes);

		//calculate LD and distance for pairs between sets of size n;
		initializelQRFromPop();
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(outfilename));
			
			int n = 500;
			
			for (int i = 0; i < n; i++) {
				double[] geno1 = projectSNP(founderGenotypes.get(i));
				int dist1 = founderGenotypes.get(i)[0];
				double dist1cm = theAGPMap.getCmFromPosition(chromosome, dist1);
				for (int j = n; j < 2*n; j++) {
					double[] geno2 = projectSNP(founderGenotypes.get(j));
					int dist2 = founderGenotypes.get(j)[0];
					double dist2cm = theAGPMap.getCmFromPosition(chromosome, dist2);
					double r2 = correctedRsquare(geno1, geno2);
					bw.write(Integer.toString(Math.abs(dist1 - dist2)));
					bw.write(",");
					bw.write(Double.toString(Math.abs(dist1cm - dist2cm)));
					bw.write(",");
					bw.write(Double.toString(r2));
					bw.newLine();
				}
			}
			bw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}		
	}
	
	public void calculateCorrectedLD_restricted() {
		//shuffle genotypes
		Collections.shuffle(founderGenotypes);

		//calculate LD and distance for pairs between sets of size n;
		initializelQRFromPop();
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(outfilename));
			int n = 500;
			int[] set = new int[2*n];
			int count = 0;
			int genoNumber = 0;
			while (count < 2*n) {
				int[] parents = int2Genotype(founderGenotypes.get(genoNumber)[1]);
				int sumOfParents = 0;
				for (int i = 0; i < 26; i++) sumOfParents += parents[i];
				if (sumOfParents > 2 && sumOfParents < 26) {
					set[count++] = genoNumber; 
				}
				genoNumber++;
			}
			
			for (int i = 0; i < n; i++) {
				int[] fg = founderGenotypes.get(set[i]);
				fg[1] = 67108863;
				double[] geno1 = projectSNP(fg);
				
				int dist1 = founderGenotypes.get(set[i])[0];
				double dist1cm = theAGPMap.getCmFromPosition(chromosome, dist1);
				for (int j = n; j < 2*n; j++) {
					fg = founderGenotypes.get(set[j]);
					fg[1] = 67108863;
					double[] geno2 = projectSNP(founderGenotypes.get(set[j]));
					
					int dist2 = founderGenotypes.get(set[j])[0];
					double dist2cm = theAGPMap.getCmFromPosition(chromosome, dist2);
					double delta = Math.abs(dist1cm - dist2cm);
					if (delta < 500) {
						double r2 = correctedRsquare(geno1, geno2);
						bw.write(Integer.toString(Math.abs(dist1 - dist2)));
						bw.write(",");
						bw.write(Double.toString(Math.abs(dist1cm - dist2cm)));
						bw.write(",");
						bw.write(Double.toString(r2));
						bw.newLine();
					}
				}
			}
			bw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}		
	}
	
	public void importSNPs() {
		founderGenotypes = new ArrayList<int[]>();
		File chrfile = new File(inputFiles[chromosome - 1]);
		Pattern sep = Pattern.compile("\\s");
		int firstCol = numberOfColumns - numberOfTaxa;
		try {
			BufferedReader br = new BufferedReader(new FileReader(chrfile));
			String inputLine;
			br.readLine();
			while ((inputLine = br.readLine()) != null) {
				String[] data = sep.split(inputLine);
				if (data.length == numberOfColumns) {
					int[] parents = new int[numberOfTaxa];
					for (int i = 0; i < numberOfTaxa; i++) parents[i] = Integer.parseInt(data[firstCol + i]);
					int geno = genotype2Int(parents);
					int pos = Integer.parseInt(data[poscolumn]);
					founderGenotypes.add(new int[]{pos, geno});
				}
			}
			
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	public int genotype2Int(int[] parents) {
		int mult = 1;
		int result = 0;
		for (int i = 0; i < numberOfTaxa; i++) {
			result += mult * parents[i];
			mult *= 2;
		}
		return result;
	}
	
	public int[] int2Genotype(int geno) {
		int mult = 1;
		int[] parents = new int[numberOfTaxa];
		for (int i = 0; i < numberOfTaxa; i++) {
			parents[i] = mult & geno;
			if (parents[i] > 0) parents[i] = 1;
			mult *= 2;
		}
		return parents;
	}
	
	public double[] projectSNP(int[] SNPinfo) {
		int pos = SNPinfo[0];
		int[] parents = int2Genotype(SNPinfo[1]);
		Object[][] markers = theAGPMap.getInterval(chromosome, pos);
		int left = AGPMap.getMarkerPosition(markers[0]);
		if (left == -1) left = 0;
		int right = AGPMap.getMarkerPosition(markers[1]);
		if (right == -1) right = ModelFitter.chromosomeLength[chromosome - 1];
		
		//proportion of distance of snp between left and right markers
		double pd = ((double)(pos - left)) / ((double)(right - left));
		
		int leftmarker = AGPMap.getMarkerNumber(markers[0]);
		if (leftmarker == -1) leftmarker = 0; //telomere
		else leftmarker = leftmarker - firstMarker + 1;
		int rightmarker = AGPMap.getMarkerNumber(markers[1]);
		if (rightmarker == -1) rightmarker = maxMarker; //telomere
		else rightmarker = rightmarker - firstMarker + 1;
		
		int nSamples = genotypes.length;
		double[] snpvalues = new double[nSamples];
		for (int i = 0; i < nSamples; i++) {
			snpvalues[i] = genotypes[i][leftmarker] * (1 - pd) + genotypes[i][rightmarker] * pd;
			snpvalues[i] = snpvalues[i] * parents[pops[i] - 1];
		}
		return snpvalues;
	}
	
	public void getAGPMap() {
		theAGPMap = new AGPMap(new File(namMapFile));
	}

	protected void importNamMarkersforMap() {
		String dir = "C:/Projects/NAM/leaf traits/markers/";
		String rils = "imputedMarkersGWAS.chr";
		String ibm = "imputedIBMMarkersGWAS.chr";
		String ending = ".082809.txt";
		
		String rilfilename = dir + rils + chromosome + ending;
		String ibmfilename = dir + ibm + chromosome + ending;
		
		try {
			BufferedReader br;
			String inputLine;
			br = new BufferedReader(new FileReader(rilfilename));
			int nRows = 0;
			String[] parsedLine = br.readLine().split("\t");
			int nCols = parsedLine.length - 1;
			firstMarker = Integer.parseInt(parsedLine[2].substring(1));
			maxMarker = nCols - 1;
			while ((inputLine = br.readLine()) != null) if (inputLine.length() > 10) nRows++;
			br.close();
						
			br = new BufferedReader(new FileReader(ibmfilename));
			br.readLine();
			while ((inputLine = br.readLine()) != null) if (inputLine.length() > 10) nRows++;
			br.close();
			
			genotypes = new double[nRows][nCols];
			pops = new int[nRows];
			br = new BufferedReader(new FileReader(rilfilename));
			br.readLine();
			int count = 0;
			while ((inputLine = br.readLine()) != null) {
				parsedLine = inputLine.split("\t");
				for (int c = 0; c < nCols; c++) genotypes[count][c] = Double.parseDouble(parsedLine[c + 1]);
				pops[count] = Integer.parseInt(parsedLine[0].substring(1,4));
				count++;
			}
			br.close();
			
			br = new BufferedReader(new FileReader(ibmfilename));
			br.readLine();
			while ((inputLine = br.readLine()) != null) {
				parsedLine = inputLine.split("\t");
				for (int i = 0; i < nCols; i++) genotypes[count][i] = Double.parseDouble(parsedLine[i + 1]);
				pops[count] = 17;
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
	
	public double rsquare(double[] a, double[] b) {
		int n = a.length;
		double suma = 0;
		double sumsqa = 0;
		double sumb = 0;
		double sumsqb = 0;
		double sumprod = 0;
		for (int i = 0; i < n; i++) {
			suma += a[i];
			sumb += b[i];
			sumsqa += a[i] * a[i];
			sumsqb += b[i] * b[i];
			sumprod += a[i] * b[i];
		}
		double cov = (sumprod - suma * sumb / n) / (n - 1);
		double vara = (sumsqa - suma * suma / n) / (n - 1);
		double varb = (sumsqb - sumb * sumb / n) / (n - 1);
		return cov * cov / vara / varb;
	}
	
	public double correctedRsquare(double[] a, double[] b) {
		DoubleMatrix2D geno = DoubleFactory2D.dense.make(26,2);
		geno.viewColumn(0).assign(mepop.getXTy(a));
		geno.viewColumn(1).assign(mepop.getXTy(b));
		DoubleMatrix2D beta = theQR.solve(geno);
		double[] ahat = mepop.getyhat(beta.viewColumn(0)).toArray();
		double[] bhat = mepop.getyhat(beta.viewColumn(1)).toArray();
		int n = a.length;
		for (int i = 0; i < n; i++) {
			a[i] -= ahat[i];
			b[i] -= bhat[i];
		}
		
		return rsquare(a,b);
	}
	
	public void initializelQRFromPop() {
		ArrayList<Integer> popList = new ArrayList<Integer>();
		for (int pop:pops) popList.add(pop);
		mepop = new ModelEffect(ModelEffect.getIntegerLevels(popList));
		theQR = new QRDecomposition(mepop.getXTX());
	}
	
	
}
