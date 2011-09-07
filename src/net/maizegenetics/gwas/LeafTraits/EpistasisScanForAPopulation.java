package net.maizegenetics.gwas.LeafTraits;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;

import net.maizegenetics.jGLiM.LinearModelUtils;
import net.maizegenetics.jGLiM.dm.SweepFast;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrix;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrixFactory;

public class EpistasisScanForAPopulation extends Thread {
	int numberOfSamples;
	final int numberOfMarkers = 1106;
	double[]pheno;
	double[][] geno;
	Random randomGenerator = new Random();
	double[][] Fval;
	EpistasisScan parent;
	String traitname;
	String outputFile;
	int pop;
	String poplabel;
	
	/**
	 * @param phenotype an array of phenotypes
	 * @param genotype	an array of marker scores. The first dimension is number of markers, the second dimension is number of samples.
	 * @param numberOfPermutations 
	 */
	public EpistasisScanForAPopulation(double[] phenotype, double[][] genotype, String phenoName, int pop) {
		numberOfSamples = phenotype.length;
		pheno = phenotype;
		geno = genotype;
		int nModels = numberOfMarkers * (numberOfMarkers - 1) / 2;
		Fval = new double[nModels][4];
		traitname = phenoName;
		this.pop = pop;
		poplabel = "0" + pop;
		poplabel = poplabel.substring(poplabel.length() - 2);
	}
	
	@Override
	public void run() {
		System.out.println("Starting analysis for pop " + pop);
		testAllMarkerPairs();
	}

	public void testAllMarkerPairs() {
		DoubleMatrix Y = DoubleMatrixFactory.DEFAULT.make(numberOfSamples, 1);
		double totalSS = 0;
		
		for (int i = 0; i < numberOfSamples; i++) {
			Y.set(i, 0, pheno[i]);
			totalSS += pheno[i] * pheno[i];
		}
		
		//do stepwise regression for the base model using 1e-3 cutoff
		boolean keepGoing = true;
		DoubleMatrix X = DoubleMatrixFactory.DEFAULT.make(numberOfSamples, 1, 1);
		while (keepGoing) {
			//add the first marker to X
			int bestmarker = 1;
			int lastColumn = X.numberOfColumns();
			X = X.concatenate(DoubleMatrixFactory.DEFAULT.make(numberOfSamples, 1, geno[0]), false);
			
			//test the first marker and set minSS equal to its residual.
			DoubleMatrix XtX = X.crossproduct(X);
			DoubleMatrix XtY = X.crossproduct(Y);
			SweepFast sf = new SweepFast(XtX, XtY, totalSS);
			sf.revg2sweep(0);
			sf.setDminFromA();
			for (int i = 1; i <= lastColumn; i++) sf.revg2sweep(i);
			double minSS = sf.getResidualSS();
			DoubleMatrix minX = X.copy();
			
			//test the rest of the markers and identify the one with the smallest residual.
			for (int m = 0; m < numberOfMarkers; m++) {
				for (int i = 0; i < numberOfSamples; i++) X.set(i,lastColumn,geno[m][i]);
				XtX = X.crossproduct(X);
				XtY = X.crossproduct(Y);
				sf = new SweepFast(XtX, XtY, totalSS);
				sf.revg2sweep(0);
				sf.setDminFromA();
				for (int i = 1; i <= lastColumn; i++) sf.revg2sweep(i);
				double resSS = sf.getResidualSS();
				if (resSS < minSS) {
					minSS = resSS;
					minX = X.copy();
					bestmarker = m;
				}
			}
			
			//if the p-value for this  marker < 1e-3 add it to the model and keepgoing
			double F,p;
			XtX = minX.crossproduct(minX);
			XtY = minX.crossproduct(Y);
			sf = new SweepFast(XtX, XtY, totalSS);
			sf.revg2sweep(0);
			sf.setDminFromA();
			for (int i = 1; i < lastColumn; i++) sf.revg2sweep(i);
			double reducedResidualSS = sf.getResidualSS();
			sf.revg2sweep(lastColumn);
			double residualSS = sf.getResidualSS();
			double errordf = numberOfSamples - minX.numberOfColumns();
			double errorMS = residualSS / errordf;
			double interactionMS = reducedResidualSS - residualSS;
			F = interactionMS / errorMS;
			try {
				p = LinearModelUtils.Ftest(F, 1, errordf);
			} catch(Exception e) {
				p = Double.NaN;
			}
			
			if (p < 1e-3) {
				X = minX;
			} else {
				System.out.println("Base model complete for pop" + pop);
				keepGoing = false;
			}
			
		}
		
		//remove the last column of X
		int n = X.numberOfColumns() - 1;
		int[] saveColumns = new int[n];
		for (int i = 0; i < n; i++) saveColumns[i] = i;
		X = X.getSelection(null, saveColumns);
		
		//define K for testing interaction
		int nXcolumns = X.numberOfColumns();
		int nKrows = nXcolumns + 3;
		DoubleMatrix K = DoubleMatrixFactory.DEFAULT.make(nKrows, 1, 0);
		K.set(nKrows - 1, 0, 1);
		DoubleMatrix X2 = DoubleMatrixFactory.DEFAULT.make(numberOfSamples, 3);
		DoubleMatrix[][] xtx = new DoubleMatrix[2][2];
		DoubleMatrix[] xty = new DoubleMatrix[2];
		xtx[0][0] = X.crossproduct(X);
		xty[0] = X.crossproduct(Y);
		
		BufferedWriter bw = null;
		try {
			String filename = outputFile + "Ftest_" + traitname + "_pop" + poplabel + ".txt";
			bw = new BufferedWriter(new FileWriter(filename));
			bw.write("Trait\tpopulation\tmarker1\tmarker2\tF\tpval\tlogp");
			bw.newLine();
		} catch (IOException e1) {
			e1.printStackTrace();
			System.exit(-1);
		}
		
		ArrayList<GLMResult> resultList = new ArrayList<GLMResult>();
		for (int m1 = 0; m1 < numberOfMarkers - 1; m1++) {
//			if (m1%1000 == 0) System.out.println("pop " + pop + ", marker 1 = " + m1);
			double[] marker1 = geno[m1];
			for (int m2 = m1 + 1; m2 < numberOfMarkers; m2++) {
				double F,p;
				double[] marker2 = geno[m2];
				for (int i = 0; i < numberOfSamples; i++) {
					X2.set(i, 0, marker1[i]);
					X2.set(i, 1, marker2[i]);
					X2.set(i, 2, marker1[i] * marker2[i]);
				}
				xtx[0][1] = X.crossproduct(X2);
				xtx[1][1] = X2.crossproduct(X2);
				xty[1] = X2.crossproduct(Y);
				SweepFast sf = new SweepFast(xtx, xty, totalSS);
				sf.revg2sweep(0);
				sf.setDminFromA();
				double reducedModeldf = 1;
				for (int i = 1; i < nKrows - 1; i++) {
					if (sf.revg2sweep(i)) reducedModeldf++;
				}
				double reducedResidualSS = sf.getResidualSS();
				int interactiondf = 0;
				if (sf.revg2sweep(nKrows - 1)) {
					interactiondf = 1;
				}
				
				if (interactiondf == 0) {
					F = Double.NaN;
					p = Double.NaN;
				} else {
					double residualSS = sf.getResidualSS();
					double errordf = numberOfSamples - reducedModeldf - interactiondf;
					double errorMS = residualSS / errordf;
					double interactionMS = reducedResidualSS - residualSS;
					F = interactionMS / errorMS;
					try {
						p = LinearModelUtils.Ftest(F, 1, errordf);
					} catch(Exception e) {
						p = Double.NaN;
					}
				}
				
				if (!Double.isNaN(p)) resultList.add(new GLMResult(m1, m2, F, p));
				
				//write the output
				StringBuilder sb = new StringBuilder(traitname);
				sb.append("\t").append(pop);
				sb.append("\t").append(m1);
				sb.append("\t").append(m2);
				sb.append("\t").append(F);
				sb.append("\t").append(p);
				sb.append("\t").append(-Math.log10(p));
				
				try {
					bw.write(sb.toString());
					bw.newLine();
				} catch (IOException e) {
					e.printStackTrace();
					System.exit(-1);
				}
			}
		}
		try {
			bw.close();
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		
		int nResults = resultList.size();
		Collections.sort(resultList);
		int maxsig = -1;
		for (int i = 0; i < nResults; i++) {
			if (resultList.get(i).p < 0.05*i/nResults) maxsig = i;
		}
		
		System.out.println("pop " + pop + ": number of significant tests = " + (maxsig + 1));
		if (maxsig >= 0) {
			String filename = outputFile + "bhresults_" + traitname + "_pop" + poplabel + ".txt"; 
			try {
				bw = new BufferedWriter(new FileWriter(filename));
				bw.write("Trait\tpopulation\tmarker1\tmarker2\tF\tpval\tlogp");
				bw.newLine();
				for (int i = 0; i <= maxsig; i++) {
					GLMResult result = resultList.get(i);
					StringBuilder sb = new StringBuilder(traitname);
					sb.append("\t").append(pop);
					sb.append("\t").append(result.m1);
					sb.append("\t").append(result.m2);
					sb.append("\t").append(result.F);
					sb.append("\t").append(result.p);
					sb.append("\t").append(-Math.log10(result.p));
					bw.write(sb.toString());
					bw.newLine();
				}
				bw.close();
			} catch (IOException e) {
				e.printStackTrace();
				System.exit(-1);
			}
		}
	}
	
	public void setOutputFile(String output) {
		outputFile = output;
	}
	
	private void shuffle(int[] index) {
		//uses Fisher-Yates shuffle to randomize the order of an index
	    // i is the number of items remaining to be shuffled.
		int n = index.length;
	    for (int i = n; i > 1; i--) {
	        // Pick a random element to swap with the i-th element.
	        int j = randomGenerator.nextInt(i);  // 0 <= j <= i-1 (0-based array)
	        // Swap array elements.
	        int tmp = index[j];
	        index[j] = index[i-1];
	        index[i-1] = tmp;
	    }
	}
	
	class GLMResult implements Comparable<GLMResult> {
		int m1;
		int m2;
		double F;
		double p;
		
		GLMResult(int m1, int m2, double F, double p) {
			this.m1 = m1;
			this.m2 = m2;
			this.F = F;
			this.p = p;
		}

		@Override
		public int compareTo(GLMResult other) {
			if (Double.isNaN(p) && Double.isNaN(other.p)) return 0;
			if (Double.isNaN(p)) return 1;
			if (Double.isNaN(other.p)) return -1;
			if (p < other.p) return -1;
			if (p > other.p) return 1;
			return 0;
		}
		
	}
	
}
