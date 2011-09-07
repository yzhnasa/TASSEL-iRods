package net.maizegenetics.gwas.LeafTraits;

import java.io.BufferedWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;

import net.maizegenetics.jGLiM.AbstractLinearModel;
import net.maizegenetics.jGLiM.CovariateModelEffect;
import net.maizegenetics.jGLiM.LinearModelWithSweep;
import net.maizegenetics.jGLiM.LinearModelforStepwiseRegression;
import net.maizegenetics.jGLiM.ModelEffect;
import net.maizegenetics.jGLiM.SimpleShuffler;

public class ModelFitterWithPermutations extends ModelFitter {
	ArrayList<double[]> permData;
	ArrayList<Integer> phenotypeIndex;
	ArrayList<Integer> populationList;
	ArrayList<Double> SST;
	double[] minp;
	double[] origdata;
	int numberOfPermutations;
	protected static HashMap<Integer, ArrayList<Integer>> popMap = null;
	boolean bootstrap = false;

	public ModelFitterWithPermutations (int chromosome, int trait) {
		super(chromosome, trait);
	}
	
	public ModelFitterWithPermutations(int chromosome, FileNames files) {
		super(chromosome, files);
		bootstrap = files.bootstrapPermutation;
		numberOfPermutations = files.iterations;
	}
	
	public void testSnps() {
		minp = new double[numberOfPermutations];
		for (int i = 0; i < numberOfPermutations; i++) minp[i] = 1;
		
 		createPopMap();
		createPermutedData();
		
		SnpData snpdata = new SnpDataImputed(chromosome, files);
		totalSnps = snpdata.getNumberOfSnps();
		int maxmarker = genotypes[0].length - 1;
		int nSamples = phenotypeIndex.size();
		double[] snpvalues = new double[nSamples];

		//create an index array into genotypes
		int[] sampleIndex = new int[nSamples];
		for (int i = 0; i < nSamples; i++) sampleIndex[i] = genotypeMap.get(residualSamples[phenotypeIndex.get(i)]);
		
		//set up model effects
		ModelEffect memean = new ModelEffect(new int[nSamples]);
		ModelEffect mepop = new ModelEffect(ModelEffect.getIntegerLevels(populationList));
		
		//set up xtx and xty arrays
		DoubleMatrix2D[][] xtx = new DoubleMatrix2D[3][3];
		xtx[0][0] = memean.getXTX();
		xtx[1][1] = mepop.getXTX();
		xtx[0][1] = memean.getX1TX2(mepop);
		xtx[1][0] = xtx[0][1].viewDice();
		DoubleMatrix2D[][] xty = new DoubleMatrix2D[3][1];
		
		int snpcount = 0;
		long time = System.currentTimeMillis();
		System.out.println("Testing snps...");

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
			
			for (int i = 0; i < nSamples; i++) {
				snpvalues[i] = genotypes[sampleIndex[i]][leftmarker] * (1 - pd) + genotypes[sampleIndex[i]][rightmarker] * pd;
				snpvalues[i] = snpvalues[i] * parents[populationList.get(i)];
			}
			
			CovariateModelEffect snpEffect = new CovariateModelEffect(snpvalues);
			xtx[2][2] = snpEffect.getXTX();
			xtx[0][2] = memean.getX1TX2(snpEffect);
			xtx[2][0] = xtx[0][2].viewDice();
			xtx[1][2] = mepop.getX1TX2(snpEffect);
			xtx[2][1] = xtx[1][2].viewDice(); 
			
			//firstpermutation
			double[] y = permData.get(0);
			xty[0][0] = DoubleFactory2D.dense.make(memean.getLevelSums(y), memean.getNumberOfLevels());
			xty[1][0] = DoubleFactory2D.dense.make(mepop.getLevelSums(y), mepop.getNumberOfLevels());
			xty[2][0] = DoubleFactory2D.dense.make(1,1, snpEffect.getSumProducts(y));
			LinearModelWithSweep lmws = new LinearModelWithSweep(xtx, xty, y);
			double F, pval, errordf;
			errordf = lmws.getErrordf();
			F = lmws.effectSS.get(2) / lmws.getErrorSS() * errordf;
			try{
				pval = AbstractLinearModel.Ftest(F, 1, errordf);
			} catch (Exception e) {
				pval = Double.NaN;
			}
			if (!Double.isNaN(pval)) {
				minp[0] = Math.min(minp[0], pval);
			}
			
			//test the snp for each of the permuted residuals
			DoubleMatrix2D G = lmws.getInverse();
			double g = 1 / G.getQuick(G.rows() - 1, G.rows() - 1);
			for (int p = 1; p < numberOfPermutations; p++) {
				y = permData.get(p);
				double ssy = 0;
				double sumy = 0;
				double sumprody = 0;
				int count = 0;
				for (double d : y) {
					ssy += d*d;
					sumy += d;
					sumprody += d*snpvalues[count++];
				}
				
				int npops = popMap.size();
				DoubleMatrix1D xtrany = DoubleFactory1D.dense.make(2 + npops);
				xtrany.setQuick(0, sumy);
				double[] popy = mepop.getLevelSums(y);
				for (int i = 0; i < npops; i++) xtrany.setQuick(i+1, popy[i]);
				xtrany.setQuick(xtrany.size()-1, sumprody);
				DoubleMatrix1D beta = G.zMult(xtrany, null);
				double modelSS = xtrany.zDotProduct(beta);
				double errorms = (ssy - modelSS) / errordf;
				
				double snpeffect = beta.getQuick(beta.size() - 1);
				F = snpeffect * snpeffect * g / errorms;
				try{
					pval = AbstractLinearModel.Ftest(F, 1, errordf);
				} catch (Exception e) {
					pval = Double.NaN;
				}
				if (!Double.isNaN(pval)) {
					minp[p] = Math.min(minp[p], pval);
				}
			}
			if (snpcount%10000 == 0) System.out.println("chr " + chromosome + "snpcount " + snpcount + ": elapsed time = " + (System.currentTimeMillis() - time)); //debug

		}
		

		//setup the output
		BufferedWriter bw = openOutputFile(files.chrmodel[chromosome - 1], true);
		for (double pval : minp) writeToOutput(Double.toString(pval), bw);
		closeOutput(bw);	
		
	}
	
	public void testSnps_alternate_method() {
		minp = new double[numberOfPermutations];
		for (int i = 0; i < numberOfPermutations; i++) minp[i] = 1;
		
		System.out.println("createPopMap"); //debug
 		createPopMap();
		System.out.println("createPermutedData"); //debug
		createPermutedData();
		
		SnpData snpdata = new SnpDataImputed(chromosome);
		totalSnps = snpdata.getNumberOfSnps();
		int maxmarker = genotypes[0].length - 1;
		int nSamples = phenotypeIndex.size();
		double[] snpvalues = new double[nSamples];

		//create an index array into genotypes
		int[] sampleIndex = new int[nSamples];
		for (int i = 0; i < nSamples; i++) sampleIndex[i] = genotypeMap.get(residualSamples[phenotypeIndex.get(i)]);
		
		//debug code
		double[] unpermutedy = new double[phenotypeIndex.size()];
		for (int i = 0; i < nSamples; i++) unpermutedy[i] = residuals[phenotypeIndex.get(i)][chromosome - 1];
		//debug code
		
		//set up model effects
		ModelEffect memean = new ModelEffect(new int[nSamples]);
		ModelEffect mepop = new ModelEffect(ModelEffect.getIntegerLevels(populationList));
		
		//set up xtx and xty arrays
		DoubleMatrix2D[][] xtx = new DoubleMatrix2D[3][3];
		xtx[0][0] = memean.getXTX();
		xtx[1][1] = mepop.getXTX();
		xtx[0][1] = memean.getX1TX2(mepop);
		xtx[1][0] = xtx[0][1].viewDice();
		DoubleMatrix2D[][] xty = new DoubleMatrix2D[3][1];
		
		int snpcount = 0;
		long time = System.currentTimeMillis();
		System.out.println("test snps"); //debug

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
			
			for (int i = 0; i < nSamples; i++) {
				snpvalues[i] = genotypes[sampleIndex[i]][leftmarker] * (1 - pd) + genotypes[sampleIndex[i]][rightmarker] * pd;
				snpvalues[i] = snpvalues[i] * parents[populationList.get(i)];
			}
			
			int[] popIndex = new int[populationList.size()];
			int popcnt = 0;
			for (Integer pop : populationList) popIndex[popcnt++] = pop.intValue();
			
			CovariateModelEffect snpEffect = new CovariateModelEffect(snpvalues);
			xtx[2][2] = snpEffect.getXTX();
			xtx[0][2] = memean.getX1TX2(snpEffect);
			xtx[2][0] = xtx[0][2].viewDice();
			xtx[1][2] = mepop.getX1TX2(snpEffect);
			xtx[2][1] = xtx[1][2].viewDice(); 
			
			
			//test the snp for each of the permuted residuals
			for (int p = 0; p < numberOfPermutations; p++) {
				double[] y;  //debug code modification
				if (p==0) y = unpermutedy;  //debug code modification
				else y = permData.get(p); //debug code modification 
				
				xty[0][0] = DoubleFactory2D.dense.make(memean.getLevelSums(y), memean.getNumberOfLevels());
				xty[1][0] = DoubleFactory2D.dense.make(mepop.getLevelSums(y), mepop.getNumberOfLevels());
				xty[2][0] = DoubleFactory2D.dense.make(1,1, snpEffect.getSumProducts(y));
				LinearModelWithSweep lmws = new LinearModelWithSweep(xtx, xty, y);
				double F, pval, errordf;
				errordf = lmws.getErrordf();
				F = lmws.effectSS.get(2) / lmws.getErrorSS() * errordf;
				try{
					pval = AbstractLinearModel.Ftest(F, 1, errordf);
				} catch (Exception e) {
					pval = Double.NaN;
				}
				if (!Double.isNaN(pval)) {
					minp[p] = Math.min(minp[p], pval);
				}
			}
			if (snpcount%10 == 0) System.out.println("snpcount " + snpcount + ": elapsed time = " + (System.currentTimeMillis() - time)); //debug

		}
		

		//setup the output
		BufferedWriter bw = openOutputFile(files.chrmodel[chromosome - 1], true);
		for (double pval : minp) writeToOutput(Double.toString(pval), bw);
		closeOutput(bw);
	}
	
	public void createPermutedData() {
		permData = new ArrayList<double[]> ();
		int npop = popMap.size();
		int nobs = 0;
		for (ArrayList<Integer> popList : popMap.values()) nobs += popList.size();
		
		phenotypeIndex = new ArrayList<Integer>();
		for (int i = 0; i < numberOfPermutations; i++) {
			permData.add(new double[nobs]);
		}
		
		int popstart = 0;
		for (ArrayList<Integer> popList : popMap.values()) {
			phenotypeIndex.addAll(popList);
			ArrayList<Integer> thispop = new ArrayList<Integer>(popList);
			int thispopsize = thispop.size();
			for (int p = 0; p < numberOfPermutations; p++) {
				Collections.shuffle(thispop);
				double[] data = permData.get(p);
				for (int i = 0; i < thispopsize; i++) {
					data[i + popstart] = residuals[thispop.get(i)][chromosome - 1];
				}
			}
			popstart += thispopsize;
		}
		
		populationList = new ArrayList<Integer>();
		int count = 0;
		for (Integer sample : phenotypeIndex) {
			populationList.add(getPopulation(residualSamples[sample]) - 1);
		}
	}
	
	public synchronized void createPopMap() {
		if (popMap == null) {
			popMap = new HashMap<Integer, ArrayList<Integer>>();
			int nSamples = residualSamples.length;
			for (int i = 0; i < nSamples; i++) {
				int ipop = getPopulation(residualSamples[i]) - 1;
				ArrayList<Integer> popList = popMap.get(ipop);
				if (popList == null) {
					popList = new ArrayList<Integer>();
					popMap.put(ipop, popList);
				}
				popList.add(i);
			}
		}
		
		if (bootstrap) {
			//replace each popList with a subsample of 150
			for (ArrayList<Integer> popList : popMap.values()) {
				int N = popList.size();
				if (N > 150) {
					Collections.shuffle(popList);
					for (int i = N-1; i >= 150; i--) popList.remove(i);
				}
				else {
					int portion = (int) (0.8 * N);
					Collections.shuffle(popList);
					int count = 0;
					for (int i = portion; i < N; i++) popList.set(i, popList.get(count++));
					for (int i = N; i < 150; i++) popList.add(popList.get(count++));
				}
			}
		}
	}

	public void setNumberOfPermutations(int n) {
		numberOfPermutations = n;
	}
}
