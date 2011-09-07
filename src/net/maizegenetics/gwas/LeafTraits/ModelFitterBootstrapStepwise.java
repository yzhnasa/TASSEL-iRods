package net.maizegenetics.gwas.LeafTraits;

import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Set;

import cern.colt.matrix.DoubleMatrix1D;
import cern.jet.random.Uniform;
import cern.jet.random.engine.MersenneTwister;
import cern.jet.random.sampling.RandomSampler;

import net.maizegenetics.jGLiM.AbstractLinearModel;
import net.maizegenetics.jGLiM.CovariateModelEffect;
import net.maizegenetics.jGLiM.LinearModelforStepwiseRegression;
import net.maizegenetics.jGLiM.ModelEffect;
import net.maizegenetics.jGLiM.Shuffler;
import net.maizegenetics.jGLiM.SimpleShuffler;

public class ModelFitterBootstrapStepwise extends ModelFitter {
	
	protected double enterLimit = 1e-6;
	protected boolean withReplacement;
	protected Uniform uniformDist;
	protected int iterations = 100;
	protected int maxsnps = 100;
	protected int startIteration = 0;
	
	protected static int[] sampleIndex = null;
	protected static HashMap<Integer, ArrayList<Integer>> popMap = null;
	protected static int[] populationIndex = null;
	protected File modelOutputFile;
	protected boolean permute;
	protected Shuffler shuffler = null;
		
	public ModelFitterBootstrapStepwise(int chromosome, FileNames files){
		super(chromosome, files);
		if (files.limits != null) enterLimit = files.limits[chromosome - 1];
		else if (!Double.isNaN(files.enterlimit)) enterLimit = files.enterlimit;
		iterations = files.iterations;
		if (files.replacement.startsWith("T")) withReplacement = true;
		else if (files.replacement.startsWith("t")) withReplacement = true;
		else withReplacement = false;
		maxsnps = files.maxsnps;
		permute = files.permute;
		startIteration = files.startIteration;
	}

	public void testSnps() {
		if (withReplacement) testSnpsWithReplacement();
		else testSnpsWithoutReplacement();
	}
	
	public void testSnpsWithReplacement() {
		SnpData snpdata;
		if (files == null) snpdata = new SnpDataImputed(chromosome); 
		else snpdata = new SnpDataImputed(chromosome, files);
		totalSnps = snpdata.getNumberOfSnps();
		int maxmarker = genotypes[0].length - 1;
		int nSamples = residualSamples.length;
		uniformDist = new Uniform(0, nSamples - 1, new MersenneTwister(new java.util.Date()));

		int[] mean = new int[nSamples];
		int[] popIndex = new int[nSamples];
		ArrayList<String> pops = new ArrayList<String>();
		int count = 0;
		for (String sample : residualSamples) {
			int ipop = getPopulation(sample) - 1;
			pops.add(SnpData.popnames[ipop]);
			popIndex[count++] = ipop;
		}

		//create an index array into genotypes
		int[] sampleIndex = new int[nSamples];
		for (int i = 0; i < nSamples; i++) sampleIndex[i] = genotypeMap.get(residualSamples[i]);
		
		//set up the mean and pop model effects
		ArrayList<ModelEffect> effects = new ArrayList<ModelEffect>();
		ModelEffect memean = new ModelEffect(mean);
		effects.add(memean);
		ModelEffect mepop = new ModelEffect(ModelEffect.getIntegerLevels(pops));
		effects.add(mepop);
		
		//open file for output
		BufferedWriter bwstep;
		if (files == null) {
			StringBuilder sb = new StringBuilder("C:/Projects/NAM/leaf traits/results/stepwise/stepwise_chr");
			sb.append(chromosome).append("_");
			sb.append(thisTraitName).append(".txt");
			bwstep = Utilities.openOutputFile(sb.toString());
		}
		else bwstep = Utilities.openOutputFile(files.chrsteps[chromosome - 1]);
		Utilities.flushToOutput("chromosome\tposition\tcM\tp-value,\titeration", bwstep);
		
		for (int iter = 0; iter < iterations; iter++) {
			//get data for this chromosome 
			double[] y = new double[nSamples];
			int[] bootstrapIndex = new int[nSamples];
			for (int i = 0; i < nSamples; i++) {
				int ndx = uniformDist.nextInt();
				y[i] = residuals[ndx][chromosome - 1];
				bootstrapIndex[i] = sampleIndex[ndx];
			}
			
			
			//create the base model
			LinearModelforStepwiseRegression lmsr = new LinearModelforStepwiseRegression(effects, y);
			
			//forward regression
			double minp;
			do {
				//create the base model
				minp = 1;
				int bestpos = -1;
				ModelEffect bestEffect = null;
				double maxF = 0;

				snpcount = 0;
				snpdata.reset();
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

					double[] snpvalues = new double[nSamples];
					for (int i = 0; i < nSamples; i++) {
						snpvalues[i] = genotypes[bootstrapIndex[i]][leftmarker] * (1 - pd) + genotypes[bootstrapIndex[i]][rightmarker] * pd;
						snpvalues[i] = snpvalues[i] * parents[popIndex[i]];
					}

					//build and solve the model
					CovariateModelEffect cme = new CovariateModelEffect(snpvalues);
					cme.setId(new Integer(pos));
					double[] Fp =  lmsr.testNewEffect(cme);

					if (!Double.isNaN(Fp[1]) && (Fp[1] < minp || (Fp[1] == minp && Fp[0] > maxF) )) {
						minp = Fp[1];
						maxF = Fp[0];
						bestEffect = cme;
						bestpos = pos;
					}
				}
				
				if (minp < enterLimit) {
					lmsr.addEffect(bestEffect);
					StringBuilder sb = new StringBuilder();
					sb.append(chromosome).append("\t").append(bestpos);
					sb.append("\t").append(theAGPMap.getCmFromPosition(chromosome, bestpos));
					sb.append("\t").append(minp);
					sb.append("\t").append(iter + startIteration);
					Utilities.flushToOutput(sb.toString(), bwstep);
					System.out.println(sb.toString());
				}
				
			} while (minp < enterLimit);

		}
	}
	
	public void testSnpsWithoutReplacement() {
		//create an index array into genotypes and a population map with a list of samples for each population
		createSampleIndex();
		createPopMap();
		Set<Integer> popset = popMap.keySet();
		int npops = popset.size();

		SnpData snpdata;
		if (files == null) snpdata = new SnpDataImputed(chromosome); 
		else snpdata = new SnpDataImputed(chromosome, files);
		totalSnps = snpdata.getNumberOfSnps();
		int maxmarker = genotypes[0].length - 1;
		int nSamples = residualSamples.length;
		int nSubsamples = npops * 150; //populations by 150 samples
		MersenneTwister mt = new MersenneTwister(new java.util.Date());
		
		//create effects, 1 for mean, npops for populations, 150/population
		Integer[] populations = new Integer[nSubsamples];
		int count = 0;
		for (Integer pop : popset) {
			for (int j = 0; j < 150; j++) populations[count++] = pop;
		}
		ModelEffect meanme = new ModelEffect(new int[nSubsamples]);
		ModelEffect popme = new ModelEffect(ModelEffect.getIntegerLevels(populations));
		
		//open file for output
		BufferedWriter bwstep;
		if (files == null) {
			StringBuilder sb = new StringBuilder("C:/Projects/NAM/leaf traits/results/stepwise/stepwise_chr");
			sb.append(chromosome).append("_");
			sb.append(thisTraitName).append(".txt");
			bwstep = Utilities.openOutputFile(sb.toString());
		}
		else bwstep = Utilities.openOutputFile(files.chrsteps[chromosome - 1]);
		Utilities.flushToOutput("chromosome\tposition\tcM\tallele\tp-value,\titeration", bwstep);
		
		//create model files
		BufferedWriter bw;
		if (files == null) {
			StringBuilder sb = new StringBuilder("C:/Projects/NAM/leaf traits/results/stepwise/");
			sb.append("effects_").append(thisTraitName);
			sb.append("_chr").append(chromosome).append(".txt");
			modelOutputFile = new File(sb.toString());
		}
		else modelOutputFile = files.chrmodel[chromosome - 1];
		bw = Utilities.openOutputFile(modelOutputFile);
		Utilities.writeToOutput("chromosome\tposition\tcM\tallele\teffect\tpvalue\titeration", bw);
		Utilities.closeOutput(bw);
		
		for (int iter = 0; iter < iterations; iter++) {
			//get data for this chromosome 
		    System.out.println("Starting iteration " + iter + "...");
			double[] y = new double[nSubsamples];
			int[] bootstrapGenotypeIndex = new int[nSubsamples];
			int[] bootstrapSampleIndex = new int[nSubsamples];
			count = 0;
			long[] sample = new long[150];
			for (Integer pop : popset) {
				ArrayList<Integer> popList = popMap.get(pop);
				int N = popList.size();
				if (N >= 150) RandomSampler.sample(150, N, 150, 0, sample, 0, mt);
				//where population size < 150, a sample of size 0.8 * the population size is created without replacement
				//the rest of the 150 is filled by resampling the selected 80%
				else {
					int subsize = N * 8 / 10;
					RandomSampler.sample(subsize, N, subsize, 0, sample, 0, mt);
					for (int i = subsize; i < 150; i++) {
						int ndx = Uniform.staticNextIntFromTo(0, subsize - 1);
						sample[i] = sample[ndx];
					}
				}
				int[] shuffledIndex = null;
				if (permute) shuffledIndex = getShuffler().getNextPermutationIndex();
				for (int i = 0; i < 150; i++) {
					int ndx = popList.get((int) sample[i]);
					if (permute) {
						int ndx2 = popList.get((int) sample[shuffledIndex[i]]);
						y[count] = residuals[ndx2][chromosome - 1];
					}
					else y[count] = residuals[ndx][chromosome - 1];
					bootstrapGenotypeIndex[count] = sampleIndex[ndx];
					bootstrapSampleIndex[count] = ndx;
					count++;
				}
				
			}
			
			//create the base model
			ArrayList<ModelEffect> modelEffects = new ArrayList<ModelEffect>();
			modelEffects.add(meanme);
			modelEffects.add(popme);
			
			LinearModelforStepwiseRegression lmsr = new LinearModelforStepwiseRegression(modelEffects, y);
			
			//forward regression
			double minp;
			do {
				//create the base model
				minp = 1;
				int bestpos = -1;
				String bestallele = "";
				ModelEffect bestEffect = null;
				double maxF = 0;

				snpcount = 0;
				snpdata.reset();
				while (snpdata.next()) {
					snpcount++;
					double[] parents = snpdata.getGenotype();
					int pos = snpdata.getPosition();
					String allele = snpdata.getAllele();
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

					double[] snpvalues = new double[nSubsamples];
					for (int i = 0; i < nSubsamples; i++) {
						snpvalues[i] = genotypes[bootstrapGenotypeIndex[i]][leftmarker] * (1 - pd) + genotypes[bootstrapGenotypeIndex[i]][rightmarker] * pd;
						snpvalues[i] = snpvalues[i] * parents[populationIndex[bootstrapSampleIndex[i]]];
					}

					//build and solve the model
					CovariateModelEffect cme = new CovariateModelEffect(snpvalues);
					cme.setId(new Object[]{new Integer(pos), allele});
					double[] Fp =  lmsr.testNewEffect(cme);

					if (!Double.isNaN(Fp[1]) && (Fp[1] < minp || (Fp[1] == minp && Fp[0] > maxF) )) {
						minp = Fp[1];
						maxF = Fp[0];
						bestEffect = cme;
						bestpos = pos;
						bestallele = allele;
					}
				}
				
				if (minp < enterLimit) {
					lmsr.addEffect(bestEffect);
					StringBuilder sb = new StringBuilder();
					sb.append(chromosome).append("\t").append(bestpos);
					sb.append("\t").append(theAGPMap.getCmFromPosition(chromosome, bestpos));
					sb.append("\t").append(bestallele);
					sb.append("\t").append(minp);
					sb.append("\t").append(iter + startIteration);
					Utilities.flushToOutput(sb.toString(), bwstep);
					System.out.println(sb.toString());
				}
				
			} while (minp < enterLimit && modelEffects.size() < maxsnps + 2);

			//code for writing to file estimates from the full model
			//write locus effect estimates to a file
			DoubleMatrix1D beta = lmsr.getLinearModel().getBeta();
			int nbeta = beta.size();
			int neffects = modelEffects.size();
			int nloci = neffects - 2;
			int firsteffect = nbeta - nloci;
			
			bw = Utilities.openOutputFile(modelOutputFile, true);
			
			double errSS = lmsr.getLinearModel().getErrorSS();
			double errdf = lmsr.getLinearModel().getErrordf();
			
			for (int i = 0; i < nloci; i++) {
				Object[] info = (Object[]) modelEffects.get(2 + i).getId();
				Integer pos = (Integer) info[0];
				String allele = (String) info[1];
				double cM = theAGPMap.getCmFromPosition(chromosome, pos);
				double estEffect = beta.getQuick(i + firsteffect);
				double[] effectSSdf = lmsr.getLinearModel().marginalEffectSSdf(2 + i);
				double F = effectSSdf[0] / effectSSdf[1] / errSS * errdf;
				double p;
				try {p = AbstractLinearModel.Ftest(F, effectSSdf[1], errdf);}
				catch (Exception e) {p = Double.NaN;}
				Utilities.writeToOutput(chromosome + "\t" + pos + "\t" + cM + "\t" + allele + "\t" + estEffect + "\t" + p+ "\t" + (iter + startIteration), bw);
			}
			Utilities.closeOutput(bw);
			
		}
	}
	
	public synchronized void createSampleIndex() {
		if (sampleIndex == null) {
			int nSamples = residualSamples.length;
			sampleIndex = new int[nSamples];
			for (int i = 0; i < nSamples; i++) sampleIndex[i] = genotypeMap.get(residualSamples[i]);
		}
	}
	
	public synchronized void createPopMap() {
		if (popMap == null) {
			popMap = new HashMap<Integer, ArrayList<Integer>>();
			int nSamples = residualSamples.length;
			populationIndex = new int[nSamples];
			for (int i = 0; i < nSamples; i++) {
				int ipop = getPopulation(residualSamples[i]) - 1;
				ArrayList<Integer> popList = popMap.get(ipop);
				if (popList == null) {
					popList = new ArrayList<Integer>();
					popMap.put(ipop, popList);
				}
				popList.add(i);
				populationIndex[i] = ipop;
			}
		}
	}

	public void setEnterLimit(double enterLimit) {
		this.enterLimit = enterLimit;
	}
	
	public Shuffler getShuffler() {
		if (shuffler == null) shuffler = new SimpleShuffler(150);
		return shuffler;
	}
}
