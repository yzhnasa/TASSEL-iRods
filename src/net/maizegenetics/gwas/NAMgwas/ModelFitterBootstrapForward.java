package net.maizegenetics.gwas.NAMgwas;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;
import java.util.TreeMap;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import net.maizegenetics.jGLiM.LinearModelUtils;
import net.maizegenetics.jGLiM.dm.CovariateModelEffect;
import net.maizegenetics.jGLiM.dm.FactorModelEffect;
import net.maizegenetics.jGLiM.dm.LinearModelForStepwiseRegression;
import net.maizegenetics.jGLiM.dm.ModelEffect;
import net.maizegenetics.jGLiM.dm.ModelEffectUtils;

public class ModelFitterBootstrapForward extends ModelFitter {
	static Random rng = new Random();
	ArrayList<int[]> populationArrays;
	int totalSubsampleSize;
	int[] popSubsampleSize;
	int[] subsample;
	double enterLimit = 1e-6;
	int maxThreads = 50;
	int numberOfThreads;
	int nSamples;
	
	public ModelFitterBootstrapForward(int chromosome, FileNames files, boolean threaded) {
		super(chromosome, files);
		if (!Double.isNaN(files.enterlimit)) enterLimit = files.enterlimit;
	}

	public ModelFitterBootstrapForward(int chromosome, FileNames files) {
		super(chromosome, files);
		if (!Double.isNaN(files.enterlimit)) enterLimit = files.enterlimit;
	}

	@Override
	protected void testSnps() {
		long start = System.currentTimeMillis();
		System.out.println("Bootstrap analysis starting at " + System.currentTimeMillis());
		buildPopulationList();
		
		//initialize output
		try {
			initializeOutput();
		} catch (IOException e) {
			e.printStackTrace();
			return;
		}
		
		if (files.randomizeSnpOrder) {
			randomSnpThreadedAnalysis();
		} else if (files.threaded) {
			numberOfThreads = Runtime.getRuntime().availableProcessors();
			threadedAnalysis();
		} else {
			nonThreadedAnalysis();
		}
		
		System.out.println("Bootstrap analysis finished, elapsed time = " + (System.currentTimeMillis() - start));
	}
	
	protected void simpleAnalysis() {
		//analyze bootstrap (subbagging) samples
		for (int i = 0; i < files.iterations; i++) {
			getSubsample();
			
			//fit forward regression
			LinearModelForStepwiseRegression lmsr = getBaseModel();
			SnpInfo nextSnp = findNextTerm(lmsr);
			while (nextSnp.p < enterLimit) {
				lmsr.addEffect(new CovariateModelEffect(nextSnp.genotype, nextSnp));
				nextSnp = findNextTerm(lmsr);
			}
			
			//sent results to output
			try {
				recordResults(i + files.startIteration, lmsr);
			} catch (IOException e) {
				String msg = "Failed to write results for iteration " + i;
				e.printStackTrace();
			}
		}
	}
	
	protected void nonThreadedAnalysis() {
		int numberOfIterations = files.iterations;
		
//		ForwardRegressionModel[] theModels = new ForwardRegressionModel[numberOfIterations];
		ArrayList<ForwardRegressionModel> theModels = new ArrayList<ForwardRegressionModel>();
		for (int iter = 0; iter < numberOfIterations; iter++) {
			//get data for this chromosome 
			double[] data = new double[totalSubsampleSize];
			getSubsample();
			for (int i = 0; i < totalSubsampleSize; i++) data[i] = residuals[subsample[i]][chromosome - 1];
			
			theModels.add(new ForwardRegressionModelSubsample(getBaseEffects(), data, enterLimit, subsample));
		}

		boolean areAnyUpdateable = true;
		
		while (areAnyUpdateable) {
			snpdata.reset();
			boolean hasnext = snpdata.next();
			ArrayList<ForwardRegressionModel> updateableModels = new ArrayList<ForwardRegressionModel>();
			for (ForwardRegressionModel frm : theModels) {
				if (frm.wasUpdated) updateableModels.add(frm);
			}
			if (updateableModels.size() > 0) areAnyUpdateable = true;
			else areAnyUpdateable = false;
			
			while (hasnext) {
				// create a list of snps to test
				LinkedList<SnpInfo> snpList = new LinkedList<SnpInfo>();
				while (hasnext && snpList.size() <=500) {
					double[] snp = projectSnp(snpdata.getGenotype(), snpdata.getPosition(), popIndex);
					snpList.add(new SnpInfo(chromosome, snpdata.getPosition(),snpdata.getAllele(), sampleArray(snp), 1, 1));
					hasnext = snpdata.next();
				}

				//create and fire the updater threads
				ForwardRegressionUpdater theUpdater;
				theUpdater = new ForwardRegressionUpdater(updateableModels, snpList, !hasnext);
				theUpdater.run();
			}
		}
		
		//write the output
		int iterationCount = files.startIteration;
		for (ForwardRegressionModel frm : theModels) {
			
			try {
				recordResults(iterationCount++, frm.lmsr);
			} catch (IOException e) {
				System.err.println("Error recording results for iteration " + (iterationCount - 1));
				e.printStackTrace();
			}
		}
	}
	
	protected void threadedAnalysis() {
		int numberOfIterations = files.iterations;
		ForwardRegressionModel[] theModels = new ForwardRegressionModel[numberOfIterations];
		ArrayList<ArrayList<ForwardRegressionModel>> theModelLists = new ArrayList<ArrayList<ForwardRegressionModel>>();
		for (int iter = 0; iter < numberOfIterations; iter++) {
			//get data for this chromosome 
			double[] data = new double[totalSubsampleSize];
			getSubsample();
			for (int i = 0; i < totalSubsampleSize; i++) data[i] = residuals[subsample[i]][chromosome - 1];
			
			theModels[iter] = new ForwardRegressionModelSubsample(getBaseEffects(), data, enterLimit, subsample);
		}


		boolean anyUpdateable = true;
		//create lists of models which need to have terms added
		//the number of lists equals the number of threads to be run
		//add snps then update the models; wash, rinse, repeat until all models are fit
		
		while (anyUpdateable) {
			System.out.println("Adding a term to the models.");
			long start = System.currentTimeMillis();
			theModelLists.clear();
			for (int i = 0; i < numberOfThreads; i++) {
				theModelLists.add(new ArrayList<ForwardRegressionModel>());
			}

			//create lists of forward regression models
			int currentList = 0;
			anyUpdateable = false;
			for (ForwardRegressionModel frm : theModels) {
				if (frm.wasUpdated) {
					theModelLists.get(currentList++).add(frm);
					if (currentList == numberOfThreads) currentList = 0;
					anyUpdateable = true;
				}
			}

			if (anyUpdateable) {
				snpdata.reset();
				boolean hasnext = snpdata.next();
				
				while (hasnext) {
					// create a list of snps to test
					LinkedList<SnpInfo> snpList = new LinkedList<SnpInfo>();
					while (hasnext && snpList.size() <=500) {
						double[] snp = projectSnp(snpdata.getGenotype(), snpdata.getPosition(), popIndex);
						snpList.add(new SnpInfo(chromosome, snpdata.getPosition(),snpdata.getAllele(), snp, 1, 1));
						hasnext = snpdata.next();
					}
					
					//create and fire the updater threads
					ForwardRegressionUpdater[] updaters = null;
					updaters = new ForwardRegressionUpdater[numberOfThreads];
					for (int i = 0; i < numberOfThreads; i++) {
						updaters[i] = new ForwardRegressionUpdater(theModelLists.get(i), snpList, !hasnext);
						updaters[i].start();
					}
					
					//wait for the threads to finish
					for (int i = 0; i < numberOfThreads; i++) {
						try {
							updaters[i].join();
						} catch (InterruptedException e) {
						}
					}
				}
			}
			System.out.println("added a term in " + (System.currentTimeMillis() - start));
		}
		//write the output
		for (int i = 0; i < numberOfIterations; i++) {
			int thisIteration = files.startIteration + i;
			try {
				recordResults(thisIteration, theModels[i].lmsr);
			} catch (IOException e) {
				System.err.println("Error recording results for iteration " + thisIteration);
				e.printStackTrace();
			}
		}
	}
	
	protected void randomSnpThreadedAnalysis() {
		//Use a thread pool. For each thread:
		//Instantiate a forward regression model. Add that and a snpdata to a thread.
		//The thread should cycle through the snps several times until no more snps are added 
		//to the model. The finished model is a completed iteration, which is recorded.
		//The iterations can be numbered and written to output in the order finished.
		
		//The thread pool
		int numberOfThreads = Runtime.getRuntime().availableProcessors();
		ExecutorService executor = Executors.newFixedThreadPool(numberOfThreads);
		
		//The loop for iterating through bootstraps
		LinkedList<ForwardRegressionFitter> modelList = new LinkedList<ForwardRegressionFitter>(); 
		for (int iter = 0; iter < files.iterations; iter++) {
			double[] data = new double[totalSubsampleSize];
			Integer[] pops = new Integer[totalSubsampleSize];
			int[] mean = new int[totalSubsampleSize];
			int[] sampleIndex = new int[totalSubsampleSize];
			getSubsample();
			for (int i = 0; i < totalSubsampleSize; i++) {
				data[i] = residuals[subsample[i]][chromosome - 1];
				pops[i] = new Integer(popIndex[subsample[i]]);
				sampleIndex[i] = subsample[i]; 
			}
			ArrayList<ModelEffect> initialEffects = new ArrayList<ModelEffect>();
			initialEffects.add(new FactorModelEffect(mean, false));
			initialEffects.add(new FactorModelEffect(ModelEffectUtils.getIntegerLevels(pops, null), true));
			ForwardRegressionModel frm = new ForwardRegressionModelSubsample(initialEffects, data, enterLimit, sampleIndex);
			modelList.add(new ForwardRegressionFitter(this, frm, snpdata.getCopy()));
		}
		System.out.println("All models have been created. Fitting will proceed.");
		List<Future<ForwardRegressionModel>> futureList = new LinkedList<Future<ForwardRegressionModel>>();
		try {
			futureList = executor.invokeAll(modelList);
		} catch (InterruptedException e1) {
			e1.printStackTrace();
		}
		executor.shutdown();
		
		//Output the results to a file
		int modelCount = files.startIteration;
		for (Future<ForwardRegressionModel> f : futureList) {
			try {
				ForwardRegressionModel frm = f.get();
				recordResults(modelCount, frm.lmsr);
				
				System.out.println("Results recorded for iteration " + modelCount);
			} catch (InterruptedException e) {
				System.err.println("Iteration " + modelCount + " failed.");
				e.printStackTrace();
			} catch (ExecutionException e) {
				System.err.println("Iteration " + modelCount + " failed.");
				e.printStackTrace();
			} catch (IOException e) {
				System.err.println("Error writing results for iteration " + modelCount);
				e.printStackTrace();
			}
			modelCount++;
		}

	}
	
	//for the next term (max modelss) need pos, allele, projection, F, p
	protected SnpInfo findNextTerm(LinearModelForStepwiseRegression lmsr) {
//		snpcount = 0;
		SnpInfo bestSnp = null;
		double bestSS = 0;
		snpdata.reset();
		while (snpdata.next()) {
//			snpcount++;
//			if (snpcount % 10000 == 0) System.out.println("Chromosome " + chromosome + ", snp " + snpcount);
			double[] parents = snpdata.getGenotype();
			int pos = snpdata.getPosition();
			
			//build and solve the model
			double[] snp = projectSnp(parents, pos, popIndex);
			double ms =  lmsr.testNewEffect(snp);
			
			if (ms > bestSS) {
				bestSS = ms;
				double[] Fp = lmsr.getFpFromModelSS(ms);
				bestSnp = new SnpInfo(chromosome, pos, snpdata.getAllele(), snp, Fp[0], Fp[1]);
			}
		}
		System.out.println(bestSnp.pos + ", " + bestSnp.allele + ", " + bestSnp.F + ", " + bestSnp.p);
		return bestSnp;
	}

	public static synchronized void shuffle(int[] array) {
	    // i is the number of items remaining to be shuffled.
		int n = array.length;
	    for (int i = n; i > 1; i--) {
	        // Pick a random element to swap with the i-th element.
	        int j = rng.nextInt(i);  // 0 <= j <= i-1 (0-based array)
	        // Swap array elements.
	        int tmp = array[j];
	        array[j] = array[i-1];
	        array[i-1] = tmp;
	    }
	}

	private void buildPopulationList() {
		TreeMap<Integer, LinkedList<Integer>> popmap = new TreeMap<Integer, LinkedList<Integer>>();
		int count = 0;
		for (int pop : popIndex) {
			LinkedList<Integer> poplist = popmap.get(pop);
			if (poplist == null) {
				poplist = new LinkedList<Integer>();
				popmap.put(pop, poplist);
			}
			poplist.add(count++);
		}
		
		populationArrays = new ArrayList<int[]>();
		int n = popmap.size();
		totalSubsampleSize = 0;
		popSubsampleSize = new int[n];
		
		Iterator<LinkedList<Integer>> mit = popmap.values().iterator();
		int pcount = 0;
		while (mit.hasNext()) {
			LinkedList<Integer> plist = mit.next();
			int m = plist.size();
			int[] samples = new int[m];
			count = 0;
			for (Integer sample : plist) {
				samples[count++] = sample.intValue();
			}
			populationArrays.add(samples);
			
			popSubsampleSize[pcount] = m * 8 / 10;
//			popSubsampleSize[i] = m; //debug check
			totalSubsampleSize += popSubsampleSize[pcount];
			pcount++;

		}
		
	}
	
	private void getSubsample() {
		int start = 0;
		int npops = populationArrays.size();
		subsample = new int[totalSubsampleSize];
		
		for (int i = 0; i < npops; i++) {
			int[] samples = populationArrays.get(i);
			shuffle(samples);
			System.arraycopy(samples, 0, subsample, start, popSubsampleSize[i]);
			start += popSubsampleSize[i];
		}
	}

	private double[] sampleArray(double[] in) {
		int nIn = in.length;
		int nOut = totalSubsampleSize;
		double[] out = new double[nOut];
		for (int i = 0; i < nOut; i++) {
			out[i] = in[subsample[i]];
		}
		return out;
	}
	
	protected LinearModelForStepwiseRegression getBaseModel() {
		ArrayList<ModelEffect> effects = getBaseEffects();
		
		//get data for this chromosome 
		double[] data = new double[totalSubsampleSize];
		for (int i = 0; i < totalSubsampleSize; i++) data[i] = residuals[subsample[i]][chromosome - 1];
		
		//create the base model
		return new LinearModelForStepwiseRegression(effects, data);
	}

	protected ArrayList<ModelEffect> getBaseEffects() {
		int[] mean = new int[totalSubsampleSize];

		//set up the mean and pop model effects
		ArrayList<ModelEffect> effects = new ArrayList<ModelEffect>();
		ModelEffect memean = new FactorModelEffect(mean, false);
		effects.add(memean);
		ArrayList<Integer> populations = new ArrayList<Integer>();
		for (int sample : subsample) populations.add(popIndex[sample]);
		ModelEffect mepop = new FactorModelEffect(ModelEffectUtils.getIntegerLevels(populations), true);
		effects.add(mepop);
		
		return effects;
	}
	
	private void initializeOutput() throws IOException {
		BufferedWriter bw = new BufferedWriter(new FileWriter(files.chrmodel[chromosome - 1]));
		bw.write("chromosome\tposition\tcM\tallele\teffect\tF\tpvalue\titeration");
		bw.newLine();
		bw.close();
	}
	
	private void recordResults(int iteration, LinearModelForStepwiseRegression lmsr) throws IOException {
		
		String tab = "\t";
		ArrayList<ModelEffect> effects = lmsr.getModelEffects();
		int nsnps = effects.size() - 2;
		double[] beta = lmsr.getLinearModel().getBeta();
		int nbeta = beta.length;
		int start = nbeta - nsnps;
		
		BufferedWriter bw = new BufferedWriter(new FileWriter(files.chrmodel[chromosome - 1], true));
		double[] errorssdf = lmsr.getLinearModel().getResidualSSdf();
		for (int i = 0; i < nsnps; i++) {
			ModelEffect thiseffect = effects.get(i + 2);
			SnpInfo snpinfo = (SnpInfo) thiseffect.getID();
			
			//calculate F and p
			double[] snpssdf = lmsr.getLinearModel().getMarginalSSdf(i + 2);
			double F = snpssdf[0] / snpssdf[1] / errorssdf[0] * errorssdf[1];
			double p;
			try {
				p = LinearModelUtils.Ftest(F, snpssdf[1], errorssdf[1]);
			} catch (Exception e) {
				p = Double.NaN;
			}
			
			StringBuilder sb = new StringBuilder();
			sb.append(chromosome);
			sb.append(tab).append(snpinfo.pos);
			sb.append(tab).append(theAGPMap.getCmFromPosition(chromosome, snpinfo.pos));
			sb.append(tab).append(snpinfo.allele);
			sb.append(tab).append(beta[start + i]);
			sb.append(tab).append(F);
			sb.append(tab).append(p);
			sb.append(tab).append(iteration);
			bw.write(sb.toString());
			bw.newLine();
		}
		bw.close();
	}

}
