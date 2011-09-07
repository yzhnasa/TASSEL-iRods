package net.maizegenetics.gwas.LeafTraits;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.LinkedList;

import net.maizegenetics.jGLiM.AbstractLinearModel;
import net.maizegenetics.jGLiM.CovariateModelEffect;
import net.maizegenetics.jGLiM.LinearModelforStepwiseRegression;
import net.maizegenetics.jGLiM.ModelEffect;
import net.maizegenetics.jGLiM.NestedCovariateModelEffect;
import cern.colt.matrix.DoubleMatrix1D;

public class ModelFitterStepwiseWithCofactors extends ModelFitter {

	double enterLimit = 1e-6;
	SnpData snpdata;
	NamMap namMap = new NamMap();
	int[] sampleIndex;
	int maxmarker;
	int step = 0;
	ArrayList<Integer> qtlList;
	Integer current_qtl = 0;
	
	public ModelFitterStepwiseWithCofactors(int chromosome, int trait) {
		super(chromosome, trait);
	}
	
	public void testSnps() {
		testSnpsByWindow();
	}
	
	public void testSnpsA() {
		snpdata = new SnpDataImputed(chromosome);
		totalSnps = snpdata.getNumberOfSnps();
		maxmarker = genotypes[0].length - 1;
		int nSamples = residualSamples.length;

		qtlList = new ArrayList<Integer>();
		qtlList.add(-1);
		qtlList.addAll(getQTLforChromosome());
		
		int[] popIndex = new int[nSamples];
		ArrayList<String> pops = new ArrayList<String>();
		int count = 0;
		for (String sample : residualSamples) {
			int ipop = getPopulation(sample) - 1;
			pops.add(SnpData.popnames[ipop]);
			popIndex[count++] = ipop;
		}

		//create an index array into genotypes
		sampleIndex = new int[nSamples];
		for (int i = 0; i < nSamples; i++) sampleIndex[i] = genotypeMap.get(residualSamples[i]);

		//open files for output
		//the file to record steps
		StringBuilder sb = new StringBuilder("C:/Projects/NAM/leaf traits/results/stepwise/stepwiseByQTL_");
		sb.append(thisTraitName);
		sb.append("_chr").append(chromosome).append(".txt");
		BufferedWriter bwstep = Utilities.openOutputFile(sb.toString());
		Utilities.flushToOutput("chromosome\tqtl\tposition\tp-value\teffect", bwstep);
		
		//the file to record final effect estimates
		sb = new StringBuilder("C:/Projects/NAM/leaf traits/results/stepwise/");
		sb.append("effectsByQTL_").append(thisTraitName);
		sb.append("_chr").append(chromosome).append(".txt");
		BufferedWriter bw = Utilities.openOutputFile(sb.toString());
		Utilities.flushToOutput("chromosome\tqtl\tposition\teffect\tpvalue", bw);

		for (Integer whichqtl : qtlList) {
			current_qtl = whichqtl;
			LinearModelforStepwiseRegression lmsr = getBaseModelWithQTL(pops, whichqtl);
			ArrayList<ModelEffect> effects = lmsr.getModelEffects();

			//forward regression
			step = 1;
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

					//build and solve the model
					CovariateModelEffect cme = new CovariateModelEffect(projectSnp(parents, pos, popIndex));
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
					snpdata.reset();
					lmsr.addEffect(bestEffect);
					
					//get the effect size
					DoubleMatrix1D beta = lmsr.getLinearModel().getBeta();
					double estEffect = beta.getQuick(beta.size() - 1);
					sb = new StringBuilder();
					sb.append(chromosome).append("\t").append(whichqtl).append("\t").append(bestpos);
					sb.append("\t").append(minp).append("\t").append(estEffect);
					Utilities.flushToOutput(sb.toString(), bwstep);
					System.out.println(sb.toString());
				}
				
				
				
				step++;
			} while (minp < enterLimit);
			step -= 2;
			
			//report effect sizes and p-values in the final model;
			DoubleMatrix1D beta = lmsr.getLinearModel().getBeta();
			int nbeta = beta.size();
			int neffects = effects.size();
			double sserr = lmsr.getLinearModel().getErrorSS();
			double dferr = lmsr.getLinearModel().getErrordf();
			
			for (int i = 0; i < step; i++) {
				int effectnumber = neffects - step + i;
				Object snpid = effects.get(effectnumber).getId();
				double estEffect = beta.getQuick(nbeta - step + i);
				double[] ssdf = lmsr.getLinearModel().marginalEffectSSdf(effectnumber);
				double F = ssdf[0] / ssdf[1] / sserr * dferr;
				double p;
				try {p = AbstractLinearModel.Ftest(F, ssdf[1], dferr);}
				catch (Exception e) {p = Double.NaN;}

				sb = new StringBuilder();
				sb.append(chromosome).append("\t").append(whichqtl).append("\t");
				sb.append(snpid.toString()).append("\t").append(estEffect).append("\t").append(p);
				Utilities.flushToOutput(sb.toString(), bw);
			}
		}
		
		Utilities.closeOutput(bwstep);
		Utilities.closeOutput(bw);
		System.out.println("Finished chromosome " + chromosome + ".");
	}
	
	public void testSnpsB() {
		snpdata = new SnpDataImputed(chromosome);
		totalSnps = snpdata.getNumberOfSnps();
		maxmarker = genotypes[0].length - 1;
		int nSamples = residualSamples.length;

		qtlList = getQTLforChromosome();
		
		int[] popIndex = new int[nSamples];
		ArrayList<String> pops = new ArrayList<String>();
		int count = 0;
		for (String sample : residualSamples) {
			int ipop = getPopulation(sample) - 1;
			pops.add(SnpData.popnames[ipop]);
			popIndex[count++] = ipop;
		}

		//create an index array into genotypes
		sampleIndex = new int[nSamples];
		for (int i = 0; i < nSamples; i++) sampleIndex[i] = genotypeMap.get(residualSamples[i]);

		//open files for output
		//the file to record steps
		StringBuilder sb = new StringBuilder("C:/Projects/NAM/leaf traits/results/stepwise/stepwiseByQTLv2_");
		sb.append(thisTraitName);
		sb.append("_chr").append(chromosome).append(".txt");
		BufferedWriter bwstep = Utilities.openOutputFile(sb.toString());
		Utilities.flushToOutput("chromosome\tqtl\tposition\tp-value\teffect", bwstep);
		
		//the file to record final effect estimates
		sb = new StringBuilder("C:/Projects/NAM/leaf traits/results/stepwise/");
		sb.append("effectsByQTLv2_").append(thisTraitName);
		sb.append("_chr").append(chromosome).append(".txt");
		BufferedWriter bw = Utilities.openOutputFile(sb.toString());
		Utilities.flushToOutput("chromosome\tqtl\tposition\teffect\tpvalue", bw);

		current_qtl = -1;
		LinearModelforStepwiseRegression lmsr = getBaseModelWithQTL(pops, current_qtl);
		ArrayList<ModelEffect> effects = lmsr.getModelEffects();
		ArrayList<ModelEffect> snpEffects = new ArrayList<ModelEffect>();
				
		//forward regression
		step = 1;
		double minp;
		do {
			//create the base model
			minp = forwardStep(lmsr, popIndex, bwstep);
			if (minp < enterLimit) {
				int neffects = effects.size();
				snpEffects.add(effects.get(neffects - 1));
			}
			step++;
		} while (minp < enterLimit);
		step -= 2;
		
		//report effect sizes and p-values in the final model;
		DoubleMatrix1D beta = lmsr.getLinearModel().getBeta();
		int nbeta = beta.size();
		int neffects = effects.size();
		double sserr = lmsr.getLinearModel().getErrorSS();
		double dferr = lmsr.getLinearModel().getErrordf();
		
		for (int i = 0; i < step; i++) {
			int effectnumber = neffects - step + i;
			Object snpid = effects.get(effectnumber).getId();
			double estEffect = beta.getQuick(nbeta - step + i);
			double[] ssdf = lmsr.getLinearModel().marginalEffectSSdf(effectnumber);
			double F = ssdf[0] / ssdf[1] / sserr * dferr;
			double p;
			try {p = AbstractLinearModel.Ftest(F, ssdf[1], dferr);}
			catch (Exception e) {p = Double.NaN;}

			sb = new StringBuilder();
			sb.append(chromosome).append("\t").append(current_qtl).append("\t");
			sb.append(snpid.toString()).append("\t").append(estEffect).append("\t").append(p);
			Utilities.flushToOutput(sb.toString(), bw);
		}

		
		//get data for this chromosome 
		double[] y = new double[nSamples];
		for (int i = 0; i < nSamples; i++) y[i] = residuals[i][chromosome - 1];

		for (Integer whichqtl : qtlList) {
			current_qtl = whichqtl;
			
			lmsr = getBaseModelWithQTL(pops, current_qtl);
			for (ModelEffect me : snpEffects) lmsr.addEffect(me);
			
			//forward regression
			step = 1;
			do {
				//create the base model
				minp = forwardStep(lmsr, popIndex, bwstep);
				step++;
			} while (minp < enterLimit);
			step -= 2;
			
			//report effect sizes and p-values in the final model;
			beta = lmsr.getLinearModel().getBeta();
			nbeta = beta.size();
			neffects = effects.size();
			sserr = lmsr.getLinearModel().getErrorSS();
			dferr = lmsr.getLinearModel().getErrordf();
			
			for (int i = 0; i < step; i++) {
				int effectnumber = neffects - step + i;
				Object snpid = effects.get(effectnumber).getId();
				double estEffect = beta.getQuick(nbeta - step + i);
				double[] ssdf = lmsr.getLinearModel().marginalEffectSSdf(effectnumber);
				double F = ssdf[0] / ssdf[1] / sserr * dferr;
				double p;
				try {p = AbstractLinearModel.Ftest(F, ssdf[1], dferr);}
				catch (Exception e) {p = Double.NaN;}

				sb = new StringBuilder();
				sb.append(chromosome).append("\t").append(whichqtl).append("\t");
				sb.append(snpid.toString()).append("\t").append(estEffect).append("\t").append(p);
				Utilities.flushToOutput(sb.toString(), bw);
			}
		}
		
		Utilities.closeOutput(bwstep);
		Utilities.closeOutput(bw);
		System.out.println("Finished chromosome " + chromosome + ".");
	}
	
	public void testSnpsByWindow() {
		double window = 10; //window is +/-window cM
		snpdata = new SnpDataImputed(chromosome);
		totalSnps = snpdata.getNumberOfSnps();
		maxmarker = genotypes[0].length - 1;
		int nSamples = residualSamples.length;

		qtlList = getQTLforChromosome();
		
		int[] popIndex = new int[nSamples];
		ArrayList<String> pops = new ArrayList<String>();
		int count = 0;
		for (String sample : residualSamples) {
			int ipop = getPopulation(sample) - 1;
			pops.add(SnpData.popnames[ipop]);
			popIndex[count++] = ipop;
		}

		//create an index array into genotypes
		sampleIndex = new int[nSamples];
		for (int i = 0; i < nSamples; i++) sampleIndex[i] = genotypeMap.get(residualSamples[i]);

		//open files for output
		//the file to record steps
//		StringBuilder sb = new StringBuilder("C:/Projects/NAM/leaf traits/results/stepwise/stepwiseByQTLWindow_");
		StringBuilder sb = new StringBuilder("C:/Projects/NAM/flowering/results/stepwiseByQTLWindow_");
		sb.append(thisTraitName);
		sb.append("_chr").append(chromosome).append(".txt");
		BufferedWriter bwstep = Utilities.openOutputFile(sb.toString());
		Utilities.flushToOutput("chromosome\tqtl\tposition\tcM\tallele\tp-value\teffect", bwstep);
		
		//the file to record final effect estimates
//		sb = new StringBuilder("C:/Projects/NAM/leaf traits/results/stepwise/");
		sb = new StringBuilder("C:/Projects/NAM/flowering/results/");
		sb.append("effectsByQTLWindow_").append(thisTraitName);
		sb.append("_chr").append(chromosome).append(".txt");
		BufferedWriter bw = Utilities.openOutputFile(sb.toString());
		Utilities.flushToOutput("chromosome\tqtl\tposition\tcM\tallele\teffect\tpvalue", bw);

		for (Integer whichqtl : qtlList) {
			current_qtl = whichqtl;
			double qtlPosition = namMap.position.get(current_qtl);
			int leftedge = theAGPMap.getPositionFromCm(chromosome, qtlPosition - window);
			int rightedge = theAGPMap.getPositionFromCm(chromosome, qtlPosition + window);
			LinearModelforStepwiseRegression lmsr = getBaseModelWithQTL(pops, whichqtl);
			ArrayList<ModelEffect> effects = lmsr.getModelEffects();

			//forward regression
			step = 1;
			double minp;
			do {
				//create the base model
				minp = 1;
				ModelEffect bestEffect = null;
				double maxF = 0;

				snpcount = 0;
				snpdata.reset();
				
				while (snpdata.next()) {
					if (snpdata.getPosition() >= leftedge && snpdata.getPosition() <= rightedge) {
						snpcount++;
						double[] parents = snpdata.getGenotype();
						int pos = snpdata.getPosition();
						SnpInfo snpinfo = new SnpInfo();
						snpinfo.chromosome = chromosome;
						snpinfo.position = pos;
						snpinfo.allele = snpdata.getAllele();
						snpinfo.scores = snpdata.getGenotype();
						
						//build and solve the model
						CovariateModelEffect cme = new CovariateModelEffect(projectSnp(parents, pos, popIndex));
						cme.setId(snpinfo);
						double[] Fp =  lmsr.testNewEffect(cme);

						if (!Double.isNaN(Fp[1]) && (Fp[1] < minp || (Fp[1] == minp && Fp[0] > maxF) )) {
							minp = Fp[1];
							maxF = Fp[0];
							bestEffect = cme;
						}
					}
				}

				if (minp < enterLimit) {
					lmsr.addEffect(bestEffect);
					
					//get the effect size
					DoubleMatrix1D beta = lmsr.getLinearModel().getBeta();
					double estEffect = beta.getQuick(beta.size() - 1);
					SnpInfo snpinfo = (SnpInfo) bestEffect.getId();
					double cmPosition = theAGPMap.getCmFromPosition(chromosome, snpinfo.position);
					sb = new StringBuilder();
					sb.append(chromosome).append("\t").append(whichqtl).append("\t").append(snpinfo.position);
					sb.append("\t").append(cmPosition);
					sb.append("\t").append(snpinfo.allele);
					sb.append("\t").append(minp).append("\t").append(estEffect);
					Utilities.flushToOutput(sb.toString(), bwstep);
					System.out.println(sb.toString());
				}
				
				step++;
			} while (minp < enterLimit);
			step -= 2;
			
			//report effect sizes and p-values in the final model;
			DoubleMatrix1D beta = lmsr.getLinearModel().getBeta();
			int nbeta = beta.size();
			int neffects = effects.size();
			double sserr = lmsr.getLinearModel().getErrorSS();
			double dferr = lmsr.getLinearModel().getErrordf();
			
			for (int i = 0; i < step; i++) {
				int effectnumber = neffects - step + i;
				SnpInfo snpid = (SnpInfo) effects.get(effectnumber).getId();
				double estEffect = beta.getQuick(nbeta - step + i);
				double[] ssdf = lmsr.getLinearModel().marginalEffectSSdf(effectnumber);
				double F = ssdf[0] / ssdf[1] / sserr * dferr;
				double p;
				try {p = AbstractLinearModel.Ftest(F, ssdf[1], dferr);}
				catch (Exception e) {p = Double.NaN;}

				double cmPosition = theAGPMap.getCmFromPosition(chromosome, snpid.position);
				
				sb = new StringBuilder();
				sb.append(chromosome).append("\t").append(whichqtl);
				sb.append("\t").append(snpid.position);
				sb.append("\t").append(cmPosition);
				sb.append("\t").append(snpid.allele);
				sb.append("\t").append(estEffect).append("\t").append(p);
				Utilities.flushToOutput(sb.toString(), bw);
			}
		}
		
		Utilities.closeOutput(bwstep);
		Utilities.closeOutput(bw);
		System.out.println("Finished chromosome " + chromosome + ".");
	}
	
	public void testMultipleSnpModelByWindow() {
		double window = 10; //window is +/-window cM
		snpdata = new SnpDataImputed(chromosome);
		totalSnps = snpdata.getNumberOfSnps();
		maxmarker = genotypes[0].length - 1;
		int nSamples = residualSamples.length;

		qtlList = getQTLforChromosome();
		
		int[] popIndex = new int[nSamples];
		ArrayList<String> pops = new ArrayList<String>();
		int count = 0;
		for (String sample : residualSamples) {
			int ipop = getPopulation(sample) - 1;
			pops.add(SnpData.popnames[ipop]);
			popIndex[count++] = ipop;
		}

		//create an index array into genotypes
		sampleIndex = new int[nSamples];
		for (int i = 0; i < nSamples; i++) sampleIndex[i] = genotypeMap.get(residualSamples[i]);

		//open files for output
		StringBuilder sb;
		
		current_qtl = 660; //on chromosome 5
		double qtlPosition = namMap.position.get(current_qtl);
		int leftedge = theAGPMap.getPositionFromCm(chromosome, qtlPosition - window);
		int rightedge = theAGPMap.getPositionFromCm(chromosome, qtlPosition + window);
		ArrayList<SnpInfo> snps = getSnpsInWindow(leftedge, rightedge);
		System.out.println("Number of snps = " + snps.size());
		
		snpcount = 0;
		long time = System.currentTimeMillis();
		int nsnps = snps.size();
		double minp = 1;
		double maxF = 0;
		SnpInfo bestsnp;
		CovariateModelEffect[] besteffects = new CovariateModelEffect[2];
		for (int i = 0; i < nsnps; i++) {
			LinearModelforStepwiseRegression lmsr = getBaseModelWithQTL(pops, current_qtl);
			SnpInfo snp1 = snps.get(i);
			CovariateModelEffect cme1 = new CovariateModelEffect(projectSnp(snp1.scores, snp1.position, popIndex));
			lmsr.addEffect(cme1);
			
			for (int j = i + 1; j < nsnps; j++) {
				snpcount++;
				SnpInfo snp2 = snps.get(j);
				CovariateModelEffect cme2 = new CovariateModelEffect(projectSnp(snp1.scores, snp1.position, popIndex));
				double[] Fp = lmsr.testNewEffect(cme2);
				if (Fp[1] < minp || (Fp[1] == minp && Fp[0] > maxF)) {
					minp = Fp[1];
					maxF = Fp[0];
					besteffects[0] = cme1;
					besteffects[1] = cme2;
				}
			}
		}
		
	}
	
	public LinearModelforStepwiseRegression getBaseModel(ArrayList<String> pops) {
		int nSamples = residualSamples.length;
		int[] mean = new int[nSamples];

		//set up the mean and pop model effects
		ArrayList<ModelEffect> effects = new ArrayList<ModelEffect>();
		ModelEffect memean = new ModelEffect(mean);
		effects.add(memean);
		ArrayList<String> popLabels = new ArrayList<String>();
		ModelEffect mepop = new ModelEffect(ModelEffect.getIntegerLevels(pops, popLabels));
		effects.add(mepop);
		
		//get data for this chromosome 
		double[] y = new double[nSamples];
		for (int i = 0; i < nSamples; i++) y[i] = residuals[i][chromosome - 1];
		
		//create the base model
		return new LinearModelforStepwiseRegression(effects, y);
	}
	
	public LinearModelforStepwiseRegression getBaseModelWithQTL(ArrayList<String> pops, int qtlToExclude) {
		int nSamples = residualSamples.length;
		int[] mean = new int[nSamples];

		//set up the mean and pop model effects
		ArrayList<ModelEffect> effects = new ArrayList<ModelEffect>();
		ModelEffect memean = new ModelEffect(mean);
		effects.add(memean);
		ArrayList<String> popLabels = new ArrayList<String>();
		ModelEffect mepop = new ModelEffect(ModelEffect.getIntegerLevels(pops, popLabels));
		effects.add(mepop);
		
		//get data for this chromosome 
		double[] y = new double[nSamples];
		for (int i = 0; i < nSamples; i++) y[i] = residuals[i][chromosome - 1];
		
		//and the QTL are ...
		
		for (Integer marker : qtlList) {
			if (marker >= 0 && marker != qtlToExclude) {
				double[] markerscores = new double[nSamples];
				for (int i = 0; i < nSamples; i++) {
					markerscores[i] = genotypes[sampleIndex[i]][marker - firstMarker + 1];
				}
				CovariateModelEffect cme = new CovariateModelEffect(markerscores);
				ModelEffect me = new NestedCovariateModelEffect(cme, mepop);
				me.setId(marker);
				effects.add(me);
			}
		}
		
		//create the base model
		return new LinearModelforStepwiseRegression(effects, y);
	}
	
	public ArrayList<Integer> getQTLforChromosome() {
//		String dir = "C:/Projects/NAM/leaf traits/models/";
//		String prefix = "StepwiseWithScan082809_leaftrait_distanceImputation_";
//		String suffix = "_model.txt";
//		String filename = dir + prefix + traitnumber + suffix;

		String filename = "C:/Projects/NAM/flowering/stepwise_markers_d2ae-4_model.txt";
		BufferedReader br = Utilities.openReader(filename);
		
		String input;
		ArrayList<Integer> qtlList = new ArrayList<Integer>();
		while ((input = Utilities.readALineFrom(br)) != null) {
			String[] parsed = input.split("\t");
			try {
				Integer marker = Integer.parseInt(parsed[1]);
				if (namMap.chromosome.get(marker) == chromosome) qtlList.add(marker);
			}
			catch(Exception e) {}
		}
		return qtlList;
	}
	
	public double[] projectSnp(double[] parents, int pos, int[] popIndex) {
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
		
		int nSamples = residualSamples.length;
		double[] snpvalues = new double[nSamples];
		for (int i = 0; i < nSamples; i++) {
			snpvalues[i] = genotypes[sampleIndex[i]][leftmarker] * (1 - pd) + genotypes[sampleIndex[i]][rightmarker] * pd;
			snpvalues[i] = snpvalues[i] * parents[popIndex[i]];
		}
		return snpvalues;
	}

	public double forwardStep(LinearModelforStepwiseRegression lmsr, int[] popIndex, BufferedWriter bw) {
		double minp = 1;
		int bestpos = -1;
		ModelEffect bestEffect = null;
		double maxF = 0;

		snpcount = 0;
		snpdata.reset();
		while (snpdata.next()) {
			snpcount++;
			double[] parents = snpdata.getGenotype();
			int pos = snpdata.getPosition();

			//build and solve the model
			CovariateModelEffect cme = new CovariateModelEffect(projectSnp(parents, pos, popIndex));
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
			snpdata.reset();
			lmsr.addEffect(bestEffect);
			
			//get the effect size
			DoubleMatrix1D beta = lmsr.getLinearModel().getBeta();
			double estEffect = beta.getQuick(beta.size() - 1);
			StringBuilder sb = new StringBuilder();
			sb.append(chromosome).append("\t").append(current_qtl).append("\t").append(bestpos);
			sb.append("\t").append(minp).append("\t").append(estEffect);
			Utilities.flushToOutput(sb.toString(), bw);
			System.out.println(sb.toString());
		}
		
		return minp;
	}
	
	public boolean backStep() {
		return false;
	}
	
	public void setEnterLimit(double limit) {enterLimit = limit;}

	public String showProgress() {
		DecimalFormat format = new DecimalFormat("#.#");
		format.setMaximumFractionDigits(2);
//		return "QTL " + current_qtl + " step " + step + ": " + format.format((double)snpcount / (double)snpdata.getNumberOfSnps() * 100) + "%";
		return "Testing snp pair " + snpcount;
	}
	
	public ArrayList<SnpInfo> getSnpsInWindow(int start, int end) {
		ArrayList<SnpInfo> snps = new ArrayList<SnpInfo>();
		snpdata.reset();
		snpdata.next();
		while (snpdata.getPosition() < start) snpdata.next();
		
		while (snpdata.getPosition() < end) {
			SnpInfo snpinfo = new SnpInfo();
			snpinfo.position = snpdata.getPosition();
			snpinfo.chromosome = snpdata.chromosome;
			snpinfo.allele = snpdata.getAllele();
			snpinfo.scores = snpdata.getGenotype();
			snps.add(snpinfo);
			snpdata.next();
		}
		return snps;
	}
	
	class SnpInfo {
		int position;
		int chromosome;
		String allele;
		double[] scores;
		
		public String toString() {
			return Integer.toString(position);
		}
	}
}
