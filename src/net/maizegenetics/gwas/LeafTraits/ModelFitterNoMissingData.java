package net.maizegenetics.gwas.LeafTraits;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;

import net.maizegenetics.jGLiM.AbstractLinearModel;
import net.maizegenetics.jGLiM.CovariateModelEffect;
import net.maizegenetics.jGLiM.LinearModelWithSweep;
import net.maizegenetics.jGLiM.LinearModelforStepwiseRegression;
import net.maizegenetics.jGLiM.ModelEffect;
import net.maizegenetics.jGLiM.NestedCovariateModelEffect;

public class ModelFitterNoMissingData extends ModelFitter {
	
	int[] sampleIndex;
	ArrayList<String> popLabels;
	SnpData snpdata;
	int maxmarker;
	ArrayList<String> pops;
	int[] popIndex;
	ArrayList<Integer> snpList = null;
	ArrayList<ModelEffect> effects;
	
	public ModelFitterNoMissingData(int chromosome, int trait) {
		super(chromosome, trait);
	}
 
	public ModelFitterNoMissingData(int chromosome, FileNames files) {
		super(chromosome, files);
	}
	
	public void testSnps() {
		testSnpsNoInteraction();
	}
	
	public void testSnpsNoInteraction() {
		snpdata = new SnpDataImputed(chromosome, files);
		totalSnps = snpdata.getNumberOfSnps();
		maxmarker = genotypes[0].length - 1;
		int nSamples = residualSamples.length;

		popIndex = new int[nSamples];
		pops = new ArrayList<String>();
		int count = 0;
		for (String sample : residualSamples) {
			int ipop = getPopulation(sample) - 1;
			pops.add(SnpData.popnames[ipop]);
			popIndex[count++] = ipop;
		}

		//create an index array into genotypes
//		sampleIndex = new int[nSamples];
//		for (int i = 0; i < nSamples; i++) sampleIndex[i] = genotypeMap.get(residualSamples[i]);
		
		//create the base model
		LinearModelforStepwiseRegression lmsr = getBaseModel();
		
		//setup the output
		BufferedWriter bw;
		if (files == null) {
//			bw = openOutputFile("C:/Projects/NAM/flowering/results/scans/d2a_test_snps_fastphase_" + thisTraitName + "_chr" + chromosome + ".txt");
			bw = openOutputFile("C:/Projects/NAM/leaf traits/results.test/chr10anglescan_leaftraits.txt");
		}
		else {
			bw = openOutputFile(files.chrsteps[chromosome - 1]);
		}
		
		StringBuilder sb = new StringBuilder("Position\tF\tp\tlog(1/p)");
		sb.append("\t").append("B97");
		sb.append("\t").append("CML103");
		sb.append("\t").append("CML228");
		sb.append("\t").append("CML247");
		sb.append("\t").append("CML277");
		sb.append("\t").append("CML322");
		sb.append("\t").append("CML333");
		sb.append("\t").append("CML52");
		sb.append("\t").append("CML69");
		sb.append("\t").append("Hp301");
		sb.append("\t").append("Il14H");
		sb.append("\t").append("Ki11");
		sb.append("\t").append("Ki3");
		sb.append("\t").append("Ky21");
		sb.append("\t").append("M162W");
		sb.append("\t").append("M37W");
		sb.append("\t").append("Mo17");
		sb.append("\t").append("Mo18W");
		sb.append("\t").append("MS71");
		sb.append("\t").append("NC350");
		sb.append("\t").append("NC358");
		sb.append("\t").append("Oh43");
		sb.append("\t").append("Oh7B");
		sb.append("\t").append("P39");
		sb.append("\t").append("Tx303");
		sb.append("\t").append("Tzi8");
		

		writeToOutput(sb.toString(), bw);

		double[] Fp = new double[]{0,1};
		while (snpdata.next()) {
			double[] parents = snpdata.getGenotype();
			int pos = snpdata.getPosition();
			
			//build and solve the model
			CovariateModelEffect cme = new CovariateModelEffect(projectSnp(parents, pos, popIndex));
			Fp =  lmsr.testNewEffect(cme);
			
			sb = new StringBuilder();
			sb.append(pos);
			sb.append("\t").append(Fp[0]);
			sb.append("\t").append(Fp[1]);
			sb.append("\t").append(Math.log10(1/Fp[1]));
			for (int i = 0; i < 26; i++) sb.append("\t").append(parents[i]);
			
			writeToOutput(sb.toString(), bw);
		}
		
		closeOutput(bw);
		
	}
	
	public void testSnpsWithInteraction() {
		snpdata = new SnpDataImputed(chromosome);
		totalSnps = snpdata.getNumberOfSnps();
		maxmarker = genotypes[0].length - 1;
		int nSamples = residualSamples.length;

		popIndex = new int[nSamples];
		pops = new ArrayList<String>();
		int count = 0;
		for (String sample : residualSamples) {
			int ipop = getPopulation(sample) - 1;
			pops.add(SnpData.popnames[ipop]);
			popIndex[count++] = ipop;
		}

		//create model effect for populations
		popLabels = new ArrayList<String>();
		ModelEffect mepop = new ModelEffect(ModelEffect.getIntegerLevels(pops, popLabels));
		
		//create an index array into genotypes
		sampleIndex = new int[nSamples];
		for (int i = 0; i < nSamples; i++) sampleIndex[i] = genotypeMap.get(residualSamples[i]);
		
		//create the base model
//		ArrayList<Integer> includedSnps = new ArrayList<Integer>();
//		includedSnps.add(21443991);
//		LinearModelforStepwiseRegression lmsr = getBaseModelWithQTL(pops, 30, includedSnps);
//		ArrayList<ModelEffect> baseEffects = new ArrayList<ModelEffect>(lmsr.getModelEffects());
		ArrayList<ModelEffect> baseEffects = new ArrayList<ModelEffect>();
		baseEffects.add(new ModelEffect(new int[nSamples]));
		baseEffects.add(new ModelEffect(ModelEffect.getIntegerLevels(pops)));
		int nterms = baseEffects.size();
		
		//get the data
		//get data for this chromosome 
		double[] y = new double[nSamples];
		for (int i = 0; i < nSamples; i++) y[i] = residuals[i][chromosome - 1];

		//setup the output
		BufferedWriter bw = openOutputFile("C:/Projects/NAM/leaf traits/results/qtl30b.test_snps_fastphase_" + thisTraitName + "_chr" + chromosome + ".txt");
		StringBuilder sb = new StringBuilder("Position\tmarkerF\tmarkerP\tlogmarkerP\tpopMarkerF\tpopMarkerP\tlogpopMarkerP");
		for (int i = 0; i < 26; i++) sb.append("\t").append(popLabels.get(i));
		writeToOutput(sb.toString(), bw);

		snpcount = 0;
		while (snpdata.next()) {
			snpcount++;
			double[] parents = snpdata.getGenotype();
			int pos = snpdata.getPosition();
			
			//build and solve the model
			CovariateModelEffect cme = new CovariateModelEffect(projectSnp(parents, pos, popIndex));
			baseEffects.add(cme);
			ModelEffect me = new NestedCovariateModelEffect(cme, mepop);
			baseEffects.add(me);
			LinearModelWithSweep lmws = new LinearModelWithSweep(baseEffects, y);
			ArrayList<Double> ss = lmws.effectSS;
			ArrayList<Double> df = lmws.effectdf;
			double errss = lmws.getErrorSS();
			double errdf = lmws.getErrordf();
			double markerss = ss.get(nterms);
			double markerdf = df.get(nterms);
			double popmarkerss = ss.get(nterms + 1);
			double popmarkerdf = df.get(nterms + 1);
			
			baseEffects.remove(nterms + 1);
			baseEffects.remove(nterms);
			
			double F, p;
			sb = new StringBuilder();
			sb.append(pos);
			
			F = markerss / markerdf / (errss + popmarkerss) * (errdf + popmarkerdf);
			try {p = AbstractLinearModel.Ftest(F, markerdf, (errdf + popmarkerdf));}
			catch (Exception e) {p = Double.NaN;}
			sb.append("\t").append(F);
			sb.append("\t").append(p);
			sb.append("\t").append(Math.log10(1/p));
			
			F = popmarkerss / popmarkerdf / errss * errdf;
			try {p = AbstractLinearModel.Ftest(F, markerdf, (errdf + popmarkerdf));}
			catch (Exception e) {p = Double.NaN;}
			sb.append("\t").append(F);
			sb.append("\t").append(p);
			sb.append("\t").append(Math.log10(1/p));
			
			for (int i = 0; i < 26; i++) sb.append("\t").append(parents[i]);
			
			writeToOutput(sb.toString(), bw);
		}
		
		closeOutput(bw);
	}
	
	//additionalSnps must be sorted in ascending order
	public LinearModelforStepwiseRegression getBaseModelWithQTL(int qtlToExclude) {
		int nSamples = residualSamples.length;
		int[] mean = new int[nSamples];
		NamMap namMap = new NamMap();
		
		//set up the mean and pop model effects
		ArrayList<ModelEffect> effects = new ArrayList<ModelEffect>();
		ModelEffect memean = new ModelEffect(mean);
		effects.add(memean);
		popLabels = new ArrayList<String>();
		ModelEffect mepop = new ModelEffect(ModelEffect.getIntegerLevels(pops, popLabels));
		effects.add(mepop);
		
		//get data for this chromosome 
		double[] y = new double[nSamples];
		for (int i = 0; i < nSamples; i++) y[i] = residuals[i][chromosome - 1];
		
		//and the QTL are ...
		ArrayList<Integer> qtlList = getQTLforChromosome();
		
		for (Integer marker : qtlList) {
			if (namMap.chromosome.get(marker) == chromosome && marker != qtlToExclude) {
				double[] markerscores = new double[nSamples];
				for (int i = 0; i < nSamples; i++) {
					markerscores[i] = genotypes[sampleIndex[i]][marker];
				}
				
				CovariateModelEffect cme = new CovariateModelEffect(markerscores);
				ModelEffect me = new NestedCovariateModelEffect(cme, mepop);
				me.setId(marker);
				effects.add(me);
			}
		}
		
		//add snps
		snpdata.reset();
		if (snpList != null && snpList.size() > 0) {
			boolean notEOF = snpdata.next();
			for (Integer snp : snpList) {
				
				while (notEOF && snpdata.getPosition() != snp) notEOF = snpdata.next();
				if (notEOF) {
					double[] parents = snpdata.getGenotype();
					int pos = snpdata.getPosition();
					CovariateModelEffect cme = new CovariateModelEffect(projectSnp(parents, pos, popIndex));
					effects.add(cme);
				}
			}
			
		}
		
		//create the base model
		return new LinearModelforStepwiseRegression(effects, y);
	}
	
	//additionalSnps must be sorted in ascending order
	public LinearModelforStepwiseRegression getBaseModel() {
		int nSamples = residualSamples.length;
		int[] mean = new int[nSamples];
		
		//set up the mean and pop model effects
		ArrayList<ModelEffect> effects = new ArrayList<ModelEffect>();
		ModelEffect memean = new ModelEffect(mean);
		effects.add(memean);
		popLabels = new ArrayList<String>();
		ModelEffect mepop = new ModelEffect(ModelEffect.getIntegerLevels(pops, popLabels));
		effects.add(mepop);
		
		//get data for this chromosome 
		double[] y = new double[nSamples];
		for (int i = 0; i < nSamples; i++) y[i] = residuals[i][chromosome - 1];
		
		//add snps
		snpdata.reset();
		if (snpList != null && snpList.size() > 0) {
			boolean notEOF = snpdata.next();
			for (Integer snp : snpList) {
				while (notEOF && snpdata.getPosition() != snp) notEOF = snpdata.next();
				if (notEOF) {
					double[] parents = snpdata.getGenotype();
					int pos = snpdata.getPosition();
					CovariateModelEffect cme = new CovariateModelEffect(projectSnp(parents, pos, popIndex));
					effects.add(cme);
				}
			}
			
		}
		
		//create the base model
		return new LinearModelforStepwiseRegression(effects, y);
	}
	
	public ArrayList<Integer> getQTLforChromosome() {
		String dir = "C:/Projects/NAM/leaf traits/models/";
		String prefix = "StepwiseWithScan082809_leaftrait_distanceImputation_";
		String suffix = "_model.txt";
		String filename = dir + prefix + traitnumber + suffix;
		
		BufferedReader br = Utilities.openReader(filename);
		
		String input;
		ArrayList<Integer> qtlList = new ArrayList<Integer>();
		while ((input = Utilities.readALineFrom(br)) != null) {
			String[] parsed = input.split("\t");
			try { qtlList.add(Integer.parseInt(parsed[1])); }
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
			snpvalues[i] = genotypes[i][leftmarker] * (1 - pd) + genotypes[i][rightmarker] * pd;
			snpvalues[i] = snpvalues[i] * parents[popIndex[i]];
		}
		return snpvalues;
	}

	public void setSnpList(ArrayList<Integer> snpList) {
		this.snpList = snpList;
		Collections.sort(this.snpList);
	}
}
