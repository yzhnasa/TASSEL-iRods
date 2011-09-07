package net.maizegenetics.gwas.LeafTraits;

import java.io.BufferedWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map.Entry;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;

import net.maizegenetics.jGLiM.AbstractLinearModel;
import net.maizegenetics.jGLiM.CovariateModelEffect;
import net.maizegenetics.jGLiM.LinearModelforStepwiseRegression;
import net.maizegenetics.jGLiM.ModelEffect;

public class ModelFitterStepwiseNoMissing extends ModelFitter {
	
	double enterLimit = 1e-6;
	
	public ModelFitterStepwiseNoMissing(int chromosome, int trait) {
		super(chromosome, trait);
	}
	
	public ModelFitterStepwiseNoMissing(int chromosome, FileNames files) {
		super(chromosome, files);
		if (!Double.isNaN(files.enterlimit)) enterLimit = files.enterlimit;
	}
	
	//fits snps for a chromosome using forward regression
	//data are the residuals from a model including QTL for other chromosomes
	
	public void testSnps() {
		if (files.analysis.equals(Analysis.ANALYSIS_HAPLOTYPE)) {
			testHaplotypes();
			return;
		}
		SnpData snpdata;
		if (files == null) snpdata = new SnpDataImputed(chromosome); 
		else snpdata = new SnpDataImputed(chromosome, files);
		totalSnps = snpdata.getNumberOfSnps();
		int maxmarker = genotypes[0].length - 1;
		int nSamples = residualSamples.length;

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
//		int[] sampleIndex = new int[nSamples];
//		for (int i = 0; i < nSamples; i++) sampleIndex[i] = genotypeMap.get(residualSamples[i]);
		
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
		LinearModelforStepwiseRegression lmsr = new LinearModelforStepwiseRegression(effects, y);
		
		//open file for output
		BufferedWriter bwstep;
		if (files == null) {
			StringBuilder sb = new StringBuilder("C:/Projects/NAM/leaf traits/results/stepwise/stepwise_chr");
			sb.append(chromosome).append("_");;
			sb.append(thisTraitName).append(".txt");
			bwstep = Utilities.openOutputFile(sb.toString());
		}
		else bwstep = Utilities.openOutputFile(files.chrsteps[chromosome - 1]);
		Utilities.flushToOutput("chromosome\tposition\tcM\tp-value", bwstep);
		
		//forward regression
		double minp;
		do {
			//create the base model
			minp = 1;
			int bestpos = -1;
			ModelEffect bestEffect = null;
			double maxF = 0;

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

				double[] snpvalues = new double[nSamples];
				for (int i = 0; i < nSamples; i++) {
//					snpvalues[i] = genotypes[sampleIndex[i]][leftmarker] * (1 - pd) + genotypes[sampleIndex[i]][rightmarker] * pd;
					snpvalues[i] = genotypes[i][leftmarker] * (1 - pd) + genotypes[i][rightmarker] * pd;
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
				snpdata.reset();
				lmsr.addEffect(bestEffect);
				StringBuilder sb = new StringBuilder();
				sb.append(chromosome).append("\t").append(bestpos);
				sb.append("\t").append(theAGPMap.getCmFromPosition(chromosome, bestpos));
				sb.append("\t").append(minp);
				Utilities.flushToOutput(sb.toString(), bwstep);
				System.out.println(sb.toString());
			}
			
		} while (minp < enterLimit);
		
		Utilities.closeOutput(bwstep);
		
		//write locus effect estimates to a file
		DoubleMatrix1D beta = lmsr.getLinearModel().getBeta();
		int nbeta = beta.size();
		int neffects = effects.size();
		int nloci = neffects - 2;
		int firsteffect = nbeta - nloci;
		
		BufferedWriter bw;
		if (files == null) {
			StringBuilder sb = new StringBuilder("C:/Projects/NAM/leaf traits/results/stepwise/");
			sb.append("effects_").append(thisTraitName);
			sb.append("_chr").append(chromosome).append(".txt");
			bw = Utilities.openOutputFile(sb.toString());
		}
		else bw = Utilities.openOutputFile(files.chrmodel[chromosome - 1]);
		
		Utilities.writeToOutput("chromosome\tposition\tcM\tallele\teffect\tpvalue", bw);
		double errSS = lmsr.getLinearModel().getErrorSS();
		double errdf = lmsr.getLinearModel().getErrordf();
		
		//the mean
		StringBuilder sb = new StringBuilder();
		sb.append("overall mean\t\t\t\t");
		sb.append(beta.getQuick(0));
		Utilities.writeToOutput(sb.toString(), bw);
		
		//the populations
		sb = new StringBuilder();
		sb.append("population\tall\t\t\t\t");
		double[] effectSSdf = lmsr.getLinearModel().marginalEffectSSdf(1);
		double F = effectSSdf[0] / effectSSdf[1] / errSS * errdf;
		double p;
		try {p = AbstractLinearModel.Ftest(F, effectSSdf[1], errdf);}
		catch (Exception e) {p = Double.NaN;}
		sb.append(p);
		Utilities.writeToOutput(sb.toString(), bw);
		
		int npops = nbeta - 1 - nloci;
		for (int i = 1; i <= npops; i++) {
			sb = new StringBuilder("population\t");
			sb.append(i).append("\t\t\t").append(beta.getQuick(i));
			Utilities.writeToOutput(sb.toString(), bw);
		}
		
		for (int i = 0; i < nloci; i++) {
			Integer pos = (Integer) effects.get(2 + i).getId();
			double cM = theAGPMap.getCmFromPosition(chromosome, pos);
			double estEffect = beta.getQuick(i + firsteffect);
			effectSSdf = lmsr.getLinearModel().marginalEffectSSdf(2 + i);
			F = effectSSdf[0] / effectSSdf[1] / errSS * errdf;
			try {p = AbstractLinearModel.Ftest(F, effectSSdf[1], errdf);}
			catch (Exception e) {p = Double.NaN;}
			Utilities.writeToOutput(chromosome + "\t" + pos + "\t" + cM + "\t" + estEffect + "\t" + p, bw);
		}
		Utilities.closeOutput(bw);
		System.out.println("Finished chromosome " + chromosome + ".");
	}

	public void testHaplotypes() {
		HapData hapdata;
		hapdata = new HapData(files, chromosome);
		totalSnps = hapdata.getNumberOfHaplotypes();
		int maxmarker = genotypes[0].length - 1;
		int nSamples = residualSamples.length;

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
		ArrayList<String> popLabels = new ArrayList<String>();
		ModelEffect mepop = new ModelEffect(ModelEffect.getIntegerLevels(pops, popLabels));
		effects.add(mepop);
		
		//get data for this chromosome 
		double[] y = new double[nSamples];
		for (int i = 0; i < nSamples; i++) y[i] = residuals[i][chromosome - 1];
		
		//create the base model
		LinearModelforStepwiseRegression lmsr = new LinearModelforStepwiseRegression(effects, y);
		
		//open file for output
		BufferedWriter bwstep;
		if (files == null) {
			StringBuilder sb = new StringBuilder("C:/Projects/NAM/leaf traits/results/stepwise/stepwise_chr");
			sb.append(chromosome).append("_");;
			sb.append(thisTraitName).append(".txt");
			bwstep = Utilities.openOutputFile(sb.toString());
		}
		else bwstep = Utilities.openOutputFile(files.chrsteps[chromosome - 1]);
		Utilities.flushToOutput("chromosome\tposition\tcM\tp-value", bwstep);
		
		//stick vgt1 in the model first
		int[] vgt1 = new int[26];
		vgt1[10] = 1;
		vgt1[16] = 1;
		vgt1[18] = 1;
		vgt1[23] = 1;
		int thispos = 130661845;
		Object[][] _markers = theAGPMap.getInterval(chromosome, thispos);
		int _left = AGPMap.getMarkerPosition(_markers[0]);
		if (_left == -1) _left = 0;
		int _right = AGPMap.getMarkerPosition(_markers[1]);
		if (_right == -1) _right = chromosomeLength[chromosome - 1];

		//proportion of distance of snp between left and right markers
		double _pd = ((double)(thispos - _left)) / ((double)(_right - _left));

		int _leftmarker = AGPMap.getMarkerNumber(_markers[0]);
		if (_leftmarker == -1) _leftmarker = 0; //telomere
		else _leftmarker = _leftmarker - firstMarker + 1;
		int _rightmarker = AGPMap.getMarkerNumber(_markers[1]);
		if (_rightmarker == -1) _rightmarker = maxmarker; //telomere
		else _rightmarker = _rightmarker - firstMarker + 1;

		HashMap<Integer, Integer> _hapmap = new HashMap<Integer, Integer>();
		int _nparents = vgt1.length;
		count = 0;
		for (int i = 0; i < _nparents; i++) {
			if (vgt1[i] > 0) {
				Integer _hap = _hapmap.get(vgt1[i]);
				if (_hap == null) _hapmap.put(vgt1[i], count++);
			}
		}
		
		int _nhaps = _hapmap.size();
		double[][] _alleles = new double[_nhaps][nSamples];
		for (int i = 0; i < nSamples; i++) {
			int ndx = vgt1[popIndex[i]];
			if (ndx > 0) {
				_alleles[_hapmap.get(ndx)][i] = genotypes[sampleIndex[i]][_leftmarker] * (1 - _pd) + genotypes[sampleIndex[i]][_rightmarker] * _pd;
			}
		}
		
		Integer[] _haplotypes = new Integer[_nhaps];
		for (Entry<Integer,Integer> ent : _hapmap.entrySet()) {
			_haplotypes[ent.getValue()] = ent.getKey();
		}
		
		MultipleCovariateModelEffect _mcme = new MultipleCovariateModelEffect(_alleles);
		_mcme.setId(new Object[]{new Integer(thispos), _haplotypes});
		lmsr.addEffect(_mcme);
		//end: stick vgt1 in the model first

		
		//forward regression
		double minp;
		do {
			//create the base model
			minp = 1;
			int bestpos = -1;
			ModelEffect bestEffect = null;
			double maxF = 0;

			snpcount = 0;
			hapdata.reset();
			while (hapdata.next()) {
				snpcount++;
				//project haplotypes
				int[] parents = hapdata.getGenotype();
				int pos = hapdata.getPosition();
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

				HashMap<Integer, Integer> hapmap = new HashMap<Integer, Integer>();
				int nparents = parents.length;
				count = 0;
				for (int i = 0; i < nparents; i++) {
					if (parents[i] > 0) {
						Integer hap = hapmap.get(parents[i]);
						if (hap == null) hapmap.put(parents[i], count++);
					}
				}
				
				int nhaps = hapmap.size();
				if (nhaps == 0) continue; //skips the rest of the loop if all haplotypes = B73
				
				double[][] alleles = new double[nhaps][nSamples];
				for (int i = 0; i < nSamples; i++) {
					int ndx = parents[popIndex[i]];
					if (ndx > 0) {
						alleles[hapmap.get(ndx)][i] = genotypes[sampleIndex[i]][leftmarker] * (1 - pd) + genotypes[sampleIndex[i]][rightmarker] * pd;
					}
				}
				
				Integer[] haplotypes = new Integer[nhaps];
				for (Entry<Integer,Integer> ent : hapmap.entrySet()) {
					haplotypes[ent.getValue()] = ent.getKey();
				}
				
				//build and solve the model
				MultipleCovariateModelEffect mcme = new MultipleCovariateModelEffect(alleles);
				mcme.setId(new Object[]{new Integer(pos), haplotypes});
				double[] Fp =  lmsr.testNewEffect(mcme);
				
				if (!Double.isNaN(Fp[1]) && (Fp[1] < minp || (Fp[1] == minp && Fp[0] > maxF) )) {
					minp = Fp[1];
					maxF = Fp[0];
					bestEffect = mcme;
					bestpos = pos;
				}
			}
			
			if (minp < enterLimit) {
				hapdata.reset();
				lmsr.addEffect(bestEffect);
				StringBuilder sb = new StringBuilder();
				sb.append(chromosome).append("\t").append(bestpos);
				sb.append("\t").append(theAGPMap.getCmFromPosition(chromosome, bestpos));
				sb.append("\t").append(minp);
				Utilities.flushToOutput(sb.toString(), bwstep);
				System.out.println(sb.toString());
			}
			
		} while (minp < enterLimit);
		
		Utilities.closeOutput(bwstep);
		
		//write locus effect estimates to a file
		BufferedWriter bw;
		if (files == null) {
			StringBuilder sb = new StringBuilder("C:/Projects/NAM/leaf traits/results/stepwise/");
			sb.append("effects_").append(thisTraitName);
			sb.append("_chr").append(chromosome).append(".txt");
			bw = Utilities.openOutputFile(sb.toString());
		}
		else bw = Utilities.openOutputFile(files.chrmodel[chromosome - 1]);
		
		Utilities.writeToOutput("chromosome\tposition\tcM\thaplotype\teffect\tpvalue", bw);
		double errSS = lmsr.getLinearModel().getErrorSS();
		double errdf = lmsr.getLinearModel().getErrordf();
		
		DoubleMatrix1D beta = lmsr.getLinearModel().getBeta();
		int nbeta = beta.size();
		int neffects = effects.size();
		int nloci = neffects - 2;
		
		//write results of overall locus test
		for (int i = 0; i < nloci; i++) {
			ModelEffect thiseffect = effects.get(i + 2);
			Object[] thisId = (Object[]) thiseffect.getId(); 
			Integer pos = (Integer) thisId[0];
			double cM = theAGPMap.getCmFromPosition(chromosome, pos);
			double[] effectSSdf = lmsr.getLinearModel().marginalEffectSSdf(2 + i);
			double F = effectSSdf[0] / effectSSdf[1] / errSS * errdf;
			double p;
			try {p = AbstractLinearModel.Ftest(F, effectSSdf[1], errdf);}
			catch (Exception e) {p = Double.NaN;}
			Utilities.writeToOutput(chromosome + "\t" + pos + "\t" + cM + "\tall\t\t" + p, bw);
		}
		
		//write haplotype effect estimates
		int currentParameter = effects.get(0).getNumberOfLevels();
		currentParameter += effects.get(1).getNumberOfLevels();
		
		DoubleMatrix2D G = lmsr.getLinearModel().getInverse();
		for (int i = 0; i < nloci; i++) {
			ModelEffect thiseffect = effects.get(i + 2);
			int levels = thiseffect.getNumberOfLevels();
			Object[] thisId = (Object[]) thiseffect.getId(); 
			Integer pos = (Integer) thisId[0];
			double cM = theAGPMap.getCmFromPosition(chromosome, pos);
			Integer[] haplotypes = (Integer[]) thisId[1];
			for (int j = 0; j < levels; j++) {
				double estEffect = beta.get(currentParameter);
				double F = estEffect * estEffect / G.getQuick(currentParameter, currentParameter) / errSS * errdf;
				double p;
				try {p = AbstractLinearModel.Ftest(F, 1, errdf);}
				catch (Exception e) {p = Double.NaN;}
				Utilities.flushToOutput(chromosome + "\t" + pos + "\t" + cM + "\t" + haplotypes[j] + "\t" + estEffect + "\t" + p, bw);
				currentParameter++;
			}
		}
		
		Utilities.closeOutput(bw);
		System.out.println("Finished chromosome " + chromosome + ".");
	}
	
	public void setEnterLimit(double limit) {enterLimit = limit;}

}
