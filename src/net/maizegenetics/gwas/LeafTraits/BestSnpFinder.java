package net.maizegenetics.gwas.LeafTraits;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import net.maizegenetics.jGLiM.CovariateModelEffect;
import net.maizegenetics.jGLiM.LinearModelforStepwiseRegression;
import net.maizegenetics.jGLiM.ModelEffect;

public class BestSnpFinder extends ModelFitter {
	double[] phenotype;
	int maxmarker;
	int[] sampleIndex;
	ArrayList<ModelEffect> baseEffects;
	SnpData snpdata;
	CovariateModelEffect bestEffect;
	double maxF;
	double minp;
	String[] sampleNames;
	int[] popIndex;
	ArrayList<Integer> pops;
	StepwiseGenomeWideSnpFinder parent;
	
	public BestSnpFinder(StepwiseGenomeWideSnpFinder parent, FileNames files, int chromosome, double[] phenotype, String[] samples, int[] popIndex) {
		super(chromosome, 0, files);
		this.parent = parent;
		this.phenotype = phenotype;
		sampleNames = samples;
		this.popIndex = popIndex;
		
		//create an index array into genotypes
		importNamMarkersforMap();
		int nSamples = phenotype.length;
		sampleIndex = new int[nSamples];
		for (int i = 0; i < nSamples; i++) sampleIndex[i] = genotypeMap.get(sampleNames[i]);
		
		snpdata = new SnpData(chromosome);
		totalSnps = snpdata.numberOfSnps;
	}

	public void run() {
		testSnps();
	}
	
	public void testSnps() {
		LinearModelforStepwiseRegression lmsr = new LinearModelforStepwiseRegression(baseEffects, phenotype);
		snpdata.reset();
		maxF = 0;
		minp = 1;
		bestEffect = null;
		snpcount = 0;
		while (snpdata.next()) {
			CovariateModelEffect cme = new CovariateModelEffect(projectSnp(snpdata.getGenotype(), snpdata.getPosition()));
			cme.setId(new int[]{chromosome,snpdata.getPosition()});
			double Fp[] = lmsr.testNewEffect(cme);
			if (Fp[1] < minp || (Fp[1] == minp && Fp[0] > maxF)) {
				minp = Fp[1];
				maxF = Fp[0];
				bestEffect = cme;
			}
			snpcount++;
		}
		parent.thisChromosomeIsFinished(chromosome);
	}
	
	public double[] projectSnp(double[] parents, int pos) {
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
		
		int nSamples = phenotype.length;
		double[] snpvalues = new double[nSamples];
		for (int i = 0; i < nSamples; i++) {
			snpvalues[i] = genotypes[sampleIndex[i]][leftmarker] * (1 - pd) + genotypes[sampleIndex[i]][rightmarker] * pd;
			snpvalues[i] = snpvalues[i] * parents[popIndex[i]];
		}
		return snpvalues;
	}

	public void setBaseEffects(ArrayList<ModelEffect> baseEffects) {
		this.baseEffects = baseEffects;
	}
	
	public int getBestPosition() {
		return ((int[]) bestEffect.getId())[1];
	}
	
	public double getBestCm() {
		int[] pos = (int[]) bestEffect.getId();
		return theAGPMap.getCmFromPosition(pos[0], pos[1]);
	}
}
