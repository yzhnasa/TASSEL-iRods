package net.maizegenetics.gwas.NAMgwas;

import java.util.ArrayList;

import net.maizegenetics.jGLiM.dm.ModelEffect;

public class ForwardRegressionModelSubsample extends ForwardRegressionModel {
	int[] sampleIndex;
	int numberOfSamples;
	
	public ForwardRegressionModelSubsample(ArrayList<ModelEffect> initialEffects, double[] phenotype, double enterlimit, int[] sampleIndex) {
		super(initialEffects, phenotype, enterlimit);
		this.sampleIndex = sampleIndex;
		numberOfSamples = sampleIndex.length;
	}

	@Override
	public void addNextSnp(SnpInfo snp) {
		double[] subgeno = new double[numberOfSamples];
		double[] geno = snp.genotype;
		for (int i = 0; i < numberOfSamples; i++) {
			subgeno[i] = geno[sampleIndex[i]];
		}
		double modelss = lmsr.testNewEffect(subgeno);
		if (modelss > bestModelSS) {
			bestModelSS = modelss;
			bestsnp = new SnpInfo(snp.chromosome, snp.pos, snp.allele, subgeno, 1, 1);
		}
	}
	
	
}
