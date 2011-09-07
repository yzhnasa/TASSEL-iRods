package net.maizegenetics.gwas.jointlinkage;

import java.util.concurrent.Callable;

import net.maizegenetics.jGLiM.dm.FactorModelEffect;
import net.maizegenetics.jGLiM.dm.NestedCovariateModelEffect;
import net.maizegenetics.jGLiM.dm.PartitionedLinearModel;

public class ForwardStepJointLinkage implements Callable<ForwardStepResult> {
	PartitionedLinearModel plm;
	SNPdata snps;
	FactorModelEffect popEffect;
	int startIndex;
	int stopIndex;
	
	public ForwardStepJointLinkage(PartitionedLinearModel plm, SNPdata snps, FactorModelEffect popEffect, int startSnp, int stopSnp) {
		this.plm = plm;
		this.snps = snps;
		this.popEffect = popEffect;
		startIndex = startSnp;
		stopIndex = stopSnp;
	}
	
	@Override
	public ForwardStepResult call() throws Exception {
		double bestss = 0;
		NestedCovariateModelEffect besteffect = null;
		
		//test all the snps. Select the one that gives the highest model SS
		for (int s = startIndex; s <= stopIndex; s++) {
			SNP snp = snps.getSnp(s);
			if (snp != null) {
				NestedCovariateModelEffect nested = new NestedCovariateModelEffect(snp.score, popEffect);
				nested.setID(snp);
				plm.testNewModelEffect(nested);
				double modelss = plm.getModelSS();
				if (modelss > bestss) {
					bestss = modelss;
					besteffect = nested;
				}
			}
		}

		plm.testNewModelEffect(besteffect);
		double[] Fp = plm.getFp();
		ForwardStepResult result = new ForwardStepResult();
		result.bestEffect = besteffect;
		result.bestModelSS = bestss;
		result.bestF = Fp[0];
		result.bestp = Fp[1];
		return result;
	}

}
