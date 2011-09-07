package net.maizegenetics.gwas.jointlinkage;

import net.maizegenetics.jGLiM.dm.NestedCovariateModelEffect;

public class ForwardStepResult {
	NestedCovariateModelEffect bestEffect;
	double bestModelSS;
	double bestF;
	double bestp;
}
