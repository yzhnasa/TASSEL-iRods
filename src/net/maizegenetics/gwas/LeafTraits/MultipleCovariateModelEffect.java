package net.maizegenetics.gwas.LeafTraits;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import net.maizegenetics.jGLiM.ModelEffect;

public class MultipleCovariateModelEffect extends ModelEffect {
	//each row is a covariate
    double[][] covariates;
    
    public MultipleCovariateModelEffect(double[][] covariates) {
        super();
        this.covariates = covariates;
        size = covariates[0].length;
        numberOfLevels = covariates.length;
    }

    public double getSumProducts(double[] x, double[] y) {
        double sp = 0;
        for (int i = 0; i < size; i++) {
            sp += x[i] * y[i];
        }
        return sp;
    }
    
    @Override
    public DoubleMatrix2D getXTCov(double[][] otherCovariates) {
        int n1 = covariates.length;
        int n2 = otherCovariates.length;
        DoubleMatrix2D XTCov = DoubleFactory2D.dense.make(n1,n2,0);
        for (int i = 0; i < n1; i++) {
        	for (int j = i; j < n2; j++) {
        		double sp = getSumProducts(covariates[i], otherCovariates[j]);
        		XTCov.setQuick(i, j, sp);
        		if(i != j) XTCov.setQuick(j, i, sp);
        	}
        }
        
        return XTCov;
    }

    @Override
    public DoubleMatrix2D getXTX() {
    	int n = covariates.length;
    	DoubleMatrix2D XTX = DoubleFactory2D.dense.make(n, n, 0);
    	for (int i = 0; i < n; i++) {
    		for (int j = i; j < n; j++) {
    			double sp = getSumProducts(covariates[i], covariates[j]);
    			XTX.setQuick(i, j, sp);
    			if (i != j) XTX.setQuick(j, i, sp);
    		}
    	}
        return XTX;
    }

    @Override
    public DoubleMatrix1D getXTy(double[] y) {
    	int n = covariates.length;
    	DoubleMatrix1D XTy = DoubleFactory1D.dense.make(n,0);
    	for (int i = 0; i < n; i++) XTy.setQuick(i, getSumProducts(covariates[i], y));
        return XTy;
    }

    @Override
    public DoubleMatrix1D getyhat(DoubleMatrix1D beta) {
    	DoubleMatrix2D cov = DoubleFactory2D.dense.make(covariates).viewDice();
        return cov.zMult(beta, null);
    }

    public double[][] getCovariates() {
        return covariates;
    }
}
