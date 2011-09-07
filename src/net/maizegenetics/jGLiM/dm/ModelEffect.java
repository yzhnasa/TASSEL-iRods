package net.maizegenetics.jGLiM.dm;

import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrix;

public interface ModelEffect {
	/**
	 * @return	an identifier for this factor
	 */
	Object getID();
	
	/**
	 * @param id	an identifier for this factor
	 */
	void setID(Object id);
	
	/**
	 * @return	the number of observations in the data set
	 */
	int getSize();
	
	/**
	 * @return	the design matrix for this factor
	 */
	DoubleMatrix getX();
	
	/**
	 * @return the crossproduct of the design matrix for this factor and itself
	 */
	DoubleMatrix getXtX();
	
	/**
	 * @param y	the dependent variable
	 * @return	the product of the transpose of the design matrix for this factor and the data, y
	 */
	DoubleMatrix getXty(double[] y);
	
	/**
	 * @param beta	 the effect estimate for each level of this factor
	 * @return	 the predicted value of this factor for each observation
	 */
	DoubleMatrix getyhat(DoubleMatrix beta);
	
	/**
	 * @param beta	 the effect estimate for each level of this factor
	 * @return the predicted value of this factor for each observation
	 */
	DoubleMatrix getyhat(double[] beta);
	
	/**
	 * @return	the number of observations for each level of this factor
	 */
	int[] getLevelCounts();
	
	/**
	 * @return	the number of levels in this factor
	 */
	int getNumberOfLevels();
	
	/**
	 * @return	a copy of this model effect
	 */
	ModelEffect getCopy();
}
