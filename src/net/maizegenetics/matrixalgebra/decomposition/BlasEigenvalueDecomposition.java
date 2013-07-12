package net.maizegenetics.matrixalgebra.decomposition;

import java.util.Arrays;

import net.maizegenetics.matrixalgebra.Matrix.BlasDoubleMatrix;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrix;

public class BlasEigenvalueDecomposition implements EigenvalueDecomposition {
	
	protected double[] eigenvalues;
	protected double[] eigenvectors;
	protected BlasDoubleMatrix bdm;
	protected int info;
	
	// use public static native int eigenValueSymmetricDecomposition(int order, double[] A, double[] eigval) from BlasDoubleMatrix
	public BlasEigenvalueDecomposition(DoubleMatrix dm) {
		bdm = (BlasDoubleMatrix) dm;
		init();
	}
	
	protected void init() {
		eigenvectors = Arrays.copyOf(bdm.getMatrix(), bdm.getSize());
		eigenvalues = new double[bdm.numberOfRows()];
		info = BlasDoubleMatrix.eigenValueSymmetricDecomposition(bdm.numberOfRows(), eigenvectors, eigenvalues);
	}
	
	@Override
	public double[] getEigenvalues() {
		return eigenvalues;
	}

	@Override
	public double getEigenvalue(int i) {
		return eigenvalues[i];
	}

	@Override
	public DoubleMatrix getEigenvectors() {
		int nrows = bdm.numberOfRows();
		return new BlasDoubleMatrix(nrows, nrows, eigenvalues, true);
	}  

	@Override
	public DoubleMatrix getEigenvalueMatrix() {
		return BlasDoubleMatrix.getDiagonalMatrix(eigenvalues);
	}

}
