package net.maizegenetics.matrixalgebra.decomposition;

import java.util.Arrays;

import org.apache.log4j.Logger;

import net.maizegenetics.matrixalgebra.Matrix.BlasDoubleMatrix;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrix;

public class BlasSingularValueDecomposition implements SingularValueDecomposition {
	private BlasDoubleMatrix bdS = null;
	private BlasDoubleMatrix bdU = null;
	private BlasDoubleMatrix bdVT = null;
	public boolean successful = false;
	private static final Logger myLogger = Logger.getLogger(BlasSingularValueDecomposition.class);
	private double tol = 1e-12;
	
	public BlasSingularValueDecomposition(BlasDoubleMatrix bdm, char jobz) {
		int nrows = bdm.numberOfRows();
		int ncols = bdm.numberOfColumns();
		int ns = Math.min(nrows, ncols);
		double[] A = bdm.getMatrixCopy();
		double[] S = new double[ns];
		int usize = nrows * nrows;
		
		double[] U;
		if (jobz == 'N' || (jobz == 'O'& nrows >= ncols)) {
			U = new double[]{0};
		} else {
			U = new double[usize];
		}
		
		double[] VT;
		int vtsize = ncols * ncols;
		if (jobz == 'N' || (jobz == 'O'& nrows < ncols)) {
			VT = new double[]{0};
		} else {
			VT = new double[vtsize];
		}
		
		int result = BlasDoubleMatrix.singularValueDecompositionDgesdd(jobz, nrows, ncols, A, nrows, S, U, nrows, VT, ncols);
		if (result == 0) {
			successful = true;
			bdS = new BlasDoubleMatrix(ns, 1, S, true);
			if (jobz == 'A') {
				bdU = new BlasDoubleMatrix(nrows, nrows, U, true);
				bdVT = new BlasDoubleMatrix(ncols, ncols, VT, true);
			} else if (jobz == 'S') {
				if (nrows < ncols) {
					bdU = new BlasDoubleMatrix(nrows, nrows, U, true);
					int[] selection = new int[nrows];
					for (int i = 0; i < nrows; i++) selection[i] = i;
					double[] subvt = BlasDoubleMatrix.getSelectionFromDoubleArray(VT, ncols, ncols, selection, null);
					bdVT = new BlasDoubleMatrix(nrows, ncols, subvt, true);
				} else if (nrows > ncols) {
					int[] selection = new int[ncols];
					for (int i = 0; i < ncols; i++) selection[i] = i;
					double[] subu = BlasDoubleMatrix.getSelectionFromDoubleArray(U, nrows, nrows, null, selection);
					bdU = new BlasDoubleMatrix(nrows, ncols, subu, true);
					bdVT = new BlasDoubleMatrix(ncols, ncols, VT, true);
				} else {
					bdU = new BlasDoubleMatrix(nrows, nrows, U, true);
					bdVT = new BlasDoubleMatrix(ncols, ncols, VT, true);
				}
			} else if (jobz == 'O') {
				if (nrows >= ncols) {
					bdVT = new BlasDoubleMatrix(ncols, ncols, VT, true);
					int[] selection = new int[ncols];
					for (int i = 0; i < ncols; i++) selection[i] = i;
					U = BlasDoubleMatrix.getSelectionFromDoubleArray(A, nrows, ncols, null, selection);
					bdU = new BlasDoubleMatrix(nrows, ncols, U, true);
				} else {
					bdU = new BlasDoubleMatrix(nrows, nrows, U, true);
					int[] selection = new int[nrows];
					for (int i = 0; i < nrows; i++) selection[i] = i;
					VT = BlasDoubleMatrix.getSelectionFromDoubleArray(A, nrows, ncols, selection, null);
					bdVT = new BlasDoubleMatrix(nrows, ncols, VT, true);
				}
			}
			
		} else {
			myLogger.error("BlasSVD failed with a return value of " + result);
		}
	}
	
	public BlasSingularValueDecomposition(BlasDoubleMatrix bdm) {
		this(bdm, 'A');
	}

	@Override
	public DoubleMatrix getU(boolean transpose) {
		if (transpose) return bdU.transpose();
		else return bdU;
	}

	@Override
	public DoubleMatrix getV(boolean transpose) {
		if (transpose) return bdVT;
		else return bdVT.transpose();
	}

	@Override
	public DoubleMatrix getS() {
		return bdS;
	}

	@Override
	public double[] getSingularValues() {
		return bdS.getMatrix();
	}

	@Override
	public int getRank() {
		int rank = 0;
		double[] sv = bdS.getMatrix();
		while (sv[rank] > tol) rank++;
		return rank;
	}

	public double getTolerance() { return tol; }
	
	public void setTolerance(double tolerance) { tol = tolerance; }
	
	public boolean wasSuccessful() { return successful; }
}
