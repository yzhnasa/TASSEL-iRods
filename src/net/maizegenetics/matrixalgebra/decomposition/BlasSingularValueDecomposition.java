package net.maizegenetics.matrixalgebra.decomposition;

import org.apache.log4j.Logger;

import net.maizegenetics.matrixalgebra.Matrix.BlasDoubleMatrix;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrix;

public class BlasSingularValueDecomposition implements
		SingularValueDecomposition {
//	use 	public static native int singularValueDecomposition(char jobu, char jobvt, int nrows, int ncols, int ns, double[] A, double[] S, double[] U, double[] VT);
//	in BlasDoubleMatrix
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
		double[] S = new double[ns];
		int urows = nrows;
		int ucols = ns;
		int usize = urows * ucols;
		double[] U = new double[usize];
		int vtrows = ns;
		int vtcols = ncols;
		int vtsize = vtrows * vtcols;
		double[] VT = new double[vtsize];

		int result = BlasDoubleMatrix.singularValueDecompositionDgesdd(jobz, nrows, ncols, bdm.getMatrix(), nrows, S, U, nrows, VT, ncols);
		if (result == 0) {
			successful = true;
			bdS = new BlasDoubleMatrix(ns, 1, S, true);
			bdU = new BlasDoubleMatrix(urows, ucols, U, true);
			bdVT = new BlasDoubleMatrix(vtrows, vtcols, VT, true);
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
}
