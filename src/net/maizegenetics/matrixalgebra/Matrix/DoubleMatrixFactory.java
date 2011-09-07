package net.maizegenetics.matrixalgebra.Matrix;

import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.SpecializedOps;

import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix2D;

public class DoubleMatrixFactory {
	public enum FactoryType {ejml, jblas, colt};
	private FactoryType myType;
	public static DoubleMatrixFactory DEFAULT = new DoubleMatrixFactory(FactoryType.ejml);
	
	public DoubleMatrixFactory(FactoryType type) {
		myType = type;
	}
	
	public static void setDefault(FactoryType type) {
		DoubleMatrixFactory.DEFAULT = new DoubleMatrixFactory(type);
	}
	
	public FactoryType getType() { return myType; }
	
	public DoubleMatrix make(int row, int col) {
		if (myType == FactoryType.ejml) return new EJMLDoubleMatrix(row, col);
		if (myType == FactoryType.colt) return new ColtDoubleMatrix(row, col);
		if (myType == FactoryType.jblas) return null;
		return null;
	}
	
	public DoubleMatrix make(int row, int col, double[] values) {
		if (myType == FactoryType.ejml) return new EJMLDoubleMatrix(row, col, values);
		if (myType == FactoryType.colt) return new ColtDoubleMatrix(row, col, values);
		if (myType == FactoryType.jblas) return null;
		return null;
	}

	public DoubleMatrix make(int row, int col, double value) {
		if (myType == FactoryType.ejml) return new EJMLDoubleMatrix(row, col, value);
		if (myType == FactoryType.colt) return new ColtDoubleMatrix(row, col, value);
		if (myType == FactoryType.jblas) return null;
		return null;
	}
	
	public DoubleMatrix make(double[][] values) {
		if (myType == FactoryType.ejml) return new EJMLDoubleMatrix(values);
		if (myType == FactoryType.colt) return new ColtDoubleMatrix(values);
		if (myType == FactoryType.jblas) return null;
		return null;
	}
	
	public DoubleMatrix identity(int n) {
		if (myType == FactoryType.ejml) return new EJMLDoubleMatrix(n);
		if (myType == FactoryType.colt) return new ColtDoubleMatrix(n);
		if (myType == FactoryType.jblas) return null;
		return null;
	}
	
	public DoubleMatrix diagonal(double[] diag) {
		if (myType == FactoryType.ejml) return new EJMLDoubleMatrix(diag);
		if (myType == FactoryType.colt) return new ColtDoubleMatrix(diag);
		if (myType == FactoryType.jblas) return null;
		return null;
	}
	
	public DoubleMatrix compose(DoubleMatrix[][] components) {
		
		if (myType == FactoryType.ejml) {
			int totalRows = 0;
			int totalCols = 0;
			int nRows = components.length;
			int nCols = components[0].length;
			for (int i = 0; i < nRows; i++) totalRows += components[i][0].numberOfRows();
			for (int i = 0; i < nCols; i++) totalCols += components[0][i].numberOfColumns();
			
			DenseMatrix64F result = new DenseMatrix64F(totalRows, totalCols);
			int startRow = 0;
			for (int r = 0; r < nRows; r++) {
				int startCol = 0;
				for (int c = 0; c < nCols; c++) {
					DenseMatrix64F dm = ((EJMLDoubleMatrix) components[r][c]).myMatrix;
					SpecializedOps.insert(dm, startRow, startCol, result);
					startCol += dm.numCols;
				}
				startRow += components[r][0].numberOfRows();
			}
			return new EJMLDoubleMatrix(result);
		}
		
		if (myType == FactoryType.colt) {
			int nRows = components.length;
			int nCols = components[0].length;
			DoubleMatrix2D[][] coltComponents = new DoubleMatrix2D[nRows][nCols];
			for (int r = 0; r < nRows; r++) {
				for (int c = 0; c < nCols; c++) {
					coltComponents[r][c] = ((ColtDoubleMatrix) components[r][c]).myMatrix;
				}
			}
			return new ColtDoubleMatrix(DoubleFactory2D.dense.compose(coltComponents));
		}
		
		if (myType == FactoryType.jblas) return null;
		return null;
	}
}
