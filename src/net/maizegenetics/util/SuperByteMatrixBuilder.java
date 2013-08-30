/*
 *  SuperByteMatrixBuilder
 */
package net.maizegenetics.util;

/**
 *
 * @author Terry Casstevens
 */
public class SuperByteMatrixBuilder {

    private SuperByteMatrixBuilder() {
    }

    public static SuperByteMatrix getInstance(int rows, int columns) {
        if ((long) rows * (long) columns > Integer.MAX_VALUE) {
            return new SuperByteMatrixMultiple(rows, columns);
        } else {
            return new SuperByteMatrixSingle(rows, columns);
        }
    }

    public static SuperByteMatrix getInstanceTranspose(int rows, int columns) {
        return new SuperByteMatrixTranspose(rows, columns);
    }

    public static SuperByteMatrix getInstanceTranspose(SuperByteMatrix matrix) {

        int numRows = matrix.getNumRows();
        int numColumns = matrix.getNumColumns();
        SuperByteMatrix result;

        if (matrix instanceof SuperByteMatrixSingle) {
            result = getInstanceTranspose(numRows, numColumns);
        } else if (matrix instanceof SuperByteMatrixTranspose) {
            result = getInstance(numRows, numColumns);
        } else {
            throw new IllegalArgumentException("SuperByteMatrixBuilder: getInstanceTranspose: Don't Know how to Transpose: " + matrix.getClass().getName());
        }
        for (int r = 0; r < numRows; r++) {
            for (int c = 0; c < numColumns; c++) {
                result.set(r, c, matrix.get(r, c));
            }
        }
        return result;
    }
}
