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

    /**
     * This returns a SuperByteMatrix designed for better performance when
     * column iteration loop inside row iteration loop.
     *
     * @param numRows number of rows
     * @param numColumns number of columns
     *
     * @return SuperByteMatrix (double dimension byte array)
     */
    public static SuperByteMatrix getInstance(int numRows, int numColumns) {
        long numElements = (long) numRows * (long) numColumns;
        if (numElements > (long) (Integer.MAX_VALUE - 10)) {
            return new SuperByteMatrixMultiple(numRows, numColumns);
        } else {
            return new SuperByteMatrixSingle(numRows, numColumns);
        }
    }

    /**
     * This returns a SuperByteMatrix designed for better performance when row
     * iteration loop inside column iteration loop.
     *
     * @param numRows number of rows
     * @param numColumns number of columns
     *
     * @return SuperByteMatrix (double dimension byte array)
     */
    public static SuperByteMatrix getInstanceTranspose(int numRows, int numColumns) {
        return new SuperByteMatrixTranspose(numRows, numColumns);
    }

    public static SuperByteMatrix getInstanceTranspose(SuperByteMatrix matrix) {

        int blockSize = 64;
        int numRows = matrix.getNumRows();
        int numColumns = matrix.getNumColumns();
        SuperByteMatrix result;

        if ((matrix instanceof SuperByteMatrixSingle) || (matrix instanceof SuperByteMatrixMultiple)) {
            result = getInstanceTranspose(numRows, numColumns);
            int currentCSize = blockSize;
            byte[][] temp = new byte[blockSize][blockSize];
            for (int bigC = 0; bigC < numColumns; bigC += blockSize) {
                if (numColumns - bigC < blockSize) {
                    currentCSize = numColumns - bigC;
                }
                int currentRSize = blockSize;
                for (int bigR = 0; bigR < numRows; bigR += blockSize) {
                    if (numRows - bigR < blockSize) {
                        currentRSize = numRows - bigR;
                    }

                    for (int r = 0; r < currentRSize; r++) {
                        for (int c = 0; c < currentCSize; c++) {
                            temp[r][c] = matrix.get(r + bigR, c + bigC);
                        }
                    }

                    for (int c = 0; c < currentCSize; c++) {
                        for (int r = 0; r < currentRSize; r++) {
                            result.set(r + bigR, c + bigC, temp[r][c]);
                        }
                    }
                }
            }
        } else if (matrix instanceof SuperByteMatrixTranspose) {
            result = getInstance(numRows, numColumns);
            int currentCSize = blockSize;
            byte[][] temp = new byte[blockSize][blockSize];
            for (int bigC = 0; bigC < numColumns; bigC += blockSize) {
                if (numColumns - bigC < blockSize) {
                    currentCSize = numColumns - bigC;
                }
                int currentRSize = blockSize;
                for (int bigR = 0; bigR < numRows; bigR += blockSize) {
                    if (numRows - bigR < blockSize) {
                        currentRSize = numRows - bigR;
                    }

                    for (int c = 0; c < currentCSize; c++) {
                        for (int r = 0; r < currentRSize; r++) {
                            temp[r][c] = matrix.get(r + bigR, c + bigC);
                        }
                    }

                    for (int r = 0; r < currentRSize; r++) {
                        for (int c = 0; c < currentCSize; c++) {
                            result.set(r + bigR, c + bigC, temp[r][c]);
                        }
                    }
                }
            }
        } else {
            throw new IllegalArgumentException("SuperByteMatrixBuilder: getInstanceTranspose: Don't Know how to Transpose: " + matrix.getClass().getName());
        }

        return result;
    }
}
