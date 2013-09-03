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
