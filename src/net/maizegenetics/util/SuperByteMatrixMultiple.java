/*
 *  SuperByteMatrixMultiple
 */
package net.maizegenetics.util;

/**
 *
 * @author Terry Casstevens
 */
public class SuperByteMatrixMultiple implements SuperByteMatrix {

    private final byte[][] myData;
    private final int myNumRows;
    private final int myNumColumns;
    private final int myNumRowsPerSingleDimArray;

    SuperByteMatrixMultiple(int rows, int columns) {

        myNumRows = rows;
        myNumColumns = columns;

        long numElements = myNumRows * myNumColumns;
        myNumRowsPerSingleDimArray = Integer.MAX_VALUE / myNumColumns;
        int numElementsPerSingleDimArray = myNumRowsPerSingleDimArray * myNumColumns;
        int numSingleDimArrays = (int) (numElements / (long) numElementsPerSingleDimArray);
        int numRemaining = (int) (numElements % (long) numElementsPerSingleDimArray);
        if (numRemaining != 0) {
            myData = new byte[numSingleDimArrays + 1][];
            for (int i = 0; i < numSingleDimArrays; i++) {
                myData[i] = new byte[numElementsPerSingleDimArray];
            }
            myData[numSingleDimArrays] = new byte[numRemaining];
        } else {
            myData = new byte[numSingleDimArrays][];
            for (int i = 0; i < numSingleDimArrays; i++) {
                myData[i] = new byte[numElementsPerSingleDimArray];
            }
        }

    }

    public void set(int row, int column, byte value) {
        myData[getFirstIndex(row)][getSecondIndex(row, column)] = value;
    }

    public byte get(int row, int column) {
        return myData[getFirstIndex(row)][getSecondIndex(row, column)];
    }

    private int getFirstIndex(int row) {
        return row / myNumRowsPerSingleDimArray;
    }

    private int getSecondIndex(int row, int column) {
        return (row % myNumRowsPerSingleDimArray) * myNumColumns + column;
    }
}
