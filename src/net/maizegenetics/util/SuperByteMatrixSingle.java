/*
 *  SuperByteMatrixSingle
 */
package net.maizegenetics.util;

/**
 *
 * @author Terry Casstevens
 */
public class SuperByteMatrixSingle implements SuperByteMatrix {

    private final byte[] myData;
    private final int myNumRows;
    private final int myNumColumns;

    SuperByteMatrixSingle(int rows, int columns) {

        myNumRows = rows;
        myNumColumns = columns;

        long numElements = myNumRows * myNumColumns;
        if (numElements > Integer.MAX_VALUE) {
            throw new IllegalArgumentException("SuperByteMatrixSingle: init: this number of rows: " + rows + "  and columns: " + columns + " is too large for SuperByteMatrixSingle.");
        }
        myData = new byte[(int) numElements];

    }

    public void set(int row, int column, byte value) {
        myData[getIndex(row, column)] = value;
    }

    public byte get(int row, int column) {
        return myData[getIndex(row, column)];
    }

    private int getIndex(int row, int column) {
        return row * myNumColumns + column;
    }
}
