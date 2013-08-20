/*
 *  SuperByteMatrix
 */
package net.maizegenetics.util;

/**
 *
 * @author Terry Casstevens
 */
public class SuperByteMatrix {

    private final byte[][] myData;
    private final int myNumRows;
    private final int myNumColumns;
    private final int myNumRowsPerSingleDimArray;

    public SuperByteMatrix(int rows, int columns) {
        
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
                myData[i] = new byte[Integer.MAX_VALUE];
            }
            myData[numSingleDimArrays] = new byte[numRemaining];
        } else {
            myData = new byte[numSingleDimArrays][];
            for (int i = 0; i < numSingleDimArrays; i++) {
                myData[i] = new byte[Integer.MAX_VALUE];
            }
        }

    }
    
    public void set (int row, int column, byte value) {
        int[] indices = getIndices(row, column);
        myData[indices[0]][indices[1]] = value;
    }
    
    public byte get (int row, int column) {
        int[] indices = getIndices(row, column);
        return myData[indices[0]][indices[1]];
    }
    
    private int[] getIndices(int row, int column) {
        int[] result = new int[2];
        result[0] = row / myNumRowsPerSingleDimArray;
        result[1] = (row % myNumRowsPerSingleDimArray) * myNumColumns + column;
        return result;
    }
}
