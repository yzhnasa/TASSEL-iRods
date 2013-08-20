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
        if (rows * columns > Integer.MAX_VALUE) {
            return new SuperByteMatrixMultiple(rows, columns);
        } else {
            return new SuperByteMatrixSingle(rows, columns);
        }
    }
}
