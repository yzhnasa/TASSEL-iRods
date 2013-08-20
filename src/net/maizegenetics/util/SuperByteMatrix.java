/*
 *  SuperByteMatrix
 */
package net.maizegenetics.util;

/**
 *
 * @author Terry Casstevens
 */
public interface SuperByteMatrix {

    public void set(int row, int column, byte value);

    public byte get(int row, int column);
}
