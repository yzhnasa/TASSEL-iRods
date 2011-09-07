// SimpleDataType.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)
// Known bugs and limitations:
package net.maizegenetics.pal.datatype;

/**
 */
public abstract class SimpleDataType implements DataType {

    private static final long serialVersionUID = 7902613264354545217L;
    private int myMaxFormatedStringLength = 1;


    public boolean isUnknownChar(final char c) {
        return isUnknownState(getState(c));
    }

    public boolean isUnknownState(final int state) {
        return ((state >= getNumStates()) || (state < 0));
    }

    public String toString() {
        return getDescription();
    }

    public void setMaxFormatedStringLength(int length) {
        myMaxFormatedStringLength = length;
    }

    public int getMaxFormatedStringLength() {
        return myMaxFormatedStringLength;
    }
}
