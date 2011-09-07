// Identifier.java
//
// (c) 1999-2000 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)
package net.maizegenetics.pal.ids;

import java.io.Serializable;

import net.maizegenetics.prefs.TasselPrefs;

/**
 * An identifier for some sampled data. This will most often be
 * for example, the accession number of a DNA sequence, or the
 * taxonomic name that the sequence represents, et cetera.
 *
 * @author terry
 */
public class Identifier implements Serializable, Comparable {

    private static final long serialVersionUID = -7873729831795750538L;
    public static final String DELIMITER = ":";
    private final String myName;
    public static Identifier ANONYMOUS = new Identifier("");

    public Identifier(String name) {
        myName = name;
    }

    public String toString() {
        return getName();
    }

    // implements Comparable interface
    public int compareTo(Object c) {
        if (c instanceof Identifier) {
            return compareTo(((Identifier) c).getFullName());
        } else {
            throw new ClassCastException();
        }
    }

    public boolean equals(Object c) {

        if (c instanceof Identifier) {
            if (TasselPrefs.getIDJoinStrict()) {
                return getFullName().equals(((Identifier) c).getFullName());
            } else {
                return compareTo(c) == 0;
            }
        } else if (c instanceof String) {
            if (TasselPrefs.getIDJoinStrict()) {
                return getFullName().equals(((String) c));
            } else {
                return compareTo((String) c) == 0;
            }
        } else {
            return false;
        }

    }

    public int compareTo(String c) {

        String[] first = myName.split(DELIMITER);
        String[] second = c.split(DELIMITER);
        int count = first.length < second.length ? first.length : second.length;
        for (int i = 0; i < count; i++) {
            int current = first[i].compareTo(second[i]);
            if (current != 0) {
                return current;
            }
        }

        if (TasselPrefs.getIDJoinStrict()) {
            if (first.length < second.length) {
                return -1;
            } else if (second.length < first.length) {
                return 1;
            } else {
                return 0;
            }
        } else {
            return 0;
        }

    }

    public String getName() {
        return getNameLevel(0);
    }

    public String getFullName() {
        return myName;
    }

    public String getFullName(String delimiter) {
        if (delimiter.equals(DELIMITER)) {
            return myName;
        }
        return myName.replaceAll(DELIMITER, delimiter);
    }

    /**
     * Returns requested level of name starting at index 0.
     * 0 will generally be most specific classification.
     *
     * @param index
     * @return Specified level.
     */
    public String getNameLevel(int index) {
        String[] temp = myName.split(DELIMITER);
        if (index < temp.length) {
            return temp[index];
        }
        return null;
    }

    /**
     * Returns name up to specified level (not including
     * specified level.  Levels start at index 0.
     *
     * @param index
     * @return name up to specified level exclusive.
     */
    public String getNameToLevel(int index) {
        return getNameToLevel(index, DELIMITER);
    }

    public String getNameToLevel(int index, String delimiter) {

        String[] temp = myName.split(DELIMITER);
        int upto = 0;
        if (index > temp.length) {
            upto = temp.length;
        } else {
            upto = index;
        }
        if (upto == 0) {
            return null;
        }

        StringBuilder builder = new StringBuilder();
        for (int i = 0; i < upto; i++) {
            if (i != 0) {
                builder.append(delimiter);
            }
            builder.append(temp[i]);
        }

        return builder.toString();
    }

    /**
     * Returns number of name levels.
     *
     * @return number of name levels.
     */
    public int getNumNameLevels() {
        return myName.split(DELIMITER).length;
    }
}

