// SimpleIdGroup.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)
package net.maizegenetics.pal.ids;

import java.io.Serializable;

import java.util.HashMap;


/**
 * Default implementation of IdGroup interface.
 * Memory-inefficient to allow fast whichIdNumber calls.
 *
 * @version $Id: SimpleIdGroup.java,v 1.6 2009/07/07 16:19:37 tcasstevens Exp $
 *
 * @author Alexei Drummond
 */
public class SimpleIdGroup implements IdGroup, Serializable {

    private Identifier[] ids;
    private HashMap myIndices;

    //
    // Serialization code
    //
    private static final long serialVersionUID = -4266575329980153075L;

    //serialver -classpath ./classes net.maizegenetics.pal.ids.SimpleIdGroup
    private void writeObject(java.io.ObjectOutputStream out) throws java.io.IOException {
        out.writeByte(1); //Version number
        out.writeObject(ids);
    }

    private void readObject(java.io.ObjectInputStream in) throws java.io.IOException, ClassNotFoundException {
        byte version = in.readByte();
        switch (version) {
            default: {
                ids = (Identifier[]) in.readObject();
                myIndices = new HashMap(ids.length);
                for (int i = 0; i < ids.length; i++) {
                    myIndices.put(ids[i].getFullName(), new Integer(i));
                }
                break;
            }
        }
    }

    /**
     * Constructor taking the size of the group.
     */
    public SimpleIdGroup(int size) {
        this(size, false);
    }

    /**
     * Constructor taking an array of strings.
     */
    public SimpleIdGroup(String[] labels) {
        this(labels.length);
        for (int i = 0; i < labels.length; i++) {
            setIdentifier(i, new Identifier(labels[i]));
        }
    }

    /**
     * Constructor taking the size of the group.
     * @param size - the number of ids
     * @param createIDs - if true creates default Identifiers.
     * Otherwise leaves blank (for user to fill in)
     */
    public SimpleIdGroup(int size, boolean createIDs) {

        ids = new Identifier[size];
        myIndices = new HashMap(size);
        if (createIDs) {
            for (int i = 0; i < size; i++) {
                setIdentifier(i, new Identifier("" + i));
            }
        }
    }

    /**
     * Constructor taking an array of identifiers.
     */
    public SimpleIdGroup(Identifier[] id) {
        this(id.length);
        for (int i = 0; i < id.length; i++) {
            setIdentifier(i, id[i]);
        }
    }

    /**
     * Constructor taking two separate id groups and merging them.
     */
    public SimpleIdGroup(IdGroup a, IdGroup b) {
        this(a.getIdCount() + b.getIdCount());

        for (int i = 0; i < a.getIdCount(); i++) {
            setIdentifier(i, a.getIdentifier(i));
        }
        for (int i = 0; i < b.getIdCount(); i++) {
            setIdentifier(i + a.getIdCount(), b.getIdentifier(i));
        }
    }

    /**
     * Impersonating Constructor.
     */
    public static SimpleIdGroup getInstance(IdGroup a) {
        Identifier[] ids = new Identifier[a.getIdCount()];
        for (int i = 0; i < a.getIdCount(); i++) {
            ids[i] = a.getIdentifier(i);
        }
        return new SimpleIdGroup(ids);
    }

    /**
     * Impersonating Constructor.
     * @param toIgnore - will ignore the identifier at the index specified by toIgnore
     */
    public SimpleIdGroup(IdGroup a, int toIgnore) {
        this((toIgnore < 0 || toIgnore > a.getIdCount() ? a.getIdCount() : a.getIdCount() - 1));
        int index = 0;
        for (int i = 0; i < a.getIdCount(); i++) {
            if (i != toIgnore) {
                setIdentifier(index++, a.getIdentifier(i));
            }
        }
    }

    /**
     * Returns the number of identifiers in this group
     */
    public int getIdCount() {
        return ids.length;
    }

    /**
     * Returns the ith identifier.
     */
    public Identifier getIdentifier(int i) {
        return ids[i];
    }

    /**
     * Convenience method to return the name of identifier i
     */
    public final String getName(int i) {
        return ids[i].getName();
    }

    /**
     * Sets the ith identifier.
     */
    public void setIdentifier(int i, Identifier id) {
        ids[i] = id;
        myIndices.put(id.getFullName(), new Integer(i));
    }

    /**
     * Return index of identifier with name or -1 if not found
     */
    public int whichIdNumber(String name) {

        Integer index = (Integer) myIndices.get(name);
        if (index != null) {
            return index.intValue();
        }

        for (int i = 0, n = ids.length; i < n; i++) {
            if (ids[i].equals(name)) {
                return i;
            }
        }

        return -1;

    }

    /**
     * Returns a string representation of this IdGroup in the form of
     * a bracketed list.
     */
    public String toString() {

        StringBuffer sb = new StringBuffer();
        sb.append("[ ");
        for (int i = 0; i < getIdCount(); i++) {
            sb.append(getIdentifier(i) + " ");
        }
        sb.append("]");
        return new String(sb);
    }

    public int whichIdNumber(Identifier id) {
        return whichIdNumber(id.getFullName());
    }
}

