// Identifier.java
//
// (c) 1999-2000 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)
package net.maizegenetics.pal.taxa;

import com.google.common.collect.ImmutableMultimap;
import net.maizegenetics.pal.util.AbstractAnnotation;
import net.maizegenetics.pal.util.GeneralAnnotation;
import net.maizegenetics.prefs.TasselPrefs;

import java.io.Serializable;

/**
 * An identifier for some sampled data. This will most often be for example, the
 * accession number of a DNA sequence, or the taxonomic name that the sequence
 * represents, et cetera.
 *
 * The generally used class for defining a taxon. Contains its name, plus a
 * series of optional annotations. Use the builder to create immutable
 * instances.
 * //TODO has a very inefficient method for storing Strings. String and/or annotation pools need to be used.
 *
 * @author Terry Casstevens
 * @author Ed Buckler
 */
public class Taxon implements Serializable, Comparable, GeneralAnnotation {

    private static final long serialVersionUID = -7873729831795750538L;
    public static final String DELIMITER = ":";
    //TODO need to change the annotation approach to the same GeneralPosition with pools
    protected final GeneralAnnotation myGA;
    protected final String myParent1;  // generally female
    protected final String myParent2;  // generally male
    protected final float myInbreedF; // inbreeding coefficient
    protected final byte mySex;  // 0=both, 1=female, 2=male
    protected final String myPedigree;
    private final String myName;
    private final String[] myNameTokens;
    public static Taxon ANONYMOUS = new Taxon("");
    private final int hashCode;

    public Taxon(String name) {
        this(new Builder(name));
    }

    private Taxon(Builder builder) {
        myName = builder.myTaxonFullName;
        myNameTokens = myName.split(DELIMITER);
        hashCode = myName.hashCode();
        this.mySex = builder.sex;
        this.myInbreedF = builder.inbreedF;
        this.myParent1 = builder.parent1;
        this.myPedigree = builder.pedigree;
        this.myGA = new AbstractAnnotation(builder.myAnnoMap);
        this.myParent2 = builder.parent2;
    }

    public static Taxon getMergedInstance(Taxon id1, Taxon id2) {
        String[] first = id1.getFullNameTokens();
        String[] second = id2.getFullNameTokens();
        int count = Math.min(first.length, second.length);
        for (int i = 0; i < count; i++) {
            if (!first[i].equals(second[i])) {
                StringBuilder builder = new StringBuilder();
                for (int x = 0; x < i; x++) {
                    if (x != 0) {
                        builder.append(DELIMITER);
                    }
                    builder.append(first[x]);
                    return new Taxon(builder.toString());
                }
            }
        }
        return id1;
    }

    public String toString() {
        return getName();
    }

    // implements Comparable interface
    @Override
    public int compareTo(Object c) {
        if (this == c) {
            return 0;
        } else if (c instanceof Taxon) {
            return compareTo(((Taxon) c).getFullNameTokens());
        } else {
            throw new ClassCastException();
        }
    }

    @Override
    public boolean equals(Object c) {

        if (this == c) {
            return true;
        }

        TasselPrefs.TASSEL_IDENTIFIER_JOIN_TYPES type = TasselPrefs.getIDJoinStrict();

        if (c instanceof Taxon) {
            if (type == TasselPrefs.TASSEL_IDENTIFIER_JOIN_TYPES.Strict) {
                return getFullName().equals(((Taxon) c).getFullName());
            } else {
                return compareTo(c) == 0;
            }
        } else if (c instanceof String) {
            if (type == TasselPrefs.TASSEL_IDENTIFIER_JOIN_TYPES.Strict) {
                return getFullName().equals(((String) c));
            } else {
                return compareTo((String) c) == 0;
            }
        } else {
            return false;
        }

    }

    public int compareTo(String c) {
        return compareTo(c.split(DELIMITER));
    }

    public int compareTo(String[] fullNameTokens) {

        TasselPrefs.TASSEL_IDENTIFIER_JOIN_TYPES type = TasselPrefs.getIDJoinStrict();

        int count = Math.min(myNameTokens.length, fullNameTokens.length);
        if (type == TasselPrefs.TASSEL_IDENTIFIER_JOIN_TYPES.NumLevels) {
            count = Math.min(count, TasselPrefs.getIDJoinNumLevels());
        }
        for (int i = 0; i < count; i++) {
            int current = myNameTokens[i].compareTo(fullNameTokens[i]);
            if (current != 0) {
                return current;
            }
        }

        if (type == TasselPrefs.TASSEL_IDENTIFIER_JOIN_TYPES.Strict) {
            if (myNameTokens.length < fullNameTokens.length) {
                return -1;
            } else if (fullNameTokens.length < myNameTokens.length) {
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

    public String[] getFullNameTokens() {
        return myNameTokens;
    }

    /**
     * Returns requested level of name starting at index 0. 0 will generally be
     * most specific classification.
     *
     * @param index
     * @return Specified level.
     */
    public String getNameLevel(int index) {
        if (index < myNameTokens.length) {
            return myNameTokens[index];
        }
        return null;
    }

    /**
     * Returns name up to specified level (not including specified level. Levels
     * start at index 0.
     *
     * @param index
     * @return name up to specified level exclusive.
     */
    public String getNameToLevel(int index) {
        return getNameToLevel(index, DELIMITER);
    }

    public String getNameToLevel(int index, String delimiter) {

        int upto = 0;
        if (index > myNameTokens.length) {
            upto = myNameTokens.length;
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
            builder.append(myNameTokens[i]);
        }

        return builder.toString();
    }

    /**
     * Returns number of name levels.
     *
     * @return number of name levels.
     */
    public int getNumNameLevels() {
        return myNameTokens.length;
    }

    @Override
    public int hashCode() {
        return hashCode;
    }

    public String getParent1() {
        return myParent1;
    }

    public String getParent2() {
        return myParent2;
    }

    public float getInbreedF() {
        return myInbreedF;
    }

    public byte getSex() {
        return mySex;
    }

    public String getPedigree() {
        return myPedigree;
    }

    @Override
    public Object[] getAnnotation(String annoName) {
        return myGA.getAnnotation(annoName);
//        switch (annoName) {  //TODO: uncomment once in Java 7
//            case "myParent1":return myLocus;
//            case "myParent2":return myPosition;
//            case "myInbreedF":return myCM;
//            case "mySex":return myStrand;
//            case "pedigree":return mySNPID;
//        }
//       }
    }

    @Override
    public String[] getTextAnnotation(String annoName) {
        return myGA.getTextAnnotation(annoName);
    }

    @Override
    public double[] getQuantAnnotation(String annoName) {
        return myGA.getQuantAnnotation(annoName);
    }

    @Override
    public String getConsensusAnnotation(String annoName) {
        return myGA.getConsensusAnnotation(annoName);
    }

    @Override
    public double getAverageAnnotation(String annoName) {
        return myGA.getAverageAnnotation(annoName);
    }

    /**
     * A builder for creating immutable Taxon instances.
     * <p> Example:
     * <pre>   {@code
     * Taxon cp= new Taxon.Builder("Z001E0001:Line:mays:Zea")
     *   .inbreedF(0.99)
     *   .parents("B73","B97")
     *   .pedigree("(B73xB97)S6I1")
     *   .build();}</pre>
     * <p>This would create an Taxon.
     */
    public static class Builder {

        // Required parameters
        private final String myTaxonFullName;
        // Optional parameters - initialized to default values
        private String parent1 = null;  //generally female
        private String parent2 = null;  //generally male
        private float inbreedF = Float.NaN;
        private byte sex = 0;  //0=both, 1=female, 2=male
        private String pedigree = null;
        private ImmutableMultimap.Builder<String, Object> myAnnoMapBld = null;
        private ImmutableMultimap<String, Object> myAnnoMap = null;

        /**
         * Constructor for Builder, requires a Taxon object
         *
         * @param aTaxon taxon object
         */
        public Builder(Taxon aTaxon) {
            myTaxonFullName = aTaxon.getFullName();
        }

        /**
         * Constructor for Builder, requires a Taxon name
         *
         * @param aTaxonName name of the taxon
         */
        public Builder(String aTaxonName) {
            myTaxonFullName = aTaxonName;
        }

        /**
         * Set sex: 0=both, 1=female, 2=male (default=0 Both)
         */
        public Builder sex(byte val) {
            sex = val;
            return this;
        }

        /**
         * Set inbreeding coefficient (default=Float.NaN)
         */
        public Builder inbreedF(float val) {
            inbreedF = val;
            return this;
        }

        /**
         * Set text definition of parents (default=null)
         */
        public Builder parents(String mom, String dad) {
            parent1 = mom;
            parent2 = dad;
            return this;
        }

        /**
         * Set text definition of pedigree (default=null)
         */
        public Builder pedigree(String val) {
            pedigree = val;
            return this;
        }

        /**
         * Add non-standard annotation
         */
        public Builder addAnno(String key, String value) {
            if (myAnnoMapBld == null) {
                myAnnoMapBld = new ImmutableMultimap.Builder();
            }
            myAnnoMapBld.put(key, value);
            return this;
        }

        /**
         * Add non-standard annotation
         */
        public Builder addAnno(String key, Number value) {
            if (myAnnoMapBld == null) {
                myAnnoMapBld = new ImmutableMultimap.Builder();
            }
            myAnnoMapBld.put(key, value);
            return this;
        }

        public Taxon build() {
            if (myAnnoMapBld != null) {
                myAnnoMap = myAnnoMapBld.build();
            }
            return new Taxon(this);
        }
    }
}
