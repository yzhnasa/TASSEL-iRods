// Identifier.java
//
// (c) 1999-2000 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)
package net.maizegenetics.pal.taxa;

import net.maizegenetics.pal.util.GeneralAnnotation;
import net.maizegenetics.pal.util.GeneralAnnotationUtils;
import net.maizegenetics.prefs.TasselPrefs;

import java.io.Serializable;
import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentMap;

/**
 * An identifier for some sampled data. This will most often be for example, the
 * accession number of a DNA sequence, or the taxonomic name that the sequence
 * represents, et cetera.
 *
 * The generally used class for defining a taxon. Contains its name, plus a
 * series of optional annotations. Use the builder to create immutable
 * instances.
 *
 * @author Terry Casstevens
 * @author Ed Buckler
 */
public class Taxon implements Serializable, Comparable, GeneralAnnotation {

    private static final long serialVersionUID = -7873729831795750538L;
    public static final String DELIMITER = ":";
    public static Taxon ANONYMOUS = new Taxon("");
    public static final String FatherKey="FATHER";
    public static final String MotherKey="MOTHER";
    public static final String PedigreeKey="PEDIGREE";
    public static final String SexKey="SEX";
    public static final String InbreedFKey="InbreedF";


    private final Map.Entry<String, String>[] myAnno;
    private final String myName;
    private final String[] myNameTokens;
    private final int hashCode;

    //since there are numerous redundant annotations and variant descriptions, this class use a annotation hash, so that
    //only the pointers are stored.  It takes a little longer to instantiate, but save 3-fold in memory.
    private static final ConcurrentMap<Map.Entry<String, String>,Map.Entry<String, String>> TAXON_ANNO_HASH = new ConcurrentHashMap<>(1_000_000);

    private static Map.Entry<String, String> getCanonicalAnnotation(String key, String value) {
        if (TAXON_ANNO_HASH.size() > 1_000_000) {
            TAXON_ANNO_HASH.clear();
        }
        Map.Entry<String, String> str= new AbstractMap.SimpleImmutableEntry(key,value);
        Map.Entry<String, String> canon = TAXON_ANNO_HASH.putIfAbsent(str, str);
        return (canon == null) ? str : canon;
    }

    public Taxon(String name) {
        this(new Builder(name));
    }

    private Taxon(Builder builder) {
        myName = builder.myTaxonFullName;
        myNameTokens = myName.split(DELIMITER);
        hashCode = myName.hashCode();
        //this looks crazy because it java doesn't support generic arrays
        myAnno=(Map.Entry<String, String>[])new Map.Entry<?,?>[builder.myAnnotations.size()];
        for (int i = 0; i < builder.myAnnotations.size(); i++) {
            myAnno[i]=builder.myAnnotations.get(i);
        }
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

    @Override
    public Object[] getAnnotation(String annoName) {
        return GeneralAnnotationUtils.getAnnotation(myAnno, annoName);
    }

    @Override
    public String[] getTextAnnotation(String annoName) {
        return GeneralAnnotationUtils.getTextAnnotation(myAnno,annoName);
    }

    @Override
    public double[] getQuantAnnotation(String annoName) {
        return GeneralAnnotationUtils.getQuantAnnotation(myAnno,annoName);
    }


    @Override
    public String getConsensusAnnotation(String annoName) {
        return GeneralAnnotationUtils.getConsensusAnnotation(myAnno,annoName);
    }

    @Override
    public double getAverageAnnotation(String annoName) {
        return GeneralAnnotationUtils.getAverageAnnotation(myAnno,annoName);
    }

    /**
     * A builder for creating immutable Taxon instances.
     * <p> Example:
     * <pre>   {@code
     * Taxon cp= new Taxon.Builder("Z001E0001:Line:mays:Zea")
     *   .inbreedF(0.99)
     *   .parents("B73","B97")
     *   .pedigree("(B73xB97)S6I1")
     *   .addAnno("Group","Dent")
     *   .build();}</pre>
     * <p>This would create an Taxon.
     */
    public static class Builder {

        // Required parameters
        private final String myTaxonFullName;
        private ArrayList<Map.Entry<String, String>> myAnnotations=new ArrayList<>(0);

        /**
         * Constructor for Builder, requires a Taxon object
         * @param aTaxon taxon object
         */
        public Builder(Taxon aTaxon) {
            myTaxonFullName = aTaxon.getFullName();
        }

        /**
         * Constructor for Builder, requires a Taxon name
         * @param aTaxonName name of the taxon
         */
        public Builder(String aTaxonName) {
            myTaxonFullName = aTaxonName;
        }

        /**Add non-standard annotation*/
        public Builder addAnno(String key, String value) {
            Map.Entry<String, String> ent=getCanonicalAnnotation(key,value);
            myAnnotations.add(ent);
            return this;
        }
        /**Add non-standard annotation*/
        public Builder addAnno(String key, Number value) {
            Map.Entry<String, String> ent=getCanonicalAnnotation(key,value.toString());
            myAnnotations.add(ent);
            return this;
        }

        /** Set sex: 0=both, 1=female, 2=male (default=0 Both) */
        public Builder sex(byte val) {return addAnno(SexKey,val);}

        /**  Set inbreeding coefficient (default=Float.NaN)*/
        public Builder inbreedF(float val) {return addAnno(InbreedFKey,val);}

        /** Set text definition of parents (default=null)*/
        public Builder parents(String mom, String dad) {
            addAnno(MotherKey,mom);
            return addAnno(FatherKey,dad);
        }

        /** Set text definition of pedigree (default=null) */
        public Builder pedigree(String val) {return addAnno(PedigreeKey,val);}


        public Taxon build() {
            return new Taxon(this);
        }
    }
}
