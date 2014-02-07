package net.maizegenetics.trait;

import net.maizegenetics.trait.MarkerPhenotypeAdapter;
import java.util.ArrayList;

public class MarkerPhenotypeAdapterUtils {

    //to prevent istantiation of this utility class
    private MarkerPhenotypeAdapterUtils() {
    }

    /**
     * Updates missing1 with missing2. Sets missing1 = missing1 || missing2.
     * @param missing1 an array of booleans
     * @param missing2 an array of booleans
     */
    public static void updateMissing(boolean[] missing1, boolean[] missing2) {
        int n = missing1.length;
        if (missing2.length != n) {
            throw new IllegalArgumentException("missing1 and missing2 must be the same length");
        }
        for (int i = 0; i < n; i++) {
            missing1[i] = missing1[i] || missing2[i];
        }
    }

    /**
     * Updates missing with data. Sets missing = missing || (data = NaN).
     * @param missing	an array of booleans
     * @param data		an array of doubles
     */
    public static void updateMissing(boolean[] missing, double[] data) {
        int n = missing.length;
        if (data.length != n) {
            throw new IllegalArgumentException("missing and data must be the same length");
        }
        for (int i = 0; i < n; i++) {
            missing[i] = missing[i] || Double.isNaN(data[i]);
        }
    }

    /**
     * Updates missing with labels. Sets missing = missing || isMissing(labels).
     * @param missing an array of booleans
     * @param labels an array of Strings
     */
    public static void updateMissing(boolean[] missing, String[] labels) {
        int n = missing.length;
        if (labels.length != n) {
            throw new IllegalArgumentException("missing and labels must be the same length");
        }
        for (int i = 0; i < n; i++) {
            missing[i] = missing[i] || isMissing(labels[i]);
        }
    }

    /**
     * @param mpa	a MarkerPhenotypeAdapter
     * @param phenotype	the phenotype for which the factor values will be returned
     * @param missing	the array that keeps track of which rows have one or more missing values
     * @return	an ArrayList of String[]. Each String[] is a list of factor values. The missing array is modified in place. Returns null if there are no covariates.
     */
    public static ArrayList<String[]> getFactorList(MarkerPhenotypeAdapter mpa, int phenotype, boolean[] missing) {
        int n = mpa.getNumberOfFactors();
        if (n == 0) {
            return null;
        }
        ArrayList<String[]> factorList = new ArrayList<String[]>(n);
        for (int i = 0; i < n; i++) {
            factorList.add(mpa.getFactorValues(phenotype, i));
            updateMissing(missing, mpa.getMissingFactors(phenotype, i));
        }
        return factorList;
    }

    /**
     * @param mpa	a MarkerPhenotypeAdapter
     * @param phenotype	the phenotype for which the factor values will be returned
     * @param missing	the array that keeps track of which rows have one or more missing values
     * @return	an ArrayList of double[]. Each double[] is a list of covariate values. The missing array is modified in place. Returns null if there are no covariates.
     */
    public static ArrayList<double[]> getCovariateList(MarkerPhenotypeAdapter mpa, int phenotype, boolean[] missing) {
        int n = mpa.getNumberOfCovariates();
        if (n == 0) {
            return null;
        }
        ArrayList<double[]> covariateList = new ArrayList<double[]>(n);
        for (int i = 0; i < n; i++) {
            covariateList.add(mpa.getCovariateValues(phenotype, i));
            updateMissing(missing, mpa.getMissingCovariates(phenotype, i));
        }
        return covariateList;
    }

    /**
     * @param doubleArray	an array of doubles
     * @return	an array of booleans equal to true if the corresponding value in doubleArray is NaN, false otherwise.
     */
    public static boolean[] whichAreMissing(double[] doubleArray) {
        int n = doubleArray.length;
        boolean[] missing = new boolean[n];
        for (int i = 0; i < n; i++) {
            missing[i] = Double.isNaN(doubleArray[i]);
        }
        return missing;
    }

    /**
     * Tests for missing values using the MarkerPhenotypeAdapterUtils isMissing() function.
     * @param stringArray	an array of String
     * @return	an array of booleans equal to true if the corresponding String equals the missing value, false otherwise.
     */
    public static boolean[] whichAreMissing(String[] stringArray) {
        int n = stringArray.length;
        boolean[] missing = new boolean[n];
        for (int i = 0; i < n; i++) {
            missing[i] = (isMissing(stringArray[i]));
        }
        return missing;
    }

    /**
     * Tests for missing values using the MarkerPhenotypeAdapterUtils isMissing() function.
     * @param stringArray	an array of Object
     * @return	an array of booleans equal to true if the corresponding Object equals the missing value, false otherwise.
     */
    public static boolean[] whichAreMissing(Object[] values) {
        int n = values.length;
        boolean[] missing = new boolean[n];
        if (values instanceof String[]) {
            String[] strvalues = (String[]) values;
            for (int i = 0; i < n; i++) {
                missing[i] = isMissing(strvalues[i]);
            }
        } else if (values instanceof Double[]) {
            Double[] dblvalues = (Double[]) values;
            for (int i = 0; i < n; i++) {
                missing[i] = isMissing(dblvalues[i]);
            }
        } else if (values instanceof Character[]) {
            Character[] chrvalues = (Character[]) values;
            for (int i = 0; i < n; i++) {
                missing[i] = isMissing(chrvalues[i]);
            }
        } else {
            for (int i = 0; i < n; i++) {
                missing[i] = (values[i] == null);
            }
        }
        return missing;
    }

    /**
     * @param values	an array of Objects
     * @return	true if any of the Objects is null or missing as tested by MarkerPhenotypeAdapterUtils isMissing() function, false otherwise.
     */
    public static boolean areAnyMissing(Object[] values) {
        if (values instanceof String[]) {
            String[] strvalues = (String[]) values;
            for (String str : strvalues) {
                if (isMissing(str)) {
                    return true;
                }
            }
            return false;
        } else if (values instanceof Double[]) {
            Double[] dblvalues = (Double[]) values;
            for (Double dbl : dblvalues) {
                if (isMissing(dbl)) {
                    return true;
                }
            }
            return false;
        } else if (values instanceof Character[]) {
            Character[] chrvalues = (Character[]) values;
            for (Character chr : chrvalues) {
                if (isMissing(chr)) {
                    return true;
                }
            }
            return false;
        } else {
            for (Object val : values) {
                if (val == null) {
                    return true;
                }
            }
            return false;
        }
    }

    /**
     * @param values	an array of doubles
     * @return	true if any of the doubles equals NaN, false otherwise.
     */
    public static boolean areAnyMissing(double[] values) {
        for (double val : values) {
            if (Double.isNaN(val)) {
                return true;
            }
        }
        return false;
    }

    /**
     * @param missing an array of booleans
     * @return	an array of the indices of non-missing elements, that is the elements of missing equal to false
     */
    public static int[] getNonMissingIndex(boolean[] missing) {
        int n = missing.length;
        int ntrue = 0;
        for (boolean miss : missing) {
            if (!miss) {
                ntrue++;
            }
        }
        int[] nonmiss = new int[ntrue];
        int count = 0;
        for (int i = 0; i < n; i++) {
            if (!missing[i]) {
                nonmiss[count++] = i;
            }
        }
        return nonmiss;
    }

    /**
     * @param missing	a boolean[] array
     * @return the number of elements of missing equal to false
     */
    public static int numberNotMissing(boolean[] missing) {
        int numberNotMissing = 0;
        for (boolean m : missing) {
            if (!m) {
                numberNotMissing++;
            }
        }
        return numberNotMissing;
    }

    /**
     * Tests for missing values. Returns true for a Double equal to Nan, a Character equal to '?', or a String equal to "?" or "N".
     * @param value	an Object
     * @return	true if equal to a missing value, false otherwise
     */
    public static boolean isMissing(Object value) {
        if (value instanceof Double) {
            return ((Double) value).isNaN();
        }
        if (value instanceof Character) {
            return value.equals('?');
        }
        if (value instanceof String) {
            String val = (String) value;
            if (val.equals("?")) {
                return true;
            }
            if (val.equals("N")) {
                return true;
            }
            return false;
        } else {
            return false;
        }
    }

    /**
     * @param value a Double
     * @return true if the Double equals NaN.
     */
    public static boolean isMissing(Double value) {
        return value.isNaN();
    }

    /**
     * @param value a String
     * @return true if the String equals "N" or "?".
     */
    public static boolean isMissing(String value) {
        return (value.equals("?") || value.equals("N"));
    }

    /**
     * @param value a Character
     * @return true if the Character equals 'N'.
     */
    public static boolean isMissing(Character value) {
        return (value.equals('?') || value.equals('N'));
    }
}
