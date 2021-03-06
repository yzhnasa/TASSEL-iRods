/*
 * TasselPrefs.java
 *
 * Created on August 5, 2007, 6:58 PM
 *
 */
package net.maizegenetics.prefs;

import java.util.HashMap;
import java.util.Map;
import java.util.prefs.Preferences;

/**
 *
 * @author terryc
 */
public class TasselPrefs {

    private static boolean PERSIST_PREFERENCES = false;
    private static Map<String, Object> TEMP_CACHED_VALUES = new HashMap<>();
    //
    // Top level preferences
    //
    public static final String TASSEL_TOP = "/tassel";
    public static final String TASSEL_SAVE_DIR = "saveDir";
    public static final String TASSEL_SAVE_DIR_DEFAULT = "";
    public static final String TASSEL_OPEN_DIR = "openDir";
    public static final String TASSEL_OPEN_DIR_DEFAULT = "";
    public static final String TASSEL_X_DIM = "xDimension";
    public static final int TASSEL_X_DIM_DEFAULT = -1;
    public static final String TASSEL_Y_DIM = "yDimension";
    public static final int TASSEL_Y_DIM_DEFAULT = -1;
    //
    // FilterAlignmentPlugin preferences
    //
    public static final String FILTER_ALIGN_PLUGIN_TOP = "/tassel/plugins/filterAlign";
    // Min. frequency for filtering sites.
    public static final String FILTER_ALIGN_PLUGIN_MIN_FREQ = "minFreq";
    public static final double FILTER_ALIGN_PLUGIN_MIN_FREQ_DEFAULT = 0.0;
    // Max. frequency for filtering sites.
    public static final String FILTER_ALIGN_PLUGIN_MAX_FREQ = "maxFreq";
    public static final double FILTER_ALIGN_PLUGIN_MAX_FREQ_DEFAULT = 1.0;
    // Min. frequency for filtering sites.
    public static final String FILTER_ALIGN_PLUGIN_MIN_COUNT = "minCount";
    public static final int FILTER_ALIGN_PLUGIN_MIN_COUNT_DEFAULT = 1;
    //
    // FilterTaxaPropertiesPlugin preferences
    //
    public static final String FILTER_TAXA_PROPS_PLUGIN_TOP = "/tassel/plugins/filterTaxaAlign";
    // Min. Not Missing Gametes Proportion
    public static final String FILTER_TAXA_PROPS_PLUGIN_MIN_NOT_MISSING = "minNotMissingFreq";
    public static final double FILTER_TAXA_PROPS_PLUGIN_MIN_NOT_MISSING_DEFAULT = 0.1;
    //Min. Heterozygotes Proportion
    public static final String FILTER_TAXA_PROPS_PLUGIN_MIN_HET = "minHetFreq";
    public static final double FILTER_TAXA_PROPS_PLUGIN_MIN_HET_DEFAULT = 0.0;
    //Max. Heterozygotes Proportion
    public static final String FILTER_TAXA_PROPS_PLUGIN_MAX_HET = "maxHetFreq";
    public static final double FILTER_TAXA_PROPS_PLUGIN_MAX_HET_DEFAULT = 1.0;
    //
    // Alignment preferences
    //
    public static final String ALIGNMENT_TOP = "/tassel/alignment";
    public static final String ALIGNMENT_RETAIN_RARE_ALLELES = "retainRareAlleles";
    public static final boolean ALIGNMENT_RETAIN_RARE_ALLELES_DEFAULT = false;

    /**
     * Creates a new instance of TasselPrefs
     */
    private TasselPrefs() {
    }

    public static boolean getPersistPreferences() {
        return PERSIST_PREFERENCES;
    }

    /**
     * Whether to Persist Preferences. Preference changes should be persisted
     * when executing GUI and set only temporarily from Command Line Flags. Also
     * getting preferences should use stored values when executing GUI. And
     * should use default values (if not temporarily set) when executing from
     * Command Line.
     *
     * @param persist whether to persist preferences
     */
    public static void setPersistPreferences(boolean persist) {
        PERSIST_PREFERENCES = persist;
    }

    public static String getPref(String path, String key, String def) {
        String pref = path + "/" + key;
        String result = (String) TEMP_CACHED_VALUES.get(pref);
        if (result != null) {
            return result;
        }
        if (PERSIST_PREFERENCES) {
            Preferences node = Preferences.userRoot();
            node = node.node(path);
            result = node.get(key, def);
        } else {
            result = def;
        }
        TEMP_CACHED_VALUES.put(pref, result);
        return result;
    }

    public static void putPref(String path, String key, String value) {
        String pref = path + "/" + key;
        TEMP_CACHED_VALUES.put(pref, value);
        if (PERSIST_PREFERENCES) {
            Preferences node = Preferences.userRoot();
            node = node.node(path);
            node.put(key, value);
        }
    }

    public static double getDoublePref(String path, String key, double def) {
        String pref = path + "/" + key;
        Double result = (Double) TEMP_CACHED_VALUES.get(pref);
        if (result != null) {
            return result;
        }
        if (PERSIST_PREFERENCES) {
            Preferences node = Preferences.userRoot();
            node = node.node(path);
            result = node.getDouble(key, def);
        } else {
            result = def;
        }
        TEMP_CACHED_VALUES.put(pref, result);
        return result;
    }

    public static void putDoublePref(String path, String key, double value) {
        String pref = path + "/" + key;
        TEMP_CACHED_VALUES.put(pref, value);
        if (PERSIST_PREFERENCES) {
            Preferences node = Preferences.userRoot();
            node = node.node(path);
            node.putDouble(key, value);
        }
    }

    public static int getIntPref(String path, String key, int def) {
        String pref = path + "/" + key;
        Integer result = (Integer) TEMP_CACHED_VALUES.get(pref);
        if (result != null) {
            return result;
        }
        if (PERSIST_PREFERENCES) {
            Preferences node = Preferences.userRoot();
            node = node.node(path);
            result = node.getInt(key, def);
        } else {
            result = def;
        }
        TEMP_CACHED_VALUES.put(pref, result);
        return result;
    }

    public static void putIntPref(String path, String key, int value) {
        String pref = path + "/" + key;
        TEMP_CACHED_VALUES.put(pref, value);
        if (PERSIST_PREFERENCES) {
            Preferences node = Preferences.userRoot();
            node = node.node(path);
            node.putInt(key, value);
        }
    }

    public static boolean getBooleanPref(String path, String key, boolean def) {
        String pref = path + "/" + key;
        Boolean result = (Boolean) TEMP_CACHED_VALUES.get(pref);
        if (result != null) {
            return result;
        }
        if (PERSIST_PREFERENCES) {
            Preferences node = Preferences.userRoot();
            node = node.node(path);
            result = node.getBoolean(key, def);
        } else {
            result = def;
        }
        TEMP_CACHED_VALUES.put(pref, result);
        return result;
    }

    public static void putBooleanPref(String path, String key, boolean value) {
        String pref = path + "/" + key;
        TEMP_CACHED_VALUES.put(pref, value);
        if (PERSIST_PREFERENCES) {
            Preferences node = Preferences.userRoot();
            node = node.node(path);
            node.putBoolean(key, value);
        }
    }

    //
    // Top level preferences
    //
    public static String getSaveDir() {
        return getPref(TASSEL_TOP, TASSEL_SAVE_DIR, TASSEL_SAVE_DIR_DEFAULT);
    }

    public static void putSaveDir(String value) {
        putPref(TASSEL_TOP, TASSEL_SAVE_DIR, value);
    }

    public static String getOpenDir() {
        return getPref(TASSEL_TOP, TASSEL_OPEN_DIR, TASSEL_OPEN_DIR_DEFAULT);
    }

    public static void putOpenDir(String value) {
        putPref(TASSEL_TOP, TASSEL_OPEN_DIR, value);
    }

    public static int getXDim() {
        return getIntPref(TASSEL_TOP, TASSEL_X_DIM, TASSEL_X_DIM_DEFAULT);
    }

    public static void putXDim(int value) {
        putIntPref(TASSEL_TOP, TASSEL_X_DIM, value);
    }

    public static int getYDim() {
        return getIntPref(TASSEL_TOP, TASSEL_Y_DIM, TASSEL_Y_DIM_DEFAULT);
    }

    public static void putYDim(int value) {
        putIntPref(TASSEL_TOP, TASSEL_Y_DIM, value);
    }

    //
    // FilterAlignmentPlugin preferences
    //
    public static double getFilterAlignPluginMinFreq() {
        return getDoublePref(FILTER_ALIGN_PLUGIN_TOP, FILTER_ALIGN_PLUGIN_MIN_FREQ, FILTER_ALIGN_PLUGIN_MIN_FREQ_DEFAULT);
    }

    public static void putFilterAlignPluginMinFreq(double value) {
        putDoublePref(FILTER_ALIGN_PLUGIN_TOP, FILTER_ALIGN_PLUGIN_MIN_FREQ, value);
    }

    public static double getFilterAlignPluginMaxFreq() {
        return getDoublePref(FILTER_ALIGN_PLUGIN_TOP, FILTER_ALIGN_PLUGIN_MAX_FREQ, FILTER_ALIGN_PLUGIN_MAX_FREQ_DEFAULT);
    }

    public static void putFilterAlignPluginMaxFreq(double value) {
        putDoublePref(FILTER_ALIGN_PLUGIN_TOP, FILTER_ALIGN_PLUGIN_MAX_FREQ, value);
    }

    public static int getFilterAlignPluginMinCount() {
        return getIntPref(FILTER_ALIGN_PLUGIN_TOP, FILTER_ALIGN_PLUGIN_MIN_COUNT, FILTER_ALIGN_PLUGIN_MIN_COUNT_DEFAULT);
    }

    public static void putFilterAlignPluginMinCount(int value) {
        putIntPref(FILTER_ALIGN_PLUGIN_TOP, FILTER_ALIGN_PLUGIN_MIN_COUNT, value);
    }

    //
    // FilterTaxaPropertiesPlugin preferences
    //
    public static double getFilterTaxaPropsMinNotMissingFreq() {
        return getDoublePref(FILTER_TAXA_PROPS_PLUGIN_TOP, FILTER_TAXA_PROPS_PLUGIN_MIN_NOT_MISSING, FILTER_TAXA_PROPS_PLUGIN_MIN_NOT_MISSING_DEFAULT);
    }

    public static void putFilterTaxaPropsMinNotMissingFreq(double value) {
        putDoublePref(FILTER_TAXA_PROPS_PLUGIN_TOP, FILTER_TAXA_PROPS_PLUGIN_MIN_NOT_MISSING, value);
    }

    public static double getFilterTaxaPropsMinHetFreq() {
        return getDoublePref(FILTER_TAXA_PROPS_PLUGIN_TOP, FILTER_TAXA_PROPS_PLUGIN_MIN_HET, FILTER_TAXA_PROPS_PLUGIN_MIN_HET_DEFAULT);
    }

    public static void putFilterTaxaPropsMinHetFreq(double value) {
        putDoublePref(FILTER_TAXA_PROPS_PLUGIN_TOP, FILTER_TAXA_PROPS_PLUGIN_MIN_HET, value);
    }

    public static double getFilterTaxaPropsMaxHetFreq() {
        return getDoublePref(FILTER_TAXA_PROPS_PLUGIN_TOP, FILTER_TAXA_PROPS_PLUGIN_MAX_HET, FILTER_TAXA_PROPS_PLUGIN_MAX_HET_DEFAULT);
    }

    public static void putFilterTaxaPropsMaxHetFreq(double value) {
        putDoublePref(FILTER_TAXA_PROPS_PLUGIN_TOP, FILTER_TAXA_PROPS_PLUGIN_MAX_HET, value);
    }

    //
    // Alignment preferences
    //
    public static boolean getAlignmentRetainRareAlleles() {
        return getBooleanPref(ALIGNMENT_TOP, ALIGNMENT_RETAIN_RARE_ALLELES, ALIGNMENT_RETAIN_RARE_ALLELES_DEFAULT);
    }

    public static void putAlignmentRetainRareAlleles(boolean value) {
        putBooleanPref(ALIGNMENT_TOP, ALIGNMENT_RETAIN_RARE_ALLELES, value);
    }
}
