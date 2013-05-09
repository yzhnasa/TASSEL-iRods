/*
 * TasselPrefs.java
 *
 * Created on August 5, 2007, 6:58 PM
 *
 */
package net.maizegenetics.prefs;

import java.util.prefs.Preferences;

/**
 *
 * @author terryc
 */
public class TasselPrefs {

    //
    // Top level preferences
    //
    public static final String TASSEL_TOP = "/tassel";
    public static final String TASSEL_SAVE_DIR = "saveDir";
    public static final String TASSEL_SAVE_DIR_DEFAULT = "";
    public static final String TASSEL_OPEN_DIR = "openDir";
    public static final String TASSEL_OPEN_DIR_DEFAULT = "";
    public static final String TASSEL_FONT_METRICS_CHAR_WIDTH = "fontMetricsCharWidth";
    public static final int TASSEL_FONT_METRICS_CHAR_WIDTH_DEFAULT = 9;

    public static enum TASSEL_IDENTIFIER_JOIN_TYPES {

        /**
         * This means taxa names must match exactly all levels. This is the same
         * as true use to be.
         */
        Strict,
        /**
         * This means taxa names must match all levels up to the least specific
         * taxa. This is the same as false use to be.
         */
        NonStrict,
        /**
         * This means taxa names must match specified number of levels.
         */
        NumLevels
    };
    public static final String TASSEL_IDENTIFIER_JOIN_STRICT = "idJoinStrict";
    public static final TASSEL_IDENTIFIER_JOIN_TYPES TASSEL_IDENTIFIER_JOIN_STRICT_DEFAULT = TASSEL_IDENTIFIER_JOIN_TYPES.NonStrict;
    private static TASSEL_IDENTIFIER_JOIN_TYPES TASSEL_IDENTIFIER_JOIN_STRICT_VALUE = initIDJoin();
    public static final String TASSEL_IDENTIFIER_JOIN_NUM_LEVELS = "idJoinNumLevels";
    public static final int TASSEL_IDENTIFIER_JOIN_NUM_LEVELS_DEFAULT = 1;
    private static int TASSEL_IDENTIFIER_JOIN_NUM_LEVELS_VALUE = getIntPref(TASSEL_TOP, TASSEL_IDENTIFIER_JOIN_NUM_LEVELS, TASSEL_IDENTIFIER_JOIN_NUM_LEVELS_DEFAULT);
    //
    // FilterAlignmentPlugin preferences
    //
    public static final String FILTER_ALIGN_PLUGIN_TOP = "/tassel/plugins/filterAlign";
    // Min. frequency for filtering sites.
    public static final String FILTER_ALIGN_PLUGIN_MIN_FREQ = "minFreq";
    public static final double FILTER_ALIGN_PLUGIN_MIN_FREQ_DEFAULT = 0.01;
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
    //Max. Heterozygotes Proportion
    public static final String FILTER_TAXA_PROPS_PLUGIN_MAX_HET = "maxHetFreq";
    public static final double FILTER_TAXA_PROPS_PLUGIN_MAX_HET_DEFAULT = 0.1;
    //
    // Alignment preferences
    //
    public static final String ALIGNMENT_TOP = "/tassel/alignment";
    public static final String ALIGNMENT_MAX_ALLELES_TO_RETAIN = "maxAllelesToRetain";
    public static final int ALIGNMENT_MAX_ALLELES_TO_RETAIN_DEFAULT = 2;
    public static final String ALIGNMENT_RETAIN_RARE_ALLELES = "retainRareAlleles";
    public static final boolean ALIGNMENT_RETAIN_RARE_ALLELES_DEFAULT = true;

    /**
     * Creates a new instance of TasselPrefs
     */
    private TasselPrefs() {
    }

    public static String getPref(String path, String key, String def) {
        Preferences node = Preferences.userRoot();
        node = node.node(path);
        return node.get(key, def);
    }

    public static void putPref(String path, String key, String value) {
        Preferences node = Preferences.userRoot();
        node = node.node(path);
        node.put(key, value);
    }

    public static double getDoublePref(String path, String key, double def) {
        Preferences node = Preferences.userRoot();
        node = node.node(path);
        return node.getDouble(key, def);
    }

    public static void putDoublePref(String path, String key, double value) {
        Preferences node = Preferences.userRoot();
        node = node.node(path);
        node.putDouble(key, value);
    }

    public static int getIntPref(String path, String key, int def) {
        Preferences node = Preferences.userRoot();
        node = node.node(path);
        return node.getInt(key, def);
    }

    public static void putIntPref(String path, String key, int value) {
        Preferences node = Preferences.userRoot();
        node = node.node(path);
        node.putInt(key, value);
    }

    public static boolean getBooleanPref(String path, String key, boolean def) {
        Preferences node = Preferences.userRoot();
        node = node.node(path);
        return node.getBoolean(key, def);
    }

    public static void putBooleanPref(String path, String key, boolean value) {
        Preferences node = Preferences.userRoot();
        node = node.node(path);
        node.putBoolean(key, value);
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

    public static int getFontMetricsCharWidth() {
        return getIntPref(TASSEL_TOP, TASSEL_FONT_METRICS_CHAR_WIDTH, TASSEL_FONT_METRICS_CHAR_WIDTH_DEFAULT);
    }

    public static void putFontMetricsCharWidth(int value) {
        putIntPref(TASSEL_TOP, TASSEL_FONT_METRICS_CHAR_WIDTH, value);
    }

    private static TASSEL_IDENTIFIER_JOIN_TYPES initIDJoin() {

        String value = getPref(TASSEL_TOP, TASSEL_IDENTIFIER_JOIN_STRICT, null);
        if (value == null) {
            value = TASSEL_IDENTIFIER_JOIN_STRICT_DEFAULT.toString();
        } else if (value.equalsIgnoreCase("true")) {
            value = TASSEL_IDENTIFIER_JOIN_TYPES.Strict.toString();
        } else if (value.equalsIgnoreCase("false")) {
            value = TASSEL_IDENTIFIER_JOIN_TYPES.NonStrict.toString();
        }

        TASSEL_IDENTIFIER_JOIN_TYPES result = null;
        try {
            result = TASSEL_IDENTIFIER_JOIN_TYPES.valueOf(value);
        } catch (Exception e) {
            result = TASSEL_IDENTIFIER_JOIN_STRICT_DEFAULT;
        }

        putPref(TASSEL_TOP, TASSEL_IDENTIFIER_JOIN_STRICT, result.toString());

        return result;

    }

    public static TASSEL_IDENTIFIER_JOIN_TYPES getIDJoinStrict() {
        // This can be called many times, so to improve performance
        // this will return value without executing system call.
        return TASSEL_IDENTIFIER_JOIN_STRICT_VALUE;
    }

    public static void putIDJoinStrict(TASSEL_IDENTIFIER_JOIN_TYPES type) {
        TASSEL_IDENTIFIER_JOIN_STRICT_VALUE = type;
        putPref(TASSEL_TOP, TASSEL_IDENTIFIER_JOIN_STRICT, type.toString());
    }

    /**
     * Sets ID Join type only for current program execution.
     *
     * @param type type
     */
    public static void putTempIDJoinStrict(TASSEL_IDENTIFIER_JOIN_TYPES type) {
        TASSEL_IDENTIFIER_JOIN_STRICT_VALUE = type;
    }

    public static int getIDJoinNumLevels() {
        return TASSEL_IDENTIFIER_JOIN_NUM_LEVELS_VALUE;
    }

    public static void putIDJoinNumLevels(int value) {
        TASSEL_IDENTIFIER_JOIN_NUM_LEVELS_VALUE = value;
        putIntPref(TASSEL_TOP, TASSEL_IDENTIFIER_JOIN_NUM_LEVELS, value);
    }

    public static void putTempIDJoinNumLevels(int value) {
        TASSEL_IDENTIFIER_JOIN_NUM_LEVELS_VALUE = value;
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

    public static double getFilterTaxaPropsMaxHetFreq() {
        return getDoublePref(FILTER_TAXA_PROPS_PLUGIN_TOP, FILTER_TAXA_PROPS_PLUGIN_MAX_HET, FILTER_TAXA_PROPS_PLUGIN_MAX_HET_DEFAULT);
    }

    public static void putFilterTaxaPropsMaxHetFreq(double value) {
        putDoublePref(FILTER_TAXA_PROPS_PLUGIN_TOP, FILTER_TAXA_PROPS_PLUGIN_MAX_HET, value);
    }

    //
    // Alignment preferences
    //
    public static int getAlignmentMaxAllelesToRetain() {
        return getIntPref(ALIGNMENT_TOP, ALIGNMENT_MAX_ALLELES_TO_RETAIN, ALIGNMENT_MAX_ALLELES_TO_RETAIN_DEFAULT);
    }

    public static void putAlignmentMaxAllelesToRetain(int value) {
        putIntPref(ALIGNMENT_TOP, ALIGNMENT_MAX_ALLELES_TO_RETAIN, value);
    }

    public static boolean getAlignmentRetainRareAlleles() {
        return getBooleanPref(ALIGNMENT_TOP, ALIGNMENT_RETAIN_RARE_ALLELES, ALIGNMENT_RETAIN_RARE_ALLELES_DEFAULT);
    }

    public static void putAlignmentRetainRareAlleles(boolean value) {
        putBooleanPref(ALIGNMENT_TOP, ALIGNMENT_RETAIN_RARE_ALLELES, value);
    }
}
