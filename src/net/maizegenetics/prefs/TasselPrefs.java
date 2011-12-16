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
    public static final String TASSEL_IDENTIFIER_JOIN_STRICT = "idJoinStrict";
    public static final boolean TASSEL_IDENTIFIER_JOIN_STRICT_DEFAULT = false;
    private static boolean TASSEL_IDENTIFIER_JOIN_STRICT_VALUE = getIDJoinStrict();
    //
    // FilterAlignmentPlugin preferences
    //
    public static final String FILTER_ALIGN_PLUGIN_TOP = "/tassel/plugins/filterAlign";
    // Min. frequency for filtering sites.
    public static final String FILTER_ALIGN_PLUGIN_MIN_FREQ = "minFreq";
    public static final double FILTER_ALIGN_PLUGIN_MIN_FREQ_DEFAULT = 0.01;
    // Min. frequency for filtering sites.
    public static final String FILTER_ALIGN_PLUGIN_MIN_COUNT = "minCount";
    public static final int FILTER_ALIGN_PLUGIN_MIN_COUNT_DEFAULT = 1;
    //
    // AlignmentPlugin Preferences
    //
    public static final String ALIGN_PLUGIN_TOP = "/tassel/plugins/alignment";
    public static final String ALIGN_PLUGIN_QUALSCORE_LOWRANGE_CUTOFF = "qualityScoreLowRangeCuttoff";
    public static final int ALIGN_PLUGIN_QUALSCORE_LOWRANGE_CUTOFF_DEFAULT = 0;
    public static final String ALIGN_PLUGIN_QUALSCORE_HIGHRANGE_CUTOFF = "qualityScoreHighRangeCuttoff";
    public static final int ALIGN_PLUGIN_QUALSCORE_HIGHRANGE_CUTOFF_DEFAULT = 20;
    public static final String ALIGN_PLUGIN_QUALSCORE_LOWRANGE_COLOR = "qualityScoreLowRangeColor";
    public static final int ALIGN_PLUGIN_QUALSCORE_LOWRANGE_COLOR_DEFAULT = -4144960; // light gray
    public static final String ALIGN_PLUGIN_QUALSCORE_MIDRANGE_COLOR = "qualityScoreMidRangeColor";
    public static final int ALIGN_PLUGIN_QUALSCORE_MIDRANGE_COLOR_DEFAULT = -256;     // yellow
    public static final String ALIGN_PLUGIN_QUALSCORE_HIGHRANGE_COLOR = "qualityScoreHighRangeColor";
    public static final int ALIGN_PLUGIN_QUALSCORE_HIGHRANGE_COLOR_DEFAULT = -1;      // white
    public static final String ALIGN_PLUGIN_SHOW_QUALSCORE = "qualityScoreShow";
    public static final boolean ALIGN_PLUGIN_SHOW_QUALSCORE_DEFAULT = true;

    /** Creates a new instance of TasselPrefs */
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

    public static boolean getIDJoinStrict() {
        // This can be called many times, so to improve performance
        // this will return value without executing system call.
        return TASSEL_IDENTIFIER_JOIN_STRICT_VALUE;
        // return getBooleanPref(TASSEL_TOP, TASSEL_IDENTIFIER_JOIN_STRICT, TASSEL_IDENTIFIER_JOIN_STRICT_DEFAULT);
    }

    public static void putIDJoinStrict(boolean value) {
        TASSEL_IDENTIFIER_JOIN_STRICT_VALUE = value;
        putBooleanPref(TASSEL_TOP, TASSEL_IDENTIFIER_JOIN_STRICT, value);
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

    public static int getFilterAlignPluginMinCount() {
        return getIntPref(FILTER_ALIGN_PLUGIN_TOP, FILTER_ALIGN_PLUGIN_MIN_COUNT, FILTER_ALIGN_PLUGIN_MIN_COUNT_DEFAULT);
    }

    public static void putFilterAlignPluginMinCount(int value) {
        putIntPref(FILTER_ALIGN_PLUGIN_TOP, FILTER_ALIGN_PLUGIN_MIN_COUNT, value);
    }

    //
    // AlignmentPlugin Preferences
    //
    public static int getAlignPluginQualscoreLowrangeCutoff() {
        return getIntPref(ALIGN_PLUGIN_TOP, ALIGN_PLUGIN_QUALSCORE_LOWRANGE_CUTOFF, ALIGN_PLUGIN_QUALSCORE_LOWRANGE_CUTOFF_DEFAULT);
    }

    public static void putAlignPluginQualscoreLowrangeCutoff(int value) {
        putIntPref(ALIGN_PLUGIN_TOP, ALIGN_PLUGIN_QUALSCORE_LOWRANGE_CUTOFF, value);
    }

    public static int getAlignPluginQualscoreHighrangeCutoff() {
        return getIntPref(ALIGN_PLUGIN_TOP, ALIGN_PLUGIN_QUALSCORE_HIGHRANGE_CUTOFF, ALIGN_PLUGIN_QUALSCORE_HIGHRANGE_CUTOFF_DEFAULT);
    }

    public static void putAlignPluginQualscoreHighrangeCutoff(int value) {
        putIntPref(ALIGN_PLUGIN_TOP, ALIGN_PLUGIN_QUALSCORE_HIGHRANGE_CUTOFF, value);
    }

    public static int getAlignPluginQualscoreLowrangeColor() {
        return getIntPref(ALIGN_PLUGIN_TOP, ALIGN_PLUGIN_QUALSCORE_LOWRANGE_COLOR, ALIGN_PLUGIN_QUALSCORE_LOWRANGE_COLOR_DEFAULT);
    }

    public static void putAlignPluginQualscoreLowrangeColor(int value) {
        putIntPref(ALIGN_PLUGIN_TOP, ALIGN_PLUGIN_QUALSCORE_LOWRANGE_COLOR, value);
    }

    public static int getAlignPluginQualscoreMidrangeColor() {
        return getIntPref(ALIGN_PLUGIN_TOP, ALIGN_PLUGIN_QUALSCORE_MIDRANGE_COLOR, ALIGN_PLUGIN_QUALSCORE_MIDRANGE_COLOR_DEFAULT);
    }

    public static void putAlignPluginQualscoreMidrangeColor(int value) {
        putIntPref(ALIGN_PLUGIN_TOP, ALIGN_PLUGIN_QUALSCORE_MIDRANGE_COLOR, value);
    }

    public static int getAlignPluginQualscoreHighrangeColor() {
        return getIntPref(ALIGN_PLUGIN_TOP, ALIGN_PLUGIN_QUALSCORE_HIGHRANGE_COLOR, ALIGN_PLUGIN_QUALSCORE_HIGHRANGE_COLOR_DEFAULT);
    }

    public static void putAlignPluginQualscoreHighrangeColor(int value) {
        putIntPref(ALIGN_PLUGIN_TOP, ALIGN_PLUGIN_QUALSCORE_HIGHRANGE_COLOR, value);
    }

    public static boolean getAlignPluginShowQualscore() {
        return getBooleanPref(ALIGN_PLUGIN_TOP, ALIGN_PLUGIN_SHOW_QUALSCORE, ALIGN_PLUGIN_SHOW_QUALSCORE_DEFAULT);
    }

    public static void putAlignPluginShowQualscore(boolean value) {
        putBooleanPref(ALIGN_PLUGIN_TOP, ALIGN_PLUGIN_SHOW_QUALSCORE, value);
    }

}
