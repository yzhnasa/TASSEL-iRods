/*
 * Utils.java
 *
 * Created on May 27, 2003, 2:03 AM
 */
package net.maizegenetics.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.InputStreamReader;

import java.net.URL;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Enumeration;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

/**
 *
 * @author  terryc
 */
public final class Utils {

    private static Collection<String> myJavaPackages = null;

    /** Creates a new instance of Utils */
    private Utils() {
    }

    /**
     * Returns the base name of a string delimited
     * with periods (i.e. Java Class).
     *
     * @param str string to parse
     *
     * @return base name
     */
    public static String getBasename(String str) {
        int index = str.lastIndexOf('.');
        index++;
        return str.substring(index);
    }

    /**
     * This returns the filename only.  Preceding
     * directories are removed and  everything after
     * last . is removed.
     *
     * @param str original filename
     * @param suffix suffix
     *
     * @return trimmed filename
     */
    public static String getFilename(String str) {

        int indexForwardSlash = str.lastIndexOf('/');
        int indexBackwardSlash = str.lastIndexOf('\\');

        int index = 0;
        if ((indexForwardSlash == -1) && (indexBackwardSlash == -1)) {
            index = 0;
        } else if (indexForwardSlash > indexBackwardSlash) {
            index = indexForwardSlash + 1;
        } else {
            index = indexBackwardSlash + 1;
        }

        String result = str.substring(index);
        if (result.indexOf('.') > 0) {
            result = result.substring(0, result.indexOf('.'));
        }

        return result;

    }

    /**
     * This returns the filename only.  Preceding
     * directories are removed and suffix.  If suffix not
     * found, then everything after last . is removed.
     *
     * @param str original filename
     * @param suffix suffix
     *
     * @return trimmed filename
     */
    public static String getFilename(String str, String suffix) {

        int indexForwardSlash = str.lastIndexOf('/');
        int indexBackwardSlash = str.lastIndexOf('\\');

        int index = 0;
        if ((indexForwardSlash == -1) && (indexBackwardSlash == -1)) {
            index = 0;
        } else if (indexForwardSlash > indexBackwardSlash) {
            index = indexForwardSlash + 1;
        } else {
            index = indexBackwardSlash + 1;
        }

        String result = str.substring(index);
        if ((suffix != null) && (result.lastIndexOf(suffix) > 0)) {
            result = result.substring(0, result.lastIndexOf(suffix));
        } else if (result.lastIndexOf('.') > 0) {
            result = result.substring(0, result.lastIndexOf('.'));
        }

        return result;

    }

    public static List<String> getFullyQualifiedClassNames(String simpleName) {

        if (myJavaPackages == null) {
            myJavaPackages = getJavaPackages();
        }

        List<String> fqns = new ArrayList<String>();
        for (String aPackage : myJavaPackages) {
            try {
                String fqn = aPackage + "." + simpleName;
                Class.forName(fqn);
                fqns.add(fqn);
            } catch (Exception e) {
                // Do Nothing
            }
        }
        return fqns;

    }

    public static Collection<String> getJavaPackages() {
        String classpath = System.getProperty("java.class.path");
        return getPackagesFromClassPath(classpath);
    }

    public static Set<String> getPackagesFromClassPath(String classpath) {
        Set<String> packages = new HashSet<String>();
        String[] paths = classpath.split(File.pathSeparator);
        for (String path : paths) {
            if (path.trim().length() == 0) {
                continue;
            } else {
                File file = new File(path);
                if (file.exists()) {
                    String childPath = file.getAbsolutePath();
                    if (childPath.endsWith(".jar")) {
                        packages.addAll(readZipFile(childPath));
                    } else {
                        packages.addAll(readDirectory(childPath));
                    }
                }
            }

        }
        return packages;
    }

    public static Set<String> readDirectory(String path) {
        Set<String> packages = new HashSet<String>();
        File file = new File(path);
        int startIndex = path.length() + 1;
        for (File child : file.listFiles()) {
            recursiveRead(child, startIndex, packages);
        }
        return packages;
    }

    public static void recursiveRead(File file, int startIndex, Set<String> packages) {
        if (!file.isFile()) {
            packages.add(file.getAbsolutePath().substring(startIndex).replace(File.separator, "."));
            for (File child : file.listFiles()) {
                recursiveRead(child, startIndex, packages);
            }
        }
    }

    public static Set<String> readZipFile(String path) {
        Set<String> packages = new HashSet<String>();
        try {
            ZipFile zFile = new ZipFile(path);
            Enumeration<? extends ZipEntry> entries = zFile.entries();
            while (entries.hasMoreElements()) {
                ZipEntry entry = entries.nextElement();
                if (!entry.isDirectory()) {
                    String dirName = new File(entry.getName()).getParent();
                    if (dirName != null) {
                        String name = dirName.replace(File.separator, ".");
                        if (name.endsWith(".")) {
                            name = name.substring(0, name.length() - 1);
                        }
                        packages.add(name);
                    }
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        return packages;
    }

    public static String shortenStrLineLen(String str, int preferredLen) {
        return shortenStrLineLen(str, preferredLen, -1);
    }

    /**
     */
    public static String shortenStrLineLen(String str, int preferredLen, int preferredLines) {

        StringBuffer buffer = new StringBuffer();

        int startIndex = 0;
        int endIndex = preferredLen;
        int strLen = str.length();
        int numLines = 0;

        while (startIndex < strLen - 1) {

            int place = str.indexOf(' ', endIndex);
            int newLine = str.indexOf('\n', startIndex);

            if ((newLine != -1) && (newLine < place)) {
                place = newLine;
            }

            String part = null;
            if (place == -1) {
                part = str.substring(startIndex);
                buffer.append(part);
                buffer.append("\n");
                break;
            } else {
                place++;
                part = str.substring(startIndex, place);
                buffer.append(part);
                buffer.append("\n");
                startIndex = place;
                endIndex = place + preferredLen;
            }

            numLines++;

            if ((preferredLines > 0) && (numLines >= preferredLines)) {
                return buffer.toString();
            }

        }

        return buffer.toString();

    }

    /**
     * Adds suffix (i.e. .txt) to end of filename if it's not
     * already there.
     *
     * @param filename filename
     * @param suffix suffix
     *
     * @return filename with suffix
     */
    public static String addSuffixIfNeeded(String filename, String suffix) {

        if (suffix.charAt(0) != '.') {
            suffix = '.' + suffix;
        }

        int periodIndex = filename.lastIndexOf(suffix);
        if (periodIndex == -1) {
            filename = filename + suffix;
        }

        return filename;

    }

    public static BufferedReader getBufferedReader(String inSourceName) {

        try {
            if (inSourceName.startsWith("http")) {
                return new BufferedReader(new InputStreamReader((new URL(inSourceName)).openStream()));
            } else {
                return new BufferedReader(new FileReader(inSourceName));
            }
        } catch (Exception e) {
            System.err.println("File IO in getBufferedReader: " + e);
        }
        return null;
    }

    public static BufferedReader getBufferedReader(String inSourceName, int bufSize) {

        try {
            if (bufSize < 1) {
                return getBufferedReader(inSourceName);
            } else {
                if (inSourceName.startsWith("http")) {
                    return new BufferedReader(new InputStreamReader((new URL(inSourceName)).openStream()), bufSize);
                } else {
                    return new BufferedReader(new FileReader(inSourceName), bufSize);
                }
            }
        } catch (Exception e) {
            System.err.println("File IO in getBufferedReader: " + e);
        }
        return null;
    }

    /**
     * Finds index of Nth occurrence of character in string.
     * 
     * @param str string
     * @param match character to match
     * @param n Nth occurrence
     * 
     * @return index 
     */
    public static int findNthOccurrenceInString(String str, char match, int n) {
        int result = str.indexOf(match);
        while (--n > 0 && result != -1) {
            result = str.indexOf(match, result + 1);
        }
        return result;
    }
    
    /**
     * Returns max heap size in MB.
     * 
     * @return max heap size 
     */
    public static long getMaxHeapSizeMB() {
        return Runtime.getRuntime().maxMemory() / 1048576l;
    }
}
