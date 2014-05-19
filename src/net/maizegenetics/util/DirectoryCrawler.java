package net.maizegenetics.util;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

/**
 * @author James Harriman
 * @author Terry Casstevens
 */
public class DirectoryCrawler {

    /**
     * The directory to search if the object is created with no parameters.
     */
    private static final String DEFAULT_DIRECTORY = ".";

    private DirectoryCrawler() {
        // utility class.
    }

    /**
     * Recursively searches the object's directory tree for a specific filename.
     *
     * @param pattern A regular expression specifying what filename to look for.
     * For example, pass the string ".*\.qseq.*" to search for .qseq files.
     * 
     * @return outputFiles A list of matching files found in the directory
     * structure.
     */
    public static File[] listFiles(String pattern, File[] inputArray) {
        List<File> outputList = null;
        for (File file : inputArray) {
            outputList = traverse(file, pattern);  //Otherwise loop over all arguments
        }
        File[] outputArray = outputList.toArray(new File[outputList.size()]);
        return (outputArray);
    }

    /**
     * @param pattern A regular expression specifying a filename to look for.
     * @param inputFile The name of a directory in which to search for files.
     * 
     * @return 
     */
    public static String[] listFileNames(String pattern, String inputFile) {
        File[] outputFiles = listFiles(pattern, inputFile);
        String[] outputNames = new String[outputFiles.length];
        for (int i = 0; i < outputFiles.length; i++) {
            outputNames[i] = outputFiles[i].getAbsolutePath();
        }
        return outputNames;
    }

    public static File[] listFiles(String pattern, File inputFile) {
        return listFiles(pattern, new File[]{inputFile});
    }

    public static File[] listFiles(String pattern) {
        return listFiles(pattern, new File[]{new File(DEFAULT_DIRECTORY)});
    }

    public static File[] listFiles(String pattern, String inputFile) {
        return listFiles(pattern, new File[]{new File(inputFile)});
    }

    public static File[] listFiles(String pattern, String[] inputFiles) {
        File[] inputArray = new File[inputFiles.length];
        for (int i = 0; i < inputFiles.length; i++) {
            inputArray[i] = new File(inputFiles[i]);
        }
        return listFiles(pattern, inputArray);
    }

    private static List<File> traverse(File file, String pattern) {
        List<File> outputList = new ArrayList<File>();
        if (file.isDirectory()) {      // If file is a directory...
            String entries[] = file.list();         // Get a list of all the entries in the directory
            if (entries != null) {   // Ensure that the directory is not empty
                for (String entry : entries) {    // Loop over all the entries
                    List<File> temp = traverse(new File(file, entry), pattern);  // Recursive call to traverse
                    outputList.addAll(temp);
                }
            }
        } else {
            if (file.getName().matches(pattern)) {
                outputList.add(file);
            }
        } //If file is a file, add to list
        return outputList;
    }
}
