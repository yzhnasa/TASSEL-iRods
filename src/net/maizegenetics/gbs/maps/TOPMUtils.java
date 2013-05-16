/*
 * TOPMUtils
 */
package net.maizegenetics.gbs.maps;

/**
 *
 * @author terry
 */
public class TOPMUtils {

    private TOPMUtils() {
        // Utility
    }

    /**
     * This reads in a TOPM file. It can be .topm.txt, .topm.bin, or .topm.h5.
     *
     * @param filename filename
     *
     * @return TOPM
     */
    public static TOPMInterface readTOPM(String filename) {

        String temp = filename.trim().toLowerCase();

        if (temp.endsWith(".topm.txt")) {
            return new TagsOnPhysicalMap(filename, false);
        } else if (temp.endsWith(".topm.bin")) {
            return new TagsOnPhysicalMap(filename, true);
        } else if (temp.endsWith(".topm.h5")) {
            return new TagsOnPhysMapHDF5(filename);
        } else {
            throw new IllegalArgumentException("TOPMUtils: readTOPM: Unknown file extension: " + filename);
        }

    }
}
