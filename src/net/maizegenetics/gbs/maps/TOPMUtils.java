/*
 * TOPMUtils
 */
package net.maizegenetics.gbs.maps;

import java.io.File;
import net.maizegenetics.util.Utils;

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

    public static void writeTOPM(TOPMInterface topm, String filename) {

        filename = Utils.addSuffixIfNeeded(filename, ".topm.h5", new String[]{".topm.bin", ".topm.txt", ".topm.h5"});

        String temp = filename.trim().toLowerCase();

        if ((topm instanceof TagsOnPhysicalMap) && (temp.endsWith(".topm.bin"))) {
            ((TagsOnPhysicalMap) topm).writeBinaryFile(new File(filename));
        } else if ((topm instanceof TagsOnPhysicalMap) && (temp.endsWith(".topm.txt"))) {
            ((TagsOnPhysicalMap) topm).writeTextFile(new File(filename));
        } else if ((topm instanceof TagsOnPhysicalMap) && (temp.endsWith(".topm.h5"))) {
            TagsOnPhysMapHDF5.createFile((TagsOnPhysicalMap) topm, filename, 1, topm.getMaxNumVariants());
        } else {
            // TBD
        }

    }
}
