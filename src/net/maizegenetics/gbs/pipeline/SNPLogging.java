/*
 * SNPLogging
 */
package net.maizegenetics.gbs.pipeline;

import java.io.BufferedWriter;
import java.io.File;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.util.Utils;

/**
 *
 * @author terry
 */
public class SNPLogging {

    private static final String HEADER = "Chr\tPosition\tAlleles\tTagLocusStart\tStrand\tPlugin\tTest\tStatus\tValue\tCuttoff\n";
    private static final String DELIMITER = "\t";
    private final BufferedWriter myWriter;

    public SNPLogging(String filename) {

        boolean exists = false;
        File file = new File(filename);
        if (file.exists()) {
            exists = true;
        }

        myWriter = Utils.getBufferedWriter(filename, true);
        if (!exists) {
            try {
                myWriter.append(HEADER);
            } catch (Exception e) {
                e.printStackTrace();
            }
        }

    }

    public void close() {
        try {
            myWriter.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void writeEntry(Alignment a, int site, String tagLocusStart, String strand, Class plugin, String test, String status, String value, String cuttoff) {
        String alleles = a.getMajorAlleleAsString(site) + "/" + a.getMinorAlleleAsString(site);
        writeEntry(a.getLocusName(site), a.getPositionInLocus(site), alleles, tagLocusStart, strand, plugin, test, status, value, cuttoff);
    }

    public void writeEntry(String chr, int position, String alleles, String tagLocusStart, String strand, Class plugin, String test, String status, String value, String cuttoff) {
        try {
            myWriter.append(chr);
            myWriter.append(DELIMITER);
            myWriter.append(String.valueOf(position));
            myWriter.append(DELIMITER);
            myWriter.append(alleles);
            myWriter.append(DELIMITER);
            myWriter.append(tagLocusStart);
            myWriter.append(DELIMITER);
            myWriter.append(strand);
            myWriter.append(DELIMITER);
            myWriter.append(plugin.getSimpleName());
            myWriter.append(DELIMITER);
            myWriter.append(test);
            myWriter.append(DELIMITER);
            myWriter.append(status);
            myWriter.append(DELIMITER);
            myWriter.append(value);
            myWriter.append(DELIMITER);
            myWriter.append(cuttoff);
            myWriter.append("\n");
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
