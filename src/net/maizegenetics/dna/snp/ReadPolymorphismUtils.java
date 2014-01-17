package net.maizegenetics.dna.snp;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.regex.Pattern;

/**
 * Methods for reading non-standard file format.
 *
 * @deprecated
 */
@Deprecated
public class ReadPolymorphismUtils {

    private static Pattern WHITESPACE = Pattern.compile("\\s+");

    //prevents instantiation
    private ReadPolymorphismUtils() {
    }

    public static GeneticMap readGeneticMapFile(String filename) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(filename));
        String inline = br.readLine();
        while (inline != null && inline.startsWith("#")) {
            inline = br.readLine();
        }
        if (!inline.trim().toLowerCase().equals("<map>")) {
            br.close();
            throw new IOException("Expected <Map> as first line of file.");
        }

        GeneticMap theMap = new GeneticMap("Genetic Map: " + filename);

        while (inline != null) {
            if (!inline.startsWith("#")) {
                theMap.addMarker(WHITESPACE.split(inline));
            }
            inline = br.readLine();
        }
        br.close();
        return theMap;
    }
}
