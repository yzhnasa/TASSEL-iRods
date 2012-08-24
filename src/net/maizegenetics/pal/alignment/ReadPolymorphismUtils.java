package net.maizegenetics.pal.alignment;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.regex.Pattern;

import net.maizegenetics.pal.ids.SimpleIdGroup;

public class ReadPolymorphismUtils {

    private static Pattern WHITESPACE = Pattern.compile("\\s+");

    //prevents instantiation
    private ReadPolymorphismUtils() {
    }

    public static Alignment readPolymorphismFile(String inFile) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(inFile));
        String inline = br.readLine();
        ArrayList<String> markerNames = new ArrayList<String>();
        ArrayList<String[]> dataList = new ArrayList<String[]>();
        int nTaxa = 0;
        int nMarkers = 0;
        float[][] scores = null;
        boolean isNumeric = false;

        if (!inline.startsWith("#") && !inline.startsWith("<")) {
            String[] data = inline.split("[:\\s]+");
            nTaxa = Integer.parseInt(data[0]);
            nMarkers = Integer.parseInt(data[1]);
            data = WHITESPACE.split(br.readLine());
            for (int i = 0; i < nMarkers; i++) {
                markerNames.add(data[i]);
            }
            for (int i = 0; i < nTaxa; i++) {
                inline = br.readLine();
                data = WHITESPACE.split(inline);
                dataList.add(data);
            }
        } else {
            while (inline != null) {
                if (!inline.startsWith("#")) {
                    if (inline.startsWith("<")) {
                        String[] data = inline.split("[<>\\s]+");
                        if (data[1].toLowerCase().startsWith("mark")) {
                            nMarkers = data.length - 2;
                            for (int i = 0; i < nMarkers; i++) {
                                markerNames.add(data[i + 2]);
                            }
                        }
                    } else {
                        dataList.add(WHITESPACE.split(inline));
                    }
                }
                inline = br.readLine();
            }
            nTaxa = dataList.size();
        }

        String[][] myData = new String[nTaxa][nMarkers];
        String[] taxa = new String[nTaxa];
        for (int t = 0; t < nTaxa; t++) {
        	String[] taxonData = dataList.get(t);
        	taxa[t] = taxonData[0];
        	for (int s = 0; s < nMarkers; s++) {
        		myData[t][s] = taxonData[s+1];
        	}
        }
        
        String[] markers = new String[nMarkers];
        markerNames.toArray(markers);
        Locus[] myLoci = new Locus[]{new Locus("Unknown", "0", 0, nMarkers, null, null)};
        return BitAlignment.getInstance(new SimpleIdGroup(taxa), myData, null, null, null, 14, myLoci, new int[]{0}, markers, true, true); 
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
