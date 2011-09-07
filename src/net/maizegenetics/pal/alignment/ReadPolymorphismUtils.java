package net.maizegenetics.pal.alignment;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.regex.Pattern;

import net.maizegenetics.pal.datatype.TextDataType;
import net.maizegenetics.pal.ids.SimpleIdGroup;

public class ReadPolymorphismUtils {
	private static Pattern WHITESPACE = Pattern.compile("\\s+");
	
	//prevents instantiation
	private ReadPolymorphismUtils(){}
	
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
			for (int i = 0; i < nMarkers; i++) markerNames.add(data[i]);
			for (int i = 0; i <nTaxa; i++) {
				inline = br.readLine();
				data = WHITESPACE.split(inline);
				dataList.add(data);
			}
		} else {
			while(inline != null) {
				if (!inline.startsWith("#")) {
					if (inline.startsWith("<")) {
						String[] data = inline.split("[<>\\s]+");
						if (data[1].toLowerCase().startsWith("mark")) {
							nMarkers = data.length - 2;
							for (int i = 0; i < nMarkers; i++) markerNames.add(data[i + 2]);
						}
					} else {
						dataList.add(WHITESPACE.split(inline));
					}
				}
				inline = br.readLine();
			}
			nTaxa = dataList.size();
		}
		
		String[] taxa = new String[nTaxa];
		//int[] pos = new int[nMarkers];
		//for (int i = 0; i < nMarkers; i++) pos[i] = -1;
		Locus locus = new Locus("Unknown", "0", 0, 0, null, null);
		String[] markers = new String[nMarkers];
		markerNames.toArray(markers);
		String[] seq = new String[nTaxa];
		TextDataType dt = new TextDataType();
		for (int i = 0; i < nTaxa; i++) {
			String[] data = dataList.get(i);
			taxa[i] = data[0];
			int n = data.length;
			if (n != nMarkers + 1) {
				StringBuilder sb = new StringBuilder("Error: not enough marker data for taxa ");
				sb.append(data[0]).append(" in line ").append(i);
				throw new IOException(sb.toString());
			}
			StringBuilder sb = new StringBuilder();
			for (int j = 0; j < nMarkers; j++) {
				char aChar;
				if (data[j+1].startsWith("?")) aChar =  TextDataType.UNKNOWN_CHARACTER;
				else aChar = dt.getCharFromTextRepresentation(data[j+1]);
				sb.append(aChar);
			}
			seq[i] = sb.toString();
		}
		
		return new SimpleAlignment(new SimpleIdGroup(taxa), seq, dt, null, null, null, locus, scores, markers, false);
	}
	
	public static GeneticMap readGeneticMapFile(String filename) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(filename));
		String inline = br.readLine();
		while (inline != null && inline.startsWith("#")) inline = br.readLine();
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