package net.maizegenetics.gwas.NAMgwas;

import java.io.BufferedReader;
import java.util.ArrayList;

public class NamMap {
	public final String mapfile = "C:/Projects/NAM/data/markers061208.txt";
	public ArrayList<String> markerName = new ArrayList<String>();
	public ArrayList<Integer> chromosome = new ArrayList<Integer>();
	public ArrayList<Double> position = new ArrayList<Double>();
	public ArrayList<Integer> markerNumber = new ArrayList<Integer>();
	public int numberOfMarkers;
	
	public NamMap() {
		BufferedReader br = Utilities.openReader(mapfile);
		Utilities.readALineFrom(br);
		String input;
		while ((input = Utilities.readALineFrom(br)) != null) {
			String[] parsed = input.split("\t");
			markerName.add(parsed[0]);
			chromosome.add(Integer.parseInt(parsed[1]));
			position.add(Double.parseDouble(parsed[2]));
			markerNumber.add(Integer.parseInt(parsed[3]));
		}
		numberOfMarkers = markerName.size();
	}
}
