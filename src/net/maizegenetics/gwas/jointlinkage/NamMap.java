package net.maizegenetics.gwas.jointlinkage;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

public class NamMap {
	public final String mapfile = "C:/Projects/NAM/data/markers061208.txt";
	public ArrayList<String> markerName = new ArrayList<String>();
	public ArrayList<Integer> chromosome = new ArrayList<Integer>();
	public ArrayList<Double> position = new ArrayList<Double>();
	public ArrayList<Integer> markerNumber = new ArrayList<Integer>();
	public int numberOfMarkers;
	
	public NamMap() {
		BufferedReader br;
		try {
			br = new BufferedReader(new FileReader(mapfile));
			br.readLine();
			String input;
			while ((input = br.readLine()) != null) {
				String[] parsed = input.split("\t");
				markerName.add(parsed[0]);
				chromosome.add(Integer.parseInt(parsed[1]));
				position.add(Double.parseDouble(parsed[2]));
				markerNumber.add(Integer.parseInt(parsed[3]));
			}
			numberOfMarkers = markerName.size();
			br.close();
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}
}
