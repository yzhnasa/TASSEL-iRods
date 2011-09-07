package net.maizegenetics.gwas.LeafTraits;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

public class AGPMap {
	static final String namMarkers = "C:/Projects/NAM/leaf traits/data/NAM_Map_20090730.txt";
	int[] chrend = new int[10];
	int[] markerPosition;
	double[] markercm;
	String[] marker;
	File mapfile;
	
	public AGPMap() {
		this(new File(namMarkers));
	}
	
	public AGPMap(File mapfile) {
		if (mapfile == null) mapfile = new File(namMarkers);
		this.mapfile = mapfile;
		loadMap();
	}
	
	private void loadMap() {
		try {
			int nRows = 0;
			BufferedReader br = new BufferedReader(new FileReader(mapfile));
			br.readLine();
			while (br.readLine() != null) nRows++;
			br.close();
			
			markerPosition = new int[nRows];
			markercm = new double[nRows];
			marker = new String[nRows];
			br = new BufferedReader(new FileReader(mapfile));
			br.readLine();
			for (int i =  0; i < nRows; i++) {
				String[] parsedLine = br.readLine().split("\t");
				int chr = Integer.parseInt(parsedLine[0]);
				chrend[chr - 1] = i;
				marker[i] = parsedLine[2];
				markercm[i] = Double.parseDouble(parsedLine[3]);
				markerPosition[i] = Integer.parseInt(parsedLine[4]);
			}
			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
			System.exit(-1);
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}

	}
	
	public Object[][] getInterval(int chromosome, int position) {
		int end = chrend[chromosome - 1];
		int start;
		if (chromosome == 1) start = 0;
		else start = chrend[chromosome - 2] + 1;
		
		if (position < markerPosition[start]) 
			return new Object[][]{{null,null}, {marker[start], markerPosition[start]}};
		if (position > markerPosition[end]) 
			return new Object[][]{{marker[end],markerPosition[end]}, {null, null}};
		if (position == markerPosition[start]) 
			return new Object[][]{{marker[start],markerPosition[start]}, {marker[start], markerPosition[start]}};
		if (position == markerPosition[end]) 
			return new Object[][]{{marker[end],markerPosition[end]}, {marker[end], markerPosition[end]}};

		int left = start;
		int right = end;
		
		while (right - left > 1) {
			int mid = left + (right - left)/2;
			if (position == markerPosition[mid]) 
				return new Object[][]{{marker[mid], markerPosition[mid]}, {marker[mid], markerPosition[mid]}};
			if (position < markerPosition[mid]) right = mid;
			else left = mid;
		}
		
		return new Object[][]{{marker[left],markerPosition[left]}, {marker[right], markerPosition[right]}};
	}
	
	public double getCmFromPosition(int chromosome, int position) {
		int end = chrend[chromosome - 1];
		int start;
		if (chromosome == 1) start = 0;
		else start = chrend[chromosome - 2] + 1;
		
		if (position < markerPosition[start]) {
			double g2pRatio = (markercm[start + 1] - markercm[start]) / (markerPosition[start + 1] - markerPosition[start]);
			return markercm[start] - (markerPosition[start] - position) * g2pRatio;
		} 
		if (position > markerPosition[end]) {
			double g2pRatio = (markercm[end] - markercm[end - 1]) / (markerPosition[end] - markerPosition[end - 1]);
			return markercm[end] + (position - markerPosition[end]) * g2pRatio;
		} 
		if (position == markerPosition[start]) return markercm[start];
		if (position == markerPosition[end]) return markercm[end];

		int left = start;
		int right = end;
		
		while (right - left > 1) {
			int mid = left + (right - left)/2;
			if (position == markerPosition[mid]) return markercm[mid];
			if (position < markerPosition[mid]) right = mid;
			else left = mid;
		}
		
		double g2pRatio = (markercm[right] - markercm[left]) / (markerPosition[right] - markerPosition[left]);
		return markercm[left] + (position - markerPosition[left]) * g2pRatio;
	}
	
	public int getPositionFromCm(int chromosome, double cM) {
		int end = chrend[chromosome - 1];
		int start;
		if (chromosome == 1) start = 0;
		else start = chrend[chromosome - 2] + 1;

		if (cM < markercm[start]) {
			double p2gRatio = (markerPosition[start + 1] - markerPosition[start]) / (markercm[start + 1] - markercm[start]);
			return markerPosition[start] - (int) Math.round((markercm[start] - cM) * p2gRatio);
		} 
		if (cM > markercm[end]) {
			double p2gRatio = (markerPosition[end] - markerPosition[end - 1]) / (markercm[end] - markercm[end - 1]);
			return markerPosition[end] + (int) Math.round((cM - markercm[end]) * p2gRatio);
		}

		if (cM == markercm[start]) return markerPosition[start];
		if (cM == markercm[end]) return markerPosition[end];
		
		int left = start;
		int right = end;
		
		while (right - left > 1) {
			int mid = left + (right - left)/2;
			if (cM == markercm[mid]) return markerPosition[mid];
			if (cM < markercm[mid]) right = mid;
			else left = mid;
		}
		
		double p2gRatio = (markerPosition[right] - markerPosition[left]) / (markercm[right] - markercm[left]);
		return markerPosition[left] + (int) Math.round((cM - markercm[left]) * p2gRatio);
	}
	
	public static int getMarkerPosition(Object[] marker){
		if (marker[1] == null) return -1;
		return ((Integer) marker[1]).intValue();
	}
	
	public static String getMarkerName(Object[] marker) {
		return (String) marker[0];
	}
	
	public static int getMarkerNumber(Object[] marker) {
		if (marker[0] == null) return -1;
		return Integer.parseInt(((String) marker[0]).substring(1));
	}
	
	
}
