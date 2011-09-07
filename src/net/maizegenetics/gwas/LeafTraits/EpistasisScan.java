package net.maizegenetics.gwas.LeafTraits;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.regex.Pattern;

public class EpistasisScan {
	int numberOfSamples;
	String[] sampleNames;
	double[][] pheno;
	int[] pop;
	double[][] geno;
	int numberOfMarkers = 1106;
	int phenoNumber;
	String phenoName;
	String outputfilename;
	
	public EpistasisScan(String[] fileNames, int whichPhenotype){
		readDataFile(fileNames[0]);
		phenoNumber = whichPhenotype;
		outputfilename = fileNames[1];
		String[] names = new String[]{"leaf length","leaf width","leaf angle","transformed leaf angle"};
		phenoName = names[whichPhenotype];
		analyzeData();
	}
	
	public static void main(String[] args) {
		try {
			String[] files = new String[]{args[0], args[1]};
			
			EpistasisScan es = new EpistasisScan(files, Integer.parseInt(args[2]) - 1);
		} catch (Exception e) {
			StringBuilder msg = new StringBuilder();
			e.printStackTrace();
			msg.append("EpistasisScan usage:\n");
			msg.append("EpistasisScan dataFileName OutputDirectoryName whichPhenotype\n");
			msg.append("All parameters are required.\n");
			msg.append("The phenotypes are\n");
			msg.append("1 - leaf length\n");
			msg.append("2 - leaf width\n");
			msg.append("3 - leaf angle\n");
			msg.append("4 - transformed angle\n");
			System.out.println(msg);
		}
		
	}
	
	public void readDataFile(String filename) {
		Pattern tab = Pattern.compile("\t");
		try {
			BufferedReader br = new BufferedReader(new FileReader(filename));
			
			//skip the first two lines then read the rest
			br.readLine();
			br.readLine();
			String inline;
			String[] data;
			ArrayList<String> linesOfData = new ArrayList<String>();
			while ((inline = br.readLine()) != null) {
				if (inline.length() > 10) linesOfData.add(inline);
			}
			numberOfSamples = linesOfData.size();
			sampleNames = new String[numberOfSamples];
			pheno = new double[4][numberOfSamples];
			pop = new int[numberOfSamples];
			geno = new double[numberOfMarkers][numberOfSamples];
			
			//then read sample,leaf_length,leaf_width,leaf_angle,transformed_angle,pop, markers 0-1105
			for (int i = 0; i < numberOfSamples; i++) {
				data = tab.split(linesOfData.get(i));
				sampleNames[i] = data[0];
				for (int j = 0; j < 4; j++) pheno[j][i] = Double.parseDouble(data[j + 1]);
				pop[i] = Integer.parseInt(data[5]);
				for (int j = 0; j < numberOfMarkers; j++) geno[j][i] = Double.parseDouble(data[j + 6]); 
			}
			
			
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	public void analyzeData() {
		ExecutorService pool = Executors.newFixedThreadPool(2);
		long start = System.currentTimeMillis();
		for (int p = 1; p <= 26; p++) {
			ArrayList<Integer> popList = new ArrayList<Integer>(200);
			for (int i = 0; i < numberOfSamples; i++) if (pop[i] == p && pheno[phenoNumber][i] > 0) popList.add(i);
			
			int n = popList.size();
			double[] popPheno = new double[n];
			double[][] popGeno = new double[numberOfMarkers][n];
			for (int i = 0; i < n; i++) {
				int whichSample = popList.get(i);
				popPheno[i] = pheno[phenoNumber][whichSample];
				for (int j = 0; j < numberOfMarkers; j++) popGeno[j][i] = geno[j][whichSample];
			}
			
			EpistasisScanForAPopulation theScan = new EpistasisScanForAPopulation(popPheno,popGeno, phenoName, p);
			theScan.setOutputFile(outputfilename);
			pool.execute(theScan);
		}
		System.out.println("All populations have started.");
		pool.shutdown();
	}
	
}
