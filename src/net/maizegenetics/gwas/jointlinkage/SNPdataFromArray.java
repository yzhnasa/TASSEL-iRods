package net.maizegenetics.gwas.jointlinkage;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.regex.Pattern;

public class SNPdataFromArray implements SNPdata {
	final Pattern tab = Pattern.compile("\t");
	double[] data;
	String[] pops;
	double[][] geno;
	int numberOfSamples;
	int numberOfMarkers = 1106;
	int startingSNP = 0;
	int endingSNP = 1105;
	int nextMarker = 0;
	String traitName;
	NamMap nammap;
	
	public SNPdataFromArray() {
		int traitcol = 2;
		traitName = "dta";
		int popcol = 4;
		int markercol = 5;
		
		nammap = new NamMap();
		ArrayList<String[]> theData = new ArrayList<String[]>();
		
		
		try {			
			BufferedReader br = new BufferedReader(new FileReader(FileNames.snpFilename));
			
			br.readLine();
			String inline;
			while((inline = br.readLine()) != null) {
				String[] data = tab.split(inline);
				try {
					double val = Double.parseDouble(data[traitcol]);
					int pop = Integer.parseInt(data[popcol]);
//					if (pop != 17) theData.add(data);
					theData.add(data);
				} catch(Exception e) {
					
				}
			}
			
			br.close();
			
			numberOfSamples = theData.size();
			data = new double[numberOfSamples];
			geno = new double[numberOfMarkers][numberOfSamples];
			pops = new String[numberOfSamples]; 
			for (int i = 0; i < numberOfSamples; i++) {
				String[] info = theData.get(i);
				data[i] = Double.parseDouble(info[traitcol]);
				pops[i] = info[popcol];
				for (int j = 0; j < numberOfMarkers; j++) {
					geno[j][i] = Double.parseDouble(info[j + markercol]);
				}
			}
			
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	@Override
	public boolean hasNext() {
		if (nextMarker < numberOfMarkers) return true;
		return false;
	}

	@Override
	public SNP next() {
		return getSnp(nextMarker++);
	}

	@Override
	public synchronized SNP getSnp(int index) {
		int chr = nammap.chromosome.get(index);
		double pos = nammap.position.get(index);
		double[] score = geno[index];
		String name = "m" + nammap.markerNumber.get(index).toString();
		SNP asnp = new SNP(chr, pos, score, name, index);
		return asnp;
	}

	@Override
	public int getNumberOfSNPs() {
		return numberOfMarkers;
	}

	@Override
	public void setStartingSNP(int start) {
		startingSNP = start;
	}

	@Override
	public void setEndingSNP(int end) {
		endingSNP = end;
	}

	@Override
	public int getChromosome(int index) {
		return 0;
	}

	@Override
	public int getPosition(int index) {
		return index + 1;
	}

	@Override
	public void resetSNPs() {
		nextMarker = 0;
	}

	@Override
	public double[] getPhenotype() {
		return data;
	}

	@Override
	public String[] getPopulations() {
		return pops;
	}

	@Override
	public String getPhenotypeName() {
		return traitName;
	}

}
