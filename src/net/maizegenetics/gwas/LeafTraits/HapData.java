package net.maizegenetics.gwas.LeafTraits;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

public class HapData {
	FileNames files;
	int chromosome;
	BufferedReader br = null;
	String[] parsedLine;
	int numberOfHaplotypes;
	
	public HapData(FileNames files, int chromosome) {
		this.files = files;
		this.chromosome = chromosome;
		countHaplotypes();
	}
	
	public void countHaplotypes() {
		reset();
		numberOfHaplotypes = 0;
		try {
			String inLine = br.readLine();
			while ((inLine = br.readLine()) != null) {
				numberOfHaplotypes ++;
			}
			br.close();
		} catch (IOException e) {
			System.err.println("Error reading haplotype file. Program ending.");
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	public int getNumberOfHaplotypes() {
		return numberOfHaplotypes;
	}
	
	public boolean next() {
		if (br == null) return false;
		try {
			String inLine = br.readLine();
			if (inLine == null) {
				br.close();
				return false;
			}
			else parsedLine = inLine.split("\t");
		} catch (IOException e) {
			e.printStackTrace();
			return false;
		}
		
		return true;
	}
	
	public int[] getGenotype() {
		int[] geno = new int[26];
		for (int i = 0; i < 26; i++) {
			geno[i] = Integer.parseInt(parsedLine[i + 4]);
		}
		return geno;
	}
	
	public int getPosition() {
		return Integer.parseInt(parsedLine[2]);
	}
	
	public void reset() {
		File file = files.snps[chromosome - 1];

		try {
			if (br != null) br.close();
			br = new BufferedReader(new FileReader(file));
			br.readLine();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
			System.exit(-1);
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		
	}
}
