package net.maizegenetics.gwas.LeafTraits;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

public class SnpData {
	int chromosome;
	BufferedReader br = null;
	String[] parsedLine;
	int numberOfSnps;
	FileNames files = null;
	
	public final static String[] popnames = new String[]{"B97","CML103","CML228","CML247","CML277","CML322","CML333","CML52R","CML69","HP301","IL14H",
			"KI11","KI3","KY21","M162W","M37W","MO17","MO18W","MS71","NC350","NC358","OH43","OH7B","P39","TX303","TZI8"};
		

	public SnpData(int chromosome) {
		this(chromosome, null);
	}
	
	public SnpData(int chromosome, FileNames files) {
		this.chromosome = chromosome;
		this.files = files;
		findTotalSnpNumber();
		reset();
	}
	
	public void findTotalSnpNumber() {
		reset();
		numberOfSnps = 0;
		try {
			while (br.readLine() != null) numberOfSnps++;
			br.close();
		} catch (IOException e) {
			e.printStackTrace();
			return;
		}
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
		
		if (IsB73Alt())  return next();
		
		return true;
	}
	
	public double[] getGenotype() {
		char ref = parsedLine[1].charAt(0);
		char alt = parsedLine[1].charAt(2);
		int count = 0;
		double[] geno = new double[26];
		for (int i = 12; i <= 38; i++) {
			if (i != 19) {
				char cval = parsedLine[i].charAt(0);
				if (cval == ref) geno[count++] = 0;
				else if (cval == alt) geno[count++] = 1;
				else geno[count++] = Double.NaN;
			}
		}
		return geno;
	}
	
	public int getPosition() {
		return Integer.parseInt(parsedLine[3]);
	}
	
	public void reset() {
		String filename = createInputFileName();

		try {
			if (br != null) br.close();
			br = new BufferedReader(new FileReader(filename));
			br.readLine();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
			System.exit(-1);
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	public String createInputFileName() {
		
		if (files != null) {
			return files.snps[chromosome - 1].getPath();
		}
		
		StringBuilder sb = new StringBuilder();
		
		//non-imputed snps
		sb.append("C:/Projects/NAM/association/0_10.snps.combined.1.13.log2ml.HP1/");
		sb.append(chromosome);
		sb.append(".snps.combined.1.13.log2ml.HP1.txt");
		
		return sb.toString();
	}
	
	public boolean IsB73Alt() {
		char alt = parsedLine[1].charAt(2);
		if (parsedLine[11].charAt(0) == alt) return true;
		return false;
	}
	
	public int getNumberOfSnps() {return numberOfSnps;}
	
	public String getAllele() {return parsedLine[1];}
}
