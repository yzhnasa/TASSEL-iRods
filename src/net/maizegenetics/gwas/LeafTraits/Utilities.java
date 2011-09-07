package net.maizegenetics.gwas.LeafTraits;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class Utilities {
	public static BufferedWriter openOutputFile(String filename) {
		try {
			return new BufferedWriter(new FileWriter(filename));
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		return null;
	}

	public static BufferedWriter openOutputFile(String filename, boolean append) {
		try {
			return new BufferedWriter(new FileWriter(filename, append));
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		return null;
	}

	public static BufferedWriter openOutputFile(File file) {
		try {
			return new BufferedWriter(new FileWriter(file));
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		return null;
	}

	public static BufferedWriter openOutputFile(File file, boolean append) {
		try {
			return new BufferedWriter(new FileWriter(file, append));
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		return null;
	}

	public static void writeToOutput(String out, BufferedWriter bw) {
		try {
			bw.write(out);
			bw.newLine();
		} catch (IOException e) {
			e.printStackTrace();
			try{bw.flush();} catch (Exception e2) {}
			System.exit(-1);
		}
	}

	public static void flushToOutput(String out, BufferedWriter bw) {
		try {
			bw.write(out);
			bw.newLine();
			bw.flush();
		} catch (IOException e) {
			e.printStackTrace();
			try{bw.flush();} catch (Exception e2) {}
			System.exit(-1);
		}
	}
	
	public static void closeOutput(BufferedWriter bw) {
		try {
			bw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public static BufferedReader openReader(String filename) {
		try {
			return new BufferedReader(new FileReader(filename));
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		return null;
	}
	
	public static String readALineFrom(BufferedReader br) {
		String result;
		try {
			result =  br.readLine();
		} catch (IOException e) {
			e.printStackTrace();
			result = null;
		}
		
		if (result == null) closeReader(br);
		return result;
	}
	
	public static void closeReader(BufferedReader br) {
		try {
			br.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public static void mergeFiles() {
		File dir = new File("C:/Projects/NAM/leaf traits/markers");
		String[] ibmfiles = new String[10];
		String[] namfiles = new String[10];
		String[] outfiles = new String[10];
		
		for (int i = 0; i < 10; i++) {
			int chr = i + 1;
			ibmfiles[i] = "imputedIBMMarkersGWAS.chr" + chr + ".082809.txt";
			namfiles[i] = "imputedMarkersGWAS.chr" + chr + ".082809.txt";
			outfiles[i] = "imputedMarkersNAMIBM.chr" + chr + ".112409.txt";
		}
		
		BufferedReader br;
		BufferedWriter bw;
		for (int i = 0; i < 10; i++) {
			try {
				br = new BufferedReader(new FileReader(new File(dir, namfiles[i])));
				bw = new BufferedWriter(new FileWriter(new File(dir, outfiles[i])));
				
				String inputline;
				while((inputline = br.readLine()) != null) {
					bw.write(inputline);
					bw.newLine();
				}
				
				br.close();
				br = new BufferedReader(new FileReader(new File(dir, ibmfiles[i])));
				br.readLine();
				while((inputline = br.readLine()) != null) {
					bw.write(inputline);
					bw.newLine();
				}
				
				br.close();
				bw.close();
				
			} catch (FileNotFoundException e) {
				e.printStackTrace();
				System.exit(-1);
			} catch (IOException e) {
				e.printStackTrace();
				System.exit(-1);
			}
		}
		
	}

}
