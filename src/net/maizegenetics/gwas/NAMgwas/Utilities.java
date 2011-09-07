package net.maizegenetics.gwas.NAMgwas;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeSet;
import java.util.regex.Pattern;

public class Utilities {
	
	public static void main(String[] args) {
//		cnvCounts();
//		Utilities u = new Utilities();
//		u.cnvBPPCounts();
//		head(new File("C:/Projects/NAM/hapmap/v2snps/cnv_ssaha_genes.dat/cnv_ssaha_genes.dat"));
//		countTermsByAnnotation();
//		extractChr9Positions();
		countLines();
	}
	
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
	
	public static void head(File file) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(file));
			for (int i = 0; i < 100; i ++) System.out.println(br.readLine());
			br.close();
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		} 
	}
	
	public static void subset(File inFile, File outFile, int nLines) {
		try{
			BufferedReader br = new BufferedReader(new FileReader(inFile));
			BufferedWriter bw = new BufferedWriter(new FileWriter(outFile));
			
			int count = 0;
			String input;
			while (count < nLines && (input = br.readLine()) != null) {
				count++;
				bw.write(input);
				bw.newLine();
			}
			br.close();
			bw.close();
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	public static void countTypes(File file) {
		TreeSet<String> annotations = new TreeSet<String>();
		Pattern tab = Pattern.compile("\t");
		try {
			BufferedReader br = new BufferedReader(new FileReader(file));
			for (int i = 0; i < 2000000; i ++) {
				String[] data = tab.split(br.readLine());
				annotations.add(data[6]);
			}
			br.close();
			System.out.println("number of annotations = " + annotations.size());
			for (String a : annotations) System.out.println(a);
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
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

	public static void convertHapmapv2(int chromosome) {
//		File inputDir = new File("Z:/SolexaAnal/hapmapV2/hp1/standardQ87Q87Union/N50IMP");
		//chr8.CSHLALLBGI.h90_f50.Q87Q87Union.imp.hmp.txt
		File outputDir = new File("C:/Projects/NAM/hapmap");
		String inName = "chr" + chromosome + ".CSHLALLBGI.h90_f50.Q87Q87Union.imp.hmp.txt";
		String outName = "chr" + chromosome + ".CSHLALLBGI.h90_f50.Q87Q87Union.gwas.hmp.txt";
		String errorFile = "chr" + chromosome + "_conversion_errors.txt";
		Pattern tab = Pattern.compile("\t");
		int[] namcols = new int[]{12,42,47,48,49,51,53,62,63,69,70,71,72,73,74,75,76,77,78,79,80,81,82,84,105,106};
		
		//input col 1 = allele, 2=chromosome, 3=position, 11=B73
		try {
			BufferedReader br = new BufferedReader(new FileReader(new File(outputDir, inName)));
			BufferedWriter err = new BufferedWriter(new FileWriter(new File(outputDir, errorFile)));
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outputDir, outName)));
			bw.write("rs#	alleles	chrom	pos	strand	assembly#	center	protLSID	assayLSID	panelLSID	QCcode	B73	B97	CML103	CML228	CML247	CML277	CML322	CML333	CML52R	CML69	HP301	IL14H	KI11	KI3	KY21	M162W	M37W	MO17	MO18W	MS71	NC350	NC358	OH43	OH7B	P39	TX303	TZI8");
			bw.newLine();
			String input;
			String[] data;
			br.readLine();
			int count = 0;
			while ((input = br.readLine()) != null) {
				if (count % 100000 == 0) {
					System.out.println("line count = " + count);
				}
				count++;
				data = tab.split(input);
//				String[] ref = data[1].split("/");
				String b73Allele = data[11];
				
				// b73 does not equal the reference allele skip
				//check to see if there is at least one non-reference allele in the nam parents
				boolean polymorphic = false;
				int namcnt = 0;
				while (polymorphic == false && namcnt < 26) {
					if (!data[namcols[namcnt]].equals(b73Allele)) polymorphic = true;
					namcnt++;
				}
				
				if (isHeterozygote(b73Allele)) {
					err.write("B73 heterozygous: " + input);
					err.newLine();
				} else if (polymorphic) {
					boolean validLine = true;
					String[] outputData = new String[38];
					for (int i = 0; i < 11; i++) outputData[i] = data[i];
					outputData[11] = "0";
					if (data[3].equals("91960685")) {
						System.out.println("stop here");
					}
					for (int i = 0; i < 26; i++) {
						String inputval = data[namcols[i]];
						if (inputval.equals(b73Allele)) {
							outputData[12 + i] = "0";
						} else if (isHeterozygote(inputval)) {
							outputData[12 + i] = "0.5";
						} else {
							outputData[12 + i] = "1";
						}
					}

					if (validLine) {
						bw.write(outputData[0]);
						for (int i = 1; i < 38; i++) {
							bw.write("\t");
							bw.write(outputData[i]);
						}
						bw.newLine();
					} else {
						err.write("Invalid SNP: " + input);
					}

				} else {
					err.write("Not polymorphic: " + input);
					err.newLine();
				}


			}

			br.close();
			err.close();
			bw.close();

		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}

	public static void processhapmap() {
		int chromosome = 8;
		File inputDir = new File("C:/Projects/NAM/hapmap");
		String inName = "chr" + chromosome + ".CSHLALLBGI.h90_f50.Q87Q87Union.imp.hmp.txt";
		
		try {
			FileInputStream fis = new FileInputStream(new File(inputDir, inName));
			FileChannel infc = fis.getChannel();
			
			int capacity = 10240;
			byte[] input = new byte[capacity];
			
		} catch (FileNotFoundException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	public static boolean isHeterozygote(String snp) {
		//valid hets are R, Y, S, W, K, M
		snp = snp.toUpperCase();
		if (snp.equals("R")) return true;
		if (snp.equals("Y")) return true;
		if (snp.equals("S")) return true;
		if (snp.equals("W")) return true;
		if (snp.equals("K")) return true;
		if (snp.equals("M")) return true;
		if (snp.equals("0")) return true;
		return false;
	}
	
	public static void convertCNVfiles() {
		Pattern tab = Pattern.compile("\t");
		for (int chr = 1; chr <= 10; chr++) {
			String inputFile = "C:/Projects/NAM/hapmap/10kb_cov_matrixes/chr." + chr + ".10kb_log2_matrix.txt";
			String outputFile = "C:/Projects/NAM/hapmap/10kb_cov_matrixes/gwas.chr." + chr + ".10kb_log2_matrix.txt";
			System.out.println("converting " + inputFile);
			try {
				BufferedReader br = new BufferedReader(new FileReader(inputFile));
				BufferedWriter bw = new BufferedWriter(new FileWriter(outputFile));
				
				bw.write("rs#	alleles	chrom	pos	strand	assembly#	center	protLSID	assayLSID	panelLSID	QCcode");
				for (int i = 0; i < 27; i++) bw.write("\tt" + i);
				bw.newLine();
				
				int[] conversion = new int[]{-1,-1,0,1,-1,-1,-1,-1,-1,-1,-1,2,3,33,38,39,40,42,44,52,53,59,60,61,62,63,64,65,66,67,68,69,70,71,72,74,95,96};
				String input;
				
				input = br.readLine();
				String[] data;
				while ((input = br.readLine()) != null) {
					data = tab.split(input);
					bw.write("CNV");
					for (int i = 1; i < 38; i++) {
						if (conversion[i] > -1) {
							bw.write("\t");
							bw.write(data[conversion[i]]);
						} else {
							bw.write("\t");
							bw.write("NA");
						}
					}
					bw.newLine();
				}
				
				br.close();
				bw.close();
			} catch (IOException e) {
				e.printStackTrace();
				System.exit(-1);
			}
			
		}
		
	}
	
	public static void renameFiles() {
		File dir = new File("C:/Projects/NAM/hapmap/results/cnv slb");
//		File dir = new File("c:/users/peter/temp");
		
		String trait = "slb";
		for (int chr = 2; chr <= 10; chr++) {
			//bootstrap.random.model.chr1.v3.txt
			//bootstrap.random.steps.chr1.v3.txt
			//genocheck.txt
			File model = new File(dir, "bootstrap.random.model.chr" + chr + ".v3.txt");
			model.renameTo(new File(dir, trait + ".bootstrap.random.model.chr" + chr + ".v3.txt"));
			
			model = new File(dir, "bootstrap.random.steps.chr" + chr + ".v3.txt");
			model.renameTo(new File(dir, trait + ".bootstrap.random.steps.chr" + chr + ".v3.txt"));
			
		}
	}

	public static void alleleCounts() {
		Pattern tab = Pattern.compile("\t");
		try {

			BufferedWriter bw = new BufferedWriter(new FileWriter("C:/Projects/NAM/hapmap/v2snps/allelecounts.txt"));
			bw.write("chromosome\tallele\tcount");
			for (int chr = 1; chr <=10; chr++) {
				HashMap<String, Integer> alleleMap = new HashMap<String, Integer>();
				String filename = "C:/Projects/NAM/hapmap/v2snps/Chr" + chr + "_nohets_zeroone_cnv.txt";
				BufferedReader br = new BufferedReader(new FileReader(filename));
				String input;
				String[] data;
				br.readLine();
				while ((input = br.readLine()) != null) {
					data = tab.split(input);
					String allele = data[1];
					Integer alleleCount = alleleMap.get(allele);
					if (alleleCount == null) alleleCount = 0;
					alleleMap.put(allele, alleleCount + 1);
				}
				br.close();
				Set<String> keys = alleleMap.keySet();
				TreeSet<String> sortedKeys = new TreeSet<String>(keys);
				bw.newLine();
				for (String key : sortedKeys) {
					StringBuffer sb = new StringBuffer();
					sb.append(chr);
					sb.append("\t").append(key);
					sb.append("\t").append(alleleMap.get(key));
					bw.write(sb.toString());
					bw.newLine();
					System.out.println(sb.toString());
				}
			}
			bw.close();
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	public static void countTermsByAnnotation() {

		class PositionInfo {
			ArrayList<String> annotations = new ArrayList<String>();
			int count = 0;
			
			@Override
			public String toString() {
				StringBuilder sb = new StringBuilder();
				sb.append(count);
				for (String anno : annotations) sb.append("\t").append(anno);
				return sb.toString();
			}
		}
		
		Pattern tab = Pattern.compile("\t");
		ArrayList<HashMap<Integer, PositionInfo>> chromosomeMaps = new ArrayList<HashMap<Integer, PositionInfo>>();
		//Object[] is count (Integer), ArrayList<String> to hold annotations
		String trait = "angle";
		
		System.out.println("Reading step terms...");
		for (int chr = 1; chr <= 10; chr++) {
			String filename = "C:/Projects/NAM/hapmap/results/minor/" + trait + "/" + trait + ".bootstrap.random.steps.chr" + chr + ".v4.txt";
			HashMap<Integer, PositionInfo> thisMap = new HashMap<Integer, PositionInfo>();
			chromosomeMaps.add(thisMap);
			try {
				BufferedReader br = new BufferedReader(new FileReader(filename));
				br.readLine();
				String input;
				while ((input = br.readLine()) != null) {
					String[] data = tab.split(input);
					Integer pos = Integer.decode(data[1]);
					PositionInfo posInfo = thisMap.get(pos);
					if (posInfo == null) {
						posInfo = new PositionInfo();
						thisMap.put(pos, posInfo);
					}
					posInfo.count++;
				}
				br.close();
			} catch(IOException e) {
				e.printStackTrace();
				System.exit(-1);
			}
		}

		//output counts here
		try{
			String filename = "C:/Projects/NAM/hapmap/results/minor/RFR_term_counts.txt";
			BufferedWriter bw = new BufferedWriter(new FileWriter(filename));
			bw.write("chromosome\tpos\tcount");
			bw.newLine();
			for (int c = 0; c < 10; c++) {
				String chr = Integer.toString(c + 1);
				HashMap<Integer, PositionInfo> thisMap = chromosomeMaps.get(c);
				TreeSet<Integer> posSet = new TreeSet<Integer>(thisMap.keySet());
				for (Integer pos : posSet) {
					PositionInfo pi = thisMap.get(pos);
					bw.write(chr + "\t" + pos + "\t" + pi.count);
					bw.newLine();
				}
			}
			bw.close();
		} catch(IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		
		String filename = "C:/Projects/NAM/hapmap/snpinfo/BGI-CSHLALL-UNION.HP1.snp_effects.txt";
		int curr_chr = 0;
		String input = "beginning";
		HashMap<String, Integer> annotationCount = new HashMap<String, Integer>();
		//record annotations for model terms
		try {
			BufferedReader br = new BufferedReader(new FileReader(filename));
			int count = -1;
			while ((input = br.readLine()) != null) {
				String[] data = tab.split(input);
				Integer pos = Integer.decode(data[1]);
				int chr = Integer.parseInt(data[0].substring(3));
				if (chr != curr_chr) {
					TreeSet<String> names = new TreeSet(annotationCount.keySet());
					for (String anno : names) {
						System.out.println("chr" + curr_chr + "\t" + anno + "\t" + annotationCount.get(anno));
					}
					curr_chr = chr;
					annotationCount.clear();
				}
				Integer cnt = annotationCount.get(data[6]);
				if (cnt == null) cnt = annotationCount.put(data[6], new Integer(1));
				else annotationCount.put(data[6], ++cnt);
				HashMap<Integer, PositionInfo> thisMap = chromosomeMaps.get(chr - 1);
				PositionInfo posInfo = thisMap.get(pos);
				if (posInfo != null) {
					posInfo.annotations.add(data[6]);
				}
			}
			TreeSet<String> names = new TreeSet<String>(annotationCount.keySet());
			for (String anno : names) {
				System.out.println("chr" + curr_chr + "\t" + anno + "\t" + annotationCount.get(anno));
			}
			
			
			br.close();
		} catch(IOException e) {
			System.err.println("failed at: " + input);
			e.printStackTrace();
			System.exit(-1);
		}

		System.out.println("Printing annotation summary...");
		
		//get set of annotations
		TreeSet<String> annotationSet = new TreeSet<String>();
		for (HashMap<Integer, PositionInfo> map : chromosomeMaps) {
			for (PositionInfo pi : map.values()) {
				annotationSet.addAll(pi.annotations);
			}
		}
		
		//count annotations by chromosome
		int nAnno = annotationSet.size();
		ArrayList<String> annotationList = new ArrayList<String>(annotationSet);
		int[][] counts = new int[nAnno][10];
		for (int chr = 0; chr < 10; chr++) {
			for (PositionInfo pi : chromosomeMaps.get(chr).values()) {
				for (String anno : pi.annotations) {
					int n = annotationList.indexOf(anno);
					counts[n][chr]++;
				}
			}
		}
		
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter("C:/Projects/NAM/hapmap/snpinfo/hp2_annotation_summary.txt"));
			bw.write("annotation\tchr1\tchr2\tchr3\tchr4\tchr5\tchr6\tchr7\tchr8\tchr9\tchr10");
			bw.newLine();
			for (int i = 0; i < nAnno; i++) {
				bw.write(annotationList.get(i));
				for (int j = 0; j < 10; j++) {
					bw.write("\t");
					bw.write(Integer.toString(counts[i][j]));
				}
				bw.newLine();
			}
			
			bw.close();
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		
		//also output annotations for each term
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter("C:/Projects/NAM/hapmap/snpinfo/hp2_model_annotations.txt"));
			bw.write("chr\tposition\tcount\tannotations");
			bw.newLine();

			for (int c = 0; c < 10; c++) {
				String chr = Integer.toString(c + 1);
				HashMap<Integer, PositionInfo> thisMap = chromosomeMaps.get(c);
				TreeSet<Integer> positions = new TreeSet<Integer>(thisMap.keySet());
				for (Integer pos : positions) {
					PositionInfo pi = thisMap.get(pos);
					bw.write(chr + "\t" + pos + "\t" + pi);
					bw.newLine();
				}
			}
			bw.close();
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		
		System.out.println("Finished");
	}
	
	public static void extractChr9Positions() {
		File inFile = new File("C:/Projects/NAM/hapmap/snpinfo/BGI-CSHLALL-UNION.HP1.snp_effects.txt");
		File outFile = new File("C:/Projects/NAM/hapmap/snpinfo/chr9sites.txt");
		Pattern tab = Pattern.compile("\t");
		
		try {
			BufferedReader br = new BufferedReader(new FileReader(inFile));
			BufferedWriter bw = new BufferedWriter(new FileWriter(outFile));
			
			String input;
			Integer prevSite = -1;
			while ((input = br.readLine()) != null) {
				String[] data = tab.split(input);
				Integer site = Integer.decode(data[1]);
				int chr = Integer.parseInt(data[0].substring(3));
				if (chr == 9 && !site.equals(prevSite)) {
					bw.write(data[0] + "\t" + data[1]);
					bw.newLine();
					prevSite = site;
				}
			}
			br.close();
			bw.close();
		} catch(IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		System.out.println("Finished");
	}
	
	public static void cnvCounts() {
		String trait = "nlb";
		String version = "v7";
		Pattern tab = Pattern.compile("\t");
		File dir = new File("C:/Projects/NAM/hapmap/results/minor_cnv_hp1/" + trait);
		for (int chr = 1; chr <=10; chr++) {
			String filename = trait + ".bootstrap.random.steps.chr" + chr + "." + version + ".txt";
			File stepfile = new File(dir, filename);
			File outfile = new File(dir, "cnv_" + trait + "_summary.txt");
			try {
				BufferedReader br = new BufferedReader(new FileReader(stepfile));
				BufferedWriter bw = new BufferedWriter(new FileWriter(outfile, true));
				br.readLine();
				String input;
				
				int total = 0;
				int nacount = 0;
				int gcount = 0;
				int hp1count = 0;
				while ((input = br.readLine()) != null) {
					String[] data = tab.split(input);
					total++;
					if (data[2].equals("NA")) nacount++;
					else if (data[2].startsWith("GRMZ")) gcount++;
					else if (data[2].startsWith("v1:")) hp1count++;
				}
				bw.write(total + "\t" + nacount + "\t" + gcount + "\t" + hp1count);
				bw.newLine();
				bw.close();
				br.close();
			} catch (IOException e) {
				e.printStackTrace();
				System.exit(-1);
			}
			
		}
	}
	
	public void cnvBPPCounts() {
		Pattern tab = Pattern.compile("\t");
		File dir = new File("C:/Projects/NAM/hapmap/results/minor_cnv/angle");
		for (int chr = 1; chr <=10; chr++) {
			String filename = "angle.bootstrap.random.steps.chr" + chr + ".v6.txt";
			File stepfile = new File(dir, filename);
			File outfile = new File(dir, "cnv_modelterms.txt");
			try {
				HashMap<Snp, Integer> SnpMap = new HashMap<Snp, Integer>(); 
				BufferedReader br = new BufferedReader(new FileReader(stepfile));
				BufferedWriter bw = new BufferedWriter(new FileWriter(outfile, true));
				br.readLine();
				String input;

				while ((input = br.readLine()) != null) {
					String[] data = tab.split(input);
					Snp snp = new Snp(Integer.parseInt(data[1]), data[2]);
					Integer count = SnpMap.get(snp);
					if (count == null) count = new Integer(0);
					count++;
					SnpMap.put(snp, count);
				}

				int total = 0;
				int nacount = 0;
				int gcount = 0;
				
				for (Entry<Snp,Integer> mapEntry : SnpMap.entrySet()) {
					int cnt = mapEntry.getValue().intValue();
					if (cnt > 4) {
						total++;
						String a = mapEntry.getKey().allele;
						if (a.equals("NA")) nacount++;
						else if (a.startsWith("GRMZ")) gcount++;
					}
				}
				
				bw.write(total + "\t" + nacount + "\t" + gcount);
				bw.newLine();
				bw.close();
				br.close();
			} catch (IOException e) {
				e.printStackTrace();
				System.exit(-1);
			}

		}
	}
	
	public static void countLines() {
		String filename = "C:/Projects/NAM/hapmap/v2snps/cnv_ssaha_genes.locs.txt/cnv_ssaha_genes.locs.txt";
		try {
			BufferedReader br = new BufferedReader(new FileReader(filename));
			int count = 0;
			while (br.readLine() != null) count++;
			System.out.println(count + " lines read.");
			br.close();
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	public class Snp {
		int pos;
		String allele;
		
		public Snp(int pos, String allele) {
			this.pos = pos;
			this.allele = allele;
		}

		@Override
		public boolean equals(Object obj) {
			if (!(obj instanceof Snp)) return false;
			Snp other = (Snp) obj;
			if (other.pos != pos) return false;
			if (other.allele.equals(allele)) return true;
			return false;
		}

		@Override
		public int hashCode() {
			return pos;
		}
		
		
	}
}
