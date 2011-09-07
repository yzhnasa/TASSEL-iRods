package net.maizegenetics.gwas.NAMgwas;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;
import java.util.regex.Pattern;


import cern.jet.stat.Probability;


public class NamGBSData {
	int chr = 0;
	final ArrayList<String> taxanames = new ArrayList<String>();
	final ArrayList<String> snplabels = new ArrayList<String>();
	final String[] headerForLabels = new String[11];
	int[] snpPositions;
	byte[][] snps;
	byte[][] imputedSnps;
	byte[] b73Allele;
	byte[] nonb73Allele;
	int numberOfTaxa, numberOfSnps, numberOfPopulations;
	
	BufferedWriter hapout;
	
	final HashMap<String, Byte> snpToByteMap = new HashMap<String, Byte>();;
	final HashMap<String, ArrayList<Integer>> popmap = new HashMap<String, ArrayList<Integer>>(); //which taxa are in each population
	final HashMap<String, boolean[]> segregatingSites = new HashMap<String, boolean[]>();
	
	static String[] byteToSnp = new String[]{"A","C","G","T","R","Y","S","W","K","M","+","0","H","B","N","-","X"};

	static final byte A = (byte) 0;
	static final byte C = (byte) 1;
	static final byte G =(byte) 2;
	static final byte T =(byte) 3;
	static final byte R =(byte) 4;		//A/G
	static final byte Y =(byte) 5; 	//C/T
	static final byte S =(byte) 6; 	//G/C
	static final byte W =(byte) 7;		//A/T
	static final byte K =(byte) 8;		//G/T
	static final byte M =(byte) 9;		//A/C	
	static final byte plus =(byte) 10;	
	static final byte zero =(byte) 11;	//+/-
	static final byte H =(byte) 12;	
	static final byte B =(byte) 13;	
	static final byte N =(byte) 14;	
	static final byte minus =(byte) 15;
	static final byte X = (byte) 16;
	
	final byte[][] hetcodes = new byte[][]{{A,M,R,W},{M,C,S,Y},{R,S,G,K},{W,Y,K,T}}; 
	


	private void loadSnpToByte() {
		snpToByteMap.put("A", A);
		snpToByteMap.put("C", C);
		snpToByteMap.put("G", G);
		snpToByteMap.put("T", T);
		snpToByteMap.put("R", R);
		snpToByteMap.put("Y", Y);
		snpToByteMap.put("S", S);
		snpToByteMap.put("W", W);
		snpToByteMap.put("K", K);
		snpToByteMap.put("M", M);
		snpToByteMap.put("+", plus);
		snpToByteMap.put("0", zero);
		snpToByteMap.put("H", H);
		snpToByteMap.put("B", B);
		snpToByteMap.put("N", N);
		snpToByteMap.put("-", minus);
	}
	
	public NamGBSData(String filename) {
		loadSnpToByte();
		ArrayList<Integer> columns = new ArrayList<Integer>();
		String[] info;
		String input;
		Pattern tab = Pattern.compile("\t");
		Pattern slash = Pattern.compile("/");
		int count;
		System.out.println("Reading the snp file...");
		try {
			BufferedReader br = new BufferedReader(new FileReader(filename));
			
			//count the number of lines in the file and get the chromosome number
			count = 1;
			br.readLine(); //the header
			info = tab.split(br.readLine());
			chr = Integer.parseInt(info[2]);
			while (br.readLine() != null) count++;
			br.close();
			numberOfSnps = count;
			
			br = new BufferedReader(new FileReader(filename));
			info = tab.split(br.readLine());
			for (int i = 0; i < 11; i++) headerForLabels[i] = info[i];
			count = 0;
			for (String str : info) {
				if (str.startsWith("Z0")) {
					taxanames.add(str);
					columns.add(count);
				}
				count++;
			}
			numberOfTaxa = taxanames.size();
			getpops();
			
			snps = new byte[numberOfTaxa][numberOfSnps];
			b73Allele = new byte[numberOfSnps];
			nonb73Allele = new byte[numberOfSnps];
			snpPositions = new int[numberOfSnps];
			int incount = 0;
			while ((input = br.readLine()) != null) {
				info = tab.split(input);
				String[] alleles = slash.split(info[1]);
				b73Allele[incount] = snpToByteMap.get(alleles[0]);
				try {nonb73Allele[incount] = snpToByteMap.get(alleles[1]);} catch(Exception e) {
					System.out.println("Error in second allele: " + input);
				}
				
				snpPositions[incount] = Integer.parseInt(info[3]);
				count = 0;
				StringBuilder sb = new StringBuilder(info[0]);
				for (int i = 1; i < 11; i++) sb.append("\t").append(info[i]);
				snplabels.add(sb.toString());
				for (Integer col:columns) {
					snps[count++][incount] = snpToByteMap.get(info[col]);
				}
				incount++;
			}
			
			br.close();
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		
		System.out.println("Imputing snps...");
		imputeABH();
		writeTaxa("C:/Projects/NAM/namgbs/datasets/cov20/cov20");
	}

	private void getpops() {
		int count = 0;
		for (String taxon:taxanames) {
			String pop = taxon.substring(0,4);
			ArrayList<Integer> poplist = popmap.get(pop);
			if (poplist == null) {
				poplist = new ArrayList<Integer>();
				popmap.put(pop, poplist);
				segregatingSites.put(pop, new boolean[numberOfSnps]);
			}
			poplist.add(count++);
		}
		numberOfPopulations = popmap.size();
		
	}
	
	public void imputeABH() {
		imputedSnps = new byte[numberOfTaxa][numberOfSnps];
		
		//count alleles by population. 
		//consider all populations with total allele count >= 15.
		int minCallsPerPopulation = 15;
		
		//a pop is non-segregating if there is only one allele or if the p(sample drawn from a binary dist. with p = 0.5)<.01
		//if all non segregating pops have the same major allele then call that B73
		//use the value 16(X) to indicated that a snp has not been imputed.
		
		//diagnostic variables
		HashMap<String, Integer> popCountMap = new HashMap<String, Integer>();
		for (String pop : popmap.keySet()) {
			popCountMap.put(pop, 0);
		}
		
		//for each snp
		System.out.println("Making initial ABH calls.");
		for (int s = 0; s < numberOfSnps; s++) {
			//count alleles in each population
			byte[][] popAlleles = new byte[numberOfPopulations][2];
			int popcount = 0;
			ArrayList<Integer> segList = new ArrayList<Integer>();
			for (Map.Entry<String, ArrayList<Integer>> ent:popmap.entrySet()) {
				boolean[] segsites = segregatingSites.get(ent.getKey());
				int[] allelecounts = new int[16];
				for (Integer taxon : ent.getValue()) {
					allelecounts[snps[taxon][s]]++;
					imputedSnps[taxon][s] = N; //initialize imputed snps to N (missing)
				}
				int high = allelecounts[0];
				int next = 0;
				byte major = A;
				byte minor = X;
				
				for (int i = 1; i < 4; i++) {
					if (allelecounts[i] > high) {
						next = high;
						high = allelecounts[i];
						minor = major;
						major = (byte) i;
					} else if (allelecounts[i] > next) {
						next = allelecounts[i];
						minor = (byte) i;
					}
				}
				
				//is this segregating?
				if (high + next < minCallsPerPopulation) { //not enough samples
					popAlleles[popcount][0] = X;
					popAlleles[popcount][1] = X;
					segsites[s] = false;
				} else {
					double pval = Probability.binomial(next, high + next, 0.5);
					if (pval < .01) {  //not segregating
						popAlleles[popcount][0] = major;
						popAlleles[popcount][1] = X;
						segsites[s] = false;
					} else { //segregating
						popAlleles[popcount][0] = major;
						popAlleles[popcount][1] = minor;
						segsites[s] = true;
					}
				}
				
				popcount++;
			}
			
			//can B73 be called?
			byte B73 = X;
			byte notB73 = X;
			boolean consistent = true;
			
			//is this site segregating or not callable in all populations?
			//if so, add to allPopSegList and do not call for now by setting consistent = false
			int nonsegpops = 0;
			int notcallable = 0;
			for (int p = 0; p < numberOfPopulations; p++) {
				if (popAlleles[p][0] == X && popAlleles[p][1] == X) notcallable++;
				if (popAlleles[p][0] == X || popAlleles[p][1] == X) nonsegpops++;
			}
			
			if (nonsegpops == 0 && notcallable < numberOfPopulations) {
				consistent = false;
			}
			
			//is the B73 allele the same in all non-seg populations. That is all non-seg populations should have the same allele.
			for (int p = 0; p < numberOfPopulations; p++) if (consistent){
				if (popAlleles[p][1] == X && popAlleles[p][0] != X) {
					if (B73 == X) B73 = popAlleles[p][0];
					else if (B73 != popAlleles[p][0]) consistent = false;
				}
			}
			
			if (consistent && B73 != X) { //find nonB73
				for (byte[] alleles : popAlleles) {
					if (alleles[0] == X) {
						//do nothing
					} else if (alleles[0] != B73) {
						if (notB73 == X) notB73 = alleles[0];
						else if (notB73 != alleles[0]) consistent = false;
					} else if (alleles[1] != X) {
						if (notB73 == X) notB73 = alleles[1];
						else if (notB73 != alleles[1]) consistent = false;
					}
				}
			}
			
			//make ABH calls only in segregating populations
			if (consistent && notB73 != X && B73 != X) {
				b73Allele[s] = B73;
				nonb73Allele[s] = notB73;
				byte het = hetcodes[B73][notB73];
				for (String pop : popmap.keySet()) {
					if (segregatingSites.get(pop)[s]) {
						popCountMap.put(pop, popCountMap.get(pop) + 1);
						ArrayList<Integer> taxaList = popmap.get(pop);
						for (Integer t:taxaList) {
							if (snps[t][s] == B73) imputedSnps[t][s] = A;
							else if (snps[t][s] == notB73) imputedSnps[t][s] = B;
							else if (snps[t][s] == het) imputedSnps[t][s] = H;
						}
					}
				}
			}
			
		} //end snp loop: for (int s = 0; s < numberOfSnps; s++)
		
		//output population counts
//		System.out.println("Initial set segregating sites\nAll,Selected");
//		for (String pop:segregatingSites.keySet()) {
//			int nseg = 0;
//			for (boolean b:segregatingSites.get(pop)) if (b) nseg++;
//			System.out.println(pop + "," + nseg + "," + popCountMap.get(pop));
//		}
		
		//check snps against surrounding haplotypes
		System.out.println("Checking snps against enclosing haplotypes.");
		
		//which sites are segregating in at least one population
		boolean[] checkSite = new boolean[numberOfSnps];
		for (int s = 0; s < numberOfSnps; s++) checkSite[s] = false;
		
		for (boolean[] segsites:segregatingSites.values()) {
			for (int s = 0; s < numberOfSnps; s++) {
				checkSite[s] = checkSite[s] || segsites[s];
			}
		}
		
		//set non-segregating sites to N in imputedSnps
		for (int s = 0; s < numberOfSnps; s++) {
			if (!checkSite[s]) {
				for (int t = 0; t < numberOfTaxa; t++) imputedSnps[t][s] = N; 
			}
		}
		
		//generate some stats at this point
//		System.out.println("pop,A,B,N,total");
//		for (String pop : popmap.keySet()) {
//			ArrayList<Integer> taxaList = popmap.get(pop);
//			int Acount = 0;
//			int Bcount = 0;
//			int Ncount = 0;
//			int total = 0;
//			for (Integer t : taxaList) {
//				for (int s = 0; s < numberOfSnps; s++) {
//					if (imputedSnps[t][s] == A) Acount++;
//					else if (imputedSnps[t][s] == B) Bcount++;
//					else if (imputedSnps[t][s] == N) Ncount++;
//					total++;
//				}
//			}
//			System.out.println(pop + "," + Acount + "," + Bcount + "," + Ncount + "," + total);
//		}
//		System.exit(0);
		
		//this checks all SNPs against the surrounding haplotype. It checks SNPs in the input data against the haplotypes in imputed SNPs
		//that have been called based on segregation ratios in individual populations. ImputedSNPs is a more stringent dataset at this point.
		//All SNPs are tested against it. Those that are consistent with the surrounding haplotype are added to ImputedSNPs.
		imputeSnpsBasedOnHaplotypes(0.85);
		
		//generate some stats at this point
		System.out.println("pop,A,B,N,total");
		for (String pop : popmap.keySet()) {
			ArrayList<Integer> taxaList = popmap.get(pop);
			int Acount = 0;
			int Bcount = 0;
			int Ncount = 0;
			int total = 0;
			for (Integer t : taxaList) {
				for (int s = 0; s < numberOfSnps; s++) {
					if (imputedSnps[t][s] == A) Acount++;
					else if (imputedSnps[t][s] == B) Bcount++;
					else if (imputedSnps[t][s] == N) Ncount++;
					total++;
				}
			}
			System.out.println(pop + "," + Acount + "," + Bcount + "," + Ncount + "," + total);
		}

//		System.exit(0);
		
		//for individual taxa, look for probable errors
		//change to missing if an A has hapsize B's on both sides change to B, likewise A's & H's
		System.out.println("Identifying errors.");
		int hapsize = 6;
		for (String popname : popmap.keySet()) {
			ArrayList<Integer> taxaList = popmap.get(popname);
			for (Integer t : taxaList) {
				byte[] taxonSnps = imputedSnps[t];

				//count non-missing snps
				int nonmissCount = 0;
				for (byte snp:taxonSnps) {
					if (snp != N) nonmissCount++;
				}
				if (nonmissCount < 10) {
					System.out.println(nonmissCount + " non-missing data points for " + taxanames.get(t) );
					for (int s = 0; s < numberOfSnps; s++) imputedSnps[t][s] = N;

				} else {
					//create an index of non-missing counts
					int[] nonmissIndex = new int[nonmissCount];
					int count = 0;
					for (int s = 0; s < numberOfSnps; s++) {
						if (taxonSnps[s] != N) nonmissIndex[count++] = s;
					}

					//count number of adjacent matching snps from left to right
					int[] leftcount = new int[nonmissCount];
					leftcount[0] = 1;
					for (int s = 1; s < nonmissCount; s++) {
						if (taxonSnps[nonmissIndex[s]] == taxonSnps[nonmissIndex[s-1]]) leftcount[s] = leftcount[s-1] + 1;
						else leftcount[s] = 1;
					}

					//count number of adjacent matching snps from right to left
					int[] rightcount = new int[nonmissCount];
					rightcount[nonmissCount - 1] = 1;
					for (int s = nonmissCount - 2; s >= 0; s--) {
						if (taxonSnps[nonmissIndex[s]] == taxonSnps[nonmissIndex[s+1]]) rightcount[s] = rightcount[s+1] + 1;
						else rightcount[s] = 1;
					}

					
					
					
					//if left and right counts = one and left - 1 >= hapsize and right + 1 >= hapsize then set to N
					//also at position s, if leftcount[s-1] = s and rightcount[s+1] >= hapsize, then set to N
					//similarly if leftcount[s-1] >= hapsize and rightcount[s+1] == nmcount - s, then set to N
					for (int s = 0; s < nonmissCount; s++) {
						if (leftcount[s] == 1 && rightcount[s] == 1) {
							int sIndex = nonmissIndex[s];
							if (s == 0) {
								if (rightcount[s+1] >= hapsize) taxonSnps[sIndex] = N;
							} else if (s <= hapsize) {
								if (rightcount[s+1] >= hapsize && leftcount[s - 1] == s) taxonSnps[sIndex] = N;
							} else if (s == nonmissCount - 1) {
								if (leftcount[s-1] >= hapsize) taxonSnps[sIndex] = N;
							} else if (s >= nonmissCount - hapsize - 1) {
								if (leftcount[s-1] >= hapsize && rightcount[s+1] == nonmissCount - s) taxonSnps[sIndex] = N;
							} else {
								if (leftcount[s-1] >= hapsize && rightcount[s+1] >= hapsize) taxonSnps[sIndex] = N;
							}
						}
						
						//if snp pairs are within 128 bp treat as a single snp
						if (s < nonmissCount - 1 && leftcount[s] == 1 && rightcount[s] == 2 && leftcount[s+1] == 2 && rightcount[s+1] == 1) {
							int sIndexLeft = nonmissIndex[s];
							int sIndexRight = nonmissIndex[s+1];
							if (snpPositions[sIndexRight] - snpPositions[sIndexLeft] < 128) { //the snps are within 128 bp of each other
								if (s == 0) {
									if (rightcount[s+2] >= hapsize) {
										taxonSnps[sIndexRight] = N;
										taxonSnps[sIndexLeft] = N;
									}
								} else if (s <= hapsize) {
									if (rightcount[s+2] >= hapsize && leftcount[s - 1] == s) {
										taxonSnps[sIndexRight] = N;
										taxonSnps[sIndexLeft] = N;
									}
								} else if (s == nonmissCount - 2) {
									if (leftcount[s-1] >= hapsize) {
										taxonSnps[sIndexRight] = N;
										taxonSnps[sIndexLeft] = N;
									}
								} else if (s >= nonmissCount - hapsize - 2) {
									if (leftcount[s-1] >= hapsize && rightcount[s+2] == nonmissCount - s - 1) {
										taxonSnps[sIndexRight] = N;
										taxonSnps[sIndexLeft] = N;
									}
								} else {
									if (leftcount[s-1] >= hapsize && rightcount[s+2] >= hapsize) {
										taxonSnps[sIndexRight] = N;
										taxonSnps[sIndexLeft] = N;
									}
								}
							}
							s++;
						}
					}
				}
			}
		}
		
	}
	
	
	/**
	 * This function finds the proportion of snps five snps on either side of snpIndex that match the A/B call of snpIndex. 
	 * Fewer than five may be used at the ends of the chromosomes.
	 * @param snpIndex	the index of the snp being tested
	 * @param popnames	names of populations
	 * @return for each population, and int array containing the number of A matches, the number of B matches and the total count
	 */
	public int[][] testEnclosingHaplotype(int snpIndex, List<String> popnames) {
		//test this snp against enclosing haplotypes for each population
		//the function returns the proportion of A/B calls for any snp that match the 5 A/B calls adjacent on either side
		
//		HashMap<String, Double> pmatch = new HashMap<String, Double>();
//		double[] pmatch = new double[popnames.size()];
		int[][] nmatch = new int[popnames.size()][3];
		int sitecount;
		int hapsize = 5;
		byte b73 = b73Allele[snpIndex];
		byte nonb73 = nonb73Allele[snpIndex];
		int popcount = 0;
		for (String popname:popnames) { //iterate through populations
			ArrayList<Integer> taxaList = popmap.get(popname);
			int matchCount = 0;
			int AmatchCount = 0;
			int BmatchCount = 0;
			int totalCount = 0; 
			for (Integer t:taxaList) { //iterate through taxa in one population
				byte snp = snps[t][snpIndex];
				byte snptype;
				if (snp == b73) snptype = A;
				else if (snp == nonb73) snptype = B;
				else snptype = N;
				if (snptype != N) {
					//count matches to the left
					sitecount = 0;
					int testsite = snpIndex - 1;
					while (sitecount < hapsize && testsite >= 0) {
						byte testsnp = imputedSnps[t][testsite--];
//						if (testsnp == A || testsnp == B) {
//							totalCount++;
//							sitecount++;
//							if (testsnp == snptype) matchCount++;
//						}
						if (testsnp == A) {
							totalCount++;
							sitecount++;
							if (testsnp == snptype) {
								matchCount++;
								AmatchCount++;
							}
						}
						else if (testsnp == B) {
							totalCount++;
							sitecount++;
							if (testsnp == snptype) {
								matchCount++;
								BmatchCount++;
							}
						}
					}
					
					//count matches to the right
					sitecount = 0;
					testsite = snpIndex + 1;
					while (sitecount < hapsize && testsite < numberOfSnps) {
						byte testsnp = imputedSnps[t][testsite++];
//						if (testsnp == A || testsnp == B) {
//							totalCount++;
//							sitecount++;
//							if (testsnp == snptype) matchCount++;
//						}
						if (testsnp == A) {
							totalCount++;
							sitecount++;
							if (testsnp == snptype) {
								matchCount++;
								AmatchCount++;
							}
						}
						else if (testsnp == B) {
							totalCount++;
							sitecount++;
							if (testsnp == snptype) {
								matchCount++;
								BmatchCount++;
							}
						}
					}
				}
			}
			
//			if (totalCount == 0) pmatch[popcount] = Double.NaN;
//			else pmatch[popcount] = ((double) matchCount) / ((double) totalCount);
			
			nmatch[popcount][0] = AmatchCount;
			nmatch[popcount][1] = BmatchCount;
			nmatch[popcount][2] = totalCount;
			
			popcount++;
		}
		
		return nmatch;
	}
	
	/*
	public void imputeSnpBasedOnHaplotype(int snpIndex, double limit, List<String> popnames) {
		double[] pmatch = testEnclosingHaplotype(snpIndex, popnames);
		int overcount = 0;
		int undercount = 0;
		for (double p:pmatch) {
			if (p > limit) overcount++;
			else if (p < (1 - limit)) undercount++;
		}
		
		if (overcount > 0 && undercount == 0) { //impute snp values in the pops with overcount>0 and set to missing in the rest
			imputeSnpValuesInSegregatingPopulations(snpIndex, pmatch, limit);
			
		} else if (overcount == 0 && undercount > 0){  //switch b73/nonb73 and impute snp values
			byte temp = b73Allele[snpIndex];
			b73Allele[snpIndex] = nonb73Allele[snpIndex];
			nonb73Allele[snpIndex] = temp;
			imputeSnpValuesInSegregatingPopulations(snpIndex, pmatch, limit);
		} else { //set all snp values to missing
			for (int t = 0; t < numberOfTaxa; t++) imputedSnps[t][snpIndex] = N;
		}
	}
	*/
	
	/**
	 * Infers whether called SNPs are either A type, B type, or cannot be called reliably. 
	 * The function sets values to A/B/N in the array imputedSnps. Only existing calls are imputed as A/B. Missing data is not imputed.
	 * @param limit the proportion of SNPs that must match for the SNP to be considered okay
	 */
	public void imputeSnpsBasedOnHaplotypes (double limit) {
		byte[][] impute2 = new byte[numberOfTaxa][numberOfSnps];
		List<String> popnames = new ArrayList<String>(popmap.keySet());
		Collections.sort(popnames);
		for (int s = 0; s < numberOfSnps; s++) {
			//pmatch contains the proportion of neigboring SNPs with A/B matching the SNP at position s in each population
			//overcount contains the number of populations in which this SNP matches the surrounding haplotype
			//undercount contains the number of populations in which this SNP would match the surrounding haplotype if the SNP call were changed,
			//indicating that the A/B designation of the SNP should be changed
//			double[] pmatch = testEnclosingHaplotype(s, popnames);
			int[][] nMatch = testEnclosingHaplotype(s, popnames);
			int n = nMatch.length;
			double[] pmatch = new double[n];
			for (int i = 0; i < n; i++) {
				if (nMatch[i][2] == 0) pmatch[i] = Double.NaN;
				else pmatch[i] = ((double)(nMatch[i][0] + nMatch[i][1])) / ((double) nMatch[i][2]);
			}
			
			int overcount = 0;
			int undercount = 0;
			for (double p:pmatch) {
				if ( !Double.isNaN(p) && p > limit ) overcount++;
				else if ( !Double.isNaN(p) && p < (1 - limit) ) undercount++;
			}
			
			//if both overcount > 0 and undercount > 0, there is a problem and the snp should not be used
			if (overcount > 0 && undercount == 0 || overcount == 0 && undercount > 0) { //impute snp values

				if (overcount == 0 && undercount > 0){  //switch b73/nonb73 
					byte temp = b73Allele[s];
					b73Allele[s] = nonb73Allele[s];
					nonb73Allele[s] = temp;
				}
				
				//impute snp values for populations in which the SNP matches the haplotype
				//and there are at least two B calls and two A calls to ensure that the population is segregating
				int popcount = 0;
				byte b73 = b73Allele[s];
				byte nonb73 = nonb73Allele[s];
				for (String popname:popnames) {
					double pval = pmatch[popcount];
					int Acount = nMatch[popcount][0];
					int Bcount = nMatch[popcount][1];
					int mincount = Math.max(2, (int)(0.05 * (Acount + Bcount) ));
					if ((pval > limit || pval < 1 - limit) && Acount >= mincount && Bcount >= mincount){
						for (Integer t:popmap.get(popname)) {
							if (snps[t][s] == b73) {
								impute2[t][s] = A;
							} else if (snps[t][s] == nonb73) {
								impute2[t][s] = B;
							} else {
								impute2[t][s] = N;
							}
						}
					} else {
						for (Integer t:popmap.get(popname)) {
							impute2[t][s] = N;
						}
					}
					popcount++;
				}

			} else { //set all snp values to missing
				for (int t = 0; t < numberOfTaxa; t++) impute2[t][s] = N;
			}
			
			
		}
		
		imputedSnps = impute2;
	}
	
	public void imputeSnpValuesInSegregatingPopulations(int snpIndex, double[] pmatch, List<String> popnames, double limit) {
		byte b73 = b73Allele[snpIndex];
		byte nonb73 = nonb73Allele[snpIndex];
		byte het = getHetCode(b73, nonb73);
		int popcount = 0;
		for (String popname:popnames) {
			double pval = pmatch[popcount];
			if (pval > limit || pval < 1 - limit){
				for (Integer t:popmap.get(popname)) {
					if (snps[t][snpIndex] == b73) {
						imputedSnps[t][snpIndex] = A;
					} else if (snps[t][snpIndex] == nonb73) {
						imputedSnps[t][snpIndex] = B;
					} else {
						imputedSnps[t][snpIndex] = N;
					}
				}
			} else {
				for (Integer t:popmap.get(popname)) {
					imputedSnps[t][snpIndex] = N;
				}
			}
			popcount++;
		}
	}
	
	public static void callHets(byte[] snps) {
		int sizecount = 1;
		int segcount = 0;
		int hetstart = 0;
		int segstart = 0;
		byte prevbyte = -1;
		int hetlimit = 11;  //the largest number of contiguous calls of the same type allowed to call a het segment
		int minsegments = 3;  //the minimum number of contiguous het segments required to call a het
		int numberOfSnps = snps.length;
		
		for (int s = 0; s < numberOfSnps; s++) {
			byte thisbyte = snps[s];
			if (thisbyte != N) {
				if ( thisbyte == prevbyte) {
					sizecount++;
				} else {
					if (sizecount > hetlimit) {  //reset everything
						if (segcount >= minsegments) {  //set bytes between start and current segment (exclusive) to H
							for (int h = hetstart; h < segstart; h++) {
								if (snps[h] != N) snps[h] = H;
							}
						}
						sizecount = 1;
						segcount = 0;
						hetstart = s;
						segstart = s;
						prevbyte = thisbyte;
					} else {  //count as a het segment, reset segstart, and keep going
						sizecount = 1;
						segcount++;
						segstart = s;
						prevbyte = thisbyte;
					}
				}
			}
		}
		
		//test for stuff at final snp
		if (sizecount > hetlimit) {
			if (segcount >= minsegments) {  //set bytes between start and current segment (exclusive) to H
				for (int h = hetstart; h < segstart; h++) {
					if (snps[h] != N) snps[h] = H;
				}
			}
		} else {
			if (segcount >= minsegments - 1) {  //set bytes between start and this site (inclusive) to H
				for (int h = hetstart; h < numberOfSnps; h++) {
					if (snps[h] != N) snps[h] = H;
				}
			}
		}
		
	}
	
	public static void callHetsFromABHFile(String inputfile, String outputfile) {
		Pattern tab = Pattern.compile("\t");
		HashMap<String, Byte> snpToByte = new HashMap<String, Byte>();
		snpToByte.put("A", A);
		snpToByte.put("H", H);
		snpToByte.put("B", B);
		snpToByte.put("N", N);
		ArrayList<String> labelList = new ArrayList<String>();
		byte snps[][];
		
		System.out.println("Reading " + inputfile + ", writing " + outputfile);
		try {
			File outFile = new File(outputfile);
//			if (outFile.exists()) {
//				String msg = "Output file already exists and will be overwritten.";
//				JOptionPane.showMessageDialog(null, msg, "Output Exists", JOptionPane.ERROR_MESSAGE);
////				return;
//			}
			BufferedWriter bw = new BufferedWriter(new FileWriter(outFile));
			
			//count lines in the input file
			BufferedReader br = new BufferedReader(new FileReader(inputfile));
			int linecount = 0;
			while (br.readLine() != null) linecount++;
			br.close();
			int numberOfSnps = linecount - 1;
			
			String input;
			String[] info;
			br = new BufferedReader(new FileReader(inputfile));
			input = br.readLine();
			info = tab.split(input);
			int numberOfTaxa = info.length - 11;
			snps = new byte[numberOfTaxa][numberOfSnps];
			bw.write(input);
			bw.write("\n");
			
			int count = 0;
			while ((input = br.readLine()) != null) {
//				if (count % 1000 == 0) System.out.println("Processing line " + count);
				info = tab.split(input);
				StringBuilder sb = new StringBuilder(info[0]);
				for (int i = 1; i < 11; i++) sb.append("\t").append(info[i]);
				labelList.add(sb.toString());

				for (int t = 0; t < numberOfTaxa; t++) {
					snps[t][count] = snpToByte.get(info[t+11]);
				}
				count++;
			}
			br.close();
			
			for (int t = 0; t < numberOfTaxa; t++) {
				callHets(snps[t]);
			}
			
			//write to output
			for (int s = 0; s < numberOfSnps; s++) {
				bw.write(labelList.get(s));
				for (int t = 0; t < numberOfTaxa; t++) {
					bw.write("\t");
					bw.write(byteToSnp[snps[t][s]]);
				}
				bw.write("\n");
			}
			bw.close();
			
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		
		System.out.println("Finished calling hets.");
	}
	
	public byte getHetCode(byte a1, byte a2) {
		if (a1 < R & a2 < R) return hetcodes[a1][a2];
		if (a1 == a2) return a1;
		if ((a1 == plus && a2 == minus) || (a1 == plus && a2 == minus)) return zero;
		return N;
	}
	
	public void writeTaxa(String filenamebase) {
		
		for (String popname:popmap.keySet()) {
			ArrayList<Integer> taxaList = popmap.get(popname);
			String filename;
			if (filenamebase == null) filename ="C:/users/peter/temp/gbs/" + popname + "_chr" + chr + "_imputedsnps.txt";
			else filename = filenamebase + "_chr" + chr + "_" + popname + ".txt";
			try {
				BufferedWriter bw = new BufferedWriter(new FileWriter(filename));
				bw.write("rs#\talleles\tchrom\tpos\tstrand\tassembly#\tcenter\tprotLSID\tassayLSID\tpanelLSID\tQCcode");
				for (Integer t:taxaList) bw.write("\t" + taxanames.get(t));
				bw.write("\n");
				for (int i = 0; i < numberOfSnps; i++) {
					Iterator<Integer> it = taxaList.iterator();
					boolean allN = true;
					while (allN && it.hasNext()) {
						if (imputedSnps[it.next()][i] != N) allN = false;
					}
					if (!allN) {
						bw.write(snplabels.get(i));
						for (Integer t:taxaList) {
							bw.write("\t");
							bw.write(byteToSnp[imputedSnps[t][i]]);
						}
						bw.write("\n");
					}
				}
				bw.close();
			} catch (IOException e) {
				e.printStackTrace();
				System.exit(-1);
			}
		}
		
		System.out.println("finished writing files");
	}
	
	public static void countClasses() {
		String filename = "C:/Projects/NAM/namgbs/allfusion_110401.LD/allLD/allfusion_110401.lh.ld.c4.dedupe.txt";
		Pattern tab = Pattern.compile("\t");
		String input;
		String[] info;
		int[][] popcounts = new int[27][2];
		int[] pop;
		
		try {
			BufferedReader br = new BufferedReader(new FileReader(filename));
			info = tab.split(br.readLine());
			int n = info.length;
			pop = new int[n - 11];
			for (int i = 11; i < n; i++) {
				pop[i - 11] = Integer.parseInt(info[i].substring(2,4));
			}
			
			while ((input = br.readLine()) != null) {
				info = tab.split(input);
				for (int i = 11; i < n; i++) {
					int taxon = i - 11;
					int popnum = pop[taxon];
					popcounts[popnum][1]++;
					if (!info[taxon].equals("N")) popcounts[popnum][0]++;
				}
			}
			br.close();
			
			for (int i = 0; i < 27; i++) {
				int totalcount = popcounts[i][1];
				int nonNcount = popcounts[i][0];
				double prop;
				if (totalcount != 0) prop = ((double) nonNcount) / ((double) totalcount);
				else prop = Double.NaN;
				
				System.out.println("pop " + i + ": " + nonNcount + ", " + totalcount + ", " + prop);
			}
			
		} catch(IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		
	}
	
	public static void imputeForHapmap(String inputfile, String outputfile) {
		String input;
		String[] info;
		Pattern tab = Pattern.compile("\t");
		String fill = "\tNA\tNA\tNA\tNA\tNA\tNA\tNA";
		
		byte A = (byte) 0;
		byte H = (byte) 1;
		byte B = (byte) 2;
		byte N = (byte) 3;
		
		HashMap<String, Byte> string2byte = new HashMap<String, Byte>();
		string2byte.put("A", A);
		string2byte.put("H", H);
		string2byte.put("B", B);
		string2byte.put("N", N);
		
		String[] b2s = new String[] {"A", "M", "C", "N"};
		
		try{
			BufferedReader br = new BufferedReader(new FileReader(inputfile));
			BufferedWriter bw = new BufferedWriter(new FileWriter(outputfile));
			
			//count lines
			int linecount = 0;
			while (br.readLine() != null) linecount++;
			br.close();
			int numberOfSnps = linecount - 1;
			
			br = new BufferedReader(new FileReader(inputfile));
			info = tab.split(br.readLine());
			int numberOfTaxa = info.length - 4;
			byte[][] snps = new byte[numberOfTaxa][numberOfSnps];
			ArrayList<String> snplabel = new ArrayList<String>();
			bw.write("rs#\talleles\tchrom\tpos\tstrand\tassembly#\tcentre\tprotLSID\tassayLSID\tpanelLSID\tQCcode");
			for (int i = 4; i < info.length; i++) bw.write("\t" + info[i]);
			bw.write("\n");
			
			System.out.println("reading data");
			linecount = 0;
			while ((input = br.readLine()) != null) {
				info = tab.split(input);
				snplabel.add(info[0] + "\t" + info[1] + "\t" + info[2] + "\t" + info[3]);
				for (int t = 0; t < numberOfTaxa; t++) {
					snps[t][linecount] = string2byte.get(info[t + 4]);
				}
				linecount++;
			}
			
			//impute snps
			System.out.println("imputing snps");
			for (int t = 0; t < numberOfTaxa; t++) {
				byte prevsnp = -1;
				int prevpos = -1;
				byte[] asnp = snps[t];
				for (int s = 0; s < numberOfSnps; s++) {
					if (asnp[s] != N) {
						if (prevsnp == -1) {
							for (int i = 0; i < s; i++) asnp[i] = asnp[s];
						} else if(asnp[s] == prevsnp) {
							for (int i = prevpos; i < s; i++) asnp[i] = prevsnp;
						}
						prevsnp = asnp[s];
						prevpos = s;
					}
				}
				if (prevpos > 0) {
					for (int i = prevpos; i < numberOfSnps; i++) asnp[i] = prevsnp;
				}
			}
			
			System.out.println("writing output");
			for (int s = 0; s < numberOfSnps; s++) {
				bw.write(snplabel.get(s));
				bw.write(fill);
				for (int t = 0; t < numberOfTaxa; t++) {
					bw.write("\t");
					bw.write(b2s[snps[t][s]]);
				}
				bw.write("\n");
			}
			
			
			br.close();
			bw.close();
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}

	public static void imputeNumericGenotype(String inputfile, String outputfile) {
		String input;
		String[] info;
		Pattern tab = Pattern.compile("\t");
		
		byte A = (byte) 0;
		byte H = (byte) 1;
		byte B = (byte) 2;
		byte N = (byte) 3;
		
		HashMap<String, Byte> string2byte = new HashMap<String, Byte>();
		string2byte.put("A", A);
		string2byte.put("H", H);
		string2byte.put("B", B);
		string2byte.put("N", N);
		
		String[] b2s = new String[] {"0", "1", "2", "N"};
		
		try{
			BufferedReader br = new BufferedReader(new FileReader(inputfile));
			BufferedWriter bw = new BufferedWriter(new FileWriter(outputfile));
			
			//count lines
			int linecount = 0;
			while (br.readLine() != null) linecount++;
			br.close();
			int numberOfSnps = linecount - 1;
			
			br = new BufferedReader(new FileReader(inputfile));
			info = tab.split(br.readLine());
			int numberOfTaxa = info.length - 4;
			byte[][] snps = new byte[numberOfTaxa][numberOfSnps];
			String[] taxa = new String[numberOfTaxa];
			for (int t = 0; t < numberOfTaxa; t++) taxa[t] = info[t + 4];
			ArrayList<String> snplabel = new ArrayList<String>();
			bw.write("<Numeric>\n");
			bw.write("<Marker>");
			for (int s = 0; s < numberOfSnps; s++) {
				bw.write("\t");
				bw.write("s" + s);
			}
			bw.write("\n");
			
			System.out.println("reading data");
			linecount = 0;
			while ((input = br.readLine()) != null) {
				info = tab.split(input);
				snplabel.add(info[0] + "\t" + info[1] + "\t" + info[2] + "\t" + info[3]);
				for (int t = 0; t < numberOfTaxa; t++) {
					snps[t][linecount] = string2byte.get(info[t + 4]);
				}
				linecount++;
			}
			
			//impute snps
			System.out.println("imputing snps");
			for (int t = 0; t < numberOfTaxa; t++) {
				byte prevsnp = -1;
				int prevpos = -1;
				byte[] asnp = snps[t];
				for (int s = 0; s < numberOfSnps; s++) {
					if (asnp[s] != N) {
						if (prevsnp == -1) {
							for (int i = 0; i < s; i++) asnp[i] = asnp[s];
						} else if(asnp[s] == prevsnp) {
							for (int i = prevpos; i < s; i++) asnp[i] = prevsnp;
						}
						prevsnp = asnp[s];
						prevpos = s;
					}
				}
				if (prevpos > 0) {
					for (int i = prevpos; i < numberOfSnps; i++) asnp[i] = prevsnp;
				}
			}
			
			System.out.println("writing output");
			
			for (int t = 0; t < numberOfTaxa; t++) {
				bw.write(taxa[t]);
				for (int s = 0; s < numberOfSnps; s++) {
					bw.write("\t");
					byte val = snps[t][s];
					bw.write(b2s[val]);
				}
				bw.write("\n");
			}
			
			
			br.close();
			bw.close();
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}

	/**
	 * compares gbs calls to 1106 marker genotyping array calls
	 * @param inputfile a hapmap formatted file of genotypes for a chromosome with hets imputed
	 * @return a HashMap with taxa names as keys and matchcount, total nam homozygote count, hetcount, totalcount as values
	 */
	public static HashMap<String, int[]> compareToNAMarray(String inputfile, int chr) {
		String arrayFileName = "C:/Projects/NAM/data/markergenotypes062508.txt";
		String mapFileName = "C:/Projects/NAM/data/markers061208agpv2.txt";
		int[] agpMarkerPositions = null;
		int[] namMarkerNumber = null;
		int numberOfMarkers = 0;
		HashMap <String, int[]> genotypeMap = new HashMap <String, int[]>();
		HashMap <String, int[]> taxaCountMap = new HashMap <String, int[]>();
		
		String input;
		String[] data;
		Pattern tab = Pattern.compile("\t");
		
		int startMarker = 0;
		
		HashMap <String, Integer> string2int = new HashMap <String, Integer>();
		string2int.put("0.0", 0);
		string2int.put("0.5", 1);
		string2int.put("1.0", 2);
		string2int.put("1.5", 3);
		string2int.put("2.0", 4);
		
		//read in the map file
		try { 
			
			//count the number of markers in this chromosome
			BufferedReader br = new BufferedReader(new FileReader(mapFileName));
			br.readLine();
			int snpcount = 0;
			int linecount = 0;
			input = br.readLine();
			data = tab.split(input);
			int chrnum = Integer.parseInt(data[1]);
			while (chrnum < chr) {
				input = br.readLine();
				data = tab.split(input);
				chrnum = Integer.parseInt(data[1]);
				linecount++;
			}
			startMarker = linecount;
			
			while (chrnum == chr) {
				snpcount++;
				input = br.readLine();
				if (input != null) {
					data = tab.split(input);
					chrnum = Integer.parseInt(data[1]);
				} else {
					chrnum = -1;
				}
				linecount++;
			}
			numberOfMarkers = snpcount;
			br.close();
			
			//dimension the marker position array and read the agpv2 positions into it
			agpMarkerPositions = new int[numberOfMarkers];
			namMarkerNumber = new int[numberOfMarkers];
			br = new BufferedReader(new FileReader(mapFileName));
			br.readLine();
			snpcount = 0;
			input = br.readLine();
			data = tab.split(input);
			chrnum = Integer.parseInt(data[1]);
			while (chrnum < chr) {
				input = br.readLine();
				data = tab.split(input);
				chrnum = Integer.parseInt(data[1]);
			}
			while (chrnum == chr) {
				agpMarkerPositions[snpcount] = Integer.parseInt(data[5]);
				namMarkerNumber[snpcount++] = Integer.parseInt(data[3]) - 1;
				input = br.readLine();
				if (input != null) {
					data = tab.split(input);
					chrnum = Integer.parseInt(data[1]);
				} else {
					chrnum = -1;
				}
			}
			br.close();
			
			
		} catch(IOException e){
			e.printStackTrace();
			System.exit(-1);
		}
		
		//read in the markergenotype data for this chromosome

		try { 
			BufferedReader br = new BufferedReader(new FileReader(arrayFileName));
			br.readLine();
			br.readLine(); //discard two header rows
			
			while ((input = br.readLine()) != null) {
				data = tab.split(input);
				String taxon = data[0];
				int[] geno = new int[1106];
				for (int m = 0; m < 1106; m++) {
					geno[m] = string2int.get(data[m + 5]);
				}
				genotypeMap.put(taxon, geno);
			}
			
			br.close();
		} catch(IOException e){
			e.printStackTrace();
			System.exit(-1);
		}
		
		//read the snp data and get counts

		try { 
			BufferedReader br = new BufferedReader(new FileReader(inputfile));
			input = br.readLine();
			data = tab.split(input);
			int numberOfTaxa = data.length - 11;
			String[] taxanames = new String[numberOfTaxa];
			for (int t = 0; t < numberOfTaxa; t++) {
				String taxon = data[t + 11];
				taxanames[t] = taxon.substring(0, taxon.indexOf(':'));
				taxaCountMap.put(taxanames[t], new int[4]);
			}
			
			while ((input = br.readLine()) != null) {
				data = tab.split(input);
				//get the position of this snp and find the flanking nam markers
				int snpPosition = Integer.parseInt(data[3]);
				int snpIndex = Arrays.binarySearch(agpMarkerPositions, snpPosition);
				int[] flankingMarkers = new int[2];
				if (snpIndex >= 0) {
					flankingMarkers[0] = snpIndex;
					flankingMarkers[1] = snpIndex;
				} else {
					int insertionPoint = - (snpIndex + 1);
					if (insertionPoint == 0) {
						flankingMarkers[0] = 0;
						flankingMarkers[1] = 0;
					} else if (insertionPoint == numberOfMarkers) {
						flankingMarkers[0] = numberOfMarkers - 1;
						flankingMarkers[1] = numberOfMarkers - 1;
					} else {
						flankingMarkers[0] = insertionPoint - 1;
						flankingMarkers[1] = insertionPoint;
					}
				}
				
				for (int t = 0; t < numberOfTaxa; t++) {
					String geno = data[t + 11];
					int[] namgeno = genotypeMap.get(taxanames[t]);
					if (namgeno != null) {
						int leftmarker = namgeno[namMarkerNumber[flankingMarkers[0]]];
						int rightmarker = namgeno[namMarkerNumber[flankingMarkers[1]]];
						int markersum = leftmarker + rightmarker;
						int[] theseCounts = taxaCountMap.get(taxanames[t]);
						if (geno.equals("A")) {
							theseCounts[3]++;
							if (markersum == 0) {
								theseCounts[0]++;
								theseCounts[1]++;
							} else if (markersum == 8) {
								theseCounts[1]++;
							}
						} else if (geno.equals("B")) {
							theseCounts[3]++;
							if (markersum == 8) {
								theseCounts[0]++;
								theseCounts[1]++;
							} else if (markersum == 0) {
								theseCounts[1]++;
							}
						} else if (geno.equals("H")) {
							theseCounts[2]++;
							theseCounts[3]++;
						} 
					}
				}
			}
			
			br.close();
		} catch(IOException e){
			e.printStackTrace();
			System.exit(-1);
		}
		
		//code to output results to a file
//		try {
//			BufferedWriter bw = new BufferedWriter(new FileWriter("c:/users/peter/temp/gbs/matchcounts.chr" + chr + ".txt"));
//			bw.write("taxon\tmatchAB\ttotalAB\thets\ttotal\n");
//			for (String t : taxaCountMap.keySet()) {
//				bw.write(t);
//				int counts[] = taxaCountMap.get(t);
//				for (int i = 0; i < 4; i++) {
//					bw.write("\t" + counts[i]);
//				}
//				bw.newLine();
//			}
//			bw.close();
//		} catch(IOException e) {
//			
//		}
		return taxaCountMap;
	}
	
	public static void countMatches() {
		int pops[] = new int[]{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,18,19,20,21,22,23,24,25,26};
		for (int pop:pops) {
			try {
				BufferedWriter bw = new BufferedWriter(new FileWriter("C:/Projects/NAM/namgbs/datasets/cov20/matchcounts.pop" + pop + ".txt"));
				bw.write("chromosome\ttaxon\tmatchAB\ttotalAB\tHcount\ttotal\tprmatch\tprhet\n");

				for (int chr = 1; chr <= 10; chr++) {
					String znum;
					if (pop < 10) znum = "Z00" + pop;
					else znum = "Z0" + pop;
					String filename = "C:/Projects/NAM/namgbs/datasets/cov20/cov20_chr" + chr + "_" + znum + ".hets.txt";
					HashMap<String, int[]> taxonCounts = compareToNAMarray(filename, chr);
					for (String t : taxonCounts.keySet()) {
						bw.write(chr + "\t" + t);
						int counts[] = taxonCounts.get(t);
						for (int i = 0; i < 4; i++) {
							bw.write("\t" + counts[i]);
						}
						if (counts[1] > 0) bw.write("\t" + (((double) counts[0])/counts[1]));
						else bw.write("\t");
						if (counts[3] > 0) bw.write("\t" + (((double) counts[2])/counts[3]));
						else bw.write("\t");
						bw.newLine();
					}

				}
				bw.close();
			} catch(IOException e) {
				e.printStackTrace();
				System.exit(-1);
			}
			System.out.println("Finished summarizing population " + pop + ".");
		}
	}
	
	public static void generateSNPdataset(int chr, String outputFileName) {
		class Snp implements Comparable<Snp> {
			String name;
			String allele;
			int chr;
			int pos;
//			double genpos;
			
			
			Snp(String name, String allele, int chr, int pos) {
				this.name = name;
				this.allele = allele;
				this.chr = chr;
				this.pos = pos;
			}

			@Override
			public int compareTo(Snp snp) {
				if (chr > snp.chr) return 1;
				if (chr < snp.chr) return -1;
				if (pos > snp.pos) return 1;
				if (pos < snp.pos) return -1;
				
				int compName = name.compareTo(snp.name);
				if (compName == 0) return allele.compareTo(snp.allele);
				return compName;
				
			}
		}
		
		HashMap<String,Integer> taxaMap = new HashMap<String,Integer>();
		HashMap<Snp, Integer> snpMap = new HashMap<Snp, Integer>();
		String input;
		String[] info;
		Pattern tab = Pattern.compile("\t");
		int numberOfTaxa;
		int numberOfSnps;
		ArrayList<String> taxaList = new ArrayList<String>();
		AGPMap agpmap = new AGPMap(true);

		//read in taxanames and delete the poor quality lines
		String taxaNameFile = "C:/Projects/NAM/NAM_map_and_genos-090921/RILs_for_NAM_Map_20071102.txt";
		try {
			BufferedReader br = new BufferedReader(new FileReader(taxaNameFile));
			info = tab.split(br.readLine());
			for (String t:info) taxaList.add(t);
			br.close();
			
			System.out.println("Taxa names from NAM list: " + taxaList.size());
			String qualityFile = "C:/Projects/NAM/namgbs/datasets/cov20/sample quality issues.txt";
			br = new BufferedReader(new FileReader(qualityFile));
			br.readLine();
			while ((input = br.readLine()) != null) {
				info = tab.split(input);
				taxaList.remove(info[0]);
			}
			System.out.println("Taxa names after removing samples with problems: " + taxaList.size());
			
			br.close();
		} catch(IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		
		//sort the taxa
		Collections.sort(taxaList);
		numberOfTaxa = taxaList.size();
		
		for (int i = 0; i < 10; i++) System.gc();
		long total = Runtime.getRuntime().totalMemory();
		long free = Runtime.getRuntime().freeMemory();
		long used = total - free;
		System.out.println("taxalist memory: " + "total = " + total + ", free = " + free + ", used = " + used);
		
		//read in the snps and sort them for faster searching
		//get a list of snp files
		System.out.println("Getting a list of all snps.");
		File snpdir = new File("C:/Projects/NAM/namgbs/datasets/cov20/");
		final String prefix = "cov20_chr" + chr + "_Z";
		File[] snpfiles = snpdir.listFiles(new FilenameFilter() {

			@Override
			public boolean accept(File dir, String filename) {
				if (filename.startsWith(prefix) && filename.endsWith(".hets.txt")) return true;
				return false;
			}
			
		});
		
		//create a sorted set of all snps (each snp contained only once)
		TreeSet<Snp> snpSet = new TreeSet<Snp>();
		int count = 0;
		for (File file:snpfiles) {
			System.out.print(" " + count++);
			if (count % 40 == 0) System.out.println();
			try {
				BufferedReader br = new BufferedReader(new FileReader(file));
				br.readLine();
				while ((input = br.readLine()) != null) {
					info = tab.split(input);
					Snp snp = new Snp(info[0], info[1], Integer.parseInt(info[2]), Integer.parseInt(info[3]));
					snpSet.add(snp);
				}
				br.close();
			} catch(IOException e) {
				e.printStackTrace();
				System.exit(-1);
			} 
		}
		System.out.println("\nFinished reading snp files.");

		ArrayList<Snp> snpList = new ArrayList<Snp>(snpSet);
		numberOfSnps = snpList.size();
		
		System.out.println("Number of snps = " + numberOfSnps);
		for (int i = 0; i < 10; i++) System.gc();
		total = Runtime.getRuntime().totalMemory();
		free = Runtime.getRuntime().freeMemory();
		used = total - free;
		System.out.println("snplist memory: " + "total = " + total + ", free = " + free + ", used = " + used);
		
		//find the start and end of each chromosome rounded to 0.2 cM
		double startgenpos = agpmap.getCmFromPosition(chr, snpList.get(0).pos);
		double endgenpos = agpmap.getCmFromPosition(chr, snpList.get(snpList.size() - 1).pos);
		
		//create a list of imputed snps at 0.2 cM between these limits (inclusive)
		ArrayList<Snp> imputedSnpList = new ArrayList<Snp>();
		double interval = 0.2;
		
		double start = Math.ceil(startgenpos / interval) * interval;
		double end = Math.floor(endgenpos / interval) * interval;
		String snpPrefix;
		if (chr < 10) {
			snpPrefix = "P0" + chr + "_";
		} else {
			snpPrefix = "P10_";
		}
		
		for (double p = start; p <= end; p += 0.2) {
			int pos = agpmap.getPositionFromCm(chr, p);
			String snpname = snpPrefix + pos;
			imputedSnpList.add(new Snp(snpname, "A/B", chr, pos)); 
		}
		
		int numberOfImputedSnps = imputedSnpList.size();
		System.out.println("Number of imputed snps =  " + numberOfImputedSnps);
		for (int i = 0; i < 10; i++) System.gc();
		total = Runtime.getRuntime().totalMemory();
		free = Runtime.getRuntime().freeMemory();
		used = total - free;
		System.out.println("imputed snplist memory: " + "total = " + total + ", free = " + free + ", used = " + used);

		//impute these snps based on known flanking snps
		
		//dimension the snp array and set to N (=3)
		byte[][] snps = new byte[numberOfSnps][numberOfTaxa];
//		float[][] snps = new float[numberOfSnps][numberOfTaxa];
		for (int s = 0; s < numberOfSnps; s++) {
			for (int t = 0; t < numberOfTaxa; t++) {
				snps[s][t] = 3;
			}
		}
		
		//read through the snp files a second time and load the snp data into the byte array
		//code bytes as A=0, H=1, B=2, N=3
		System.out.println("Reading snp files and loading the snp array:");

		HashMap<String,Byte> snp2byte = new HashMap<String,Byte>();
		snp2byte.put("A", (byte) 0);
		snp2byte.put("H", (byte) 1);
		snp2byte.put("B", (byte) 2);
		snp2byte.put("N", (byte) 3);
		
		
		for (File file:snpfiles) {
			try {
				System.out.println("Processing " + file.getName());
				
				BufferedReader br = new BufferedReader(new FileReader(file));
				info = tab.split(br.readLine());
				int ntaxa = info.length - 11;
				String[] taxanames = new String[ntaxa];
				for (int t = 0; t < ntaxa; t++) {
					taxanames[t] = info[t + 11].substring(0, info[t+11].indexOf(':'));
				}
				
				while ((input = br.readLine()) != null) {
					info = tab.split(input);
					Snp snp = new Snp(info[0], info[1], Integer.parseInt(info[2]), Integer.parseInt(info[3]));
					int ndxSnp = Collections.binarySearch(snpList, snp);
					for (int t = 0; t < ntaxa; t++) {
						int ndxTaxon = Collections.binarySearch(taxaList, taxanames[t]);
						if (ndxTaxon > -1) {
							snps[ndxSnp][ndxTaxon] = snp2byte.get(info[t + 11]);
						}
					}
				}
				br.close();
			} catch(IOException e) {
				e.printStackTrace();
				System.exit(-1);
			} 
		}
		
		for (int i = 0; i < 10; i++) System.gc();
		total = Runtime.getRuntime().totalMemory();
		free = Runtime.getRuntime().freeMemory();
		used = total - free;
		System.out.println("after loading all snps memory: " + "total = " + total + ", free = " + free + ", used = " + used);

		//impute snps
		System.out.println("Creating set of imputed snps.");
		
		float[][] imputedSnps = new float[numberOfImputedSnps][numberOfTaxa];
		for (int t = 0; t < numberOfTaxa; t++) {
			for (int s = 0; s < numberOfImputedSnps; s++) {
				Snp snp = imputedSnpList.get(s);
				int snpIndex = Collections.binarySearch(snpList, snp);
				int endndx = snpList.size() - 1;
				int rightIndex = snpIndex;
				int leftIndex = snpIndex;
				if (snpIndex < 0) {
					rightIndex = -(snpIndex + 1);
					leftIndex = rightIndex - 1;
					if (leftIndex < 0) leftIndex = 0;
					if (rightIndex > endndx) rightIndex = endndx;
				}				
				while (snps[rightIndex][t] == 3 && rightIndex < endndx) rightIndex++;
				while (snps[leftIndex][t] == 3 && leftIndex > 0) leftIndex--;
				if (snps[rightIndex][t] == 3 || rightIndex == leftIndex) imputedSnps[s][t] = snps[leftIndex][t];
				else if (snps[leftIndex][t] == 3) imputedSnps[s][t] = snps[rightIndex][t];
				else {
					Snp leftsnp = snpList.get(leftIndex);
					Snp rightsnp = snpList.get(rightIndex);
					float pd = ((float)(snp.pos - leftsnp.pos)) / ((float)(rightsnp.pos - leftsnp.pos));
					imputedSnps[s][t] = ((float) snps[leftIndex][t]) * (1 - pd) + ((float) snps[rightIndex][t]) * pd;
				}

			}
		}
		
		//write all imputed snp data to a single file
		System.out.println("Writing output to " + outputFileName);
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(outputFileName));
			StringBuilder sb = new StringBuilder("Snp\tallele\tchr\tpos\tcm");
			for (int t = 0; t < numberOfTaxa; t++) {
				sb.append("\t").append(taxaList.get(t));
			}
			bw.write(sb.toString());
			bw.newLine();
			
			for (int s = 0; s < numberOfImputedSnps; s++) {
				Snp snp = imputedSnpList.get(s);
				sb = new StringBuilder(snp.name);
				sb.append("\t").append(snp.allele);
				sb.append("\t").append(snp.chr);
				sb.append("\t").append(snp.pos);
				sb.append("\t").append(agpmap.getCmFromPosition(snp.chr, snp.pos));
				for (int t = 0; t < numberOfTaxa; t++) {
					sb.append("\t");
					sb.append(imputedSnps[s][t]);
				}
				bw.write(sb.toString());
				bw.newLine();
			}
			bw.close();
		} catch(IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		System.out.println("Finished.");
	}
	
	public static void createImputedSNPDatasets() {
		for (int c = 10; c <= 10; c++) {
			String filename = "C:/Projects/NAM/namgbs/datasets/cov20/imputedSnpsChr" + c + ".txt";
			generateSNPdataset(c, filename);
		}
	}
	
	public static void concatenateSomeFiles() {
		File[] input = new File[10];
		
		for (int chr = 1; chr <= 10; chr++) {
			input[chr - 1] = new File("C:/Projects/NAM/namgbs/datasets/cov20/imputedSnpsChr" + chr + ".txt");
		}
		File output = new File("C:/Projects/NAM/namgbs/datasets/cov20/imputedSnpsAllChr.txt");
		concatenateOutput(input, output);
	}
	
	public static void concatenateOutput(File[] inputFiles, File outputFile) {
		//assumes that the first line consists of column headers and that the headers are the same for all files
		int numberOfInputFiles = inputFiles.length;
		byte lf = (byte) '\n';
		
		try {
			FileOutputStream fos = new FileOutputStream(outputFile);
			FileChannel outfc = fos.getChannel();
			
			System.out.println("Processing " + inputFiles[0].getPath());
			FileInputStream fis = new FileInputStream(inputFiles[0]);
			FileChannel infc = fis.getChannel();
			int bufferCapacity = 100000;
			ByteBuffer bb = ByteBuffer.allocate(bufferCapacity);
			bb.clear();
			
			while (infc.read(bb) > 0) {
				bb.flip();
				outfc.write(bb);
				bb.clear();
			}
			infc.close();
			
			for (int f = 1; f < numberOfInputFiles; f++) {
				System.out.println("Processing " + inputFiles[f].getPath());
				fis = new FileInputStream(inputFiles[f]);
				infc = fis.getChannel();
				bb.clear();
				int bytesread = infc.read(bb);
				bb.flip();
				byte b = bb.get();
				while(b != lf) {
					b = bb.get();
				}
				outfc.write(bb);
				bb.clear();
				while (infc.read(bb) > 0) {
					bb.flip();
					outfc.write(bb);
					bb.clear();
				}
				infc.close();
			}

			outfc.close();
			
		} catch(IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	
}
