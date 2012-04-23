package net.maizegenetics.gwas.imputation;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.regex.Pattern;

import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.FilterAlignment;
import net.maizegenetics.pal.alignment.ImportUtils;
import net.maizegenetics.pal.alignment.MutableNucleotideAlignment;
import net.maizegenetics.pal.alignment.NucleotideAlignmentConstants;
import net.maizegenetics.pal.alignment.SBitAlignment;
import net.maizegenetics.pal.alignment.TBitAlignment;
import net.maizegenetics.pal.distance.DistanceMatrix;
import net.maizegenetics.pal.distance.IBSDistanceMatrix;
import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.IdGroupUtils;
import net.maizegenetics.pal.ids.SimpleIdGroup;
import net.maizegenetics.pal.tree.Tree;
import net.maizegenetics.pal.tree.TreeClusters;
import net.maizegenetics.pal.tree.UPGMATree;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.BitUtil;
import net.maizegenetics.util.OpenBitSet;
import net.maizegenetics.util.ProgressListener;

public class NucleotideImputor {
	/*
	 * Imputations steps for individual populations:
	 * 1. identify parent haplotypes in a window
	 * 2. add snps in LD with the parent haplotypes
	 * 3. identify parent alleles
	 * 4. score progeny for parent alleles
	 * 5. use Viterbi algorithm
	 * 
	 * */
	
	Pattern tab = Pattern.compile("\t");
	SBitAlignment gbsSnps;
	TBitAlignment gbsSnpsByTaxa;
	HashMap<String, Population> familyMap = new HashMap<String, Population>();
	String baseOutFilename;
	String assembly = "NA";
	int minAlleleCount = 20;
	
	static final byte AA = NucleotideAlignmentConstants.getNucleotideDiploidByte("AA");
	static final byte CC = NucleotideAlignmentConstants.getNucleotideDiploidByte("CC");
	static final byte GG = NucleotideAlignmentConstants.getNucleotideDiploidByte("GG");
	static final byte TT = NucleotideAlignmentConstants.getNucleotideDiploidByte("TT");
	static final byte AC = NucleotideAlignmentConstants.getNucleotideDiploidByte("AC");
	static final byte AG = NucleotideAlignmentConstants.getNucleotideDiploidByte("AC");
	static final byte AT = NucleotideAlignmentConstants.getNucleotideDiploidByte("AC");
	static final byte CG = NucleotideAlignmentConstants.getNucleotideDiploidByte("AC");
	static final byte CT = NucleotideAlignmentConstants.getNucleotideDiploidByte("AC");
	static final byte GT = NucleotideAlignmentConstants.getNucleotideDiploidByte("AC");
	
	static final byte[] byteval = new byte[] {AA,CC,GG,TT,AC};
	final static HashMap<Byte, Integer> genotypeMap = new HashMap<Byte, Integer>();
	
	static{
		genotypeMap.put(AA, 0);
		genotypeMap.put(CC, 1);
		genotypeMap.put(GG, 2);
		genotypeMap.put(TT, 3);
		genotypeMap.put(AC, 4);
	}
	
	static final byte[][] genoval = new byte[][]{{AA,AC,AG,AT},{AC,CC,CG,CT},{AG,CG,GG,GT},{AT,CT,GT,TT}};
	
	class Population {
		ArrayList<String> members;
		BitSet index;
		String parent1;
		String parent2;
		double contribution1;
		double contribution2;
		int Fgen;
		double inbredCoef;
	}
	
	public static void main(String[] args) {
		NucleotideImputor ni = new NucleotideImputor();
		ni.importChromosome("10");
	}
	
	public void importChromosome(String chr) {
		System.out.println("Imputing nucleotides for chromosome " + chr);
		String inputFile = "/Volumes/Macintosh HD 2/data/namgbs/genos_20120110/merged_nam_ibm/NAM_IBM_282_Zea20120110_scv10mF8maf002_mgs_E1pLD5kpUn.c" + chr + ".hmp.txt";
		try {
			BufferedReader br = new BufferedReader(new FileReader(inputFile));
			br.readLine();
			String[] info = br.readLine().split("\t", 12);
			assembly = info[5];
			br.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		int n = inputFile.indexOf(".hmp.txt");
		baseOutFilename = inputFile.substring(0, n);
		gbsSnps = (SBitAlignment) ImportUtils.readFromHapmap(inputFile, null);
		
		System.out.println("Importing population information...");
		importPopulationInformation("/Volumes/Macintosh HD 2/data/namgbs/genos_20120110/merged_nam_ibm/namibm.pedigree.info.txt");
		
		System.out.println("Scoring parents...");
		scoreParentsForAPopulation("1");
		
	}
	
	public void importPopulationInformation(String filename) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(filename));
			br.readLine();
			String input;
			while ((input = br.readLine()) != null) {
				String[] info = tab.split(input);
				Population family = familyMap.get(info[0]);
				if (family == null) {
					family = new Population ();
					family.members = new ArrayList<String>();
					family.members.add(info[2]);  //add parents to family members
					family.members.add(info[3]);
					family.members.add(info[1]);
					family.parent1 = info[2];
					family.parent2 = info[3];
					family.contribution1 = Double.parseDouble(info[4]);
					family.contribution2 = Double.parseDouble(info[5]);
					familyMap.put(info[0], family);
				}
				else family.members.add(info[1]);
			}
			br.close();
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
		
		createPopulationIndex();
	}
	
	public void createPopulationIndex() {
		int ntaxa = gbsSnps.getSequenceCount();
		HashMap<String, Integer> taxaMap = new HashMap<String, Integer>();
		for (int t = 0; t < ntaxa; t++) taxaMap.put(gbsSnps.getTaxaName(t), t);
		
		for (Population pop:familyMap.values()) {
			OpenBitSet mask = new OpenBitSet(ntaxa);
			for (String name:pop.members) {
				try {
					Integer ndx = taxaMap.get(name);
					if (ndx != null) mask.set(taxaMap.get(name));
				} catch(Exception e) {
					System.out.println("Error at name = " + name);
					e.printStackTrace();
					System.exit(-1);
				}
			}
			pop.index = mask;
		}
	}
	
	public double computeR(int site1, int site2, Alignment a) {
		int s1Count = 0;
		int s2Count = 0;
		int prodCount = 0;
		int totalCount = 0;
		
		long[] m11 = a.getAllelePresenceForAllTaxa(site1, 0).getBits();
		long[] m12 = a.getAllelePresenceForAllTaxa(site1, 1).getBits();
		long[] m21 = a.getAllelePresenceForAllTaxa(site2, 0).getBits();
		long[] m22 = a.getAllelePresenceForAllTaxa(site2, 1).getBits();
		int n = m11.length;
		for (int i = 0; i < n; i++) {

			long valid = (m11[i] ^ m12[i]) & (m21[i] ^ m22[i]); //only count non-het & non-missing
			long s1major = m11[i] & valid;
			long s2major = m21[i] & valid;
			long s1s2major = s1major & s2major & valid;
			s1Count += BitUtil.pop(s1major);
			s2Count += BitUtil.pop(s2major);
			prodCount += BitUtil.pop(s1s2major);
			totalCount += BitUtil.pop(valid);
		}
		
		if (totalCount < 2) return Double.NaN;
		//if major site one is x=1, minor site one is x = 0 and for site 2 y = 1 or 0
		// r = [sum(xy) - sum(x)sum(y)/N] / sqrt[(sum(x) - sum(x)*sum(x)/N) * ((sum(y) - sum(y)*sum(y)/N)]
		//and sum(x) - sum(x)*sum(x)/N = sum(x)(N - sum(x))/N
		// because sum(x^2) = sum(x)
		double num = ((double) prodCount - ((double) s1Count * s2Count) / ((double) totalCount));
		double denom = ((double) (s1Count * (totalCount - s1Count))) / ((double) totalCount);
		denom *= ((double) (s2Count * (totalCount - s2Count))) / ((double) totalCount);
		if (denom == 0) return Double.NaN;
		return  num / Math.sqrt(denom);
	}
	
	/*
	 * Strategy a single population
	 * 1. find parental haplotypes within populations
	 * 	a. remove monomorphic snps
	 * 	a. choose a window
	 * 	a. create an alignment for this window and these taxa
	 * 	b. within the window, cluster snps (UPGMA)
	 * 	c. cut the tree at 0.3
	 * 	d. use snps in largest cluster to define parental haplotypes 
	 * 	e. extend haplotypes to ends of chromosome
	 * 2. assign parental haplotypes to specific parents
	 * 3. score progeny
	 * 
	 * */
	public void scoreParentsForAPopulation(String family) {
		Population pop = familyMap.get(family);
		String[] ids = new String[pop.members.size()];
		pop.members.toArray(ids);
		Alignment popAlignment = FilterAlignment.getInstance(gbsSnps, new SimpleIdGroup(ids), false);
		
		//which sites are polymorphic? minor allele count > 2 and exceed the minimum allele count
		int nsites = popAlignment.getSiteCount();
		OpenBitSet polybits = new OpenBitSet(nsites);
		for (int s = 0; s < nsites; s++) {
			int[][] freq = popAlignment.getAllelesSortedByFrequency(s);
			if (freq[1].length > 1 && freq[1][1] > 2) {
				int alleleCount = freq[1][0] + freq[1][1];
				if (alleleCount >= minAlleleCount) polybits.fastSet(s);
			}
		}
		int totalpoly = (int) polybits.cardinality();
		System.out.println("polymorphic sites = " + totalpoly);
		
		//define a window
		int windowSize = 100;
		
		//find windowSize polymorphic snps centered on the midpoint
		int mid = totalpoly / 2;
		int polyCount = 0;
		int start = 0;
		while (polyCount < mid) { //find the midpoint of the polymorphic snps
			if (polybits.get(start++)) polyCount++;
		}
		
		polyCount = 0;
		while (polyCount < windowSize) { //find the starting point of the window
			if (polybits.get(start--)) polyCount++;
		}
		start++;
		
		int[] snpIds = new int[windowSize];
		polyCount = 0;
		int snpid = start;
		while (polyCount < windowSize) {
			if (polybits.get(snpid)) {
				snpIds[polyCount++] = snpid; 
			}
			snpid++;
		}
		
		FilterAlignment filteredPopAlignment = FilterAlignment.getInstance(popAlignment, snpIds);
		
		
		//cluster polymorphic snps within the window by creating a UPGMA tree (cluster on snps)
		SBitAlignment haplotypeAlignment = SBitAlignment.getInstance(filteredPopAlignment);
		UPGMATree myTree = new UPGMATree(snpDistance(haplotypeAlignment));
		
		//cut the tree
		TreeClusters clusterMaker = new TreeClusters(myTree);
		int[] groups = clusterMaker.getGroups(0.3);
		
		//find the biggest group
		int maxGroup = 0;
		for (int grp:groups) maxGroup = Math.max(maxGroup, grp);
		int ngroups = maxGroup + 1;
		int[] groupCount = new int[ngroups];
		for (int grp:groups) groupCount[grp]++;
		int biggestGroup = 0;
		for (int i = 1; i < ngroups; i++) {
			if (groupCount[i] > groupCount[biggestGroup]) biggestGroup = i;
		}
		
		//some diagnostics??
		 //List groups
		for (int i = 0; i < ngroups; i++) {
			if (groupCount[i] > 5) System.out.println("Snp group " + i + " has " + groupCount[i] + " members.");
		}

		//this group of Snps is the starting point
		OpenBitSet ldbits = new OpenBitSet(nsites);
		int ngrp = groupCount[biggestGroup];
		int[] ldSnps = new int[ngrp];
		int grpCount = 0;
		for (int i = 0; i < windowSize; i++) {
			if (groups[i] == biggestGroup) ldSnps[grpCount++] = snpIds[i];
			ldbits.fastSet(snpIds[i]);
		}
		
		//cluster taxa for these snps to find parental haplotypes (cluster on taxa)
		filteredPopAlignment = FilterAlignment.getInstance(popAlignment, ldSnps);
		IBSDistanceMatrix dm = new IBSDistanceMatrix(SBitAlignment.getInstance(filteredPopAlignment));
		estimateMissingDistances(dm);
		myTree = new UPGMATree(dm);
		clusterMaker = new TreeClusters(myTree);
		groups = clusterMaker.getGroups(0.3);
		
		//parents should be in different clusters, there should only be two major clusters
		//find biggest groups
		maxGroup = 0;
		for (int grp:groups) maxGroup = Math.max(maxGroup, grp);
		ngroups = maxGroup + 1;
		groupCount = new int[ngroups];
		for (int grp:groups) groupCount[grp]++;
		
		 //List groups
		for (int i = 0; i < ngroups; i++) {
			if (groupCount[i] > 5) System.out.println("Taxa group " + i + " has " + groupCount[i] + " members.");
		}
		
		//use haplotypes to score parental type
		SBitAlignment sbitPopAlignment = SBitAlignment.getInstance(popAlignment);
		MutableNucleotideAlignment parentAlignment = MutableNucleotideAlignment.getInstance(sbitPopAlignment);
		
		//find the biggest and next biggest groups
		int majorGroup = 0;
		int minorGroup = 1;
		for (int i = 1; i < ngroups; i++) {
			if (groupCount[i] > groupCount[majorGroup]) {
				minorGroup = majorGroup;
				majorGroup = i;
			} else if (groupCount[i] > groupCount[minorGroup]) {
				minorGroup = i;
			}
		}
		
		//create an alignment for each cluster
		IdGroup allTaxa = filteredPopAlignment.getIdGroup();
		int ntaxa = allTaxa.getIdCount();
		boolean[] isInMajor = new boolean[ntaxa];
		boolean[] isInMinor = new boolean[ntaxa];
		for (int i = 0; i < ntaxa; i++) {
			isInMajor[i] = false;
			isInMinor[i] = false;
			if (groups[i] == majorGroup) isInMajor[i] = true; 
			else if (groups[i] == minorGroup) isInMinor[i] = true;
		}
		IdGroup majorTaxa = IdGroupUtils.idGroupSubset(allTaxa, isInMajor);
		IdGroup minorTaxa = IdGroupUtils.idGroupSubset(allTaxa, isInMinor);
		if (majorTaxa.whichIdNumber(pop.parent1) > -1) {
			//nothing to do
		} else if (minorTaxa.whichIdNumber(pop.parent1) > -1) {
			//swap groups
			IdGroup temp = majorTaxa;
			majorTaxa = minorTaxa;
			minorTaxa = temp;
		} else {
			System.out.println("Parent 1 not in a taxa cluster, parent assignment will be arbitrary");
		} 
		
		Alignment aAlignment = FilterAlignment.getInstance(filteredPopAlignment, majorTaxa);
		Alignment cAlignment = FilterAlignment.getInstance(filteredPopAlignment, minorTaxa);
		
		//set first parent to AA, second parent to CC for ldSnps (snps used to form taxa clusters)
		byte AA = NucleotideAlignmentConstants.getNucleotideDiploidByte("AA");
		byte AC = NucleotideAlignmentConstants.getNucleotideDiploidByte("AC");
		byte CC = NucleotideAlignmentConstants.getNucleotideDiploidByte("CC");
		byte missing = NucleotideAlignmentConstants.getNucleotideDiploidByte("NN");
		for (int i = 0; i < ldSnps.length; i++) {
			int snp = ldSnps[i];
			byte parentA = aAlignment.getMajorAllele(snp);
			byte parentC = cAlignment.getMajorAllele(snp);
			
			for (int t = 0; t < ntaxa; t++) {
				byte[] taxon = popAlignment.getBaseArray(t, snp);
				if (taxon[0] == taxon[1]) {
					if (taxon[0] == parentA) parentAlignment.setBase(t, snp, AA);
					else if (taxon[0] == parentC) parentAlignment.setBase(t, snp, CC);
					else parentAlignment.setBase(t, snp, missing);
				} else if (taxon[0] == parentA) {
					if (taxon[1] == parentC) parentAlignment.setBase(t, snp, AC);
					else parentAlignment.setBase(t, snp, missing);
				} else if (taxon[0] == parentC) {
					if (taxon[1] == parentA) parentAlignment.setBase(t, snp, AC);
					else parentAlignment.setBase(t, snp, missing);
				} else {
					parentAlignment.setBase(t, snp, missing);
				}
			}
		}
		
		
		//extend haplotypes
		//add snps in ld
		double minr = 0.9;
		int testSize = 25;
		
		//add snps from middle to start
		LinkedList<Integer> testSnps = new LinkedList<Integer>();
		for (int i = 0; i < testSize; i++) testSnps.add(snpIds[i]);
		for (int testsnp = snpIds[0] - 1; testsnp >= 0; testsnp--) {
			if (polybits.get(testsnp)) {
				double sumr = 0;
				double count = 0;
				for (Integer snp:testSnps) {
					double r = computeR(testsnp, snp, sbitPopAlignment);
					if (!Double.isNaN(r)) {
						sumr += Math.abs(r);
						count++;
					}
				}
				if (count > 0.9) sumr /= count;
				if (sumr > minr && count > 5) {
					ldbits.fastSet(testsnp);
					testSnps.remove();
					testSnps.add(testsnp);
				}
			}
		}
		
		//add snps from middle to end
		testSnps.clear();
		for (int i = testSize; i >= 1; i--) testSnps.add(snpIds[windowSize - i]);
		for (int testsnp = snpIds[windowSize - 1] + 1; testsnp < nsites; testsnp++) {
			if (polybits.get(testsnp)) {
				double sumr = 0;
				double count = 0;
				for (Integer snp:testSnps) {
					double r = computeR(testsnp, snp, sbitPopAlignment);
					if (!Double.isNaN(r)) {
						sumr += Math.abs(r);
						count++;
					}
				}
				if (count > 0.9) sumr /= count;
				if (sumr > minr && count > 5) {
					ldbits.fastSet(testsnp);
					testSnps.remove();
					testSnps.add(testsnp);
				}
			}
		}
		
		//set ldbits to true for snps in the LD set
		System.out.println("Snps added to the LD set = " + ldbits.cardinality());
		int nsnps = (int) ldbits.cardinality(); //nsnps is now the number of snps that will be used
		int[] snpIndex = new int[nsnps];
		int ldCount = 0;
		int n = popAlignment.getSiteCount();
		for (int s = 0; s < n; s++) {
			if (ldbits.fastGet(s)) snpIndex[ldCount++] = s;
		}
		
		//get the parental haplotypes
		SBitAlignment ldSnpAlignment = SBitAlignment.getInstance(FilterAlignment.getInstance(popAlignment, snpIndex));
		BitSet parent0 = new OpenBitSet(nsnps); //records whether a parent carries the major allele
		parent0.fastSet(0); //arbitrarily say that parent 0 is the parent carrying the major allele at site 0.
		
		BitSet[][] mySites = new BitSet[2][2];
		int[][] matches = new int[2][2];
		mySites[0][0] = ldSnpAlignment.getAllelePresenceForAllTaxa(0, 0);
		mySites[0][1] = ldSnpAlignment.getAllelePresenceForAllTaxa(0, 1);
		for (int s = 1; s < nsnps; s++) {
			mySites[1][0] = ldSnpAlignment.getAllelePresenceForAllTaxa(s, 0);
			mySites[1][1] = ldSnpAlignment.getAllelePresenceForAllTaxa(s, 1);
			for (int i = 0; i < 2; i++) {
				for (int j = 0; j < 2; j++) {
					matches[0][0] = (int) OpenBitSet.intersectionCount(mySites[0][i], mySites[1][j]);
				}
			}
			
			if (matches[0][0] + matches[1][1] > matches[0][1] + matches[1][0]) {
				//if major alleles alleles at site s associate with major alleles at site s-1, parent values for s should match s-1.
				if (parent0.fastGet(s - 1)) parent0.fastSet(s);
			} else {
				if (!parent0.fastGet(s - 1)) parent0.fastSet(s);
			}
			
			mySites[0][0] = mySites[1][0];
			mySites[0][1] = mySites[1][1];
		}
		
		//get parental haplotypes
		//cluster taxa using first
		
		//score the progeny as parent 0(2), parent 1(0), het(1), or missing(-1), i.e. count of parent 0 alleles
		nsites = snpIndex.length;
		ntaxa = popAlignment.getSequenceCount();
		byte[][] scores = new byte[nsites][ntaxa];
		for (int i = 0; i < nsites; i++) { //set to missing
			for (int j = 0; j < ntaxa; j++) {
				scores[i][j] = -1;
			}
		}
		
		TBitAlignment taxaAlignment = TBitAlignment.getInstance(ldSnpAlignment);
		for (int t = 0; t < ntaxa; t++) {
			BitSet major = taxaAlignment.getAllelePresenceForAllSites(t, 0);
			BitSet minor = taxaAlignment.getAllelePresenceForAllSites(t, 1);
			for (int s = 0; s < nsites; s++) {
				boolean ismajor = major.fastGet(s);
				boolean isminor = minor.fastGet(s);
				if (ismajor && isminor) scores[s][t] = 1;
				else if (!ismajor & !isminor) scores[s][t] = -1;
				else if (parent0.fastGet(s)) {
					if (ismajor) scores[s][t] = 2;
					else scores[s][t] = 0;
				} else {
					if (isminor) scores[s][t] = 2;
					else scores[s][t] = 0;
				}
			}
		}
		
		writeScoresNumeric(ldSnpAlignment, scores, family);
		
		System.out.println("Finished.");
	}
	
	private boolean addSnpToTarget(int snp, LinkedList<Integer> testSnps, MutableNucleotideAlignment targetAlignment) {
		
		int ntaxa = targetAlignment.getSequenceCount();
		byte[] snpvals = new byte[ntaxa];
		for (int t = 0; t < ntaxa; t++) {
			snpvals[t] = targetAlignment.getBase(t, snp);
		}
		
		int[] acount = new int[5];
		int[] ccount = new int[5];
		for (Integer testsnp:testSnps) {
			for (int t = 0; t < ntaxa; t++) {
				if (snpvals[t] == AA) {
					Integer ndx = genotypeMap.get(targetAlignment.getBase(t, testsnp));
					if (ndx != null) acount[ndx]++;
				} else if (snpvals[t] == CC) {
					Integer ndx = genotypeMap.get(targetAlignment.getBase(t, testsnp));
					if (ndx != null) ccount[ndx]++;
				}
			}
		}
		
		int maxa = 0;
		int maxc = 0;
		for (int i = 1; i < 4; i++) {
			if (acount[i] > acount[maxa]) maxa = i;
			if (ccount[i] > ccount[maxc]) maxc = i;
		}
		
		return false;
		
	}
	
	private double calculateR(int[][] counts) {
		int N = 0;
		for (int[] row:counts) {
			for (int cell:row) {
				N += cell;
			}
		}
		
		double sumx = counts[0][0] + counts[0][1];
		double sumy = counts[0][0] + counts[1][0];
		double sumxy = counts[0][0];
		
		return (sumxy - sumx * sumy / N) / Math.sqrt( (sumx - sumx * sumx / N) * (sumy - sumy * sumy / N) );
		
	}
	
	public void imputeUsingViterbi(byte[][] scores, Alignment a) {
		
	}
	
	public void writeScoresHapmap(Alignment a, byte[][] scores, String family) {
		String outfile = baseOutFilename + "_family" + family + ".hmp.txt";
		int ntaxa = a.getSequenceCount();
		int nsites = a.getSiteCount();
		String[] code = new String[]{"A","M","C"};
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(outfile));
			bw.write("rs#\talleles\tchrom\tpos\tstrand\tassembly#\tcenter\tprotLSID\tassayLSID\tpanelLSID\tQCcode");
			
			for (int t = 0; t < ntaxa; t++) {
				bw.write("\t");
				bw.write(a.getTaxaName(t));
			}
			
			bw.newLine();
			
			for (int s = 0; s < nsites; s++) {
				bw.write(a.getSNPID(s));
				bw.write("\t");
				bw.write(a.getMajorAlleleAsString(s));
				bw.write("/");
				bw.write(a.getMinorAlleleAsString(s));
				bw.write("\t");
				bw.write(a.getLocus(s).getChromosomeName());
				bw.write("\t");
				bw.write(Integer.toString(a.getPositionInLocus(s)));
				bw.write("\tNA\t");
				bw.write(assembly);
				bw.write("\tNA\tNA\tNA\tNA\tNA");
				for (int t = 0; t < ntaxa; t++) {
					bw.write("\t");
					if (scores[s][t] == -1) bw.write("N");
					else bw.write(code[scores[s][t]]);
				}
				bw.newLine();
			}
			
			bw.close(); 
		} catch (IOException e) {
			e.printStackTrace();
			
		}
	}
	
	public void writeScoresNumeric(Alignment a, byte[][] scores, String family) {
		String outfile = baseOutFilename + "_family" + family + ".numeric.txt";
		int ntaxa = a.getSequenceCount();
		int nsites = a.getSiteCount();
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(outfile));
			bw.write("chr\tpos\talleles");
			for (int t = 0; t < ntaxa; t++) {
				bw.write("\t");
				bw.write(a.getTaxaName(t));
			}
			
			bw.newLine();
			
			for (int s = 0; s < nsites; s++) {
				bw.write(a.getLocus(s).getChromosomeName());
				bw.write("\t");
				bw.write(Integer.toString(a.getPositionInLocus(s)));
				bw.write("\t");
				bw.write(a.getMajorAlleleAsString(s));
				bw.write("/");
				bw.write(a.getMinorAlleleAsString(s));
				for (int t = 0; t < ntaxa; t++) {
					bw.write("\t");
					bw.write(Integer.toString(scores[s][t]));
				}
				bw.newLine();
			}
			
			bw.close(); 
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public DistanceMatrix snpDistance(Alignment a) {
		int nsnps = a.getSiteCount();
		SimpleIdGroup snpIds = new SimpleIdGroup(nsnps, true);
		double[][] distance = new double[nsnps][nsnps];
		double sum = 0;
		int count = 0;
		for (int i = 0; i < nsnps; i++) {
			for (int j = i; j < nsnps; j++) {
				distance[i][j] = distance[j][i] = 1 - Math.abs(computeR(i, j, a));
				if (!Double.isNaN(distance[i][j])) {
					sum += distance[i][j];
					count++;
				}
			}
		}
		
		//set missing to average
		double avg = sum/count;
		for (int i = 0; i < nsnps; i++) {
			for (int j = i; j < nsnps; j++) {
				if (Double.isNaN(distance[i][j])) {
					distance[i][j] = distance[j][i] = avg;
				}
			}
		}
		
		return new DistanceMatrix(distance, snpIds);
	}
	
	public void estimateMissingDistances(DistanceMatrix dm) {
		int nsize = dm.getSize();
		
		//average distance
		double totalDistance = 0;
		int count = 0;
		for (int i = 0; i < nsize; i++) {
			for (int j = i + 1; j < nsize; j++) {
				double distance = dm.getDistance(i, j);
				if (!Double.isNaN(distance)) {
					totalDistance += dm.getDistance(i, j);
					count++;
				}
			}
		}
		double avgDist = totalDistance / count;
		
		for (int i = 0; i < nsize; i++) {
			if ( Double.isNaN(dm.getDistance(i,i)) ) dm.setDistance(i, i, 0);
			for (int j = i + 1; j < nsize; j++) {
				if ( Double.isNaN(dm.getDistance(i,j)) ) {
					dm.setDistance(i, j, avgDist);
				}
			}
		}
	}
}
