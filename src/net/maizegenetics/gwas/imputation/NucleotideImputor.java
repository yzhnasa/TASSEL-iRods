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
	 * Imputations steps:
	 * 1. score parental alleles
	 * 2. check LD; eliminate SNPs not in LD with neighbors (within a population)
	 * 3. use Viterbi algorithm
	 * 
	 * 
	 * */
	
	Pattern tab = Pattern.compile("\t");
	SBitAlignment gbsSnps;
	TBitAlignment gbsSnpsByTaxa;
	HashMap<String, Population> familyMap = new HashMap<String, Population>();
	String baseOutFilename;
	String assembly = "NA";
	
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
		
		importPopulationInformation("/Volumes/Macintosh HD 2/data/namgbs/genos_20120110/merged_nam_ibm/namibm.pedigree.info.txt");
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
		
		//which sites are polymorphic? minor allele count > 1
		int nsites = popAlignment.getSiteCount();
		OpenBitSet polybits = new OpenBitSet(nsites);
		for (int s = 0; s < nsites; s++) {
			int[][] freq = popAlignment.getAllelesSortedByFrequency(s);
			if (freq[1].length > 1 && freq[1][1] > 1) {
				polybits.fastSet(s);
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
		
		//assign to parents
		//are the parents in the tree? which are they?
		int[] parents = new int[2];
		parents[0] = parents[1] = -1;
		int ntaxa = myTree.getIdCount();
		int tcount = 0;
		while (tcount < ntaxa && (parents[0] == -1 || parents[1] == -1)) {
			if (myTree.getIdentifier(tcount).compareTo(pop.parent1) == 0) parents[0] = tcount;
			else if (myTree.getIdentifier(tcount).compareTo(pop.parent2) == 0) parents[1] = tcount;
			tcount++;
		}
		
		int[] parentGroup = new int[2];
		parentGroup[0] = groups[parents[0]];
		parentGroup[1] = groups[parents[1]];
		System.out.println("parent 1, " + pop.parent1 + ", is in group " + parentGroup[0] + ".");
		System.out.println("parent 2, " + pop.parent2 + ", is in group " + parentGroup[1] + ".");
		
		//extend haplotypes
		//add snps in ld
		double minr = 0.9;
		int testSize = 25;
		
		//add snps from middle to start
		SBitAlignment sbitPopAlignment = SBitAlignment.getInstance(popAlignment);
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
		int nsnps = (int) ldbits.cardinality();
		int[] snpIndex = new int[nsnps];
		int ldCount = 0;
		for (int s = 0; s < nsnps; s++) {
			if (ldbits.fastGet(s)) snpIndex[ldCount++] = s;
		}
		
		//get the parental haplotypes
		filteredPopAlignment = FilterAlignment.getInstance(popAlignment, snpIndex);
		IdGroup allTaxa = popAlignment.getIdGroup();
		ntaxa = allTaxa.getIdCount();
		boolean[] p0group = new boolean[ntaxa];
		boolean[] p1group = new boolean[ntaxa]; 
		for (int t = 0; t < ntaxa; t++) {
			if (groups[t] == parentGroup[1]) p1group[t] = true;
			else p1group[t] = false;
		}
		IdGroup p1IdGroup = IdGroupUtils.idGroupSubset(allTaxa, p1group);
		Alignment p1Alignment = FilterAlignment.getInstance(filteredPopAlignment, p1IdGroup);
		
		//score the progeny as parent 0(0), parent 1(2), het(1), or missing(-1)
		nsites = snpIndex.length;
		byte[][] scores = new byte[nsites][ntaxa];
		for (int s = 0; s < nsites; s++) {
			byte p1 = p1Alignment.getMajorAllele(s);
			byte p0;
			byte major = filteredPopAlignment.getMajorAllele(s);
			byte minor = filteredPopAlignment.getMinorAllele(s);
			if (p1 == major) p0 = minor;
			else p0 = major;
			for ( int t = 0; t < ntaxa; t++) {
				if (filteredPopAlignment.isHeterozygous(t, s)) {
					scores[s][t]= 1;
				} else {
					byte[] alleles = filteredPopAlignment.getBaseArray(t, s);
					if (alleles[0] == Alignment.UNKNOWN_ALLELE || alleles[1] == Alignment.UNKNOWN_ALLELE) scores[s][t] = -1;
					else if (filteredPopAlignment.isHeterozygous(t, s)) scores[s][t]= 1;
					else if (alleles[0] == p1) scores[s][t] = 2;
					else if (alleles[0] == p0) scores[s][t] = 0;
					else scores[s][t] = -1;
				}
			}
		}
		
		
		
	}
	
	public void imputeUsingViterbi(byte[][] scores, Alignment a) {
		
	}
	
	public void writeScoresHapmap(Alignment a, byte[][] scores, String family) {
		String outfile = baseOutFilename + "family:" + family + ".hmp.txt";
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
				bw.write(a.getPositionInLocus(s));
				bw.write("\tNA\t");
				bw.write(assembly);
				bw.write("\tNA\tNA\tNA\tNA\tNA");
				for (int t = 0; t < ntaxa; t++) {
					bw.write("\t");
					if (scores[s][t] == -1) bw.write("N");
					else bw.write(code[scores[s][t]]);
				}
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
