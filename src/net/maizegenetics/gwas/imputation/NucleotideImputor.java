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
import net.maizegenetics.baseplugins.ConvertSBitTBitPlugin;

import net.maizegenetics.baseplugins.TreeDisplayPlugin;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.ExportUtils;
import net.maizegenetics.pal.alignment.FilterAlignment;
import net.maizegenetics.pal.alignment.ImportUtils;
import net.maizegenetics.pal.alignment.Locus;
import net.maizegenetics.pal.alignment.MutableNucleotideAlignment;
import net.maizegenetics.pal.alignment.MutableSingleEncodeAlignment;
import net.maizegenetics.pal.alignment.NucleotideAlignmentConstants;
import net.maizegenetics.pal.alignment.BitAlignment;
import net.maizegenetics.pal.distance.DistanceMatrix;
import net.maizegenetics.pal.distance.IBSDistanceMatrix;
import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.Identifier;
import net.maizegenetics.pal.ids.SimpleIdGroup;
import net.maizegenetics.pal.tree.Tree;
import net.maizegenetics.pal.tree.TreeClusters;
import net.maizegenetics.pal.tree.UPGMATree;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.BitUtil;
import net.maizegenetics.util.OpenBitSet;

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
	Alignment gbsSnps;
	Alignment gbsSnpsByTaxa;
	ArrayList<PopulationData> familyList;
	String baseOutFilename;
	String assembly = "NA";
	int minAlleleCount = 2;
	double expectedProbHet = 0.0625;
	int maxcaCount = 0;
	double minrForLD = 0.4;
	double cutHeightForSnpClusters = 0.3;
	int windowSize = 100; //window size for calling snps
	
	static final byte AA = NucleotideAlignmentConstants.getNucleotideDiploidByte("AA");
	static final byte CC = NucleotideAlignmentConstants.getNucleotideDiploidByte("CC");
	static final byte GG = NucleotideAlignmentConstants.getNucleotideDiploidByte("GG");
	static final byte TT = NucleotideAlignmentConstants.getNucleotideDiploidByte("TT");
	static final byte AC = NucleotideAlignmentConstants.getNucleotideDiploidByte("AC");
	static final byte AG = NucleotideAlignmentConstants.getNucleotideDiploidByte("AG");
	static final byte AT = NucleotideAlignmentConstants.getNucleotideDiploidByte("AT");
	static final byte CG = NucleotideAlignmentConstants.getNucleotideDiploidByte("CG");
	static final byte CT = NucleotideAlignmentConstants.getNucleotideDiploidByte("CT");
	static final byte GT = NucleotideAlignmentConstants.getNucleotideDiploidByte("GT");
	static final byte NN = NucleotideAlignmentConstants.getNucleotideDiploidByte("NN");
	static final byte CA = NucleotideAlignmentConstants.getNucleotideDiploidByte("CA");
	
	static final byte[] byteval = new byte[] {AA,CC,GG,TT,AC};
	final static HashMap<Byte, Integer> genotypeMap = new HashMap<Byte, Integer>();
	
	static{
		genotypeMap.put(AA, 0);
		genotypeMap.put(CC, 1);
		genotypeMap.put(GG, 2);
		genotypeMap.put(TT, 3);
		genotypeMap.put(AC, 4);
		genotypeMap.put(CA, 4);
	}
	
	static final byte[][] genoval = new byte[][]{{AA,AC,AG,AT},{AC,CC,CG,CT},{AG,CG,GG,GT},{AT,CT,GT,TT}};
	
//	class Population {
//		PopulationData data = new PopulationData();
//	}
	
	public static void main(String[] args) {
		NucleotideImputor ni = new NucleotideImputor();
//		ni.processMaizeChromosome(1);
		ni.processCimmyt("1");
//		ni.imputeParentCalls();
//		ni.extractDataForAFamily();
	}
	
	public static void processRiceChromosomes() {
			NucleotideImputor ni = new NucleotideImputor();
			ni.importRiceChromosome(1);
	}
	
	public void importRiceChromosome(int chr) {
		String in = "/Volumes/Macintosh HD 2/data/chen/MAB_hmp_all_chr/MAB1_ct5.chr" + chr + ".full.hmp.txt";
		String out = "/Volumes/Macintosh HD 2/data/chen/MAB_hmp_all_chr/MAB1_ct5.chr" + chr + ".HMMimputed.hmp.txt";
		String filled = "/Volumes/Macintosh HD 2/data/chen/MAB_hmp_all_chr/MAB1_ct5.chr" + chr + ".HMMimputed.filled.hmp.txt";
		String ped = "/Volumes/Macintosh HD 2/data/chen/MAB_hmp_all_chr/MAB1_pedigree.txt";

		System.out.println("Imputing nucleotides from " + in);
		try {
			BufferedReader br = new BufferedReader(new FileReader(in));
			br.readLine();
			String[] info = br.readLine().split("\t", 12);
			assembly = info[5];
			br.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		gbsSnps = ImportUtils.readFromHapmap(in, true, null);
		
		System.out.println("Importing population information...");
		familyList = PopulationData.readPedigreeFile(ped);

		System.out.println("Scoring parents...");
		MutableNucleotideAlignment result = scoreParentsForAPopulation(familyList.get(0));
		ExportUtils.writeToHapmap(result, true, out, '\t', null);

		//set missing values to flanking values, when flanking markers are identical
		int ntaxa = result.getSequenceCount();
		int nsites = result.getSiteCount();
		int prevsite = -1;
		byte prevValue = -1;
		for (int t = 0; t < ntaxa; t++) {
			for (int s = 0; s < nsites; s++) {
				byte val = result.getBase(t, s);
				if (val != NN) {
					if (prevsite == -1) {
						prevsite = s;
						prevValue = val;
					} else if(val == prevValue) {
						for (int site = prevsite + 1; site < s; site++) {
							result.setBase(t, site, prevValue);
							prevsite = s;
						}
					} else {
						prevsite = s;
						prevValue = val;
					}
				}
			}
		}

		result.clean();
		ExportUtils.writeToHapmap(result, true, filled, '\t', null);

	}

	public void processMaizeChromosome(int chr) {
		String inHapmap = "/Volumes/Macintosh HD 2/data/namgbs/genos_20120110/merged_nam_ibm/NAM_IBM_282_Zea20120110_scv10mF8maf002_mgs_E1pLD5kpUn.c" + chr + ".hmp.txt";
		String outHapmap = "/Volumes/Macintosh HD 2/data/namgbs/genos_20120110/merged_nam_ibm/NAM_IBM_282_Zea20120110_scv10mF8maf002_mgs_E1pLD5kpUn.c" + chr + ".hmpImputed.txt";
//		String outFilledAC = "";
//		String outImputedLinkage = "";
		String ped = "/Volumes/Macintosh HD 2/data/namgbs/genos_20120110/merged_nam_ibm/namibm.pedigree.info.txt";

		gbsSnps = ImportUtils.readFromHapmap(inHapmap, true, null);
		
		System.out.println("Importing population information...");
		familyList = PopulationData.readPedigreeFile(ped);

		System.out.println("Scoring parents...");
		int familyCount = 0;
		for (PopulationData popdata : familyList) {
			System.out.println("Processing family " + popdata.name + "in NucleotideImputor.processMaizeChromosome()");
			String[] ids = new String[popdata.members.size()];
			popdata.members.toArray(ids);
			Alignment a =  FilterAlignment.getInstance(gbsSnps, new SimpleIdGroup(ids), false);
//			imputeSnpsForPopulation(a, pop);
			Alignment tba = callParentAlleles(a, popdata);
			String parentout = "/Volumes/Macintosh HD 2/data/namgbs/genos_20120110/merged_nam_ibm/NAM_IBM_282_Zea20120110_scv10mF8maf002_mgs_E1pLD5kpUn.c" + chr + "." + popdata.name + "parent.calls.hmpImputed.txt";
			ExportUtils.writeToHapmap(tba, false, parentout, '\t', null);
			
//			String out = "/Volumes/Macintosh HD 2/data/namgbs/genos_20120110/merged_nam_ibm/NAM_IBM_282_Zea20120110_scv10mF8maf002_mgs_E1pLD5kpUn.c" + chr + "." + family + ".hmpImputed.txt";
//			ExportUtils.writeToHapmap(pop.align, false, out, '\t', null);
		}
		
		//update the original alignment
		MutableNucleotideAlignment imputedGbsSnps = MutableNucleotideAlignment.getInstance(gbsSnps);
		for (PopulationData pop : familyList) {
			//TODO finish this part
		}
		
		String out = "";
		ExportUtils.writeToHapmap(imputedGbsSnps, false, out, '\t', null);
	}
	
	public void processCimmyt(String chr) {
//		String inHapmap = "/Volumes/Macintosh HD 2/data/cimmyt/cimmyt_imp_20120220/cimmyt_chr" + chr + "_imp_20120220.hmp.txt";
		String inHapmap = "/Volumes/Macintosh HD 2/data/cimmyt/cimmyt_unimp_20120220/cimmyt_chr" + chr + "_20120220.hmp.txt";
		String ped = "/Volumes/Macintosh HD 2/data/cimmyt/CIMMYT.3pop.peds.txt";
		gbsSnps = ImportUtils.readFromHapmap(inHapmap, true, null);
		
		familyList = PopulationData.readPedigreeFile(ped);

		for (PopulationData popdata:familyList) {
			System.out.println("Processing family " + popdata.name + "in NucleotideImputor.processMaizeChromosome()");
			String[] ids = new String[popdata.members.size()];
			popdata.members.toArray(ids);
			Alignment a =  FilterAlignment.getInstance(gbsSnps, new SimpleIdGroup(ids), false);
			Alignment parentCalls = callParentAlleles(a, popdata);
//			TBitAlignment parentCalls = callParentAllelesFromParents(a, pop);
			System.out.println("Parents called on " + parentCalls.getSiteCount() + " sites for family " + popdata.name);
//			String out = "/Volumes/Macintosh HD 2/data/cimmyt/cimmyt_unimp_20120220/cimmyt_chr1_20120220.parentsfromparents." + family + ".hmp.txt";
//			ExportUtils.writeToHapmap(parentCalls, false, out, '\t', null);
			Alignment va = imputeUsingViterbiFiveState(parentCalls);
			
			String out = "/Volumes/Macintosh HD 2/data/cimmyt/cimmyt_unimp_20120220/cimmyt_chr" + chr + "_20120220.imputed.family." + popdata.name + ".hmp.txt";
			ExportUtils.writeToHapmap(va, false, out, '\t', null);
		}

	}
	
	public void imputeParentCalls() {
		System.out.println("Imputing calls from parental calls...");
		String[] inFiles = new String[]{
			"/Volumes/Macintosh HD 2/data/cimmyt/cimmyt_unimp_20120220/cimmyt_chr1_20120220.parents.1.hmp.txt",
			"/Volumes/Macintosh HD 2/data/cimmyt/cimmyt_unimp_20120220/cimmyt_chr1_20120220.parents.2.hmp.txt",
			"/Volumes/Macintosh HD 2/data/cimmyt/cimmyt_unimp_20120220/cimmyt_chr1_20120220.parents.3.hmp.txt"
		};
		String[] outFiles = new String[]{
				"/Volumes/Macintosh HD 2/data/cimmyt/cimmyt_unimp_20120220/cimmyt_chr1_20120220.imputed.1.hmp.txt",
				"/Volumes/Macintosh HD 2/data/cimmyt/cimmyt_unimp_20120220/cimmyt_chr1_20120220.imputed.2.hmp.txt",
				"/Volumes/Macintosh HD 2/data/cimmyt/cimmyt_unimp_20120220/cimmyt_chr1_20120220.imputed.3.hmp.txt"
			};
		
		int n = inFiles.length;
		for (int i = 0; i < n; i++) {
			Alignment tba = ImportUtils.readFromHapmap(inFiles[i], false, null);
			//TBitAlignment tba = TBitAlignment.getInstance(a);
			ExportUtils.writeToHapmap(imputeUsingViterbiFiveState(tba), false, outFiles[i], '\t', null);
		}
		System.out.println("Finished.");
	}
	
	
	public void processWheatSomething() {
		String inHapmap = "/Volumes/Macintosh HD 2/data/wheat/wheat.f2.genotypes.txt";
		String ped = "/Volumes/Macintosh HD 2/data/wheat/wheat.ped.txt";
		gbsSnps = ImportUtils.readFromHapmap(inHapmap, true, null);
		
		familyList = PopulationData.readPedigreeFile(ped);
		
		PopulationData popdata = familyList.get(0);
		String[] ids = new String[popdata.members.size()];
		popdata.members.toArray(ids);
		Alignment a =  FilterAlignment.getInstance(gbsSnps, new SimpleIdGroup(ids), false);

		Locus[] theLoci = a.getLoci();
		int[] offsets = a.getLociOffsets();
		int nsites = a.getSiteCount();
		int nloci = theLoci.length;
		Alignment[] theAlignments = new Alignment[nloci];
		
		for (int locusnumber = 0; locusnumber < nloci; locusnumber++) {
			int start = offsets[locusnumber];
			int end = nsites - 1;
			if (locusnumber < nloci - 1) end = offsets[locusnumber + 1];
			
			Alignment locusAlignment = FilterAlignment.getInstance(a, start, end);
			theAlignments[locusnumber] = callParentAlleles(locusAlignment, popdata);
			System.out.println("Finished calling parental alleles for locus " + theLoci[locusnumber].getName());
		}
		
		Alignment combinedAlignments = MutableSingleEncodeAlignment.getInstance(theAlignments);
		String out = "/Volumes/Macintosh HD 2/data/wheat/family1.parent.encoding.hmp.txt";
		ExportUtils.writeToHapmap(combinedAlignments, false, out, '\t', null);
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
	 * Strategy for a single population
	 * 1. find parental haplotypes within populations
	 * 	a. remove monomorphic snps
	 * 	a. choose a window
	 * 	a. create an alignment for this window and these taxa
	 * 	b. within the window, cluster snps (UPGMA)
	 * 	c. cut the tree at 0.3
	 * 	d. use snps in largest cluster to define parental haplotypes
	 *  e. cluster on taxa to find parental haplotypes 
	 * 	e. extend haplotypes to ends of chromosome
	 * 2. assign parental haplotypes to specific parents
	 * 3. score progeny
	 * 
	 * */
	public MutableNucleotideAlignment scoreParentsForAPopulation(PopulationData popdata) {
		String[] ids = new String[popdata.members.size()];
		popdata.members.toArray(ids);
		Alignment popAlignment = FilterAlignment.getInstance(gbsSnps, new SimpleIdGroup(ids), false);
		
		//which sites are polymorphic? minor allele count > 2 and exceed the minimum allele count
		int nsites = popAlignment.getSiteCount();
		int ntaxa = popAlignment.getSequenceCount();
		OpenBitSet polybits = new OpenBitSet(nsites);
		for (int s = 0; s < nsites; s++) {
			int[][] freq = popAlignment.getAllelesSortedByFrequency(s);
			if (freq[1].length > 1 && freq[1][1] > 2) {
				int alleleCount = freq[1][0] + freq[1][1];
				if (alleleCount >= minAlleleCount) polybits.fastSet(s);
			}
		}
		int totalpoly = (int) polybits.cardinality();
		System.out.println("polymorphic sites = " + totalpoly + ", total sites = " + popAlignment.getSiteCount());
		
		//define parameters
		int windowSize = 100;
		FilterAlignment filteredPopAlignment;
		UPGMATree myTree;
		TreeClusters clusterMaker;
		int[] ldSnps = null;
		int[] snpIds = null;
		OpenBitSet ldbits = null;
		
		int ntrials = 20;
		for (int trial = 0; trial < ntrials; trial++) {
			//find windowSize polymorphic snps centered on the midpoint
			int mid = totalpoly / (ntrials + 1) * (trial + 1);
			int polyCount = 0;
			int start = 0;
			while (polyCount < mid) { //find the midpoint of the polymorphic snps
				if (polybits.get(start++)) polyCount++;
			}

			polyCount = 0;
			while (polyCount < (windowSize/2)) { //find the starting point of the window
				if (polybits.get(start--)) polyCount++;
			}
			start++;

			snpIds = new int[windowSize];
			polyCount = 0;
			int snpid = start;
			while (polyCount < windowSize && snpid < nsites) {
				if (polybits.get(snpid)) {
					snpIds[polyCount++] = snpid; 
				}
				snpid++;
			}
			
			if (polyCount < windowSize) snpIds = Arrays.copyOf(snpIds, polyCount);

			//create a filtered alignment containing only the test snps
			filteredPopAlignment = FilterAlignment.getInstance(popAlignment, snpIds);

			//cluster polymorphic snps within the window by creating a UPGMA tree (cluster on snps)
			Alignment haplotypeAlignment = BitAlignment.getInstance(filteredPopAlignment, true);
			myTree = new UPGMATree(snpDistance(haplotypeAlignment));

			//display tree for debugging
			TreeDisplayPlugin tdp = new TreeDisplayPlugin(null, true);
			tdp.performFunction(new DataSet(new Datum("Snp Tree", myTree, "Snp Tree"), null));

			//cut the tree
			clusterMaker = new TreeClusters(myTree);
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
				if (groupCount[i] > 5) System.out.println("For trial " + trial + " SNP group " + i + " has " + groupCount[i] + " members.");
			}

			//sort the groups
			int[] order = ImputationUtils.reverseOrder(groupCount);
			
			
			
			//this group of Snps is the starting point
			ldbits = new OpenBitSet(nsites);
			int selectedGroup = order[0];
			int ngrp = groupCount[selectedGroup];
			ldSnps = new int[ngrp];
			int grpCount = 0;
			for (int i = 0; i < groups.length; i++) {
				if (groups[i] == selectedGroup) {
					//the tree nodes are not in the same order as the alignment, so have to do the conversion
					int snpnumber = Integer.parseInt(myTree.getIdentifier(i).getFullName());
					ldSnps[grpCount++] = snpIds[snpnumber];
					ldbits.fastSet(snpIds[snpnumber]);
				}
			}
		}
		
		//cluster taxa for these snps to find parental haplotypes (cluster on taxa)
		filteredPopAlignment = FilterAlignment.getInstance(popAlignment, ldSnps);
		IBSDistanceMatrix dm = new IBSDistanceMatrix(BitAlignment.getInstance(filteredPopAlignment, true));
		estimateMissingDistances(dm);
		myTree = new UPGMATree(dm);
		clusterMaker = new TreeClusters(myTree);
		
		int majorCount = ntaxa;
		int minorCount = 0;
		int ngroups = 1;
		int[] groups = new int[0];
		int[] groupCount = null;
		int majorGroup = 0;
		int minorGroup = 1;
		
		while (majorCount > ntaxa / 2 && minorCount < 10) {
			ngroups++;
			groups = clusterMaker.getGroups(ngroups);
			groupCount = new int[ngroups];
			for (int gr : groups) groupCount[gr]++;
			
			for (int i = 1; i < ngroups; i++) {
				if (groupCount[i] > groupCount[majorGroup]) {
					minorGroup = majorGroup;
					majorGroup = i;
				} else if (groupCount[i] > groupCount[minorGroup]) minorGroup = i;
			}
			
			majorCount = groupCount[majorGroup];
			minorCount = groupCount[minorGroup];
		}
		
		//parents should be in different clusters, there should only be two major clusters
		
		 //List groups
		for (int i = 0; i < ngroups; i++) {
			if (groupCount[i] > 5) System.out.println("Taxa group " + i + " has " + groupCount[i] + " members.");
		}
		
		//display tree for debugging
		TreeDisplayPlugin tdp = new TreeDisplayPlugin(null, true);
		tdp.performFunction(new DataSet(new Datum("Snp Tree", myTree, "Snp Tree"), null));

		//use haplotypes to score parental type
		
		//create major and minor id groups
		String[] majorids = new String[groupCount[majorGroup]];
		String[] minorids = new String[groupCount[minorGroup]];
		majorCount = 0;
		minorCount = 0;
		for (int i = 0; i < groups.length; i++) {
			if (groups[i] == majorGroup) majorids[majorCount++] = myTree.getIdentifier(i).getFullName();
			else if (groups[i] == minorGroup) minorids[minorCount++] = myTree.getIdentifier(i).getFullName();
		}
		IdGroup majorTaxa = new SimpleIdGroup(majorids);
		IdGroup minorTaxa =  new SimpleIdGroup(minorids);
		
		//create an alignment for each cluster
		if (majorTaxa.whichIdNumber(popdata.parent1) > -1) {
			//nothing to do
		} else if (minorTaxa.whichIdNumber(popdata.parent1) > -1) {
			//swap groups
			IdGroup temp = majorTaxa;
			majorTaxa = minorTaxa;
			minorTaxa = temp;
		} else if(majorTaxa.whichIdNumber(popdata.parent2) > -1) {	
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
		Alignment sbitPopAlignment = BitAlignment.getInstance(popAlignment, true);
		MutableNucleotideAlignment parentAlignment = MutableNucleotideAlignment.getInstance(sbitPopAlignment);
		System.out.println("snps in parent Alignment = " + parentAlignment.getSiteCount());
		ntaxa = parentAlignment.getSequenceCount();
		
		for (int i = 0; i < ldSnps.length; i++) {
			int snp = ldSnps[i];
			byte parentA = aAlignment.getMajorAllele(i);
			byte parentC = cAlignment.getMajorAllele(i);
			
			for (int t = 0; t < ntaxa; t++) {
				byte[] taxon = popAlignment.getBaseArray(t, snp);
				if (taxon[0] == taxon[1]) {
					if (taxon[0] == parentA) parentAlignment.setBase(t, snp, AA);
					else if (taxon[0] == parentC) parentAlignment.setBase(t, snp, CC);
					else parentAlignment.setBase(t, snp, NN);
				} else if (taxon[0] == parentA) {
					if (taxon[1] == parentC) parentAlignment.setBase(t, snp, AC);
					else parentAlignment.setBase(t, snp, NN);
				} else if (taxon[0] == parentC) {
					if (taxon[1] == parentA) parentAlignment.setBase(t, snp, AC);
					else parentAlignment.setBase(t, snp, NN);
				} else {
					parentAlignment.setBase(t, snp, NN);
				}
			}
		}
		
		
		//extend haplotypes
		//add snps in ld
		int testSize = 25;
		double minr = 0.4;
		
		//add snps from middle to start; test only polymorphic snps
		LinkedList<Integer> testSnps = new LinkedList<Integer>();
		for (int i = testSize - 1; i >= 0; i--) testSnps.add(snpIds[i]);
		for (int snp = snpIds[0] - 1; snp >= 0; snp--) {
			if (polybits.fastGet(snp)) {
				byte[] ac = recodeParentalSnps(snp, testSnps, parentAlignment, minr);
				if (ac != null) {
					ldbits.fastSet(snp);
					testSnps.add(snp);
					testSnps.remove();
				}
			}
		}
		
		//add snps from middle to end
		testSnps.clear();
		int n = snpIds.length;
		for (int i = n - testSize; i < n; i++) testSnps.add(snpIds[i]);
		for (int snp = snpIds[n - 1] + 1; snp < nsites; snp++) {
			if (polybits.fastGet(snp)) {
				byte[] ac = recodeParentalSnps(snp, testSnps, parentAlignment, minr);
				if (ac != null) {
					ldbits.fastSet(snp);
					testSnps.add(snp);
					testSnps.remove();
				}
			}
		}
		
		parentAlignment.clean();
		n = (int) ldbits.size();
		int[] retainedSites = new int[(int) ldbits.cardinality()];
		int snpcount = 0;
		nsites = parentAlignment.getSiteCount();
		for (int i = 0; i < n; i++) {
			if (ldbits.fastGet(i)) retainedSites[snpcount++] = i;
		}
		
		FilterAlignment ldAlignment = FilterAlignment.getInstance(parentAlignment, retainedSites);
		System.out.println("Number of sites in the ldAlignment = " + ldAlignment.getSiteCount());
		Alignment tba = BitAlignment.getInstance(ldAlignment, false);
		
		System.out.println("Starting Viterbi algorithm...");
		
		MutableNucleotideAlignment a = imputeUsingViterbiFiveState(tba);
		
		return a;
	}
	
	public void imputeSnpsForPopulation(Alignment input, PopulationData popdata) {
		BitSet polybits = whichSitesArePolymorphic(input);
		int[] coreSnps = findCoreSnps(input, polybits, 5);
		String parentA = popdata.parent1;
		String parentC = popdata.parent2;
		
		OpenBitSet ldbits = new OpenBitSet(input.getSiteCount());
		int n = coreSnps.length;
		for (int i = 0; i < n; i++) ldbits.fastSet(coreSnps[i]);
		
		IdGroup[] taxaGroup =  findTaxaGroups(input, coreSnps);
		
		//create an alignment for each cluster
		IdGroup parentAGroup;
		IdGroup parentCGroup;
		if (taxaGroup[0].whichIdNumber(parentA) > -1) {
			parentAGroup = taxaGroup[0];
			parentCGroup = taxaGroup[1];
		} else if (taxaGroup[1].whichIdNumber(parentA) > -1) {
			parentAGroup = taxaGroup[1];
			parentCGroup = taxaGroup[0];
		} else if(taxaGroup[0].whichIdNumber(parentC) > -1) {	
			parentAGroup = taxaGroup[1];
			parentCGroup = taxaGroup[0];
		} else {
			parentAGroup = taxaGroup[0];
			parentCGroup = taxaGroup[1];
		} 
		
		Alignment aAlignment = FilterAlignment.getInstance(input, parentAGroup);
		Alignment cAlignment = FilterAlignment.getInstance(input, parentCGroup);
		byte[] Asnp = new byte[input.getSiteCount()];
		byte[] Csnp = new byte[input.getSiteCount()];
				
		//set first parent to AA, second parent to CC for snps used to form taxa clusters
		Alignment sbitPopAlignment = BitAlignment.getInstance(input, true);
		MutableNucleotideAlignment parentAlignment = MutableNucleotideAlignment.getInstance(sbitPopAlignment);
		System.out.println("snps in parent Alignment = " + parentAlignment.getSiteCount());
		int ntaxa = parentAlignment.getSequenceCount();
		
		for (int i = 0; i < coreSnps.length; i++) {
			int snp = coreSnps[i];
			byte alleleA = aAlignment.getMajorAllele(snp);
			byte alleleC = cAlignment.getMajorAllele(snp);
			Asnp[snp] = alleleA;
			Csnp[snp] = alleleC;
			 
			for (int t = 0; t < ntaxa; t++) {
				byte[] taxon = input.getBaseArray(t, snp);
				if (taxon[0] == taxon[1]) {
					if (taxon[0] == alleleA) parentAlignment.setBase(t, snp, AA);
					else if (taxon[0] == alleleC) parentAlignment.setBase(t, snp, CC);
					else parentAlignment.setBase(t, snp, NN);
				} else if (taxon[0] == alleleA) {
					if (taxon[1] == alleleC) parentAlignment.setBase(t, snp, AC);
					else parentAlignment.setBase(t, snp, NN);
				} else if (taxon[0] == alleleC) {
					if (taxon[1] == alleleA) parentAlignment.setBase(t, snp, AC);
					else parentAlignment.setBase(t, snp, NN);
				} else {
					parentAlignment.setBase(t, snp, NN);
				}
			}
		}

		//extend haplotypes
		//add snps in ld
		int testSize = 25;
		double minr = minrForLD;
		
		//add snps from middle to start; test only polymorphic snps
		LinkedList<Integer> testSnps = new LinkedList<Integer>();
		for (int i = testSize - 1; i >= 0; i--) testSnps.add(coreSnps[i]);
		for (int snp = coreSnps[0] - 1; snp >= 0; snp--) {
			if (polybits.fastGet(snp)) {
				byte[] ac = recodeParentalSnps(snp, testSnps, parentAlignment, minr);
				if (ac != null) {
					ldbits.fastSet(snp);
					testSnps.add(snp);
					testSnps.remove();
					Asnp[snp] = ac[0];
					Csnp[snp] = ac[1];
				}
			}
		}
		
		//add snps from middle to end
		testSnps.clear();
		n = coreSnps.length;
		int nsites = parentAlignment.getSiteCount();
		for (int i = n - testSize; i < n; i++) testSnps.add(coreSnps[i]);
		for (int snp = coreSnps[n - 1] + 1; snp < nsites; snp++) {
			if (polybits.fastGet(snp)) {
				byte[] ac = recodeParentalSnps(snp, testSnps, parentAlignment, minr);
				if (ac != null) {
					ldbits.fastSet(snp);
					testSnps.add(snp);
					testSnps.remove();
					Asnp[snp] = ac[0];
					Csnp[snp] = ac[1];
				}
			}
		}
		
		parentAlignment.clean();
		n = (int) ldbits.size();
		int nRetained = (int) ldbits.cardinality();
		popdata.alleleA = new byte[nRetained];
		popdata.alleleC = new byte[nRetained];
		int[] retainedSites = new int[nRetained];
		int snpcount = 0;
		for (int i = 0; i < n; i++) {
			if (ldbits.fastGet(i)) {
				popdata.alleleA[snpcount] = Asnp[i];
				popdata.alleleC[snpcount] = Csnp[i];
				retainedSites[snpcount++] = i;
			}
		}
		
		FilterAlignment ldAlignment = FilterAlignment.getInstance(parentAlignment, retainedSites);
		popdata.snpIndex = ldbits;
		
		Alignment tba = BitAlignment.getInstance(ldAlignment, false);
		System.out.println("number of original sites = " + input.getSiteCount() + ", number of polymorphic sites = " + polybits.cardinality() + ", number of ld sites = " + tba.getSiteCount());
		
		//debug export parent alignment
//		ExportUtils.writeToHapmap(tba, false, "/Volumes/Macintosh HD 2/temp/pop14.txt", '\t', null);
		
		System.out.println("Starting Viterbi algorithm...");
		popdata.imputed = imputeUsingViterbiFiveState(tba);
		 
		System.out.println("Finished.");
	}
	
	public Alignment callParentAlleles(Alignment input, PopulationData popdata) {
		BitSet polybits = whichSitesArePolymorphic(input);
		int[] coreSnps = findCoreSnps(input, polybits, 20);
		String parentA = popdata.parent1;
		String parentC = popdata.parent2;
		
		OpenBitSet ldbits = new OpenBitSet(input.getSiteCount());
		int n = coreSnps.length;
		for (int i = 0; i < n; i++) ldbits.fastSet(coreSnps[i]);
		
		IdGroup[] taxaGroup =  findTaxaGroups(input, coreSnps);
		
		//create an alignment for each cluster
		IdGroup parentAGroup;
		IdGroup parentCGroup;
		if (taxaGroup[0].whichIdNumber(parentA) > -1) {
			parentAGroup = taxaGroup[0];
			parentCGroup = taxaGroup[1];
		} else if (taxaGroup[1].whichIdNumber(parentA) > -1) {
			parentAGroup = taxaGroup[1];
			parentCGroup = taxaGroup[0];
		} else if(taxaGroup[0].whichIdNumber(parentC) > -1) {	
			parentAGroup = taxaGroup[1];
			parentCGroup = taxaGroup[0];
		} else {
			parentAGroup = taxaGroup[0];
			parentCGroup = taxaGroup[1];
		} 
		
		Alignment aAlignment = FilterAlignment.getInstance(input, parentAGroup);
		Alignment cAlignment = FilterAlignment.getInstance(input, parentCGroup);
		byte[] Asnp = new byte[input.getSiteCount()];
		byte[] Csnp = new byte[input.getSiteCount()];
				
		//set first parent to AA, second parent to CC for snps used to form taxa clusters
		Alignment sbitPopAlignment = BitAlignment.getInstance(input, true);
		MutableNucleotideAlignment parentAlignment = MutableNucleotideAlignment.getInstance(sbitPopAlignment);
		System.out.println("snps in parent Alignment = " + parentAlignment.getSiteCount());
		int ntaxa = parentAlignment.getSequenceCount();
		
		for (int i = 0; i < coreSnps.length; i++) {
			int snp = coreSnps[i];
			byte alleleA = aAlignment.getMajorAllele(snp);
			byte alleleC = cAlignment.getMajorAllele(snp);
			Asnp[snp] = alleleA;
			Csnp[snp] = alleleC;
			 
			for (int t = 0; t < ntaxa; t++) {
				byte[] taxon = input.getBaseArray(t, snp);
				if (taxon[0] == taxon[1]) {
					if (taxon[0] == alleleA) parentAlignment.setBase(t, snp, AA);
					else if (taxon[0] == alleleC) parentAlignment.setBase(t, snp, CC);
					else parentAlignment.setBase(t, snp, NN);
				} else if (taxon[0] == alleleA) {
					if (taxon[1] == alleleC) parentAlignment.setBase(t, snp, AC);
					else parentAlignment.setBase(t, snp, NN);
				} else if (taxon[0] == alleleC) {
					if (taxon[1] == alleleA) parentAlignment.setBase(t, snp, AC);
					else parentAlignment.setBase(t, snp, NN);
				} else {
					parentAlignment.setBase(t, snp, NN);
				}
			}
		}

		//extend haplotypes
		//add snps in ld
		int testSize = 25;
		double minr = minrForLD;
		
		//add snps from middle to start; test only polymorphic snps
		LinkedList<Integer> testSnps = new LinkedList<Integer>();
		for (int i = testSize - 1; i >= 0; i--) testSnps.add(coreSnps[i]);
		for (int snp = coreSnps[0] - 1; snp >= 0; snp--) {
			if (polybits.fastGet(snp)) {
				byte[] ac = recodeParentalSnps(snp, testSnps, parentAlignment, minr);
				if (ac != null) {
					ldbits.fastSet(snp);
					testSnps.add(snp);
					testSnps.remove();
					Asnp[snp] = ac[0];
					Csnp[snp] = ac[1];
				}
			}
		}
		
		//add snps from middle to end
		testSnps.clear();
		n = coreSnps.length;
		int nsites = parentAlignment.getSiteCount();
		for (int i = n - testSize; i < n; i++) testSnps.add(coreSnps[i]);
		for (int snp = coreSnps[n - 1] + 1; snp < nsites; snp++) {
			if (polybits.fastGet(snp)) {
				byte[] ac = recodeParentalSnps(snp, testSnps, parentAlignment, minr);
				if (ac != null) {
					ldbits.fastSet(snp);
					testSnps.add(snp);
					testSnps.remove();
					Asnp[snp] = ac[0];
					Csnp[snp] = ac[1];
				}
			}
		}
		
		parentAlignment.clean();
		n = (int) ldbits.size();
		int nRetained = (int) ldbits.cardinality();
		popdata.alleleA = new byte[nRetained];
		popdata.alleleC = new byte[nRetained];
		int[] retainedSites = new int[nRetained];
		int snpcount = 0;
		for (int i = 0; i < n; i++) {
			if (ldbits.fastGet(i)) {
				popdata.alleleA[snpcount] = Asnp[i];
				popdata.alleleC[snpcount] = Csnp[i];
				retainedSites[snpcount++] = i;
			}
		}
		
		FilterAlignment ldAlignment = FilterAlignment.getInstance(parentAlignment, retainedSites);
		popdata.snpIndex = ldbits;
		System.out.println("number of original sites = " + input.getSiteCount() + ", number of polymorphic sites = " + polybits.cardinality() + ", number of ld sites = " + ldAlignment.getSiteCount());

		return BitAlignment.getInstance(ldAlignment, false);
	}
	
	public BitSet whichSitesArePolymorphic(Alignment a) {
		//which sites are polymorphic? minor allele count > 2 and exceed the minimum allele count
		int nsites = a.getSiteCount();
		OpenBitSet polybits = new OpenBitSet(nsites);
		for (int s = 0; s < nsites; s++) {
			int[][] freq = a.getAllelesSortedByFrequency(s);
			if (freq[1].length > 1 && freq[1][1] > 2) {
				int alleleCount = freq[1][0] + freq[1][1];
				if (alleleCount >= minAlleleCount) polybits.fastSet(s);
			}
		}
		return polybits;
	}
	
	/**
	 * This function finds a set of snps within a window of the specified size (100) that are in LD with each other. It trys multiple windows and uses the
	 * window that yields the largest number of snps. 
	 * @param a	the input alignment
	 * @param polybits	a BitSet corresponding to SNPs in a, set if a snp is polymorphic. Only polymorphic SNPs will be considered.
	 * @param numberToTry the number of windows to try. The function will use the window returning the largest set of SNPs.
	 * @return indices of the core snps. 
	 */
	public int[] findCoreSnps(Alignment a, BitSet polybits, int numberToTry) {
		//define a window
		int totalpoly = (int) polybits.cardinality();
		int[][] snpSets = new int[numberToTry][];
		
		//find windowSize polymorphic snps centered on the midpoint
		int snpInterval = totalpoly / (numberToTry + 1);
		int start = - windowSize / 2;
		int snpCount = 0;
		int polyCount = 0;
		int nsites = a.getSiteCount();
		for (int setnum = 0; setnum < numberToTry; setnum++) {
			start += snpInterval;
			if (start < 0) start = 0;
			while (polyCount < start) {
				if (polybits.fastGet(snpCount)) polyCount++;
				snpCount++;
			}
			
			int[] snpIds = new int[windowSize];
			int windowCount = 0;
			while (windowCount < windowSize && snpCount < nsites) {
				if (polybits.fastGet(snpCount)) {
					snpIds[windowCount++] = snpCount;
					polyCount++;
				}
				snpCount++;
			}
			
			//adjust the size of the array if all the snps were used before the array was filled
			if (windowCount < windowSize) snpIds = Arrays.copyOf(snpIds, windowCount);
			
			//create a filtered alignment containing only the test snps
			FilterAlignment filteredPopAlignment = FilterAlignment.getInstance(a, snpIds);
			
			//cluster polymorphic snps within the window by creating a UPGMA tree (cluster on snps)
			Alignment haplotypeAlignment = BitAlignment.getInstance(filteredPopAlignment, true);
			UPGMATree myTree = new UPGMATree(snpDistance(haplotypeAlignment));
			
			//debug - display the tree 
//			TreeDisplayPlugin tdp = new TreeDisplayPlugin(null, true);
//			tdp.performFunction(new DataSet(new Datum("Snp Tree", myTree, "Snp Tree"), null));

			//cut the tree to create two parent groups
			TreeClusters clusterMaker = new TreeClusters(myTree);
			int[] groups = clusterMaker.getGroups(cutHeightForSnpClusters);
			
			//find the biggest group
			int maxGroup = 0;
			for (int grp:groups) maxGroup = Math.max(maxGroup, grp);
			int ngroups = maxGroup + 1;
			int[] groupCount = new int[ngroups];
			for (int grp:groups) groupCount[grp]++;
			int[]groupIndex = ImputationUtils.reverseOrder(groupCount);
			
			snpSets[setnum] = new int[ groupCount[groupIndex[0]] ];
			int count = 0;
			for (int i = 0; i < snpIds.length; i++) {
				if (groups[i] == groupIndex[0]) {
					int snpIndex = Integer.parseInt(myTree.getIdentifier(i).getFullName());
					snpSets[setnum][count++] = snpIds[snpIndex];
				}
			}
			Arrays.sort(snpSets[setnum]);
		}
		
		int bestSet = 0;
		for (int i = 1; i < numberToTry; i++) {
			if (snpSets[i].length > snpSets[bestSet].length) bestSet = i;
		}
		return snpSets[bestSet];
	}
	
	public IdGroup[] findTaxaGroups(Alignment a, int[] coreSnps) {
		
		//cluster taxa for these snps to find parental haplotypes (cluster on taxa)
		IBSDistanceMatrix dm = new IBSDistanceMatrix(BitAlignment.getInstance(FilterAlignment.getInstance(a, coreSnps), true));
		estimateMissingDistances(dm);
		Tree myTree = new UPGMATree(dm);
		TreeClusters clusterMaker = new TreeClusters(myTree);
		
		int ntaxa = a.getSequenceCount();
		int majorCount = ntaxa;
		int minorCount = 0;
		int ngroups = 1;
		int[] groups = new int[0];
		int[] groupCount = null;
		int majorGroup = 0;
		int minorGroup = 1;
		
		while (majorCount > ntaxa / 2 && minorCount < 10) {
			ngroups++;
			groups = clusterMaker.getGroups(ngroups);
			groupCount = new int[ngroups];
			for (int gr : groups) groupCount[gr]++;
			
			for (int i = 1; i < ngroups; i++) {
				if (groupCount[i] > groupCount[majorGroup]) {
					minorGroup = majorGroup;
					majorGroup = i;
				} else if (groupCount[i] > groupCount[minorGroup]) minorGroup = i;
			}
			
			majorCount = groupCount[majorGroup];
			minorCount = groupCount[minorGroup];
		}

		//debug - display the tree 
//		TreeDisplayPlugin tdp = new TreeDisplayPlugin(null, true);
//		tdp.performFunction(new DataSet(new Datum("Snp Tree", myTree, "Snp Tree"), null));

		 //List groups
		for (int i = 0; i < ngroups; i++) {
			if (groupCount[i] > 5) System.out.println("Taxa group " + i + " has " + groupCount[i] + " members.");
		}
		
		//create major and minor id groups
		String[] majorids = new String[groupCount[majorGroup]];
		String[] minorids = new String[groupCount[minorGroup]];
		majorCount = 0;
		minorCount = 0;
		for (int i = 0; i < groups.length; i++) {
			if (groups[i] == majorGroup) majorids[majorCount++] = myTree.getIdentifier(i).getFullName();
			else if (groups[i] == minorGroup) minorids[minorCount++] = myTree.getIdentifier(i).getFullName();
		}
		IdGroup majorTaxa = new SimpleIdGroup(majorids);
		IdGroup minorTaxa =  new SimpleIdGroup(minorids);
		return new IdGroup[]{majorTaxa,minorTaxa};
	}
	
	//returns a byte array containing the A allele as element 0 and the C allele as element 1, or null if the snp is not in LD
	private byte[] recodeParentalSnps(int snp, LinkedList<Integer> testSnps, MutableNucleotideAlignment snpAlignment, double minr) {
		int ntaxa = snpAlignment.getSequenceCount();
		byte[] snpvals = new byte[ntaxa];
		for (int t = 0; t < ntaxa; t++) {
			snpvals[t] = snpAlignment.getBase(t, snp);
		}
		
		int[] acount = new int[5];
		int[] ccount = new int[5];
		for (int t = 0; t < ntaxa; t++) {
			Integer ndx = genotypeMap.get(snpvals[t]);
			if (ndx != null) {
				int indx = ndx.intValue();
				for (Integer testsnp:testSnps) {
					byte testval = snpAlignment.getBase(t, testsnp);
					if (testval == AA) acount[ndx]++;
					else if (testval == CC) ccount[ndx]++;
				}
			}
		}
		
		//calculate r
		int maxa = 0;
		int maxc = 0;
		for (int i = 1; i < 4; i++) {
			if (acount[i] > acount[maxa]) maxa = i;
			if (ccount[i] > ccount[maxc]) maxc = i;
		}
		
		int[][] counts = new int[2][2];
		counts[0][0] = acount[maxa];
		counts[0][1] = acount[maxc];
		counts[1][0] = ccount[maxa];
		counts[1][1] = ccount[maxc];
		
		double r = Math.abs(calculateR(counts));
//		if (maxa == maxc) {
//			maxcaCount++;
//			System.out.print("");
//		}
		
		//if abs(r) > 0.8, recode the snp
		if ( maxa != maxc && r >= minr) {
			byte hetval = genoval[maxa][maxc];
			for (int t = 0; t < ntaxa; t++) {
				byte val = snpvals[t];
				if (val == byteval[maxa]) snpAlignment.setBase(t, snp, AA);
				else if (val == byteval[maxc]) snpAlignment.setBase(t, snp, CC);
				else if (val == hetval) snpAlignment.setBase(t, snp, AC);
				else snpAlignment.setBase(t, snp, NN);
			}
			return new byte[]{(byte) maxa, (byte) maxc};
		}
		
		return null;
		
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
	
	public MutableNucleotideAlignment imputeUsingViterbiFiveState(Alignment a) {
                a = ConvertSBitTBitPlugin.convertAlignment(a, ConvertSBitTBitPlugin.CONVERT_TYPE.tbit, null);
		//states are in {all A; 3A:1C; 1A:1C, 1A:3C; all C}
		//obs are in {A, C, M}, where M is heterozygote A/C
		int maxIterations = 50;
		HashMap<Byte, Byte> obsMap = new HashMap<Byte, Byte>();
		obsMap.put(AA, (byte) 0);
		obsMap.put(AC, (byte) 1);
		obsMap.put(CA, (byte) 1);
		obsMap.put(CC, (byte) 2);
		
		int ntaxa = a.getSequenceCount();
		int nsites = a.getSiteCount();
		
		//initialize the transition matrix
		double[][] transition = new double[][] {
				{.999,.0001,.0003,.0001,.0005},
				{.0002,.999,.00005,.00005,.0002},
				{.0002,.00005,.999,.00005,.0002},
				{.0002,.00005,.00005,.999,.0002},
				{.0005,.0001,.0003,.0001,.999}
		};
		
		TransitionProbability tp = new TransitionProbability();
		tp.setTransitionProbability(transition);
		int chrlength = a.getPositionInLocus(nsites - 1) - a.getPositionInLocus(0);
		tp.setAverageSegmentLength( chrlength / nsites );
		
		
		//initialize the emission matrix, states (5) in rows, observations (3) in columns
		double[][] emission = new double[][] {
				{.98,.001,.001},
				{.6,.2,.2},
				{.4,.2,.4},
				{.2,.2,.6},
				{.001,.001,.98}
		};
		
		EmissionProbability ep = new EmissionProbability();
		ep.setEmissionProbability(emission);
		
		//set up indices to non-missing data
		ArrayList<BitSet> notMissingIndex = new ArrayList<BitSet>();
		int[] notMissingCount = new int[ntaxa];
		ArrayList<byte[]> nonMissingObs = new ArrayList<byte[]>();
		ArrayList<int[]> snpPositions = new ArrayList<int[]>();
		
		for (int t = 0; t < ntaxa; t++) {
			long[] bits = a.getAllelePresenceForAllSites(t, 0).getBits();
			BitSet notMiss = new OpenBitSet(bits, bits.length);
			notMiss.or(a.getAllelePresenceForAllSites(t, 1));
			notMissingIndex.add(notMiss);
			notMissingCount[t] = (int) notMiss.cardinality();
		}
		
		for (int t = 0; t < ntaxa; t++) {
			byte[] obs = new byte[notMissingCount[t]];
			int[] pos = new int[notMissingCount[t]];
			nonMissingObs.add(obs);
			snpPositions.add(pos);
			BitSet isNotMissing = notMissingIndex.get(t);
			int nmcount = 0;
			for (int s = 0; s < nsites; s++) {
				byte base = a.getBase(t, s);
				if (isNotMissing.fastGet(s) && obsMap.get(a.getBase(t, s)) == null) {
					System.out.println("null from " + Byte.toString(base));
				}
				if (isNotMissing.fastGet(s)) {
					obs[nmcount] = obsMap.get(a.getBase(t, s));
					pos[nmcount++] = a.getPositionInLocus(s);
				}
				
			}
			
		}
		
		double phom = (1 - expectedProbHet) / 2;
		double[] pTrue = new double[]{phom, .25*expectedProbHet ,.5 * expectedProbHet, .25*expectedProbHet, phom};
		
		//iterate
		ArrayList<byte[]> bestStates = new ArrayList<byte[]>();
		int[][] previousStateCount = new int[5][3];
		int iter = 0;
		boolean hasNotConverged = true;
		while (iter < maxIterations && hasNotConverged) {
			//apply Viterbi
			System.out.println("Iteration " + iter++);
			bestStates.clear();
			for (int t = 0; t < ntaxa; t++) {
				tp.setPositions(snpPositions.get(t));
				int nobs = notMissingCount[t];
				if (nobs >= 20) {
					ViterbiAlgorithm va = new ViterbiAlgorithm(nonMissingObs.get(t), tp, ep, pTrue);
					va.calculate();
					bestStates.add(va.getMostProbableStateSequence());
				} else { //do not impute if obs < 20
					System.out.println("Fewer then 20 observations for " + a.getTaxaName(t));
					byte[] states = new byte[nobs];
					byte[] obs = nonMissingObs.get(t);
					for (int i = 0; i < nobs; i++) {
						if (obs[i] == AA) states[i] = 0;
						else if (obs[i] == CC) states[i] = 4;
						else states[i] = 2;
					}
					bestStates.add(states);
				}
			}
			
			//re-estimate transition probabilities
			int[][] transitionCounts = new int[5][5];
			double[][] transitionProb = new double[5][5];
			for (int t = 0; t < ntaxa; t++) {
				byte[] states = bestStates.get(t);
				for (int s = 1; s < notMissingCount[t]; s++) {
					transitionCounts[states[s-1]][states[s]]++;
				}
			}
			
			//transition is prob(state2 | state1) = count(cell)/count(row)
			for (int row = 0; row < 5; row++) {
				double rowsum = 0;
				for (int col = 0; col < 5; col++) rowsum += transitionCounts[row][col];
				for (int col = 0; col < 5; col++) transitionProb[row][col] = ((double) transitionCounts[row][col]) / rowsum;
			}
			tp.setTransitionCounts(transitionCounts, chrlength, ntaxa);
			
			//print transition probabilities
			System.out.println("Transition probabilities");
			for (double[] row:transitionProb) {
				for (double cell:row) {
					System.out.print(cell  + " ");
				}
				System.out.println();
			}
			System.out.println();
			
			//re-estimate emission probabilities
			int[][] emissionCounts = new int[5][3];
			double[][] emissionProb = new double[5][3];
			for (int t = 0; t < ntaxa; t++) {
				byte[] obs = nonMissingObs.get(t);
				byte[] states = bestStates.get(t);
				for (int s = 0; s < notMissingCount[t]; s++) {
					emissionCounts[states[s]][obs[s]]++;
				}
			}
			
			//print observation/state counts
			System.out.println("Imputation counts, rows=states, columns=observations");
			for (int[] row:emissionCounts) {
				for (int cell:row) {
					System.out.print(cell  + " ");
				}
				System.out.println();
			}
			System.out.println();
			
			//check to see if there is a change in the observation/state counts
			hasNotConverged = false;
			for (int r = 0; r < 5; r++) {
				for (int c = 0; c < 3; c++) {
					if (previousStateCount[r][c] != emissionCounts[r][c]) {
						hasNotConverged = true;
						previousStateCount[r][c] = emissionCounts[r][c];
					}
				}
			}
				
			
			//emission is prob(obs | state) = count(cell)/count(row)
			for (int row = 0; row < 5; row++) {
				double rowsum = 0;
				for (int col = 0; col < 3; col++) rowsum += emissionCounts[row][col];
				for (int col = 0; col < 3; col++) emissionProb[row][col] = ((double) emissionCounts[row][col]) / rowsum;
			}
			ep.setEmissionProbability(emissionProb);
			
			//print emission probabilities
			System.out.println("Emission probabilities");
			for (double[] row:emissionProb) {
				for (double cell:row) {
					System.out.print(cell  + " ");
				}
				System.out.println();
			}
			System.out.println();
		}
		
		MutableNucleotideAlignment result = MutableNucleotideAlignment.getInstance(a);
		nsites = result.getSiteCount();
		for (int t = 0; t < ntaxa; t++) {
			BitSet hasData = notMissingIndex.get(t);
			byte[] states = bestStates.get(t);
			int stateCount = 0;
			for (int s = 0; s < nsites; s++) {
				if (hasData.fastGet(s)) {
					if (states[stateCount] == 0) result.setBase(t, s, AA);
					else if (states[stateCount] < 4) result.setBase(t, s, AC);
					else if (states[stateCount] == 4) result.setBase(t, s, CC);
					stateCount++;
				}
			}
		}
		
		result.clean();
		
		return result;
	}
	
	public void writeScoresHapmap(Alignment a, String family) {
		String outfile = baseOutFilename + "_family" + family + ".hmp.txt";
		ExportUtils.writeToHapmap(a, true, baseOutFilename, '\t', null);
	}
	
	public void writeScoresNumeric(Alignment a, String family) {
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
					byte val = a.getBase(t, s);
					if (val == AA) bw.write("0");
					else if (val == CC) bw.write("2");
					else if (val == AC) bw.write("1");
					else bw.write("-1");
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
//				distance[i][j] = distance[j][i] = 1 - Math.abs(computeR(i, j, a));
				double r = computeR(i, j, a);
				distance[i][j] = distance[j][i] = 1 - r*r;
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
	
	public DistanceMatrix snpDistance2(Alignment a) {
                a = ConvertSBitTBitPlugin.convertAlignment(a, ConvertSBitTBitPlugin.CONVERT_TYPE.sbit, null);
		
		int nsnps = a.getSiteCount();
		SimpleIdGroup snpIds = new SimpleIdGroup(nsnps, true);
		double[][] distance = new double[nsnps][nsnps];
		double sum = 0;
		int count = 0;
		for (int i = 0; i < nsnps; i++) {
			long[][][] snps = new long[2][2][];
			snps[0][0] = a.getAllelePresenceForAllTaxa(i, 0).getBits();
			snps[0][1] = a.getAllelePresenceForAllTaxa(i, 1).getBits();
			int n = snps[0][0].length;
			distance[i][i] = 0;
			for (int j = i + 1; j < nsnps; j++) {
				snps[1][0] = a.getAllelePresenceForAllTaxa(j, 0).getBits();
				snps[1][1] = a.getAllelePresenceForAllTaxa(j, 1).getBits();
				long[] bothHom = new long[n];
				for (int m = 0; m < n; m++) bothHom[m] = (snps[0][0][m]^snps[0][1][m]) & (snps[1][0][m]^snps[1][1][m]);
				long[] same = new long[n];
				for (int m = 0; m < n; m++) same[m] = bothHom[m] & ~(snps[0][0][m]^snps[1][0][m]);
				distance[i][j] = distance[j][i] = 1 - ((double) BitUtil.pop_array(same, 0, n)) / ((double) BitUtil.pop_array(bothHom, 0, n));
			}
		}
		
		double avgDistance = 0;
		count = 0;
		for (double[] row : distance) {
			for (double cell : row) {
				if (!Double.isNaN(cell)) {
					avgDistance += cell;
					count++;
				}
			}
		}
		
		avgDistance /= count;
		for (double[] row : distance) {
			for (double cell : row) {
				if (!Double.isNaN(cell)) {
					cell = avgDistance;
				}
			}
		}
		
		return new DistanceMatrix(distance, new SimpleIdGroup(nsnps, true));
		
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

	public void updateSnpAlignment(MutableNucleotideAlignment mna, PopulationData popdata) {
		//set monomorphic sites to major (or only allele) (or not)
		//set polymorphic sites consistent with flanking markers if equal, unchanged otherwise
		//do not change sites that are not clearly monomorhpic or polymorphic
		
		fillGapsInAlignment(popdata);
		
		//map population taxa to mna
		int ntaxa = popdata.imputed.getSequenceCount();
		int[] taxaIds = new int[ntaxa];
		for (int t = 0; t < ntaxa; t++) {
			Identifier ident = popdata.imputed.getIdGroup().getIdentifier(t);
			if (popdata.parent1.equals(ident.getFullName()) || popdata.parent2.equals(ident.getFullName())) {
				taxaIds[t] = -1;
			} else {
				taxaIds[t] = mna.getIdGroup().whichIdNumber(ident);
			}
		}
		
		BitSet isPopSnp = popdata.snpIndex;
		
		int nsites = mna.getSiteCount();
		int popSnpCount = 0;
		for (int s = 0; s < nsites; s++) {
			if (isPopSnp.fastGet(s)) {
				int Acall = popdata.alleleA[popSnpCount];
				int Ccall = popdata.alleleC[popSnpCount];
				byte AAcall = (byte) ((Acall << 4) | Acall);
				byte CCcall = (byte) ((Ccall << 4) | Ccall);
				byte ACcall = (byte) ((Acall << 4) | Ccall);
				for (int t = 0; t < ntaxa; t++) {
					if (taxaIds[t] > -1) {
						byte parentCall = popdata.imputed.getBase(t, popSnpCount);
						if (parentCall == AA) {
							mna.setBase(taxaIds[t], s, AAcall);
						} else if (parentCall == CC) {
							mna.setBase(taxaIds[t], s, CCcall);
						} else if (parentCall == AC || parentCall == CA) {
							mna.setBase(taxaIds[t], s, ACcall);
						} else {
							mna.setBase(taxaIds[t], s, NN);
						}
					}
				}
				popSnpCount++;
			} else { //if the site is monomorphic fill in all the genotypes
				//do nothing for now.
			}
		}
		
	}
	
	public void fillGapsInAlignment(PopulationData popdata) {
		MutableNucleotideAlignment a;
		if (popdata.imputed instanceof MutableNucleotideAlignment) {
			a = (MutableNucleotideAlignment) popdata.imputed;
		} else {
			a = MutableNucleotideAlignment.getInstance(popdata.imputed);
			popdata.imputed = a;
		}
		
		int ntaxa = a.getSequenceCount();
		int nsites = a.getSiteCount();
		int prevsite = -1;
		byte prevValue = -1;
		for (int t = 0; t < ntaxa; t++) {
			for (int s = 0; s < nsites; s++) {
				byte val = a.getBase(t, s);
				if (val != NN) {
					if (prevsite == -1) {
						prevsite = s;
						prevValue = val;
					} else if(val == prevValue) {
						for (int site = prevsite + 1; site < s; site++) {
							a.setBase(t, site, prevValue);
							prevsite = s;
						}
					} else {
						prevsite = s;
						prevValue = val;
					}
				}
			}
		}

		a.clean();
				
	}
	
	public Alignment callParentAllelesFromParents(Alignment a, PopulationData popdata) {
		int minAlleleCount = 3;
		
		//which are the parents?
		int p1 = a.getIdGroup().whichIdNumber(popdata.parent1);
		int p2 = a.getIdGroup().whichIdNumber(popdata.parent2);
		
		if (p1 == -1 || p2 == -1) {
			if (p2 > -1) System.out.println("In callParentAllelesFromImputedParents, Parent 1 : " + popdata.parent1 + " has no data. Parent 2 does : " + popdata.parent2);
			else if (p1 > -1) System.out.println("In callParentAllelesFromImputedParents, Parent 2 : " + popdata.parent2 + " has no data. Parent 1 does : " + popdata.parent1);
			else System.out.println("In callParentAllelesFromImputedParents, Parent 1 : " + popdata.parent1 + " andParent 2 : " + popdata.parent2 + " have no data.");
			return null;
		}
		
		//convert a to an sbit alignment
		Alignment sba = ConvertSBitTBitPlugin.convertAlignment(a, ConvertSBitTBitPlugin.CONVERT_TYPE.sbit, null);
		//if (a instanceof SBitAlignment) {
		//	sba = (SBitAlignment) a;
		//} else {
		//	sba = SBitAlignment.getInstance(a);
		//}
		
		int nsites = sba.getSiteCount();
		OpenBitSet ispoly = new OpenBitSet(nsites);
		
		for (int s = 0; s < nsites; s++) {
			BitSet major = sba.getAllelePresenceForAllTaxa(s, 0);
			BitSet minor = sba.getAllelePresenceForAllTaxa(s, 1);
			
			//are parents different and not heterozygous?
			if ( (major.fastGet(p1)^minor.fastGet(p1)) && (major.fastGet(p1)^major.fastGet(p2)) && (minor.fastGet(p1)^minor.fastGet(p2)) ) {
				int majorCount = (int) OpenBitSet.andNotCount(major, minor);
				int minorCount = (int) OpenBitSet.andNotCount(minor, major);
				if (majorCount >= minAlleleCount && minorCount >= minAlleleCount) ispoly.fastSet(s);
			}
		}
		
		int npoly = (int) ispoly.cardinality();
		int[] snpids = new int[npoly];
		
		int polycount = 0;
		for (int s = 0; s < nsites; s++) {
			if (ispoly.fastGet(s)) snpids[polycount++] = s;
		}
		
		FilterAlignment fa = FilterAlignment.getInstance(sba, snpids);
		
		MutableNucleotideAlignment mna = MutableNucleotideAlignment.getInstance(BitAlignment.getInstance(fa, true));
		
		fa = null;
		sba = null;
		
		int ntaxa = mna.getSequenceCount();
		nsites = mna.getSiteCount();
		for (int s = 0; s < nsites; s++) {
			byte genop1 = mna.getBase(p1, s);
			byte genop2 = mna.getBase(p2, s);
			byte het = (byte) ((genop1 & 0x0F) | (genop2 & 0xF0));
			byte revhet = (byte) ((genop1 & 0xF0) | (genop2 & 0x0F));
			for (int t = 0; t < ntaxa; t++) {
				byte bval = mna.getBase(t, s);
				if (bval == genop1) mna.setBase(t, s, AA);
				else if (bval == genop2) mna.setBase(t, s, CC);
				else if (bval == het) mna.setBase(t, s, AC);
				else if (bval == revhet) mna.setBase(t, s, AC);
				else mna.setBase(t, s, NN);
			}
		}
		
		mna.clean();
		return BitAlignment.getInstance(mna, false);
	}
	
	public void extractDataForAFamily() {
		String genoFile = "/Volumes/Macintosh HD 2/data/cimmyt/cimmyt_imp_20120220/cimmyt_chr1_imp_20120220.hmp.txt";
		String pedFile = "/Volumes/Macintosh HD 2/data/cimmyt/CIMMYT.3pop.peds.txt";
		
		Alignment a = ImportUtils.readFromHapmap(genoFile, null);
		familyList = PopulationData.readPedigreeFile(pedFile);
		
		for (PopulationData family:familyList) {
			System.out.println("Extracting data for family " + family.name);
			int n = family.members.size();
			String[] ids = new String[n];
			family.members.toArray(ids);
			Alignment b = FilterAlignment.getInstance(a, new SimpleIdGroup(ids));
			String out = "/Volumes/Macintosh HD 2/data/cimmyt/cimmyt_imp_20120220/cimmyt_chr1_imp_20120220.family." + family.name + ".hmp.txt";
			ExportUtils.writeToHapmap(b, false, out, '\t', null);
		}
	}
}


