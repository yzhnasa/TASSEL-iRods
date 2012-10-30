package net.maizegenetics.gwas.imputation;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Random;
import java.util.Set;
import java.util.regex.Pattern;

import net.maizegenetics.baseplugins.ConvertSBitTBitPlugin;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.BitAlignment;
import net.maizegenetics.pal.alignment.FilterAlignment;
import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.IdGroupUtils;
import net.maizegenetics.util.BitSet;

public class ImputationUtils {
	public static Pattern tab = Pattern.compile("\t");
	
	public static int[] order(int[] array) {
		class SortElement implements Comparable<SortElement> {
			int val;
			int ndx;
			
			SortElement(int x, int index) {
				val = x;
				ndx = index;
			}
			
			@Override
			public int compareTo(SortElement se) {
				return val - se.val;
			}
		}

		int n = array.length;
		SortElement[] sortArray = new SortElement[n];
		for (int i = 0; i < n; i++) {
			sortArray[i] = new SortElement(array[i], i);
		}
		
		Arrays.sort(sortArray);
		int[] order = new int[n];
		for (int i = 0; i < n; i++) {
			order[i] = sortArray[i].ndx;
		}
		return order;
	}
	
	public static int[] reverseOrder(int[] array) {
		class SortElement implements Comparable<SortElement> {
			int val;
			int ndx;
			
			SortElement(int x, int index) {
				val = x;
				ndx = index;
			}
			
			@Override
			public int compareTo(SortElement se) {
				return se.val - val;
			}
		}

		int n = array.length;
		SortElement[] sortArray = new SortElement[n];
		for (int i = 0; i < n; i++) {
			sortArray[i] = new SortElement(array[i], i);
		}
		
		Arrays.sort(sortArray);
		int[] order = new int[n];
		for (int i = 0; i < n; i++) {
			order[i] = sortArray[i].ndx;
		}
		return order;
	}
	
	public static Alignment[] getTwoClusters(Alignment a, int[] parentIndex) {
		int maxiter = 5;
		Alignment tb = ConvertSBitTBitPlugin.convertAlignment(a, ConvertSBitTBitPlugin.CONVERT_TYPE.tbit, null);
		//if (a instanceof TBitAlignment) tb = (TBitAlignment) a;
		//else tb = TBitAlignment.getInstance(a);
		
		//if the parents are in the data set use these as seeds
		//if one parent is in the dataset pick the taxon farthest from it as the other seed
		//if neither parent is in the dataset choose random seeds
		int ntaxa = tb.getSequenceCount();
		int nsnps = tb.getSiteCount();
		int seed1 = parentIndex[0];
		int seed2 = parentIndex[1];
		float[] loc1, loc2;
		Random rand = new Random();
		if (seed1 == -1) {
			if (seed2 == -1) {
				//both parents are not in the data set
				seed1 = rand.nextInt(ntaxa);
				loc1 = snpsAsFloatVector(new BitSet[]{tb.getAllelePresenceForAllSites(seed1, 0), tb.getAllelePresenceForAllSites(seed1, 1)}, nsnps);
				while (seed2 == -1 || seed1 == seed2) {
					seed2 = rand.nextInt(ntaxa);
				}
				loc2 = snpsAsFloatVector(new BitSet[]{tb.getAllelePresenceForAllSites(seed2, 0), tb.getAllelePresenceForAllSites(seed2, 1)}, nsnps);
			} else {
				//parent2 is in the data set
				loc1 = snpsAsFloatVector(new BitSet[]{tb.getAllelePresenceForAllSites(0, 0), tb.getAllelePresenceForAllSites(0, 1)}, nsnps);
				loc2 = snpsAsFloatVector(new BitSet[]{tb.getAllelePresenceForAllSites(seed2, 0), tb.getAllelePresenceForAllSites(seed2, 1)}, nsnps);
				seed1 = 0;
				float prevdist = getManhattanDistance(loc2, loc1, nsnps);
				for (int t = 1; t < ntaxa; t++) {
					float[] tloc = snpsAsFloatVector(new BitSet[]{tb.getAllelePresenceForAllSites(t, 0), tb.getAllelePresenceForAllSites(t, 1)}, nsnps);
					float dist = getManhattanDistance(loc2, tloc, nsnps);
					if (dist > prevdist) {
						prevdist = dist;
						loc1 = tloc;
						seed1 = t;
					}
				}
			}
		} else if (seed2 == -1) {
			//parent1 is in the data set
			loc1 = snpsAsFloatVector(new BitSet[]{tb.getAllelePresenceForAllSites(seed1, 0), tb.getAllelePresenceForAllSites(seed1, 1)}, nsnps);
			loc2 = snpsAsFloatVector(new BitSet[]{tb.getAllelePresenceForAllSites(0, 0), tb.getAllelePresenceForAllSites(0, 1)}, nsnps);
			seed2 = 0;
			float prevdist = getManhattanDistance(loc1, loc2, nsnps);
			for (int t = 1; t < ntaxa; t++) {
				float[] tloc = snpsAsFloatVector(new BitSet[]{tb.getAllelePresenceForAllSites(t, 0), tb.getAllelePresenceForAllSites(t, 1)}, nsnps);
				float dist = getManhattanDistance(loc1, tloc, nsnps);
				if (dist > prevdist) {
					prevdist = dist;
					loc2 = tloc;
					seed2 = t;
				}
			}
		} else {
			//both parents are in the data set
			loc1 = snpsAsFloatVector(new BitSet[]{tb.getAllelePresenceForAllSites(seed1, 0), tb.getAllelePresenceForAllSites(seed1, 1)}, nsnps);
			loc2 = snpsAsFloatVector(new BitSet[]{tb.getAllelePresenceForAllSites(seed2, 0), tb.getAllelePresenceForAllSites(seed2, 1)}, nsnps);
		}
		
		int[] size1 = new int[nsnps];
		int[] size2 = new int[nsnps];
		for (int i = 0; i < nsnps; i++) {
			if (loc1[i] >= 0) size1[i] = 1;
			if (loc2[i] >= 0) size2[i] = 1;
		}
		boolean[] isInCluster1 = new boolean[ntaxa];
		isInCluster1[seed1] = true;
		isInCluster1[seed2] = false;
		
		
		//do initial cluster assignment
		for (int t = 0; t < ntaxa; t++) {
			if (t != seed1 && t != seed2) {
				float[] tloc = snpsAsFloatVector(new BitSet[]{tb.getAllelePresenceForAllSites(t, 0), tb.getAllelePresenceForAllSites(t, 1)}, nsnps);
				float dist1 = getManhattanDistance(loc1, tloc, nsnps);
				float dist2 = getManhattanDistance(loc2, tloc, nsnps);
				if (dist1 <= dist2) {
					isInCluster1[t] = true;
					loc1 = getMeanLocation(loc1, size1, tloc, true, nsnps);
				} else {
					isInCluster1[t] = false;
					loc2 = getMeanLocation(loc2, size2, tloc, true, nsnps);
				}
			}
		}
		
		//update cluster membership until there are no changes or for the maximum number of iterations
		for (int iter = 0; iter < maxiter; iter++) {
			boolean noChanges = true;
			for (int t = 0; t < ntaxa; t++) {
				float[] tloc = snpsAsFloatVector(new BitSet[]{tb.getAllelePresenceForAllSites(t, 0), tb.getAllelePresenceForAllSites(t, 1)}, nsnps);
				float dist1 = getManhattanDistance(loc1, tloc, nsnps);
				float dist2 = getManhattanDistance(loc2, tloc, nsnps);
				if (dist1 <= dist2 && isInCluster1[t] == false) {
					isInCluster1[t] = true;
					loc1 = getMeanLocation(loc1, size1, tloc, true, nsnps);
					loc2 = getMeanLocation(loc2, size2, tloc, false, nsnps);
					noChanges = false;
				} else if (dist1 > dist2 && isInCluster1[t] == true){
					isInCluster1[t] = false;
					loc1 = getMeanLocation(loc1, size1, tloc, false, nsnps);
					loc2 = getMeanLocation(loc2, size2, tloc, true, nsnps);
					noChanges = false;
				}
			}

			if (noChanges) break;
		}
		
		System.out.println("distance between clusters = " + getManhattanDistance(loc1, loc2, nsnps));
		
		//make alignments based on the clusters
		boolean[] isInCluster2 = new boolean[ntaxa];
		for (int t = 0; t < ntaxa; t++) isInCluster2[t] = !isInCluster1[t];
		IdGroup id1 = IdGroupUtils.idGroupSubset(tb.getIdGroup(), isInCluster1);
		IdGroup id2 = IdGroupUtils.idGroupSubset(tb.getIdGroup(), isInCluster2);
		
		Alignment a1 = FilterAlignment.getInstance(tb, id1);
		Alignment a2 = FilterAlignment.getInstance(tb, id2);
		
		return new Alignment[]{a1, a2};
	}
	
	public static Alignment[] getTwoClusters(Alignment inputAlignment, int minGametesPerTaxon) {
		
		//filter out low coverage taxa
		int ntaxa = inputAlignment.getSequenceCount();
		boolean[] include = new boolean[ntaxa];
		
		for (int t = 0; t < ntaxa; t++) {
			if (inputAlignment.getTotalGametesNotMissingForTaxon(t) >= minGametesPerTaxon) include[t] = true;
			else include[t] = false;
		}
		
		Alignment fa = FilterAlignment.getInstance(inputAlignment, IdGroupUtils.idGroupSubset(inputAlignment.getIdGroup(), include));
		Alignment myAlignment = BitAlignment.getInstance(fa, false);
		int ntrials = 5;
		int maxiter = 5;
		
		//if the parents are in the data set use these as seeds
		//if one parent is in the dataset pick the taxon farthest from it as the other seed
		//if neither parent is in the dataset choose random seeds
		ntaxa = myAlignment.getSequenceCount();
		int nsnps = myAlignment.getSiteCount();
		boolean[][] isInCluster1 = new boolean[ntrials][ntaxa];
		int bestTrial = -1;
		float maxDistance = 0;
		
		Random rand = new Random();
		
		float[][] taxaLocs = new float[ntaxa][nsnps];
		
		myAlignment.optimizeForTaxa(null);
		for (int t = 0; t < ntaxa; t++) {
			taxaLocs[t] = snpsAsFloatVector(new BitSet[]{myAlignment.getAllelePresenceForAllSites(t, 0), myAlignment.getAllelePresenceForAllSites(t, 1)}, nsnps);
		}
		
		for (int trial = 0; trial < ntrials; trial++) {
			int seed1 = rand.nextInt(ntaxa);
			int seed2 = -1;
			while (seed2 == -1 || seed1 == seed2) {
				seed2 = rand.nextInt(ntaxa);
			}

			isInCluster1[trial][seed1] = true;
			isInCluster1[trial][seed2] = false;
			
			
			//do initial cluster assignment
			for (int t = 0; t < ntaxa; t++) {
				if (t != seed1 && t != seed2) {
					float dist1 = getManhattanDistance(taxaLocs[seed1], taxaLocs[t], nsnps);
					float dist2 = getManhattanDistance(taxaLocs[seed2], taxaLocs[t], nsnps);
					if (dist1 < dist2) {
						isInCluster1[trial][t] = true;
					} else if (dist1 > dist2){
						isInCluster1[trial][t] = false;
					} else if (rand.nextDouble() > 0.5) {
						isInCluster1[trial][t] = true;
					} else {
						isInCluster1[trial][t] = false;
					}
				}
			}
			
			//update cluster membership until there are no changes or for the maximum number of iterations
			float[][] meanLocs = new float[2][];
			boolean badclusters = false;
			for (int iter = 0; iter < maxiter; iter++) {
				boolean noChanges = true;
				
				int nCluster1 = 0;
				int nCluster2 = 0;
				for (int t = 0; t < ntaxa; t++) {
					if (isInCluster1[trial][t]) nCluster1++;
					else nCluster2++;
				}
				
				if (nCluster1 == 0 || nCluster2 == 0) {
					badclusters = true;
					break;
				}
				
				float[][] cluster1Locs = new float[nCluster1][];
				float[][] cluster2Locs = new float[nCluster2][];
				int countCluster1 = 0;
				int countCluster2 = 0;
				for (int t = 0; t < ntaxa; t++) {
					if (isInCluster1[trial][t]) cluster1Locs[countCluster1++] = taxaLocs[t]; 
					else cluster2Locs[countCluster2++] = taxaLocs[t];
				}
				
				meanLocs = new float[2][];
				meanLocs[0] = getMeanLocation(cluster1Locs);
				meanLocs[1] = getMeanLocation(cluster2Locs);
				for (int t = 0; t < ntaxa; t++) {
					float[] tloc = snpsAsFloatVector(new BitSet[]{myAlignment.getAllelePresenceForAllSites(t, 0), myAlignment.getAllelePresenceForAllSites(t, 1)}, nsnps);
					float dist1 = getManhattanDistance(meanLocs[0], tloc, nsnps);
					float dist2 = getManhattanDistance(meanLocs[1], tloc, nsnps);
					if (dist1 < dist2 && isInCluster1[trial][t] == false) {
						isInCluster1[trial][t] = true;
						noChanges = false;
					} else if (dist1 > dist2 && isInCluster1[trial][t] == true){
						isInCluster1[trial][t] = false;
						noChanges = false;
					}
				}

				if (noChanges) break;
			}
			
			if (badclusters == true) {
				System.out.println("Trial " + trial + ": bad clustering, no distance could be calculated");
			} else {
				float distanceBetweenClusters = getManhattanDistance(meanLocs[0], meanLocs[1], nsnps);
				if (distanceBetweenClusters > maxDistance) {
					maxDistance = distanceBetweenClusters;
					bestTrial = trial;
				}
				System.out.println("Trial " + trial + ": distance between clusters = " + distanceBetweenClusters);
			}
		}

		
		//make alignments based on the clusters
		boolean[] isInCluster2 = new boolean[ntaxa];
		for (int t = 0; t < ntaxa; t++) isInCluster2[t] = !isInCluster1[bestTrial][t];
		IdGroup id1 = IdGroupUtils.idGroupSubset(myAlignment.getIdGroup(), isInCluster1[bestTrial]);
		IdGroup id2 = IdGroupUtils.idGroupSubset(myAlignment.getIdGroup(), isInCluster2);
		
		Alignment a1 = FilterAlignment.getInstance(myAlignment, id1);
		Alignment a2 = FilterAlignment.getInstance(myAlignment, id2);
		
		return new Alignment[]{a1, a2};
	}
	
	public static float[] snpsAsFloatVector(BitSet[] alleles, int nsnps) {
		float[] result = new float[nsnps];
		for (int s = 0; s < nsnps; s++) {
			if (alleles[0].fastGet(s)) {
				result[s] = 2;
				if (alleles[1].fastGet(s)) result[s] = 1;
			} else {
				if (alleles[1].fastGet(s)) result[s] = 0;
//				else result[s] = -1;
				else result[s] = 1;
			}
		}
		return result;
	}
	
	public static float getManhattanDistance(float[] loc, float[] t, int nsnps) {
		float d = 0;
		int nsites = 0;
		for (int s = 0; s < nsnps; s++) {
			if (loc[s] >= 0 && t[s] >= 0) {
				d += Math.abs(loc[s] - t[s]);
				nsites++;
			}
		}
		return d / nsites;
	}
	
	public static float[] getMeanLocation(float[] loc, int[] size, float[] t, boolean add, int nsnps) {
		float[] result = new float[nsnps];
		if (add) {
			for (int s = 0; s < nsnps; s++) {
				if (t[s] >= 0) {
					if (size[s] > 0) {
						result[s] = (loc[s] * size[s] + t[s]) / ((float) (size[s] + 1));
						size[s]++;
					} else {
						result[s] = t[s];
						size[s] = 1;
					}
				} else {
					result[s] = loc[s];
				}
			}
		} else {
			for (int s = 0; s < nsnps; s++) {
				if (t[s] >= 0) {
					if (size[s] > 1) {
						result[s] = (loc[s] * size[s] - t[s]) / ((float) (size[s] - 1));
						size[s]--;
					} else if (size[s] == 1){
						result[s] = 0;
						size[s] = 0;
					}
				} else {
					result[s] = loc[s];
				}
			}
		}
		
		return result;
	}
	
	public static float[] getMeanLocation(float[][] locs) {
		int nsites = locs[0].length;
		int ntaxa = locs.length;
		float[] result = new float[nsites];
		
		for (int s = 0; s < nsites; s++) {
			int count = 0;
			float sum = 0;
			for (int t = 0; t < ntaxa; t++) {
				if (!Float.isNaN(locs[t][s])) {
					count++;
					sum += locs[t][s];
				}
				if (count > 0) result[s] = sum / count;
				else result[s] = Float.NaN;
			}
		}
		return result;
	}
	
	public static void printAlleleStats(Alignment a, String name) {
		int monoCount = 0;
		int polyCount = 0;
		int[] binCount = new int[21];
		int nsites = a.getSiteCount();
		for (int s = 0; s < nsites; s++) {
			if (a.getMajorAlleleFrequency(s) > 0.75) monoCount++;
			else {
				polyCount++;
				int bin = (int) Math.floor(20 * a.getMajorAlleleFrequency(s));
				binCount[bin]++;
			}
		}
		System.out.println(name);
		System.out.println("mono count = " + monoCount + ", poly count = " + polyCount);
		System.out.print("bins: ");
		for (int i = 0; i < 20; i++) System.out.print(" " + binCount[i]);
		System.out.println();
		System.out.println();
	}
	
	public static void mergeNonconsensusFiles(String dir, String match, String outfileName) {
		File[] mergeFiles = filterFiles(dir, match);
		
		int nfiles = mergeFiles.length;
		String[] colLabel = new String[nfiles];
		HashMap<String, String[]> taxonMap = new HashMap<String, String[]>();
		
		String input;
		String[] info;
		for (int f = 0; f < nfiles; f++) {
			System.out.println("processing" + mergeFiles[f].getName());
			try {
				BufferedReader br = new BufferedReader(new FileReader(mergeFiles[f]));
				br.readLine();
				input = br.readLine();
				info = tab.split(input);
				colLabel[f] = info[1];
				while (input != null) {
					info = tab.split(input);
					String[] values = taxonMap.get(info[0]);
					if (values == null) {
						values = new String[nfiles];
						for (int i = 0; i < nfiles; i++) values[i] = "";
						taxonMap.put(info[0], values);
					}
					values[f] = info[2];
					input =  br.readLine();
				}
				
				br.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		
		File outfile = new File(dir, outfileName);
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(outfile));
			StringBuilder sb = new StringBuilder("Taxon");
			for (int i = 0; i < nfiles; i++) {
				sb.append("\t").append(colLabel[i]);
			}
			bw.write(sb.toString());
			bw.newLine();
			LinkedList<String> taxaList = new LinkedList<String>(taxonMap.keySet());
			Collections.sort(taxaList);
			for (String taxon : taxaList) {
				sb = new StringBuilder(taxon);
				String[] values = taxonMap.get(taxon);
				for (int f = 0; f < nfiles; f++) sb.append("\t").append(values[f]);
				bw.write(sb.toString());
				bw.newLine();
			}
			bw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
	public static File[] filterFiles(String dir, String match) {
		File matchdir = new File(dir);
		final String pattern = new String(match);
		File[] filteredFiles = matchdir.listFiles(new FilenameFilter() {
			
			@Override
			public boolean accept(File dir, String name) {
				if (name.matches(pattern)) return true;
				return false;
			}
		});
		return filteredFiles;
	}
	
	public static void mergeFiles(File[] mergeFiles, int idcol, int datacol, int[] colOrder, String outfile) {
		int nfiles = mergeFiles.length;
		String input;
		String[] info;
		if (colOrder == null) colOrder = new int[]{0,2,3,4,5,6,7,8,9,1};
		HashMap<String, String[]> taxonMap = new HashMap<String, String[]>();
		
		for (int f = 0; f < nfiles; f++) {
			System.out.println("processing" + mergeFiles[f].getName());
			try {
				BufferedReader br = new BufferedReader(new FileReader(mergeFiles[f]));
				br.readLine();
				input = br.readLine();
				info = tab.split(input);
				while (input != null) {
					info = tab.split(input);
					String[] values = taxonMap.get(info[idcol]);
					if (values == null) {
						values = new String[nfiles];
						for (int i = 0; i < nfiles; i++) values[i] = "";
						taxonMap.put(info[idcol], values);
					}
					values[f] = info[datacol];
					input =  br.readLine();
				}
				
				br.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(outfile));
			StringBuilder sb = new StringBuilder("Taxon");
			for (int i = 0; i < nfiles; i++) {
				sb.append("\t").append(i + 1);
			}
			sb.append("\t").append("average");
			bw.write(sb.toString());
			bw.newLine();
			LinkedList<String> taxaList = new LinkedList<String>(taxonMap.keySet());
			Collections.sort(taxaList);
			for (String taxon : taxaList) {
				sb = new StringBuilder(taxon);
				String[] values = taxonMap.get(taxon);
				double sum = 0;
				double count = 0;
				for (int f = 0; f < nfiles; f++) {
					sb.append("\t").append(values[f]);
					try {
						sum += Double.parseDouble(values[f]);
						count++;
					} catch (NumberFormatException e) {
						//do nothing
					}
				}
				sb.append("\t").append(sum/count);
				bw.write(sb.toString());
				bw.newLine();
			}
			bw.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
}
