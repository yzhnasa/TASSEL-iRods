package net.maizegenetics.gwas.imputation;

import java.util.Arrays;
import java.util.LinkedList;
import java.util.Random;

import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.FilterAlignment;
import net.maizegenetics.pal.alignment.SBitAlignment;
import net.maizegenetics.pal.alignment.TBitAlignment;
import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.IdGroupUtils;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;

public class ImputationUtils {

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
		TBitAlignment tb;
		if (a instanceof TBitAlignment) tb = (TBitAlignment) a;
		else tb = TBitAlignment.getInstance(a);
		
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
	
	public static float[] snpsAsFloatVector(BitSet[] alleles, int nsnps) {
		float[] result = new float[nsnps];
		for (int s = 0; s < nsnps; s++) {
			if (alleles[0].fastGet(s)) {
				result[s] = 2;
				if (alleles[1].fastGet(s)) result[s] -= 1;
			} else {
				if (alleles[1].fastGet(s)) result[s] = 0;
				else result[s] = -1;
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
}
