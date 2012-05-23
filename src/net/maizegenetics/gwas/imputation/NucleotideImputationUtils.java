package net.maizegenetics.gwas.imputation;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedList;

import org.apache.log4j.Logger;
import org.apache.log4j.xml.DOMConfigurator;

import net.maizegenetics.baseplugins.TreeDisplayPlugin;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.FilterAlignment;
import net.maizegenetics.pal.alignment.MutableNucleotideAlignment;
import net.maizegenetics.pal.alignment.NucleotideAlignmentConstants;
import net.maizegenetics.pal.alignment.SBitAlignment;
import net.maizegenetics.pal.alignment.TBitAlignment;
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

public class NucleotideImputationUtils {
	private static final Logger myLogger = Logger.getLogger(NucleotideImputationUtils.class);
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
	
	//prevents instantiation of this class
	private NucleotideImputationUtils() {}
	
	public static void callParentAlleles(PopulationData popdata, int minAlleleCount, int windowSize, int numberToTry, double cutHeightSnps, double minR) {
		BitSet polybits = whichSitesArePolymorphic(popdata.original, minAlleleCount);
		int[] coreSnps = findCoreSnps(popdata.original, polybits, windowSize, numberToTry, cutHeightSnps);
		String parentA = popdata.parent1;
		String parentC = popdata.parent2;
		
		OpenBitSet ldbits = new OpenBitSet(popdata.original.getSiteCount());
		int n = coreSnps.length;
		for (int i = 0; i < n; i++) ldbits.fastSet(coreSnps[i]);
		
		IdGroup[] taxaGroup =  findTaxaGroups(popdata.original, coreSnps);
		
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
		
		Alignment aAlignment = FilterAlignment.getInstance(popdata.original, parentAGroup);
		Alignment cAlignment = FilterAlignment.getInstance(popdata.original, parentCGroup);
		byte[] Asnp = new byte[popdata.original.getSiteCount()];
		byte[] Csnp = new byte[popdata.original.getSiteCount()];
				
		//set first parent to AA, second parent to CC for snps used to form taxa clusters
		SBitAlignment sbitPopAlignment = SBitAlignment.getInstance(popdata.original);
		MutableNucleotideAlignment parentAlignment = MutableNucleotideAlignment.getInstance(sbitPopAlignment);
		myLogger.info("snps in parent Alignment = " + parentAlignment.getSiteCount());
		int ntaxa = parentAlignment.getSequenceCount();
		
		for (int i = 0; i < coreSnps.length; i++) {
			int snp = coreSnps[i];
			//debug
			int[][] acounts = aAlignment.getAllelesSortedByFrequency(snp);
			int[][] ccounts = cAlignment.getAllelesSortedByFrequency(snp);

			byte alleleA = aAlignment.getMajorAllele(snp);
			byte alleleC = cAlignment.getMajorAllele(snp);
			Asnp[snp] = alleleA;
			Csnp[snp] = alleleC;
			
			for (int t = 0; t < ntaxa; t++) {
				byte[] taxon = popdata.original.getBaseArray(t, snp);
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
		
		//add snps from middle to start; test only polymorphic snps
		LinkedList<Integer> testSnps = new LinkedList<Integer>();
		for (int i = testSize - 1; i >= 0; i--) testSnps.add(coreSnps[i]);
		for (int snp = coreSnps[0] - 1; snp >= 0; snp--) {
			if (polybits.fastGet(snp)) {
				byte[] ac = recodeParentalSnps(snp, testSnps, parentAlignment, minR);
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
				byte[] ac = recodeParentalSnps(snp, testSnps, parentAlignment, minR);
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
		myLogger.info("number of original sites = " + popdata.original.getSiteCount() + ", number of polymorphic sites = " + polybits.cardinality() + ", number of ld sites = " + ldAlignment.getSiteCount());

		popdata.imputed = TBitAlignment.getInstance(ldAlignment);
	}

	public static BitSet whichSitesArePolymorphic(Alignment a, int minAlleleCount) {
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
	
	//returns a byte array containing the A allele as element 0 and the C allele as element 1, or null if the snp is not in LD
	public static byte[] recodeParentalSnps(int snp, LinkedList<Integer> testSnps, MutableNucleotideAlignment snpAlignment, double minr) {
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
		int N = 0;
		N += acount[maxa];
		N += acount[maxc];
		N += ccount[maxa];
		N += ccount[maxc];
		
		double sumx = acount[maxa] + acount[maxc];
		double sumy = acount[maxa] + ccount[maxa];
		double sumxy = acount[maxa];
		
		double r = (sumxy - sumx * sumy / N) / Math.sqrt( (sumx - sumx * sumx / N) * (sumy - sumy * sumy / N) );
		r = Math.abs(r);
		
		//if abs(r) > minr, recode the snp
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

	/**
	 * This function finds a set of snps within a window of the specified size (100) that are in LD with each other. It trys multiple windows and uses the
	 * window that yields the largest number of snps. 
	 * @param a	the input alignment
	 * @param polybits	a BitSet corresponding to SNPs in a, set if a snp is polymorphic. Only polymorphic SNPs will be considered.
	 * @param numberToTry the number of windows to try. The function will use the window returning the largest set of SNPs.
	 * @return indices of the core snps. 
	 */
	public static int[] findCoreSnps(Alignment a, BitSet polybits, int windowSize, int numberToTry, double cutHeightForSnpClusters) {
		
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
			SBitAlignment haplotypeAlignment = SBitAlignment.getInstance(filteredPopAlignment);
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
	
	public static IdGroup[] findTaxaGroups(Alignment a, int[] coreSnps) {
		
		//cluster taxa for these snps to find parental haplotypes (cluster on taxa)
		
		IBSDistanceMatrix dm = new IBSDistanceMatrix(SBitAlignment.getInstance(FilterAlignment.getInstance(a, coreSnps)));
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
			if (groupCount[i] > 5) myLogger.info("Taxa group " + i + " has " + groupCount[i] + " members.");
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
	
	public static void estimateMissingDistances(DistanceMatrix dm) {
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

	public static DistanceMatrix snpDistance(Alignment a) {
		
		int nsnps = a.getSiteCount();
		SimpleIdGroup snpIds = new SimpleIdGroup(nsnps, true);
		double[][] distance = new double[nsnps][nsnps];
		double sum = 0;
		int count = 0;
		for (int i = 0; i < nsnps; i++) {
			for (int j = i; j < nsnps; j++) {
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
	
	public static double computeR(int site1, int site2, Alignment a) {
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
		
		//Explanation of method:
		//  if major site one is x=1, minor site one is x = 0 and for site 2 y = 1 or 0
		//  r = [sum(xy) - sum(x)sum(y)/N] / sqrt[(sum(x) - sum(x)*sum(x)/N) * ((sum(y) - sum(y)*sum(y)/N)]
		//  and sum(x) - sum(x)*sum(x)/N = sum(x)(N - sum(x))/N
		//  because sum(x^2) = sum(x)
		double num = ((double) prodCount - ((double) s1Count * s2Count) / ((double) totalCount));
		double denom = ((double) (s1Count * (totalCount - s1Count))) / ((double) totalCount);
		denom *= ((double) (s2Count * (totalCount - s2Count))) / ((double) totalCount);
		if (denom == 0) return Double.NaN;
		return  num / Math.sqrt(denom);
	}

	public static MutableNucleotideAlignment imputeUsingViterbiFiveState(TBitAlignment a, double probHeterozygous, String familyName) {
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
					myLogger.info("null from " + Byte.toString(base));
				}
				if (isNotMissing.fastGet(s)) {
					obs[nmcount] = obsMap.get(a.getBase(t, s));
					pos[nmcount++] = a.getPositionInLocus(s);
				}
				
			}
			
		}
		
		double phom = (1 - probHeterozygous) / 2;
		double[] pTrue = new double[]{phom, .25*probHeterozygous ,.5 * probHeterozygous, .25*probHeterozygous, phom};
		
		//iterate
		ArrayList<byte[]> bestStates = new ArrayList<byte[]>();
		int[][] previousStateCount = new int[5][3];
		int iter = 0;
		boolean hasNotConverged = true;
		while (iter < maxIterations && hasNotConverged) {
			//apply Viterbi
			myLogger.info("Iteration " + iter++ + " for " + familyName);
			bestStates.clear();
			for (int t = 0; t < ntaxa; t++) {
				tp.setPositions(snpPositions.get(t));
				int nobs = notMissingCount[t];
				if (nobs >= 20) {
					ViterbiAlgorithm va = new ViterbiAlgorithm(nonMissingObs.get(t), tp, ep, pTrue);
					va.calculate();
					bestStates.add(va.getMostProbableStateSequence());
				} else { //do not impute if obs < 20
					myLogger.info("Fewer then 20 observations for " + a.getTaxaName(t));
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
			
			//if the model has converged  or if the max iterations has been reached print tables
			if (!hasNotConverged || iter == maxIterations) {
				StringBuilder sb = new StringBuilder("Family ");
				sb.append(familyName).append(", chromosome ").append(a.getLocusName(0));
				if (iter < maxIterations) {
					sb.append(": EM algorithm converged at iteration ").append(iter).append(".\n");
				} else {
					sb.append(": EM algorithm failed to converge after ").append(iter).append(" iterations.\n");
				}
				
				//print transition counts
				sb = new StringBuilder("Transition counts:\n");
				for (int[] row:transitionCounts) {
					for (int cell:row) {
						sb.append(cell).append("\t");
					}
					sb.append("\n");
				}
				sb.append("\n");
				myLogger.info(sb.toString());
				
				//print transition probabilities
				sb = new StringBuilder("Transition probabilities:\n");
				for (double[] row:transitionProb) {
					for (double cell:row) {
						sb.append(cell).append("\t");
					}
					sb.append("\n");
				}
				sb.append("\n");
				myLogger.info(sb.toString());
				
				//print observation/state counts
				sb = new StringBuilder("Imputation counts, rows=states, columns=observations:\n");
				for (int[] row:emissionCounts) {
					for (int cell:row) {
						sb.append(cell).append("\t");
					}
					sb.append("\n");
				}
				sb.append("\n");
				myLogger.info(sb.toString());
				
				//print emission probabilities
				sb = new StringBuilder("Emission probabilities:\n");
				for (double[] row:emissionProb) {
					for (double cell:row) {
						sb.append(cell).append("\t");
					}
					sb.append("\n");
				}
				sb.append("\n");
				myLogger.info(sb.toString());

			}
			
			
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

	public static void fillGapsInAlignment(PopulationData popdata) {
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

	public static void updateSnpAlignment(PopulationData popdata) {
		//set monomorphic sites to major (or only allele) (or not)
		//set polymorphic sites consistent with flanking markers if equal, unchanged otherwise
		//do not change sites that are not clearly monomorhpic or polymorphic

		MutableNucleotideAlignment mna;
		if (popdata.original instanceof MutableNucleotideAlignment) mna = (MutableNucleotideAlignment) popdata.original;
		else mna = MutableNucleotideAlignment.getInstance(popdata.original);
		
		BitSet isPopSnp = popdata.snpIndex;
		
		int nsites = mna.getSiteCount();
		int ntaxa = mna.getSequenceCount();
		int popSnpCount = 0;
		for (int s = 0; s < nsites; s++) {
			if (isPopSnp.fastGet(s)) {
				int Acall = popdata.alleleA[popSnpCount];
				int Ccall = popdata.alleleC[popSnpCount];
				byte AAcall = (byte) ((Acall << 4) | Acall);
				byte CCcall = (byte) ((Ccall << 4) | Ccall);
				byte ACcall = (byte) ((Acall << 4) | Ccall);
				for (int t = 0; t < ntaxa; t++) {
					byte parentCall = popdata.imputed.getBase(t, popSnpCount);
					if (parentCall == AA) {
						mna.setBase(t, s, AAcall);
					} else if (parentCall == CC) {
						mna.setBase(t, s, CCcall);
					} else if (parentCall == AC || parentCall == CA) {
						mna.setBase(t, s, ACcall);
					} else {
						mna.setBase(t, s, NN);
					}
				}
				popSnpCount++;
			} else { //if the site is monomorphic fill in all the genotypes
				//do nothing for now.
			}
		}
		
		mna.clean();
		popdata.original = mna;
	}

}
