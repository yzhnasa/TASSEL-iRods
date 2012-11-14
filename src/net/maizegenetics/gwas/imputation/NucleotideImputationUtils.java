package net.maizegenetics.gwas.imputation;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import net.maizegenetics.baseplugins.ConvertSBitTBitPlugin;

import org.apache.log4j.Logger;

import net.maizegenetics.baseplugins.TreeDisplayPlugin;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.FilterAlignment;
import net.maizegenetics.pal.alignment.Locus;
import net.maizegenetics.pal.alignment.MutableNucleotideAlignment;
import net.maizegenetics.pal.alignment.MutableSingleEncodeAlignment;
import net.maizegenetics.pal.alignment.NucleotideAlignmentConstants;
import net.maizegenetics.pal.alignment.BitAlignment;
import net.maizegenetics.pal.distance.DistanceMatrix;
import net.maizegenetics.pal.distance.IBSDistanceMatrix;
import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.IdGroupUtils;
import net.maizegenetics.pal.ids.SimpleIdGroup;
import net.maizegenetics.pal.math.GammaFunction;
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
		
		//debug
		examineTaxaClusters(popdata.original, polybits);
		
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
		//SBitAlignment sbitPopAlignment = SBitAlignment.getInstance(popdata.original);
                Alignment sbitPopAlignment = ConvertSBitTBitPlugin.convertAlignment(popdata.original, ConvertSBitTBitPlugin.CONVERT_TYPE.sbit, null);
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

		//popdata.imputed = TBitAlignment.getInstance(ldAlignment);
                popdata.imputed = ConvertSBitTBitPlugin.convertAlignment(ldAlignment, ConvertSBitTBitPlugin.CONVERT_TYPE.tbit, null);
	}

	public static void callParentAllelesByWindow(PopulationData popdata, double maxMissing, double minMaf, int windowSize, double minR) {
		
		BitSet polybits;
		double segratio = popdata.contribution1;
		if (segratio == 0.5 || segratio == 0.25 || segratio == 0.75) {
			polybits = whichSitesSegregateCorrectly(popdata.original, maxMissing, segratio);
		} else {
			polybits = whichSitesArePolymorphic(popdata.original, maxMissing, minMaf);
		}
		myLogger.info("polybits cardinality = " + polybits.cardinality());

		OpenBitSet filteredBits = whichSnpsAreFromSameTag(popdata.original, 0.8);
		filteredBits.and(polybits);
		System.out.println("filteredBits.cardinality = " + filteredBits.cardinality());
		
		BitSet ldFilteredBits;
		if (minR > 0) {
			int halfWindow = windowSize / 2;
			ldFilteredBits = ldfilter(popdata.original, halfWindow, minR, filteredBits);
		} else {
			ldFilteredBits = filteredBits;
		}
		myLogger.info("ldFilteredBits.cardinality = " + ldFilteredBits.cardinality());
		
		//get het mask
//		double phet = 0.075;
//		double phet = (1 - popdata.inbredCoef)/2;
//		BitSet[] hetMask = hetMasker(popdata.original, phet);

		int nsites = popdata.original.getSiteCount();
		int ntaxa = popdata.original.getSequenceCount();
		popdata.alleleA = new byte[nsites];
		popdata.alleleC = new byte[nsites];
		popdata.snpIndex = new OpenBitSet(nsites);
		for (int s = 0; s < nsites; s++) {
			popdata.alleleA[s] = Alignment.UNKNOWN_ALLELE;
			popdata.alleleC[s] = Alignment.UNKNOWN_ALLELE;
		}
		
		int[] parentIndex = new int[2];
		parentIndex[0] = popdata.original.getIdGroup().whichIdNumber(popdata.parent1);
		parentIndex[1] = popdata.original.getIdGroup().whichIdNumber(popdata.parent2);
	
		//iterate through windows
		Alignment[] prevAlignment = null;
		int[][] snpIndices = getWindows(filteredBits, windowSize);
		boolean append = false;
		
		int nWindows = snpIndices.length;
		for (int w = 0; w < nWindows; w++) {
//		for (int[] snpIndex : snpIndices) {
			int[] snpIndex;
			if (append) {
				int n1 = snpIndices[w-1].length;
				int n2 = snpIndices[w].length;
				snpIndex = new int[n1 + n2];
				System.arraycopy(snpIndices[w-1], 0, snpIndex, 0, n1);
				System.arraycopy(snpIndices[w], 0, snpIndex, n1, n2);
				append = false;
			} else {
				snpIndex = snpIndices[w];
			}
			
			//mask the hets
//			MutableNucleotideAlignment mna = MutableNucleotideAlignment.getInstance(FilterAlignment.getInstance(popdata.original, snpIndex));
//			int n = snpIndex.length;
//			byte missing = NucleotideAlignmentConstants.getNucleotideDiploidByte('N');
//			for (int i = 0; i < n; i++) {
//				for (int t = 0; t < ntaxa; t++) {
//					if (hetMask[t].fastGet(snpIndex[i])) mna.setBase(t, i, missing);
//				}
//			}
//			mna.clean();
			
			Alignment windowAlignment = FilterAlignment.getInstance(popdata.original, snpIndex);
			LinkedList<Integer> snpList = new LinkedList<Integer>(); //snpList is a list of snps (indices) in this window
			for (int s:snpIndex) snpList.add(s);
			
			Alignment[] taxaAlignments = getTaxaGroupAlignments(windowAlignment, parentIndex, snpList);
			
			if (taxaAlignments == null) {
				append = true;
			} else {
				//are groups in this alignment correlated with groups in the previous alignment
				double r = 0;
				if (prevAlignment != null) {
					r = getIdCorrelation(new IdGroup[][] {{prevAlignment[0].getIdGroup(), prevAlignment[1].getIdGroup()},{taxaAlignments[0].getIdGroup(), taxaAlignments[1].getIdGroup()}});
					myLogger.info("For " + popdata.name + " the window starting at " + popdata.original.getSNPID(snpIndex[0]) + ", r = " + r + " , # of snps in alignment = " + snpList.size());
				} else {
					myLogger.info("For " + popdata.name + " the window starting at " + popdata.original.getSNPID(snpIndex[0]) + ", # of snps in alignment = " + snpList.size());
				}
				
				checkAlignmentOrderIgnoringParents(taxaAlignments, popdata, r); //if r is negative switch alignment order (which will result in a positive r) 
				
				//debug -check upgma tree
//				int[] selectSnps = new int[snpList.size()];
//				int cnt = 0;
//				for (Integer s : snpList) selectSnps[cnt++] = s;
//				SBitAlignment sba = SBitAlignment.getInstance(FilterAlignment.getInstance(popdata.original, selectSnps));
//				IBSDistanceMatrix dm = new IBSDistanceMatrix(sba);
//				estimateMissingDistances(dm);
//				Tree myTree = new UPGMATree(dm);
//				TreeDisplayPlugin tdp = new TreeDisplayPlugin(null, true);
//				tdp.performFunction(new DataSet(new Datum("Snp Tree", myTree, "Snp Tree"), null));
				
				prevAlignment = taxaAlignments;
				callParentAllelesUsingTaxaGroups(popdata, taxaAlignments, snpList);
			}
			
		}
		
		myLogger.info("number of called snps = " + popdata.snpIndex.cardinality());
		
		//create the imputed array with A/C calls
		int nsnps = (int) popdata.snpIndex.cardinality();
		ntaxa = popdata.original.getSequenceCount();
		nsites = popdata.original.getSiteCount();
		int[] snpIndex = new int[nsnps];
		int snpcount = 0;
		for (int s = 0; s < nsites; s++) {
			if (popdata.snpIndex.fastGet(s)) snpIndex[snpcount++] = s;
		}
		
		snpIndex = Arrays.copyOf(snpIndex, snpcount);
		Alignment target = FilterAlignment.getInstance(popdata.original, snpIndex);
		MutableNucleotideAlignment mna = MutableNucleotideAlignment.getInstance(target);
		
		nsnps = snpIndex.length;
		for (int s = 0; s < nsnps; s++) {
			byte Aallele = popdata.alleleA[snpIndex[s]];
			byte Callele = popdata.alleleC[snpIndex[s]];
			byte genotypeA = (byte) (Aallele << 4 | Aallele);
			byte genotypeC = (byte) (Callele << 4 | Callele);
			byte het1 = (byte) (Aallele << 4 | Callele);
			byte het2 = (byte) (Callele << 4 | Aallele);
			for (int t = 0; t < ntaxa; t++) {
				byte val = mna.getBase(t, s);
				if (val == genotypeA) {
					mna.setBase(t, s, AA);
				} else if (val == genotypeC) {
					mna.setBase(t, s, CC);
				} else if (val == het1 || val == het2) {
					mna.setBase(t, s, AC);
				} else {
					mna.setBase(t, s, NN);
				}
			}
		}
		mna.clean();
		popdata.imputed = BitAlignment.getInstance(mna, true); 
	}

	public static void callParentAllelesByWindowForBackcrosses(PopulationData popdata, double maxMissing, double minMaf, int windowSize, double minR) {
		
		BitSet polybits = whichSitesSegregateCorrectly(popdata.original, maxMissing, .25);
		myLogger.info("polybits cardinality = " + polybits.cardinality());

		OpenBitSet filteredBits = whichSnpsAreFromSameTag(popdata.original, 0.8);
		filteredBits.and(polybits);
		System.out.println("filteredBits.cardinality = " + filteredBits.cardinality());
		
		BitSet ldFilteredBits;
		if (minR > 0) {
			int halfWindow = windowSize / 2;
			ldFilteredBits = ldfilter(popdata.original, halfWindow, minR, filteredBits);
		} else {
			ldFilteredBits = filteredBits;
		}
		myLogger.info("ldFilteredBits.cardinality = " + ldFilteredBits.cardinality());
		
		int nsites = popdata.original.getSiteCount();
		int ntaxa = popdata.original.getSequenceCount();
		popdata.alleleA = new byte[nsites];
		popdata.alleleC = new byte[nsites];
		popdata.snpIndex = new OpenBitSet(nsites);
		for (int s = 0; s < nsites; s++) {
			if(ldFilteredBits.fastGet(s)) {
				popdata.alleleA[s] = popdata.original.getMajorAllele(s);
				popdata.alleleC[s] = popdata.original.getMinorAllele(s);
				popdata.snpIndex.fastSet(s);
			} 
		}
		
		
		myLogger.info("number of called snps = " + popdata.snpIndex.cardinality());
		
		//create the imputed array with A/C calls
		int nsnps = (int) popdata.snpIndex.cardinality();
		ntaxa = popdata.original.getSequenceCount();
		nsites = popdata.original.getSiteCount();
		int[] snpIndex = new int[nsnps];
		int snpcount = 0;
		for (int s = 0; s < nsites; s++) {
			if (popdata.snpIndex.fastGet(s)) snpIndex[snpcount++] = s;
		}
		
		snpIndex = Arrays.copyOf(snpIndex, snpcount);
		Alignment target = FilterAlignment.getInstance(popdata.original, snpIndex);
		MutableNucleotideAlignment mna = MutableNucleotideAlignment.getInstance(target);
		
		nsnps = snpIndex.length;
		for (int s = 0; s < nsnps; s++) {
			byte Aallele = popdata.alleleA[snpIndex[s]];
			byte Callele = popdata.alleleC[snpIndex[s]];
			byte genotypeA = (byte) (Aallele << 4 | Aallele);
			byte genotypeC = (byte) (Callele << 4 | Callele);
			byte het1 = (byte) (Aallele << 4 | Callele);
			byte het2 = (byte) (Callele << 4 | Aallele);
			for (int t = 0; t < ntaxa; t++) {
				byte val = mna.getBase(t, s);
				if (val == genotypeA) {
					mna.setBase(t, s, AA);
				} else if (val == genotypeC) {
					mna.setBase(t, s, CC);
				} else if (val == het1 || val == het2) {
					mna.setBase(t, s, AC);
				} else {
					mna.setBase(t, s, NN);
				}
			}
		}
		mna.clean();
		popdata.imputed = BitAlignment.getInstance(mna, true); 
	}
	
	public static void checkAlignmentOrder(Alignment[] alignments, PopulationData family, double r) {
		boolean swapAlignments = false;
		boolean parentsInSameGroup = false;
		boolean parentsInWrongGroups = false;
		double minR = -0.05;
		
		int p1group, p2group;
		
		//if parent 1 is in alignment 0, p1group = 0
		//if parent1 is in alignment 1, p1group = 1
		//if parent 1 is not in either alignment p1group = -1
		//likewise for parent 2
		if (alignments[0].getIdGroup().whichIdNumber(family.parent1) > -1) p1group = 0;
		else if (alignments[1].getIdGroup().whichIdNumber(family.parent1) > -1) p1group = 1;
		else p1group = -1;

		if (alignments[1].getIdGroup().whichIdNumber(family.parent2) > -1) p2group = 1;
		else if (alignments[0].getIdGroup().whichIdNumber(family.parent2) > -1) p2group = 0;
		else p2group = -1;

		//find out if parents are either in the same group or if they are in groups opposite expected
		//parent 1 is expected to be in group 0 since it was used to seed group 0
		if (p1group == 0) {
			if (p2group == 0) {
				parentsInSameGroup = true;
			} else if (p2group == 1) {
				if (r < 0) parentsInWrongGroups = true;
			} else {
				if (r < 0) parentsInWrongGroups = true;
			}
		} else if (p1group == 1) {
			if (p2group == 0) {
				if (r > 0) parentsInWrongGroups = true;
			} else if (p2group == 1) {
				parentsInSameGroup = true;
			} else {
				if (r > 0) parentsInWrongGroups = true;
			}
		} else {
			if (p2group == 0) {
				if (r > 0) parentsInWrongGroups = true;
			} else if (p2group == 1) {
				if (r < 0) parentsInWrongGroups = true;
			} else {
				//do nothing
			}
		}
		
		//r is the correlation between taxa assignments in this window and the previous one
		//if r is negative then parental assignments should be swapped, which will make it positive
		//this can happen when the parents are not in the data set and the assignment of parents to group is arbitrary
		if (r < minR) swapAlignments = true;
		if (swapAlignments) {
			Alignment temp = alignments[0];
			alignments[0] = alignments[1];
			alignments[1] = temp;
		}
		
		if (parentsInSameGroup) {
			myLogger.warn("Both parents in the same group for family " + family.name + " at " + alignments[0].getSNPID(0));
		}
		
		if (parentsInWrongGroups) {
			myLogger.warn("Parents in unexpected group for family " + family.name + " at " + alignments[0].getSNPID(0));
		}
	}
	
	public static void checkAlignmentOrderIgnoringParents(Alignment[] alignments, PopulationData family, double r) {
		boolean swapAlignments = false;
		double minR = -0.05;
		
		int p1group, p2group;
		
		//r is the correlation between taxa assignments in this window and the previous one
		//if r is negative then parental assignments should be swapped, which will make it positive
		//this can happen when the parents are not in the data set and the assignment of parents to group is arbitrary
		if (r < minR) swapAlignments = true;
		if (swapAlignments) {
			Alignment temp = alignments[0];
			alignments[0] = alignments[1];
			alignments[1] = temp;
		}
		
	}
	
	public static int[][] getWindows(BitSet ispoly, int windowSize) {
		int npoly = (int) ispoly.cardinality();
		int nsnps = (int) ispoly.size();
		int nwindows = npoly/windowSize;
		int remainder = npoly % windowSize;
		if (remainder > windowSize/2) nwindows++; //round up
		int[][] windows = new int[nwindows][];
		int setsize = npoly/nwindows;
		
		int windowCount = 0;
		int snpCount = 0;
		int polyCount = 0;
		while (snpCount < nsnps && windowCount < nwindows) {
			int numberLeft = npoly - polyCount;
			if (numberLeft < setsize * 2) setsize = numberLeft;
			int[] set = new int[setsize];
			int setcount = 0;
			while (setcount < setsize && snpCount < nsnps) {
				if (ispoly.fastGet(snpCount)) {
					set[setcount++] = snpCount;
					polyCount++;
				}
				snpCount++;
			}
			windows[windowCount++] = set;
		}
		
		return windows;
	}
	
	/**
	 * @param family	a PopulationData object containing information for this family
	 * @param taxaGroups	an array of two alignments corresponding to two clusters of taxa
	 * @param snpList	the list of snps to be called
	 */
	public static void callParentAllelesUsingTaxaGroups(PopulationData family, Alignment[] taxaGroups, LinkedList<Integer> snpList) {
		int nsnps = taxaGroups[0].getSiteCount();
		Iterator<Integer> snpit = snpList.iterator();
		for ( int s = 0; s < nsnps; s++) {
			byte[] major = new byte[2];
			major[0] = taxaGroups[0].getMajorAllele(s);
			major[1] = taxaGroups[1].getMajorAllele(s);
			Integer snpIndex = snpit.next();
			if(major[0] != Alignment.UNKNOWN_ALLELE && major[1] != Alignment.UNKNOWN_ALLELE && major[0] != major[1]) {
				family.alleleA[snpIndex] = major[0];
				family.alleleC[snpIndex] = major[1];
				family.snpIndex.fastSet(snpIndex);
			}
		}
	}
	
	public static double getIdCorrelation(IdGroup[][] id) {
		double[][] counts = new double[2][2];
		counts[0][0] = IdGroupUtils.getCommonIds(id[0][0], id[1][0]).getIdCount();
		counts[0][1] = IdGroupUtils.getCommonIds(id[0][0], id[1][1]).getIdCount();
		counts[1][0] = IdGroupUtils.getCommonIds(id[0][1], id[1][0]).getIdCount();
		counts[1][1] = IdGroupUtils.getCommonIds(id[0][1], id[1][1]).getIdCount();
		double num = counts[0][0] * counts[1][1] - counts[0][1] * counts[1][0];
		double p1 = counts[0][0] + counts[0][1];
		double q1 = counts[1][0] + counts[1][1];
		double p2 =  counts[0][0] + counts[1][0];
		double q2 =  counts[0][1] + counts[1][1];
		return num / Math.sqrt(p1 * q1 * p2 * q2);
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
	
	public static BitSet whichSitesArePolymorphic(Alignment a, double maxMissing, double minMaf) {
		//which sites are polymorphic? minor allele count > 2 and exceed the minimum allele count
		int nsites = a.getSiteCount();
		int ntaxa = a.getSequenceCount();
		double totalgametes = 2 * ntaxa;
		OpenBitSet polybits = new OpenBitSet(nsites);
		for (int s = 0; s < nsites; s++) {
			int[][] freq = a.getAllelesSortedByFrequency(s);
			int ngametes = a.getTotalGametesNotMissing(s);
			double pMissing = (totalgametes - ngametes) / totalgametes;
			if (freq[1].length > 1 && freq[1][1] > 2 && pMissing <= maxMissing && a.getMinorAlleleFrequency(s) > minMaf) {
				polybits.fastSet(s);
			}
		}
		return polybits;
	}
	
	public static BitSet whichSitesSegregateCorrectly(Alignment a, double maxMissing, double ratio) {
		int nsites = a.getSiteCount();
		int ntaxa = a.getSequenceCount();
		double totalgametes = 2 * ntaxa;
		OpenBitSet polybits = new OpenBitSet(nsites);
		for (int s = 0; s < nsites; s++) {
			int[][] freq = a.getAllelesSortedByFrequency(s);
			int ngametes = a.getTotalGametesNotMissing(s);
			double pMissing = (totalgametes - ngametes) / totalgametes;
			if (freq[1].length > 1 && pMissing <= maxMissing) {
				int Mj = freq[1][0];
				int Mn = freq[1][1];
				double pmono = binomialProbability(Mj + Mn, Mn, 0.002);
				double pquarter = binomialProbability(Mj + Mn, Mn, 0.25);
				double phalf = binomialProbability(Mj + Mn, Mn, 0.5);
				if (ratio == 0.25 || ratio == 0.75) {
					if (pquarter / (pmono + phalf) > 4) polybits.fastSet(s);
				} else {
					if (phalf / (pmono + pquarter) > 4) polybits.fastSet(s);
				}
			}
		}
		return polybits;

	}
	
	private static double binomialProbability(int trials, int successes, double pSuccess) {
		double n = trials;
		double k = successes;

		double logprob = GammaFunction.lnGamma(n + 1.0) -
		GammaFunction.lnGamma(k + 1.0) - GammaFunction.lnGamma(n - k + 1.0) + k * Math.log(pSuccess) + (n - k) * Math.log(1 - pSuccess);

		return Math.exp(logprob);
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
	
	public static IdGroup[] findTaxaGroups(Alignment a, int[] coreSnps) {
		
		//cluster taxa for these snps to find parental haplotypes (cluster on taxa)
		
		//IBSDistanceMatrix dm = new IBSDistanceMatrix(SBitAlignment.getInstance(FilterAlignment.getInstance(a, coreSnps)));
                IBSDistanceMatrix dm = new IBSDistanceMatrix(ConvertSBitTBitPlugin.convertAlignment(FilterAlignment.getInstance(a, coreSnps), ConvertSBitTBitPlugin.CONVERT_TYPE.sbit, null));
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
		TreeDisplayPlugin tdp = new TreeDisplayPlugin(null, true);
		tdp.performFunction(new DataSet(new Datum("Snp Tree", myTree, "Snp Tree"), null));

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
	
	public static Alignment[] getTaxaGroupAlignments(Alignment a, int[] parentIndex, LinkedList<Integer> snpIndices) {
		
		//cluster taxa for these snps to find parental haplotypes (cluster on taxa)
		Alignment[] taxaClusters = ImputationUtils.getTwoClusters(a, 20);
		LinkedList<Integer> originalList = new LinkedList<Integer>(snpIndices);
		int nsites = a.getSiteCount();
		boolean[] include = new boolean[nsites];
		int[] includedSnps = new int[nsites];
		int snpcount = 0;
		for (int s = 0; s < nsites; s++) {
			Integer snpIndex = snpIndices.remove();
			if ( taxaClusters[0].getMajorAllele(s) != taxaClusters[1].getMajorAllele(s) ) {
//				if ( taxaClusters[0].getMajorAllele(s) != Alignment.UNKNOWN_ALLELE && taxaClusters[1].getMajorAllele(s) != Alignment.UNKNOWN_ALLELE && 
//						taxaClusters[0].getMajorAlleleFrequency(s) > .6 &&  taxaClusters[1].getMajorAlleleFrequency(s) > .6) {
//					include[s] = true;
//					includedSnps[snpcount++] = s;
//					snpIndices.add(snpIndex);
//				} else include[s] = false;
				include[s] = true;
				includedSnps[snpcount++] = s;
				snpIndices.add(snpIndex);
			} else {
//				System.out.println("alleles equal at " + s);
//				include[s] = false;
				include[s] = false;
			}
		}
		
		if (snpcount > 5) {
			includedSnps = Arrays.copyOf(includedSnps, snpcount);
			if (snpcount == nsites) return taxaClusters;
			else return ImputationUtils.getTwoClusters(FilterAlignment.getInstance(a, includedSnps), 20);
		} else {
			snpIndices.clear();
			snpIndices.addAll(originalList);
			return taxaClusters;
		}
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

	public static MutableNucleotideAlignment imputeUsingViterbiFiveState(Alignment a, double probHeterozygous, String familyName) {
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

			//debug - print observation/state counts
//			StringBuilder strb = new StringBuilder("Imputation counts, rows=states, columns=observations:\n");
//			for (int[] row:emissionCounts) {
//				for (int cell:row) {
//					strb.append(cell).append("\t");
//				}
//				strb.append("\n");
//			}
//			strb.append("\n");
//			myLogger.info(strb.toString());

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
			double[] rowSums = new double[5];
			double total = 0;
			for (int row = 0; row < 5; row++) {
				double rowsum = 0;
				for (int col = 0; col < 3; col++) rowsum += emissionCounts[row][col];
				for (int col = 0; col < 3; col++) emissionProb[row][col] = ((double) emissionCounts[row][col]) / rowsum;
				rowSums[row] = rowsum;
				total += rowsum;
			}
			ep.setEmissionProbability(emissionProb);
			
			//re-estimate pTrue
			for (int i = 0; i < 5; i++) {
				pTrue[i
				      ] = rowSums[i] / total;
			}
			
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
				sb = new StringBuilder("Transition counts from row to column:\n");
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
		for (int t = 0; t < ntaxa; t++) {
			int prevsite = -1;
			byte prevValue = -1;
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
		popdata.imputed = a;		
	}

	public static Alignment convertParentCallsToNucleotides(PopulationData popdata) {
		//set monomorphic sites to major (or only allele) (or not)
		//set polymorphic sites consistent with flanking markers if equal, unchanged otherwise
		//do not change sites that are not clearly monomorhpic or polymorphic

		MutableNucleotideAlignment mna = MutableNucleotideAlignment.getInstance(popdata.imputed);
		BitSet isPopSnp = popdata.snpIndex;
		if (isPopSnp.cardinality() != mna.getSiteCount()) myLogger.info("size of imputed snps not equal to snpIndex cardinality in convertParentCallsToNucleotides.");
		int nsites = (int) isPopSnp.capacity();
		double ngametes = nsites * 2;
		int ntaxa = mna.getSequenceCount();
		int imputedSnpCount = 0;
		
		for (int s = 0; s < nsites; s++) {
			if (isPopSnp.fastGet(s)) {
				int Acall = popdata.alleleA[s];
				int Ccall = popdata.alleleC[s];
				byte AAcall = (byte) ((Acall << 4) | Acall);
				byte CCcall = (byte) ((Ccall << 4) | Ccall);
				byte ACcall = (byte) ((Acall << 4) | Ccall);
				for (int t = 0; t < ntaxa; t++) {
					byte parentCall = popdata.imputed.getBase(t, imputedSnpCount);
					if (parentCall == AA) {
						mna.setBase(t, imputedSnpCount, AAcall);
					} else if (parentCall == CC) {
						mna.setBase(t, imputedSnpCount, CCcall);
					} else if (parentCall == AC || parentCall == CA) {
						mna.setBase(t, imputedSnpCount, ACcall);
					} else {
						mna.setBase(t, imputedSnpCount, NN);
					}
				}
				imputedSnpCount++;
			}
		}
		
		mna.clean();
		myLogger.info("Original alignment updated for family " + popdata.name + " chromosome " + popdata.original.getLocusName(0) + ".\n");
		return mna;
	}

	public static void examineTaxaClusters(Alignment a, BitSet polybits) {
		int nsnps = a.getSiteCount();
		int ntaxa = a.getSequenceCount();
		int sitecount = 500;
		int window = 200;
		while (sitecount < nsnps) {
			int[] snpndx = new int[window];
			int snpcount = 0;
			while (snpcount < window && sitecount < nsnps) {
				if (polybits.fastGet(sitecount)) snpndx[snpcount++] = sitecount;
				sitecount++;
			}
			if (sitecount < nsnps) {
				Alignment subAlignment = BitAlignment.getInstance(FilterAlignment.getInstance(a, snpndx), true);
				IBSDistanceMatrix dm = new IBSDistanceMatrix(subAlignment);
				estimateMissingDistances(dm);

				Tree myTree = new UPGMATree(dm);
				TreeClusters tc = new TreeClusters(myTree);
				
				int[] groups = null;
				int[] order = null;
				int ngrp = 2;
				while (true) {
					groups = tc.getGroups(ngrp);
					int[] grpSize = new int[ngrp];
					for (int g:groups) grpSize[g]++;
					order = ImputationUtils.reverseOrder(grpSize);
					if (((double) grpSize[order[0]]) / ((double) grpSize[order[1]]) < 2.0) {
						String[] taxaA = new String[grpSize[order[0]]];
						String[] taxaB = new String[grpSize[order[1]]];
						int cntA = 0;
						int cntB = 0;
						for (int t = 0; t < ntaxa; t++) {
							String taxon = myTree.getIdentifier(t).getFullName();
							if (groups[t] == order[0]) taxaA[cntA++] = taxon;
							else if (groups[t] == order[1]) taxaB[cntB++] = taxon;
						}
						Alignment alignA = FilterAlignment.getInstance(subAlignment, new SimpleIdGroup(taxaA));
						Alignment alignB = FilterAlignment.getInstance(subAlignment, new SimpleIdGroup(taxaB));
						boolean[] include = new boolean[window];
						for (int s = 0; s < window; s++) {
							if (alignA.getMajorAllele(s) != alignB.getMajorAllele(s)) {
								if ( ((double) alignA.getMajorAlleleCount(s))/((double) alignA.getMinorAlleleCount(s)) > 2.0 && ((double) alignB.getMajorAlleleCount(s))/((double) alignB.getMinorAlleleCount(s)) > 2.0) {
									include[s] = true;
								} else include[s] = false;
							} else {
								System.out.println("alleles equal at " + s);
								include[s] = false;
							}
						}
						
						int ngoodsnps = 0;
						for (boolean b:include) if (b) ngoodsnps++;
						
						int[] goodSnpIndex = new int[ngoodsnps];
						int cnt = 0;
						for (int s = 0; s < window; s++) {
							if (include[s]) goodSnpIndex[cnt++] = s;
						}
						
						IBSDistanceMatrix dm2 = new IBSDistanceMatrix(BitAlignment.getInstance(FilterAlignment.getInstance(subAlignment, goodSnpIndex), true));
						estimateMissingDistances(dm2);
						Tree thisTree = new UPGMATree(dm2);
						
						//display the tree 
						TreeDisplayPlugin tdp = new TreeDisplayPlugin(null, true);
						tdp.performFunction(new DataSet(new Datum("Snp Tree", thisTree, "Snp Tree"), null));
						
						
						System.out.println("n good snps = " + ngoodsnps);
						break;
					}
//					else if (ngrp > 2) {
//						if (((double) grpSize[order[0]]) / ((double) (grpSize[order[1]] + grpSize[order[2]])) < 1.5) {
//							
//						}
//					}
					ngrp++;
					
					if (ngrp > 20) break;
				}
				
				//display the tree 
				TreeDisplayPlugin tdp = new TreeDisplayPlugin(null, true);
				tdp.performFunction(new DataSet(new Datum("Snp Tree", myTree, "Snp Tree"), null));
			}
		}
		
	}
	
	/**
	 * SNPs on the same tag will have correlated errors. While alignments do not have information about which SNPs come from the same tag, SNPs from one tag will be <64 bp distant. 
	 * They will also be highly correlated. This function tests removes the second of any pair of SNPs that could come from the same tag.
	 * @param alignIn	the input Alignment
	 * @return	an alignment with the one of any correlated pairs of SNPs removed
	 */
	public static Alignment removeSNPsFromSameTag(Alignment alignIn, double minRsq) {
		Alignment sba = ConvertSBitTBitPlugin.convertAlignment(alignIn, ConvertSBitTBitPlugin.CONVERT_TYPE.sbit, null);
		//if (alignIn instanceof SBitAlignment) sba = (SBitAlignment) alignIn;
		//else sba = SBitAlignment.getInstance(alignIn);
		
		int nsites = sba.getSiteCount();
		int ntaxa = sba.getSequenceCount();
		int firstSite = 0;
		int[] sitesSelected = new int[nsites];
		int selectCount = 0;
		sitesSelected[selectCount++] = 0;
		String firstSnpLocus = sba.getLocus(0).getName();
		int firstSnpPos = sba.getPositionInLocus(0);
		while (firstSite < nsites - 1) {
			int nextSite = firstSite + 1;
			int nextSnpPos = sba.getPositionInLocus(nextSite);
			String nextSnpLocus = sba.getLocus(nextSite).getName();
			while (firstSnpLocus.equals(nextSnpLocus) && nextSnpPos - firstSnpPos < 64) {
				//calculate r^2 between snps
	            BitSet rMj = sba.getAllelePresenceForAllTaxa(firstSite, 0);
	            BitSet rMn = sba.getAllelePresenceForAllTaxa(firstSite, 1);
	            BitSet cMj = sba.getAllelePresenceForAllTaxa(nextSite, 0);
	            BitSet cMn = sba.getAllelePresenceForAllTaxa(nextSite, 1);
	            int n = 0;
	            int[][] contig = new int[2][2];
	            n += contig[0][0] = (int) OpenBitSet.intersectionCount(rMj, cMj);
	            n += contig[1][0] = (int) OpenBitSet.intersectionCount(rMn, cMj);
	            n += contig[0][1] = (int) OpenBitSet.intersectionCount(rMj, cMn);
	            n += contig[1][1] = (int) OpenBitSet.intersectionCount(rMn, cMn);
				
				double rsq = calculateRSqr(contig[0][0], contig[0][1], contig[1][0], contig[1][1], 2);
				if (Double.isNaN(rsq) || rsq >= minRsq) sitesSelected[selectCount++] = nextSite;
				nextSite++;
				nextSnpPos = sba.getPositionInLocus(nextSite);
				nextSnpLocus = sba.getLocus(nextSite).getName();
			}
			firstSite = nextSite;
			firstSnpLocus = nextSnpLocus;
			firstSnpPos = nextSnpPos;
			sitesSelected[selectCount++] = firstSite;
		}
		
		return FilterAlignment.getInstance(sba, sitesSelected);
	}
	
	public static OpenBitSet whichSnpsAreFromSameTag(Alignment alignIn, double minRsq) {
		Alignment sba = ConvertSBitTBitPlugin.convertAlignment(alignIn, ConvertSBitTBitPlugin.CONVERT_TYPE.sbit, null);
		//if (alignIn instanceof SBitAlignment) sba = (SBitAlignment) alignIn;
		//else sba = SBitAlignment.getInstance(alignIn);
		
		int nsites = sba.getSiteCount();
		int firstSite = 0;
		OpenBitSet isSelected = new OpenBitSet(nsites);
		isSelected.fastSet(0);
		Locus firstSnpLocus = sba.getLocus(0);
		int firstSnpPos = sba.getPositionInLocus(0);
		while (firstSite < nsites - 1) {
			int nextSite = firstSite + 1;
			int nextSnpPos = sba.getPositionInLocus(nextSite);
			Locus nextSnpLocus = sba.getLocus(nextSite);
			while (firstSnpLocus.equals(nextSnpLocus) && nextSnpPos - firstSnpPos < 64) {
				//calculate r^2 between snps
	            BitSet rMj = sba.getAllelePresenceForAllTaxa(firstSite, 0);
	            BitSet rMn = sba.getAllelePresenceForAllTaxa(firstSite, 1);
	            BitSet cMj = sba.getAllelePresenceForAllTaxa(nextSite, 0);
	            BitSet cMn = sba.getAllelePresenceForAllTaxa(nextSite, 1);
	            int[][] contig = new int[2][2];
	            contig[0][0] = (int) OpenBitSet.intersectionCount(rMj, cMj);
	            contig[1][0] = (int) OpenBitSet.intersectionCount(rMn, cMj);
	            contig[0][1] = (int) OpenBitSet.intersectionCount(rMj, cMn);
	            contig[1][1] = (int) OpenBitSet.intersectionCount(rMn, cMn);
				
				double rsq = calculateRSqr(contig[0][0], contig[0][1], contig[1][0], contig[1][1], 2);
				//if rsq cannot be calculated or rsq is less than the minimum rsq for a snp to be considered highly correlated, select this snp
				if (Double.isNaN(rsq) || rsq < minRsq) {
					isSelected.fastSet(nextSite);
					break;
				}
				nextSite++;
				if (nextSite >= nsites) break;
				nextSnpPos = sba.getPositionInLocus(nextSite);
				nextSnpLocus = sba.getLocus(nextSite);
			}
			firstSite = nextSite;
			if (firstSite < nsites) {
				firstSnpLocus = nextSnpLocus;
				firstSnpPos = nextSnpPos;
				isSelected.fastSet(firstSite);
			}
		}
		
		return isSelected;
	}
	
    static double calculateRSqr(int countAB, int countAb, int countaB, int countab, int minTaxaForEstimate) {
        //this is the Hill & Robertson measure as used in Awadella Science 1999 286:2524
        double freqA, freqB, rsqr, nonmissingSampleSize;
        nonmissingSampleSize = countAB + countAb + countaB + countab;
        if (nonmissingSampleSize < minTaxaForEstimate) {
            return Double.NaN;
        }
        freqA = (double) (countAB + countAb) / nonmissingSampleSize;
        freqB = (double) (countAB + countaB) / nonmissingSampleSize;

        //Through missing data & incomplete datasets some alleles can be fixed this returns missing value
        if ((freqA == 0) || (freqB == 0) || (freqA == 1) || (freqB == 1)) {
            return Double.NaN;
        }

        rsqr = ((double) countAB / nonmissingSampleSize) * ((double) countab / nonmissingSampleSize);
        rsqr -= ((double) countaB / nonmissingSampleSize) * ((double) countAb / nonmissingSampleSize);
        rsqr *= rsqr;
        rsqr /= freqA * (1 - freqA) * freqB * (1 - freqB);
        return rsqr;
    }

    public static OpenBitSet filterSnpsOnLDandDistance(Alignment alignIn) {
    	return null;
    }

    public static BitSet[] hetMasker(Alignment a, double estimatedHetFraction) {
    	//set up Viterbi algorithm
    	//states in hom, het
    	//observations = hom, het, missing
    	
    	TransitionProbability tp = new TransitionProbability();
    	int nsites = a.getSiteCount();
    	int ntaxa = a.getSequenceCount();
    	int chrlen = a.getPositionInLocus(nsites - 1) - a.getPositionInLocus(0);
    	
    	double phet = estimatedHetFraction;
    	int totalTransitions = (nsites - 1) * ntaxa /10;
    	int hetHet = (int) Math.floor(phet*totalTransitions);
    	int hetHom = 2 * ntaxa;
    	int[][] transCount = new int[][]{{totalTransitions - hetHet - 2*hetHom, hetHom},{hetHom, hetHet}};
    	
    	tp.setTransitionCounts(transCount, chrlen, ntaxa);
    	tp.setPositions(a.getPhysicalPositions());
    	
    	//count number of het loci
    	int hetCount = 0;
    	int nonMissingCount = 0;
    	for (int s = 0; s < nsites; s++) {
    		hetCount += a.getHeterozygousCount(s);
    		nonMissingCount += a.getTotalGametesNotMissing(s) / 2;
    	}
    	
    	double estimatedPhet = ((double) hetCount) / ((double) nonMissingCount); 
    	double expectedPhet = 0.08;
    	double hetGivenHet = Math.min(.9, estimatedPhet/expectedPhet);
    	
    	double[] initialProb = new double[]{1 - expectedPhet, expectedPhet};
    	
    	try {
    		a.optimizeForTaxa(null);
    	} catch (Exception e) {
    		a = BitAlignment.getInstance(a, false);
    	}
    	
    	BitSet[] taxaStates = new BitSet[ntaxa];
    	for (int t = 0; t < ntaxa; t++) {
    		BitSet major = a.getAllelePresenceForAllSites(t, 0);
    		BitSet minor = a.getAllelePresenceForAllSites(t, 1);
    		OpenBitSet notMissing = new OpenBitSet(major.getBits(), major.getNumWords());
    		notMissing.union(minor);
    		OpenBitSet het = new OpenBitSet(major.getBits(), major.getNumWords());
    		het.intersect(minor);
    		
    		double nNotMissing = notMissing.cardinality();
    		byte[] obs = new byte[nsites];
    		for (int s = 0; s < nsites; s++) {
    			if (notMissing.fastGet(s)) {
    				if (het.fastGet(s)) obs[s] = 1;
    				else obs[s] = 0;
    			} else {
    				obs[s] = 2;
    			}
    		}
    		
    		EmissionProbability ep = new EmissionProbability();
    		
    		//states rows, observations columns
    		double[][] probMatrix = new double[2][3];
    		
    		// homozygous state
    		probMatrix[0][2] = (nsites - nNotMissing) / ((double) nsites); //observe missing
    		probMatrix[0][0] = .998 * (1 - probMatrix[0][2]); //observe hom
    		probMatrix[0][1] = .002 * (1 - probMatrix[0][2]);//observe het
    		//heterozygous state
    		probMatrix[1][2] = probMatrix[0][2]; //observe missing
    		probMatrix[1][0] = (1 - hetGivenHet) * (1 - probMatrix[0][2]); //observe hom
    		probMatrix[1][1] = hetGivenHet * (1 - probMatrix[0][2]);//observe het
    		
    		ep.setEmissionProbability(probMatrix);
    		
    		ViterbiAlgorithm va = new ViterbiAlgorithm(obs, tp, ep, initialProb);
    		va.calculate();
    		byte[] states = va.getMostProbableStateSequence();
    		OpenBitSet bitStates = new OpenBitSet(nsites);
    		for (int s = 0; s < nsites; s++) if (states[s] == 1) bitStates.fastSet(s);
    		taxaStates[t] = bitStates;
    	}
    	
    	return taxaStates;
    }
    
    public static BitSet ldfilter(Alignment a, int window, double minR) {
    	try {
    		a.optimizeForSites(null);
    	} catch(Exception e) {
    		a = BitAlignment.getInstance(a, true);
    	}

    	int nsites = a.getSiteCount();
    	OpenBitSet ldbits = new OpenBitSet(nsites);
    	for (int s = 0; s < nsites; s++) {
    		double avgr = neighborLD(a, s, window);
    		if (avgr >= minR) ldbits.fastSet(s);
    	}
    	return ldbits;
    }
    
    public static BitSet ldfilter(Alignment a, int window, double minR, BitSet filterBits) {
    	try {
    		a.optimizeForSites(null);
    	} catch(Exception e) {
    		a = BitAlignment.getInstance(a, true);
    	}
    	
    	int nsites = a.getSiteCount();
    	OpenBitSet ldbits = new OpenBitSet(nsites);
    	for (int s = 0; s < nsites; s++) {
    		if (filterBits.fastGet(s)) {
        		double avgr = neighborLD(a, s, window, filterBits);
        		if (avgr >= minR) ldbits.fastSet(s);
    		}
    	}
    	return ldbits;
    }
    
    public static double neighborLD(Alignment a, int snp, int window) {
    	double sumr = 0;
    	double count = 0;
    	int nsites = a.getSiteCount();
    	for (int i = 0; i < window; i++) {
    		int site = snp - window;
    		if (site > -1) {
    			sumr += Math.abs(computeR(snp, snp + i, a));
    			count++;
    		}
    		site = snp + window;
    		if (site < nsites) {
        		sumr += Math.abs(computeR(snp, snp - i, a));
    			count++;
    		}
    		
    	}
    	return sumr / count;
    }
    
    public static double neighborLD(Alignment a, int snp, int window, BitSet filterBits) {
    	double sumr = 0;
    	int nsites = a.getSiteCount();
    	
    	//test next <window> sites or to the end of the chromosome
    	int site = snp + 1;
    	int siteCount = 0;
    	int sitesTestedCount = 0;
    	while (siteCount < window && site < nsites) {
    		if (filterBits.fastGet(site)) {
    			sumr += Math.abs(computeR(snp, site, a));
    			siteCount++;
    			sitesTestedCount++;
    		}
    		site++;
    	}
    	
    	//test previous <window> sites or to the beginning of the chromosome
    	site = snp - 1;
    	siteCount = 0;
    	while (siteCount < window && site >= 0) {
    		if (filterBits.fastGet(site)) {
    			sumr += Math.abs(computeR(snp, site, a));
    			siteCount++;
    			sitesTestedCount++;
    		}
    		site--;
    	}
    	
    	return sumr / sitesTestedCount;
    }

    public MutableNucleotideAlignment imputeUsingPhasedViterbi(PopulationData family, double probHeterozygous, String familyName) {
    	//call parent alleles by window first for intial estimate
    	//impute haplotypes using Viterbi without iteration
    	
    	PhasedEmissionProbability pep = new PhasedEmissionProbability();
    	
    	//calculate parent allele probability given data
    	//P(A=N1|y) = P(y|A=N1)/[P(y|A=N1) + P(y|A=N2)]
    	int nsites = family.imputed.getSiteCount();
    	int ntaxa = family.imputed.getSequenceCount();
    	double[][][] alleleProb = new double[nsites][4][4]; //1st dim is site (node), 2nd dim is haplotype, 3rd dim is nucleotide (A,C,G,T)
    	for (int s = 0; s < nsites; s++) {
    		family.imputed.getSNPID(s);
    		
    		for (int t = 0; t < ntaxa; t++) {
    			
    		}
    	}
    	pep.setParentAlleleProbability(null);
    	
    	return null;
    }
}


