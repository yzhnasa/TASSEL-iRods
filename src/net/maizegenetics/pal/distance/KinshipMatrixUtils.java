package net.maizegenetics.pal.distance;


import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.datatype.DataType;
import net.maizegenetics.pal.datatype.IUPACNucleotides;
import net.maizegenetics.pal.datatype.TextDataType;

public class KinshipMatrixUtils {
	
	/**
	 * Calculates an IBS kinship matrix. Each element equals the probability that an allele drawn at random from a site 
	 * in one genotype has the same state as an allele drawn from the same site in a second genotype. 
	 * For two genotypes, only sites with no missing data for either genotype are considered. 
	 * If there is missing data the resulting kinship matrix might not be non-negative definite.
	 * @param genotypes an Alignment
	 * @return an IBS kinship matrix for this alignment
	 */
	public static DistanceMatrix getIBSMatrix(Alignment genotypes) {
		int ntaxa = genotypes.getSequenceCount();
		int nsites = genotypes.getSiteCount();
		DataType dt = genotypes.getDataType();
		
		//get ploidy and appropriate divisor for the matrix
		int ploidy = getPloidy(genotypes);
		double divisor = ploidy * ploidy;
		
		int size = ntaxa * (ntaxa + 1) / 2;
		int[] nonmissing = new int[size];
		int[] ibs = new int[size];
		char[] sites = new char[ntaxa];
		
		for (int s = 0; s < nsites; s++) {
			for (int t = 0; t < ntaxa; t++) sites[t] = genotypes.getBaseChar(t, s);
			int count = 0;
			for (int i = 0; i < ntaxa; i++) {
				char sitescore = sites[i];
				if (!dt.isUnknownChar(sitescore)) {
					nonmissing[count]++;
					ibs[count] += dt.getDiploidIdentity(sitescore, sitescore);
					count++;
					for (int j = i + 1; j < ntaxa; j++) {
						if (!dt.isUnknownChar(sites[j])) {
							nonmissing[count]++;
							ibs[count] += dt.getDiploidIdentity(sitescore, sites[j]);
							count++;
						} else {
							count++;
						}
					}
				} else {
					count += ntaxa - i;
				}
			}
		}
		
		double[][] ibsmatrix = new double[ntaxa][ntaxa];
		int count = 0;
		for (int i = 0; i < ntaxa; i++) {
			ibsmatrix[i][i] = ((double) ibs[count]) / ((double) nonmissing[count]) / divisor;
			count++;
			for (int j = i + 1; j < ntaxa; j++) {
				double value = ((double) ibs[count]) / ((double) nonmissing[count]) / divisor;
				ibsmatrix[i][j] = value;
				ibsmatrix[j][i] = value;
				count++;
			}
		}
		
		return new DistanceMatrix(ibsmatrix, genotypes.getIdGroup());
	}
	
	public static int getPloidy(Alignment a) {
		DataType dt = a.getDataType();
		if (dt instanceof IUPACNucleotides) return 2;
		if (dt instanceof TextDataType) {
			String geno = dt.getFormattedString(a.getBaseChar(0, 0));
			String[] alleles = geno.split(":");
			return alleles.length;
		}
		return -1;
	}
	
}
