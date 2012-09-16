package net.maizegenetics.gwas.imputation;

import java.awt.Frame;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

import javax.swing.ImageIcon;

import org.apache.log4j.Logger;

import net.maizegenetics.baseplugins.FixedEffectLMPlugin;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.FilterAlignment;
import net.maizegenetics.pal.ids.Identifier;
import net.maizegenetics.pal.ids.SimpleIdGroup;
import net.maizegenetics.pal.popgen.LinkageDisequilibrium;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.BitUtil;
import net.maizegenetics.util.OpenBitSet;

public class QualityChecksPlugin extends AbstractPlugin {
    private static final Logger myLogger = Logger.getLogger(QualityChecksPlugin.class);
	private String pedigreeFile;
	private int windowSizeForR2;
	private double minNonMissingProportionForTaxon = 0.05;
	private double minNonMissingProportionForSNP = 0.05;
	
	public enum checkType {AVERAGE_R2, NONCONSENSUS_PROPORTION, NONCONSENSUS_SITES};
	private ArrayList<checkType> analysisList = new ArrayList<checkType>();
	
	public QualityChecksPlugin(Frame parentFrame) {
		super(parentFrame, false);
	}
	
	@Override
	public DataSet performFunction(DataSet input) {
		List<Datum> datumList = input.getDataOfType(Alignment.class);
		
		for (Datum datum : datumList) {
			Alignment anAlignment = (Alignment) datum.getData();
			anAlignment = preFilterAlignment(anAlignment);
			for (checkType analysis : analysisList) {
				switch (analysis) {
				case AVERAGE_R2:
					calculateAverageR2ForSnps(anAlignment);
					break;
				case NONCONSENSUS_PROPORTION:
					calculateProportionNonConsensusPerTaxon(anAlignment);
					break;
				case NONCONSENSUS_SITES:
					saveNonConsensusSites(anAlignment);
					break;
				}
			}
		}
		return null;
	}
	
	public Alignment preFilterAlignment(Alignment align) {
		int ntaxa = align.getSequenceCount();
		int nsites = align.getSiteCount();
		int ngametes = 2 * nsites;
		int minTaxaGametes = (int) Math.ceil(ngametes * minNonMissingProportionForTaxon);
		int minSiteGametes = (int) Math.ceil(ngametes * minNonMissingProportionForSNP);
		
		//create list of taxa with too much missing data
		LinkedList<Identifier> taxaDiscardList = new LinkedList<Identifier>();
		for (int t = 0; t < ntaxa; t++) {
			if (align.getTotalGametesNotMissingForTaxon(t) < minTaxaGametes) taxaDiscardList.add(align.getIdGroup().getIdentifier(t));
		}
		if (taxaDiscardList.size() > 0) {
			myLogger.info("\nThe following taxa will not be included in the analysis because the proportion of nonMissing data is below " + minNonMissingProportionForTaxon + ":\n");
			for (Identifier id:taxaDiscardList) myLogger.info(id.getFullName());
			myLogger.info("\n");
			
			Identifier[] ids = new Identifier[taxaDiscardList.size()];
			taxaDiscardList.toArray(ids);
			
			align = FilterAlignment.getInstanceRemoveIDs(align, new SimpleIdGroup(ids));
		}
		
		//number of non-missing values per site
		int[] sitesToKeep = new int[nsites];
		int nsitesKept = 0;
		for (int s = 0; s < nsites; s++) {
			if (align.getTotalGametesNotMissing(s) >= minSiteGametes) sitesToKeep[nsitesKept++] = s;
		}
		
		if (nsitesKept < nsites) {
			myLogger.info(nsitesKept + " sites had more than " + minSiteGametes + " and were retained.");
			sitesToKeep = Arrays.copyOf(sitesToKeep, nsitesKept);
			align = FilterAlignment.getInstance(align, sitesToKeep);
		}
		
		return align;
	}
	
	private double[] calculateAverageR2ForSnps(Alignment align) {
		align.optimizeForSites(null);
		int nsites = align.getSiteCount();
		int ntaxa = align.getSequenceCount();
		double[] avgRsq = new double[nsites];
		
		for (int s = 0; s < nsites; s++) {
			int start = Math.max(s - windowSizeForR2, 0);
			int end = Math.min(nsites - 1, s + windowSizeForR2);
			double sum = 0;
			double count = 0;
            BitSet sMj = align.getAllelePresenceForAllTaxa(s, 0);
            BitSet sMn = align.getAllelePresenceForAllTaxa(s, 1);

			for (int i = start; i <= end; i++) {
				if (i != s) {
					int[][] contig = new int[2][2];
		            BitSet iMj = align.getAllelePresenceForAllTaxa(i, 0);
		            BitSet iMn = align.getAllelePresenceForAllTaxa(i, 1);
		            contig[0][0] = (int) OpenBitSet.intersectionCount(sMj, iMj);
		            contig[1][0] = (int) OpenBitSet.intersectionCount(sMn, iMj);
		            contig[0][1] = (int) OpenBitSet.intersectionCount(sMj, iMn);
		            contig[1][1] = (int) OpenBitSet.intersectionCount(sMn, iMn);
		            double rsq = calculateRSqr(contig[0][0], contig[1][0], contig[0][1], contig[1][1], 4);
		            if (!Double.isNaN(rsq)) {
		            	sum += rsq;
		            	count++;
		            }
				}
			}
			
			if (count > 0) {
				avgRsq[s] = sum / count;
			} else {
				avgRsq[s] = Double.NaN;
			}
			
		}
		
		return avgRsq;
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

	private double[] calculateProportionNonConsensusPerTaxon(Alignment align) {
		double maxMaf = 0.05;
		int ntaxa = align.getSequenceCount();
		int nsites = align.getSiteCount();
		
		OpenBitSet lowmaf = new OpenBitSet(nsites);
		for (int s = 0; s < nsites; s++) {
			if (align.getMinorAlleleFrequency(s) < maxMaf) lowmaf.set(s);
		}
		
		align.optimizeForTaxa(null);
		double[] proportionMinor = new double[ntaxa];
		double nmono = lowmaf.cardinality();
		for (int t = 0; t < ntaxa; t++) {
			long minorCount = OpenBitSet.intersectionCount(lowmaf, align.getAllelePresenceForAllSites(t, 1));
			proportionMinor[t] = minorCount / nmono;
		}
		
		return proportionMinor;
	}
	
	private void saveNonConsensusSites(Alignment align) {
		
	}
	
	public void addAnalysis(checkType analysis) {
		analysisList.add(analysis);
	}
	
	@Override
	public void setParameters(String[] args) {
		// TODO Auto-generated method stub
		super.setParameters(args);
	}

	@Override
	public ImageIcon getIcon() {
		return null;
	}

	@Override
	public String getButtonName() {
		return null;
	}

	@Override
	public String getToolTipText() {
		return null;
	}

	
}
