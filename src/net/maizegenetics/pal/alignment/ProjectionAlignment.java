/*
 * ProjectionAlignment
 */
package net.maizegenetics.pal.alignment;

import java.util.Arrays;

import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.OpenBitSet;

/**
 * This class projects high Density markers on large group of taxa through 
 * a look up table system.  The lookup table generally needs be built through
 * some imputation approach.
 * @author ed, terry
 */
public class ProjectionAlignment extends AbstractAlignment {

    private int[][] mySiteBreaks;  //temporary - saving not needed
    private int[][] myHDTaxa;  //taxa ids should be saved
    private int[][] myPosBreaks;  //positions should be saved
    private final SBitAlignment myBaseAlignment;  //high density marker alignment that is being projected.

    public ProjectionAlignment(Alignment hdAlign, IdGroup ldIDGroup) {
        super(ldIDGroup, hdAlign.getAlleleEncodings());
        if (hdAlign instanceof SBitAlignment) {
            myBaseAlignment = (SBitAlignment) hdAlign;
        } else {
            myBaseAlignment = SBitAlignment.getInstance(hdAlign);
        }
        mySiteBreaks = new int[getSequenceCount()][];
        myPosBreaks = new int[getSequenceCount()][];
        myHDTaxa = new int[getSequenceCount()][];
    }

    public void setCompositionOfTaxon(int taxon, int[] posBreaks, int[] hdTaxa) {
        if (posBreaks.length != hdTaxa.length) {
            throw new IllegalArgumentException("ProjectionAlignment: setCompositionOfTaxon: number of positions should equal number taxa.");
        }
        myPosBreaks[taxon] = posBreaks;
        myHDTaxa[taxon] = hdTaxa;
        mySiteBreaks[taxon] = new int[myPosBreaks[taxon].length];
        for (int i = 0; i < myPosBreaks[taxon].length; i++) {
            int site = myBaseAlignment.getSiteOfPhysicalPosition(posBreaks[i], null);
            if (site < 0) {
                site = -(site + 1);
            }
            this.mySiteBreaks[taxon][i] = site;
        }
    }

    public void setCompositionOfTaxon(String taxa, int[] posBreaks, int[] hdTaxa) {
        int taxon = getIdGroup().whichIdNumber(taxa);
        setCompositionOfTaxon(taxon, posBreaks, hdTaxa);
    }
    
    public String getCompositionOfTaxon(int taxon) {
        StringBuilder sb=new StringBuilder(this.getIdGroup().getIdentifier(taxon).getFullName()+"\t");
        for (int i = 0; i < myPosBreaks[taxon].length; i++) {
            sb.append(myPosBreaks[taxon][i]+":");
            sb.append(mySiteBreaks[taxon][i]+":");
            sb.append(myHDTaxa[taxon][i]+"\t");
        }
        return sb.toString();
    }
    
    public void reportPAComposition() {
        for (int i = 0; i < this.getSequenceCount(); i++) {
            System.out.println(getCompositionOfTaxon(i));
        }
    }

    private int translateTaxon(int taxon, int site) {
        if(mySiteBreaks[taxon]==null) { 
            if(site%100==0) System.out.printf("Taxon null:%d site:%d %n",taxon,site);
            return 0;
        }
        int b = Arrays.binarySearch(mySiteBreaks[taxon], site);
        if (b<0) {
            b=-(b+2);  //this will not work if it does not start with zero.
        }
 //       if(myHDTaxa[taxon]==null) return 0;//this should be a missing taxon.
        return myHDTaxa[taxon][b];
    }

    @Override
    public String getBaseAsString(int taxon, int site) {
        return myBaseAlignment.getBaseAsString(translateTaxon(taxon, site), site);
    }

    @Override
    public String getDiploidAsString(int site, byte value) {
        return myBaseAlignment.getDiploidAsString(site, value);
    }

    @Override
    /**
     * This is the slow implementation of this.  Most of these should be buffered bit
     * sets.  Note the imputation of the taxon is likely to be the the same for 64 or more sites
     * in a row (potentially, 10,000s of sites in many cases).
     */
    public byte getBase(int taxon, int site) {
        return myBaseAlignment.getBase(translateTaxon(taxon, site), site);
    }

    @Override
    public boolean isSBitFriendly() {
        return false;
    }

    @Override
    public boolean isTBitFriendly() {
        return false;
    }

    @Override
    public byte[] getBaseArray(int taxon, int site) {
        return myBaseAlignment.getBaseArray(translateTaxon(taxon, site), site);
    }

    @Override
    public String[] getBaseAsStringArray(int taxon, int site) {
        return myBaseAlignment.getBaseAsStringArray(translateTaxon(taxon, site), site);
    }

    @Override
    public byte[] getBaseRow(int taxon) {

        int numBreaks = mySiteBreaks[taxon].length;
        byte[] result = new byte[myNumSites];
        for (int i = 0; i < numBreaks - 1; i++) {
            int hdTaxon = myHDTaxa[taxon][i];
            for (int j = mySiteBreaks[taxon][i]; j < mySiteBreaks[taxon][i + 1]; j++) {
                result[j] = myBaseAlignment.getBase(hdTaxon, j);
            }
        }

        int hdTaxon = myHDTaxa[taxon][numBreaks - 1];
        for (int j = mySiteBreaks[taxon][numBreaks - 1], n = getSiteCount(); j < n; j++) {
            result[j] = myBaseAlignment.getBase(hdTaxon, j);
        }

        return result;
    }

    @Override
    // TERRY - This Could be Optimized like getBaseRow()
    public byte[] getBaseRange(int taxon, int startSite, int endSite) {

        byte[] result = new byte[endSite - startSite];
        for (int i = startSite; i < endSite; i++) {
            result[i] = getBase(taxon, i);
        }
        return result;

    }

    @Override
    public byte getBase(int taxon, Locus locus, int physicalPosition) {
        return getBase(taxon, getSiteOfPhysicalPosition(physicalPosition, locus));
    }

    @Override
    public BitSet getAllelePresenceForAllSites(int taxon, int alleleNumber) {
        throw new UnsupportedOperationException();
    }

    @Override
    public BitSet getAllelePresenceForAllTaxa(int site, int alleleNumber) {
        BitSet baseBitSet = myBaseAlignment.getAllelePresenceForAllTaxa(site, alleleNumber);
        BitSet result = new OpenBitSet(getSequenceCount());
        for (int i = 0, n = getSequenceCount(); i < n; i++) {
            int index = translateTaxon(i, site);
            if (baseBitSet.fastGet(index)) {
                result.fastSet(i);
            }
        }
        return result;
    }

    @Override
    public long[] getAllelePresenceForSitesBlock(int taxon, int alleleNumber, int startBlock, int endBlock) {
        throw new UnsupportedOperationException();
    }

    @Override
    public int getIndelSize(int site) {
        return myBaseAlignment.getIndelSize(site);
    }

    @Override
    public boolean isIndel(int site) {
        return isIndel(site);
    }

    @Override
    public float getSiteScore(int seq, int site) {
        return myBaseAlignment.getSiteScore(translateTaxon(seq, site), site);
    }

    @Override
    public boolean hasSiteScores() {
        return false;
    }

    @Override
    public SITE_SCORE_TYPE getSiteScoreType() {
        return Alignment.SITE_SCORE_TYPE.None;
    }

    @Override
    public boolean hasReference() {
        return myBaseAlignment.hasReference();
    }

    @Override
    public byte getReferenceAllele(int site) {
        return myBaseAlignment.getReferenceAllele(site);
    }

    @Override
    public byte[] getReference() {
        return myBaseAlignment.getReference();
    }

    @Override
    public byte[] getReference(int startSite, int endSite) {
        return myBaseAlignment.getReference(startSite, endSite);
    }

    @Override
    public boolean isPhased() {
        return false;
    }

    @Override
    public boolean isPositiveStrand(int site) {
        return myBaseAlignment.isPositiveStrand(site);
    }

    @Override
    public String getGenomeAssembly() {
        return myBaseAlignment.getGenomeAssembly();
    }

    @Override
    public GeneticMap getGeneticMap() {
        return myBaseAlignment.getGeneticMap();
    }

    @Override
    public int[] getPhysicalPositions() {
        return myBaseAlignment.getPhysicalPositions();
    }

    @Override
    public int getSiteCount() {
        return myBaseAlignment.getSiteCount();
    }

    @Override
    public int getPositionInLocus(int site) {
        return myBaseAlignment.getPositionInLocus(site);
    }

    @Override
    public Locus getLocus(int site) {
        return myBaseAlignment.getLocus(site);
    }

    @Override
    public Locus[] getLoci() {
        return myBaseAlignment.getLoci();
    }

    @Override
    public int getNumLoci() {
        return myBaseAlignment.getNumLoci();
    }

    @Override
    public int[] getLociOffsets() {
        return myBaseAlignment.getLociOffsets();
    }

    @Override
    public int getLocusSiteCount(Locus locus) {
        return myBaseAlignment.getLocusSiteCount(locus);
    }

    @Override
    public String[] getSNPIDs() {
        return myBaseAlignment.getSNPIDs();
    }

    @Override
    public String getSNPID(int site) {
        return myBaseAlignment.getSNPID(site);
    }

    @Override
    public int getSiteOfPhysicalPosition(int physicalPosition, Locus locus) {
        return myBaseAlignment.getSiteOfPhysicalPosition(physicalPosition, locus);
    }

    @Override
    public byte getPositionType(int site) {
        return myBaseAlignment.getPositionType(site);
    }

    @Override
    public byte[] getPositionTypes() {
        return myBaseAlignment.getPositionTypes();
    }

    @Override
    public boolean retainsRareAlleles() {
        return myBaseAlignment.retainsRareAlleles();
    }

    @Override
    public String[][] getAlleleEncodings() {
        return myBaseAlignment.getAlleleEncodings();
    }

    @Override
    public String[] getAlleleEncodings(int site) {
        return myBaseAlignment.getAlleleEncodings(site);
    }

    @Override
    public String getBaseAsString(int site, byte value) {
        return myBaseAlignment.getBaseAsString(site, value);
    }

    @Override
    public int getMaxNumAlleles() {
        return myBaseAlignment.getMaxNumAlleles();
    }

    @Override
    public int getTotalGametesNotMissingForTaxon(int taxon) {

        int numBreaks = mySiteBreaks[taxon].length;
        int result = 0;
        for (int i = 0; i < numBreaks - 1; i++) {
            int hdTaxon = myHDTaxa[taxon][i];
            for (int j = mySiteBreaks[taxon][i]; j < mySiteBreaks[taxon][i + 1]; j++) {
                byte[] current = myBaseAlignment.getBaseArray(hdTaxon, j);
                if (current[0] != Alignment.UNKNOWN_ALLELE) {
                    result++;
                }
                if (current[1] != Alignment.UNKNOWN_ALLELE) {
                    result++;
                }
            }
        }

        int hdTaxon = myHDTaxa[taxon][numBreaks - 1];
        for (int j = mySiteBreaks[taxon][numBreaks - 1], n = getSiteCount(); j < n; j++) {
            byte[] current = myBaseAlignment.getBaseArray(hdTaxon, j);
            if (current[0] != Alignment.UNKNOWN_ALLELE) {
                result++;
            }
            if (current[1] != Alignment.UNKNOWN_ALLELE) {
                result++;
            }
        }

        return result;

    }

    // TERRY - This Needs Work.
    public Object[][] getMajorMinorCounts() {

        if (myAlleleStates.length != 1) {
            return new Object[0][0];
        }

        long[][] counts = new long[16][16];

        if (myMaxNumAlleles >= 2) {
            for (int site = 0; site < myNumSites; site++) {
                byte indexI = myAlleles[site][0];
                byte indexJ = myAlleles[site][1];
                if (indexJ == UNKNOWN_ALLELE) {
                    indexJ = indexI;
                }
                counts[indexI][indexJ]++;
            }
        } else {
            for (int site = 0; site < myNumSites; site++) {
                byte indexI = myAlleles[site][0];
                counts[indexI][indexI]++;
            }
        }

        int numAlleles = 0;
        for (byte x = 0; x < 16; x++) {
            for (byte y = 0; y < 16; y++) {
                if (counts[x][y] != 0) {
                    numAlleles++;
                }
            }
        }

        Object[][] result = new Object[2][numAlleles];
        int nextResult = 0;
        for (byte x = 0; x < 16; x++) {
            for (byte y = 0; y < 16; y++) {
                if (counts[x][y] != 0) {
                    result[0][nextResult] = getBaseAsString(0, x) + ":" + getBaseAsString(0, y);
                    result[1][nextResult++] = counts[x][y];
                }
            }
        }

        boolean change = true;
        while (change) {

            change = false;

            for (int k = 0; k < numAlleles - 1; k++) {

                if ((Long) result[1][k] < (Long) result[1][k + 1]) {

                    Object temp = result[0][k];
                    result[0][k] = result[0][k + 1];
                    result[0][k + 1] = temp;

                    Object tempCount = result[1][k];
                    result[1][k] = result[1][k + 1];
                    result[1][k + 1] = tempCount;

                    change = true;
                }
            }

        }

        return result;
    }

    @Override
    public int[][] getAllelesSortedByFrequency(int site) {

        int maxNumAlleles = myBaseAlignment.getMaxNumAlleles();
        int numDataRows = myBaseAlignment.getTotalNumAlleles();
        BitSet[] data = new BitSet[numDataRows];
        for (int i = 0; i < numDataRows; i++) {
            data[i] = getAllelePresenceForAllTaxa(site, i);
        }
        // TERRY - need to think about this.
        byte[] alleles = myBaseAlignment.myAlleles[site];

        int[] counts = new int[16];
        for (int i = 0; i < numDataRows; i++) {
            byte indexI;
            if ((retainsRareAlleles()) && (i == maxNumAlleles)) {
                indexI = Alignment.RARE_ALLELE;
            } else {
                indexI = alleles[i];
            }
            counts[indexI] += (int) data[i].cardinality() * 2;
            for (int j = i + 1; j < numDataRows; j++) {
                byte indexJ;
                if ((retainsRareAlleles()) && (j == maxNumAlleles)) {
                    indexJ = Alignment.RARE_ALLELE;
                } else {
                    indexJ = alleles[j];
                }
                int ijHet = (int) OpenBitSet.intersectionCount(data[i], data[j]);
                counts[indexI] -= ijHet;
                counts[indexJ] -= ijHet;
            }
        }

        int numAlleles = 0;
        for (byte x = 0; x < Alignment.UNKNOWN_ALLELE; x++) {
            if (counts[x] != 0) {
                numAlleles++;
            }
        }

        int current = 0;
        int[][] result = new int[2][numAlleles];
        for (byte x = 0; x < Alignment.UNKNOWN_ALLELE; x++) {
            if (counts[x] != 0) {
                result[0][current] = x;
                result[1][current++] = counts[x];
            }
        }

        boolean change = true;
        while (change) {

            change = false;

            for (int k = 0; k < numAlleles - 1; k++) {

                if (result[1][k] < result[1][k + 1]) {

                    int temp = result[0][k];
                    result[0][k] = result[0][k + 1];
                    result[0][k + 1] = temp;

                    int tempCount = result[1][k];
                    result[1][k] = result[1][k + 1];
                    result[1][k + 1] = tempCount;

                    change = true;
                }
            }

        }

        return result;

    }

    @Override
    public int getTotalNumAlleles() {
        return myBaseAlignment.getTotalNumAlleles();
    }
}
