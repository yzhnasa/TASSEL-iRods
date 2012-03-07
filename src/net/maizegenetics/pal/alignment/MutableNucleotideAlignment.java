/*
 * MutableNucleotideAlignment
 */
package net.maizegenetics.pal.alignment;

import cern.colt.GenericSorting;
import cern.colt.Swapper;
import cern.colt.function.IntComparator;

import java.util.Arrays;
import java.util.HashMap;

import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.Identifier;
import net.maizegenetics.pal.ids.SimpleIdGroup;

/**
 *
 * @author edbuckler
 */
public class MutableNucleotideAlignment extends AbstractAlignment implements MutableAlignment {

    private byte[][] seq;
    private Locus[] myLoci;
    private int[] lociForSite;
    private int[] variableSites;
    private byte[] strand;
    private byte[] sitePrefix;
    private int siteNumber;
    private int nextFreeSite = 0;
    protected byte[] majorAlleles, minorAlleles;
    private HashMap<String, Integer> locusNameToLociIndex;

    public MutableNucleotideAlignment(IdGroup taxa, int maxSites, Locus[] loci) {
        super(taxa);
        myLoci = loci;
        siteNumber = maxSites;
        seq = new byte[getSequenceCount()][siteNumber];
        variableSites = new int[siteNumber];
        strand = new byte[siteNumber];
        lociForSite = new int[siteNumber];
        sitePrefix = new byte[siteNumber];
        for (int t = 0; t < getSequenceCount(); t++) {
            Arrays.fill(seq[t], Alignment.UNKNOWN_DIPLOID_ALLELE);
        }
        Arrays.fill(variableSites, Integer.MAX_VALUE);
        Arrays.fill(strand, Byte.MAX_VALUE);
        Arrays.fill(lociForSite, Integer.MAX_VALUE);
        Arrays.fill(sitePrefix, (byte) 'S');
        locusNameToLociIndex = new HashMap<String, Integer>();
        for (int i = 0; i < loci.length; i++) {
            locusNameToLociIndex.put(loci[i].getName(), i);
        }
        initMajorMinorAlleles();
    }

    public MutableNucleotideAlignment(String[] taxaNames, int maxSites, Locus[] loci) {
        this(new SimpleIdGroup(taxaNames), maxSites, loci);
    }

    public MutableNucleotideAlignment(Alignment a) {
        this(a.getIdGroup(), a.getSiteCount(), a.getLoci());
        for (int s = 0; s < a.getSiteCount(); s++) {
            setLocusOfSite(s, a.getLocusName(s));
            setPositionOfSite(s, a.getPositionInLocus(s));
            setSitePrefix(s, (byte) a.getSNPID(s).charAt(0));
            for (int i = 0; i < a.getSequenceCount(); i++) {
                setBase(i, s, a.getBase(i, s));
            }
        }
        initMajorMinorAlleles();
    }

    public void sortSiteByPhysicalPosition() {
        System.out.println("initPhysicalSort");
        Swapper swapperPos = new Swapper() {

            public void swap(int a, int b) {
                int it;
                it = lociForSite[a];
                lociForSite[a] = lociForSite[b];
                lociForSite[b] = it;
                byte bt;
                bt = strand[a];
                strand[a] = strand[b];
                strand[b] = bt;
                bt = sitePrefix[a];
                sitePrefix[a] = sitePrefix[b];
                sitePrefix[b] = bt;
                for (int t = 0; t < getSequenceCount(); t++) {
                    bt = getBase(t, a);
                    setBase(t, a, getBase(t, b));
                    setBase(t, b, bt);
                }
                it = variableSites[a];
                variableSites[a] = variableSites[b];
                variableSites[b] = it;
            }
        };
        IntComparator compPos = new IntComparator() {

            public int compare(int a, int b) {
                if (lociForSite[a] < lociForSite[b]) {
                    return -1;
                }
                if (lociForSite[a] > lociForSite[b]) {
                    return 1;
                }
                if (variableSites[a] < variableSites[b]) {
                    return -1;
                }
                if (variableSites[a] > variableSites[b]) {
                    return 1;
                }
                if (strand[a] < strand[b]) {
                    return -1;
                }
                if (strand[a] > strand[b]) {
                    return 1;
                }
                return 0;
            }
        };
        System.out.println("Alignment sort begin.");
        GenericSorting.quickSort(0, this.getSiteCount(), compPos, swapperPos);
        System.out.println("Alignment sort end.");
        for (int i = 0; i < siteNumber; i++) {
            if (variableSites[i] == Integer.MAX_VALUE) {
                this.nextFreeSite = i;
                break;
            }
        }
        initMajorMinorAlleles();
    }

    private void initMajorMinorAlleles() {
        majorAlleles = new byte[getSiteCount()];
        minorAlleles = new byte[getSiteCount()];
        for (int i = 0; i < getSiteCount(); i++) {
            byte[] alleles = getAlleles(i);
            majorAlleles[i] = (alleles.length > 0) ? alleles[0] : Alignment.UNKNOWN_ALLELE;
            minorAlleles[i] = (alleles.length > 1) ? alleles[1] : Alignment.UNKNOWN_ALLELE;
        }
    }

    public void clearSite(int site) {
        lociForSite[site] = Integer.MAX_VALUE;
        variableSites[site] = Integer.MAX_VALUE;
        strand[site] = Byte.MAX_VALUE;
        for (int t = 0; t < getSequenceCount(); t++) {
            setBase(t, site, Alignment.UNKNOWN_DIPLOID_ALLELE);
        }
    }

    @Override
    public byte getBase(int taxon, int site) {
        return seq[taxon][site];
    }

    public void setPositionOfSite(int site, int position) {
        variableSites[site] = position;
        if (site >= nextFreeSite) {
            nextFreeSite = site + 1;
        }
    }

    public void setStrandOfSite(int site, byte strand) {
        this.strand[site] = strand;
        if (site >= nextFreeSite) {
            nextFreeSite = site + 1;
        }
    }

    void setLocusIndexOfSite(int site, int locusIndex) {
        lociForSite[site] = locusIndex;
        if (site >= nextFreeSite) {
            nextFreeSite = site + 1;
        }
    }

    public void setLocusOfSite(int site, String locusName) {
        setLocusIndexOfSite(site, locusNameToLociIndex.get(locusName));
        if (site >= nextFreeSite) {
            nextFreeSite = site + 1;
        }
    }

    public byte getSitePrefix(int site) {
        return sitePrefix[site];
    }

    public void setSitePrefix(int site, byte prefix) {
        this.sitePrefix[site] = prefix;
    }

    public int getNextFreeSite() {
        return nextFreeSite;
    }

    @Override
    public String[] getSNPIDs() {
        int sites = getSiteCount();
        String[] SNPids = new String[sites];
        for (int i = 0; i < sites; i++) {
            SNPids[i] = getSNPID(i);
        }
        return SNPids;
    }

    @Override
    public byte getMajorAllele(int site) {
        return majorAlleles[site];
    }

    @Override
    public byte getMinorAllele(int site) {
        return minorAlleles[site];
    }

    @Override
    public String getSNPID(int site) {
        return ((char) getSitePrefix(site)) + getLocus(site).getChromosomeName() + "_" + variableSites[site];
    }

    @Override
    public int getSiteCount() {
        return nextFreeSite;
    }

    @Override
    public int getLocusSiteCount(Locus locus) {
        return getSiteCount();
    }

    @Override
    public int getPositionInLocus(int site) {
        return variableSites[site];
    }

    @Override
    public int getSiteOfPhysicalPosition(int physicalPosition, Locus locus) {
        try {
            return Arrays.binarySearch(variableSites, physicalPosition);
        } catch (Exception e) {
            return physicalPosition;
        }
    }

    @Override
    public byte getPositionType(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public byte[] getPositionTypes() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public Locus getLocus(int site) {
        return myLoci[lociForSite[site]];
    }

    @Override
    public Locus[] getLoci() {
        return myLoci;
    }

    @Override
    public int getNumLoci() {
        return myLoci.length;
    }

    /**
     * Return sorted list of alleles from highest frequency to lowest
     * at given site in alignment. Assumes IUPAC and diploid.
     * Resulting double dimension array
     * holds alleles (actually chars) in result[0].  And the counts
     * are in result[1]. All alleles are counted once.  No support
     * for higher ploidys.
     *
     * @param a alignment
     * @param site site
     * @param includeGAP whether to include GAP
     * @return sorted list of alleles and counts
     */
    @Override
    public int[][] getAllelesSortedByFrequency(int site) {
        byte[] nuc = {'A', 'C', 'G', 'T', '-', '+'};  // note that 'N' is not included as an allele
        int[] nucIndex = new int[Byte.MAX_VALUE];
        for (int i = 0; i < nuc.length; i++) {
            nucIndex[nuc[i]] = i;
        }
        int[] cnt = new int[nuc.length];
        for (int i = 0; i < getSequenceCount(); i++) {
            byte b[] = this.getBaseArray(i, site);
            if (b[0] != Alignment.UNKNOWN_ALLELE) {
                cnt[nucIndex[b[0]]]++;
            }
            if (b[1] != Alignment.UNKNOWN_ALLELE) {
                cnt[nucIndex[b[1]]]++;
            }
        }
        int[][] result = new int[2][nuc.length]; // result[0]=allele; result[1]=count   
        for (int i = 0; i < nuc.length; i++) {
            result[0][i] = nuc[i];
            result[1][i] = cnt[i];
        }
        boolean change = true;
        while (change) { // sort the alleles by descending frequency
            change = false;
            for (int k = 0; k < nuc.length - 1; k++) {
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
    public String[][] getAlleleEncodings() {
        return NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES;
    }

    //
    // Mutable methods
    //
    @Override
    public void setBase(int taxon, int site, byte base) {
        seq[taxon][site] = base;
    }

    public void setBaseRange(int taxon, int startSite, byte[] newBases) {
        throw new UnsupportedOperationException();
    }

    public void addSite(int site) {
        throw new UnsupportedOperationException();
    }

    public void removeSite(int site) {
        throw new UnsupportedOperationException();
    }

    public void addTaxon(Identifier id) {
        throw new UnsupportedOperationException();
    }

    public void removeTaxon(int taxon) {
        throw new UnsupportedOperationException();
    }

    public void clean() {
        throw new UnsupportedOperationException();
    }

    public boolean isDirty() {
        throw new UnsupportedOperationException();
    }
    
    @Override
    public boolean isSBitFriendly() {
        return false;
    }
    
    @Override
    public boolean isTBitFriendly() {
        return false;
    }
}
