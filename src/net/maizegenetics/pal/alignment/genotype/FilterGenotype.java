/*
 *  FilterGenotype
 */
package net.maizegenetics.pal.alignment.genotype;

import net.maizegenetics.pal.alignment.Alignment;

/**
 *
 * @author Terry Casstevens
 */
class FilterGenotype extends AbstractGenotype {

    private final Genotype myBaseGenotype;
    private final boolean myIsTaxaFilter;
    private final boolean myIsSiteFilter;
    private final boolean myIsSiteFilterByRange;
    private final int[] myTaxaRedirect;
    private final int[] mySiteRedirect;
    private final int myRangeStart;
    private final int myRangeEnd;

    FilterGenotype(Genotype genotype, int numTaxa, int[] taxaRedirect, int numSites, int rangeStart, int rangeEnd) {

        super(numTaxa, numSites, genotype.isPhased(), null, genotype.getMaxNumAlleles());
        myBaseGenotype = genotype;

        myTaxaRedirect = taxaRedirect;
        if (taxaRedirect == null) {
            myIsTaxaFilter = false;
        } else {
            myIsTaxaFilter = true;
        }

        myIsSiteFilter = false;
        mySiteRedirect = null;

        myRangeStart = rangeStart;
        myRangeEnd = rangeEnd;
        if ((rangeStart == -1) || (rangeEnd == -1)) {
            myIsSiteFilterByRange = false;
        } else {
            myIsSiteFilterByRange = true;
        }

    }

    FilterGenotype(Genotype genotype, int numTaxa, int[] taxaRedirect, int numSites, int[] siteRedirect) {

        super(numTaxa, numSites, genotype.isPhased(), null, genotype.getMaxNumAlleles());
        myBaseGenotype = genotype;

        myTaxaRedirect = taxaRedirect;
        if (taxaRedirect == null) {
            myIsTaxaFilter = false;
        } else {
            myIsTaxaFilter = true;
        }

        myRangeStart = -1;
        myRangeEnd = -1;
        myIsSiteFilterByRange = false;

        mySiteRedirect = siteRedirect;
        if (siteRedirect == null) {
            myIsSiteFilter = false;
        } else {
            myIsSiteFilter = true;
        }

    }

    @Override
    public byte getBase(int taxon, int site) {
        int taxaIndex = translateTaxon(taxon);
        if (taxaIndex == -1) {
            return Alignment.UNKNOWN_ALLELE;
        } else {
            return myBaseGenotype.getBase(taxaIndex, translateSite(site));
        }
    }

    @Override
    public void transposeData(boolean siteInnerLoop) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String[][] getAlleleEncodings() {
        String[][] encodings = myBaseGenotype.getAlleleEncodings();
        if (encodings.length == 1) {
            return encodings;
        } else if ((myIsSiteFilterByRange) || (myIsSiteFilter)) {
            int numSites = getSiteCount();
            String[][] result = new String[numSites][];
            for (int i = 0; i < numSites; i++) {
                result[i] = getAlleleEncodings(i);
            }
            return result;
        } else {
            return encodings;
        }
    }

    @Override
    public String[] getAlleleEncodings(int site) {
        return myBaseGenotype.getAlleleEncodings(translateSite(site));
    }

    @Override
    public int getMaxNumAlleles() {
        return myBaseGenotype.getMaxNumAlleles();
    }

    public int translateSite(int site) {

        if (myIsSiteFilterByRange) {
            return site + myRangeStart;
        } else if (myIsSiteFilter) {
            return mySiteRedirect[site];
        } else {
            return site;
        }

    }

    public int translateTaxon(int taxon) {

        if (myIsTaxaFilter) {
            return myTaxaRedirect[taxon];
        } else {
            return taxon;
        }

    }
}
