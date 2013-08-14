/*
 *  ByteGenotype
 */
package net.maizegenetics.pal.alignment.genotype;

import net.maizegenetics.pal.alignment.AlignmentNew;
import net.maizegenetics.pal.alignment.AlignmentUtils;
import net.maizegenetics.util.ProgressListener;
import org.apache.log4j.Logger;

/**
 *
 * @author Terry Casstevens
 */
public class ByteGenotype implements Genotype {

    private static final Logger myLogger = Logger.getLogger(ByteGenotype.class);
    private final byte[][] myGenotype;
    private final int myTaxaCount;
    private final int mySiteCount;

    ByteGenotype(byte[][] genotype) {
        myTaxaCount = genotype.length;
        mySiteCount = genotype[0].length;
        myGenotype = new byte[myTaxaCount][mySiteCount];
        for (int t = 0; t < myTaxaCount; t++) {
            System.arraycopy(genotype[t], 0, myGenotype[t], 0, mySiteCount);
        }
    }

    @Override
    public byte getBase(int taxon, int site) {
        return myGenotype[taxon][site];
    }

    @Override
    public byte[] getBaseArray(int taxon, int site) {
        return AlignmentUtils.getDiploidValues(myGenotype[taxon][site]);
    }

    @Override
    public byte[] getBaseRange(int taxon, int startSite, int endSite) {
        byte[] result = new byte[endSite - startSite];
        for (int i = startSite; i < endSite; i++) {
            result[i - startSite] = getBase(taxon, i);
        }
        return result;
    }

    @Override
    public byte[] getBaseRow(int taxon) {
        byte[] result = new byte[mySiteCount];
        for (int i = 0; i < mySiteCount; i++) {
            result[i] = getBase(taxon, i);
        }
        return result;
    }

    @Override
    public String getBaseAsString(int taxon, int site) {
        String[][] alleleStates = getAlleleEncodings();
        byte[] temp = getBaseArray(taxon, site);
        return alleleStates[0][temp[0]] + ":" + alleleStates[0][temp[1]];
    }

    @Override
    public String getBaseAsStringRange(int taxon, int startSite, int endSite) {
        StringBuilder builder = new StringBuilder();
        for (int i = startSite; i < endSite; i++) {
            if (i != startSite) {
                builder.append(";");
            }
            builder.append(getBaseAsString(taxon, i));
        }
        return builder.toString();
    }

    @Override
    public String getBaseAsStringRow(int taxon) {
        return getBaseAsStringRange(taxon, 0, mySiteCount);
    }

    @Override
    public String[] getBaseAsStringArray(int taxon, int site) {
        String[][] alleleStates = getAlleleEncodings();
        byte[] temp = getBaseArray(taxon, site);
        return new String[]{alleleStates[0][temp[0]], alleleStates[0][temp[1]]};
    }

    @Override
    public boolean isHeterozygous(int taxon, int site) {
        byte[] values = getBaseArray(taxon, site);
        if (values[0] == values[1]) {
            return false;
        } else {
            return true;
        }
    }

    @Override
    public int getHeterozygousCount(int site) {
        int result = 0;
        for (int i = 0, n = myTaxaCount; i < n; i++) {
            if (isHeterozygous(i, site)) {
                result++;
            }
        }
        return result;
    }

    @Override
    public boolean isPolymorphic(int site) {

        byte first = AlignmentNew.UNKNOWN_ALLELE;
        for (int i = 0, n = myTaxaCount; i < n; i++) {
            byte[] current = getBaseArray(i, site);
            if (current[0] != AlignmentNew.UNKNOWN_ALLELE) {
                if (first == AlignmentNew.UNKNOWN_ALLELE) {
                    first = current[0];
                } else if (first != current[0]) {
                    return true;
                }
            }
            if (current[1] != AlignmentNew.UNKNOWN_ALLELE) {
                if (first == AlignmentNew.UNKNOWN_ALLELE) {
                    first = current[1];
                } else if (first != current[1]) {
                    return true;
                }
            }
        }

        return false;

    }

    @Override
    public boolean isAllPolymorphic() {

        for (int i = 0, n = mySiteCount; i < n; i++) {
            if (!isPolymorphic(i)) {
                return false;
            }
        }

        return true;

    }

    @Override
    public boolean isPhased() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public boolean retainsRareAlleles() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String[][] getAlleleEncodings() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String[] getAlleleEncodings(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String getBaseAsString(int site, byte value) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String getDiploidAsString(int site, byte value) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int getMaxNumAlleles() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int getTotalNumAlleles() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int getTotalGametesNotMissing(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int getTotalNotMissing(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int getMinorAlleleCount(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int getMajorAlleleCount(int site) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public Object[][] getDiploidCounts() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public Object[][] getMajorMinorCounts() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int getTotalGametesNotMissingForTaxon(int taxon) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int getHeterozygousCountForTaxon(int taxon) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public int getTotalNotMissingForTaxon(int taxon) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public boolean isSBitFriendly() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public boolean isTBitFriendly() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public void optimizeForTaxa(ProgressListener listener) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public void optimizeForSites(ProgressListener listener) {
        throw new UnsupportedOperationException("Not supported yet.");
    }
}
