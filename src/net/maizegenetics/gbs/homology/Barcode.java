package net.maizegenetics.gbs.homology;

import java.util.Arrays;

import net.maizegenetics.gbs.util.BaseEncoder;

/**
 * Contains methods and information on a barcode used in encoding tags.
 * 
 * @author edbuckler
 */
public class Barcode implements Comparable<Barcode> {

    String barcodeS;
    String[] overhangS;
    String taxaName, flowcell, lane;
    long[] barOverLong;
    int barOverLength, barLength;

    public Barcode(String barcodeS, String[] overhangSunsort, String taxa, String flowcell, String lane) {
        this.barcodeS = barcodeS;
        Arrays.sort(overhangSunsort);
        this.overhangS = overhangSunsort;
        this.flowcell = flowcell;
        this.lane = lane;
        this.taxaName = taxa;
        barOverLong = new long[overhangS.length];
        for (int i = 0; i < overhangS.length; i++) {
            barOverLong[i] = BaseEncoder.getLongFromSeq(barcodeS + overhangS[i]);
        }
        barOverLength = barcodeS.length() + overhangS[0].length();
        barLength = barcodeS.length();
    }

    public int compareSequence(long queryLong, int maxDivCheck) {
        int div = barOverLength;
        for (long targetLong : barOverLong) {
            int c = BaseEncoder.seqDifferencesForSubset(targetLong, queryLong, barOverLength, maxDivCheck);
            if (c < div) {
                div = c;
            }
        }
        return div;
    }

    @Override
    public int compareTo(Barcode anotherBarcode) {
        if (this.barOverLong[0] < anotherBarcode.barOverLong[0]) {
            return -1;
        }
        if (this.barOverLong[0] > anotherBarcode.barOverLong[0]) {
            return 1;
        }
        return 0;
    }

    public String getTaxaName() {
        return taxaName;
    }
}
