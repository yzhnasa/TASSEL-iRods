/*
 * ReadBarcodeResult
 */
package net.maizegenetics.gbs.homology;

import net.maizegenetics.gbs.util.BaseEncoder;

/**
 * Container class for returning the results of parsed barcoded sequencing read.
 * 
 * @author Fei Lu
 */
public class PEReadBarcodeResult {
    long[] readF;
    short lengthF;
    long[] readB;
    short lengthB;
    String taxonName;

    public PEReadBarcodeResult(ShortReadBarcodeResult rF, ShortReadBarcodeResult rB) {
        this.taxonName = rF.taxonName;
        this.readF = rF.read;
        this.lengthF = rF.length;
        this.readB = rB.read;
        this.lengthB = rB.length;
        this.orderReadFReadB();
    }

    public short getLengthF() {
        return lengthF;
    }

    public short getLengthB() {
        return lengthB;
    }
    
    public long[] getReadF() {
        return readF;
    }
    
    public long[] getReadB() {
        return readB;
    }

    public String getTaxonName() {
        return taxonName;
    }
    
    public int getTagLengthInLong () {
        return readF.length;
    }
    
    private boolean orderReadFReadB () {
        if (this.compareReadFReadB() == 1) {
            this.switchReadFReadB();
            return true;
        }
        return false;
    }
    
    private void switchReadFReadB () {
        for (int i = 0; i < readF.length; i++) {
            long temp = readF[i];
            readF[i] = readB[i];
            readB[i] = temp;  
        }
        short tem = lengthF;
        lengthF = lengthB;
        lengthB = tem;
    }
    
    private int compareReadFReadB () {
        int tagLengthInLong = readF.length;
        for (int i = 0; i < tagLengthInLong; i++) {
            if (readF[i] < readB[i]) return -1;
            if (readF[i] > readB[i]) return 1;
        }
        return 0;
    }
}
