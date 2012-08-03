package net.maizegenetics.gbs.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;

import java.util.Arrays;

/**
 * User: ed
 * Date: Jan 22, 2007
 * Time: 9:43:05 PM
 */
public class BaseEncoder {
    //todo count the fraction of MAGIs hit

    public static final int chunkSize = 32;
    public static final int chunkSizeForInt = 16;
    public static final char[] bases = {'A', 'C', 'G', 'T'};

    private BaseEncoder() {
    }

    public static long[] readTag(File inputFile, int skipLines, int maxTags) {
        int c = 0;
        long[] tempTags = new long[maxTags];
        String temp, tagNumText;
        try {
            FileReader fr = new FileReader(inputFile);
            BufferedReader br = new BufferedReader(fr);
            for (int i = 0; i < skipLines; i++) {
                br.readLine();
            }
            while (br.ready()) {
                temp = br.readLine();
                if (temp.indexOf(' ') < 0) {
                    tempTags[c] = Long.parseLong(temp);
                } else {
                    tagNumText = temp.substring(0, temp.indexOf(' '));
                    tempTags[c] = Long.parseLong(tagNumText);
                }
                c++;
                if (c % 100000 < 1) {
                    System.out.println("read number = " + c);
                }
            }
            fr.close();
            System.out.println("read number = " + c);
        } catch (Exception e) {
            System.out.println("Catch c=" + c + " e=" + e);
        }
        return Arrays.copyOfRange(tempTags, 0, c);
    }

    public static long[] uniqueReads(long[] allSeqs) {
        int uniqueCnt = 1, c = 0;
        for (int i = 0; i < allSeqs.length - 1; i++) {
            if (allSeqs[i] != allSeqs[i + 1]) {
                uniqueCnt++;
            }
        }
        long[] uniSeqs = new long[uniqueCnt];
        // System.out.println("uniqueCnt = " + uniqueCnt);
        for (int i = 1; i < allSeqs.length - 1; i++) {
            // System.out.println("c = " + c + " i="+i);
            if (allSeqs[i] != allSeqs[i + 1]) {
                uniSeqs[c] = allSeqs[i];
                c++;
            }
        }
        return uniSeqs;
    }

    public static char getCharBase(byte base) {
        char c = 'N';
        switch (base) {
            case 0:
                return 'A';  //00
            case 1:
                return 'C';  //01
            case 2:
                return 'G';  //10
            case 3:
                return 'T';  //11
        }
        return c;
    }

    public static long getLongFromSeq(String seq) {
        int seqLength = seq.length();
        long v = 0;
        for (int i = 0; i < seqLength; i++) {
            switch (seq.charAt(i)) {
                case 'A':
                case 'a':
                    v = v << 2;
                    break;
                case 'C':
                case 'c':
                    v = (v << 2) + (byte) 1;
                    break;
                case 'G':
                case 'g':
                    v = (v << 2) + (byte) 2;
                    break;
                case 'T':
                case 't':
                    v = (v << 2) + (byte) 3;
                    break;
                default:
                    return -1;
            }
        }
        if (seqLength == chunkSize) {
            return v;
        }
        if (seqLength > chunkSize) {
            return -1;
        }
        v = (v << (2 * (chunkSize - seqLength))); //if shorter fill with AAAA
        return v;
    }

    /**
     * @param seq A String containing a DNA sequence.
     * @return result A array of Long ints containing the binary representation of the sequence.
     * null if sequence length is not a multiple of BaseEncoder.chunksize.
     */
    public static long[] getLongArrayFromSeq(String seq) {
        if (seq.length() % chunkSize != 0) {
            return null;
        }
        long[] result = new long[seq.length() / chunkSize];
        for (int i = 0; i < result.length; i++) {
            result[i] = getLongFromSeq(seq.substring(i * chunkSize, (i + 1) * chunkSize));
        }
        return result;
    }

    //polyA is used represent unknown, but reverse complement will change it to polyT
    //which does not mean the same
    //sometimes it is best to reverseComplement by text below
    public static long getReverseComplement(long seq, byte len) {
        // if(seq==-1) return -1;
        long rev = 0;
        // byte b=0;
        long mask = 3;
        seq = ~seq;
        for (int i = 0; i < len; i++) {
            rev = (rev << 2) + (seq & mask);
            seq = seq >> 2;
            // System.out.println("v = " + v);
        }
        return rev;
    }

    public static long getReverseComplement(long seq) {
        return getReverseComplement(seq, (byte) chunkSize);
    }

    public static long[] getReverseComplement(long[] seq) {
        long[] rev = new long[seq.length];
        for (int i = 0; i < rev.length; i++) {
            rev[i] = getReverseComplement(seq[seq.length - i - 1], (byte) chunkSize);
        }
        return rev;
    }

    public static String getReverseComplement(String seq) {
        StringBuilder sb = new StringBuilder(seq.length());
        for (int i = seq.length() - 1; i >= 0; i--) {
            sb.append(getComplementBase(seq.charAt(i)));
        }
        return sb.toString();
    }

    public static char getComplementBase(char base) {
        switch (base) {
            case 'A':
                return 'T';
            case 'C':
                return 'G';
            case 'G':
                return 'C';
            case 'T':
                return 'A';
        }
        return 'N';
    }

    public static byte[] getByteSeqFromLong(long val) {
        byte[] b = new byte[chunkSize];
        long mask = 3;
        for (int i = 0; i < chunkSize; i++) {
            b[chunkSize - i - 1] = (byte) (val & mask);
            val = val >> 2;
        }
        return b;
    }

    public static byte[] getByteSeqFromLong(long[] valA) {
        byte[] b = new byte[chunkSize * valA.length];
        long mask = 3;
        long val;
        for (int j = 0; j < valA.length; j++) {
            val = valA[j];
            for (int i = 0; i < chunkSize; i++) {
                b[(j * chunkSize) + chunkSize - i - 1] = (byte) (val & mask);
                val = val >> 2;
            }
        }
        return b;
    }

    public static long getLongSeqFromByteArray(byte[] b) {
        //the byte array must be in 0-3 coding for A, C, G, T
        long v = 0;
        if (b.length != chunkSize) {
            return -1;
        }
        for (int i = 0; i < b.length; i++) {
            v = (v << 2) + b[i];
        }
        return v;
    }

    public static String getSequenceFromLong(long val, byte len) {
        StringBuilder seq = new StringBuilder(chunkSize + 4);
        long mask = 3;
        for (int i = 0; i < len; i++) {
            byte base = (byte) (val & mask);
            seq.insert(0, bases[base]);
            val = val >> 2;
        }
        return seq.toString();
    }

    public static String getSequenceFromLong(long[] val) {
        StringBuilder seq = new StringBuilder();
        for (long v : val) {
            seq.append(getSequenceFromLong(v));
        }
        return seq.toString();
    }

    public static int[] getIntFromLong(long val) {
        int[] ival = new int[2];
        ival[0] = (int) (val >> chunkSize);
        ival[1] = (int) (val);
        return ival;
    }

    public static String getSequenceFromInt(int val) {
        StringBuilder seq = new StringBuilder(chunkSizeForInt + 1);
        long mask = 3;
        for (int i = 0; i < chunkSizeForInt; i++) {
            byte base = (byte) (val & mask);
            seq.insert(0, bases[base]);
            val = val >> 2;
        }
        return seq.toString();
    }

    public static String[] createDegenerateTags(String site, boolean onlyMethylMutations) {
        String[] ps = new String[site.length() * 4];
        StringBuilder sb;
        int numDegenTags = 0;

        if (onlyMethylMutations) {
            ps[numDegenTags++] = site;
            for (int i = 0; i < site.length(); i++) {
                if (site.charAt(i) == 'C') {
                    sb = new StringBuilder(site);
                    ps[numDegenTags++] = sb.replace(i, i + 1, "T").toString();
                }
                if (site.charAt(i) == 'G') {
                    sb = new StringBuilder(site);
                    ps[numDegenTags++] = sb.replace(i, i + 1, "A").toString();
                }
            }
        } else {
            ps[numDegenTags++] = site;
            for (int i = 0; i < site.length(); i++) {
                for (int j = 0; j < bases.length; j++) {
                    if (site.charAt(i) != bases[j]) {
                        sb = new StringBuilder(site);
                        ps[numDegenTags++] = sb.replace(i, i + 1, "" + bases[j]).toString();
                    }
                }
            }
        }
        return Arrays.copyOf(ps, numDegenTags);
    }

    public static int getAvgQuality(String quality) {
        int sum = 0;
        for (int i = 0; i < quality.length(); i++) {
            sum += (int) quality.charAt(i) - 64;
        }
        int avg = sum / quality.length();
        if ((avg < 50) && (avg > -1)) {
            return avg;
        }
        return -1;
    }

    public static int getFirstLowQualityPos(String quality, int minQual) {
        int qualInt = 0;
        for (int i = 0; i < quality.length(); i++) {
            qualInt = (int) quality.charAt(i) - 64;
            if (qualInt < minQual) {
                return i;
            }
        }
        return quality.length();
    }

    public static String getSequenceFromLong(long val) {
        return getSequenceFromLong(val, (byte) chunkSize);
    }

    public static byte seqDifferences(long seq1, long seq2, int maxDivergence) {
        long mask = 3;
        byte cnt = 0;
        long diff = seq1 ^ seq2;
        for (int x = 0; x < chunkSize && cnt <= maxDivergence; x++) {
            if ((diff & mask) > 0) {
                cnt++;
            }
            diff = diff >> 2;
            // System.out.println("v = " + v);
        }
        if (cnt > maxDivergence) {
            cnt = (byte) chunkSize;
        }
        // if(x<(chunkSize-1)) cnt=(byte)chunkSize;  //if didn't get to the end of the sequence set to maximum
        return cnt;
    }

    public static byte seqDifferences(long seq1, long seq2) {
        long mask = 3;
        byte cnt = 0;
        long diff = seq1 ^ seq2;
        for (int x = 0; x < chunkSize; x++) {
            if ((diff & mask) > 0) {
                cnt++;
            }
            diff = diff >> 2;
            // System.out.println("v = " + v);
        }
        return cnt;
    }

    public static byte seqDifferencesForSubset(long seq1, long seq2, int lengthOfComp, int maxDivergence) {
        long mask = 3;
        byte cnt = 0;
        long diff = seq1 ^ seq2;
        diff = diff >> (2 * (chunkSize - lengthOfComp));  //shift to 5' end of sequence
        for (int x = 0; x < lengthOfComp && cnt < maxDivergence; x++) {
            if ((diff & mask) > 0) {
                cnt++;
            }
            diff = diff >> 2;
        }
        return cnt;
    }

    public static byte seqDifferencesWithGapWrong(long seq1, long seq2) {
        //allow a single gap, identify gap by two adjacent mismatches
        long mask = 3;
        byte cnt = 0;
        long diff = seq1 ^ seq2;
        for (int x = 0; x < chunkSize; x++) {
            if ((diff & mask) > 0) {
                cnt++;
            }
            diff = diff >> 2;
            // System.out.println("v = " + v);
        }
        return cnt;
    }

    public static String removePolyAFromEnd(String s) {
        int index = s.length() - 1;
        while (s.charAt(index) == 'A') {
            index--;
            if (index < 1) {
                return null;
            }
        }
        return s.substring(0, index + 1);
    }
}
