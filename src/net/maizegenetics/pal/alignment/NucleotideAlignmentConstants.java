/*
 * NucleotideAlignmentConstants
 */
package net.maizegenetics.pal.alignment;

import java.util.HashMap;
import java.util.Map;

/**
 *
 * @author terry
 */
public final class NucleotideAlignmentConstants {

    public static byte A_ALLELE = (byte) 0x0;
    public static byte C_ALLELE = (byte) 0x1;
    public static byte G_ALLELE = (byte) 0x2;
    public static byte T_ALLELE = (byte) 0x3;
    public static byte INSERT_ALLELE = (byte) 0x4;
    public static byte GAP_ALLELE = 0x5;
    public static byte GAP_DIPLOID_ALLELE = (byte) 0x55;
    public static byte A_DIPLOID_ALLELE = (byte) 0x00;
    public static byte C_DIPLOID_ALLELE = (byte) 0x11;
    public static byte M_DIPLOID_ALLELE = (byte) 0x01;
    public static String[][] NUCLEOTIDE_ALLELES = new String[][]{{"A", "C", "G", "T", "+", "-",
            Alignment.UNKNOWN_ALLELE_STR, Alignment.UNKNOWN_ALLELE_STR, Alignment.UNKNOWN_ALLELE_STR, Alignment.UNKNOWN_ALLELE_STR, Alignment.UNKNOWN_ALLELE_STR,
            Alignment.UNKNOWN_ALLELE_STR, Alignment.UNKNOWN_ALLELE_STR, Alignment.UNKNOWN_ALLELE_STR, Alignment.RARE_ALLELE_STR, Alignment.UNKNOWN_ALLELE_STR}};
    /**
     * Number of nucleotide states excluding rare and unknown.
     */
    public static final int NUMBER_NUCLEOTIDE_ALLELES = 6;
    private static final Map NUCLEOTIDE_DIPLOID_HASH = new HashMap<String, Byte>();

    static {
        NUCLEOTIDE_DIPLOID_HASH.put("AA", (byte) 0x00);
        NUCLEOTIDE_DIPLOID_HASH.put("AC", (byte) 0x01);
        NUCLEOTIDE_DIPLOID_HASH.put("AG", (byte) 0x02);
        NUCLEOTIDE_DIPLOID_HASH.put("AT", (byte) 0x03);
        NUCLEOTIDE_DIPLOID_HASH.put("A+", (byte) 0x04);
        NUCLEOTIDE_DIPLOID_HASH.put("A-", (byte) 0x05);
        NUCLEOTIDE_DIPLOID_HASH.put("AN", (byte) 0x0F);
        NUCLEOTIDE_DIPLOID_HASH.put("AX", (byte) 0x0F);
        NUCLEOTIDE_DIPLOID_HASH.put("AZ", (byte) 0x0E);

        NUCLEOTIDE_DIPLOID_HASH.put("CA", (byte) 0x10);
        NUCLEOTIDE_DIPLOID_HASH.put("CC", (byte) 0x11);
        NUCLEOTIDE_DIPLOID_HASH.put("CG", (byte) 0x12);
        NUCLEOTIDE_DIPLOID_HASH.put("CT", (byte) 0x13);
        NUCLEOTIDE_DIPLOID_HASH.put("C+", (byte) 0x14);
        NUCLEOTIDE_DIPLOID_HASH.put("C-", (byte) 0x15);
        NUCLEOTIDE_DIPLOID_HASH.put("CN", (byte) 0x1F);
        NUCLEOTIDE_DIPLOID_HASH.put("CX", (byte) 0x1F);
        NUCLEOTIDE_DIPLOID_HASH.put("CZ", (byte) 0x1E);

        NUCLEOTIDE_DIPLOID_HASH.put("GA", (byte) 0x20);
        NUCLEOTIDE_DIPLOID_HASH.put("GC", (byte) 0x21);
        NUCLEOTIDE_DIPLOID_HASH.put("GG", (byte) 0x22);
        NUCLEOTIDE_DIPLOID_HASH.put("GT", (byte) 0x23);
        NUCLEOTIDE_DIPLOID_HASH.put("G+", (byte) 0x24);
        NUCLEOTIDE_DIPLOID_HASH.put("G-", (byte) 0x25);
        NUCLEOTIDE_DIPLOID_HASH.put("GN", (byte) 0x2F);
        NUCLEOTIDE_DIPLOID_HASH.put("GX", (byte) 0x2F);
        NUCLEOTIDE_DIPLOID_HASH.put("GZ", (byte) 0x2E);

        NUCLEOTIDE_DIPLOID_HASH.put("TA", (byte) 0x30);
        NUCLEOTIDE_DIPLOID_HASH.put("TC", (byte) 0x31);
        NUCLEOTIDE_DIPLOID_HASH.put("TG", (byte) 0x32);
        NUCLEOTIDE_DIPLOID_HASH.put("TT", (byte) 0x33);
        NUCLEOTIDE_DIPLOID_HASH.put("T+", (byte) 0x34);
        NUCLEOTIDE_DIPLOID_HASH.put("T-", (byte) 0x35);
        NUCLEOTIDE_DIPLOID_HASH.put("TN", (byte) 0x3F);
        NUCLEOTIDE_DIPLOID_HASH.put("TX", (byte) 0x3F);
        NUCLEOTIDE_DIPLOID_HASH.put("TZ", (byte) 0x3E);

        NUCLEOTIDE_DIPLOID_HASH.put("+A", (byte) 0x40);
        NUCLEOTIDE_DIPLOID_HASH.put("+C", (byte) 0x41);
        NUCLEOTIDE_DIPLOID_HASH.put("+G", (byte) 0x42);
        NUCLEOTIDE_DIPLOID_HASH.put("+T", (byte) 0x43);
        NUCLEOTIDE_DIPLOID_HASH.put("++", (byte) 0x44);
        NUCLEOTIDE_DIPLOID_HASH.put("+-", (byte) 0x45);
        NUCLEOTIDE_DIPLOID_HASH.put("+N", (byte) 0x4F);
        NUCLEOTIDE_DIPLOID_HASH.put("+X", (byte) 0x4F);
        NUCLEOTIDE_DIPLOID_HASH.put("+Z", (byte) 0x4E);

        NUCLEOTIDE_DIPLOID_HASH.put("-A", (byte) 0x50);
        NUCLEOTIDE_DIPLOID_HASH.put("-C", (byte) 0x51);
        NUCLEOTIDE_DIPLOID_HASH.put("-G", (byte) 0x52);
        NUCLEOTIDE_DIPLOID_HASH.put("-T", (byte) 0x53);
        NUCLEOTIDE_DIPLOID_HASH.put("-+", (byte) 0x54);
        NUCLEOTIDE_DIPLOID_HASH.put("--", (byte) 0x55);
        NUCLEOTIDE_DIPLOID_HASH.put("-N", (byte) 0x5F);
        NUCLEOTIDE_DIPLOID_HASH.put("-X", (byte) 0x5F);
        NUCLEOTIDE_DIPLOID_HASH.put("-Z", (byte) 0x5E);

        NUCLEOTIDE_DIPLOID_HASH.put("NA", (byte) 0xF0);
        NUCLEOTIDE_DIPLOID_HASH.put("NC", (byte) 0xF1);
        NUCLEOTIDE_DIPLOID_HASH.put("NG", (byte) 0xF2);
        NUCLEOTIDE_DIPLOID_HASH.put("NT", (byte) 0xF3);
        NUCLEOTIDE_DIPLOID_HASH.put("N+", (byte) 0xF4);
        NUCLEOTIDE_DIPLOID_HASH.put("N-", (byte) 0xF5);
        NUCLEOTIDE_DIPLOID_HASH.put("NN", (byte) 0xFF);
        NUCLEOTIDE_DIPLOID_HASH.put("NX", (byte) 0xFF);
        NUCLEOTIDE_DIPLOID_HASH.put("NZ", (byte) 0xFE);

        NUCLEOTIDE_DIPLOID_HASH.put("XA", (byte) 0xF0);
        NUCLEOTIDE_DIPLOID_HASH.put("XC", (byte) 0xF1);
        NUCLEOTIDE_DIPLOID_HASH.put("XG", (byte) 0xF2);
        NUCLEOTIDE_DIPLOID_HASH.put("XT", (byte) 0xF3);
        NUCLEOTIDE_DIPLOID_HASH.put("X+", (byte) 0xF4);
        NUCLEOTIDE_DIPLOID_HASH.put("X-", (byte) 0xF5);
        NUCLEOTIDE_DIPLOID_HASH.put("XN", (byte) 0xFF);
        NUCLEOTIDE_DIPLOID_HASH.put("XX", (byte) 0xFF);
        NUCLEOTIDE_DIPLOID_HASH.put("XZ", (byte) 0xFE);

        NUCLEOTIDE_DIPLOID_HASH.put("ZA", (byte) 0xE0);
        NUCLEOTIDE_DIPLOID_HASH.put("ZC", (byte) 0xE1);
        NUCLEOTIDE_DIPLOID_HASH.put("ZG", (byte) 0xE2);
        NUCLEOTIDE_DIPLOID_HASH.put("ZT", (byte) 0xE3);
        NUCLEOTIDE_DIPLOID_HASH.put("Z+", (byte) 0xE4);
        NUCLEOTIDE_DIPLOID_HASH.put("Z-", (byte) 0xE5);
        NUCLEOTIDE_DIPLOID_HASH.put("ZN", (byte) 0xEF);
        NUCLEOTIDE_DIPLOID_HASH.put("ZX", (byte) 0xEF);
        NUCLEOTIDE_DIPLOID_HASH.put("ZZ", (byte) 0xEE);

        NUCLEOTIDE_DIPLOID_HASH.put("A", (byte) 0x00); // AA
        NUCLEOTIDE_DIPLOID_HASH.put("C", (byte) 0x11); // CC
        NUCLEOTIDE_DIPLOID_HASH.put("G", (byte) 0x22); // GG
        NUCLEOTIDE_DIPLOID_HASH.put("T", (byte) 0x33); // TT
        NUCLEOTIDE_DIPLOID_HASH.put("+", (byte) 0x44); // ++
        NUCLEOTIDE_DIPLOID_HASH.put("-", (byte) 0x55); // --
        NUCLEOTIDE_DIPLOID_HASH.put("Z", (byte) 0xEE); // ZZ
        NUCLEOTIDE_DIPLOID_HASH.put("N", (byte) 0xFF); // NN
        NUCLEOTIDE_DIPLOID_HASH.put("X", (byte) 0xFF); // NN

        NUCLEOTIDE_DIPLOID_HASH.put("R", (byte) 0x02); // AG
        NUCLEOTIDE_DIPLOID_HASH.put("Y", (byte) 0x13); // CT
        NUCLEOTIDE_DIPLOID_HASH.put("S", (byte) 0x21); // GC
        NUCLEOTIDE_DIPLOID_HASH.put("W", (byte) 0x03); // AT
        NUCLEOTIDE_DIPLOID_HASH.put("K", (byte) 0x23); // GT
        NUCLEOTIDE_DIPLOID_HASH.put("M", (byte) 0x01); // AC
        NUCLEOTIDE_DIPLOID_HASH.put("0", (byte) 0x54); // -+
    }
    private static final Map<Byte, String> NUCLEOTIDE_IUPAC_HASH = new HashMap<Byte, String>();

    static {
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x00, "A"); // AA
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x01, "M"); // AC
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x02, "R"); // AG
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x03, "W"); // AT
        // NUCLEOTIDE_IUPAC_HASH.put((byte) 0x04, "?"); // A+
        // NUCLEOTIDE_IUPAC_HASH.put((byte) 0x05, "?"); // A-
        // NUCLEOTIDE_IUPAC_HASH.put((byte) 0x0E, "?"); // AZ
        // NUCLEOTIDE_IUPAC_HASH.put((byte) 0x0F, "?"); // AN

        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x10, "M"); // CA
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x11, "C"); // CC
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x12, "S"); // CG
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x13, "Y"); // CT
        // NUCLEOTIDE_IUPAC_HASH.put((byte) 0x14, "?"); // C+
        // NUCLEOTIDE_IUPAC_HASH.put((byte) 0x15, "?"); // C-
        // NUCLEOTIDE_IUPAC_HASH.put((byte) 0x1E, "?"); // CZ
        // NUCLEOTIDE_IUPAC_HASH.put((byte) 0x1F, "?"); // CN

        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x20, "R"); // GA
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x21, "S"); // GC
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x22, "G"); // GG
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x23, "K"); // GT
        // NUCLEOTIDE_IUPAC_HASH.put((byte) 0x24, "?"); // G+
        // NUCLEOTIDE_IUPAC_HASH.put((byte) 0x25, "?"); // G-
        // NUCLEOTIDE_IUPAC_HASH.put((byte) 0x2E, "?"); // GZ
        // NUCLEOTIDE_IUPAC_HASH.put((byte) 0x2F, "?"); // GN

        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x30, "W"); // TA
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x31, "Y"); // TC
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x32, "K"); // TG
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x33, "T"); // TT
        // NUCLEOTIDE_IUPAC_HASH.put((byte) 0x34, "?"); // T+
        // NUCLEOTIDE_IUPAC_HASH.put((byte) 0x35, "?"); // T-
        // NUCLEOTIDE_IUPAC_HASH.put((byte) 0x3E, "?"); // TZ
        // NUCLEOTIDE_IUPAC_HASH.put((byte) 0x3F, "?"); // TN

        // NUCLEOTIDE_IUPAC_HASH.put((byte) 0x40, "?"); // +A
        // NUCLEOTIDE_IUPAC_HASH.put((byte) 0x41, "?"); // +C
        // NUCLEOTIDE_IUPAC_HASH.put((byte) 0x42, "?"); // +G
        // NUCLEOTIDE_IUPAC_HASH.put((byte) 0x43, "?"); // +T
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x44, "+"); // ++
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x45, "0"); // +-
        // NUCLEOTIDE_IUPAC_HASH.put((byte) 0x4E, "?"); // +Z
        // NUCLEOTIDE_IUPAC_HASH.put((byte) 0x4F, "?"); // +N

        // NUCLEOTIDE_IUPAC_HASH.put((byte) 0x50, "?"); // -A
        // NUCLEOTIDE_IUPAC_HASH.put((byte) 0x51, "?"); // -C
        // NUCLEOTIDE_IUPAC_HASH.put((byte) 0x52, "?"); // -G
        // NUCLEOTIDE_IUPAC_HASH.put((byte) 0x53, "?"); // -T
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x54, "0"); // -+
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0x55, "-"); // --
        // NUCLEOTIDE_IUPAC_HASH.put((byte) 0x5E, "?"); // -Z
        // NUCLEOTIDE_IUPAC_HASH.put((byte) 0x5F, "?"); // -N

        // NUCLEOTIDE_IUPAC_HASH.put((byte) 0xE0, "?"); // ZA
        // NUCLEOTIDE_IUPAC_HASH.put((byte) 0xE1, "?"); // ZC
        // NUCLEOTIDE_IUPAC_HASH.put((byte) 0xE2, "?"); // ZG
        // NUCLEOTIDE_IUPAC_HASH.put((byte) 0xE3, "?"); // ZT
        // NUCLEOTIDE_IUPAC_HASH.put((byte) 0xE4, "?"); // Z+
        // NUCLEOTIDE_IUPAC_HASH.put((byte) 0xE5, "?"); // Z-
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0xEE, "Z"); // ZZ
        // NUCLEOTIDE_IUPAC_HASH.put((byte) 0xEF, "?"); // ZN

        // NUCLEOTIDE_IUPAC_HASH.put((byte) 0xF0, "?"); // NA
        // NUCLEOTIDE_IUPAC_HASH.put((byte) 0xF1, "?"); // NC
        // NUCLEOTIDE_IUPAC_HASH.put((byte) 0xF2, "?"); // NG
        // NUCLEOTIDE_IUPAC_HASH.put((byte) 0xF3, "?"); // NT
        // NUCLEOTIDE_IUPAC_HASH.put((byte) 0xF4, "?"); // N+
        // NUCLEOTIDE_IUPAC_HASH.put((byte) 0xF5, "?"); // N-
        // NUCLEOTIDE_IUPAC_HASH.put((byte) 0xFE, "?"); // NZ
        NUCLEOTIDE_IUPAC_HASH.put((byte) 0xFF, "N"); // NN

    }

    private NucleotideAlignmentConstants() {
        // do not instantiate
    }

    /**
     * Returns diploid byte value for given nucleotide value. First four bits
     * contain first allele value. And second four bits contain second allele
     * value.
     *
     * @param value
     *
     * @return nucleotide diploid allele byte value
     */
    public static byte getNucleotideDiploidByte(String value) {
        try {
            return ((Byte) NUCLEOTIDE_DIPLOID_HASH.get(value)).byteValue();
        } catch (NullPointerException e) {
            throw new IllegalArgumentException("NucleotideAlignmentConstants: getNucleotideDiploidByte: unknown allele value: " + value);
        }
    }

    public static byte getNucleotideDiploidByte(char value) {
        try {
            return ((Byte) NUCLEOTIDE_DIPLOID_HASH.get(String.valueOf(value))).byteValue();
        } catch (NullPointerException e) {
            throw new IllegalArgumentException("NucleotideAlignmentConstants: getNucleotideDiploidByte: unknown allele value: " + value);
        }
    }

    public static String getNucleotideIUPAC(byte value) {
        try {
            return ((String) NUCLEOTIDE_IUPAC_HASH.get(value));
        } catch (NullPointerException e) {
            throw new IllegalArgumentException("NucleotideAlignmentConstants: getNucleotideIUPAC: unknown allele value: " + value);
        }
    }

    public static byte getNucleotideComplement(byte geno) {

        if (geno == A_ALLELE) {
            return T_ALLELE;
        } else if (geno == T_ALLELE) {
            return A_ALLELE;
        } else if (geno == C_ALLELE) {
            return G_ALLELE;
        } else if (geno == G_ALLELE) {
            return C_ALLELE;
        } else {
            return geno;
        }

    }

    public static byte getNucleotideDiploidComplement(byte diploidAllele) {

        byte first = (byte) ((diploidAllele >>> 4) & 0xf);
        byte second = (byte) (diploidAllele & 0xf);
        first = getNucleotideComplement(first);
        second = getNucleotideComplement(second);
        return (byte) ((first << 4) | second);

    }

    public static boolean isNucleotideEncodings(String[][] alleleStates) {

        boolean isNucleotide = false;
        if (alleleStates.length == 1) {
            isNucleotide = true;
            if (alleleStates[0].length == NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES[0].length) {
                for (int i = 0; i < alleleStates.length; i++) {
                    if (!alleleStates[0][i].equals(NucleotideAlignmentConstants.NUCLEOTIDE_ALLELES[0][i])) {
                        isNucleotide = false;
                    }
                }
            }

        }

        return isNucleotide;

    }
}
