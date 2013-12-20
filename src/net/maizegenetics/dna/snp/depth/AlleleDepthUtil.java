/*
 *  AlleleDepthUtil
 */
package net.maizegenetics.dna.snp.depth;

import java.util.Arrays;

/**
 *
 * @author Terry Casstevens
 * @author Robert Bukowski
 */
public class AlleleDepthUtil {

    private static final double LOG_BASE = 1.0746;  // LOG_BASE^128 = 10,482
    private static final double R_LOG_CONV = 1.0 / Math.log(LOG_BASE);
    private static final double LOG_CONV = 1.0 / R_LOG_CONV;
    private static final int[] BYTE_TO_INT = new int[256];
    private static final int MAX_ACC_DEPTH = 182;
    private static final int MIN_ACC_BYTE = 127 - MAX_ACC_DEPTH;
    private static final int OFFSET = 126;
    private static final double ADJ = 0.5;

    static {
        Arrays.fill(BYTE_TO_INT, -1);
        for (int i = 0; i < 256; i++) {
            BYTE_TO_INT[i] = decode((byte) i);
        }
    }

    private AlleleDepthUtil() {
        // utility
    }

    public static byte depthIntToByte(int depth) {
        byte bdepth;
        int itd;

        if (depth <= 127) {
            itd = depth;
        } else if (depth <= MAX_ACC_DEPTH) {
            itd = 127 - depth;
        } else {
            itd = (int) (-R_LOG_CONV * Math.log(depth - OFFSET));
            if (itd < -128) {
                itd = -128;
            }
        }
        bdepth = (byte) itd;

        return bdepth;
    }

    public static int depthByteToInt(byte depth) {
        return BYTE_TO_INT[depth & 0xFF];
    }

    private static int decode(byte bdepth) {
        int depth;
        if (bdepth >= 0) {
            depth = bdepth;
        } else if (bdepth >= MIN_ACC_BYTE) {
            depth = 127 - bdepth;
        } else {
            depth = OFFSET + (int) (Math.exp(-LOG_CONV * (bdepth - ADJ)));
        }

        return depth;
    }

}
