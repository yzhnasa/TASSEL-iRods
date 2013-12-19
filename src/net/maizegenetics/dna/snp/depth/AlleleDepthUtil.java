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

    private static final double LOG_BASE = 1.0746;  // LOG_BASE^128 = 10,356
    private static final double LOG_CONV = 1.0 / Math.log(LOG_BASE);
    private static final double DECODE_LOG_CONV = Math.log(LOG_BASE);
    private static final int[] BYTE_TO_INT = new int[256];

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
        if (depth > 127) {
            depth = (int) (-LOG_CONV * Math.log(depth));
            if (depth < -128) {
                depth = -128;
            }
        }

        return (byte) depth;
    }

    public static int depthByteToInt(byte depth) {
        return BYTE_TO_INT[depth & 0xFF];
    }

    private static int decode(byte depth) {
        if (depth >= 0) {
            return depth;
        } else {
            return (int) (Math.exp(-DECODE_LOG_CONV * (depth - 0.5)));
        }
    }

}
