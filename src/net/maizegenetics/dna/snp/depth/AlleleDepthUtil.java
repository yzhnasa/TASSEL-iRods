/*
 *  AlleleDepthUtil
 */
package net.maizegenetics.dna.snp.depth;

/**
 *
 * @author Terry Casstevens
 * @author Robert Bukowski
 */
public class AlleleDepthUtil {

    private static double LOG_BASE = 1.0746;  // logbase^128 = 10,000
    private static double LOG_CONV = 1.0 / Math.log(LOG_BASE);

    private AlleleDepthUtil() {
        // utility
    }

    public static byte depthIntToByte(int depth) {
        if (depth > 127) {
            depth = (int) (-LOG_CONV * Math.log(depth));
            if (depth < -128) {
                depth = -127;
            }
        }

        return (byte) depth;
    }

    public static int depthByteToInt(byte depth) {
        if (depth < 128) {
            return depth;
        } else {
            // TODO: Need calculation for higher depths.
            return 127;
        }
    }

}
