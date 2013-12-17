/*
 *  AlleleDepthUtil
 */
package net.maizegenetics.dna.snp.depth;

/**
 *
 * @author Terry Casstevens
 */
public class AlleleDepthUtil {

    private AlleleDepthUtil() {
        // utility
    }

    public static byte depthIntToByte(int depth) {
        if (depth < 128) {
            return (byte) depth;
        } else {
            // TODO: Need calculation for higher depths
            return (byte) 127;
        }
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
