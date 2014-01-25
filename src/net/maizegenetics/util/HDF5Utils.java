package net.maizegenetics.util;

import ch.systemsx.cisd.base.mdarray.MDArray;
import ch.systemsx.cisd.hdf5.IHDF5Writer;

import java.util.Arrays;

/**
 * @author Ed Buckler
 * @author Terry Casstevens
 */
public final class HDF5Utils {

    private HDF5Utils() {
    }

    public static void writeHDF5GenotypesMaxNumAlleles(IHDF5Writer h5w, int maxNumAlleles) {
        h5w.setIntAttribute(HDF5Constants.GENOTYPES_ATTRIBUTES_PATH, HDF5Constants.GENOTYPES_MAX_NUM_ALLELES, maxNumAlleles);
    }

    public static void writeHDF5GenotypesRetainRareAlleles(IHDF5Writer h5w, boolean retainRareAlleles) {
        h5w.setBooleanAttribute(HDF5Constants.GENOTYPES_ATTRIBUTES_PATH, HDF5Constants.GENOTYPES_RETAIN_RARE_ALLELES, retainRareAlleles);
    }

    public static void writeHDF5GenotypesNumTaxa(IHDF5Writer h5w, int numTaxa) {
        h5w.setIntAttribute(HDF5Constants.GENOTYPES_ATTRIBUTES_PATH, HDF5Constants.GENOTYPES_MAX_NUM_ALLELES, numTaxa);
    }

    public static void writeHDF5GenotypesScoreType(IHDF5Writer h5w, String scoreType) {
        h5w.setStringAttribute(HDF5Constants.GENOTYPES_ATTRIBUTES_PATH, HDF5Constants.GENOTYPES_MAX_NUM_ALLELES, scoreType);
    }

    public static void writeHDF5GenotypesAlleleStates(IHDF5Writer h5w, String[][] aEncodings) {
        int numEncodings = aEncodings.length;
        int numStates = aEncodings[0].length;
        MDArray<String> alleleEncodings = new MDArray<>(String.class, new int[]{numEncodings, numStates});
        for (int s = 0; s < numEncodings; s++) {
            for (int x = 0; x < numStates; x++) {
                alleleEncodings.set(aEncodings[s][x], s, x);
            }
        }
        h5w.createStringMDArray(HDF5Constants.GENOTYPES_ALLELE_STATES, 100, new int[]{numEncodings, numStates});
        h5w.writeStringMDArray(HDF5Constants.GENOTYPES_ALLELE_STATES, alleleEncodings);
    }

    public static void writeHDF5GenotypesCalls(IHDF5Writer h5w, String taxon, byte[] calls) {
        writeHDF5EntireArray(HDF5Constants.getGenotypesCallsPath(taxon), h5w, calls.length, 1 << 16, calls);
    }

    /**
     *
     * @param objectPath
     * @param myWriter
     * @param objMaxLength
     * @param blockSize
     * @param val
     */
    public static void writeHDF5EntireArray(String objectPath, IHDF5Writer myWriter, int objMaxLength, int blockSize, Object val) {
        int blocks = ((objMaxLength - 1) / blockSize) + 1;
        for (int block = 0; block < blocks; block++) {
            int startPos = block * blockSize;
            int length = ((objMaxLength - startPos) > blockSize) ? blockSize : objMaxLength - startPos;
            if (val instanceof byte[][]) {
                byte[][] oval = (byte[][]) val;
                byte[][] sval = new byte[oval.length][length];
                for (int j = 0; j < oval.length; j++) {
                    sval[j] = Arrays.copyOfRange(oval[j], startPos, startPos + length);
                }
                writeHDF5Block(objectPath, myWriter, blockSize, block, sval);
            } else if (val instanceof int[][]) {
                int[][] oval = (int[][]) val;
                int[][] sval = new int[oval.length][length];
                for (int j = 0; j < oval.length; j++) {
                    sval[j] = Arrays.copyOfRange(oval[j], startPos, startPos + length);
                }
                writeHDF5Block(objectPath, myWriter, blockSize, block, sval);
            } else if (val instanceof byte[]) {
                writeHDF5Block(objectPath, myWriter, blockSize, block, Arrays.copyOfRange((byte[]) val, startPos, startPos + length));
            } else if (val instanceof float[]) {
                writeHDF5Block(objectPath, myWriter, blockSize, block, Arrays.copyOfRange((float[]) val, startPos, startPos + length));
            } else if (val instanceof int[]) {
                writeHDF5Block(objectPath, myWriter, blockSize, block, Arrays.copyOfRange((int[]) val, startPos, startPos + length));
            } else if (val instanceof String[]) {
                writeHDF5Block(objectPath, myWriter, blockSize, block, Arrays.copyOfRange((String[]) val, startPos, startPos + length));
            }
        }
    }

    /**
     *
     * @param objectPath
     * @param myWriter
     * @param blockSize
     * @param block
     * @param val
     */
    public static void writeHDF5Block(String objectPath, IHDF5Writer myWriter, int blockSize, int block, Object val) {
        int startPos = block * blockSize;
        if (val instanceof byte[]) {
            byte[] fval = (byte[]) val;
            myWriter.writeByteArrayBlockWithOffset(objectPath, fval, fval.length, (long) startPos);
        } else if (val instanceof float[]) {
            float[] fval = (float[]) val;
            myWriter.writeFloatArrayBlockWithOffset(objectPath, fval, fval.length, (long) startPos);
        } else if (val instanceof int[]) {
            int[] fval = (int[]) val;
            myWriter.writeIntArrayBlockWithOffset(objectPath, fval, fval.length, (long) startPos);
        } else if (val instanceof int[][]) {
            int[][] ival = (int[][]) val;
            myWriter.writeIntMatrixBlockWithOffset(objectPath, ival, ival.length, ival[0].length, 0l, (long) startPos);
        } else if (val instanceof byte[][]) {
            byte[][] bval = (byte[][]) val;
            myWriter.writeByteMatrixBlockWithOffset(objectPath, bval, bval.length, bval[0].length, 0l, (long) startPos);
        } else if (val instanceof String[]) {
            String[] sval = (String[]) val;
            myWriter.writeStringArrayBlockWithOffset(objectPath, sval, sval.length, (long) startPos);
        }
    }

}
