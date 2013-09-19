package net.maizegenetics.util;

import ch.systemsx.cisd.hdf5.IHDF5Writer;

import java.util.Arrays;

/**
 * Created with IntelliJ IDEA.
 * User: edbuckler
 * Date: 9/19/13
 * Time: 4:28 PM
 * To change this template use File | Settings | File Templates.
 */
public final class HDF5Utils {
    private HDF5Utils() {}

    /**
     * Needs to go the a JHDF5 Utils.
     */
    public static void writeHDF5EntireArray(String objectPath, IHDF5Writer myWriter, int objMaxLength, int blockSize, Object val) {
        int blocks=((objMaxLength-1)/blockSize)+1;
        for (int block = 0; block < blocks; block++) {
            int startPos=block*blockSize;
            int length=((objMaxLength-startPos)>blockSize)?blockSize:objMaxLength-startPos;
            if(val instanceof byte[][]) {
                byte[][] oval=(byte[][])val;
                byte[][] sval=new byte[oval.length][length];
                for (int j = 0; j < oval.length; j++) {
                    sval[j]= Arrays.copyOfRange(oval[j], startPos, startPos + length);
                }
                writeHDF5Block(objectPath,  myWriter, blockSize, block, sval);
            }  else if(val instanceof int[][]) {
                int[][] oval=(int[][])val;
                int[][] sval=new int[oval.length][length];
                for (int j = 0; j < oval.length; j++) {
                    sval[j]=Arrays.copyOfRange(oval[j], startPos, startPos+length);
                }
                writeHDF5Block(objectPath,  myWriter, blockSize, block, sval);
            } else if(val instanceof byte[]) {
                writeHDF5Block(objectPath,  myWriter, blockSize, block, Arrays.copyOfRange((byte[])val, startPos, startPos+length));
            } else if(val instanceof float[]) {
                writeHDF5Block(objectPath,  myWriter, blockSize, block, Arrays.copyOfRange((float[])val, startPos, startPos+length));
            } else if(val instanceof int[]) {
                writeHDF5Block(objectPath,  myWriter, blockSize, block, Arrays.copyOfRange((int[])val, startPos, startPos+length));
            } else if(val instanceof String[]) {
                writeHDF5Block(objectPath,  myWriter, blockSize, block, Arrays.copyOfRange((String[])val, startPos, startPos+length));
            }
        }
    }

    /**
     * Needs to go the a JHDF5 Utils.
     */
    public static void writeHDF5Block(String objectPath, IHDF5Writer myWriter, int blockSize, int block, Object val) {
        int startPos=block*blockSize;
        if(val instanceof byte[]) {
            byte[] fval=(byte[])val;
            myWriter.writeByteArrayBlockWithOffset(objectPath, fval, fval.length, (long)startPos);
        } else if(val instanceof float[]) {
            float[] fval=(float[])val;
            myWriter.writeFloatArrayBlockWithOffset(objectPath, fval, fval.length, (long)startPos);
        } else if(val instanceof int[]) {
            int[] fval=(int[])val;
            myWriter.writeIntArrayBlockWithOffset(objectPath, fval, fval.length, (long)startPos);
        } else if(val instanceof int[][]) {
            int[][] ival=(int[][])val;
            myWriter.writeIntMatrixBlockWithOffset(objectPath, ival, ival.length, ival[0].length, 0l, (long)startPos);
        } else if(val instanceof byte[][]) {
            byte[][] bval=(byte[][])val;
            myWriter.writeByteMatrixBlockWithOffset(objectPath, bval, bval.length, bval[0].length, 0l, (long)startPos);
        } else if(val instanceof String[]) {
            String[] sval=(String[])val;
            myWriter.writeStringArrayBlockWithOffset(objectPath, sval, sval.length, (long)startPos);
        }
    }

}
