package net.maizegenetics.util;


import java.io.*;
import java.security.MessageDigest;

/**
 * User: dkroon
 * Date: 4/9/13
 * Time: 11:47 AM
 * Adapted from http://www.rgagnon.com/javadetails/java-0416.html  ("Real's HowTo")
 */
public class CheckSum {

    /**
     *
     * @param filename
     * @param protocol  Current natively supported protocols should include MD5, MD2, SHA, SHA1, SHA-256, SHA-384, SHA-512
     * @return
     * @throws Exception
     */
    public static byte[] createChecksum(File filename, String protocol) throws Exception
    {
        InputStream fis =  new FileInputStream(filename);

        byte[] buffer = new byte[1024];
        MessageDigest complete = MessageDigest.getInstance(protocol);
        int numRead;
        do {
            numRead = fis.read(buffer);
            if (numRead > 0) {
                complete.update(buffer, 0, numRead);
            }
        } while (numRead != -1);
        fis.close();
        return complete.digest();
    }


    private static String getMD5ChecksumByFile(File aFile)  {

        String result = "";
        try{
            byte[] b = createChecksum(aFile, "MD5");

            for (int i=0; i < b.length; i++) {
                result += Integer.toString( ( b[i] & 0xff ) + 0x100, 16).substring( 1 );
            }
        }catch(Exception e){ /* ignore */ }
        return result;
    }


    public static String getMD5Checksum(File aFile) {

        if(aFile.isDirectory()){
            StringBuffer sb = new StringBuffer();
                File[] dirFiles = aFile.listFiles();
                for(int i=0; i<dirFiles.length; i++){
                    if(dirFiles[i].isDirectory()){
                        sb.append(getMD5Checksum(dirFiles[i]));
                    }else{
                         sb.append(getAnnotatedMD5ChecksumByFile(dirFiles[i]));
                    }
                }
            return sb.toString();
        }else{
            return getAnnotatedMD5ChecksumByFile(aFile);
        }
    }

    private static String getAnnotatedMD5ChecksumByFile(File fileIn ){
        String checksum = getMD5ChecksumByFile(fileIn);

        // get canonical path if possible so symbolic links are resolved
        String path = null;
        try{
            path = fileIn.getCanonicalPath();
        }catch(IOException ioe){ /* ignore */ }
        if(path == null) path = fileIn.getAbsolutePath();

        return path + " checksum: " + checksum;
    }
}
