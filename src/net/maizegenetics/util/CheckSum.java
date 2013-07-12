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
        byte[] b = null;
        try{
         b = createChecksum(aFile, "MD5");
        }catch(Exception e){ /* ignore */ } 
        String result = null;
        if(b != null){
            for (int i=0; i < b.length; i++) {
                result += Integer.toString( ( b[i] & 0xff ) + 0x100, 16).substring( 1 );
            }
        }
        return result;
    }


    public static String getMD5Checksum(File aFile) {
        if(aFile.isDirectory()){
            StringBuffer sb = new StringBuffer();
                File[] topmFiles = aFile.listFiles();
                for(int i=0; i<topmFiles.length; i++){
                    if(topmFiles[i].isDirectory()){ sb.append(getMD5Checksum(topmFiles[i])); }
                    String checksum = null;
                    try{
                        checksum = getMD5ChecksumByFile(topmFiles[i]);
                    }catch(Exception e){ /* ignore */  }

                    sb.append(topmFiles[i].getAbsolutePath() + " checksum: " + checksum + "\n");
                }
            return sb.toString();
        }else{
            return getMD5ChecksumByFile(aFile);
        }
    }
}
