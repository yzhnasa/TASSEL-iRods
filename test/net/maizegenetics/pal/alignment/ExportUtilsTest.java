/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.pal.alignment;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Ignore;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author yz79
 */
public class ExportUtilsTest {

    public static char delimChar = '\t';
    public static boolean isDiploid = false;
    public static String myHapmapTestFileName = "Temporary Hapmap Output File.hmp.txt";
    public static String myFlapjackTestFileName = "Temporary Flapjack Output File";
    public static String myPlinkTestFileName = "Temporary Plink Output File";
    public static String myZipTestFileName = "Temporary Zip Output File.BLOB.zip";

    public ExportUtilsTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
        Alignment hapmapAlign = ImportUtils.readFromHapmap("test/datafiles/mdp_genotype.hmp.txt");
        ExportUtils.writeToHapmap(hapmapAlign, isDiploid, myHapmapTestFileName, delimChar);
        Alignment flapjackAlign = ImportUtils.readFromFlapjack("test/datafiles/mdp_genotype.flpjk.geno", "test/datafiles/mdp_genotype.flpjk.map");
        ExportUtils.writeToFlapjack(flapjackAlign, myFlapjackTestFileName, delimChar);
        Alignment plinkAlign = ImportUtils.readFromPLINK("test/datafiles/mdp_genotype.plk.ped", "test/datafiles/mdp_genotype.plk.map");
        ExportUtils.writeToPlink(plinkAlign, myPlinkTestFileName, delimChar);

    }

    @AfterClass
    public static void tearDownClass() throws Exception {
        // Deleting temporary Hapmap file
//        File myHapmapFile = new File(myHapmapTestFileName);
//        myHapmapFile.delete();
        // Deleting temporary flapjack files
//        File myFlapjackGenoFile = new File(myFlapjackTestFileName + "flpjk.geno");
//        File myFlapjackMapFile = new File(myFlapjackTestFileName + "flpjak.map");
//        myFlapjackGenoFile.delete();
//        myFlapjackGenoFile.delete();
        // Deleting temporary plink files
//        File myPlinkPedFile = new File(myPlinkTestFileName + "plk.ped");
//        File myPlinkMapFile = new File(myPlinkTestFileName + "plk.map");
//        myPlinkPedFile.delete();
//        myPlinkMapFile.delete();
        // Deleting temporary Zip file
//        File myZipFile = new File(myZipTestFileName);
//        myZipFile.delete();
    }

    @Before
    public void setUp() {
    }

    @After
    public void tearDown() {
    }

    /**
     * Test of writeToHapmap method, of class ExportUtils.
     */
    @Test
    public void testWriteToHapmap() {
        System.out.println("JUnit test: writeToHapmap");
        String filename = "test/datafiles/mdp_genotype.hmp.txt";
        Alignment expResult = ImportUtils.readFromHapmap(myHapmapTestFileName);
        Alignment result = ImportUtils.readFromHapmap(filename);
        assertEquals(expResult.getSNPIDs(), result.getSNPIDs());
        assertEquals(expResult.getSiteCount(), result.getSiteCount());
        assertEquals(expResult.getSequenceCount(), result.getSequenceCount());
        for (int i = 0; i < expResult.getSiteCount(); i++) {
            assertEquals(expResult.getPositionInLocus(i), result.getPositionInLocus(i));
            for (int j = 0; j < expResult.getSequenceCount(); j++) {
                assertEquals(expResult.getBase(j, i), result.getBase(j, i));
            }
        }
    }

    /**
     * Test of writeToLZMA method, of class ExportUtils.
     */
    @Ignore
    @Test
    public void testWriteToLZMA() {
        System.out.println("writeToLZMA");
        Pack1Alignment alignment = null;
        String filenameRoot = "";
        ExportUtils.writeToLZMA(alignment, filenameRoot);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of writeToZip method, of class ExportUtils.
     */
    @Ignore
    @Test
    public void testWriteToZip() {
        System.out.println("writeToZip");
        Alignment alignment = null;
        String filenameRoot = "";
        ExportUtils.writeToZip(alignment, filenameRoot);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of writeToGZIP method, of class ExportUtils.
     */
    @Ignore
    @Test
    public void testWriteToGZIP() {
        System.out.println("writeToGZIP");
        Alignment alignment = null;
        String fileName = "";
        ExportUtils.writeToGZIP(alignment, fileName);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of writeToPlink method, of class ExportUtils.
     */
    @Test
    public void testWriteToPlink() {
        System.out.println("JUnit test: writeToPlink");
        String PEDfileName = "test/datafiles/mdp_genotype.plk.ped";
        String MAPfileName = "test/datafiles/mdp_genotype.plk.map";
        Alignment expResult = ImportUtils.readFromPLINK(myPlinkTestFileName + ".plk.ped", myPlinkTestFileName + ".plk.map");
        Alignment result = ImportUtils.readFromPLINK(PEDfileName, MAPfileName);
        assertEquals(expResult.getSNPIDs(), result.getSNPIDs());
        assertEquals(expResult.getSiteCount(), result.getSiteCount());
        assertEquals(expResult.getSequenceCount(), result.getSequenceCount());
        for (int i = 0; i < expResult.getSiteCount(); i++) {
            assertEquals(expResult.getPositionInLocus(i), result.getPositionInLocus(i));
            for (int j = 0; j < expResult.getSequenceCount(); j++) {
                assertEquals(expResult.getBase(j, i), result.getBase(j, i));
            }
        }
    }

    /**
     * Test of writeToFlapjack method, of class ExportUtils.
     */
    @Test
    public void testWriteToFlapjack() {
        System.out.println("JUnit test: writeToFlapjack");
        String genotypeFileName = "test/datafiles/mdp_genotype.flpjk.geno";
        String MAPfileName = "test/datafiles/mdp_genotype.flpjk.map";
        Alignment expResult = ImportUtils.readFromFlapjack(myFlapjackTestFileName + ".flpjk.geno", myFlapjackTestFileName + ".flpjk.map");
        Alignment result = ImportUtils.readFromFlapjack(genotypeFileName, MAPfileName);
        assertEquals(expResult.getSNPIDs(), result.getSNPIDs());
        assertEquals(expResult.getSiteCount(), result.getSiteCount());
        assertEquals(expResult.getSequenceCount(), result.getSequenceCount());
        for (int i = 0; i < expResult.getSiteCount(); i++) {
            assertEquals(expResult.getPositionInLocus(i), result.getPositionInLocus(i));
            for (int j = 0; j < expResult.getSequenceCount(); j++) {
                assertEquals(expResult.getBase(j, i), result.getBase(j, i));
            }
        }
    }

}