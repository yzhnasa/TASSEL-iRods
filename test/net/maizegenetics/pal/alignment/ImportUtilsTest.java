/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.pal.alignment;

import java.io.File;
import java.util.ArrayList;
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
public class ImportUtilsTest {

    public static char delimChar = '\t';
    public static boolean isDiploid = false;
    public static String myHapmapTestFileName = "Temporary Hapmap Output File.hmp.txt";
    public static String myFlapjackTestFileName = "Temporary Flapjack Output File";
    public static String myPlinkTestFileName = "Temporary Plink Output File";
    public static String myZipTestFileName = "Temporary Zip Output File.BLOB.zip";

    public ImportUtilsTest() {
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
     * Test of readFromHapmap method, of class ImportUtils.
     */
    @Test
    public void testReadFromHapmap_String_String() {
        System.out.println("JUnit test: readFromHapmap_String_String");
        String filename = "test/datafiles/mdp_genotype.hmp.txt";
        String chrom = "3";
        Alignment expResult = ImportUtils.readFromHapmap(myHapmapTestFileName, chrom);
        Alignment result = ImportUtils.readFromHapmap(filename, chrom);
        assertEquals(expResult.getSNPIDs(), result.getSNPIDs());
        assertEquals(expResult.getSiteCount(), result.getSiteCount());
        assertEquals(expResult.getSequenceCount(), result.getSequenceCount());
        for (int i = 0; i < expResult.getSiteCount(); i++) {
            assertEquals(expResult.getPositionInLocus(i), result.getPositionInLocus(i));
            for (int j = 0; j < expResult.getSequenceCount(); j++) {
                assertEquals(expResult.getBase(j, i), result.getBase(j, i));
            }
        }
        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
    }

    /**
     * Test of readFromHapmap method, of class ImportUtils.
     */
    @Test
    public void testReadFromHapmap_String() {
        System.out.println("JUnit test: readFromHapmap_String");
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
        // TODO review the generated test code and remove the default call to fail.
    }

    /**
     * Test of readFromLZMA method, of class ImportUtils.
     */
    @Ignore
    @Test
    public void testReadFromLZMA() {
        System.out.println("JUnit test: readFromLZMA");
        String filenameRoot = "";
        Pack1Alignment expResult = null;
        Pack1Alignment result = ImportUtils.readFromLZMA(filenameRoot);
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of readFromZip method, of class ImportUtils.
     */
    @Ignore
    @Test
    public void testReadFromZip() {
        System.out.println("JUnit test: readFromZip");
        String filenameRoot = "test/datafiles/all_8pos";
        Alignment expResult = ImportUtils.readFromZip(myZipTestFileName);
        Alignment result = ImportUtils.readFromZip(filenameRoot);
        assertEquals(expResult.getSNPIDs(), result.getSNPIDs());
        assertEquals(expResult.getSiteCount(), result.getSiteCount());
        assertEquals(expResult.getSequenceCount(), result.getSequenceCount());
        for (int i = 0; i < expResult.getSiteCount(); i++) {
            assertEquals(expResult.getPositionInLocus(i), result.getPositionInLocus(i));
            for (int j = 0; j < expResult.getSequenceCount(); j++) {
                assertEquals(expResult.getBase(j, i), result.getBase(j, i));
            }
        }
        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
    }

    /**
     * Test of readFromGZIP method, of class ImportUtils.
     */
    @Ignore
    @Test
    public void testReadFromGZIP() {
        System.out.println("JUnit test: readFromGZIP");
        String fileName = "";
        Alignment expResult = null;
        Alignment result = ImportUtils.readFromGZIP(fileName);
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of readFromPLINK method, of class ImportUtils.
     */
    @Test
    public void testReadFromPLINK_3args() {
        System.out.println("JUnit test: readFromPLINK_3args");
        String PEDfileName = "test/datafiles/mdp_genotype.plk.ped";
        String MAPfileName = "test/datafiles/mdp_genotype.plk.map";
        String chrom = "3";
        Alignment expResult = ImportUtils.readFromPLINK(myPlinkTestFileName + ".plk.ped", myPlinkTestFileName + ".plk.map", chrom);
        Alignment result = ImportUtils.readFromPLINK(PEDfileName, MAPfileName, chrom);
        assertEquals(expResult.getSNPIDs(), result.getSNPIDs());
        assertEquals(expResult.getSiteCount(), result.getSiteCount());
        assertEquals(expResult.getSequenceCount(), result.getSequenceCount());
        for (int i = 0; i < expResult.getSiteCount(); i++) {
            assertEquals(expResult.getPositionInLocus(i), result.getPositionInLocus(i));
            for (int j = 0; j < expResult.getSequenceCount(); j++) {
                assertEquals(expResult.getBase(j, i), result.getBase(j, i));
            }
        }
        // TODO review the generated test code and remove the default call to fail.
    }

    /**
     * Test of readFromPLINK method, of class ImportUtils.
     */
    @Test
    public void testReadFromPLINK_String_String() {
        System.out.println("JUnit test: readFromPLINK_String_String");
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
        // TODO review the generated test code and remove the default call to fail.
    }

    /**
     * Test of readFromFlapjack method, of class ImportUtils.
     */
    @Test
    public void testReadFromFlapjack_3args() {
        System.out.println("JUnit test: readFromFlapjack_3args");
        String genotypeFileName = "test/datafiles/mdp_genotype.flpjk.geno";
        String MAPfileName = "test/datafiles/mdp_genotype.flpjk.map";
        String chrom = "3";
        Alignment expResult = ImportUtils.readFromFlapjack(myFlapjackTestFileName + ".flpjk.geno", myFlapjackTestFileName + ".flpjk.map", chrom);
        Alignment result = ImportUtils.readFromFlapjack(genotypeFileName, MAPfileName, chrom);
        assertEquals(expResult.getSNPIDs(), result.getSNPIDs());
        assertEquals(expResult.getSiteCount(), result.getSiteCount());
        assertEquals(expResult.getSequenceCount(), result.getSequenceCount());
        for (int i = 0; i < expResult.getSiteCount(); i++) {
            assertEquals(expResult.getPositionInLocus(i), result.getPositionInLocus(i));
            for (int j = 0; j < expResult.getSequenceCount(); j++) {
                assertEquals(expResult.getBase(j, i), result.getBase(j, i));
            }
        }
        // TODO review the generated test code and remove the default call to fail.
    }

    /**
     * Test of readFromFlapjack method, of class ImportUtils.
     */
    @Test
    public void testReadFromFlapjack_String_String() {
        System.out.println("JUnit test: readFromFlapjack_String_String");
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
        // TODO review the generated test code and remove the default call to fail.
    }

    /**
     * Test of createPack1AlignmentFromFile method, of class ImportUtils.
     */
    @Ignore
    @Test
    public void testCreatePack1AlignmentFromFile_ArrayList() {
        System.out.println("JUnit test: createPack1AlignmentFromFile_ArrayList");
        ArrayList<String> fileList = null;
        Alignment expResult = null;
        Alignment result = ImportUtils.createPack1AlignmentFromFile(fileList);
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of createPack1AlignmentFromFile method, of class ImportUtils.
     */
    @Ignore
    @Test
    public void testCreatePack1AlignmentFromFile_3args() {
        System.out.println("JUnit test: createPack1AlignmentFromFile_3args");
        String alleleFileName = "";
        String coverageFileName = "";
        String refseqFileName = "";
        Alignment expResult = null;
        Alignment result = ImportUtils.createPack1AlignmentFromFile(alleleFileName, coverageFileName, refseqFileName);
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

}