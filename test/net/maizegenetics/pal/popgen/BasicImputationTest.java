/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.pal.popgen;

import java.util.Arrays;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.ImportUtils;
import net.maizegenetics.pal.alignment.Pack1Alignment;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author ed
 */
public class BasicImputationTest {
    String infile="test/datafiles/mdp_genotype.hmp.txt";
//    String basePath="/maizegenetics/test/datafiles/";  //there got to be someway to do it relative
//    String infile=basePath+"all_8pos.zip";
    Pack1Alignment p1a;

    public BasicImputationTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
    }

    @Before
    public void setUp() {
        p1a = (Pack1Alignment)ImportUtils.readFromHapmap(infile,"1");
        System.out.println("setUp");
    }

    @After
    public void tearDown() {
    }

    /**
     * Test of imputeBySite method, of class BasicImputation.
     */
    @Test
    public void testImputeBySite() {
        System.out.println("imputeBySite");
        Pack1Alignment align = p1a;
        int minLength = 30;
        int mismatchNum = 0;
        Pack1Alignment expResult = null;
        Pack1Alignment result = null;
        Pack1Alignment p2a=BasicImputation.imputeBySite(p1a, minLength, mismatchNum);
        assertEquals('T', p2a.getBaseChar(0, 14));
        // TODO review the generated test code and remove the default call to fail.
      //  fail("The test case is a prototype.");
    }

    /**
     * Test of maxHaplotypeLength method, of class BasicImputation.
     */
    @Test
    public void testMaxHaplotypeLength() {
        System.out.println("maxHaplotypeLength");
        Alignment align = p1a;
        int initialSite = 0;
        int taxa1 = 0;
        int taxa2 = 1;
        int mismatch = 0;
        int[] expResult = {1, 0, 1};
        int[] result = BasicImputation.maxHaplotypeLength(align, initialSite, taxa1, taxa2, mismatch);
        System.out.println(Arrays.toString(result));
        for(int i=0; i<result.length; i++) {
            assertEquals(expResult[i], result[i]);
        }
        // TODO review the generated test code and remove the default call to fail.
      //  fail("The test case is a prototype.");
    }

    /**
     * Test of maxHaplotypeLengthMatrix method, of class BasicImputation.
     */
    @Test
    public void testMaxHaplotypeLengthMatrix() {
        System.out.println("maxHaplotypeLengthMatrix");
        Alignment align = p1a;
        int initialSite = 0;
        int maxMismatch = 0;
        boolean maskKnown = false;
        int[] expResult = {0, 2, 7, 2, 1, 3, 1, 1, 7, 2, 7, 1, 2, 3, 3, 3, 1, -1, 1, 2, 3, 3, 3, 1, 1, 1, 3, 3, 7, 2, 3, 2, 1, 1, 2, 3, 2, 2, 3, 3, 3, 3, 1, 2, 2, 1, 2, 2, 1, 1, 2, 1, 126, 1, 3, 3, 2, 2, 1, 2, 2, 3, 1, 1, 1, 2, 2, 1, 1, 2, 3, 1, 3, 2, 2, 2, 2, 1, 2, 3, 2, -6, 3, 3, 2, 1, 2, 3, -1, 3, 3, 1, 2, 2, 2, 2, 1, 2, 3, 2, 3, 3, 1, 3, 2, 1, 2, 1, 2, 2, 2, 6, 2, 2, 1, 2, 1, -2, 1, 3, 2, 2, 3, 2, 1, 3, 7, 7, 2, 1, 3, 3, 1, 1, 2, 1, 3, 1, 1, 1, 2, 7, 1, 1, 3, 1, 7, 7, 7, 2, 2, 3, 2, 3, 2, 3, 3, 1, 7, 1, 2, 1, 2, 1, 2, 2, 3, 3, 3, 2, 7, 2, 1, 1, 3, -6, 1, 3, 7, 5, 7, 1, 1, 2, 2, 2, 2, 1, 2, 7, 2, 2, 2, 2, -1, 3, 2, 3, 3, 3, 3, 7, 7, 3, 3, 3, 2, 2, 2, 3, 2, 2, 2, 2, 3, 2, 2, 3, 1, 2, 3, 3, 2, 3, 1, 2, 2, 2, 2, 2, 3, 1, 2, 2, 2, 2, 2, 2, 2, 3, 1, 3, 2, 1, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 2, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 2, 3, 2, 2, 6, 2, 1, 2, 2};
        int[][] result = BasicImputation.maxHaplotypeLengthMatrix(align, initialSite, maxMismatch, maskKnown);
        System.out.println(Arrays.toString(result[0]));
        for(int i=0; i<result[0].length; i++) {
            assertEquals(expResult[i], result[0][i]);
        }
    }

}