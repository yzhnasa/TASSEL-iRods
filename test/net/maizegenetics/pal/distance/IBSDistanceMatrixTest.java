/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.pal.distance;

import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.BitAlignment;
import net.maizegenetics.pal.alignment.Locus;
import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.SimpleIdGroup;
import net.maizegenetics.prefs.TasselPrefs;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author edbuckler
 */
public class IBSDistanceMatrixTest {
    
    public IBSDistanceMatrixTest() {
    }
    
    @BeforeClass
    public static void setUpClass() {
    }
    
    @AfterClass
    public static void tearDownClass() {
    }
    
    @Before
    public void setUp() {
    }
    
    @After
    public void tearDown() {
    }
    
    private Alignment createAlignment(String[] taxaNames, String[] sequenceArray) {
        if(taxaNames==null) {
            taxaNames=new String[sequenceArray.length];
            for (int i = 0; i < sequenceArray.length; i++) {
                taxaNames[i]="Taxa"+i;   
            }
        }
        IdGroup idGroup = new SimpleIdGroup(taxaNames);

        Locus unknown = new Locus("Unknown", "0", 0, sequenceArray[0].length(), null, null);
        return BitAlignment.getNucleotideInstance(idGroup, sequenceArray, null, null, null, 
                TasselPrefs.getAlignmentMaxAllelesToRetain(), new Locus[]{unknown}, new int[]{0}, 
                null, TasselPrefs.getAlignmentRetainRareAlleles(), false);
    }

    /**
     * Test of computeHetBitDistances method, of class IBSDistanceMatrix.
     */
    @Test
    public void testComputeHetBitDistances_3args() {
        System.out.println("computeHetBitDistances");
        int taxon1 = 0;
        int taxon2 = 1;
        String[] seqs={"AAA","AAA"};
        Alignment theTBA = createAlignment(null,seqs);
        double result = IBSDistanceMatrix.computeHetBitDistances(theTBA, taxon1, taxon2);
        assertEquals(0.0, result, 0.00001);
        seqs[1]="CAA";
        theTBA = createAlignment(null,seqs);
        result = IBSDistanceMatrix.computeHetBitDistances(theTBA, taxon1, taxon2);
        assertEquals(0.3333333333, result, 0.00001);
        seqs[1]="CGT";
        theTBA = createAlignment(null,seqs);
        result = IBSDistanceMatrix.computeHetBitDistances(theTBA, taxon1, taxon2);
        assertEquals(1.0, result, 0.00001);
        seqs[1]="AAR";
        theTBA = createAlignment(null,seqs);
        result = IBSDistanceMatrix.computeHetBitDistances(theTBA, taxon1, taxon2);
        assertEquals(1.0/6.0, result, 0.00001);
        seqs[0]="RMR";
        seqs[1]="RMR";
        theTBA = createAlignment(null,seqs);
        result = IBSDistanceMatrix.computeHetBitDistances(theTBA, taxon1, taxon2);
        assertEquals(0.5, result, 0.00001);
    }

    /**
     * Test of computeHetBitDistances method, of class IBSDistanceMatrix.
     */
    @Test
    public void testComputeHetBitDistances_5args_1() {
        System.out.println("computeHetBitDistances");
        Alignment theTBA = null;
        int taxon1 = 0;
        int taxon2 = 1;
        int minSitesCompared = 0;
        boolean isTrueIBS = false;
        double expResult = 0.0;
        double result = IBSDistanceMatrix.computeHetBitDistances(theTBA, taxon1, taxon2, minSitesCompared, isTrueIBS);
        assertEquals(expResult, result, 0.0);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of computeHetBitDistances method, of class IBSDistanceMatrix.
     */
    @Test
    public void testComputeHetBitDistances_5args_2() {
        System.out.println("computeHetBitDistances");
        long[] iMj = null;
        long[] iMn = null;
        long[] jMj = null;
        long[] jMn = null;
        int minSitesCompared = 0;
        double expResult = 0.0;
        double result = IBSDistanceMatrix.computeHetBitDistances(iMj, iMn, jMj, jMn, minSitesCompared);
        assertEquals(expResult, result, 0.0);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of getAverageTotalSites method, of class IBSDistanceMatrix.
     */
    @Test
    public void testGetAverageTotalSites() {
        System.out.println("getAverageTotalSites");
        IBSDistanceMatrix instance = null;
        double expResult = 0.0;
        double result = instance.getAverageTotalSites();
        assertEquals(expResult, result, 0.0);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of getDistance method, of class IBSDistanceMatrix.
     */
    @Test
    public void testGetDistance() {
        System.out.println("getDistance");
        IBSDistanceMatrix instance = null;
        double[][] expResult = null;
        double[][] result = instance.getDistance();
        assertArrayEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of toString method, of class IBSDistanceMatrix.
     */
    @Test
    public void testToString_int() {
        System.out.println("toString");
        int d = 0;
        IBSDistanceMatrix instance = null;
        String expResult = "";
        String result = instance.toString(d);
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of toString method, of class IBSDistanceMatrix.
     */
    @Test
    public void testToString_0args() {
        System.out.println("toString");
        IBSDistanceMatrix instance = null;
        String expResult = "";
        String result = instance.toString();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of fireProgress method, of class IBSDistanceMatrix.
     */
    @Test
    public void testFireProgress() {
        System.out.println("fireProgress");
        int percent = 0;
        IBSDistanceMatrix instance = null;
        instance.fireProgress(percent);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }

    /**
     * Test of isTrueIBS method, of class IBSDistanceMatrix.
     */
    @Test
    public void testIsTrueIBS() {
        System.out.println("isTrueIBS");
        IBSDistanceMatrix instance = null;
        boolean expResult = false;
        boolean result = instance.isTrueIBS();
        assertEquals(expResult, result);
        // TODO review the generated test code and remove the default call to fail.
        fail("The test case is a prototype.");
    }
}
