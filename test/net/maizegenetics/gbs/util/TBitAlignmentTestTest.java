/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.gbs.util;

import net.maizegenetics.pal.datatype.DataType;
import net.maizegenetics.pal.alignment.AlignmentUtils;
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
 * @author edbuckler
 */
public class TBitAlignmentTestTest {
    String infile="test/datafiles/mdp_genotype.hmp.txt";
    Pack1Alignment p1a;
    TBitAlignmentTest theTBA;

    public TBitAlignmentTestTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
    }

    @Before
    public void setUp() {
     //   String infile="/Users/edbuckler/SolexaAnal/GBS/hmp/HS55K_110215.hmp.txt";
        p1a = (Pack1Alignment)ImportUtils.readFromHapmap(infile,"1");
    //    p1a=(Pack1Alignment)ImportUtils.createPack1AlignmentFromFile(infile2,"","");
        theTBA=new TBitAlignmentTest(p1a);
        long time=System.currentTimeMillis();
        int count=0;
        for (int i = 0; i < 200000000; i++) {
            //if(p1a.getBase(i%16, i%512)==DataType.UNKNOWN_BYTE) count++;
            if(theTBA.getBase(i%16, i%512)==DataType.UNKNOWN_BYTE) count++;
        }
        System.out.println("Time: "+(System.currentTimeMillis()-time)+" count:"+count);
        System.out.println(AlignmentUtils.getSequenceString(p1a, 0));
        System.out.println(AlignmentUtils.getSequenceString(theTBA, 0));
        System.out.println("setUp");
        time=System.currentTimeMillis();
        count=0;
        int same=0;
        for (int i = 0; i < p1a.getSequenceCount(); i++) {
            for (int j = 0; j < i; j++) {
                for (int k = 0; k < p1a.getSiteCount(); k++) {
                    if((p1a.getBase(i, k)==DataType.UNKNOWN_BYTE)||(p1a.getBase(j, k)==DataType.UNKNOWN_BYTE)) continue;
                    if(p1a.getBase(i, k)==p1a.getBase(j, k)) same++;
                    count++;
                }
            }
        }
        System.out.println("Basic Get Base - Time: "+(System.currentTimeMillis()-time)+" count:"+count+" same:"+same);
        time=System.currentTimeMillis();
        count=0;
        same=0;
        for (int i = 0; i < p1a.getSequenceCount(); i++) {
            for (int j = 0; j < i; j++) {
                byte[] b1=p1a.getBaseRange(i, 0, p1a.getSiteCount()-1);
                byte[] b2=p1a.getBaseRange(j, 0, p1a.getSiteCount()-1);
                for (int k = 0; k < p1a.getSiteCount(); k++) {
                    if((b1[k]==DataType.UNKNOWN_BYTE)||(b2[k]==DataType.UNKNOWN_BYTE)) continue;
                    if(b1[k]==b2[k]) same++;
                    count++;
                }
            }
        }
        System.out.println("Basic Get Base Range - Time: "+(System.currentTimeMillis()-time)+" count:"+count+" same:"+same);
        time=System.currentTimeMillis();
        count=0;
        same=0;
        for (int i = 0; i < p1a.getSequenceCount(); i++) {
            for (int j = 0; j < i; j++) {
                int[][] s=new int[2][2];
                for (int k = 0; k < 2; k++) {
                    for (int m = 0; m < 2; m++) {
                        //this could be made even faster by combining the next three into on function
//                        OpenBitSet b1=theTBA.getTaxaBitsWithClone(i, k);
//                        OpenBitSet b2=theTBA.getTaxaBitsWithClone(j, m);
//                        b1.and(b2);
//                        s[k][m] =(int)b1.cardinality();
                        s[k][m] =(int)OpenBitSet.intersectionCount(theTBA.getTaxaBitsNoClone(i, k), theTBA.getTaxaBitsWithClone(j, m));
                        count+=s[k][m];
                    }
                }
                same+=s[0][0]+s[1][1];
            }
        }
        System.out.println("Bit Intersect - Time: "+(System.currentTimeMillis()-time)+" count:"+count+" same:"+same);

    }

    @After
    public void tearDown() {
    }

    /**
     * Test of getBaseChar method, of class TBitAlignmentTest.
     */
    @Test
    public void testGetBaseChar() {
        System.out.println("getDataChar");
        char result = p1a.getBaseChar(0, 2);
        assertEquals('G', result);
        result = p1a.getBaseChar(210, 3);
        assertEquals('W', result);
    }
//
//    /**
//     * Test of getBase method, of class TBitAlignmentTest.
//     */
//    @Test
//    public void testGetBase_int_int() {
//        System.out.println("getBase");
//        int taxon = 0;
//        int site = 0;
//        TBitAlignmentTest instance = null;
//        byte expResult = 0;
//        byte result = instance.getBase(taxon, site);
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of getBase method, of class TBitAlignmentTest.
//     */
//    @Test
//    public void testGetBase_3args() {
//        System.out.println("getBase");
//        int taxon = 0;
//        int site = 0;
//        int allele = 0;
//        TBitAlignmentTest instance = null;
//        byte expResult = 0;
//        byte result = instance.getBase(taxon, site, allele);
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of getBase method, of class TBitAlignmentTest.
//     */
//    @Test
//    public void testGetBase_4args() {
//        System.out.println("getBase");
//        int taxon = 0;
//        Locus locus = null;
//        int physicalPosition = 0;
//        int allele = 0;
//        TBitAlignmentTest instance = null;
//        byte expResult = 0;
//        byte result = instance.getBase(taxon, locus, physicalPosition, allele);
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of getSNPIDs method, of class TBitAlignmentTest.
//     */
//    @Test
//    public void testGetSNPIDs() {
//        System.out.println("getSNPIDs");
//        TBitAlignmentTest instance = null;
//        String[] expResult = null;
//        String[] result = instance.getSNPIDs();
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of getSNPID method, of class TBitAlignmentTest.
//     */
//    @Test
//    public void testGetSNPID() {
//        System.out.println("getSNPID");
//        int site = 0;
//        TBitAlignmentTest instance = null;
//        String expResult = "";
//        String result = instance.getSNPID(site);
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of getSiteCount method, of class TBitAlignmentTest.
//     */
//    @Test
//    public void testGetSiteCount() {
//        System.out.println("getSiteCount");
//        TBitAlignmentTest instance = null;
//        int expResult = 0;
//        int result = instance.getSiteCount();
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of getLocusSiteCount method, of class TBitAlignmentTest.
//     */
//    @Test
//    public void testGetLocusSiteCount() {
//        System.out.println("getLocusSiteCount");
//        Locus locus = null;
//        TBitAlignmentTest instance = null;
//        int expResult = 0;
//        int result = instance.getLocusSiteCount(locus);
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of getPositionInLocus method, of class TBitAlignmentTest.
//     */
//    @Test
//    public void testGetPositionInLocus() {
//        System.out.println("getPositionInLocus");
//        int site = 0;
//        TBitAlignmentTest instance = null;
//        int expResult = 0;
//        int result = instance.getPositionInLocus(site);
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of getSiteOfPhysicalPosition method, of class TBitAlignmentTest.
//     */
//    @Test
//    public void testGetSiteOfPhysicalPosition() {
//        System.out.println("getSiteOfPhysicalPosition");
//        int physicalPosition = 0;
//        Locus locus = null;
//        TBitAlignmentTest instance = null;
//        int expResult = 0;
//        int result = instance.getSiteOfPhysicalPosition(physicalPosition, locus);
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of getPositionType method, of class TBitAlignmentTest.
//     */
//    @Test
//    public void testGetPositionType() {
//        System.out.println("getPositionType");
//        int site = 0;
//        TBitAlignmentTest instance = null;
//        byte expResult = 0;
//        byte result = instance.getPositionType(site);
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of getPositionTypes method, of class TBitAlignmentTest.
//     */
//    @Test
//    public void testGetPositionTypes() {
//        System.out.println("getPositionTypes");
//        TBitAlignmentTest instance = null;
//        byte[] expResult = null;
//        byte[] result = instance.getPositionTypes();
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of getLocus method, of class TBitAlignmentTest.
//     */
//    @Test
//    public void testGetLocus() {
//        System.out.println("getLocus");
//        int site = 0;
//        TBitAlignmentTest instance = null;
//        Locus expResult = null;
//        Locus result = instance.getLocus(site);
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of getLoci method, of class TBitAlignmentTest.
//     */
//    @Test
//    public void testGetLoci() {
//        System.out.println("getLoci");
//        TBitAlignmentTest instance = null;
//        Locus[] expResult = null;
//        Locus[] result = instance.getLoci();
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
//
//    /**
//     * Test of getNumLoci method, of class TBitAlignmentTest.
//     */
//    @Test
//    public void testGetNumLoci() {
//        System.out.println("getNumLoci");
//        TBitAlignmentTest instance = null;
//        int expResult = 0;
//        int result = instance.getNumLoci();
//        assertEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }

}