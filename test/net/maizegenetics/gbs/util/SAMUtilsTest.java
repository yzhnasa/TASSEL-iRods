/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.gbs.util;

import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author jvh39
 */
public class SAMUtilsTest {

    public SAMUtilsTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
    }

    @Test
    public void testGetVariants() {
    }

    @Test
    public void testEndCoordinate() {
        assertEquals(
            SAMUtils.endCoordinate("64M", 100), //CIGAR code and start coordinate
            new Long(163) //Expected end coordinate
        );

        assertEquals(
            SAMUtils.endCoordinate("17M1D2M1I20M1D24M", 222278120),
            new Long(222278184)
        );
        assertEquals(
            SAMUtils.endCoordinate("22M1I29M3D11M1S", 218794214),
            new Long(218794279)
        );

        assertEquals(
            SAMUtils.endCoordinate("6M1I1M2D56M", 209120270),
            new Long(209120334)
        );

        assertEquals(
            SAMUtils.endCoordinate("35M1I24M1D3M1S",138877212),
            new Long(138877275)
        );
    }

    @Test
    public void testParseCIGAR() {
    }

    @Test
    public void testParseMDField() {
    }

}