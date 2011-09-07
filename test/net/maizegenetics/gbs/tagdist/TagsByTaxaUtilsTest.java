/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.gbs.tagdist;

import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.*;
import java.io.BufferedReader;
import java.io.FileReader;

/** @author jvh39 */
public class TagsByTaxaUtilsTest {

    public TagsByTaxaUtilsTest() {
    }

    @Before
    public void setUp() {
    }



    /*@Test
    public void testBitDistFile() {
        ArrayList<String> mergedTaxaNameList = new ArrayList<String>();
        TreeSet<String> duplicateTaxa = new TreeSet<String>();
        String testFileName = "/media/CA568FCE568FB9A9/tbt/Maize_282110315.tbt.bit";
        int expectedTaxaNumber = 0;
        int expectedTagNumber = 0;
        int expectedTagLengthInLong = 0;
        String currTaxonName;

        try {
            DataInputStream rw = new DataInputStream(
                    new BufferedInputStream(
                    new FileInputStream(
                    new File(testFileName)),
                    4000000));

            expectedTagNumber = rw.readInt(); //Read header
            expectedTagLengthInLong = rw.readInt();
            expectedTaxaNumber = rw.readInt();
            String[] taxonNames = new String[expectedTaxaNumber];

            //Add taxon names to TreeSet.  "Set" objects are checked before adding
            //new elements so that they never contain duplicates.
            for (int t = 0; t < expectedTaxaNumber; t++) {
                currTaxonName = rw.readUTF();
                duplicateTaxa.add(currTaxonName);
                taxonNames[t] = currTaxonName;
            }
            assertEquals(expectedTaxaNumber, taxonNames.length);

            long[] currTag = new long[expectedTagLengthInLong]; 
            int longsInBitset = OpenBitSet.bits2words(expectedTaxaNumber);
            long[] distInLong = new long[longsInBitset]; //Array for current bitset
            byte expectedTagLength;
            String distInText;

            //  Begin loop of PAIN
            for (int i = 0; i < expectedTagNumber; i++) {
                    if(i%1000==0){System.out.println("Checked "+i+" tags.");}
                for (int j = 0; j < expectedTagLengthInLong; j++) {currTag[j] = rw.readLong();}//Read tag sequence
                expectedTagLength = rw.readByte();                 //Read tag length
                for (int j = 0; j < longsInBitset; j++) {distInLong[j] = rw.readLong();}  //Read distribution into array
                OpenBitSet obs = new OpenBitSet(distInLong, expectedTaxaNumber); //Convert array to bitset
                assertNotNull(obs);
                for(int v=0; v<expectedTaxaNumber; v++){
                    if (duplicateTaxa.contains(taxonNames[v])) {
                        System.out.println(taxonNames[v]+" is a duplicate.  Its distribution string is:");
                            System.out.println(TagsByTaxaUtils.bitsetToString(obs));
                    }
                }
            }
                rw.close();
        } catch (Exception e) {
            System.out.println("Caught exception while merging output file: " + e);
            e.printStackTrace();
        }

    }*/
}
