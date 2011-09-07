package net.maizegenetics.pal.alignment;


import junit.framework.Assert;

import net.maizegenetics.pal.datatype.DataType;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

public class HetFunctionTest {
    Pack1Alignment p1a;
    Alignment sa;

	@Before
	public void setUp() throws Exception {
        p1a = (Pack1Alignment)ImportUtils.readFromHapmap("test/datafiles/mdp_genotype.hmp.txt","1");
		String markerFile = "test/datafiles/diploid_SSR.txt";
		sa = ReadPolymorphismUtils.readPolymorphismFile(markerFile);
	}

	@After
	public void tearDown() throws Exception {
	}

	@Test
	public void testHetFunctionsForTextData() {
		
		try {
			int[] index = new int[]{0,20,21};
			int n = index.length;
			DataType dt = sa.getDataType();
			
			//test isHeterozygous
			for (int j = 0; j<n; j++) {
				for (int i =0; i <6; i++) {
					if (j == 2 && i == 3) {
						String msg = "isHeterozygous failed at j=2, i=3";
						Assert.assertTrue(msg, dt.isHeterozygote(sa.getBaseChar(index[j],i)));
					}
					else {
						String msg = "is Heterozygous failed at i=" + i + ", j=" + j;
						Assert.assertFalse(msg, dt.isHeterozygote(sa.getBaseChar(index[j],i)));
					}
				}
			}
			
			//test isDiploid identity
			char[] c = new char[10];
			c[0] = sa.getBaseChar(0, 0); //123:123
			c[1] = sa.getBaseChar(3, 0); //119:119
			c[2] = sa.getBaseChar(6, 2); //?:?
			c[3] = sa.getBaseChar(7, 0); //119:119
			c[4] = sa.getBaseChar(18, 3); //165:165
			c[5] = sa.getBaseChar(21, 3); //165:200
			c[6] = sa.getBaseChar(0, 3); //200:200
			c[7] = sa.getBaseChar(51, 3); //165:183
			c[8] = sa.getBaseChar(1, 3); //166:183
			
			
			//test unequal homozygotes
			Assert.assertEquals("text diploid identity test 1 failed", 0, dt.getDiploidIdentity(c[0], c[1]));
			//test equal homozygotes
			Assert.assertEquals("text diploid identity test 2 failed", 4, dt.getDiploidIdentity(c[1], c[3]));
			//test one genotype missing
			Assert.assertEquals("text diploid identity test 3 failed", -1, dt.getDiploidIdentity(c[0], c[2]));
			//test equal heterozygotes
			Assert.assertEquals("text diploid identity test 4 failed", 2, dt.getDiploidIdentity(c[5], c[5]));
			//test heterozygotes with one allele equal
			Assert.assertEquals("text diploid identity test 5 failed", 1, dt.getDiploidIdentity(c[5], c[7]));
			//test unequal heterozygotes
			Assert.assertEquals("text diploid identity test 6 failed", 0, dt.getDiploidIdentity(c[5], c[8]));
			//test homozygote vs heterozygote, one allele the same
			Assert.assertEquals("text diploid identity test 7 failed", 2, dt.getDiploidIdentity(c[5], c[6]));
			//test homozygote vs heterozygote unequal
			Assert.assertEquals("text diploid identity test 8 failed", 0, dt.getDiploidIdentity(c[0], c[7]));
			
			
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	@Test
	public void testIUPACHetFunctions() {
		DataType dt = p1a.getDataType();
		int nsites = p1a.getSiteCount();
		int nseq = p1a.getSequenceCount();
		
		//test isHet
		for (int i = 0; i < nseq; i++) {
			for (int j = 0; j < nsites; j++){
				char base = p1a.getBaseChar(i, j);
				Assert.assertTrue("IUPAC isHet failed at seq=" + i + ", site=" + j, dt.isHeterozygote(base)==confirmHet(base));
			}
		}
		
		//test get identity
		//bases = {'A', 'C', 'G', 'T', 'R', 'Y', 'S', 'W', 'K', 'M', '+', '0', 'H', 'V', 'N', '-'};
		//K = G/T, M = A/C, R = A/G, S = G/C, W = A/T, Y = C/T
		
		byte[] bases = GdpdmBLOBUtils.bases;
		//homozygotes equal
		Assert.assertEquals("IUPAC diploid identity test 1 failed", 4, dt.getDiploidIdentity(bases[0], bases[0]));
		Assert.assertEquals("IUPAC diploid identity test 2 failed", 4, dt.getDiploidIdentity(bases[1], bases[1]));
		Assert.assertEquals("IUPAC diploid identity test 3 failed", 4, dt.getDiploidIdentity(bases[2], bases[2]));
		Assert.assertEquals("IUPAC diploid identity test 4 failed", 4, dt.getDiploidIdentity(bases[3], bases[3]));
		Assert.assertEquals("IUPAC diploid identity test 5 failed", 4, dt.getDiploidIdentity(bases[10], bases[10]));
		Assert.assertEquals("IUPAC diploid identity test 6 failed", 4, dt.getDiploidIdentity(bases[15], bases[15]));

		//heterozygotes equal
		Assert.assertEquals("IUPAC diploid identity test 7 failed", 2, dt.getDiploidIdentity(bases[4], bases[4]));
		Assert.assertEquals("IUPAC diploid identity test 8 failed", 2, dt.getDiploidIdentity(bases[5], bases[5]));
		Assert.assertEquals("IUPAC diploid identity test 9 failed", 2, dt.getDiploidIdentity(bases[6], bases[6]));
		Assert.assertEquals("IUPAC diploid identity test 10 failed", 2, dt.getDiploidIdentity(bases[7], bases[7]));
		Assert.assertEquals("IUPAC diploid identity test 11 failed", 2, dt.getDiploidIdentity(bases[8], bases[8]));
		Assert.assertEquals("IUPAC diploid identity test 12 failed", 2, dt.getDiploidIdentity(bases[9], bases[9]));
		Assert.assertEquals("IUPAC diploid identity test 13 failed", 2, dt.getDiploidIdentity(bases[11], bases[11]));

		//heterozygotes that share one allele
		Assert.assertEquals("IUPAC diploid identity test 14 failed", 1, dt.getDiploidIdentity(bases[4], bases[6]));
		Assert.assertEquals("IUPAC diploid identity test 15 failed", 1, dt.getDiploidIdentity(bases[6], bases[8]));
		Assert.assertEquals("IUPAC diploid identity test 16 failed", 1, dt.getDiploidIdentity(bases[7], bases[5]));
		Assert.assertEquals("IUPAC diploid identity test 17 failed", 1, dt.getDiploidIdentity(bases[8], bases[5]));

		//heterozygotes unequal
		Assert.assertEquals("IUPAC diploid identity test 18 failed", 0, dt.getDiploidIdentity(bases[4], bases[5]));
		Assert.assertEquals("IUPAC diploid identity test 19 failed", 0, dt.getDiploidIdentity(bases[6], bases[7]));
		Assert.assertEquals("IUPAC diploid identity test 21 failed", 0, dt.getDiploidIdentity(bases[8], bases[9]));
		
		//homozygotes unequal
		Assert.assertEquals("IUPAC diploid identity test 22 failed", 0, dt.getDiploidIdentity(bases[0], bases[1]));
		Assert.assertEquals("IUPAC diploid identity test 23 failed", 0, dt.getDiploidIdentity(bases[2], bases[3]));

	}
	
	boolean confirmHet(char base) {
		if (base=='R' || base=='Y'|| base=='S' || base=='W'|| base=='K'|| base=='M'|| base=='0') return true;
		return false;
	}
}
