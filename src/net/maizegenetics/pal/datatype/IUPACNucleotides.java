// IUPACNucleotides.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)
package net.maizegenetics.pal.datatype;

import java.util.Arrays;

/**
 * implements DataType for nucleotides with ambiguous characters
 *
 * @version $Id: IUPACNucleotides.java,v 1.1 2007/01/12 03:26:13 tcasstevens Exp $
 *
 * @author Alexei Drummond
 */
public class IUPACNucleotides extends SimpleDataType {

    Class stateSize=Byte.class;
    private static final char[] DNA_CONVERSION_TABLE = {'A', 'C', 'G', 'T', 'K', 'M', 'R', 'S', 'W', 'Y', '+', '0', 'X', 'X', UNKNOWN_CHARACTER, GAP_CHARACTER};
    private static final byte[] HETS_BYTES = {'K', 'M', 'R', 'S', 'W', 'Y', '0'};
    private static final char[] HETS_CHAR = {'K', 'M', 'R', 'S', 'W', 'Y', '0'};
    public static final byte[] AIUPAC = {(byte)'A',(byte)'A'};
    public static final byte[] CIUPAC = {(byte)'C',(byte)'C'};
    public static final byte[] GIUPAC = {(byte)'G',(byte)'G'};
    public static final byte[] TIUPAC = {(byte)'T',(byte)'T'};
    public static final byte[] RIUPAC = {(byte)'A', (byte)'G'};
    public static final byte[] YIUPAC = {(byte)'C', (byte)'T'};
    public static final byte[] SIUPAC = {(byte)'G', (byte)'C'};
    public static final byte[] WIUPAC = {(byte)'A', (byte)'T'};
    public static final byte[] KIUPAC = {(byte)'G', (byte)'T'};
    public static final byte[] MIUPAC = {(byte)'A', (byte)'C'};
    public static final byte[] BIUPAC = {(byte)'+',(byte)'+'};
    public static final byte[] DIUPAC = {(byte)'+',(byte)'-'};
    public static final byte[] NIUPAC = {(byte)'N',(byte)'N'};
    public static final byte[] gapIUPAC = {(byte)'-',(byte)'-'};
    
    /* Translation table for diploid identity
		0	1	2	3	4	5	6	7	8	9	10	11	12	13	14	15
		A/A	C/C	G/G	T/T	A/G	C/T	G/C	A/T	G/T	A/C	+/+	+/-	RES	RES	n	-/-
0	A/A	4	0	0	0	2	0	0	2	0	2	0	0	-1	-1	0	0
1	C/C	0	4	0	0	0	2	2	0	0	2	0	0	-1	-1	0	0
2	G/G	0	0	4	0	2	0	2	0	2	0	0	0	-1	-1	0	0
3	T/T	0	0	0	4	0	2	0	2	2	0	0	0	-1	-1	0	0
4	A/G	2	0	2	0	2	0	1	1	1	1	0	0	-1	-1	0	0
5	C/T	0	2	0	2	0	2	1	1	1	1	0	0	-1	-1	0	0
6	G/C	0	2	2	0	1	1	2	0	1	1	0	0	-1	-1	0	0
7	A/T	2	0	0	2	1	1	0	2	1	1	0	0	-1	-1	0	0
8	G/T	0	0	2	2	1	1	1	1	2	0	0	0	-1	-1	0	0
9	A/C	2	2	0	0	1	1	1	1	0	2	0	0	-1	-1	0	0
10	+/+	0	0	0	0	0	0	0	0	0	0	4	2	-1	-1	0	0
11	+/-	0	0	0	0	0	0	0	0	0	0	2	2	-1	-1	0	2
12	RES	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1
13	RES	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1
14	n	0	0	0	0	0	0	0	0	0	0	0	0	-1	-1	0	0
15	-/-	0	0	0	0	0	0	0	0	0	0	0	2	-1	-1	0	4

    */
    
    private static final int[][] diploidIDTable = {
    	{4,0,0,0,2,0,0,2,0,2,0,0,-1,-1,0,0},
    	{0,4,0,0,0,2,2,0,0,2,0,0,-1,-1,0,0},
    	{0,0,4,0,2,0,2,0,2,0,0,0,-1,-1,0,0},
    	{0,0,0,4,0,2,0,2,2,0,0,0,-1,-1,0,0},
    	{2,0,2,0,2,0,1,1,1,1,0,0,-1,-1,0,0},
    	{0,2,0,2,0,2,1,1,1,1,0,0,-1,-1,0,0},
    	{0,2,2,0,1,1,2,0,1,1,0,0,-1,-1,0,0},
    	{2,0,0,2,1,1,0,2,1,1,0,0,-1,-1,0,0},
    	{0,0,2,2,1,1,1,1,2,0,0,0,-1,-1,0,0},
    	{2,2,0,0,1,1,1,1,0,2,0,0,-1,-1,0,0},
    	{0,0,0,0,0,0,0,0,0,0,4,2,-1,-1,0,0},
    	{0,0,0,0,0,0,0,0,0,0,2,2,-1,-1,0,2},
    	{-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    	{-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    	{0,0,0,0,0,0,0,0,0,0,0,0,-1,-1,0,0},
    	{0,0,0,0,0,0,0,0,0,0,0,2,-1,-1,0,4}
    	};
    
    //Must stay after static CONVERSION_TABLE stuff!
    public static final IUPACNucleotides DEFAULT_INSTANCE = new IUPACNucleotides();
    public static final IUPACNucleotides DNA_INSTANCE = new IUPACNucleotides();

    char[] conversionTable_;

    //
    // Serialization Stuff
    //
    private static final long serialVersionUID = 8863411606027017687L;

  

    public IUPACNucleotides() {
    }


    /**
     * Get number of states.
     */
    public int getNumStates() {
        return 16;
    }

    public final int getState(final char c) {
        switch (c) {
            case 'A':
                return 0;
            case 'C':
                return 1;
            case 'G':
                return 2;
            case 'T':
                return 3;
            case 'U':
                return 3;
            case 'R':
                return 4;
            case 'Y':
                return 5;
            case 'S':
                return 6;
            case 'W':
                return 7;
            case 'K':
                return 8;
            case 'M':
                return 9;
            case '+':
                return 10;
            case '0':
                return 11;
            case 'X':
                return 12;
//            case 'X':
//                return 13;
            case UNKNOWN_CHARACTER:
                return 14;
            case GAP_CHARACTER:
                return 15;

            case 'a':
                return 0;
            case 'c':
                return 1;
            case 'g':
                return 2;
            case 't':
                return 3;
            case 'u':
                return 3;
            case 'r':
                return 4;
            case 'y':
                return 5;
            case 's':
                return 6;
            case 'w':
                return 7;
            case 'k':
                return 8;
            case 'm':
                return 9;
            case 'x':
                return 12;
//            case 'v':
//                return 13;
            case 'n':
                return 14;
        }
        return -1;
    }

    // Get character corresponding to a given state
    public char getChar(final int state) {
        if (state < conversionTable_.length && state >= 0) {
            return conversionTable_[state];
        }
        return UNKNOWN_CHARACTER;
    }

    // String describing the data type
    public String getDescription() {
        return "IUPACNucleotide";
    }

//==========================================================
//================ Static Utility Methods ===================
//==========================================================
    private static final int[] IUPAC_CONV = {
        -1, //Nothing... eh?
        0 /* A */, 1 /* C */, 5 /* CT */, 2 /* G */,
        6 /* GC */, 7 /* AT */, 13 /* Reserved */, 3 /* T */,
        8 /* GT */, 9 /* AC */, 12 /* Reserved */, 4 /* AG */,
        11 /* 0/0 (indel het.) */, 10 /* +/+ (insertion homoz) */, 14 /* Unknown */};
    
    public static final int getIUPACState(boolean maybeA, boolean maybeC, boolean maybeG, boolean maybeT) {
        int index = 0;
        if (maybeA) {
            index += 1;
        }
        if (maybeC) {
            index += 2;
        }
        if (maybeG) {
            index += 4;
        }
        if (maybeT) {
            index += 8;
        }
        return IUPAC_CONV[index];
    }

      /**
     * Conversion of SNP Values to a half byte (4bit) encoding.
     * This uses the standard IUPAC, except gap is 15 or Hex F
     * possible TODO: currently assumes snpValues characters are capital letters
     * mistaking one character for two is not possible as standard ASCII values
     * go up to 126 and lowest two character summation in use is 132
     * @param snpValue
     * @return two bit encoding of base in a byte
     */
    public static final byte getDegerateSNPByteFromTwoSNPs(byte[] snpValue) {
        int valueSum = 0;
        for(int i = 0; i < snpValue.length; i++) {
            valueSum += snpValue[i];
        }
        return getDegerateSNPByteFromTwoSNPs(valueSum);
    }



    private static final byte getDegerateSNPByteFromTwoSNPs(int valueSum) {
        switch (valueSum) {
            case (int)DataType.UNKNOWN_BYTE:
                return DataType.UNKNOWN_BYTE;
            case 65: //'A'
                return 'A';
            case 67: //'C'
                return 'C';
            case 71: //'G'
                return 'G';
            case 84: //'T'
                return 'T';
            case (int)DataType.GAP_BYTE:
                return DataType.GAP_BYTE;
            case 90: //'-' + '-'
                return DataType.GAP_BYTE;
            case 136: //'A+G'
                return 'R';
            case 151: //'C+T'
                return 'Y';
            case 138: //'G+C'
                return 'S';
            case 149: //'A+T'
                return 'W';
            case 155: //'G+T'
                return 'K';
            case 132: //'A+C'
                return 'M';
            case 43: //'+'
                return '+';
            case 86: //'+' + '+'
                return '+';
            case 48: //'0'
                return '0';
            case 96: //'0' + '0'
                return '0';
            case 88: //'+' + '-'
                return '0';
            case 110: //'A' + '-'
                return '0';
            case 112: //'C' + '-'
                return '0';
            case 116: //'G' + '-'
                return '0';
            case 129: //'T' + '-'
                return '0';
            case 130: //'A+A'
                return 'A';
            case 134: //'C+C'
                return 'C';
            case 142: //'G+G'
                return 'G';
            case 168: //'T+T'
                return 'T';
            default:
                return DataType.UNKNOWN_BYTE;
        }
    }
    
    /**
     * Conversion of ASCII byte encoding into ASCII bytes of the diploid.
     *  
     * For example (byte)'A' return {(byte)'A',(byte)'A'}
     * @param halfByte
     * @return character array of length 1 to 3 based on halfByte
     */
    public static byte[] getDiploidValueFromIUPACCode(byte degBase) {
        switch (degBase) {
            case DataType.UNKNOWN_BYTE:
                return NIUPAC;
            case (byte)'A':
                return AIUPAC;
            case (byte)'C':
                return CIUPAC;
            case (byte)'G':
                return GIUPAC;
            case (byte)'T':
                return TIUPAC;
            case (byte)'-':
                return gapIUPAC;
            case (byte)'R':
                return RIUPAC;
            case (byte)'Y':
                return YIUPAC;
            case (byte)'S':
                return SIUPAC;
            case (byte)'W':
                return WIUPAC;
            case (byte)'K':
                return KIUPAC;
            case (byte)'M':
                return MIUPAC;
            case (byte)'+':
                return BIUPAC;
            case (byte)'0':
                return DIUPAC;
            default:
                return NIUPAC;
        }
    }



    public static final byte getDegerateSNPByteFromTwoSNPs(byte snp1, byte snp2) {
        return getDegerateSNPByteFromTwoSNPs(snp1+snp2);
    }

    public String toString() {
        return getDescription();
    }

    public boolean hasGap() {
        return true;
    }

    /**
     * @return true if this character is a '.' or a '_'
     */
    public boolean isGapChar(final char c) {
        if (c == GAP_CHARACTER) {
            return true;
        } else {
            return false;
        }
    }

    /**
     * @return true if state is gap state (-2), false other wise
     */
    public boolean isGapState(final int state) {
        return state == 15;
    }

    public String getFormattedString(char c) {
        return String.valueOf(c);
    }


	@Override
	public boolean isHeterozygote(byte b) {
		for (byte het : HETS_BYTES) if (b==het) return true;
		return false;
	}


	@Override
	public boolean isHeterozygote(char c) {
		for (char het : HETS_CHAR) if (c==het) return true;
		return false;
	}


	@Override
	public int getDiploidIdentity(byte b1, byte b2) {
		return diploidIDTable[getState((char)b1)][getState((char)b2)];
	}


	@Override
	public int getDiploidIdentity(char c1, char c2) {
		return diploidIDTable[getState(c1)][getState(c2)];
	}
    
}
