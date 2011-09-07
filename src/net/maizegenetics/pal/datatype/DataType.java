// DataType.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)
// Known bugs and limitations:
// - all states must have a non-negative value 0..getNumStates()-1
// - ? (unknown state) has value getNumStates()
package net.maizegenetics.pal.datatype;

import java.io.Serializable;

/**
 * interface for sequence data types
 * History: 21 March 2003, Added gap stuff, to counter frustration and not being
 * able to differentiat unknowns from gaps. Gap characters should still be treated
 * as unknowns (for compatibility), but a data type should be able to identify
 * a gap from other unknowns.
 *
 * @version $Id: DataType.java,v 1.1 2007/01/12 03:26:13 tcasstevens Exp $
 *
 * @author Korbinian Strimmer
 * @author Alexei Drummond
 */
public interface DataType extends Serializable {

    public static final char UNKNOWN_CHARACTER = 'N';
    public static final byte UNKNOWN_BYTE = (byte) UNKNOWN_CHARACTER;
    public static final char GAP_CHARACTER = '-';
    public static final byte GAP_BYTE = (byte) '-';

    /**
     * Get number of unique states
     *
     * @return number of unique states
     */
    public int getNumStates();

    /**
     * Get state corresponding to a character
     *
     * @param c character
     *
     * @return state
     */
    public int getState(char c);

    /**
     * Get character corresponding to a given state
     *
     * @param state state
     *
     * return corresponding character
     */
    public char getChar(int state);

    /**
     * Get the  String version of a particular character (eg a -> A)
     * This is important for datatypes best presented by multiple characters
     * eg. TextDataType
     */
    public String getFormattedString(char c);

    /**
     * Description of data type
     *
     * @return string describing the data type
     */
    public String getDescription();

    /**
     * @return true if this state is an unknown state
     * (the same as check if a state is >= the number of states... but neater)
     */
    public boolean isUnknownState(int state);

    /**
     * @return true if this character is a gap
     */
    public boolean isUnknownChar(char c);

    /**
     * @return true if this data type supports having a gap character
     */
    public boolean hasGap();

    /**
     * @return true if this data type interprets c as a gap
     */
    public boolean isGapChar(char c);

    /**
     * @return true if this data type interprets state as a gap state
     */
    public boolean isGapState(int state);

    /**
     * For dealing with readable datatype with different string lengths.  This is most important
     * for Numeric and TextDataType
     */
    public void setMaxFormatedStringLength(int length);

    /**
     * For dealing with readable datatype with different string lengths.  This is most important
     * for Numeric and TextDataType
     */
    public int getMaxFormatedStringLength();
    
    /**
     * @param b	a genotype for a site
     * @return	true if the genotype is heterozygous, false otherwise
     */
    public boolean isHeterozygote(byte b);
    
    /**
     * @param c	a genotype for a site
     * @return	true if the genotype is heterozygous, false otherwise
     */
    public boolean isHeterozygote(char c);
    
    /**
     * Counts the number of pairwise allele matches. For identical homozygotes, this equals n * n, where n is the ploidy level. 
     * The probability that alleles drawn from the two loci are identical in state, P(IBS), equals this count divided by n*n.
     * If either genotype is unknown returns -1. 
     * @param b1	the genotype for the first site
     * @param b2	the genotype for the second site
     * @return	a count of the number of times pairs of alleles match examining all possible pairwise combinations
     */
    public int getDiploidIdentity(byte b1, byte b2);
    
    /**
     * Counts the number of pairwise allele matches. For identical homozygotes, this equals n * n, where n is the ploidy level. 
     * The probability that alleles drawn from the two loci are identical in state, P(IBS), equals this count divided by n*n. 
     * If either genotype is unknown returns -1. 
     * @param c1	the genotype for the first site
     * @param c2	the genotype for the second site
     * @return	P(IBS) scaled to return an integer
     */
    public int getDiploidIdentity(char c1, char c2);
}
