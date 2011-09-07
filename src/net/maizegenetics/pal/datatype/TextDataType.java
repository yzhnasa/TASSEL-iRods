package net.maizegenetics.pal.datatype;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;

/**
 * User: ed
 * Date: Feb 3, 2005
 * Time: 1:01:55 PM
 */
public class TextDataType implements DataType {

    private static final long serialVersionUID = 7902613264354545217L;
    public static final TextDataType DEFAULT_INSTANCE = new TextDataType();
    public static final int UNKNOWN_STATE = UNKNOWN_CHARACTER - 65;
    public static final String UNKNOWN_STRING = "N";
    private final ArrayList<String> myStateLookup;
    private int maxStringLength = 0;
    
    /**
     * Empty constructor
     */
    public TextDataType() {
        myStateLookup = new ArrayList<String>();
    }

    public TextDataType(Collection<String> states) {
        myStateLookup = new ArrayList<String>(states);
        Iterator<String> sit = states.iterator();
        while (sit.hasNext()) maxStringLength = Math.max(maxStringLength, sit.next().length());
    }

    public TextDataType(String[] states) {

        myStateLookup = new ArrayList<String>(states.length);
        for (int i = 0; i < states.length; i++) {
            myStateLookup.add(states[i]);
            maxStringLength = Math.max(maxStringLength, states[i].length());
        }
    }

    @Override
	public boolean isUnknownState(int state) {
    	return state == UNKNOWN_STATE;
	}

	@Override
	public boolean isUnknownChar(char c) {
		return c == UNKNOWN_CHARACTER;
	}

	@Override
	public void setMaxFormatedStringLength(int length) {
		maxStringLength = length;
	}

	@Override
	public int getMaxFormatedStringLength() {
		return maxStringLength;
	}

	public void addStates(Collection<String> states) {
        if (states != null) {
            myStateLookup.addAll(states);
        }
    }

    // Get number of bases
    public int getNumStates() {
        return myStateLookup.size();
    }

    public String getFormattedString(char c) {
    	return getTextRepresentationFromState(getState(c));
    }

    /**
     * Get state corresponding to a character.
     * 
     * @param c character
     * @return state
     */
    public int getState(final char c) {
        if (c == UNKNOWN_CHARACTER) {
            return UNKNOWN_STATE;
        }
        return c - 65;
    }

    /**
     * Get character corresponding to a given state
     */
    public char getChar(final int state) {
        if (state >= myStateLookup.size() || state < 0) {
            return UNKNOWN_CHARACTER;
        }
        return (char) (state + 65);
    }

    /**
     * Returns a unique ascii character for any given numeric size
     */
    public final char getCharFromTextRepresentation(String s) {
    	return getChar(getStateFromTextRepresentation(s));
    }

    /**
     * Returns a unique ascii character for any given numeric size
     */
    public final int getStateFromTextRepresentation(String s) {
    	
        int index = myStateLookup.indexOf(s);
        if (index < 0) {
            //this next lines prevents UNKNOWN_CHARACTER from being used
            if (myStateLookup.size() == UNKNOWN_STATE) {
                myStateLookup.add(UNKNOWN_STRING);
            }
        	
            myStateLookup.add(s);
            index = myStateLookup.indexOf(s);
            if (s.length() > getMaxFormatedStringLength()) {
                setMaxFormatedStringLength(s.length());
            }
        }
        return (index);
    }

    /**
     * Returns the actual string/sequence by location within the data structure
     */
    public final String getTextRepresentationFromState(int index) {
    	if (index==UNKNOWN_STATE) return "N";
        return (String) myStateLookup.get(index);
    }

    /** String describing the data type */
    public String getDescription() {
        return "Text";
    }

    public boolean hasGap() {
        return false;
    }

    public boolean isGapChar(final char c) {
        throw new UnsupportedOperationException();
    }

    public boolean isGapState(final int state) {
        throw new UnsupportedOperationException();
    }

	@Override
	public boolean isHeterozygote(byte state) {
		//this implementation assumes that b is the state not its char representation
		if (isUnknownState(state)) return false;
		String[] genotype = myStateLookup.get(state).split(":");
		int n = genotype.length;
		for (int i = 1; i < n; i++) if (!genotype[i].equals(genotype[0])) return true;
		return false;
	}

	@Override
	public boolean isHeterozygote(char c) {
		if (isUnknownChar(c)) return false;
		String[] genotype = getFormattedString(c).split(":");
		int n = genotype.length;
		for (int i = 1; i < n; i++) if (!genotype[i].equals(genotype[0])) return true;
		return false;
	}

	@Override
	public int getDiploidIdentity(byte state1, byte state2) {
		//this implementation assumes that b is the state not its char representation
		if (isUnknownState(state1) || isUnknownState(state2)) return -1;
		String[] genotype1 = myStateLookup.get(state1).split(":");
		String[] genotype2 = myStateLookup.get(state2).split(":");
		int n = genotype1.length;
		int matches = 0;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (genotype1[i].equals(genotype2[j])) matches++;
			}
		}
		if (n == 1) return 4*matches;
		return matches;
	}

	@Override
	public int getDiploidIdentity(char c1, char c2) {
		if (isUnknownChar(c1) || isUnknownChar(c2)) return -1;
		String[] genotype1 = getFormattedString(c1).split(":");
		String[] genotype2 = getFormattedString(c2).split(":");
		int n = genotype1.length;
		int matches = 0;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (genotype1[i].equals(genotype2[j])) matches++;
			}
		}
		if (n == 1) return 4*matches;
		return matches;
	}
}
