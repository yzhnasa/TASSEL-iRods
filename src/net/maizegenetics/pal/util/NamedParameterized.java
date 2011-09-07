// Parameterized.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)


package net.maizegenetics.pal.util;




/**
 * interface for class with (optimizable) named parameters
 *
 * @version $Id: NamedParameterized.java,v 1.1 2007/01/12 03:26:18 tcasstevens Exp $
 *
 * @author Alexei Drummond
 */
public interface NamedParameterized extends Parameterized {
	
	/**
	 * @return a short identifier for this parameter type. Should be the same for 
	 * all instances of a given class!
	 */
	String getParameterName(int i);
}
