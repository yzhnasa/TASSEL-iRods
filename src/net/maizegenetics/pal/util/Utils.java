// Utils.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)

package net.maizegenetics.pal.util;

import net.maizegenetics.pal.math.MultivariateFunction;
import net.maizegenetics.pal.math.OrthogonalHints;

/**
 * Provides some miscellaneous methods.
 *
 * @version $Id: Utils.java,v 1.1 2007/01/12 03:26:18 tcasstevens Exp $
 *
 * @author Matthew Goode
 */
public class Utils {

	/** Clones an array of doubles
		* @return null if input is null, otherwise return complete copy.
		*/
	public static final double[] getCopy(double[] array) {
		if(array == null) { return null; }
		double[] copy = new double[array.length];
		System.arraycopy(array,0,copy,0,array.length);
		return copy;
	}
	/** Clones an array of doubles from index start (inclusive) to index end (exclusive)
		* @return null if input is null
		*/
	public static final double[] getCopy(double[] array, int start, int end) {
		if(array == null) { return null; }
		double[] copy = new double[end-start];
		System.arraycopy(array,start,copy,0,copy.length);
		return copy;
	}
	/** Clones an array of doubles from index start (inclusive) to end
		* @return null if input is null
		*/
	public static final double[] getCopy(double[] array, int start) {
		return getCopy(array,start,array.length);
	}

	/** Clones an array of bytes
		* @return null if input is null, otherwise return complete copy.
		*/
	public static final byte[] getCopy(byte[] array) {
		if(array == null) {	return null; }
		byte[] copy = new byte[array.length];
		System.arraycopy(array,0,copy,0,array.length);
		return copy;
	}

	/**
	 * Clones an array of doubles
	 * @return null if input is null, otherwise return complete copy.
	 */
	public static final double[][] getCopy(double[][] array) {
		if(array == null) { return null; }
		double[][] copy = new double[array.length][];
		for(int i = 0 ; i < copy.length ; i++) {
			copy[i] = new double[array[i].length];
			System.arraycopy(array[i],0,copy[i],0,array[i].length);
		}
		return copy;
	}
	/**
	 * Clones an array of doubles
	 * @return null if input is null, otherwise return complete copy.
	 */
	public static final double[][][] getCopy(double[][][] array) {
		if(array == null) { return null; }
		double[][][] copy = new double[array.length][][];
		for(int i = 0 ; i < copy.length ; i++) {
			copy[i] = getCopy(array[i]);
		}
		return copy;
	}
	/**
	 * Clones an array of bytes
	 * @return null if input is null, otherwise return complete copy.
	 */
	public static final byte[][] getCopy(byte[][] array) {
		if(array == null) { return null; }
		byte[][] copy = new byte[array.length][];
		for(int i = 0 ; i < copy.length ; i++) {
			copy[i] = new byte[array[i].length];
			System.arraycopy(array[i],0,copy[i],0,array[i].length);
		}
		return copy;
	}
	/**
	 * Clones an array of booleans
	 * @return null if input is null, otherwise return complete copy.
	 */
	public static final boolean[][] getCopy(boolean[][] array) {
		if(array == null) { return null; }
		boolean[][] copy = new boolean[array.length][];
		for(int i = 0 ; i < copy.length ; i++) {
			copy[i] = new boolean[array[i].length];
			System.arraycopy(array[i],0,copy[i],0,array[i].length);
		}
		return copy;
	}
	/**
	 * Clones an array of ints
	 * @return null if input is null, otherwise return complete copy.
	 */
	public static final int[] getCopy(int[] array) {
		if(array == null) { return null; }
		int[] copy = new int[array.length];
		System.arraycopy(array,0,copy,0,array.length);
		return copy;
	}
	/**
	 * Clones an array of ints
	 * @param startingIndex, starts copying from this index
	 * @return null if input is null, otherwise return complete copy.
	 */
	public static final int[] getCopy(int[] array, int startingIndex) {
		if(array == null) { return null; }
		int[] copy = new int[array.length-startingIndex];
		System.arraycopy(array,startingIndex,copy,0,array.length-startingIndex);
		return copy;
	}

	/** Copies all of source into dest - assumes dest to be large enough */
	public static final void copy(double[][] source, double[][] dest) {
		for(int i = 0 ; i < source.length ; i++) {
			System.arraycopy(source[i],0,dest[i],0,source[i].length);
		}
	}
	/**
	 * A simple toString method for an array of doubles.
	 * No fancy formating.
	 * Puts spaces between each value
	 * @param number number of elements to process starting from first element
	 */
	public static final String toString(double[] array, int number) {
		StringBuffer sb = new StringBuffer(array.length*7);
		for(int i = 0 ; i < number ; i++) {
			sb.append(array[i]);
			sb.append(' ');
		}
		return sb.toString();
	}

	/**
	 * A simple toString method for an array of objects.
	 * No fancy formating.
	 * Puts spaces between each value
	 * @param number number of elements to process starting from first element
	 */
	public static final String toString(Object[] array, int number) {
		StringBuffer sb = new StringBuffer(array.length*7);
		for(int i = 0 ; i < number ; i++) {
			sb.append(array[i]);
			sb.append(' ');
		}
		return sb.toString();
	}
	/**
	 * A simple toString method for an array of objects.
	 * No fancy formating.
	 * Puts user defined string between each value
	 * @param number number of elements to process starting from first element
	 */
	public static final String toString(Object[] array, String divider) {
		return toString(array,divider,array.length);
	}
	/**
	 * A simple toString method for an array of objects.
	 * No fancy formating.
	 * Puts user defined string between each value
	 * @param number number of elements to process starting from first element
	 */
	public static final String toString(Object[] array, String divider, int number) {
		StringBuffer sb = new StringBuffer(array.length*7);
		for(int i = 0 ; i < number ; i++) {
			sb.append(array[i]);
			if(i!=number-1) {
				sb.append(divider);
			}
		}
		return sb.toString();
	}
	/**
	 * A simple toString method for an array of doubles.
	 * No fancy formating.
	 * Puts spaces between each value
	 */
	public static final String toString(Object[] array) {
		return toString(array,array.length);
	}
	/**
	 * A simple toString method for an array of doubles.
	 * No fancy formating.
	 * Puts spaces between each value
	 */
	public static final String toString(double[] array) {
		return toString(array,array.length);
	}
	/**
	 * A simple toString method for an array of ints.
	 * No fancy formating.
	 * Puts spaces between each value
	 */
	public static final String toString(int[] array) {
		return toString(array,array.length);
	}
	public static final String toString(int[] array, int number) {
		StringBuffer sb = new StringBuffer(array.length*7);
		for(int i = 0 ; i < number ; i++) {
			sb.append(array[i]);
			sb.append(' ');
		}
		return sb.toString();
	}

	/**
	 * A simple toString method for an array of doubles.
	 * No fancy formating.
	 * Puts spaces between each value
	 */
	public static final String toString(double[][] array) {
		String ss = "";
		for(int i = 0 ; i < array.length ; i++) {
			ss+= i+":"+toString(array[i])+'\n';
		}
		return ss;
	}
	/**
	 * A simple toString method for an array of ints.
	 * No fancy formating.
	 * Puts spaces between each value
	 */
	public static final String toString(int[][] array) {
		String ss = "";
		for(int i = 0 ; i < array.length ; i++) {
			ss+= i+":"+toString(array[i])+'\n';
		}
		return ss;
	}
	/**
	 * Creates an interface between a parameterised object to allow it to act as
	 * a multivariate minimum.
	 */
	public static final MultivariateFunction combineMultivariateFunction(MultivariateFunction base, Parameterized[] additionalParameters) {
		return new CombineMultiParam(base,additionalParameters);
	}

	private static class CombineMultiParam implements MultivariateFunction {
		Parameterized[] additionalParameters_;
		MultivariateFunction base_;

		double[] baseArgumentStorage_;
		double[] lowerBounds_;
		double[] upperBounds_;

		public CombineMultiParam(MultivariateFunction base, Parameterized[] additionalParameters) {
			this.additionalParameters_ = additionalParameters;
			this.base_ = base;
			int numberOfArguments = base_.getNumArguments();
			baseArgumentStorage_ = new double[numberOfArguments];

			for(int i = 0 ; i < additionalParameters.length ; i++) {
				numberOfArguments+=additionalParameters[i].getNumParameters();
			}
			lowerBounds_ = new double[numberOfArguments];
			upperBounds_ = new double[numberOfArguments];
			int argumentNumber = 0;
			for(int i = 0 ; i < base_.getNumArguments() ; i++) {
				lowerBounds_[argumentNumber] = base_.getLowerBound(i);
				upperBounds_[argumentNumber] = base_.getUpperBound(i);
				argumentNumber++;
			}
			for(int i = 0 ; i < additionalParameters.length ; i++) {
				for(int j = 0 ; j < additionalParameters[i].getNumParameters() ; j++) {
					lowerBounds_[argumentNumber] = additionalParameters[i].getLowerLimit(j);
					upperBounds_[argumentNumber] = additionalParameters[i].getUpperLimit(j);
					argumentNumber++;
				}
			}

		}


		public double evaluate(double[] argument) {
			int argumentNumber = 0;
			for(int i = 0 ; i < baseArgumentStorage_.length ; i++) {
				baseArgumentStorage_[i] = argument[argumentNumber];
				argumentNumber++;
			}
			for(int i = 0 ; i < additionalParameters_.length ; i++) {
				int numParam = additionalParameters_[i].getNumParameters();
				for(int j = 0 ; j < numParam ;	j++) {
					additionalParameters_[i].setParameter(argument[argumentNumber], j);
					argumentNumber++;
				}
			}
			return base_.evaluate(baseArgumentStorage_);
		}

		public int getNumArguments() {
			return lowerBounds_.length;
		}
		public double getLowerBound(int n) {
			return lowerBounds_[n];
		}

		public double getUpperBound(int n) {
			return upperBounds_[n];
		}

		/**
		 * @note NEEDS TO BE IMPLEMENTED CORRECTLY
		 * @return null
		 */
		public OrthogonalHints getOrthogonalHints() { return null; }
	}

}

