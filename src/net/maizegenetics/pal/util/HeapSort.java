// HeapSort.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)


package net.maizegenetics.pal.util;


/**
 * sorts numbers and comparable objects by treating contents of array as a binary tree.
 * KNOWN BUGS: There is a horrible amount of code duplication here!
 *
 * @version $Id: HeapSort.java,v 1.2 2009/06/09 18:46:00 tcasstevens Exp $
 *
 * @author Alexei Drummond
 * @author Korbinian Strimmer
 */
public class HeapSort {

	//
	// Public stuff
	//

	/**
	 * Sorts an array of indices into an array of doubles
	 * into increasing order.
	 */
	public static void sort(double[] array, int[] indices)
	{

		// ensures we are starting with valid indices
		for (int i = 0; i < indices.length; i++)
		{
			indices[i] = i;
		}

		int temp;
		int j, n = array.length;

		// turn input array into a heap
		for (j = n/2; j > 0; j--) {
			adjust(array, indices, j, n);
		}

		// remove largest elements and put them at the end
		// of the unsorted region until you are finished
		for (j = n-1; j > 0; j--) {
			temp = indices[0];
			indices[0] = indices[j];
			indices[j] = temp;
			adjust(array, indices, 1, j);
		}
	}


	/**
	 * helps sort an array of indices into an array of doubles.
	 * Assumes that array[lower+1] through to array[upper] is
	 * already in heap form and then puts array[lower] to
	 * array[upper] in heap form.
	 */
	private static void adjust(double[] array, int[] indices, int lower, int upper) {

		int j, k;
		int temp;

		j = lower;
		k = lower * 2;

		while (k <= upper)
		{
			if ((k < upper) && (array[indices[k-1]] < array[indices[k]]))
			{
				k += 1;
			}
			if (array[indices[j-1]] < array[indices[k-1]])
			{
				temp = indices[j-1];
				indices[j-1] = indices[k-1];
				indices[k-1] = temp;
			}
			j = k;
			k *= 2;
		}
	}
}

