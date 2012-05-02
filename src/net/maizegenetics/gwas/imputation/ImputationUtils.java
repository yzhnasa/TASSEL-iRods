package net.maizegenetics.gwas.imputation;

import java.util.Arrays;
import java.util.LinkedList;

public class ImputationUtils {

	public static int[] order(int[] array) {
		class SortElement implements Comparable<SortElement> {
			int val;
			int ndx;
			
			SortElement(int x, int index) {
				val = x;
				ndx = index;
			}
			
			@Override
			public int compareTo(SortElement se) {
				return val - se.val;
			}
		}

		int n = array.length;
		SortElement[] sortArray = new SortElement[n];
		for (int i = 0; i < n; i++) {
			sortArray[i] = new SortElement(array[i], i);
		}
		
		Arrays.sort(sortArray);
		int[] order = new int[n];
		for (int i = 0; i < n; i++) {
			order[i] = sortArray[i].ndx;
		}
		return order;
	}
	
	public static int[] reverseOrder(int[] array) {
		class SortElement implements Comparable<SortElement> {
			int val;
			int ndx;
			
			SortElement(int x, int index) {
				val = x;
				ndx = index;
			}
			
			@Override
			public int compareTo(SortElement se) {
				return se.val - val;
			}
		}

		int n = array.length;
		SortElement[] sortArray = new SortElement[n];
		for (int i = 0; i < n; i++) {
			sortArray[i] = new SortElement(array[i], i);
		}
		
		Arrays.sort(sortArray);
		int[] order = new int[n];
		for (int i = 0; i < n; i++) {
			order[i] = sortArray[i].ndx;
		}
		return order;
	}
}
