package net.maizegenetics.util;

import cern.jet.random.engine.MersenneTwister;

import java.util.Arrays;
import java.util.Comparator;

/**
 * Created using IntelliJ IDEA.
 * Author: Peter Bradbury
 * Date: Jun 15, 2004
 * Time: 10:01:34 PM
 */
public class ArrayRandomizer {
    //fields
    MersenneTwister randomizer;
    int[][] index;
    Object[] obj;
    int nObjects;
    Comparator comp;

    //constructors
    public ArrayRandomizer(Object[] objArray) {
        obj = objArray;
        randomizer = new MersenneTwister(new java.util.Date());
        nObjects = obj.length;
        index = new int[nObjects][2];
        for (int i = 0; i < nObjects; i++) {
            index[i][0] = i;
        }

        comp = new Comparator() {
            public int compare(Object o1, Object o2) {
                int[] a = (int[]) o1;
                int[] b = (int[]) o2;
                if (a[1] > b[1]) return 1;
                if (a[1] == b[1]) return 0;
                return -1;
            }
        };
    }

    //methods
    public static void main(String[] args) {
        String[][] strArray = {{"one", "two", "three", "four", "five"},{"one", "two", "three", "four", "five"}};
        Arrays.sort(strArray[0]);
        ArrayRandomizer ar = new ArrayRandomizer(strArray[1]);
        ar.shuffleExistingArray();

        for (int i = 0; i < 5; i++) {
            System.out.println(strArray[0][i] + ", " + strArray[1][i]);
        }
        System.out.println("--------------------------------");
        System.out.println("one, " + strArray[1][Arrays.binarySearch(strArray[0], "one")]);
        System.out.println("two, " + strArray[1][Arrays.binarySearch(strArray[0], "two")]);
        System.out.println("three, " + strArray[1][Arrays.binarySearch(strArray[0], "three")]);
        System.out.println("four, " + strArray[1][Arrays.binarySearch(strArray[0], "four")]);
        System.out.println("five, " + strArray[1][Arrays.binarySearch(strArray[0], "five")]);


    }

    public Object[] getNewRandomArray() {
        Object[] A = new Object[nObjects];
        shuffleIndex();
        for (int i = 0; i < nObjects; i++) {
            A[i] = obj[index[i][0]];
        }
        return A;
    }

    public void shuffleExistingArray() {
        Object[] tempArray = getNewRandomArray();
        for (int i = 0; i < nObjects; i++) {
            obj[i] = tempArray[i];
        }

    }

    private void shuffleIndex() {
        for (int i = 0; i < nObjects; i++) {
            index[i][1] = randomizer.nextInt();
        }
        Arrays.sort(index, comp);
    }

}
