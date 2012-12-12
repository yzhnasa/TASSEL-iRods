/*
 * Locus
 */
package net.maizegenetics.pal.alignment;

import java.io.Serializable;

import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

/**
 *
 * @author terry
 */
public class Locus implements Serializable {

    private static final long serialVersionUID = -5197800047652332969L;
    public static Locus UNKNOWN = new Locus("Unknown", "0", 0, 0, null, null);
    private final String myName;
    private final String myChromosome;
    private final int myStart;
    private final int myEnd;
    private final Map<String, Integer> myFeatures = new HashMap();

    public Locus(String name, String chromosome, int start, int end, String[] featureNames, int[] featurePositions) {

        myName = name;
        myChromosome = chromosome;
        myStart = start;
        myEnd = end;

        if ((featureNames == null) || (featurePositions == null)) {
            return;
        }
        if ((featureNames != null) && (featurePositions != null) && (featureNames.length != featurePositions.length)) {
            throw new IllegalArgumentException("Locus: init: number of feature names should equals number of feature positions.");
        }

        for (int i = 0, n = featureNames.length; i < n; i++) {
            myFeatures.put(featureNames[i], featurePositions[i]);
        }

    }

    public Locus(String name) {
        this(name, name, -1, -1, null, null);
    }

    public String getName() {
        return myName;
    }

    public int getStart() {
        return myStart;
    }

    public int getEnd() {
        return myEnd;
    }

    public String getChromosomeName() {
        return myChromosome;
    }

    public boolean isChromosome() {
        if (myChromosome == null) {
            return false;
        }
        return myChromosome.equals(myName);
    }

    public Map getFeatures() {
        return Collections.unmodifiableMap(myFeatures);
    }

    @Override
    public String toString() {
        return getName();
    }

    @Override
    public boolean equals(Object obj) {

        if (obj == this) {
            return true;
        }

        if (!(obj instanceof Locus)) {
            return false;
        }
        Locus other = (Locus) obj;

        if (!myName.equals(other.getName())) {
            return false;
        }

        if (!myChromosome.equals(other.getChromosomeName())) {
            return false;
        }

        if (myStart != other.getStart()) {
            return false;
        }

        if (myEnd != other.getEnd()) {
            return false;
        }

        Map<String, Integer> otherFeatures = other.getFeatures();
        if (myFeatures.size() != otherFeatures.size()) {
            return false;
        }

        Iterator itr = myFeatures.keySet().iterator();
        while (itr.hasNext()) {
            String key = (String) itr.next();
            Integer value = myFeatures.get(key);
            Integer otherValue = otherFeatures.get(key);
            if ((otherValue == null) && (value != otherValue)) {
                return false;
            }
        }

        return true;

    }
}
