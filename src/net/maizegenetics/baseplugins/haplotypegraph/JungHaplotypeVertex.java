package net.maizegenetics.baseplugins.haplotypegraph;

import edu.uci.ics.jung.graph.impl.SparseVertex;

/**
 * Created using IntelliJ IDEA.
 * Author: Peter Bradbury
 * Date: Aug 17, 2004
 * Time: 1:50:27 PM
 */
public class JungHaplotypeVertex extends SparseVertex {
    //fields
    private String toStringKey = null;
    //constructors
    
    //methods
    public String toString() {
        if (toStringKey == null) {
            return super.toString();
        }

        return getUserDatum(toStringKey).toString();
    }

    public void setToStringKey(String key) {
        toStringKey = key;
    }
}
