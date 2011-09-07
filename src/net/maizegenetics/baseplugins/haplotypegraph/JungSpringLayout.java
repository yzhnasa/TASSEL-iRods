package net.maizegenetics.baseplugins.haplotypegraph;

import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.Vertex;
import edu.uci.ics.jung.visualization.SpringLayout;

import java.util.Iterator;

/**
 * Created using IntelliJ IDEA.
 * Author: Peter Bradbury
 * Date: Aug 16, 2004
 * Time: 11:13:01 AM
 */
public class JungSpringLayout extends SpringLayout {
    //fields
   
    //constructors
    JungSpringLayout(Graph g, LengthFunction f) {
        super(g,f);
    }

    //methods

    /**
     * returns the vertex for which x,y falls within the actual vertex symbol, null otherwise
     * @param x
     * @param y
     * @return
     */
    public Vertex getVertex(double x, double y) {
        //these should match the values for vertex size in the renderer (SettableRenderer)
        int maxXdist = 7;
        int maxYdist = 5;

        Vertex v = null;
        Iterator iter = getVisibleVertices().iterator();

        while (iter.hasNext()) {
            v = (Vertex) iter.next();
            if (Math.abs(getX(v) - x) <= maxXdist &&  Math.abs(getY(v) - y) <= maxYdist) {
                return v;
            }
        }

        return null;

    }
}
