package net.maizegenetics.baseplugins.haplotypegraph;

import edu.uci.ics.jung.graph.Vertex;
import edu.uci.ics.jung.graph.decorators.EdgeStringer;
import edu.uci.ics.jung.graph.decorators.StringLabeller;
import edu.uci.ics.jung.visualization.graphdraw.SettableRenderer;

import java.awt.*;

/**
 * Created using IntelliJ IDEA.
 * Author: Peter Bradbury
 * Date: Aug 17, 2004
 * Time: 4:45:33 PM
 */
public class JungHaplotypeRenderer extends SettableRenderer {
    //fields
    
    //constructors
    public JungHaplotypeRenderer(StringLabeller sl, EdgeStringer el) {
        super(sl, el);
    }
    //methods
    public void paintVertex(Graphics g, Vertex v, int x, int y) {
        String label = getLabel(v);
        //if (mDrawLightBoxes) {
        //    paintLightVertex(g, v, x, y, label);
        //    return;
        //}

        Color fg = (vertexColorFunction == null) ? vertexFGColor
                : vertexColorFunction.getForeColor(v);

        if (vertexColorFunction == null) {
            if (isPicked(v)) {
                g.setColor(vertexPickedColor);
            } else
                g.setColor(vertexBGColor);
        } else {
            g.setColor(vertexColorFunction.getBackColor(v));
        }

        int ht = g.getFontMetrics().getHeight() + 4;
        int descent = g.getFontMetrics().getDescent();
        int width = g.getFontMetrics().stringWidth(label) + 4;
        int top = y - ht/2;
        int left = x - width/2;

        g.fillRect(left , top, width, ht);
        g.setColor(fg);
        g.drawString(label, left + 3, y + ht/2 - descent - 2);
        g.drawRect(left, top, width, ht);

    }
}
