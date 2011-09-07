package net.maizegenetics.baseplugins.haplotypegraph;

import edu.uci.ics.jung.graph.Vertex;
import edu.uci.ics.jung.visualization.GraphMouseListener;
import edu.uci.ics.jung.visualization.VisualizationViewer;

import javax.swing.*;
import java.awt.event.InputEvent;
import java.awt.event.MouseEvent;
import java.util.ArrayList;

/**
 * Created using IntelliJ IDEA.
 * Author: Peter Bradbury
 * Date: Aug 16, 2004
 * Time: 2:41:45 PM
 */
public class JungGraphMouseListener implements GraphMouseListener {
    //fields
    JungHaplotypeGraph jhg;
    VisualizationViewer vv;


    //constructors
    public JungGraphMouseListener(JungHaplotypeGraph jhg, VisualizationViewer vv) {
        this.jhg = jhg;
        this.vv = vv;
    }

    //methods
    public void graphClicked(Vertex v, MouseEvent me) {

    }

    public void graphPressed(Vertex v, MouseEvent me) {
        int btn = me.getModifiersEx();
        if (btn == InputEvent.BUTTON3_DOWN_MASK ||
                btn == InputEvent.BUTTON1_DOWN_MASK + InputEvent.CTRL_DOWN_MASK ||
                btn == InputEvent.BUTTON1_DOWN_MASK + InputEvent.BUTTON3_DOWN_MASK + InputEvent.META_DOWN_MASK ) {

            if (v != null) {
                HaplotypeVertex hv = jhg.getHaplotypeVertexFromJungVertex(v);
                JPopupMenu popup = new JPopupMenu();
                popup.add(hv.getSequence());
                popup.addSeparator();
                ArrayList names = hv.getNames();
                for (int i = 0; i < names.size(); i++)
                    popup.add((String) names.get(i));
                popup.pack();
                popup.show(vv, me.getX(), me.getY());
            }

        }
    }

    public void graphReleased(Vertex v, MouseEvent me) {

    }
}
