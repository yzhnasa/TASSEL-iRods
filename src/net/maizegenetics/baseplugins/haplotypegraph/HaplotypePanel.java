package net.maizegenetics.baseplugins.haplotypegraph;

import net.maizegenetics.pal.alignment.Alignment;
import org.apache.batik.dom.GenericDOMImplementation;
import org.apache.batik.svggen.SVGGraphics2D;
import org.w3c.dom.DOMImplementation;
import org.w3c.dom.Document;

import javax.swing.*;
import java.awt.*;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.io.File;
import java.io.FileOutputStream;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.util.ArrayList;

/**
 * Created by IntelliJ IDEA.
 * User: ajf25
 * Date: Jul 6, 2004
 * Time: 8:43:08 AM
 * To change this template use File | Settings | File Templates.
 */
public class HaplotypePanel extends JPanel implements MouseListener, MouseMotionListener {

    public static final int MINIMUM_SPANNING_TREE = 1;
    public static final int MINIMUM_SPANNING_TREE_SAME = 2;
    public static final int MINIMUM_SPANNING_TREE_THRESHOLD = 3;
    int drawType;
    int drawThreshold;
    HaplotypeGraph haplotypeGraph;

    public HaplotypePanel(File file) {
        super();
        haplotypeGraph = new HaplotypeGraph(file);
        addMouseListener(this);
        addMouseMotionListener(this);
    }

    public HaplotypePanel(Alignment aa) {
        super();
        haplotypeGraph = new HaplotypeGraph(aa);
        addMouseListener(this);
        addMouseMotionListener(this);
    }

    public void setDrawType(int drawType) {
        if (drawType == MINIMUM_SPANNING_TREE || drawType == MINIMUM_SPANNING_TREE_SAME || drawType == MINIMUM_SPANNING_TREE_THRESHOLD) {
            this.drawType = drawType;
        } else {
            this.drawType = MINIMUM_SPANNING_TREE;
        }
    }

    public void setDrawThreshold(int drawThreshold) {
        this.drawThreshold = drawThreshold;
    }

    public void paint(Graphics g) {
        super.paint(g);
        ((Graphics2D) g).setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        haplotypeGraph.paintGraph(g, drawType, drawThreshold, getSize());
    }

    public void svgSave(File saveFile) {
        if (saveFile == null) {
            return;
        }
        // Get a DOMImplementation
        DOMImplementation domImpl =
                GenericDOMImplementation.getDOMImplementation();

        // Create an instance of org.w3c.dom.Document
        Document document = domImpl.createDocument(null, "svg", null);

        // Create an instance of the SVG Generator
        SVGGraphics2D svgGenerator = new SVGGraphics2D(document);

        // Ask the test to render into the SVG Graphics2D implementation
        //      TestSVGGen test = new TestSVGGen();
        this.paint(svgGenerator);
        //this.paint(svgGenerator);

        // Finally, stream out SVG to the standard output using UTF-8
        // character to byte encoding
        boolean useCSS = true; // we want to use CSS style attribute
        try {
            FileOutputStream fos = new FileOutputStream(saveFile);
            Writer out = new OutputStreamWriter(fos, "UTF-8");
            svgGenerator.stream(out, useCSS);
            fos.flush();
            fos.close();
        } catch (Exception ee) {
            System.out.println("Error with svg button" + ee);
        }
    }

    // implementation of the interfaces
    public void mouseClicked(MouseEvent e) {
        if (e.getModifiers() == Event.META_MASK) {
            // display properties about the vertex just right clicked on
            // get(0) will give sequence, get(1) will give arraylist of names
            ArrayList info = haplotypeGraph.getVertexInfoAt(e.getX(), e.getY());
            if (info != null) {
                JPopupMenu popup = new JPopupMenu();
                popup.setInvoker(this);
                popup.add((String) info.get(0));
                popup.addSeparator();
                ArrayList names = (ArrayList) info.get(1);
                for (int i = 0; i < names.size(); i++) {
                    popup.add((String) names.get(i));
                }
                this.add(popup);
                popup.show(this, e.getX(), e.getY());
                return;
            }
            info = haplotypeGraph.getEdgeInfoAt(e.getX(), e.getY());
            if (info != null) {
                JPopupMenu popup = new JPopupMenu();
                popup.setInvoker(this);
                popup.add((String) info.get(0));
                popup.addSeparator();
                popup.add((String) info.get(1));
                popup.addSeparator();
                popup.add("Difference is " + (String) info.get(2));
                this.add(popup);
                popup.show(this, e.getX(), e.getY());
                return;
            }
        }
    }

    public void mouseEntered(MouseEvent e) {
    }

    public void mouseExited(MouseEvent e) {
    }

    public void mouseReleased(MouseEvent e) {
        haplotypeGraph.clearSelectedVertex();
        repaint();
    }

    public void mousePressed(MouseEvent e) {
        if (e.getModifiers() != Event.META_MASK) {
            haplotypeGraph.selectVertexAt(e.getX(), e.getY());
        }
        repaint();
    }

    public void mouseDragged(MouseEvent e) {
        haplotypeGraph.updateSelectedVertex(e.getX(), e.getY());
        repaint();
    }

    public void mouseMoved(MouseEvent e) {
        repaint();
    }
}