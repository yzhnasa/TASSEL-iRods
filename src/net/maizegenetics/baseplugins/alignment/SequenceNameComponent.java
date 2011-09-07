package net.maizegenetics.baseplugins.alignment;


import javax.swing.*;
import java.awt.*;

//import com.neogenesis.pfaat.colorscheme.ColorScheme;


/**
 * <code>Component</code> for displaying a single sequence name with
 * line annotation names..
 *
 * @author $Author: tcasstevens $
 * @version $Revision: 1.1 $, $Date: 2007/08/07 21:13:05 $ */
public class SequenceNameComponent extends JPanel implements SequenceListener
{
    // underlying sequence
    private Sequence sequence = null;
    private int xOffset, yOffset;
    FontMetrics theFM;

    public SequenceNameComponent(Sequence sequence, FontMetrics theFM, int maxNameLengthInChar)
    {
        super();
        this.theFM=theFM;
        xOffset=theFM.charWidth('G');
        yOffset=theFM.getHeight();
 //       Font font = theFM.getFont();
 //       Font font = theFM.getFont();
        setFont(theFM.getFont());
        setBackground(Color.lightGray);
        setSequence(sequence);
        Dimension d = new Dimension((maxNameLengthInChar+2)*xOffset,yOffset);
        setPreferredSize(d);
        setMinimumSize(d);     //
        setMaximumSize(d);   //
        setSize(d);
//        repaint();
    }

    // set the underlying sequence
    public void setSequence(Sequence sequence)
    {
        this.sequence = sequence;
        revalidate();
        repaint();
    }

    public Sequence getSequence()
    {
        return this.sequence;
    }

    // SequenceListener interface
    public void sequenceAAChanged(Sequence aaseq) {}

    public void sequenceNameChanged(Sequence aaseq, String old_name)
    {
        if (aaseq != sequence)
            throw new RuntimeException("bound to incorrect Sequence");
        revalidate();
        repaint();
    }

    public void sequenceAnnotationChanged(Sequence aaseq) {}
    public void sequenceGroupChanged(Sequence aaseq) {}
    public void sequenceColumnAnnotationsChanged(Sequence aaseq, int column) {}

    public void sequenceLineAnnotationsChanged(Sequence aaseq)
    {
        if (aaseq != sequence)
            throw new RuntimeException("bound to incorrect Sequence");
    }

    public void sequenceColorChanged(Sequence aaseq)
    {
        if (aaseq != sequence)
        {
            throw new RuntimeException("bound to incorrect Sequence");
        }
        revalidate();
        repaint();
    }

    // DisplayPropertiesListener interface
    public void displayAnnViewChanged(Sequence seq, boolean show)
    {
        if (seq == sequence)
        {
            revalidate();
            repaint();
        }
    }

    public void displaySeqSelectChanged(Sequence seq, boolean select)
    {
        if (seq == sequence)
        {
            revalidate();
            repaint();
        }
    }
/*
    // default size
    public Dimension getMinimumSize()
    {
        return getPreferredSize();
    }

    public Dimension getMaximumSize()
    {
        Dimension d = getPreferredSize();
        d.width = Integer.MAX_VALUE;
        return d;
    }
*/
    public void paint(Graphics g)
    {
         Dimension d = getSize();
        super.paint(g);
        d = getSize();         //ed
 //       System.out.println("SequenceNameComponent="+d.getHeight());     //ed

        Font font = theFM.getFont();
        g.setFont(font);
        g.setColor(Color.BLACK);
        g.drawString(sequence.getName(), 0, 0+d.height-3);
    }

}
