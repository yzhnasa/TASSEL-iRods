package net.maizegenetics.baseplugins.alignment;

import javax.swing.*;
import java.awt.*;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;

/**
 * Created by IntelliJ IDEA.
 * User: ajf25
 * Date: Jun 3, 2004
 * Time: 10:19:15 AM
 * To change this template use File | Settings | File Templates.
 */
public class AlignmentNamePanel extends JPanel implements MouseListener
{
    protected AlignmentFrame owner;
    
    // vars used for MouseListener Interface
    private boolean mouse_pressed = false;
    private int mouse_row_pos;
    private int mouse_ann_pos;
    private int last_selected_row = -1;
    private JPanel spacePanel;
    private int xOffset, yOffset, maxNameLengthInChar;
    FontMetrics theFM;

    public AlignmentNamePanel(AlignmentFrame owner, Sequence[] seq, FontMetrics theFM)
    {
        super();
        this.owner = owner;
        this.theFM=theFM;
        xOffset=theFM.charWidth('G');
        yOffset=theFM.getHeight();
        maxNameLengthInChar=getMaxNameLength(seq);
        setBackground(Color.lightGray);
        setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));

        // initialize with a few empty jlabels to correctly align names with their sequences
/*        add(new JLabel(" "));
        add(new JLabel(" "));
        add(new JLabel(" "));
  */
        spacePanel=new JPanel();
        add(spacePanel);
        Dimension d=new Dimension((maxNameLengthInChar+2)*xOffset,owner.headerOffset);
        spacePanel.setPreferredSize(d);
        spacePanel.setMinimumSize(d);
        spacePanel.setMaximumSize(d);

        for(int i = 0; i < seq.length; i++)
        {
            add(new SequenceNameComponent(seq[i], theFM, maxNameLengthInChar));
        }
        addMouseListener(this);
    }
    private int getMaxNameLength(Sequence[] seqs) {
        int maxLength=0;
        for(int i = 0; i < seqs.length; i++)
        {   if(seqs[i].getName().length()>maxLength) {maxLength=seqs[i].getName().length();}}
        return maxLength;
    }
    public AlignmentNamePanel(AlignmentFrame owner, Sequence[] seq, StringBuffer sb, FontMetrics theFM)
    {
        super();
        this.owner = owner;
        this.theFM=theFM;
        xOffset=theFM.charWidth('G');
        yOffset=theFM.getHeight();
        maxNameLengthInChar=getMaxNameLength(seq);
        setBackground(Color.white);
        setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));

        spacePanel=new JPanel();
        add(spacePanel);
        Dimension d=new Dimension((maxNameLengthInChar+2)*xOffset,owner.headerOffset);
        spacePanel.setPreferredSize(d);
        spacePanel.setMinimumSize(d);
        spacePanel.setMaximumSize(d);

/*        // these labels are used to align the names with the sequences when there is extra heading information
        add(new JLabel(" "));
        add(new JLabel(" "));
        add(new JLabel(" "));

        // initialize with a few empty jlabels to correctly align names with their sequences
        add(new JLabel(" "));
        add(new JLabel(" "));
        add(new JLabel(" "));   */
        for(int i = 0; i < seq.length; i++)
        {
            add(new SequenceNameComponent(seq[i], theFM, maxNameLengthInChar));
        }
        addMouseListener(this);
    }

    public void paint(Graphics g) {
        super.paint(g);    //To change body of overridden methods use File | Settings | File Templates
/*        Font font = theFM.getFont();
        FontMetrics fm=g.getFontMetrics(font);
        int h=fm.getHeight();
        Dimension d=new Dimension(10,5*h);
        spacePanel.setPreferredSize(d);
        spacePanel.setMinimumSize(d);
        spacePanel.setMaximumSize(d);       */
    }

    public Dimension getMinimumSize()
    {
        return getPreferredSize();
    }

    public Dimension getMaximumSize()
    {
        return getPreferredSize();
    }

    public Dimension getPreferredSize()
    {
        return getLayout().minimumLayoutSize(this);
    }

    // MouseListener Interface
    public void mouseClicked(MouseEvent e) {}
    public void mouseEntered(MouseEvent e) {}
    public void mouseExited(MouseEvent e) {}

    public void mousePressed(MouseEvent e)
    {
        mouse_row_pos = findRow(e.getX(), e.getY());
        mouse_pressed = mouse_row_pos >= 0;
        if (!mouse_pressed)
            return;
    }

    public void mouseReleased(MouseEvent e) {}

    private int findRow(int x, int y)
    {
        // return the row located at coordinates x,y
        Component c = findComponentAt(x, y);
        Component[] comps = getComponents();
        int idx;

        for (idx = comps.length - 1; idx >= 0; idx--)
        {
            if (comps[idx] == c)
                break;
        }
        return idx;
    }

}
