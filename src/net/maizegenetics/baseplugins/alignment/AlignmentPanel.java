package net.maizegenetics.baseplugins.alignment;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;

import java.util.Enumeration;
import java.util.Hashtable;

import net.maizegenetics.pal.datatype.DataType;
import net.maizegenetics.prefs.TasselPrefs;

/**
 * Created by IntelliJ IDEA.
 * User: ajf25
 * Date: Jun 3, 2004
 * Time: 10:49:42 AM
 * To change this template use File | Settings | File Templates.
 */
public class AlignmentPanel extends JPanel implements MouseListener, MouseMotionListener {

    private AlignmentFrame owner;
    private Sequence[] sequences;
    // maxSeqLen helps determine the dimensions of all the sequenceComponent panels
    private int maxSeqLength;
    // used for mouse listening for the selection and dragging of sequences
    private boolean mouse_pressed = false;
    private int mouse_row_pos;
    private int mouse_col_pos;
    // for matching and saving sequences
    // anyonomous action listeners are added to the match and save
    JPopupMenu mainPopupMenu;
    JMenuItem saveMainMenuItem;
    JCheckBoxMenuItem matchCheckBoxMenuItem;
    JCheckBoxMenuItem colorTopRowMenuItem;
    JCheckBoxMenuItem colorConsensusMenuItem;
    JCheckBoxMenuItem colorQualityScoreMenuItem;
    private int xOffset,  yOffset;
    FontMetrics theFM;

    public AlignmentPanel(AlignmentFrame owner, Sequence[] seq, FontMetrics theFM) {
        super();
        this.owner = owner;
        this.sequences = seq;
        this.theFM = theFM;
        xOffset = theFM.charWidth('G');
        yOffset = theFM.getHeight();
        setBackground(Color.white);
        setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));

        recalculateMaxSeqLength();

        // initialize with a few empty jlabels in order to align the names and sequences
//        add(new JLabel(" "));
        for (int i = 0; i < sequences.length; i++) {
            SequenceComponent sc = new SequenceComponent(sequences[i], theFM);
//            Preferences.getInstance().addIObserver(sc);
            add(sc);
        }

        // make all the sequence components the same size
        resizeSeqComponents(maxSeqLength);

        // set up the matching and saving menu for sequences
        mainPopupMenu = new JPopupMenu();
        mainPopupMenu.setInvoker(this);
        saveMainMenuItem = new JMenuItem();
        matchCheckBoxMenuItem = new JCheckBoxMenuItem();
        colorTopRowMenuItem = new JCheckBoxMenuItem();
        colorConsensusMenuItem = new JCheckBoxMenuItem();
        colorQualityScoreMenuItem = new JCheckBoxMenuItem();

        saveMainMenuItem.setText("Save");
        matchCheckBoxMenuItem.setText("Match");
        colorTopRowMenuItem.setText("Color Based on Top Row");
        colorConsensusMenuItem.setText("Color Consensus");
        colorQualityScoreMenuItem.setText("Color Based on Quality Scores");
        colorQualityScoreMenuItem.setSelected(TasselPrefs.getAlignPluginShowQualscore());

        matchCheckBoxMenuItem.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                matchCheckBoxMenuItem_actionPerformed(e);
            }
        });
        colorTopRowMenuItem.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                colorTopRowMenuItem_actionPerformed(e);
            }
        });
        colorConsensusMenuItem.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                colorConsensusMenuItem_actionPerformed(e);
            }
        });

        colorQualityScoreMenuItem.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent e) {
                colorQualityScoreMenuItem_actionPerformed(e);
            }
        });

        mainPopupMenu.add(matchCheckBoxMenuItem);
        mainPopupMenu.add(saveMainMenuItem);
        mainPopupMenu.add(colorTopRowMenuItem);
        mainPopupMenu.add(colorConsensusMenuItem);
        mainPopupMenu.add(colorQualityScoreMenuItem);

        addMouseListener(this);
        addMouseMotionListener(this);
    }

    // mouseListener and mouseMotionListener Interface implementations
    public void mouseClicked(MouseEvent e) {
        // For this example event, we are checking for right-mouse click.
        if (e.getModifiers() == Event.META_MASK) {
            // JPopupMenu must be added to the component whose event is chosen.
            this.add(mainPopupMenu);
            // Make the jPopupMenu visible relative to the current mouse position in the container.
            mainPopupMenu.show(this, e.getX(), e.getY());
        }
    }

    public void mouseEntered(MouseEvent e) {
    }

    public void mouseExited(MouseEvent e) {
    }

    public void mousePressed(MouseEvent e) {
        mouse_row_pos = findRow(e.getX(), e.getY());
        // if valid row grab the component there
        if (mouse_row_pos >= 0) {
            SequenceComponent sc;
            try {
                sc = (SequenceComponent) getComponent(mouse_row_pos);
            } catch (ClassCastException cce) {
                mouse_pressed = false;
                return;
            }
            Point p = sc.getLocation();

            mouse_col_pos = sc.findColumn(e.getX() - p.x, e.getY() - p.y);
        } else {
            mouse_col_pos = -1;
        }
        mouse_pressed = mouse_row_pos >= 0 && mouse_col_pos >= 0;

        // if pressed set the start for highlighting
        if (mouse_pressed && !e.isShiftDown() && !e.isControlDown()) {
            SequenceComponent sc = (SequenceComponent) getComponent(mouse_row_pos);
            sc.startHighlight(mouse_col_pos);
            sc.endHighlight(mouse_col_pos);
        }
    }

    public void mouseReleased(MouseEvent e) {
        mouse_pressed = false;
        resetAllHighlights();
    }

    public void mouseDragged(MouseEvent e) {
        // member of a sequence to realign
        if (mouse_pressed) {
            // get the row we are dragging in
            int idx = findRow(e.getX(), e.getY());

            // dragging of rows
            if (e.isShiftDown()) {
                // if in matching view do not allow the dragging of rows
                // allow dragging if in consensus coloring view
                if (matchCheckBoxMenuItem.isSelected() || colorTopRowMenuItem.isSelected()) {
                    return;
                }
                int rowShift = idx - mouse_row_pos;
                this.owner.swapSequences(mouse_row_pos, mouse_row_pos + rowShift);
                mouse_row_pos = idx;

                // push back the changes into the AnnotationAlignment
                this.owner.updateAnnotationAlignment();
            } // drag within a sequence
            else if (idx >= 1 && idx == mouse_row_pos && e.isControlDown()) {
                SequenceComponent sc = (SequenceComponent) getComponent(mouse_row_pos);
                Point p = sc.getLocation();
                int col = sc.findColumn(e.getX() - p.x, e.getY() - p.y);
                int row = findRow(e.getX(), e.getY());
                // if we have dragged into a new column
                if (mouse_col_pos != col) {
                    int shift = col - mouse_col_pos;

                    // if did not shift any cols, or if dragged out of the row/sequence
                    if (shift == 0 || mouse_row_pos != row) {
                        return;
                    }
                    sc.shiftSequence(mouse_col_pos, shift);

                    recalculateMaxSeqLength();
                    resizeSeqComponents(maxSeqLength);

                    mouse_col_pos = col;

                    // push back the changes into the AnnotationAlignment
                    this.owner.updateAnnotationAlignment();
                }
            } // highlighting text
            else if (idx >= 1 && idx == mouse_row_pos) {
                SequenceComponent sc = (SequenceComponent) getComponent(mouse_row_pos);
                Point p = sc.getLocation();
                int col = sc.findColumn(e.getX() - p.x, e.getY() - p.y);
                sc.endHighlight(col);
            }
        }
        createMatchedSequence();
        createColorTopSequence();
        createColorConsensus();
        repaint();
    }

    public void mouseMoved(MouseEvent e) {
        repaint();
    }

    // return the index of the component located at x,y within the array returned by getComponents()
    private int findRow(int x, int y) {
        Component c = findComponentAt(x, y);
        Component[] comps = getComponents();
        int idx;

        for (idx = comps.length - 1; idx >= 0; idx--) {
            if (comps[idx] == c) {
                //System.out.println(c);    //****
                break;
            }
        }
        return idx;
    }

    private void recalculateMaxSeqLength() {
        int maxLen = -1;
        for (int i = 0; i < sequences.length; i++) {
            Sequence tempSeq = sequences[i];
            int tempLen = tempSeq.length();
            if (tempLen > maxLen) {
                maxLen = tempLen;
            }
        }
        // using max of maxLen and 200 so that the component will always
        // be at least as long as the screen, this will keep the sequences
        // aligned with the column headings
        maxSeqLength = Math.max(maxLen + 7, 200);
    }

    // resize all sequence components to the newSize in characters not pixels
    // mult the newSize by the xOffset to get the new dimension in pixels
    private void resizeSeqComponents(int newSize) {
        maxSeqLength = newSize;
        Component[] comps = getComponents();
        for (int i = 0; i < comps.length; i++) {
            // if the component is a sequence component resize it
            if (comps[i] instanceof SequenceComponent) {
                Dimension d = ((SequenceComponent) comps[i]).getSize();
                ((SequenceComponent) comps[i]).setSize((maxSeqLength) * xOffset, d.height);
                d = comps[i].getSize();
                if (d.width / xOffset != newSize) {
                    System.out.println("Resize failed on component " + i);    //****
                }
                ((SequenceComponent) comps[i]).setPreferredSize(d);
                ((SequenceComponent) comps[i]).setMinimumSize(d);
                ((SequenceComponent) comps[i]).setMaximumSize(d);
            }
        }
    }

    // reset all the highlighting of all the sequences so nothing is highlighted
    public void resetAllHighlights() {
        Component[] comps = getComponents();
        for (int i = 0; i < comps.length; i++) {
            // if the component is a sequence component resize it
            if (comps[i] instanceof SequenceComponent) {
                ((SequenceComponent) comps[i]).resetHighlight();
            }
        }
    }

    public void matchCheckBoxMenuItem_actionPerformed(ActionEvent e) {
        colorTopRowMenuItem.setSelected(false);
        colorConsensusMenuItem.setSelected(false);
        colorQualityScoreMenuItem.setSelected(false);
        setUseQualityScore(false);
        createMatchedSequence();
        createColorTopSequence();
        createColorConsensus();
        repaint();
    }

    public void colorTopRowMenuItem_actionPerformed(ActionEvent e) {
        matchCheckBoxMenuItem.setSelected(false);
        colorConsensusMenuItem.setSelected(false);
        colorQualityScoreMenuItem.setSelected(false);
        setUseQualityScore(false);
        createMatchedSequence();
        createColorTopSequence();
        createColorConsensus();
        repaint();
    }

    public void colorConsensusMenuItem_actionPerformed(ActionEvent e) {
        colorTopRowMenuItem.setSelected(false);
        matchCheckBoxMenuItem.setSelected(false);
        colorQualityScoreMenuItem.setSelected(false);
        setUseQualityScore(false);
        createMatchedSequence();
        createColorTopSequence();
        createColorConsensus();
        repaint();
    }

    public void colorQualityScoreMenuItem_actionPerformed(ActionEvent e) {
        colorTopRowMenuItem.setSelected(false);
        matchCheckBoxMenuItem.setSelected(false);
        colorConsensusMenuItem.setSelected(false);
        setUseQualityScore(colorQualityScoreMenuItem.isSelected());
        repaint();
    }

    private void setUseQualityScore(boolean useQualScore) {
        for (int i = 1; i < sequences.length; i++) {
            sequences[i].setUseQualityScore(useQualScore);
        }

    }
    // will create the matched sequence for all sequences based on the first sequence

    private void createMatchedSequence() {
        // run through all the sequences comparing to the first sequence
        // if a base pair is the same display a "."
        // if different display the base pair
        String baseSeq = sequences[0].getSequence();
        String compareSeq;
        StringBuffer compareMatchSeq;
        for (int i = 1; i < sequences.length; i++) {
            // set whether or not to use the matched sequence
            sequences[i].useMatchedSequence(matchCheckBoxMenuItem.isSelected());
            compareSeq = sequences[i].getSequence();
            compareMatchSeq = new StringBuffer("");

            // run through all the base pairs in baseSeq comparing to compareSeq
            int max = Math.max(baseSeq.length(), compareSeq.length());
            int min = Math.min(baseSeq.length(), compareSeq.length());
            for (int j = 0; j < max; j++) {
                // make sure in the bounds of the the two sequences
                if (j < min) {
                    if (baseSeq.charAt(j) == compareSeq.charAt(j)) {
                        compareMatchSeq.append('.');
                    } else {
                        compareMatchSeq.append(compareSeq.charAt(j));
                    }
                } else if (j < compareSeq.length()) {
                    compareMatchSeq.append(compareSeq.charAt(j));
                }
            }
            sequences[i].setMatchedSequence(compareMatchSeq.toString());
        }
    }

    private void createColorTopSequence() {
        // create the top sequence. if there is a % the character following the % should be colored as it does not match
        String baseSeq = sequences[0].getSequence();
        String compareSeq;
        StringBuffer compareColorTopSeq;
        for (int i = 1; i < sequences.length; i++) {
            // set whether or not to use the colorTopSequence
            sequences[i].useColorTopSequence(colorTopRowMenuItem.isSelected());
            compareSeq = sequences[i].getSequence();
            compareColorTopSeq = new StringBuffer("");

            // run through all the base pairs in baseSeq comparing to compareSeq
            int max = Math.max(baseSeq.length(), compareSeq.length());
            int min = Math.min(baseSeq.length(), compareSeq.length());
            for (int j = 0; j < max; j++) {
                // make sure in the bounds of the two sequences
                if (j < min) {
                    if (baseSeq.charAt(j) == compareSeq.charAt(j)) {
                        compareColorTopSeq.append(baseSeq.charAt(j));
                    } else {
                        compareColorTopSeq.append('%');
                        compareColorTopSeq.append(compareSeq.charAt(j));
                    }
                } else if (j < compareSeq.length()) {
                    compareColorTopSeq.append('%');
                    compareColorTopSeq.append(compareSeq.charAt(j));
                }
            }
            sequences[i].setColorTopSequence(compareColorTopSeq.toString());
        }
    }

    public void createColorConsensus() {
        // first create the consensus sequence
        // then create the color consensus sequence for each sequence
        // if the color consensus sequence will contain a '%' if the sequence differs from the consensus sequence
        String consensus = "";
        StringBuffer consensusSB = new StringBuffer("");
        // table will store (key, Count) where key is a member of the sequence (a character of the sequence)
        Hashtable table = new Hashtable();
        Sequence sequence;
        String seq;
        for (int i = 0; i < maxSeqLength; i++) {
            for (int j = 0; j < sequences.length; j++) {
                sequence = sequences[j];
                seq = sequence.getSequence();
                // make sure this sequence is long enough compared to i
                if (seq.length() > i) {
                    char c = seq.charAt(i);
                    if (table.containsKey(c + "")) {
                        Integer oldCount = (Integer) table.remove(c + "");
                        table.put(c + "", new Integer(oldCount.intValue() + 1));
                    } else {
                        table.put(c + "", new Integer(1));
                    }
                }
            }

            // enumerate the hashtable storing the greatest value and the corresponding key
            Enumeration keys = table.keys();
            String key = "";
            String maxKey = "";
            int count = 0;
            int maxCount = -1;
            while (keys.hasMoreElements()) {
                key = (String) keys.nextElement();
                count = ((Integer) table.get(key)).intValue();
                // see if have a new majority consensus character, but do not allow '-' or '?' to be the consensus
                if (count > maxCount && count > 0 && !key.equals("" + DataType.GAP_CHARACTER) && !key.equals("" + DataType.UNKNOWN_CHARACTER)) {
                    maxCount = count;
                    maxKey = key;
                }
            }
            if (maxCount == -1 || maxKey.equals("-") || maxKey.equals("?")) {
                consensusSB.append('?');
            } else {
                consensusSB.append(maxKey.charAt(0));
            }
            table = new Hashtable();
        }
        // grab the consensus Sequence
        consensus = consensusSB.toString().toUpperCase();

        // run through all sequences and compare to consensus
        String compareSeq;
        StringBuffer compareColorCon;
        for (int i = 0; i < sequences.length; i++) {
            // set whether or not to use the colorTopSequence
            sequences[i].useColorConsensus(colorConsensusMenuItem.isSelected());
            compareSeq = sequences[i].getSequence().toUpperCase();
            compareColorCon = new StringBuffer("");

            // run through all the base pairs in baseSeq comparing to compareSeq
            int max = Math.max(consensus.length(), compareSeq.length());
            int min = Math.min(consensus.length(), compareSeq.length());
            for (int j = 0; j < max; j++) {
                // make sure in the bounds of the two sequences
                if (j < min) {
                    if ((consensus.charAt(j) == '?') || (compareSeq.charAt(j) == 'N') || (compareSeq.charAt(j) == '?')) {
                        // compareColorCon.append('%');
                        compareColorCon.append(compareSeq.charAt(j));
                    } else if (consensus.charAt(j) == compareSeq.charAt(j)) {
                        compareColorCon.append(consensus.charAt(j));
                    } else {
                        compareColorCon.append('%');
                        compareColorCon.append(compareSeq.charAt(j));
                    }
                } else if (j < compareSeq.length()) {
                    compareColorCon.append('%');
                    compareColorCon.append(compareSeq.charAt(j));
                }
            }
            sequences[i].setColorConsensus(compareColorCon.toString());
        }
    }
}
