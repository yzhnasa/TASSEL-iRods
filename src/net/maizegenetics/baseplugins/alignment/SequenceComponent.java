package net.maizegenetics.baseplugins.alignment;

import javax.swing.*;

import java.awt.*;
import java.awt.datatransfer.Clipboard;
import java.awt.datatransfer.ClipboardOwner;
import java.awt.datatransfer.StringSelection;
import java.awt.datatransfer.Transferable;

import java.util.Arrays;
import java.util.LinkedList;

import net.maizegenetics.pal.datatype.DataType;
import net.maizegenetics.prefs.TasselPrefs;

/**
 * User: ajf25
 * Date: Jun 7, 2004
 * Time: 10:30:30 AM
 */
public class SequenceComponent extends JPanel implements ClipboardOwner {

    private Sequence sequence;
    // start and end index will indicate which portion of the
    // sequence will be highlighted
    private int startIndex,  endIndex;
    private int xOffset,  yOffset;
    FontMetrics theFM;

    public SequenceComponent(Sequence sequence, FontMetrics theFM) {
        super();
        this.theFM = theFM;
        xOffset = theFM.charWidth('G');
        yOffset = theFM.getHeight();
        this.sequence = sequence;
        this.startIndex = -1;
        this.endIndex = -1;
        setSize(sequence.getSequence().length() * xOffset, yOffset);
        Dimension d = getSize();
        setPreferredSize(d);
        setMinimumSize(d);
        setMaximumSize(d);

        PreferencesDialog.addSequenceComponent(this);
    }

    public Sequence getSequence() {
        return this.sequence;
    }

    public void setSequence(Sequence seq) {
        this.sequence = seq;
    }

    public int seqLength() {
        return sequence.length();
    }

    public int findColumn(int x, int y) {
        int idx = x / xOffset;
        return idx;
    }

    public void resetHighlight() {
        startIndex = -1;
        endIndex = -1;
    }

    // sets the startIndex for highlighting
    public void startHighlight(int start) {
        startIndex = start;
    }

    // sets the endIndex for highlighting and copies highlight to clipboard
    public void endHighlight(int end) {
        endIndex = end;
        String highlight = "";
        if (startIndex <= endIndex && startIndex != -1) {
            highlight = ((String) sequence.getSequence()).substring(startIndex, endIndex + 1);
        } else if (endIndex < startIndex) {
            highlight = ((String) sequence.getSequence()).substring(endIndex, startIndex + 1);
        }

        StringSelection ss = new StringSelection(highlight);
        Clipboard cb = Toolkit.getDefaultToolkit().getSystemClipboard();
        cb.setContents(ss, this);
    }

    // implementing interface ClipboardOwner
    public void lostOwnership(Clipboard clipboard, Transferable contents) {
    }

    /**
     * shift the sequence of shiftAmt from sequence[shiftStart]
     * insert -'s if shiftAmt is > 0
     * remove -'s if shiftAmt is < 0
     * if shiftAmt is < 0 and there are no -'s remove don't alter sequence
     * @param shiftStart - the index to start shifting from
     * @param shiftAmt - the number of positions to shift the sequence
     */
    public void shiftSequence(int shiftStart, int shiftAmt) {
        // only shift if not viewing matched sequences, shiftAmt non zero, and shiftstart non negative
        if (shiftAmt == 0 || shiftStart < 0) {
            return;
        }
        // inserting -'s
        if (shiftAmt > 0) {
            String seq = sequence.getSequence();

            int[] qualityScores = sequence.getQualityScores();
            Integer[] qualSc = new Integer[qualityScores.length];
            for (int i = 0; i < qualityScores.length; i++) {
                qualSc[i] = new Integer(qualityScores[i]);
            }
            LinkedList ll = new LinkedList(Arrays.asList(qualSc));
            StringBuffer sb = new StringBuffer(seq);
            for (int i = 0; i < shiftAmt; i++) {
                sb.insert(shiftStart, "-");
                ll.add(shiftStart, "-1");
            }

            sequence.setSequence(sb.toString());

            // if/when TASSEL requires JRE 1.5 or greater, this code can be simplified to rely on autoboxing...

            Object lastItem = ll.getLast();

            int length = ll.indexOf(lastItem);
            qualityScores = new int[length + 1000];
            for (int i = 0; i < length; i++) {
                qualityScores[i] = qualSc[i].intValue();
            }

            sequence.setQualityScores(qualityScores);
            return;
        }
        // removing -'s
        if (shiftAmt < 0) {
            String seq = sequence.getSequence();
            StringBuffer sb = new StringBuffer(seq);
            int[] qualityScores = sequence.getQualityScores();
            Integer[] qualSc = new Integer[qualityScores.length];
            for (int i = 0; i < qualityScores.length; i++) {
                qualSc[i] = new Integer(qualityScores[i]);
            }
            LinkedList ll = new LinkedList(Arrays.asList(qualSc));

            for (int i = shiftAmt; i < 0; i++) {
                // if trying to start the shift from the first col return because cannot shift in this direction
                if (shiftStart <= 0) {
                    return;
                }
                // if there is a '-' to remove remove it
                if (sb.charAt(shiftStart - 1) == DataType.GAP_CHARACTER) {
                    sb.deleteCharAt(shiftStart - 1);
                    ll.remove(shiftStart - 1);

                    shiftStart--;
                }
            }
            seq = sb.toString();
            sequence.setSequence(seq);

            // if/when TASSEL requires JRE 1.5 or greater, this code can be simplified to rely on autoboxing...
            int length = ll.size();
            qualityScores = new int[length];
            for (int i = 0; i < length; i++) {
                qualityScores[i] = qualSc[i].intValue();
            }

            sequence.setQualityScores(qualityScores);
            return;
        }
    }

    public void paint(Graphics g) {
        // help indicate if we found an escape character '%' - increment i at the end of the loop
        boolean flag = false;
        Dimension d = getSize();
//       System.out.println("SequenceComponen paintt="+d.getHeight());
        Font font = theFM.getFont();
        g.setFont(font);
        g.setColor(Color.BLACK);

        // index of how far into the panel to start drawing the letter
        int xindex = 0;
        String seq = "";
        int[] qualityScores = {};
        if (sequence.displayQualityScore()) {
//            setQualityScoreValues();
            seq = sequence.getSequence();
        } else {
            if (sequence.displayMatchedSequence()) {
                seq = sequence.getMatchedSequence();
            } else if (sequence.displayColorTopSequence()) {
                seq = sequence.getColorTopSequence();
            } else if (sequence.displayColorConsensus()) {
                seq = sequence.getColorConsensus();
            } else {
                seq = sequence.getSequence();
            }
        }
        for (int i = 0; i < seq.length(); i++) {
            int qualScore = sequence.getQualityScore(i);

            if (seq.charAt(i) == '%') {
                flag = true;
            }
            Rectangle rect = g.getClipBounds();
            double rectWidth = rect.getWidth();
            // added extra set of clipping bounds, because hitClip would return true in cases it really should not have
            // when viewing the end of the sequence, alot of hits would occur for several hundred chars not visible on the screen
            if (g.hitClip(xindex, 0, xindex + xOffset, d.height) && ((xindex + xOffset) > rect.getX()) && (xindex < (rect.getX() + rectWidth))) {
                // take care of the highlighting part, draw black and then white text
                if (startIndex <= endIndex) {
                    if (i >= startIndex && i <= endIndex) {
                        g.setColor(Color.black);
                        g.fillRect(xindex, 0, xOffset, d.height);
                        g.setColor(Color.white);

                        drawSequence(g, seq, i, xindex, d, qualScore);

                        g.setColor(Color.black);
                    } else {
                        drawSequence(g, seq, i, xindex, d, qualScore);
                    }
                } else if (endIndex < startIndex) {
                    if (i >= endIndex && i <= startIndex) {
                        g.setColor(Color.black);
                        g.fillRect(xindex, 0, xOffset, d.height);
                        g.setColor(Color.white);

                        drawSequence(g, seq, i, xindex, d, qualScore);

                        g.setColor(Color.black);
                    } else {
                        drawSequence(g, seq, i, xindex, d, qualScore);
                    }
                } else {
                    drawSequence(g, seq, i, xindex, d, qualScore);
                }
            }
            xindex += xOffset;
            if (flag) {
                i++;
                flag = false;
            }
        }
    }

    public void drawSequence(Graphics g, String seq, int i, int xindex, Dimension d, int qualityScore) {
        if (sequence.displayQualityScore()) {
            Color qualColor = Color.WHITE;

            if (qualityScore < TasselPrefs.getAlignPluginQualscoreLowrangeCutoff()) {
                qualColor = new Color(TasselPrefs.getAlignPluginQualscoreLowrangeColor());
            } else if (qualityScore >= TasselPrefs.getAlignPluginQualscoreHighrangeCutoff()) {
                qualColor = new Color(TasselPrefs.getAlignPluginQualscoreHighrangeColor());
            } else {
                qualColor = new Color(TasselPrefs.getAlignPluginQualscoreMidrangeColor());
            }

            g.setColor(qualColor);
            g.fillRect(xindex, 0, xOffset, d.height);
            g.setColor(Color.black);
            g.drawString(seq.substring(i, i + 1), xindex + 1, 0 + d.height - 3);
        } else {
            if (seq.charAt(i) == DataType.GAP_CHARACTER) {
                g.drawString(seq.substring(i, i + 1), xindex + 2, 0 + d.height - 3);
            } else if (seq.charAt(i) == '.') {
                g.drawString(seq.substring(i, i + 1), xindex + 1, 0 + d.height - 3);
            } // using color sequence '%' means it disagrees with the sequence comparing to
            else if (seq.charAt(i) == '%') {
                i = i + 1;
                // draw the red background
                g.setColor(Color.red);
                g.fillRect(xindex, 0, xOffset, d.height);
                g.setColor(Color.black);
                g.drawString(seq.substring(i, i + 1), xindex, 0 + d.height - 3);
            } else {
                g.drawString(seq.substring(i, i + 1), xindex, 0 + d.height - 3);
            }
        }
    }

    /**
     */
    public void update() {
        repaint();
    }
}
