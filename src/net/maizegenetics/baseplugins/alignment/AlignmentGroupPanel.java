package net.maizegenetics.baseplugins.alignment;

import net.maizegenetics.pal.alignment.Alignment;

import javax.swing.*;
import java.awt.*;

/**
 * Created by IntelliJ IDEA.
 * User: ajf25
 * Date: Jun 3, 2004
 * Time: 10:50:14 AM
 * To change this template use File | Settings | File Templates.
 */
public class AlignmentGroupPanel extends JPanel {

    private AlignmentFrame owner;
    private Sequence[] sequences;
    private int maxSeqLength;

    // indicates if this panel should be drawn for viewing haplotypes or sequences
    private boolean forHaplotypes;
    Alignment aa;

    // offset to space out the dashes and numbers correctly
    private int xOffset;
    FontMetrics theFM;

    public AlignmentGroupPanel(AlignmentFrame owner, Sequence[] seqs, FontMetrics theFM) {
        super();
        this.owner = owner;
        this.sequences = seqs;
        this.theFM = theFM;
        xOffset = theFM.charWidth('G');
        forHaplotypes = false;
        setBackground(Color.white);

        recalculateMaxSeqLength();
        setSize((maxSeqLength) * xOffset, owner.headerOffset);
        Dimension d = getSize();
        setPreferredSize(d);
        setMinimumSize(d);
    }

    public AlignmentGroupPanel(AlignmentFrame owner, Sequence[] seqs, StringBuffer sb, Alignment aa, FontMetrics theFM) {
        super();
        this.aa = aa;
        this.owner = owner;
        this.sequences = seqs;
        this.theFM = theFM;
        xOffset = theFM.charWidth('G');
        forHaplotypes = true;
        setBackground(Color.white);

        recalculateMaxSeqLength();
        //todo:
        setSize((maxSeqLength) * xOffset, owner.headerOffset);
        Dimension d = getSize();
        setPreferredSize(d);
    }

    public void paint(Graphics g) {
        recalculateMaxSeqLength();
        if (forHaplotypes) {
            paintFH(g);
        } else {
            paintFNH(g);
        }
    }

    private void paintFH(Graphics g) {
        setSize((maxSeqLength) * xOffset, owner.headerOffset);
        Dimension dim = getSize();
        setPreferredSize(dim);

        Font font = theFM.getFont();
        g.setFont(font);
        g.setColor(Color.LIGHT_GRAY);
        g.fillRect(0, 0, dim.width, dim.height);

        g.setColor(Color.black);

        Dimension d = getSize();

        String[] headers = getHeaders();
        int yOffset;

        for (int i = 0; i < 5; i++) {
            yOffset = -62 + 15 * i;
            if (yOffset == -2) {
                yOffset = 0;
            }
            g.drawString(headers[i], 0, 0 + d.height + yOffset);
        }

    }

    private void paintFNH(Graphics g) {
        setSize((maxSeqLength) * xOffset, owner.headerOffset);
        Dimension dim = getSize();
        setPreferredSize(dim);

        g.setColor(Color.lightGray);
        g.fillRect(0, 0, dim.width, dim.height);

        g.setColor(Color.black);

        int xindex = 0;
        Dimension d = getSize();
        int evenTen = maxSeqLength + 10 - (maxSeqLength % 10);
        for (int i = 1; i <= evenTen; i++) {
            if (g.hitClip(xindex, 0, xindex + xOffset, d.height)) {
                if (i % 5 == 0) {
                    g.drawString("|", xindex + 3, 0 + d.height);
                }
                // see if we are going to be displaying numbers in the 10-thousands
                if ((i + 4) / 10000 > 0 && ((i + 4) % 10) == 0) {
                    int temp = (i + 4) / 10000;
                    String tenthousands = "" + temp;
                    //System.out.println("marking place " + (i+4) + " and tenthousands is " + temp);    //****
                    g.drawString(tenthousands, xindex, 0 + d.height - 15 - 2);
                } // see if we are going to be displaying numbers in the thousands
                else if ((i + 3) / 1000 > 0 && ((i + 3) % 10) == 0) {
                    int temp = (i + 3) / 10000;
                    temp = (i + 3) - temp * 10000;
                    temp = temp / 1000;
                    String thousands = "" + temp;
                    //System.out.println("marking place " + (i+3) + " and thousands is " + temp);    //****
                    g.drawString(thousands, xindex, 0 + d.height - 15 - 2);
                } // see if we are going to be displaying numbers in the hundreds
                else if ((i + 2) / 100 > 0 && ((i + 2) % 10) == 0) {
                    // remove the tenthousands from the number
                    int tenthou = (i + 2) / 10000;
                    tenthou = (i + 2) - tenthou * 10000;
                    // remove the thousands from the number
                    int temp = (tenthou) / 1000;
                    temp = (tenthou) - temp * 1000;
                    temp = temp / 100;
                    String hundreds = "" + temp;
                    //System.out.println("marking place " + (i+2) + " and hundreds is " + temp);    //****
                    g.drawString(hundreds, xindex, 0 + d.height - 15 - 2);
                } // see if we are going to be displaying numbers in the tens
                else if ((i + 1) / 10 > 0 && ((i + 1) % 10) == 0) {
                    // remove the ten thousands
                    int tenthou = (i + 1) / 10000;
                    tenthou = (i + 1) - tenthou * 10000;
                    // remove the thousands from the number
                    int temp = (tenthou) / 1000;
                    temp = (tenthou) - temp * 1000;
                    // remove the hundreds
                    int temp2 = temp / 100;
                    temp2 = temp - temp2 * 100;
                    temp2 = temp2 / 10;
                    //System.out.println("marking place " + (i+1) + " and tens is " + temp2);    //****
                    String tens = "" + temp2;
                    g.drawString(tens, xindex, 0 + d.height - 15 - 2);
                } // dipslay the 0 in the ones place
                else if (i % 10 == 0) {
                    g.drawString("0", xindex, 0 + d.height - 15 - 2);
                }
            }
            xindex += xOffset;
        }
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
        maxSeqLength = Math.max(maxLen + 7, 150);
    //System.out.println("Newly calculated maxSeqLength is " + maxSeqLength);    //****
    }

    private String[] getHeaders() {
        int sites = aa.getSiteCount();
        char heads[][] = new char[5][sites];
        String[] headerString = new String[5];

        //initialize char to blank
        for (int i = 0; i < 5; i++) {
            for (int j = 0; j < sites; j++) {
                heads[i][j] = ' ';
            }
        }

        //use locus position for headers
        if (aa.getLocusName(0).equals("Unknown")) {
            for (int j = 0; j < sites; j++) {
                String pos = Integer.toString(aa.getPositionInLocus(j));
                int firstChar = 5 - pos.length();
                for (int i = 0; i < pos.length(); i++) {
                    heads[i + firstChar][j] = pos.charAt(i);
                }
            }
        } //use locus name for headers
        else {
            for (int j = 0; j < sites; j++) {
                String name = aa.getLocusName(j);
                int strLength = Math.min(5, name.length());
                for (int i = 0; i < strLength; i++) {
                    heads[i][j] = name.charAt(i);
                }
            }
        }

        for (int i = 0; i < 5; i++) {
            headerString[i] = new String(heads[i]);
        }
        return headerString;
    }
}
