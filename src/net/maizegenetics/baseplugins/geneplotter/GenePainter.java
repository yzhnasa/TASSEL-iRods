/*
 * TASSEL - Trait Analysis by a aSSociation Evolution & Linkage
 * Copyright (C) 2003 Ed Buckler
 *
 * This software evaluates linkage disequilibrium nucletide diversity and 
 * associations. For more information visit http://www.maizegenetics.net
 *
 * This software is distributed under GNU general public license and without
 * any warranty ot technical support.
 *
 * You can redistribute and/or modify it under the terms of GNU General 
 * public license. 
 *
 */
package net.maizegenetics.baseplugins.geneplotter;

import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.PositionType;

import javax.swing.*;
import java.awt.*;

/**
 * Title:        TASSEL
 * Description:  A java program to deal with diversity
 * Copyright:    Copyright (c) 2000
 * Company:      USDA-ARS/NCSU
 * @author Ed Buckler
 * @version 1.0
 */
public class GenePainter extends JPanel {

    BorderLayout borderLayout1 = new BorderLayout();
    Alignment aa;
    int[][] exonBoundaries;
    int geneLengthInBases;
    int geneLengthInPixels;
    double xScaleFactor;
    float scale = 0;
    Color ExonColor = Color.blue;
    Color LineColor = Color.blue;
    int ScaleInc = 1000;
    int width, exonwidth = 1, linewidth = 1, refaxis;
    boolean baseindexflag = true;

    public GenePainter(Alignment aa) {
        this.aa = aa;
        findBoundaries();
        try {
            jbInit();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
        repaint();
    }

    void jbInit() throws Exception {
        this.setBackground(Color.lightGray);
        this.setLayout(borderLayout1);
    }

    private void findBoundaries() {
        geneLengthInBases = aa.getSiteCount();
        int exoncnt = 0;
        boolean inIntron = true;
        byte baseType;
        for (int i = 0; i < geneLengthInBases; i++) {
            baseType = aa.getPositionType(i);
            if (PositionType.isExon(baseType) && (inIntron == true)) {
                exoncnt++;
                inIntron = false;
            } else if (!PositionType.isExon(baseType)) {
                inIntron = true;
            }
        }
        exonBoundaries = new int[exoncnt][2];
        exoncnt = 0;
        inIntron = true;
        for (int i = 0; i < geneLengthInBases; i++) {
            baseType = aa.getPositionType(i);
            if (PositionType.isExon(baseType) && (inIntron == true)) {
                exonBoundaries[exoncnt][0] = i;
                inIntron = false;
            } else if (!PositionType.isExon(baseType) && (inIntron == false)) {
                exonBoundaries[exoncnt][1] = i - 1;
                exoncnt++;
                inIntron = true;
            } else if (!PositionType.isExon(baseType)) {
                inIntron = true;
            }
        }
        if (!inIntron) // if the gene end in an exon then close things off
        {
            exonBoundaries[exoncnt][1] = geneLengthInBases - 1;
        }
    }

    public void setExonColor(Color theColor) {
        ExonColor = theColor;
    }

    public void setLineWidth(int width) {
        linewidth = width;
    }

    public void setExonWidth(int width) {
        exonwidth = width;
    }

    public void setLineColor(Color theColor) {
        LineColor = theColor;
    }

    public void setBaseNumbers(boolean ff) {
        baseindexflag = ff;
    }

    public void setBaseNumberInc(int bi) {
        ScaleInc = (bi + 1) * 500;
    }

    protected void paintComponent(Graphics g) {
        /**todo: Override this javax.swing.JComponent method*/
        Dimension d = this.getSize();
        int midH = d.height / 2;
        int boxH = d.height / 20;
        g.setColor(LineColor);
        if (linewidth == 1) {
            g.drawLine(30, midH, d.width, midH);
        } else {
            g.fillRect(30, midH - ((linewidth - 1) * 10), d.width, (linewidth - 1) * 20);
        }
        refaxis = midH + linewidth * 15;
        scale = (float) (d.width - 30) / geneLengthInBases;

        g.setColor(ExonColor);
        for (int i = 0; i < exonBoundaries.length; i++) {
            g.fillRect((int) (exonBoundaries[i][0] * scale) + 30, midH - boxH * exonwidth, (int) (exonBoundaries[i][1] * scale) - (int) (exonBoundaries[i][0] * scale), 2 * exonwidth * boxH);
        }
        if (refaxis < midH + boxH * exonwidth) {
            refaxis = midH + boxH * exonwidth;
        }

        if (baseindexflag) {
            g.setColor(Color.black);
            g.drawLine(30, 0, 30, d.height);
            g.drawLine(0, refaxis + 20, d.width, refaxis + 20);
            for (int i = 0; i < geneLengthInBases; i += ScaleInc) {
                g.drawLine((int) (i * scale) + 30, refaxis + 20, (int) (i * scale) + 30, refaxis + 30);
                g.drawString("" + i, (int) (i * scale) + 30, refaxis + 40);
            }
        }
    }

    public void paint(Graphics g) {
        paintComponent(g);
    }

    public Dimension getPreferredSize() {
        /**todo: Override this java.awt.Container method*/
        return new Dimension(500, 200);
    }
}