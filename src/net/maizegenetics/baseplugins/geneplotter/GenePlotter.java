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
public class GenePlotter extends JPanel {

    BorderLayout borderLayout1 = new BorderLayout();
    Alignment aa;
    int[][] exonBoundaries;
    int geneLengthInBases;
    int geneLengthInPixels;
    double xScaleFactor;
    private Color LineColor = Color.blue;
    int scaleinc = 1000, glwidth = 1, numdata;
    float maxX, maxY;
    boolean baseindexflag = true;
    float Xscale, Yscale;
    double[][] numData;
    float threshold;
    char icon = '*';
    boolean iconflag, lineflag = true, thresholdflag;
    String[] numericColumnHeaders;
    String displayString;
    int referenceDataIndex;

    public GenePlotter(Alignment aa, double[][] numData, String[] numericColumnHeaders, int index) {
        this.aa = aa;
        this.numData = numData;
        this.numericColumnHeaders = numericColumnHeaders;
        this.referenceDataIndex = index;
        try {
            jbInit();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    void jbInit() throws Exception {
        this.setLayout(borderLayout1);
    }

    public void isDrawLine(boolean ff) {
        lineflag = ff;
    }

    public void setLineColor(Color theColor) {
        if (theColor != null) {
            LineColor = theColor;
        }
    }

    public void setLineWidth(int gw) {
        glwidth = gw;
    }

    public void isBaseNumbers(boolean ff) {
        baseindexflag = ff;
    }

    public void setBaseNumberInc(int bi) {
        scaleinc = (bi + 1) * 500;
    }

    public void SetNumData(int num) {
        numdata = num;
    }

    public void setThreshold(boolean ff) {
        thresholdflag = ff;
    }

    public void SetThresholdData(float val) {
        threshold = val;
    }

    public void setIcon(boolean ff) {
        iconflag = ff;
    }

    public void setIconStyle(int style) {
        switch (style) {
            case 0:
                icon = 'o';
                break;
            case 1:
                icon = '*';
                break;
            case 2:
                icon = '+';
                break;
            case 3:
                icon = '.';
                break;
            default:
                icon = '*';
                break;
        }
    }

    protected void paintComponent(Graphics g) {
        geneLengthInBases = numData[referenceDataIndex].length;  //aa.getSiteCount();
        /**todo: Override this javax.swing.JComponent method*/
        Dimension d = this.getSize();
        int width = d.width - 30;
        int height = d.height - 30;
        g.setColor(Color.black);
        g.drawLine(0, height, d.width, height);
        g.drawLine(30, 0, 30, d.height);
        maxY = (float) numData[numdata][0];
        maxX = (float) numData[referenceDataIndex][0];

        for (int i = 0; i < numData[referenceDataIndex].length; i++) {
            if (maxX < numData[referenceDataIndex][i]) {
                maxX = (float) numData[referenceDataIndex][i];
            }
        }
        Xscale = (float) (width / maxX);

        for (int i = 0; i < numData[numdata].length; i++) {
            if (maxY < numData[numdata][i]) {
                maxY = (float) numData[numdata][i];
            }
        }
        Yscale = (float) (height / maxY);

        for (int i = 0; i < numData[numdata].length; i++) {
            if ((lineflag) && (i > 0)) {
                g.setColor(LineColor);
                //g.drawLine((int)(i*Xscale)+30, height-(int)(numData[numdata][i]*Yscale), (int)((i+1)*Xscale)+30, height-(int)(numData[numdata][i+1]*Yscale));
                g.drawLine((int) (numData[referenceDataIndex][i - 1] * Xscale) + 30, height - (int) (numData[numdata][i - 1] * Yscale), (int) (numData[referenceDataIndex][i] * Xscale) + 30, height - (int) (numData[numdata][i] * Yscale));
            }
            if (iconflag) {
                g.setColor(Color.blue);
                g.drawString("" + icon, (int) (numData[referenceDataIndex][i] * Xscale) + 30, height - (int) (numData[numdata][i] * Yscale));
            }
            if (thresholdflag) {
                if (numData[numdata][i] <= threshold) {
                    g.setColor(Color.red);
                    g.fillRect((int) (numData[referenceDataIndex][i] * Xscale) + 30, height - 50, 5, 5);
                }
            }
        }
        scaleinc = 25;
        if (baseindexflag) {
            g.setColor(Color.black);
            //Xscale = (float) width/geneLengthInBases;
            for (int i = 0; i < maxX; i += scaleinc) {
                g.drawLine((int) (i * Xscale) + 30, height, (int) (i * Xscale) + 30, height + 10);
                g.drawString("" + i, (int) (i * Xscale) + 30, height + 20);
            //g.drawLine((int)(numData[referenceDataIndex][i]*Xscale)+30, height,(int)(numData[referenceDataIndex][i]*Xscale)+30,height+10);
            //g.drawString(""+(int)(numData[referenceDataIndex][i]),(int)(i*Xscale)+30,height+20);
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