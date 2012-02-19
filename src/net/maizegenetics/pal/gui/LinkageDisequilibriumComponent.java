// LinkageDisequilibriumComponent.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)
package net.maizegenetics.pal.gui;

import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.popgen.LinkageDisequilibrium;

import java.awt.*;
import javax.swing.JComponent;
import java.text.DecimalFormat;

/**
 * An AWT Component for displaying information on linkage disequilibrium.
 *
 * Nice schematics are produced if an annotation alignment is used to construct
 * LinkageDisequilibrium.  It can portray things both on the gene and chromosomal
 * scale.
 *
 *
 * @author Ed Buckler
 * @version $Id: LinkageDisequilibriumComponent.java
 */
public class LinkageDisequilibriumComponent extends JComponent {

    public final static int P_VALUE = 0;
    public final static int DPRIME = 1;
    public final static int RSQUARE = 2;
    double minimumChromosomeLength = 10;
    LinkageDisequilibrium theLD;
    Alignment theAA;
    boolean includeBlockSchematic, chromosomalScale;
    boolean includeLabels = true;
    int totalVariableSites, totalLoci, totalChromosomes, totalIntervals, totalBlocks;
    double[] startPos, endPos; //These are the relative positions of the polymorphisms
    double[] blockBeginPos, blockEndPos;
    String[] blockNames;
    int[] xPos, yPos, xEndPos;  //these hold positions of the upper left corners for each site
    int[] blockBeginX, blockEndX;//These are the absolute positions of the genes & chromosomes
    int ih, iw;
    double totalUnits;
    double[] blockStart, blockEnd;
    //this will range from 0 to 1
    String upperLabel, lowerLabel;
    double[][] diseq;
    Color theColor = new Color(0, 0, 0);
    int distanceBetweenGraphAndGene = 40;
    int hoff = 70, h2off = 70, voff = 20;
    boolean probability = true, upperProb = false, lowerProb = true;
    double normalizer;
    //viewer attribute variables
    int myWindowSize;
    int myWindowX;
    int myWindowY;

    public LinkageDisequilibriumComponent(LinkageDisequilibrium theLD, boolean includeBlockSchematic, boolean chromosomalScale, int windowSize, int windowX, int windowY) {
        this.theLD = theLD;
        theAA = theLD.getAlignment();
        this.includeBlockSchematic = includeBlockSchematic;
        this.chromosomalScale = chromosomalScale;
        int numSites = theLD.getSiteCount();
        myWindowSize = windowSize;
        myWindowX = windowX;
        myWindowY = windowY;
        this.diseq = new double[windowSize][windowSize];
        for (int x = 0; x < myWindowSize; x++) {
            for (int y = 0; y < myWindowSize; y++) {
                diseq[x][y] = Double.NaN;
            }
        }
        setUpperCorner(RSQUARE);
        setLowerCorner(P_VALUE);
        totalVariableSites = theLD.getSiteCount();
        if (theAA != null) {
            countGenesAndChromosomes();
            calculateStartAndEndPositions();
        } else {
            includeBlockSchematic = false;
        }
        xPos = new int[windowSize + 1];
        yPos = new int[windowSize + 1];
        xEndPos = new int[windowSize + 1];
        try {
            jbInit();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    /**
     * This sets a new viewable size
     */
    public void setWindowSize(int newSize, int ldMeasureLower, int ldMeasureUpper) {
        System.out.println("new window size: " + newSize);
        myWindowSize = newSize;
        diseq = new double[myWindowSize][myWindowSize];

        setLowerCorner(ldMeasureLower);
        setUpperCorner(ldMeasureUpper);

        if (theAA != null) {
            countGenesAndChromosomes();
            calculateStartAndEndPositions();
        }

        xPos = new int[myWindowSize + 1];
        yPos = new int[myWindowSize + 1];
        xEndPos = new int[myWindowSize + 1];
    }

    /**
     * This sets a new X position
     */
    public void setWindowX(int newX, int ldMeasureLower, int ldMeasureUpper) {
        myWindowX = newX;
        calculateStartAndEndPositions();
        setLowerCorner(ldMeasureLower);
        setUpperCorner(ldMeasureUpper);
    }

    /**
     * This sets a new Y position
     */
    public void setWindowY(int newY, int ldMeasureLower, int ldMeasureUpper) {
        myWindowY = newY;
        calculateStartAndEndPositions();
        setLowerCorner(ldMeasureLower);
        setUpperCorner(ldMeasureUpper);
    }

    public void setLowerCorner(int ldMeasure) {

        for (int r = 0; r < myWindowSize; r++) {
            for (int c = Math.max(r + myWindowX - myWindowY, 0); c < myWindowSize; c++) {
                diseq[r][c] = Double.NaN;
            }
        }

        switch (ldMeasure) {
            case P_VALUE: {
                for (int z = 0, n = theLD.getRowCount(); z < n; z++) {
                    int r = theLD.getX(z) - myWindowX + (myWindowSize / 2);
                    int c = theLD.getY(z) - myWindowY + (myWindowSize / 2);
                    if ((r >= 0) && (r < myWindowSize) && (c >= Math.max(r + myWindowX - myWindowY, 0)) && (c < myWindowSize)) {
                        diseq[r][c] = theLD.getP(z);
                    }
                }
                lowerLabel = "P value";
                break;
            }
            case DPRIME: {
                for (int z = 0, n = theLD.getRowCount(); z < n; z++) {
                    int r = theLD.getX(z) - myWindowX + (myWindowSize / 2);
                    int c = theLD.getY(z) - myWindowY + (myWindowSize / 2);
                    if ((r >= 0) && (r < myWindowSize) && (c >= Math.max(r + myWindowX - myWindowY, 0)) && (c < myWindowSize)) {
                        diseq[r][c] = theLD.getDPrime(z);
                    }
                }
                lowerLabel = "D'";
                break;
            }
            case RSQUARE: {
                for (int z = 0, n = theLD.getRowCount(); z < n; z++) {
                    int r = theLD.getX(z) - myWindowX + (myWindowSize / 2);
                    int c = theLD.getY(z) - myWindowY + (myWindowSize / 2);
                    if ((r >= 0) && (r < myWindowSize) && (c >= Math.max(r + myWindowX - myWindowY, 0)) && (c < myWindowSize)) {
                        diseq[r][c] = theLD.getRSqr(z);
                    }
                }
                lowerLabel = "R^2";
                break;
            }
        }

        lowerProb = (ldMeasure == P_VALUE) ? true : false;
    }

    public void setUpperCorner(int ldMeasure) {

        for (int c = 0; c < myWindowSize; c++) {
            for (int r = Math.max(c + myWindowY - myWindowX, 0); r < myWindowSize; r++) {
                diseq[r][c] = Double.NaN;
            }
        }

        switch (ldMeasure) {
            case P_VALUE: {
                for (int z = 0, n = theLD.getRowCount(); z < n; z++) {
                    int r = theLD.getX(z) - myWindowY + (myWindowSize / 2);
                    int c = theLD.getY(z) - myWindowX + (myWindowSize / 2);
                    if ((r >= 0) && (r < myWindowSize) && (c >= Math.max(r + myWindowY - myWindowX, 0)) && (c < myWindowSize)) {
                        diseq[c][r] = theLD.getP(z);
                    }
                }
                upperLabel = "P value";
                break;
            }
            case DPRIME: {
                for (int z = 0, n = theLD.getRowCount(); z < n; z++) {
                    int r = theLD.getX(z) - myWindowY + (myWindowSize / 2);
                    int c = theLD.getY(z) - myWindowX + (myWindowSize / 2);
                    if ((r >= 0) && (r < myWindowSize) && (c >= Math.max(r + myWindowY - myWindowX, 0)) && (c < myWindowSize)) {
                        diseq[c][r] = theLD.getDPrime(z);
                    }
                }
                upperLabel = "D'";
                break;
            }
            case RSQUARE: {
                for (int z = 0, n = theLD.getRowCount(); z < n; z++) {
                    int r = theLD.getX(z) - myWindowY + (myWindowSize / 2);
                    int c = theLD.getY(z) - myWindowX + (myWindowSize / 2);
                    if ((r >= 0) && (r < myWindowSize) && (c >= Math.max(r + myWindowY - myWindowX, 0)) && (c < myWindowSize)) {
                        diseq[c][r] = theLD.getRSqr(z);
                    }
                }
                upperLabel = "R^2";
                break;
            }
        }

        upperProb = (ldMeasure == P_VALUE) ? true : false;

    }

    /**
     * This sets the scale of the LD view, either sites are organized by chromosomes if
     * chromosomalScale is true, otherwise they are organized by genes
     */
    public void setScaleOfView(boolean chromosomalScale) {
        this.chromosomalScale = chromosomalScale;
        countGenesAndChromosomes();
        calculateStartAndEndPositions();
    }

    /**
     * This sets whether a schematic is displayed.  If true a schematic of genes or
     * chromosomes is displayed, otherwise no schematic is displayed
     */
    public void setShowSchematic(boolean includeBlockSchematic) {
        if (theAA == null) {
            return;  //if there is no annotation don't produce the schematic
        }
        this.includeBlockSchematic = includeBlockSchematic;
        countGenesAndChromosomes();
        calculateStartAndEndPositions();
    }

    public void setShowLabels(boolean includeLabels) {
        if (this.includeLabels == includeLabels) {
            return;
        }
        this.includeLabels = includeLabels;
        countGenesAndChromosomes();
        calculateStartAndEndPositions();
    }

    /**
     * this counts the number of separate blocks
     * if on chromosomal scale then chromosomes are counted otherwise only loci are counted
     * It then deteremines the total span in terms of cM or bases depending on scale
     */
    private void countGenesAndChromosomes() {
        totalLoci = totalChromosomes = 0;
        normalizer = theAA.getPositionInLocus(0) - 0.5;
        if (normalizer < 0) {
            normalizer = 0;
        }
        int currc = -1999;
        String currLocus = "";
        for (int r = 0; r < totalVariableSites; r++) {
            if (!currLocus.equals(theAA.getLocusName(r))) {
                totalLoci++;
                currLocus = theAA.getLocusName(r);
            }
        }

        //the number of separate totalBlocks
        totalBlocks = (chromosomalScale) ? totalChromosomes : totalLoci;
        if (totalBlocks == 0) {
            totalBlocks = 1;
        }
        System.out.println("totalBlocks: " + totalBlocks);
        blockStart = new double[totalBlocks];
        blockEnd = new double[totalBlocks];
        blockNames = new String[totalBlocks];

        for (int i = 0; i < totalChromosomes; i++) {
            blockStart[i] = 999999;
            blockEnd[i] = -999999;
        }
        int c = -1;
        currLocus = "unknown locus";
        currc = -1999;
        for (int r = 0; r < totalVariableSites; r++) {
            if (!currLocus.equals(theAA.getLocusName(r))) {
                c++;
                currLocus = theAA.getLocusName(r);
                blockNames[c] = currLocus;
            }
            if (blockStart[c] > theAA.getPositionInLocus(r)) {
                blockStart[c] = theAA.getPositionInLocus(r);
            }
            if (blockEnd[c] < theAA.getPositionInLocus(r)) {
                blockEnd[c] = theAA.getPositionInLocus(r);
            }
        }
        totalUnits = 0.5f;
        for (int i = 0; i < totalBlocks; i++) {
            if ((chromosomalScale) && ((blockEnd[i] - blockStart[i]) < minimumChromosomeLength)) {
                blockEnd[i] = blockStart[i] + minimumChromosomeLength;
            } else if ((blockEnd[i] - blockStart[i]) < 1) {
                blockEnd[i] = blockStart[i] + 1;
            }
            totalUnits += blockEnd[i] - blockStart[i];
        }//+1;}
    }

    /**
     * this determines to relative positions of the sites and cartoons (everything ranges from 0..1)
     *
     */
    private void calculateStartAndEndPositions() {
        //This will determine were all the relative positions of the sites go
        double proportionPerPolymorphism, proportionPerUnit = 0.0f;
        if (includeBlockSchematic) {
            totalIntervals = myWindowSize + totalBlocks - 1;
            proportionPerPolymorphism = 1 / (double) myWindowSize;
            proportionPerUnit = (proportionPerPolymorphism * myWindowSize) / totalUnits;
            blockBeginPos = new double[totalBlocks];    //These hold the start and end points of the genes
            blockEndPos = new double[totalBlocks];
        } else {
            totalIntervals = myWindowSize;
            proportionPerPolymorphism = 1 / (double) totalIntervals;
        }
        startPos = new double[myWindowSize];
        endPos = new double[myWindowSize];

        // int r,b=0,currg=-1999,currc=-1999;
        startPos[0] = 0;
        endPos[0] = 0;
        double currStartBase = 0, currEndBase = 0, geneToChromosomeSpace = 0;
        for (int r = myWindowX - (myWindowSize / 2); r < myWindowX + (myWindowSize / 2 + (myWindowSize % 2)); r++) {
            startPos[r - myWindowX + (myWindowSize / 2)] = currStartBase;
            currStartBase += proportionPerPolymorphism;
        }  //end of going through sites
        if (includeBlockSchematic) {
            currStartBase = 0;
            for (int b = 0; b < totalBlocks; b++) {
                blockBeginPos[b] = b / (double) totalBlocks;
                blockEndPos[b] = (b + 1) / (double) totalBlocks;
            }
            int currB = 0;
            for (int i = 1; i < myWindowX - (myWindowSize / 2) + 1; i++) {
                if (!theAA.getLocusName(i).equals(theAA.getLocusName(i - 1))) {
                    currB++;
                }
            }
            endPos[0] = blockBeginPos[currB] + (theAA.getPositionInLocus(myWindowX - (myWindowSize / 2)) - blockStart[currB]) / (blockEnd[currB] - blockStart[currB]) / totalBlocks;
            for (int r = myWindowX - (myWindowSize / 2) + 1; r < myWindowX + (myWindowSize / 2 + (myWindowSize % 2)); r++) {
                if (!theAA.getLocusName(r).equals(theAA.getLocusName(r - 1))) {
                    currB++;
                }
                endPos[r - myWindowX + (myWindowSize / 2)] = blockBeginPos[currB] + (theAA.getPositionInLocus(r) - blockStart[currB]) / (blockEnd[currB] - blockStart[currB]) / totalBlocks;
            }
        }
    }

    private void jbInit() throws Exception {
        this.setBackground(Color.red);
        this.setSize(400, 400);
        // this.setPreferredSize(new Dimension(400, 400));
        // this.setLayout(borderLayout1);
    }

    private Color getMagnitudeColor(int r, int c) {
        if (r + myWindowX == c + myWindowY) {
            return theColor.getHSBColor(0.999f, (float) diseq[r][c], 1f);
        }
        if (Double.isNaN(diseq[r][c])) {
            return theColor.lightGray;
        }
        if (diseq[r][c] > 0.999) {
            return theColor.getHSBColor(1f, 1f, 1f);
        }
        if (diseq[r][c] < -998.0) {
            return theColor.lightGray;
        }
        if ((((float) diseq[r][c]) + (1.0000f - ((float) diseq[r][c])) / 2) < 0.52f) {
            return theColor.getHSBColor(((float) diseq[r][c]) + (1.0000f - ((float) diseq[r][c])) / 2 - .5f, ((float) diseq[r][c]) + (1.0000f - ((float) diseq[r][c])) / 2 - .5f, 1f);
        } else {
            return theColor.getHSBColor(((float) diseq[r][c]) + (1.0000f - ((float) diseq[r][c])) / 2, ((float) diseq[r][c]) + (1.0000f - ((float) diseq[r][c])) / 2, 1f);
        }
    }

    private Color getProbabilityColor(int r, int c) {
        double p1 = 0.01, p2 = 0.001, p3 = 0.0001;
        if (Double.isNaN(diseq[r][c])) {
            return theColor.lightGray;
        }
        if (diseq[r][c] < -998.0) {
            return theColor.lightGray;
        }
        if (diseq[r][c] > p1) {
            return theColor.white;
        }
        if (diseq[r][c] > p2) {
            return theColor.blue;
        }
        if (diseq[r][c] > p3) {
            return theColor.green;
        }
        return theColor.red;
    }

    private void addPolymorphismLabels(Graphics g, int ih) {

        String s;
        g.setFont(new java.awt.Font("Dialog", 0, 9));
        g.setColor(theColor.black);
        String locus = theAA.getLocusName(myWindowY - (myWindowSize / 2));
        int jump = 0;
        for (int c = myWindowY - (myWindowSize / 2); c < myWindowY + (myWindowSize / 2 + (myWindowSize % 2)); c++) {
            if (locus.equals(theAA.getLocusName(c + jump))) {
                s = theAA.getLocusName(c + jump) + "s" + theAA.getPositionInLocus(c + jump);
            } else {
                s = "";
                locus = theAA.getLocusName(c + jump);
                jump--;
            }
            g.drawString(s, 4, yPos[c - myWindowY + (myWindowSize / 2)] + ih - 1);
        }
    }

    /**
     * This converts all those relative positions to real coordinates based on the size of the component
     */
    private void calculateCoordinates(Graphics gr) {
        Dimension d = this.getSize();
        double iwf, ihf, xSize, ySize;
        ySize = d.height - voff - distanceBetweenGraphAndGene;
        ihf = ySize / (double) myWindowSize;
        xSize = d.width - hoff - h2off;
        iwf = xSize / (double) myWindowSize;
        ih = (int) Math.round(ihf);
        iw = (int) Math.round(iwf);
        for (int r = 0; r < myWindowSize; r++) {
            xPos[r] = (int) ((startPos[r] * xSize) + (double) hoff);
            yPos[r] = (int) ((startPos[r] * ySize) + (double) voff);
        }  //end of going through sites
        xPos[myWindowSize] = (int) d.width - h2off;
        yPos[myWindowSize] = (int) ySize + voff;
        if (includeBlockSchematic) {
            for (int r = 0; r < myWindowSize; r++) {
                xEndPos[r] = (int) Math.round((endPos[r] * xSize) + hoff);
            }  //end of going through sites
            blockBeginX = new int[totalBlocks];
            blockEndX = new int[totalBlocks];
            for (int b = 0; b < totalBlocks; b++) {
                blockBeginX[b] = (int) Math.round((blockBeginPos[b] * xSize) + hoff);
                blockEndX[b] = (int) Math.round((blockEndPos[b] * xSize) + hoff);
            }
        }
    }

    protected void paintComponent(Graphics g) {
        if (diseq == null) {
            return;
        }
        Dimension d = this.getSize();
        calculateCoordinates(g);
        g.setColor(theColor.white);
        g.fillRect(0, 0, d.width, d.height);
        System.out.println("UpperProb=" + upperProb + "  LowerProb=" + lowerProb);
        g.setColor(theColor.darkGray);
        g.fillRect(xPos[0], yPos[0], xPos[myWindowSize] - xPos[0], yPos[myWindowSize] - yPos[0] + 2);
        String xLocus = theAA.getLocusName(myWindowX - (myWindowSize / 2));
        int xJump = 0;
        for (int r = 0; r < myWindowSize; r++) {
            if (!xLocus.equals(theAA.getLocusName(r + xJump + myWindowX - (myWindowSize / 2)))) {
                xLocus = theAA.getLocusName(r + xJump + myWindowX - (myWindowSize / 2));
                g.setColor(theColor.darkGray);
                for (int c = 0; c < myWindowSize; c++) {
                    g.fillRect(xPos[r], yPos[c], iw + 1, ih + 1);
                }
                xJump--;
            } else {
                String yLocus = theAA.getLocusName(myWindowY - (myWindowSize / 2));
                int yJump = 0;
                for (int c = 0; c < myWindowSize; c++) {
                    if (!yLocus.equals(theAA.getLocusName(c + yJump + myWindowY - (myWindowSize / 2)))) {
                        yLocus = theAA.getLocusName(c + yJump + myWindowY - (myWindowSize / 2));
                        g.setColor(theColor.darkGray);
                        yJump--;
                    } else if (((c + yJump + myWindowY < r + xJump + myWindowX) && (upperProb == true)) || ((c + yJump + myWindowY > r + xJump + myWindowX) && (lowerProb == true))) {
                        g.setColor(getProbabilityColor(r + xJump, c + yJump));
                    } else if (r + xJump + myWindowX == c + yJump + myWindowY) {
                        g.setColor(theColor.black);
                    } else {
                        g.setColor(getMagnitudeColor(r + xJump, c + yJump));
                    }
                    g.fillRect(xPos[r], yPos[c], iw + 1, ih + 1);
                }
            }
        }

        if (includeLabels) {
            addPolymorphismLabels(g, ih);
        }

        if (includeBlockSchematic) {
            addGenePicture(g, ih, iw);
        }
        addLegend(g);
    }

    public void paint(Graphics g) {
        paintComponent(g);
    }

    private void addLegend(Graphics g) {
        Dimension d = this.getSize();
        int localX = d.width - h2off + 10;
        int mid = d.height / 2;
        g.setColor(Color.black);
        g.drawString(upperLabel, localX, 10);
        addLegendGraph(g, upperProb, localX, 20, mid - 10);
        g.setColor(Color.black);
        g.drawString(lowerLabel, localX, mid + 10);
        addLegendGraph(g, lowerProb, localX, mid + 20, d.height - 10);
    }

    private void addLegendGraph(Graphics g, boolean prob, int xStart, int yStart, int yEnd) {
        DecimalFormat dF; //=new DecimalFormat("0.0000");
        int yInc, currY = yStart;
        int barWidth = 10;
        if (prob) {
            yInc = (yEnd - yStart) / 4;
            g.setColor(theColor.white);
            g.fillRect(xStart, currY, barWidth, yInc);
            g.setColor(Color.black);
            g.drawRect(xStart, currY, barWidth, yInc);
            g.drawString(">0.01", xStart + barWidth + 5, currY + 10);
            currY += yInc;
            g.setColor(theColor.blue);
            g.fillRect(xStart, currY, barWidth, yInc);
            g.setColor(Color.black);
            g.drawRect(xStart, currY, barWidth, yInc);
            g.drawString("<0.01", xStart + barWidth + 5, currY + 10);
            currY += yInc;
            g.setColor(theColor.green);
            g.fillRect(xStart, currY, barWidth, yInc);
            g.setColor(Color.black);
            g.drawRect(xStart, currY, barWidth, yInc);
            g.drawString("<0.001", xStart + barWidth + 5, currY + 10);
            currY += yInc;
            g.setColor(theColor.red);
            g.fillRect(xStart, currY, barWidth, yInc);
            g.setColor(Color.black);
            g.drawRect(xStart, currY, barWidth, yInc);
            g.drawString("<0.0001", xStart + barWidth + 5, currY + 10);
        } else {
            yInc = (yEnd - yStart) / 11;
            dF = new DecimalFormat("0.00");
            for (double d = 1.0000f; d >= 0.5; d -= 0.05) {
                g.setColor(theColor.getHSBColor((float) d, (float) d, 1f));
                g.fillRect(xStart, currY, barWidth, yInc);
                g.setColor(Color.black);
                g.drawRect(xStart, currY, barWidth, yInc);
                g.drawString(dF.format(d - (1.0001f - d)), xStart + barWidth + 5, currY + 10);
                currY += yInc;
            }
            g.setColor(theColor.getHSBColor((float) 0, (float) 0, 1f));
            g.fillRect(xStart, currY, barWidth, yInc);
            g.setColor(Color.black);
            g.drawRect(xStart, currY, barWidth, yInc);
            g.drawString(dF.format(0.00), xStart + barWidth + 5, currY + 10);
        }
    }

    private void addGenePicture(Graphics g, int ih, int iw) {
        //This will add the gene picture to the left of the polymorphisms
        int yOfLinkBlock, yOfGene, yOfGeneLabel;//,totalBases,spacer, cpos;
        int halfIW = iw / 2;
        // MultiAlleleSiteCharacteristic theMSC, lastMSC;
        Dimension d = this.getSize();
        yOfLinkBlock = yPos[myWindowSize];
        yOfGene = yOfLinkBlock + (distanceBetweenGraphAndGene / 2);
        yOfGeneLabel = yOfLinkBlock + (int) (0.8f * (double) distanceBetweenGraphAndGene);
        String locus = theAA.getLocusName(myWindowX - (myWindowSize / 2));
        int jump = 0;
        for (int r = 0; r < myWindowSize; r++) {
            if (locus.equals(theAA.getLocusName(r + jump + myWindowX - (myWindowSize / 2)))) {
                g.drawLine(xPos[r] + halfIW, yOfLinkBlock + 1, xEndPos[r], yOfGene);
            } else {
                locus = theAA.getLocusName(r + jump + myWindowX - (myWindowSize / 2));
                jump--;
            }
        }  //end of going through sites
        for (int b = 0; b < totalBlocks; b++) {
            g.setColor(iterColor(b));
            g.drawLine(blockBeginX[b], yOfGene, blockEndX[b], yOfGene);
            g.drawLine(blockBeginX[b], yOfGene + 1, blockEndX[b], yOfGene + 1);
            g.drawString(blockNames[b], blockBeginX[b], yOfGeneLabel);
        }
    }

    private Color iterColor(int iter) {
        Color newColor;
        if (iter % 5 == 0) {
            newColor = theColor.blue;
        } else if (iter % 5 == 1) {
            newColor = theColor.red;
        } else if (iter % 5 == 2) {
            newColor = theColor.green;
        } else if (iter % 5 == 3) {
            newColor = theColor.cyan;
        } else {
            newColor = theColor.orange;
        }
        return newColor;
    }
}