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
public class LinkageDisequilibriumComponent extends Component {

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
    //hoff is on the left side for site labels
    //h2off is on the right side for legends
    boolean probability = true, upperProb = false, lowerProb = true;
    // boolean genesOrChromo=true;  //true if display genes , false if display chromosomes

    double normalizer;

    public LinkageDisequilibriumComponent(LinkageDisequilibrium theLD, boolean includeBlockSchematic, boolean chromosomalScale) {
        this.theLD = theLD;
        theAA = theLD.getAnnotatedAlignment();
        this.includeBlockSchematic = includeBlockSchematic;
        this.chromosomalScale = chromosomalScale;
        this.diseq = new double[theLD.getSiteCount()][theLD.getSiteCount()];
        setUpperCorner(RSQUARE);
        setLowerCorner(P_VALUE);
        totalVariableSites = theLD.getSiteCount();
        if (theAA != null) {
            countGenesAndChromosomes();
            calculateStartAndEndPositions();
        } else {
            includeBlockSchematic = false;
        }
        xPos = new int[theLD.getSiteCount() + 1];
        yPos = new int[theLD.getSiteCount() + 1];
        xEndPos = new int[theLD.getSiteCount() + 1];
        try {
            jbInit();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    /**
     * This determines what is displayed in the lower left corner.
     * Options are: P_VALUE, DPRIME, and RSQUARE
     */
    public void setLowerCorner(int ldMeasure) {
        for (int r = 0; r < theLD.getSiteCount(); r++) {
            for (int c = r; c < theLD.getSiteCount(); c++) {
                switch (ldMeasure) {
                    case P_VALUE: {
                        diseq[r][c] = theLD.getP(r, c);
                        lowerLabel = "P value";
                        break;
                    }
                    case DPRIME: {
                        diseq[r][c] = theLD.getDPrime(r, c);
                        lowerLabel = "D'";
                        break;
                    }
                    case RSQUARE: {
                        diseq[r][c] = theLD.getRSqr(r, c);
                        lowerLabel = "R^2";
                        break;
                    }
                }
            }
        }
        lowerProb = (ldMeasure == P_VALUE) ? true : false;
    }

    /**
     * This determines what is displayed in the upper right corner.
     * Options are: P_VALUE, DPRIME, and RSQUARE
     */
    public void setUpperCorner(int ldMeasure) {
        for (int c = 0; c < theLD.getSiteCount(); c++) {
            for (int r = c; r < theLD.getSiteCount(); r++) {
                switch (ldMeasure) {
                    case P_VALUE: {
                        diseq[r][c] = theLD.getP(r, c);
                        upperLabel = "P value";
                        break;
                    }
                    case DPRIME: {
                        diseq[r][c] = theLD.getDPrime(r, c);
                        upperLabel = "D'";
                        break;
                    }
                    case RSQUARE: {
                        diseq[r][c] = theLD.getRSqr(r, c);
                        upperLabel = "R^2";
                        break;
                    }
                }
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
        for (int r = 0; r < totalVariableSites; r++) //sites and chromosomes need to be sorted
        {
            //if (theAA.getChromosome(r) != currc) {
            //    totalChromosomes++;
            //    currc = theAA.getChromosome(r);
            //}
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
            //if (chromosomalScale) {
            //    if (theAA.getChromosome(r) != currc) {
            //        c++;
            //        currc = theAA.getChromosome(r);
            //        blockNames[c] = "Chr." + currc;
            //    }
            //    if (blockStart[c] > theAA.getChromosomePosition(r)) {
            //        blockStart[c] = theAA.getChromosomePosition(r);
            //    }
            //    if (blockEnd[c] < theAA.getChromosomePosition(r)) {
            //        blockEnd[c] = theAA.getChromosomePosition(r);
            //    }
            //} else {
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
            //}
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
    void calculateStartAndEndPositions() {
        //This will determine were all the relative positions of the sites go
        double proportionPerPolymorphism, proportionPerUnit = 0.0f;
        if (includeBlockSchematic) {
            totalIntervals = totalVariableSites + totalBlocks - 1;
            proportionPerPolymorphism = 1 / (double) totalIntervals;
            proportionPerUnit = (proportionPerPolymorphism * totalVariableSites) / totalUnits;
            blockBeginPos = new double[totalBlocks];    //These hold the start and end points of the genes
            blockEndPos = new double[totalBlocks];
        } else {
            totalIntervals = totalVariableSites;
            proportionPerPolymorphism = 1 / (double) totalIntervals;
        }
        startPos = new double[totalVariableSites];
        endPos = new double[totalVariableSites];

        // int r,b=0,currg=-1999,currc=-1999;
        startPos[0] = 0;
        endPos[0] = 0;
        double currStartBase = 0, currEndBase = 0, geneToChromosomeSpace = 0;
        for (int r = 0; r < totalVariableSites; r++) {
            //if ((chromosomalScale) && (r > 0) && (includeBlockSchematic) && (theAA.getChromosome(r) != theAA.getChromosome(r - 1))) //transition between chromosomes if on chromosomal scale
            //{
            //    currStartBase += proportionPerPolymorphism;
            //}
            if ((!chromosomalScale) && (r > 0) && (includeBlockSchematic) && (!theAA.getLocusName(r).equals(theAA.getLocusName(r - 1)))) //transition between loci if not at chromosomal scale
            {
                currStartBase += proportionPerPolymorphism;
            }
            startPos[r] = currStartBase;
            currStartBase += proportionPerPolymorphism;
        }  //end of going through sites
        if (includeBlockSchematic) {
            currStartBase = 0;
            for (int b = 0; b < totalBlocks; b++) {
                blockBeginPos[b] = currStartBase;
                blockEndPos[b] = blockBeginPos[b] + ((blockEnd[b] - blockStart[b]) * proportionPerUnit);
                currStartBase = blockEndPos[b] + proportionPerPolymorphism;
            }
            int currB = 0;
            //if (chromosomalScale) {
            //    endPos[0] = blockBeginPos[0] + ((theAA.getChromosomePosition(0) - blockStart[0]) * proportionPerUnit);
            //} else {
                endPos[0] = blockBeginPos[0] + ((theAA.getPositionInLocus(0) - blockStart[0] - normalizer) * (1/((1/proportionPerUnit) - normalizer)));
            //}
            for (int r = 1; r < totalVariableSites; r++) {
                //if (chromosomalScale) {
                //    if (theAA.getChromosome(r) != theAA.getChromosome(r - 1)) {
                //        currB++;
                //    }
                //    endPos[r] = blockBeginPos[currB] + ((theAA.getChromosomePosition(r) - blockStart[currB]) * proportionPerUnit);
                //} else {
                    if (!theAA.getLocusName(r).equals(theAA.getLocusName(r - 1))) {
                        currB++;
                    }
                    //endPos[r] = ((theAA.getPositionInLocus(r) - blockStart[currB]) * proportionPerUnit);
                    endPos[r] = blockBeginPos[currB] + ((theAA.getPositionInLocus(r) - blockStart[currB] - normalizer) * (1/((1/proportionPerUnit) - normalizer)));
                //}
            }
        }
    // blockBeginPos[b]=currEndBase;
    // blockEndPos[b]=endPos[r-1];
    }

    private void jbInit() throws Exception {
        this.setBackground(Color.red);
        this.setSize(400, 400);
    // this.setPreferredSize(new Dimension(400, 400));
    // this.setLayout(borderLayout1);
    }

    private Color getMagnitudeColor(int r, int c) {
        if (r == c) {
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
        if ((((float) diseq[r][c]) + (1.0000f - ((float) diseq[r][c]))/2) < 0.52f) {
            return theColor.getHSBColor(((float) diseq[r][c]) + (1.0000f - ((float) diseq[r][c]))/2 - .5f, ((float) diseq[r][c]) + (1.0000f - ((float) diseq[r][c]))/2 - .5f, 1f);
        }
        else {
            return theColor.getHSBColor(((float) diseq[r][c]) + (1.0000f - ((float) diseq[r][c]))/2, ((float) diseq[r][c]) + (1.0000f - ((float) diseq[r][c]))/2, 1f);
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
        for (int r = 0; r < totalVariableSites; r++) {
            //if (chromosomalScale) {
            //    s = theAA.getChromosome(r) + "c" + Math.round(theAA.getChromosomePosition(r));
            //} else {
                s = theAA.getLocusName(r) + "s" + theAA.getPositionInLocus(r);
            //}
            g.drawString(s, 4, yPos[r] + ih - 1);
        }
    }

    /**
     * This converts all those relative positions to real coordinates based on the size of the component
     */
    private void calculateCoordinates(Graphics gr) {
        Dimension d = this.getSize();
        double iwf, ihf, xSize, ySize;
        ySize = d.height - voff - distanceBetweenGraphAndGene;
        ihf = ySize / (double) totalIntervals;
        xSize = d.width - hoff - h2off;
        iwf = xSize / (double) totalIntervals;
        ih = (int) Math.round(ihf);
        iw = (int) Math.round(iwf);
        for (int r = 0; r < totalVariableSites; r++) {
            xPos[r] = (int) ((startPos[r] * xSize) + (double) hoff);
            yPos[r] = (int) ((startPos[r] * ySize) + (double) voff);
        //xEndPos[r]=Math.round((endPos[r]*xSize)+hoff);
        }  //end of going through sites
        xPos[totalVariableSites] = (int) d.width - h2off;
        yPos[totalVariableSites] = (int) ySize + voff;
        if (includeBlockSchematic) {
            for (int r = 0; r < totalVariableSites; r++) {
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
        // super.paintComponent(g);
        int hue;
        Dimension d = this.getSize();
        calculateCoordinates(g);
        g.setColor(theColor.white);
        g.fillRect(0, 0, d.width, d.height);
        System.out.println("UpperProb=" + upperProb + "  LowerProb=" + lowerProb);
        g.setColor(theColor.darkGray);
        g.fillRect(xPos[0], yPos[0], xPos[totalVariableSites] - xPos[0], yPos[totalVariableSites] - yPos[0] + 2);
        for (int r = 0; r < totalVariableSites; r++) {
            for (int c = 0; c < totalVariableSites; c++) {
                if (((c < r) && (upperProb == true)) || ((c > r) && (lowerProb == true))) {
                    g.setColor(getProbabilityColor(r, c));
                } else if (r == c) {
                    g.setColor(theColor.black);
                } else {
                    g.setColor(getMagnitudeColor(r, c));
                }
                g.fillRect(xPos[r], yPos[c], iw + 1, ih + 1);
            }
        }

        // Removed grid lines because the cover too much
        // on large graphs.  -terryc
        /*
        g.setColor(theColor.darkGray);
        for(int r=0; r<totalVariableSites; r++) {
        g.drawLine(xPos[r], yPos[0], xPos[r], yPos[totalVariableSites]);
        g.drawLine(xPos[0], yPos[r], xPos[totalVariableSites], yPos[r]);
        }
         */

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
        g.drawString("Upper " + upperLabel, localX, 10);
        addLegendGraph(g, upperProb, localX, 20, mid - 10);
        g.setColor(Color.black);
        g.drawString("Lower " + lowerLabel, localX, mid + 10);
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
//                if (d >= 0.55) {
                    g.setColor(theColor.getHSBColor((float) d, (float) d, 1f));
//                }
//                else {
//                    g.setColor(theColor.getHSBColor((float) 0, (float) 0, 1f));
//                }
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
        yOfLinkBlock = yPos[totalVariableSites];
        yOfGene = yOfLinkBlock + (distanceBetweenGraphAndGene / 2);
        yOfGeneLabel = yOfLinkBlock + (int) (0.8f * (double) distanceBetweenGraphAndGene);

        for (int r = 0; r < totalVariableSites; r++) {
            g.drawLine(xPos[r] + halfIW, yOfLinkBlock + 1, xEndPos[r], yOfGene);
        }  //end of going through sites
        g.setColor(theColor.blue);
        for (int b = 0; b < totalBlocks; b++) {
            g.drawLine(blockBeginX[b], yOfGene, blockEndX[b], yOfGene);
            g.drawLine(blockBeginX[b], yOfGene + 1, blockEndX[b], yOfGene + 1);
            g.drawString(blockNames[b], blockBeginX[b], yOfGeneLabel);
        }
    }
}