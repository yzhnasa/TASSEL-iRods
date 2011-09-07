package net.maizegenetics.baseplugins.haplotypegraph;

/**
 * Created by IntelliJ IDEA.
 * User: ajf25
 * Date: Jul 6, 2004
 * Time: 1:04:32 PM
 * To change this template use File | Settings | File Templates.
 */
class HaplotypeEdge
{
    String key;
    String seqStart;
    String seqEnd;
    int cost;
    boolean isMinimum;
    boolean isDrawn;
    // x,y locations for drawing
    int startX;
    int startY;
    int endX;
    int endY;
    int thickness;
    // ideal length of this edge
    int L;

    HaplotypeEdge(String sequenceStart, String sequenceEnd, int cost)
    {
        this.key = sequenceStart + " " + sequenceEnd;
        this.seqStart = sequenceStart;
        this.seqEnd = sequenceEnd;
        this.cost = cost;
        isMinimum = false;
        isDrawn = false;
        startX = 0;
        startY = 0;
        endX = 0;
        endY = 0;
        thickness = 3;
        L = cost*50;
    }

    String getKey()
    {
        return this.key;
    }

    boolean isMinimum()
    {
        return this.isMinimum;
    }

    void setMinimum(boolean set)
    {
        this.isMinimum = set;
    }

    boolean isDrawn()
    {
        return this.isDrawn;
    }

    void setDrawn(boolean set)
    {
        this.isDrawn = set;
    }

    void moveLocationOf(String sequence, int x, int y)
    {
        if( seqStart.equalsIgnoreCase(sequence) )
        {
            this.startX = x;
            this.startY = y;
        }
        else if( seqEnd.equalsIgnoreCase(sequence) )
        {
            this.endX = x;
            this.endY = y;
        }
    }

    int getCost()
    {
        return this.cost;
    }

    boolean isStart(String seq)
    {
        return this.seqStart.equalsIgnoreCase(seq);
    }

    boolean isEnd(String seq)
    {
        return this.seqEnd.equalsIgnoreCase(seq);
    }

    String getStart()
    {
        return this.seqStart;
    }

    String getEnd()
    {
        return this.seqEnd;
    }

    // [0] is startX and [1] is startY
    int[] getStartLocation()
    {
        int[] startLoc = new int[2];
        startLoc[0] = startX;
        startLoc[1] = startY;
        return startLoc;
    }

    // [0] is endX and [1] is endY
    int[] getEndLocation()
    {
        int[] endLoc = new int[2];
        endLoc[0] = endX;
        endLoc[1] = endY;
        return endLoc;
    }

}