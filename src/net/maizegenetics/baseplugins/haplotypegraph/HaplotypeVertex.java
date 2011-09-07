package net.maizegenetics.baseplugins.haplotypegraph;

import java.awt.*;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.Hashtable;

/**
 * Created by IntelliJ IDEA.
 * User: ajf25
 * Date: Jul 6, 2004
 * Time: 12:58:38 PM
 * To change this template use File | Settings | File Templates.
 */
class HaplotypeVertex
{
    int x;
    int y;
    int size;
    // the cost of the edge that connected this vertex to the min spanning tree
    int minCost;
    Hashtable edges;
    ArrayList names;
    String sequence;
    // indicates if the vertex is being dragged by the user
    boolean isSelected = false;

    //todo - take into account the sequence names when determing the color of the vertex
    Color color;

    HaplotypeVertex(ArrayList names, String sequence)
    {
        this.names = names;
        this.sequence = sequence;
        edges = new Hashtable();
        x = 0;
        y = 0;
        size = 15;
        minCost = -1;
        color = Color.blue;
    }

    void addEdge(HaplotypeEdge edge)
    {
        edges.put(edge.getKey(), edge);
    }

    void moveLocation(int x, int y)
    {
        this.x = x;
        this. y = y;
    }

    String getSequence()
    {
        return this.sequence;
    }

    ArrayList getNames()
    {
        return this.names;
    }

    Hashtable getEdges()
    {
        return this.edges;
    }

    int getX()
    {
        return this.x;
    }

    void setX(int x)
    {
        this.x = x;
    }

    int getY()
    {
        return this.y;
    }

    void setY(int y)
    {
        this.y = y;
    }

    int getSize()
    {
        return this.size;
    }

    // x is in [0] and y is in [1]
    int[] getLocation()
    {
        int[] location = new int[2];
        location[0] = x;
        location[1] = y;
        return location;
    }

    // x is in [0] and y is in [1]
    void setLocation(int[] location)
    {
        this.x = location[0];
        this.y = location[1];
    }

    void setMinimumCost(int cost)
    {
        this.minCost = cost;
    }

    int getMinimumCost()
    {
        return this.minCost;
    }

    boolean isSelected()
    {
        return this.isSelected;
    }

    void setSelected(boolean set)
    {
        this.isSelected = set;
    }

    Color getColor()
    {
        return this.color;
    }

    // update all edge positions based on the current position of this vertex
    void updateEdgePositions()
    {
        Enumeration edgeKeys = edges.keys();
        while(edgeKeys.hasMoreElements())
        {
            String edgeKey = (String)edgeKeys.nextElement();
            HaplotypeEdge edge = (HaplotypeEdge)edges.get(edgeKey);
            edge.moveLocationOf(sequence, x, y);
        }
    }

}