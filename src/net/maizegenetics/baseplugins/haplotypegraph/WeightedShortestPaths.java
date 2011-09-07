package net.maizegenetics.baseplugins.haplotypegraph;

import org._3pq.jgrapht.alg.DijkstraShortestPath;
import org._3pq.jgrapht.edge.UndirectedWeightedEdge;
import org._3pq.jgrapht.graph.SimpleWeightedGraph;

import java.util.Enumeration;
import java.util.Hashtable;
import java.util.List;
import java.util.ListIterator;

/**
 * Created by IntelliJ IDEA.
 * User: ajf25
 * Date: Jul 8, 2004
 * Time: 11:15:56 AM
 * To change this template use File | Settings | File Templates.
 */
// only stores the weightedshortestpaths of the drawn edges, when add new edges create new instance of this class
public class WeightedShortestPaths
{
    Hashtable distances;

    public WeightedShortestPaths(HaplotypeGraph graph)
    {
        distances = new Hashtable();
        SimpleWeightedGraph swg = new SimpleWeightedGraph();

        // run through the graph adding all vertices to swg
        Enumeration vertexKeys = graph.vertices.keys();
        while(vertexKeys.hasMoreElements())
        {
            String vertexKey = (String)vertexKeys.nextElement();
            swg.addVertex( (HaplotypeVertex)graph.vertices.get(vertexKey) );
        }

        // run through the graph adding all the edges that are drawn to swg
        Enumeration edgeKeys = graph.allEdges.keys();
        while(edgeKeys.hasMoreElements())
        {
            String edgeKey = (String)edgeKeys.nextElement();
            HaplotypeEdge edge = (HaplotypeEdge)graph.allEdges.get(edgeKey);
            if( edge.isDrawn() )
            {
                UndirectedWeightedEdge uwe = new UndirectedWeightedEdge( graph.verticesBySequence.get(edge.seqStart), graph.verticesBySequence.get(edge.seqEnd), edge.cost );
                swg.addEdge(uwe);
            }
        }

        // run through all combinations of vertex pairs finding the shortest path distance
        for(int i = 0; i < graph.vertexCount; i++)
            for(int j = i+1; j < graph.vertexCount; j++)
            {
                HaplotypeVertex v1 = (HaplotypeVertex)graph.vertices.get(i+"");
                HaplotypeVertex v2 = (HaplotypeVertex)graph.vertices.get(j+"");
                List list = DijkstraShortestPath.findPathBetween(swg, v1, v2);
                // no path between the vertices
                if( list == null )
                {
                    distances.put( v1.getSequence() + " " + v2.getSequence(), new Integer(-1) );
                    distances.put( v2.getSequence() + " " + v1.getSequence(), new Integer(-1) );
                }
                // there is a path between them
                else
                {
                    int dist = 0;
                    ListIterator iter = list.listIterator();
                    while(iter.hasNext())
                    {
                        UndirectedWeightedEdge uwe = (UndirectedWeightedEdge)iter.next();
                        dist += uwe.getWeight();
                    }
                    // since undirected shortest path will be symmetric
                    distances.put( v1.getSequence() + " " + v2.getSequence(), new Integer(dist) );
                    distances.put( v2.getSequence() + " " + v1.getSequence(), new Integer(dist) );
                }
            }
    }

    public int getDistance(HaplotypeVertex v1, HaplotypeVertex v2)
    {
        return ((Integer)distances.get( v1.getSequence() + " " + v2.getSequence() )).intValue();
    }

}