package net.maizegenetics.baseplugins.haplotypegraph;

import edu.uci.ics.jung.graph.Edge;
import edu.uci.ics.jung.graph.Vertex;
import edu.uci.ics.jung.graph.decorators.EdgeWeightLabeller;
import edu.uci.ics.jung.graph.impl.UndirectedSparseGraph;
import edu.uci.ics.jung.utils.GraphUtils;
import edu.uci.ics.jung.utils.UserData;
import net.maizegenetics.pal.alignment.Alignment;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.*;

/**
 * Created using IntelliJ IDEA.
 * Author: Peter Bradbury
 * Date: Aug 13, 2004
 * Time: 9:51:44 AM
 */
public class JungHaplotypeGraph extends UndirectedSparseGraph {

    Hashtable vertices;
    Hashtable verticesBySequence;
    int vertexCount;
    Hashtable allEdges;
    private ArrayList seqNames;
    private ArrayList seqs;
    private int lastDrawnType;
    private int lastDrawnThreshold;
    private HaplotypeVertex selectedVertex;
    private HashMap jungVertexToSeq;

    JungHaplotypeGraph(File file) {
        super();

        vertices = new Hashtable();
        verticesBySequence = new Hashtable();
        allEdges = new Hashtable();
        jungVertexToSeq = new HashMap();
        seqNames = new ArrayList();
        seqs = new ArrayList();
        lastDrawnType = -1;
        lastDrawnThreshold = -1;
        selectedVertex = null;
        try {
            BufferedReader br = new BufferedReader(new FileReader(file));
            String line = "";
            StringTokenizer st;
            String token = "";
            while ((line = br.readLine()) != null) {
                st = new StringTokenizer(line);
                token = st.nextToken();
                seqNames.add(token);
                token = st.nextToken();
                seqs.add(token);
            }
            buildGraph();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    JungHaplotypeGraph(Alignment aa) {
        super();
        vertices = new Hashtable();
        verticesBySequence = new Hashtable();
        allEdges = new Hashtable();
        jungVertexToSeq = new HashMap();
        seqNames = new ArrayList();
        seqs = new ArrayList();
        lastDrawnType = -1;
        lastDrawnThreshold = -1;
        selectedVertex = null;
        int aaCount = aa.getIdGroup().getIdCount();
        for (int i = 0; i < aaCount; i++) {
            seqNames.add(aa.getIdGroup().getIdentifier(i).toString());
            seqs.add(aa.getAlignedSequenceString(i).trim());
        }
        buildGraph();
    }

    void selectVertexAt(int x, int y) {
        // enumerate the vertices and find the first one located at x, y and set as the selectedVertex for dragging
        Enumeration vertexKeys = vertices.keys();
        while (vertexKeys.hasMoreElements()) {
            String vertexKey = (String) vertexKeys.nextElement();
            HaplotypeVertex vertex = (HaplotypeVertex) vertices.get(vertexKey);
            int verX = vertex.getX();
            int verY = vertex.getY();
            int verSize = vertex.getSize();

            int distance = (int) Math.sqrt((verX - x) * (verX - x) + (verY - y) * (verY - y));
            if (distance <= verSize / 2) {
                selectedVertex = vertex;
                selectedVertex.setSelected(true);
                return;
            }
        }
    }

    // get(0) is the sequence get(1) is the names
    ArrayList getVertexInfoAt(int x, int y) {
        // enumerate the vertices and find the first one located at x, y and return the names and sequence
        Enumeration vertexKeys = vertices.keys();
        while (vertexKeys.hasMoreElements()) {
            String vertexKey = (String) vertexKeys.nextElement();
            HaplotypeVertex vertex = (HaplotypeVertex) vertices.get(vertexKey);
            int verX = vertex.getX();
            int verY = vertex.getY();
            int verSize = vertex.getSize();

            int distance = (int) Math.sqrt((verX - x) * (verX - x) + (verY - y) * (verY - y));
            // if distance is within the radius of the vertex
            if (distance <= verSize / 2) {
                ArrayList info = new ArrayList();
                info.add(vertex.getSequence());
                info.add(vertex.getNames());
                return info;
            }
        }
        return null;
    }

    // get(0) is the start sequence, get(1) is the end sequence, get(2) is the cost/difference
    ArrayList getEdgeInfoAt(int x, int y) {
        // enumerate the edges and find the first one located at x, y and return the start sequence, end sequence, cost
        Enumeration edgeKeys = allEdges.keys();
        while (edgeKeys.hasMoreElements()) {
            String edgeKey = (String) edgeKeys.nextElement();
            HaplotypeEdge edge = (HaplotypeEdge) allEdges.get(edgeKey);

            // only check edges that are drawn
            if (edge.isDrawn) {
                int x1 = edge.startX;
                int y1 = edge.startY;
                int x2 = edge.endX;
                int y2 = edge.endY;

                // distance point to a line equations from: http://astronomy.swin.edu.au/~pbourke/geometry/pointline/
                double uNumerator = (x - x1) * (x2 - x1) + (y - y1) * (y2 - y1);
                double uDenominator = (x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1);
                double u = uNumerator / uDenominator;

                // make sure it is on the line segment we are looking at and not on a line based on some other two points
                if (u >= 0.0 && u <= 1.0) {
                    int xTan = (int) (x1 + u * (x2 - x1));
                    int yTan = (int) (y1 + u * (y2 - y1));

                    int distance = (int) Math.sqrt((xTan - x) * (xTan - x) + (yTan - y) * (yTan - y));
                    if (distance <= 5) {
                        ArrayList info = new ArrayList();
                        info.add(edge.seqStart);
                        info.add(edge.seqEnd);
                        info.add(edge.cost + "");
                        return info;
                    }
                }
            } // end isDrawn test on edge
        } // end of while loop
        return null;
    }

    void updateSelectedVertex(int x, int y) {
        if (selectedVertex != null) {
            int[] newLoc = new int[2];
            newLoc[0] = x;
            newLoc[1] = y;
            selectedVertex.setLocation(newLoc);

            // update all edges on the vertex's new position
            Hashtable edges = selectedVertex.getEdges();
            Enumeration edgeKeys = edges.keys();
            while (edgeKeys.hasMoreElements()) {
                String edgeKey = (String) edgeKeys.nextElement();
                HaplotypeEdge edge = (HaplotypeEdge) edges.get(edgeKey);
                edge.moveLocationOf(selectedVertex.getSequence(), x, y);
            }
        }
    }

    void clearSelectedVertex() {
        if (selectedVertex == null) {
            return;
        }
        selectedVertex.setSelected(false);
        selectedVertex = null;
    }

    private void buildGraph() {
        vertexCount = 0;
        // this will contain the indices into the names and seqs lists of haplotypes that have '?'s in them
        ArrayList unknownHaplotypes = new ArrayList();
        // run through the seqs ArrayList grouping same seqs together
        for (int i = 0; i < seqs.size(); i++) {
            String seq = (String) seqs.get(i);

            // if this seq does not contain ? process and compare all other seqs to it
            // todo: decide how to handle haplotypes that have a '?' in them place in group they match aside from the '?', if more than 1 group fits that criteria place in all the groups and mark it differently
            if (seq.indexOf("?") < 0) {
                ArrayList sameSeqNames = new ArrayList();
                for (int j = i; j < seqs.size(); j++) {
                    String newSeq = (String) seqs.get(j);

                    // if the sequences are the same add the new sequences name to the arraylist of sameseqnames
                    // and remove the sequence and seqname from the class variable arraylists, do not want to process
                    // same sequence more than once
                    if (seq.equalsIgnoreCase(newSeq)) {
                        String newName = (String) seqNames.get(j);
                        sameSeqNames.add(newName);

                        // if i and j are not the same remove j
                        // this prevents the accidental jumping over of a sequence
                        // (removing i everything gets bumped back so the next sequence we want to look at is at i,
                        // but the loop increments i and we would end up jumping that sequence
                        if (i != j) {
                            seqNames.remove(j);
                            seqs.remove(j);
                            // avoid jumping over one sequence
                            j--;
                        }
                    }
                }
                HaplotypeVertex vertex = new HaplotypeVertex(sameSeqNames, seq);
                vertices.put(vertexCount + "", vertex);
                verticesBySequence.put(vertex.getSequence(), vertex);
                vertexCount++;
            } else {
                unknownHaplotypes.add(new Integer(i));
            }
        } // end of building all the vertices of the graph of none '?' containing sequences

        // run through the unknownHaplotypes placing the one's with '?'s into the group(s) they probably belong to
        for (int i = 0; i < unknownHaplotypes.size(); i++) {
            // vertexMatchKeys will store the keys of the vertices this unknownSequence match
            // if the size of this arraylist is > 1, need to indicate that this sequence is in more than vertex
            ArrayList vertexMatchKeys = new ArrayList();
            String unknownSequence = (String) seqs.get(((Integer) unknownHaplotypes.get(i)).intValue());
            String unknownName = (String) seqNames.get(((Integer) unknownHaplotypes.get(i)).intValue());
            Enumeration vertexKeys = vertices.keys();
            while (vertexKeys.hasMoreElements()) {
                String vertexKey = (String) vertexKeys.nextElement();
                HaplotypeVertex vertex = (HaplotypeVertex) vertices.get(vertexKey);
                String vertexSequence = vertex.getSequence();

                // do character comparasion between vertexSequence and unknownSequence
                // if the only difference is the '?' locations place unknownSequence into this vertex
                boolean addToVertex = true;
                if (unknownSequence.length() == vertexSequence.length()) {
                    for (int index = 0; index < unknownSequence.length(); index++) {
                        if (unknownSequence.charAt(index) != vertexSequence.charAt(index) && unknownSequence.charAt(index) != '?') {
                            addToVertex = false;
                        }
                    }
                    if (addToVertex) {
                        vertexMatchKeys.add(vertexKey);
                        ArrayList vertexNames = vertex.getNames();
                        vertexNames.add(unknownName + " " + unknownSequence);
                    }
                }
            }
            // if this unknown sequence was added to more than 1 vertex, need to indicate such
            if (vertexMatchKeys.size() > 1) {
                for (int j = 0; j < vertexMatchKeys.size(); j++) {
                    String key = (String) vertexMatchKeys.get(j);
                    HaplotypeVertex vertex = (HaplotypeVertex) vertices.get(key);
                    ArrayList vertexNames = vertex.getNames();
                    int location = vertexNames.indexOf(unknownName + " " + unknownSequence);
                    vertexNames.remove(location);
                    vertexNames.add(unknownName + " " + unknownSequence + " part of other groups");
                }
            }
            // if the unknownSequence did not match any known sequences comment on it for now will fix soon
            if (vertexMatchKeys.size() <= 0) {
                //System.out.println("an unknownSequence did not match any known sequences");    //****
                // create a HaplotypeVertex for this sequence
                ArrayList name = new ArrayList();
                name.add(unknownName);
                HaplotypeVertex vertex = new HaplotypeVertex(name, unknownSequence);
                vertices.put(vertexCount + "", vertex);
                verticesBySequence.put(vertex.getSequence(), vertex);
                vertexCount++;
            }
        }

        // grab each vertex and compare to everyother vertex updating adjacencyMatrix with the correct edge costs
        for (int startVertex = 0; startVertex < vertexCount; startVertex++) {
            HaplotypeVertex startVertexObject = (HaplotypeVertex) vertices.get(startVertex + "");
            String startVertexSeq = (String) startVertexObject.getSequence();

            // run through all other vertices comparing and calculating costs
            for (int endVertex = startVertex + 1; endVertex < vertexCount; endVertex++) {
                HaplotypeVertex endVertexObject = (HaplotypeVertex) vertices.get(endVertex + "");
                String endVertexSeq = (String) endVertexObject.getSequence();

                // all sequences will be the same length
                int diffCount = 0;
                for (int charLoc = 0; charLoc < startVertexSeq.length(); charLoc++) {
                    // only count as a difference if not a '?'
                    if (startVertexSeq.charAt(charLoc) != endVertexSeq.charAt(charLoc) && startVertexSeq.charAt(charLoc) != '?' && endVertexSeq.charAt(charLoc) != '?') {
                        diffCount++;
                    }
                }
                // undirected edges
                HaplotypeEdge edge = new HaplotypeEdge(startVertexSeq, endVertexSeq, diffCount);
                startVertexObject.addEdge(edge);
                endVertexObject.addEdge(edge);
                allEdges.put(edge.getKey(), edge);
            }
        } // end creating edges and adding to corresponding vertices
        createMinimumSpanningTree();
        lastDrawnType = HaplotypePanel.MINIMUM_SPANNING_TREE_SAME;
        updateTree();
        loadVerticesAndEdges();
    }

    // used for debugging purposes, prints out each vertex's key and sequence
    private void printVertices() {
        Enumeration vertexKeys = vertices.keys();
        while (vertexKeys.hasMoreElements()) {
            String vertexKey = (String) vertexKeys.nextElement();
            HaplotypeVertex vertex = (HaplotypeVertex) vertices.get(vertexKey);
            System.out.println("key is: " + vertexKey + " : seq is: " + vertex.getSequence());    //****
        }
    }

    private void createMinimumSpanningTree() {
        Hashtable group = new Hashtable();
        if (vertexCount < 0) {
            return;
        }
        HaplotypeVertex vertex = (HaplotypeVertex) vertices.get("0");
        group.put(vertex.getSequence(), vertex);

        // while not every vertex has been reached
        while (group.size() != vertices.size()) {
            // run through group and find the shortest edge cost leaving group
            Enumeration vertexKeys = group.keys();
            HaplotypeVertex minVertex = null;
            HaplotypeEdge minEdge = null;
            while (vertexKeys.hasMoreElements()) {
                String vertexKey = (String) vertexKeys.nextElement();
                vertex = (HaplotypeVertex) group.get(vertexKey);
                Hashtable edges = vertex.getEdges();
                Enumeration edgeKeys = edges.keys();
                // run through this vertex's edges finding the min that leaves group
                while (edgeKeys.hasMoreElements()) {
                    String edgeKey = (String) edgeKeys.nextElement();
                    HaplotypeEdge tempEdge = (HaplotypeEdge) edges.get(edgeKey);
                    // if this vertex is the start of the tempEdge and the end of the tempEdge is not in group
                    if (tempEdge.isStart(vertex.getSequence()) && !group.containsKey(tempEdge.getEnd())) {
                        // compare tempEdges cost to minEdge's cost, if less we have a new min edge
                        if (minEdge == null) {
                            minEdge = tempEdge;
                            minVertex = (HaplotypeVertex) verticesBySequence.get(tempEdge.getEnd());
                        } else if (tempEdge.getCost() < minEdge.getCost()) {
                            minEdge = tempEdge;
                            minVertex = (HaplotypeVertex) verticesBySequence.get(tempEdge.getEnd());
                        }
                    } // if this vertex is the end of the tempEdge and the start of the tempEdge is not in group
                    else if (tempEdge.isEnd(vertex.getSequence()) && !group.containsKey(tempEdge.getStart())) {
                        // compare tempEdges cost to minEdge's cost, if less we have a new min edge
                        if (minEdge == null) {
                            minEdge = tempEdge;
                            minVertex = (HaplotypeVertex) verticesBySequence.get(tempEdge.getStart());
                            ;
                        } else if (tempEdge.getCost() < minEdge.getCost()) {
                            minEdge = tempEdge;
                            minVertex = (HaplotypeVertex) verticesBySequence.get(tempEdge.getStart());
                            ;
                        }
                    }
                }
            }
            // if minEdge or minVertex is null we have a problem
            if (minEdge == null || minVertex == null) {
                System.out.println("VERY BAD MIN EDGE OR MIN VERTEX IS NULL SHOULD NOT HAPPEN");    //****
            }
            minEdge.setMinimum(true);
            minEdge.setDrawn(true);
            minVertex.setMinimumCost(minEdge.getCost());
            group.put(minVertex.getSequence(), minVertex);
        }
    }

    private void updateTree() {
        // enumerate the vertices hashtable changing the isDrawn in the edges based on the drawType
        Enumeration vertexKeys = vertices.keys();
        while (vertexKeys.hasMoreElements()) {
            String vertexKey = (String) vertexKeys.nextElement();
            HaplotypeVertex vertex = (HaplotypeVertex) vertices.get(vertexKey);

            // grab the edges of this vertex
            Hashtable edges = vertex.getEdges();
            // enumerate the edges and update based on drawType
            Enumeration edgeKeys = edges.keys();
            while (edgeKeys.hasMoreElements()) {
                String edgeKey = (String) edgeKeys.nextElement();
                HaplotypeEdge edge = (HaplotypeEdge) edges.get(edgeKey);
                if (lastDrawnType == HaplotypePanel.MINIMUM_SPANNING_TREE) {
                    // only draw edges who are in the minimum spanning tree
                    if (!edge.isMinimum()) {
                        edge.setDrawn(false);
                    } else {
                        edge.setDrawn(true);
                    }
                } else if (lastDrawnType == HaplotypePanel.MINIMUM_SPANNING_TREE_SAME) {
                    // only draw edges that are minimum or have same cost as minimum
                    if (!edge.isMinimum && edge.getCost() == vertex.getMinimumCost()) {
                        edge.setDrawn(true);
                    } else if (!edge.isMinimum) {
                        edge.setDrawn(false);
                    } else if (edge.isMinimum) {
                        edge.setDrawn(true);
                    }
                } else if (lastDrawnType == HaplotypePanel.MINIMUM_SPANNING_TREE_THRESHOLD) {
                    // only draw edges that have cost <= lastDrawnThreshold
                    if (edge.getCost() <= (lastDrawnThreshold)) {
                        edge.setDrawn(true);
                    } else {
                        edge.setDrawn(false);
                    }
                }
            }
        }
    }

    protected void initialize() {
        super.initialize();

    }

    /**
     * loads the graph vertice and edges from the Haplotype into Jung
     */
    public void loadVerticesAndEdges() {
        //create a hash map of jung vertices keyed on sequence and add the jung vertices to this graph
        this.removeAllEdges();
        this.removeAllVertices();
        HashMap jungVertices = new HashMap();
        jungVertexToSeq.clear();

        Enumeration enumVertex = vertices.elements();
        while (enumVertex.hasMoreElements()) {
            HaplotypeVertex haplotypeVertex = (HaplotypeVertex) enumVertex.nextElement();
            //add jung vertices to graph using addVertex and to jungVertices with sequence as key
            JungHaplotypeVertex v = new JungHaplotypeVertex();
            jungVertices.put(haplotypeVertex.getSequence(), v);
            v.setUserDatum("freq", new Integer(haplotypeVertex.getNames().size()), UserData.SHARED);
            v.setUserDatum("blank", "  ", UserData.SHARED);
            v.setToStringKey("freq");
            jungVertexToSeq.put(v, haplotypeVertex.getSequence());
            this.addVertex(v);
        }

        //add edges to the graph
        EdgeWeightLabeller ewl = EdgeWeightLabeller.getLabeller(this);
        Enumeration enumEdges = allEdges.elements();
        while (enumEdges.hasMoreElements()) {
            HaplotypeEdge haplotypeEdge = (HaplotypeEdge) enumEdges.nextElement();
            if (haplotypeEdge.isDrawn()) {
                Vertex v1 = (Vertex) jungVertices.get(haplotypeEdge.seqStart);
                Vertex v2 = (Vertex) jungVertices.get(haplotypeEdge.seqEnd);
                Edge e = GraphUtils.addEdge(this, v1, v2);
                ewl.setWeight(e, haplotypeEdge.getCost());
            }
        }
    }

    public HaplotypeVertex getHaplotypeVertexFromJungVertex(Vertex v) {
        String seq = (String) jungVertexToSeq.get(v);
        return (HaplotypeVertex) verticesBySequence.get(seq);
    }
}
