package net.maizegenetics.baseplugins.haplotypegraph;

import java.awt.*;
import java.util.Enumeration;

/**
 * Created by IntelliJ IDEA.
 * User: ajf25
 * Date: Jul 8, 2004
 * Time: 10:55:06 AM
 * To change this template use File | Settings | File Templates.
 */
// only stores the weightedshortestpaths of the drawn edges, when add new edges create new instance of this class
class KamadaKawai
{
    HaplotypeGraph graph;
    HaplotypeVertex[] vertexArray;

    final float EPSILON = 0.1f;

	int currentIteration;
    int maxIterations = 2000;

    // arbitrary const number for the spring force
	final double Kspring = 150000;
    // arbitrary const number for the repulsion force
    final double Krepulse= 1000;
	// distance matrix
    int[] dm;

    // used to pull the graph towards the center of the screen
	boolean adjustForGravity = true;
	boolean exchangeVertices = true;

    WeightedShortestPaths weightedShortestPaths;

    // The diameter of the visible graph. In other words, length of
    // the longest shortest path between any two vertices of the visible graph.
	int diameter;

    KamadaKawai(HaplotypeGraph graph)
    {
        this.graph = graph;
        currentIteration = 0;

        // populate vertexArray and edgeArray
        vertexArray = new HaplotypeVertex[ graph.vertices.size() ];
        dm = new int[ graph.vertices.size()*graph.vertices.size() ];
        int index = 0;
        Enumeration vertexKeys = graph.vertices.keys();
        while(vertexKeys.hasMoreElements())
        {
            String vertexKey = (String)vertexKeys.nextElement();
            vertexArray[index] = (HaplotypeVertex)graph.vertices.get(vertexKey);
            index++;
        }

        // populate the weightedShortestPaths
        weightedShortestPaths = new WeightedShortestPaths(this.graph);

        // calculate the diameter of the drawn graph - the longest shortest path in the graph
        diameter = 0;
		for (int i = 0; i < vertexArray.length - 1; i++)
			for (int j = i + 1; j < vertexArray.length; j++)
            {
				int dist = weightedShortestPaths.getDistance(vertexArray[i], vertexArray[j]);
				if (dist > diameter)
					diameter = dist;
			}
        // initialize the distance matrix dm this is done separately because the graph is not necessarily connected
        for (int i = 0; i < vertexArray.length - 1; i++)
			for (int j = i + 1; j < vertexArray.length; j++)
            {
				int dist = weightedShortestPaths.getDistance(vertexArray[i], vertexArray[j]);
				if (dist < 0)
                {
                    dm[i*vertexArray.length+j] = diameter+1;
				    dm[j*vertexArray.length+i] = diameter+1;
                }
                else
                {
                    dm[i*vertexArray.length+j] = dist;
				    dm[j*vertexArray.length+i] = dist;
                }
			}
    }

    // makes 1 step towards the minimum energy state of the system
    public void advancePositions( Dimension panelSize )
    {
		currentIteration++;
		double energy = calcEnergy();

		int n = vertexArray.length;
        if (n == 0)
            return;

		double maxDeltaM = 0;
        // the index of the node having max deltaM
		int pm = -1;
		for (int i = 0; i < n; i++)
        {
            if ( vertexArray[i].isSelected() )
                continue;
			double deltam = calcDeltaM(i);
			if (maxDeltaM < deltam)
            {
				maxDeltaM = deltam;
				pm = i;
			}
		}
		if (pm == -1)
            return;

        for (int i = 0; i < 100; i++)
        {
			int[] dxy = calcDeltaXY(pm);
            //System.out.println("in advance positions dxy[0]: " + dxy[0] + " dxy[1]: " + dxy[1]);
            //System.out.println("in advance vertexpos: " + vertexArray[pm].getX() + " " + vertexArray[pm].getY());    //****
            dxy[0] += vertexArray[pm].getX();
            dxy[1] += vertexArray[pm].getY();

            // do boundary check to make sure still on the panel
            if( dxy[0] < 0 )
                dxy[0] = 0;
            if( dxy[0] > panelSize.width )
                dxy[0] = panelSize.width;
            if( dxy[1] < 0 )
                dxy[1] = 0;
            if( dxy[1] > panelSize.height )
                dxy[1] = panelSize.height;

			vertexArray[pm].setLocation(dxy);
            vertexArray[pm].updateEdgePositions();
			double deltam = calcDeltaM(pm);
            if (deltam < EPSILON)
                break;
		}

		if (adjustForGravity)
			adjustForGravity(panelSize);

		if (exchangeVertices && maxDeltaM < EPSILON)
        {
            energy = calcEnergy();
			for (int i = 0; i < n - 1; i++)
            {
                if ( vertexArray[i].isSelected() )
                    continue;
				for (int j = i + 1; j < n; j++)
                {
                    if ( vertexArray[j].isSelected() )
                        continue;
					double xenergy = calcEnergyIfExchanged(i, j);
					if (energy > xenergy)
                    {
                        int sx = vertexArray[i].getX();
                        int sy = vertexArray[i].getY();
                        vertexArray[i].setX( vertexArray[j].getX() );
                        vertexArray[i].setY( vertexArray[j].getY() );
                        vertexArray[i].updateEdgePositions();

                        vertexArray[j].setX(sx);
                        vertexArray[j].setY(sy);
                        vertexArray[j].updateEdgePositions();
						return;
					}
				}
			}
		}
	}


    // Shift all vertices so that the center of gravity is located at the center of the screen.
	public void adjustForGravity(Dimension panelSize)
    {
		double height = panelSize.getHeight();
		double width = panelSize.getWidth();
		double gx = 0;
		double gy = 0;
		for (int i = 0; i < vertexArray.length; i++)
        {
			gx += vertexArray[i].getX();
			gy += vertexArray[i].getY();
		}
		gx /= vertexArray.length;
		gy /= vertexArray.length;
		double diffx = width / 2 - gx;
		double diffy = height / 2 - gy;
		for (int i = 0; i < vertexArray.length; i++)
        {
            if( vertexArray[i].isSelected() )
                continue;
            int[] newLoc = new int[2];
            newLoc[0] = (int)diffx + vertexArray[i].getX();
            newLoc[1] = (int)diffy + vertexArray[i].getY();
            // do boundary checking to make sure stays within the panel
            if( newLoc[0] < 0 )
                newLoc[0] = 0;
            if( newLoc[0] > panelSize.width )
                newLoc[0] = panelSize.width;
            if( newLoc[1] < 0 )
                newLoc[1] = 0;
            if( newLoc[1] > panelSize.height )
                newLoc[1] = panelSize.height;

			vertexArray[i].setLocation(newLoc);
            vertexArray[i].updateEdgePositions();
		}
	}

    // Determines a step to new position of the vertex m.
	private int[] calcDeltaXY(int m)
    {
		double dE_dxm = 0;
		double dE_dym = 0;
		double d2E_d2xm = 0;
		double d2E_dxmdym = 0;
		double d2E_dymdxm = 0;
		double d2E_d2ym = 0;

        // attempt at repulsion
        double deltaXrepulse = 0;
        double deltaYrepulse = 0;

		for (int i = 0; i < vertexArray.length; i++)
        {
			if (i != m)
            {
                // dealing with distance from vertex[m] and vertex[i]
				int dist = dm[m*vertexArray.length + i];

                HaplotypeEdge edge = (HaplotypeEdge)graph.allEdges.get( vertexArray[m].getSequence() + " " + vertexArray[i].getSequence() );
                if(edge == null)
                    edge = (HaplotypeEdge)graph.allEdges.get( vertexArray[i].getSequence() + " " + vertexArray[m].getSequence() );
                if(edge == null)
                    System.out.println("in calcDeltaXY edge from node " + m + " " + i + " was not found");

                // make sure the edge is drawn if not continue
                if(edge.isDrawn())
                {
                    int l_mi = edge.L * dist;
                    double k_mi = Kspring / (dist * dist);
                    int dx = vertexArray[m].getX() - vertexArray[i].getX();
                    int dy = vertexArray[m].getY() - vertexArray[i].getY();
                    double d = Math.sqrt(dx * dx + dy * dy);
                    double ddd = d * d * d;

                    dE_dxm += k_mi * (1 - l_mi / d) * dx;
                    dE_dym += k_mi * (1 - l_mi / d) * dy;
                    d2E_d2xm += k_mi * (1 - l_mi * dy * dy / ddd);
                    d2E_dxmdym += k_mi * l_mi * dx * dy / ddd;
                                                //d2E_dymdxm += k_mi * l_mi * dy * dx / ddd;
                    d2E_d2ym += k_mi * (1 - l_mi * dx * dx / ddd);

                }
                // attempt to add in repuslsionary forces between nodes
                int dx = vertexArray[m].getX() - vertexArray[i].getX();
                int dy = vertexArray[m].getY() - vertexArray[i].getY();
                double d = Math.sqrt(dx * dx + dy * dy);
                /*
                dE_dxm += (Krepulse/(3*d*dx))*100;
                dE_dym += (Krepulse/(3*d*dx))*100;
                d2E_d2xm += ((Krepulse*d*d*dx*dx)/(3*d*dx))*100;
                d2E_d2ym += ((Krepulse*d*d*dy*dy)/(3*d*dy))*100;
                d2E_dxmdym += ((Krepulse*d)/(3*dx*dy))*100;
                */
                /*
                dE_dxm += Krepulse*d/dx;
                dE_dym += Krepulse*d/dy;
                d2E_d2xm += Krepulse*dx/d;
                d2E_d2ym += Krepulse*dy/d;
                d2E_dxmdym += Krepulse*dy/(d*dx);
                */
                double repulseForce = Krepulse/(d*d);
                deltaXrepulse += dx*repulseForce;
                deltaYrepulse += dy*repulseForce;

			}
		}
		// d2E_dymdxm equals to d2E_dxmdym.
		d2E_dymdxm = d2E_dxmdym;

		double denomi = d2E_d2xm * d2E_d2ym - d2E_dxmdym * d2E_dymdxm;
		double deltaX = (d2E_dxmdym * dE_dym - d2E_d2ym * dE_dxm) / denomi;
		double deltaY = (d2E_dymdxm * dE_dxm - d2E_d2xm * dE_dym) / denomi;
        deltaX += deltaXrepulse*0.1;
        deltaY += deltaYrepulse*0.1;
        //System.out.println("dx: " + deltaX + " dy: " + deltaY);    //****
        //System.out.println("dx: " + (int)deltaX + " dy: " + (int)deltaY);    //****
		return new int[]{(int)deltaX, (int)deltaY};
	}

    // Calculates the gradient of energy function at the vertex m.
	private double calcDeltaM(int m)
    {
		double dEdxm = 0;
		double dEdym = 0;
		for (int i = 0; i < vertexArray.length; i++)
        {
			if (i != m)
            {
				double dist = dm[m*vertexArray.length + i];

                HaplotypeEdge edge = (HaplotypeEdge)graph.allEdges.get( vertexArray[m].getSequence() + " " + vertexArray[i].getSequence() );
                if(edge == null)
                    edge = (HaplotypeEdge)graph.allEdges.get( vertexArray[i].getSequence() + " " + vertexArray[m].getSequence() );
                if(edge == null)
                    System.out.println("in calcDeltaXY edge from node " + m + " " + i + " was not found");

                // make sure the edge is drawn if not continue
                if(edge.isDrawn())
                {
                    double l_mi = edge.L * dist;

                    double k_mi = Kspring / (dist * dist);

                    double dx = vertexArray[m].getX() - vertexArray[i].getX();
                    double dy = vertexArray[m].getY() - vertexArray[i].getY();
                    double d = Math.sqrt(dx * dx + dy * dy);

                    double common = k_mi * (1 - l_mi / d);
                    dEdxm += common * dx;
                    dEdym += common * dy;
                    //System.out.println("adding a spring force of: " + common*dx);    //****
                }
                // attempt to add in repuslsionary forces between nodes
                double dx = vertexArray[m].getX() - vertexArray[i].getX();
                double dy = vertexArray[m].getY() - vertexArray[i].getY();
                double d = Math.sqrt(dx * dx + dy * dy);

                //dEdxm += Krepulse/(3*d*dx);
                //dEdym += Krepulse/(3*d*dx);
                //dEdxm += Krepulse*d/dx;
                //dEdym += Krepulse*d/dy;
                //dEdxm += 1/dx*Krepulse/(d*d)*0.01;
                //dEdym += 1/dy*Krepulse/(d*d)*0.01;
                //System.out.println("adding a repulse force of: " + Math.abs(Krepulse/(3*d*dx)));    //****
			}
		}
		return Math.sqrt(dEdxm * dEdxm + dEdym * dEdym);
	}

    // Calculates the energy function E.
	private double calcEnergy()
    {
		double energy = 0;
		for (int i = 0; i < vertexArray.length - 1; i++)
        {
			for (int j = i + 1; j < vertexArray.length; j++)
            {
				double dist = dm[i*vertexArray.length + i];

                HaplotypeEdge edge = (HaplotypeEdge)graph.allEdges.get( vertexArray[i].getSequence() + " " + vertexArray[j].getSequence() );
                if(edge == null)
                    edge = (HaplotypeEdge)graph.allEdges.get( vertexArray[j].getSequence() + " " + vertexArray[i].getSequence() );
                if(edge == null)
                    System.out.println("in calcDeltaXY edge from node " + i + " " + j + " was not found");

                // make sure the edge is drawn if not continue
                if(edge.isDrawn())
                {
                    double l_ij = edge.L * dist;

                    double k_ij = Kspring / (dist * dist);
                    double dx = vertexArray[i].getX() - vertexArray[j].getX();
                    double dy = vertexArray[i].getY() - vertexArray[j].getY();
                    double d = Math.sqrt(dx * dx + dy * dy);

                    energy += k_ij / 2 * (dx * dx + dy * dy + l_ij * l_ij - 2 * l_ij * d);
                }
                // attempt to add in repuslsionary forces between nodes
                double dx = vertexArray[i].getX() - vertexArray[j].getX();
                double dy = vertexArray[i].getY() - vertexArray[j].getY();
                //energy += Math.abs((Krepulse*dx + Krepulse*dy)/(dist*dist*dist));
                energy += Krepulse/(dist*dist);
			}
		}
		return energy;
	}

    // Calculates the energy function E as if positions of the specified vertices are exchanged.
	private double calcEnergyIfExchanged(int p, int q)
    {
		if (p >= q)
			throw new RuntimeException("p should be < q");
		double energy = 0;		// < 0
		for (int i = 0; i < vertexArray.length - 1; i++)
			for (int j = i + 1; j < vertexArray.length; j++)
            {
				int ii = i;
				int jj = j;
				if (i == p) ii = q;
				if (j == q) jj = p;

				double dist = dm[j*vertexArray.length + i];
                HaplotypeEdge edge = (HaplotypeEdge)graph.allEdges.get( vertexArray[i].getSequence() + " " + vertexArray[j].getSequence() );
                if(edge == null)
                    edge = (HaplotypeEdge)graph.allEdges.get( vertexArray[j].getSequence() + " " + vertexArray[i].getSequence() );
                if(edge == null)
                    System.out.println("in calcDeltaXY edge from node " + i + " " + j + " was not found");

                // make sure the edge is drawn if not continue
                if(edge.isDrawn())
                {
                    double l_ij = edge.L * dist;

                    double k_ij = Kspring / (dist * dist);
                    double dx = vertexArray[ii].getX() - vertexArray[jj].getX();
                    double dy = vertexArray[ii].getY() - vertexArray[jj].getY();
                    double d = Math.sqrt(dx * dx + dy * dy);

                    energy += k_ij / 2 * (dx * dx + dy * dy + l_ij * l_ij - 2 * l_ij * d);
                }
                // attempt to add in repuslsionary forces between nodes
                double dx = vertexArray[ii].getX() - vertexArray[jj].getX();
                double dy = vertexArray[ii].getY() - vertexArray[jj].getY();
                //energy += Math.abs((Krepulse*dx + Krepulse*dy)/(dist*dist*dist));
                energy += Krepulse/(dist*dist);
			}
		return energy;
	}

}
