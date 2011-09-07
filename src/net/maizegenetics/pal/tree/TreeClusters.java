package net.maizegenetics.pal.tree;

import java.util.Arrays;

import net.maizegenetics.pal.ids.Identifier;

/**
 * @author Peter Bradbury 2/23/2010
 * For a Tree and a given number of clusters, this class finds the members of those clusters.
 */
public class TreeClusters {
	Tree theTree;
	int nTaxa;
	double[] sortedHeight;
	int nNodes;
	
	public TreeClusters(Tree theTree) {
		this.theTree = theTree;
		nTaxa = theTree.getIdCount();
		nNodes = theTree.getInternalNodeCount();
		sortedHeight = new double[nNodes];
		for (int n = 0; n < nNodes; n++) {
			sortedHeight[n] = theTree.getInternalNode(n).getNodeHeight();
		}
		Arrays.sort(sortedHeight);
	}
	
	public int[] getGroups (int numberOfGroups) {
		int[] groups = new int[nTaxa];
		if (numberOfGroups == 1) {
			for (int t = 0; t < nTaxa; t++) groups[t] = 0;
		} else if (numberOfGroups == nTaxa) {
			for (int t = 0; t <nTaxa; t++) groups[t] = t;
		} else {
			double maxHeight = sortedHeight[nNodes - numberOfGroups];
			for (int t = 0; t <nTaxa; t++) groups[t] = -1;
			int group = 0;
			for (int t = 0; t < nTaxa; t++) {
				if (groups[t] == -1) {
					Node aNode = theTree.getExternalNode(t);
					while (aNode.getParent().getNodeHeight() <= maxHeight) {
						aNode = aNode.getParent();
					}
					setNodeToGroup(aNode, group, groups);
					group++;
				}
			}
		}
		return groups;
	}
	
	public void setNodeToGroup(Node aNode, int group, int[] groups) {
		if (aNode.isLeaf()) {
			Identifier nodeId = aNode.getIdentifier();
			int index = theTree.whichIdNumber(nodeId);
			groups[index] = group;
		}
		int nChildren = aNode.getChildCount();
		for (int c = 0; c < nChildren; c++) {
			setNodeToGroup(aNode.getChild(c), group, groups);
		}
	}
	
	public Identifier getTaxon(int whichTaxon) { return theTree.getIdentifier(whichTaxon);}
	public Tree getTree() {return theTree;}
}
