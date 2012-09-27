package net.maizegenetics.baseplugins;

import java.net.URL;
import java.util.List;

import javax.swing.ImageIcon;
import javax.swing.JOptionPane;

import net.maizegenetics.pal.distance.DistanceMatrix;
import net.maizegenetics.pal.tree.Tree;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;

public class ArchaeopteryxPlugin extends AbstractPlugin {

	@Override
	public DataSet performFunction(DataSet input) {
		List<Datum> treeInList = input.getDataOfType(Tree.class);
        if (treeInList.size() != 1) {
            String message = "Invalid selection.  Please select one tree.";
            if (isInteractive()) {
                JOptionPane.showMessageDialog(getParentFrame(), message);
            } else {
                System.out.println(message);
            }
            return null;
        }

        Tree myTree = (Tree) treeInList.get(0).getData();
        
        
		return null;
	}

	@Override
	public ImageIcon getIcon() {
        URL imageURL = TreeDisplayPlugin.class.getResource("images/Tree.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
	}

	@Override
	public String getButtonName() {
		return "Tree2";
	}

	@Override
	public String getToolTipText() {
		return "Tree View using Archaeopteryx viewer";
	}

}
