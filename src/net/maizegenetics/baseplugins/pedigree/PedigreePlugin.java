/*
 * PedigreePlugin.java
 *
 * Created on December 22, 2006, 5:02 PM
 *
 */

package net.maizegenetics.baseplugins.pedigree;


import gov.usda.gdpc.browser.Browser;
import gov.usda.gdpc.pedigree.PedigreeViewer;

import javax.swing.ImageIcon;

import javax.swing.JMenu;
import javax.swing.JPanel;

import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;





/**
 *
 * @author terryc
 */
public class PedigreePlugin extends AbstractPlugin{
    
    private final PedigreeViewer myViewer;
    
    
    /** Creates a new instance of PedigreePlugin */
    public PedigreePlugin() {
        
        // This really should be the same one used by Tassel
        Browser browser = Browser.createInstanceForTASSEL();
        
        myViewer = new PedigreeViewer(browser);
        
    }
    
    
    public DataSet performFunction (DataSet input) {
        // Need to work on this.
        return null;
    }
    
    
    public JPanel getPanel() {
        return myViewer;
    }
    
    
    public JMenu getMenu() {
        return myViewer.getPedigreeMenu();
    }

    public String getToolTipText() {
        return "Pedigree Viewer";
    }

    public ImageIcon getIcon() {
        return null;
    }

    public String getButtonName() {
        return "Pedigree Viewer";
    }
    
}
