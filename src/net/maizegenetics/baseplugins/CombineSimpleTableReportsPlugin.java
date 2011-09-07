/*
 * CombineSimpleTableReportsPlugin.java
 *
 * Created on January 31, 2007, 10:10 AM
 *
 */

package net.maizegenetics.baseplugins;


import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import javax.swing.ImageIcon;

import net.maizegenetics.pal.report.SimpleTableReport;

import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;

/**
 *
 * @author terryc
 */
public class CombineSimpleTableReportsPlugin extends AbstractPlugin{
    
    private final List <SimpleTableReport> myList = new ArrayList();
    
    
    /**
     * Creates a new instance of CombineSimpleTableReportsPlugin
     */
    public CombineSimpleTableReportsPlugin() {
        super(null, false);
    }
    
    public DataSet performFunction(DataSet input) {
        
        addReportsToList(input);
        
        if (myList.size() == 0) {
            return null;
        }
        
        SimpleTableReport combinedReport = (SimpleTableReport) SimpleTableReport.getInstance(myList);
        myList.clear();
        DataSet result = new DataSet(new Datum("SimpleTableReport", combinedReport, null), this);
        fireDataSetReturned(result);

        return result;
        
    }
    
    public String getToolTipText() {
        return "";
    }
    
    public ImageIcon getIcon() {
        return null;
    }
    
    public String getButtonName() {
        return "";
    }
    
    
    private void addReportsToList(DataSet input) {
        
        if (input == null) {
            return;
        }
        
        List data = input.getDataOfType(SimpleTableReport.class);
        if ((data != null) && (data.size() != 0)) {
            Iterator itr = data.iterator();
            while (itr.hasNext()) {
                Datum current = (Datum)itr.next();
                myList.add((SimpleTableReport) current.getData());
            }
        }
        
    }
    
}
