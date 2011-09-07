/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.baseplugins.chart;

import java.util.Arrays;
import java.util.Hashtable;
import net.maizegenetics.pal.report.TableReport;
import org.jfree.data.xy.DefaultTableXYDataset;

/**
 *
 * @author yz79
 */
public class TableReportManhattanDataset extends DefaultTableXYDataset {
    double[][] theData;
    String[] seriesNames;
    String xName;
    int numberYAxes;

    String[] myChromNames;
    Object[] myColumnNames;
    double[] myPValues;
    double[] myLogPValues;
    int[] myPositions;
    String[] myMarkers;
    Hashtable myLookupTable;

    int myPValueColumnIndex = -1;
    int myChromColumnIndex = -1;
    int myPositionColumnIndex = -1;
    int myMarkerColumnIndex = -1;
    int myNumRows;

    // 1 = skip first row
    // 0 = don't skip
    int skipFirstRow;

    public TableReportManhattanDataset(TableReport theTable) {
        numberYAxes = 0;
        setTableReport(theTable);
    }

    public int getItemCount(int parm1) {
        return theData.length;
        //throw new java.lang.UnsupportedOperationException("Method getItemCount() not yet implemented.");
    }


    public Number getX(int series, int item) {
        Double x = new Double(theData[item][0]);
        return x;
        //    throw new java.lang.UnsupportedOperationException("Method getXValue() not yet implemented.");
    }

    public int getSeriesCount() {
        return numberYAxes;
        //throw new java.lang.UnsupportedOperationException("Method getSeriesCount() not yet implemented.");
    }


    public Number getY(int series, int item) {
        Double y = new Double(theData[item][1+series]);
        return y;
        //    throw new java.lang.UnsupportedOperationException("Method getYValue() not yet implemented.");
    }
    public String getSeriesName(int series) {
        /**current*/
        return seriesNames[series];
        //    throw new java.lang.UnsupportedOperationException("Method getSeriesName() not yet implemented.");
    }

    public String getSeriesKey(int series) {
        /**current*/
        return seriesNames[series];
        //    throw new java.lang.UnsupportedOperationException("Method getSeriesName() not yet implemented.");
    }

    public String getXName() {
        /**current*/
        return xName;
        //    throw new java.lang.UnsupportedOperationException("Method getSeriesName() not yet implemented.");
    }

    private void setPValueColumnIndex() {
        for (int i = 0; i < myColumnNames.length; i++) {
            if (myColumnNames[i].equals("p") || myColumnNames[i].equals("marker_p")) {
                myPValueColumnIndex = i;
                return;
            }
        }
        throw new IllegalArgumentException("No P-values in selected data");
    }

    private void setChromColumnIndex() {
        for (int i = 0; i < myColumnNames.length; i++) {
            if (myColumnNames[i].equals("Locus")) {
                myChromColumnIndex = i;
                return;
            }
        }
        throw new IllegalArgumentException("No Chromosome names in selected data");
    }

    private void setMarkerColumnIndex() {
        for (int i = 0; i < myColumnNames.length; i++) {
            if (myColumnNames[i].equals("Marker")) {
                myMarkerColumnIndex = i;
                return;
            }
        }
    }

    // MLM analysis produces extra blank first row, method to compensate
    private void setNumRows(TableReport myTableReport) {
        if (((String)myTableReport.getValueAt(0, myChromColumnIndex)).equals("")) {
            skipFirstRow = 1;
            myNumRows = myTableReport.getRowCount() - 1;
        }
        else {
            skipFirstRow = 0;
            myNumRows = myTableReport.getRowCount();
        }
    }

    private void setPositionColumnIndex() {
        for (int i = 0; i < myColumnNames.length; i++) {
            if (myColumnNames[i].equals("Locus_pos") || myColumnNames[i].equals("Site")) {
                myPositionColumnIndex = i;
                return;
            }
        }
        throw new IllegalArgumentException("No positions in selected data");
    }

    private void setChromNames(TableReport myTableReport) {
        String currentChrom = "";
        for (int i = 0; i < myChromNames.length; i++) {
            myChromNames[i]  = ((String)myTableReport.getValueAt(i + skipFirstRow, myChromColumnIndex));
            if (!currentChrom.equals(myChromNames[i])) {
                numberYAxes++;
                currentChrom = myChromNames[i];
            }
        }
    }

    private void setMarkers(TableReport myTableReport) {
        for (int i = 0; i < myMarkers.length; i++) {
            myMarkers[i] = ((String)myTableReport.getValueAt(i + skipFirstRow, myMarkerColumnIndex));
        }
    }

    private void setPValues(TableReport myTableReport) {
        for (int i = 0; i < myPValues.length; i++) {
            myPValues[i] = ((Double)myTableReport.getValueAt(i + skipFirstRow, myPValueColumnIndex)).doubleValue();
            myLookupTable.put(myPValues[i], i);
        }
    }

    public String[] getMarkers() {
        return myMarkers;
    }

    public String[] getChroms() {
        return myChromNames;
    }

    private void setPositions(TableReport myTableReport) {
        int offset = 0;
        int previousPosition = 0;
        int currentPosition = 0;
        // GLM positions formated as int
        if (myColumnNames[myPositionColumnIndex].equals("Locus_pos")) {
            myPositions[0] = ((Integer)myTableReport.getValueAt(skipFirstRow, myPositionColumnIndex)).intValue();
            previousPosition = myPositions[0];
            for (int i = 1; i < myPValues.length; i++) {
                currentPosition = ((Integer)myTableReport.getValueAt(i + skipFirstRow, myPositionColumnIndex)).intValue();
                myPositions[i] = currentPosition + offset;
                if (currentPosition < previousPosition) {
                    offset = offset + previousPosition;
                }
                previousPosition = currentPosition;
            }
        }
        // MLM positions formatted as string
        else if (myColumnNames[myPositionColumnIndex].equals("Site")) {
            myPositions[0] = Integer.valueOf((String)myTableReport.getValueAt(skipFirstRow, myPositionColumnIndex));
            previousPosition = myPositions[0];
            for (int i = 1; i < myPValues.length; i++) {
                currentPosition = Integer.valueOf((String)myTableReport.getValueAt(i + skipFirstRow, myPositionColumnIndex));
                myPositions[i] = currentPosition + offset;
                if (currentPosition < previousPosition) {
                    offset = offset + previousPosition;
                }
                previousPosition = currentPosition;
            }
        }
    }

    private void setLogPValues() {
        for (int i = 0; i < myLogPValues.length; i++) {
            myLogPValues[i] = -Math.log10(myPValues[i]);
        }
    }

    public void setTableReport(TableReport theTable) {
        myColumnNames = theTable.getTableColumnNames();
        setPValueColumnIndex();
        setChromColumnIndex();
        setPositionColumnIndex();
        setMarkerColumnIndex();
        setNumRows(theTable);
        myPValues = new double[myNumRows];
        myLogPValues = new double[myNumRows];
        myChromNames = new String[myNumRows];
        myPositions = new int[myNumRows];
        myMarkers = new String[myNumRows];
        myLookupTable = new Hashtable(myNumRows);
        setPValues(theTable);
        setLogPValues();
        setChromNames(theTable);
        setPositions(theTable);
        setMarkers(theTable);
        seriesNames = new String[numberYAxes];
        seriesNames[0] = myChromNames[0];

        String currentChrom = myChromNames[0];
        int chromIndex = 1;
        theData = new double[myNumRows][1 + numberYAxes];
        for (int i = 0; i < theData.length; i++) {
            for (int j = 0; j < theData[0].length; j++) {
                theData[i][j] = Double.NaN;
            }
        }
        for (int i = 0; i < myNumRows; i++) {
            try {
                theData[i][0] = myPositions[i];
                if (!currentChrom.equals(myChromNames[i])) {
                    chromIndex++;
                    currentChrom = myChromNames[i];
                    seriesNames[chromIndex - 1] = currentChrom;
                }
                theData[i][chromIndex] = myLogPValues[i];
            }
            catch (NumberFormatException ex) {
                System.out.println("throw new NumberFormatException();");
            }
        }
        xName = "Position";
    }
}