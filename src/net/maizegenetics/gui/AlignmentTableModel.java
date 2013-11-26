/*
 * AlignmentTableModel.java
 *
 */
package net.maizegenetics.gui;

import net.maizegenetics.dna.snp.Alignment;
import net.maizegenetics.taxa.TaxaList;

import javax.swing.*;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.table.AbstractTableModel;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author  terryc
 */
public class AlignmentTableModel extends AbstractTableModel implements ChangeListener {

    private static final NumberFormat NUMBER_FORMAT = NumberFormat.getPercentInstance();

    public enum COLUMN_NAME_TYPE {

        physicalPosition, siteNumber, locus, alleles, siteName
    };
    private COLUMN_NAME_TYPE myColumnNameType = COLUMN_NAME_TYPE.physicalPosition;
    private boolean myIsPhysicalPosition = true;
    private final Alignment myAlignment;
    // Left and Right variables
    private int myHorizontalPageSize = 0;
    private int myHorizontalCenter = 0;
    private int myHorizontalStart = 0;
    private int myHorizontalEnd = 0;

    public AlignmentTableModel(Alignment alignment, int horizontalPageSize) {

        if (alignment == null) {
            throw new IllegalArgumentException("AlignmentTableModel: init: alignment can not be null.");
        }

        myAlignment = alignment;

        myHorizontalCenter = myAlignment.numberOfSites() / 2;

        setHorizontalPageSize(horizontalPageSize);

    }

    public AlignmentTableModel(Alignment alignment) {
        this(alignment, 100);
    }

    // Return values appropriate for the visible table part
    public int getRowCount() {
        return myAlignment.numberOfTaxa();
    }

    public int getColumnCount() {
        return myHorizontalPageSize;
    }

    // Only works on the visible part of the table
    public Object getValueAt(int row, int col) {

        Object result = null;
        int realColumn = 0;

        try {
            realColumn = col + myHorizontalStart;
            result = myAlignment.genotypeAsString(row, realColumn);
        } catch (Exception e) {
            e.printStackTrace();
            System.out.println("row: " + row + "   col: " + col + "   realColumn: " + realColumn);
            throw new IllegalStateException("AlignmentTableModel: getValueAt: row: " + row + "   col: " + col + "   realColumn: " + realColumn);
        }

        return result;

    }

    public int getRealColumnIndex(int col) {
        return col + myHorizontalStart;
    }

    public Object getRealValueAt(int row, int col) {
        return myAlignment.genotype(row, col);
    }

    public String getColumnName(int col) {

        int realColumn = col + myHorizontalStart;

        if (myColumnNameType == COLUMN_NAME_TYPE.physicalPosition) {
            return String.valueOf(realColumn) + ": " + String.valueOf(myAlignment.chromosomalPosition(realColumn));
        } else if (myColumnNameType == COLUMN_NAME_TYPE.siteNumber) {
            return String.valueOf(realColumn) + ": " + String.valueOf(myAlignment.chromosomalPosition(realColumn));
        } else if (myColumnNameType == COLUMN_NAME_TYPE.locus) {
            return myAlignment.getChromosomeName(realColumn);
        } else if (myColumnNameType == COLUMN_NAME_TYPE.siteName) {
            return myAlignment.siteName(realColumn);
        } else if (myColumnNameType == COLUMN_NAME_TYPE.alleles) {
            int[][] alleles = myAlignment.allelesSortedByFrequency(realColumn);
            int numAlleles = alleles[0].length;
            double total = 0.0;

            for (int i = 0; i < numAlleles; i++) {
                total = total + alleles[1][i];
            }

            StringBuilder builder = new StringBuilder();
            for (int i = 0; i < numAlleles; i++) {
                if (i != 0) {
                    builder.append("; ");
                }
                builder.append(NUMBER_FORMAT.format((double) alleles[1][i] / total));
                builder.append(myAlignment.genotypeAsString(realColumn, (byte) alleles[0][i]));
            }

            return builder.toString();
        }

        return String.valueOf(myAlignment.chromosomalPosition(realColumn));

    }

    public void setColumnNameType(COLUMN_NAME_TYPE type) {
        myColumnNameType = type;
        if (myColumnNameType == COLUMN_NAME_TYPE.physicalPosition) {
            myIsPhysicalPosition = true;
        } else if (myColumnNameType == COLUMN_NAME_TYPE.siteNumber) {
            myIsPhysicalPosition = false;
        }
        fireTableStructureChanged();
    }

    public boolean isPhysicalPosition() {
        return myIsPhysicalPosition;
    }

    public COLUMN_NAME_TYPE getColumnNameType() {
        return myColumnNameType;
    }

    public int getRealColumnCount() {
        return myAlignment.numberOfSites();
    }

    public int getHorizontalPageSize() {
        return myHorizontalPageSize;
    }

    public void setHorizontalPageSize(int horizontalSize) {

        if (horizontalSize == myHorizontalPageSize) {
            return;
        }

        if (horizontalSize > getRealColumnCount()) {
            myHorizontalPageSize = getRealColumnCount();
        } else {
            myHorizontalPageSize = horizontalSize;
        }

        adjustPositionInternal(myHorizontalCenter);

    }

    public int getHorizontalCenter() {
        return myHorizontalCenter;
    }

    private void adjustPositionInternal(int position) {

        if (position < 0) {
            position = 0;
        } else if (position > getRealColumnCount() - 1) {
            position = getRealColumnCount() - 1;
        }

        myHorizontalCenter = position;

        int start = position - myHorizontalPageSize / 2;
        if (start < 0) {
            start = 0;
        }

        int end = start + myHorizontalPageSize - 1;
        if (end >= getRealColumnCount()) {
            end = getRealColumnCount() - 1;
            start = end - myHorizontalPageSize + 1;
        }

        myHorizontalStart = start;
        myHorizontalEnd = end;

        fireTableStructureChanged();

    }

    public void adjustPosition(int position) {

        if (isPhysicalPosition()) {
            position = myAlignment.siteOfPhysicalPosition(position, myAlignment.getChromosome(0));
            if (position < 0) {
                position = -position;
            }
        }

        adjustPositionInternal(position);

    }

    public void adjustPositionToSite(int site) {
        adjustPositionInternal(site);
    }

    public void adjustPositionToCenter() {
        adjustPositionInternal(myAlignment.numberOfSites() / 2);
    }

    /**
     * Resets the table backing this matrix table model to
     * an empty table.
     */
    public void resetTable() {
    }

    /**
     * Always returns false.
     */
    public boolean isCellEditable(int rowIndex, int columnIndex) {
        return false;
    }

    /**
     * No operation.
     */
    public void setValueAt(Object aValue, int rowIndex, int columnIndex) {
        // NO OPERATION
    }

    public List getRowHeaders() {

        List result = new ArrayList();

        TaxaList idGroup = myAlignment.taxa();
        for (int i = 0, n = idGroup.getTaxaCount(); i < n; i++) {
            result.add(idGroup.get(i));
        }

        return result;

    }

    public void fireTableChanged() {
        fireTableStructureChanged();
    }

    public void stateChanged(ChangeEvent e) {

        JSlider source = (JSlider) e.getSource();
        int position = (int) source.getValue();
        adjustPosition(position);

    }
}
