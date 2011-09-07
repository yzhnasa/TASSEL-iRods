/*
 * TableReportTableModel.java
 *
 */
package net.maizegenetics.gui;

import java.awt.BorderLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.SwingConstants;
import javax.swing.table.AbstractTableModel;

import net.maizegenetics.pal.report.TableReport;

/**
 *
 * @author  terryc
 */
public class TableReportTableModel extends AbstractTableModel {

    private TableReport myTable = null;
    private Object[] myColumnHeadings = null;
    // Up and Down variables
    private int myPageSizeMax = 0;
    private int myPageSize = 0;
    private int myPageOffset = 0;
    private int myLastPageSize = 0;
    // Left and Right variables
    private int myHorizontalPageSizeMax = 0;
    private int myHorizontalPageSize = 0;
    private int myHorizontalPageOffset = 0;
    private int myHorizontalLastPageSize = 0;

    public TableReportTableModel(TableReport table, int pageVerticalSize, int pageHorizontalSize) {
        if (table == null) {
            throw new IllegalArgumentException("TableReportTableModel: init: table can not be null.");
        }

        myPageSizeMax = pageVerticalSize;
        myHorizontalPageSizeMax = pageHorizontalSize;

        myTable = table;
        myColumnHeadings = myTable.getTableColumnNames();
        myLastPageSize = myTable.getRowCount() % myPageSizeMax;
        myHorizontalLastPageSize = myTable.getColumnCount() % myHorizontalPageSizeMax;

        updateVerticalPageSize();
        updateHorizontalPageSize();

    }

    public TableReportTableModel(TableReport table) {
        this(table, 100, 100);
    }

    // Return values appropriate for the visible table part
    public int getRowCount() {
        return myPageSize;
    }

    public int getColumnCount() {
        return myHorizontalPageSize;
    }

    // Only works on the visible part of the table
    public Object getValueAt(int row, int col) {

        Object result = null;
        int realRow = 0;
        int realColumn = 0;

        try {
            realRow = row + (myPageOffset * myPageSizeMax);
            realColumn = col + (myHorizontalPageOffset * myHorizontalPageSizeMax);
            result = myTable.getValueAt(realRow, realColumn);
        } catch (RuntimeException e) {
            e.printStackTrace();
            System.out.println("row: " + row + "   col: " + col + "   realRow: " + realRow);
            System.out.println("pageOffset: " + myPageOffset + "   pageSize: " + myPageSize);
            System.out.println("Table column count: " + myTable.getColumnCount() + "   Table row count: " + myTable.getRowCount());
            throw e;
        }

        return result;

    }

    public Object getRealValueAt(int row, int col) {
        Object[] temp = myTable.getRow(row);
        return temp[col];
    }

    public String getColumnName(int col) {
        int realColumn = col + (myHorizontalPageOffset * myHorizontalPageSizeMax);
        return myColumnHeadings[realColumn].toString();
    }

    public int getPageOffset() {
        return myPageOffset;
    }

    public int getHorizontalPageOffset() {
        return myHorizontalPageOffset;
    }

    public int getPageCount() {
        return (int) Math.ceil((double) getRealRowCount() / myPageSizeMax);
    }

    public int getHorizontalPageCount() {
        return (int) Math.ceil((double) getRealColumnCount() / myHorizontalPageSizeMax);
    }

    public int getRealRowCount() {
        return myTable.getRowCount();
    }

    public int getRealColumnCount() {
        return myTable.getColumnCount();
    }

    public int getVerticalPageSize() {
        return myPageSize;
    }

    public int getHorizontalPageSize() {
        return myHorizontalPageSize;
    }

    public void setVerticalPageSize(int verticalSize) {

        if (verticalSize == myPageSizeMax) {
            return;
        }

        int oldPageSize = myPageSize;
        myPageSizeMax = verticalSize;
        myLastPageSize = myTable.getRowCount() % myPageSizeMax;
        updateVerticalPageSize();
        if (myPageSize < oldPageSize) {
            fireTableRowsDeleted(myPageSize, oldPageSize - 1);
        } else {
            fireTableRowsInserted(oldPageSize, myPageSize - 1);
        }

    }

    public void setHorizontalPageSize(int horizontalSize) {

        if (horizontalSize == myHorizontalPageSizeMax) {
            return;
        }

        int oldPageSize = myHorizontalPageSize;
        myHorizontalPageSizeMax = horizontalSize;
        myHorizontalLastPageSize = myTable.getColumnCount() % myHorizontalPageSizeMax;
        updateHorizontalPageSize();
        if (myHorizontalPageSize < oldPageSize) {
            fireTableStructureChanged();
        } else {
            fireTableStructureChanged();
        }

    }

    public void pageDown() {

        if (myPageOffset < getPageCount() - 1) {
            myPageOffset++;
            updateVerticalPageSize();
            fireTableDataChanged();
        }

    }

    public void pageUp() {

        if (myPageOffset > 0) {
            myPageOffset--;
            updateVerticalPageSize();
            fireTableDataChanged();
        }

    }

    public void pageLeft() {

        if (myHorizontalPageOffset > 0) {
            myHorizontalPageOffset--;
            updateHorizontalPageSize();
            fireTableStructureChanged();
        }

    }

    public void pageRight() {

        if (myHorizontalPageOffset < getHorizontalPageCount() - 1) {
            myHorizontalPageOffset++;
            updateHorizontalPageSize();
            fireTableStructureChanged();
        }

    }

    private void updateVerticalPageSize() {

        if (myPageOffset == (getPageCount() - 1)) {
            myPageSize = myLastPageSize;
        } else {
            myPageSize = myPageSizeMax;
        }

        System.out.println("updateVerticalPageSize: page: " + (myPageOffset + 1) + " of " + getPageCount() + "   page size: " + myPageSize);

    }

    private void updateHorizontalPageSize() {

        if (myHorizontalPageOffset == (getHorizontalPageCount() - 1)) {
            myHorizontalPageSize = myHorizontalLastPageSize;
        } else {
            myHorizontalPageSize = myHorizontalPageSizeMax;
        }

        System.out.println("updateHorizontalPageSize: page: " + (myHorizontalPageOffset + 1) + " of " + getHorizontalPageCount() + "   page size: " + myHorizontalPageSize);

    }

    public static JPanel createPagingControls(final TableReportTableModel model) {

        JPanel result = new JPanel(new BorderLayout());

        final JButton upButton = new JButton(new ArrowIcon(ArrowIcon.UP));
        upButton.setEnabled(false);
        final JButton downButton = new JButton(new ArrowIcon(ArrowIcon.DOWN));

        final JButton leftButton = new JButton(new ArrowIcon(ArrowIcon.LEFT));
        leftButton.setEnabled(false);
        final JButton rightButton = new JButton(new ArrowIcon(ArrowIcon.RIGHT));

        upButton.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent ae) {
                model.pageUp();

                if (model.getPageOffset() == 0) {
                    upButton.setEnabled(false);
                } else {
                    upButton.setEnabled(true);
                }

                if (model.getPageOffset() == (model.getPageCount() - 1)) {
                    downButton.setEnabled(false);
                } else {
                    downButton.setEnabled(true);
                }

            }
        });

        downButton.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent ae) {
                model.pageDown();

                if (model.getPageOffset() == 0) {
                    upButton.setEnabled(false);
                } else {
                    upButton.setEnabled(true);
                }

                if (model.getPageOffset() == (model.getPageCount() - 1)) {
                    downButton.setEnabled(false);
                } else {
                    downButton.setEnabled(true);
                }

            }
        });

        leftButton.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent ae) {
                model.pageLeft();

                if (model.getHorizontalPageOffset() == 0) {
                    leftButton.setEnabled(false);
                } else {
                    leftButton.setEnabled(true);
                }

                if (model.getHorizontalPageOffset() == (model.getHorizontalPageCount() - 1)) {
                    rightButton.setEnabled(false);
                } else {
                    rightButton.setEnabled(true);
                }

            }
        });

        rightButton.addActionListener(new ActionListener() {

            public void actionPerformed(ActionEvent ae) {
                model.pageRight();

                if (model.getHorizontalPageOffset() == 0) {
                    leftButton.setEnabled(false);
                } else {
                    leftButton.setEnabled(true);
                }

                if (model.getHorizontalPageOffset() == (model.getHorizontalPageCount() - 1)) {
                    rightButton.setEnabled(false);
                } else {
                    rightButton.setEnabled(true);
                }

            }
        });

        result.add(upButton, BorderLayout.NORTH);
        result.add(downButton, BorderLayout.SOUTH);
        result.add(new JLabel("Paging Controls", SwingConstants.CENTER));
        result.add(leftButton, BorderLayout.WEST);
        result.add(rightButton, BorderLayout.EAST);

        return result;

    }

    /**
     * Resets the table backing this matrix table model to
     * an empty table.
     */
    public void resetTable() {
    }

    public Object getColumnObject(int columnIndex) {
        return myColumnHeadings[columnIndex];
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

    public void fireTableChanged() {
        fireTableStructureChanged();
    }

    public Class<?> getColumnClass(int columnIndex) {
        return getValueAt(0, columnIndex).getClass();
    }
}
