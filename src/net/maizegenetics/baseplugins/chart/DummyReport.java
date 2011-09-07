package net.maizegenetics.baseplugins.chart;

import java.util.Random;

import net.maizegenetics.pal.report.AbstractTableReport;
import net.maizegenetics.pal.report.TableReport;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2004</p>
 * <p>Company: </p>
 * @author not attributable
 * @version 1.0
 */
public class DummyReport extends AbstractTableReport implements TableReport {

    static Random random = new Random();
    private Object[][] myData;
    private String[] myColumnHeadings;

    public DummyReport() {
        init();
    }

    //Implementation of TableReport Interface
    public Object[] getTableColumnNames() {
        return myColumnHeadings;
    }

    public Object[][] getTableData() {
        return myData;
    }

    /**
     * Returns specified row.
     *
     * @param row row number
     *
     * @return row
     */
    public Object[] getRow(int row) {
        return myData[row];
    }

    private void init() {

        myColumnHeadings = new String[]{"Site_Type", "StartSite", "XGauss", "YGAUSS"};

        Object[][] data;
        java.text.NumberFormat nf = new java.text.DecimalFormat();
        nf.setMaximumFractionDigits(8);
        int basicCols = 4, labelOffset;
        data = new String[500][basicCols];
        for (int i = 0; i < 500; i++) {
            labelOffset = 0;
            data[i][labelOffset++] = "Hap" + (char) ((i % 4) + 65);
            data[i][labelOffset++] = "" + i;
            data[i][labelOffset++] = "" + nf.format(random.nextGaussian());
            data[i][labelOffset++] = "" + nf.format(random.nextGaussian() + .3);
        }
        labelOffset = 0;
        data[25][labelOffset++] = "HapX" + (char) ((25 % 4) + 65);
        data[25][labelOffset++] = "BLAH";
        data[25][labelOffset++] = "NaN";
        data[25][labelOffset++] = "BLAH";

        myData = data;
    }

    public String getTableTitle() {
        return "Diversity estimates";
    }

    public int getRowCount() {
        return myData.length;
    }

    public int getElementCount() {
        return getRowCount() * getColumnCount();
    }

    public int getColumnCount() {
        return myColumnHeadings.length;
    }
}
