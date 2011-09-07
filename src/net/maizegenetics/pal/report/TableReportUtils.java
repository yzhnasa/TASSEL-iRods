package net.maizegenetics.pal.report;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;

/**
 * Created by IntelliJ IDEA.
 * User: Ed
 * Date: Dec 7, 2004
 * Time: 2:10:50 PM
 * To change this template use File | Settings | File Templates.
 */
public class TableReportUtils {

    public static String toDelimitedString(TableReport theTableSource, String delimit, String headerName, String headerValue) {
        Object[] colNames = theTableSource.getTableColumnNames();
        StringBuilder sb = new StringBuilder();

        int cols;
        cols = colNames.length;
        int rows = theTableSource.getRowCount();
        if (headerName != null) {
            sb.append(headerName + delimit);
        }
        for (int j = 0; j < cols; j++) {
            sb.append(colNames[j]);
            if (j < (cols - 1)) {
                sb.append(delimit);
            }
        }
        sb.append("\n");
        if (headerValue != null) {
            sb.append(headerValue + delimit);
        }
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                sb.append(theTableSource.getValueAt(i, j));
                if (j < (cols - 1)) {
                    sb.append(delimit);
                }
            }
            sb.append("\n");
        }
        return sb.toString();
    }

    public static String toDelimitedString(TableReport theTableSource, String delimit) {
        return toDelimitedString(theTableSource, delimit, null, null);
    }

    public static void saveDelimitedTableReport(TableReport theTableSource, String delimit, File saveFile) {

        if (saveFile == null) {
            return;
        }
        FileWriter fw = null;
        BufferedWriter bw = null;
        try {

            fw = new FileWriter(saveFile);
            bw = new BufferedWriter(fw);

            Object[] colNames = theTableSource.getTableColumnNames();
            for (int j = 0; j < colNames.length; j++) {
                if (j != 0) {
                    bw.write(delimit);
                }
                bw.write(colNames[j].toString());
            }
            bw.write("\n");

            for (int r = 0, n = theTableSource.getRowCount(); r < n; r++) {
                Object[] theRow = theTableSource.getRow(r);
                for (int i = 0; i < theRow.length; i++) {
                    if (i != 0) {
                        bw.write(delimit);
                    }
                    bw.write(theRow[i].toString());
                }
                bw.write("\n");
            }

        } catch (Exception e) {
            System.out.println("TableReportUtils: writeReport: problem writing file: " + e.getMessage());
        } finally {
            try {
                bw.close();
                fw.close();
            } catch (Exception e) {
                // do nothing
            }
        }

    }
}