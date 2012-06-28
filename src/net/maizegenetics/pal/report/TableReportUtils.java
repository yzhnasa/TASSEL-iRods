package net.maizegenetics.pal.report;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;

import net.maizegenetics.util.DoubleFormat;

/**
 * Created by IntelliJ IDEA.
 * User: Ed
 * Date: Dec 7, 2004
 * Time: 2:10:50 PM
 * To change this template use File | Settings | File Templates.
 */
public class TableReportUtils {

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
                    if (theRow[i] == null) {
                        // do nothing
                    } else if (theRow[i] instanceof Double) {
                        bw.write(DoubleFormat.format((Double) theRow[i]));
                    } else {
                        bw.write(theRow[i].toString());
                    }
                }
                bw.write("\n");
            }

        } catch (Exception e) {
            e.printStackTrace();
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
