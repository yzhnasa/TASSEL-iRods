package net.maizegenetics.util;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.regex.Pattern;

import net.maizegenetics.util.DoubleFormat;
import net.maizegenetics.util.ExceptionUtils;
import net.maizegenetics.util.Utils;

import org.apache.log4j.Logger;

/**
 * @author terry
 */
public class TableReportUtils {

    private static final Logger myLogger = Logger.getLogger(TableReportUtils.class);

    public static void saveDelimitedTableReport(TableReport theTableSource, String delimit, File saveFile) {

        if (saveFile == null) {
            return;
        }

        BufferedWriter bw = null;
        try {

            bw = Utils.getBufferedWriter(saveFile);

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
            } catch (Exception e) {
                // do nothing
            }
        }

    }

    public static TableReport readDelimitedTableReport(String saveFile, String delimit) {

        myLogger.info("readDelimitedTableReport: Reading: " + saveFile);
        int numLines = Utils.getNumberLines(saveFile) - 1;
        myLogger.info("readDelimitedTableReport: Num Lines (Not including header): " + numLines);

        Pattern delimitPattern = Pattern.compile(delimit);
        BufferedReader br = null;
        try {
            br = Utils.getBufferedReader(saveFile);
            String[] columnHeaders = delimitPattern.split(br.readLine().trim());

            int numThreads = Runtime.getRuntime().availableProcessors();
            ExecutorService pool = Executors.newFixedThreadPool(numThreads);
            String[][] data = new String[numLines][];
            int maxNumLinesPerThread = 100000;
            for (int i = 0; i < numLines; i = i + maxNumLinesPerThread) {
                int numLinesForThread = Math.min(maxNumLinesPerThread, numLines - i);
                String[] lines = new String[numLinesForThread];
                for (int j = 0; j < numLinesForThread; j++) {
                    lines[j] = br.readLine().trim();
                }
                pool.execute(new SplitTableReportString(data, i, lines, delimitPattern));
            }
            pool.shutdown();
            if (!pool.awaitTermination(6000, TimeUnit.SECONDS)) {
                throw new IllegalStateException("TableReportUtils: readDelimitedTableReport: processing threads timed out.");
            }
            return new SimpleTableReport(saveFile, columnHeaders, data);
        } catch (Exception e) {
            e.printStackTrace();
            throw new IllegalArgumentException("Problem creating TableReport: " + saveFile + ": " + ExceptionUtils.getExceptionCauses(e));
        } finally {
            try {
                br.close();
            } catch (Exception ex) {
                // do nothing
            }
        }

    }

    private static class SplitTableReportString implements Runnable {

        private final String[][] myData;
        private int myLineNum;
        private final String[] myLines;
        private final Pattern myPattern;

        public SplitTableReportString(String[][] data, int lineNum, String[] lines, Pattern pattern) {
            myData = data;
            myLineNum = lineNum;
            myLines = lines;
            myPattern = pattern;
        }

        @Override
        public void run() {
            for (int i = 0, n = myLines.length; i < n; i++) {
                myData[myLineNum++] = myPattern.split(myLines[i]);
            }
        }
    }
}
