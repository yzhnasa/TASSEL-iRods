/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package net.maizegenetics.baseplugins.chart;

import java.awt.BorderLayout;
import java.awt.Color;
import javax.swing.JComponent;
import net.maizegenetics.pal.report.TableReport;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.ChartRenderingInfo;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.labels.StandardXYToolTipGenerator;
import org.jfree.chart.labels.XYToolTipGenerator;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.StandardXYItemRenderer;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.data.function.Function2D;
import org.jfree.data.function.LineFunction2D;
import org.jfree.data.general.DatasetUtilities;
import org.jfree.data.xy.XYDataset;

/**
 *
 * @author yz79
 */
public class XYScatterAndLinePanel extends BasicChartPanel {

    ChartPanel myChartPanel;
    TableReportQQDataset dataset;

    TableReport myTableReport;

    public XYScatterAndLinePanel(TableReport theTable) {
        myTableReport = theTable;
        try {
            dataset = new TableReportQQDataset(theTable);
            chart = createChart(dataset);
            myChartPanel = new ChartPanel(chart);
            myChartPanel.setPreferredSize(new java.awt.Dimension(500, 270));
            myTableReport = theTable;
            jbInit();
        }
        catch(Exception ex) {
            ex.printStackTrace();
        }
    }

    private void jbInit() throws Exception {
        this.add(myChartPanel, BorderLayout.CENTER);
    }


    public JFreeChart createChart(TableReportQQDataset dataset) {
        String name="Please select numeric variables";
        String xName="X";
        String y1Name="Y";
        String y2Name="Y2";
        if(dataset!=null) {
            xName=dataset.getXName();
            y1Name=dataset.getSeriesName(0);
            name=xName+" vs. "+y1Name;
            if(dataset.getSeriesCount()==2) {
                y2Name=dataset.getSeriesName(1);
                name+=" and "+y2Name;
            }
        }
        chart = ChartFactory.createScatterPlot(
        name,
        xName,y1Name,
        dataset,
        PlotOrientation.VERTICAL,
        true,
        true,
        false
        );
        chart.getXYPlot().setForegroundAlpha(0.75f);
        chart.getXYPlot().getRenderer().setToolTipGenerator(new XYAndLineToolTipGenerator());
        createLine(chart.getXYPlot(), dataset, 0, Color.BLACK);
        return chart;
    }

    private void createLine(XYPlot plot, XYDataset data, int series, Color theColor) {
        Function2D curve = new LineFunction2D(0, 1);
        double max=DatasetUtilities.findMaximumDomainValue(data).doubleValue();
        double min=DatasetUtilities.findMinimumDomainValue(data).doubleValue();
        XYDataset regressionData = DatasetUtilities.sampleFunction2D(
            curve,
            min,max, 2,
            "Expected Values"
        );
        int datasetCount=plot.getDatasetCount();
        plot.setDataset(datasetCount, regressionData);
        XYItemRenderer renderer2 = new StandardXYItemRenderer(StandardXYItemRenderer.LINES);
        renderer2.setSeriesPaint(0, theColor);
        plot.setRenderer(datasetCount, renderer2);
        setAxis(plot, max, min);
    }

    private void setAxis(XYPlot plot, double domainMax, double domainMin) {
        ValueAxis xAxis = plot.getDomainAxis();
        ValueAxis yAxis = plot.getRangeAxis();
        double xMax = domainMax;
        double xMin = domainMin;
        double yMax = yAxis.getUpperBound();
        double yMin = yAxis.getLowerBound();
        xAxis.setRange(0, xMax);
        yAxis.setRange(0, yMax);
        plot.setDomainAxis(xAxis);
        plot.setRangeAxis(yAxis);
    }

//    private void setToolTips(XYPlot plot, XYDataset data) {
//        XYAndLineToolTipGenerator tipGenerator = new XYAndLineToolTipGenerator();
//        XYItemRenderer itemRenderer = new StandardXYItemRenderer(StandardXYItemRenderer.SHAPES);
//        itemRenderer.setToolTipGenerator(tipGenerator);
//        plot.setRenderer(itemRenderer);
//    }

    public JComponent getMainComponent() {
        return myChartPanel;
    }

}