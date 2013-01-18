package nutils;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.io.File;
import java.io.IOException;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYDotRenderer;
import org.jfree.chart.title.LegendTitle;
import org.jfree.chart.title.TextTitle;
import org.jfree.data.xy.DefaultXYDataset;

public class GraphUtils {

	// ========================================================================
	/** Saves a JFreeChart as an image, given certain parameters. */
	public static void saveChartAsPNG(String outFilenameRoot, JFreeChart theChart, int chartWidth, int chartHeight) {
		try {
			//ChartUtilities.saveChartAsJPEG(theFile, theChart, chartWidth, chartHeight);
			File theFilePNG = new File(outFilenameRoot + ".png");
			ChartUtilities.saveChartAsPNG(theFilePNG, theChart, chartWidth, chartHeight);
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(-1);
		}
	}
	
	// ========================================================================
	public static XYDotRenderer getXYDotRendererHelper(int size) {
		XYDotRenderer xyDotRend = new XYDotRenderer();
		xyDotRend.setDotWidth(size);
		xyDotRend.setDotHeight(size);
		return xyDotRend;
	}

	// ========================================================================
	public static void graphDistProbabilities(double[][] distProb, String chartTitleStr, String outFilenamePrefix) {
		DefaultXYDataset xyDataset = new DefaultXYDataset();
		Dimension chartDimension = new Dimension(800, 800);
		xyDataset.addSeries("Graph",  distProb);		
			
		boolean showLegend = true;
		
		JFreeChart theChart = 
			ChartFactory.createScatterPlot(chartTitleStr, "X-coordinate", "Probability", 
					xyDataset, PlotOrientation.VERTICAL, showLegend, false, false);
		XYPlot xyPlot = (XYPlot) theChart.getPlot();
		//xyPlot.getDomainAxis().setTickLabelsVisible(false);
			
		XYDotRenderer xyDotRend = getXYDotRendererHelper(4);
		xyDotRend.setSeriesPaint(0, Color.blue);
		xyDotRend.setSeriesPaint(1, Color.green);
		xyPlot.setRenderer(xyDotRend);
		xyPlot.setBackgroundPaint(Color.white);
		xyPlot.setDomainGridlinePaint(Color.gray);
		xyPlot.setRangeGridlinePaint(Color.gray);
		
		// Set the chart title font
		TextTitle title = theChart.getTitle();
		Font titleFont = new Font("Arial", Font.BOLD, 44);
		title.setFont(titleFont);
		
		// Set the chart legend font
		LegendTitle legendTitle = theChart.getLegend();
		Font legendFont = new Font("Arial", Font.BOLD, 24);
		if (showLegend) { legendTitle.setItemFont(legendFont); }
		
		// Set the chart range axis font
		Font rangeAxisLabelFont = new Font("Arial", Font.BOLD, 32);
		Font rangeAxisTickFont = new Font("Arial", Font.BOLD, 22);
		xyPlot.getRangeAxis().setLabelFont(rangeAxisLabelFont);		
		xyPlot.getRangeAxis().setTickLabelFont(rangeAxisTickFont);
		xyPlot.getDomainAxis().setLabelFont(rangeAxisLabelFont);
		xyPlot.getDomainAxis().setTickLabelFont(rangeAxisTickFont);
	
		// TODO - MAKE ALL		
		saveChartAsPNG(outFilenamePrefix + ".ChartDist", theChart, chartDimension.width, chartDimension.height);
	}
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
