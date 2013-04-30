package lohcate.clustering;

import genomeEnums.Chrom;

import java.awt.Color;
import java.awt.Font;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Rectangle2D;
import java.awt.geom.Ellipse2D.Double;
import java.io.File;
import java.util.ArrayList;
import java.util.EnumMap;
import java.util.Set;

import lohcateEnums.ClusterType;
import lohcateEnums.ColorPastel;
import nutils.ArrayUtils;
import nutils.ArrayUtils.ParallelArrayDoubleDynamic;
import nutils.EnumMapSafe;
import nutils.GraphUtils;
import nutils.PrimitiveWrapper;
import nutils.PrimitiveWrapper.WInteger;
import nutils.counter.DynamicBucketCounter;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.chart.renderer.xy.XYShapeRenderer;
import org.jfree.chart.title.LegendTitle;
import org.jfree.data.xy.DefaultXYDataset;
import org.jfree.data.xy.XYDataset;

import com.carrotsearch.hppc.DoubleArrayList;

public class ClusteringPlotting {

	// ========================================================================
	/** Plots the VAF (variant allele frequency) of the normal tissue comapred to the tumor tissue. */
	public static void plotVAFComparison(XYDataset xyDataset, String outFilenameRoot, String sampleName) {
		String xAxisLabel = "VAF Tumor";
		String yAxisLabel = "VAF Normal";
		String title = "VAF Comparison: " + sampleName;
		
		JFreeChart theChart = 
				ChartFactory.createScatterPlot(title, xAxisLabel, yAxisLabel, xyDataset, PlotOrientation.VERTICAL, true, false, false);
		XYPlot xyPlot = (XYPlot) theChart.getPlot();
		
		XYItemRenderer itemRenderer = getXYItemRendererHelper(5);		
		setSeriesPaintPerCluster(itemRenderer);
		xyPlot.setRenderer(itemRenderer);		
		
		xyPlot.setBackgroundPaint(ColorPastel.Gray_15.getColor());
		xyPlot.setDomainGridlinePaint(Color.white);
		xyPlot.setRangeGridlinePaint(Color.white);		
		
		LegendTitle legendTitle = theChart.getLegend();
		legendTitle.setID("Clusters");
		legendTitle.setItemFont(new Font("Arial", Font.BOLD, 20));
		
		Font rangeAxisLabelFont = new Font("Arial", Font.BOLD, 20);
		Font rangeAxisTickFont = new Font("Arial", Font.BOLD, 20);
		xyPlot.getRangeAxis().setRange(0, 1.02);
		xyPlot.getDomainAxis().setRange(0, 1.02);
		xyPlot.getRangeAxis().setLabelFont(rangeAxisLabelFont);		
		xyPlot.getRangeAxis().setTickLabelFont(rangeAxisTickFont);
		xyPlot.getDomainAxis().setLabelFont(rangeAxisLabelFont);
		xyPlot.getDomainAxis().setTickLabelFont(rangeAxisTickFont);		
	
		// Now write the plot		
		int width  = 900;
		int height = 900;
		GraphUtils.saveChartAsPNG(outFilenameRoot, theChart, width, height);
	}
	
	// ========================================================================
	public static<E extends Enum<E>> void plotEventsByCoordinateAcrossSamples(EnumMapSafe<E, ParallelArrayDoubleDynamic> eventsByCoordinate, String outDir) {
		DefaultXYDataset xyDataset = new DefaultXYDataset();
		
		Set<E> enumKeys = eventsByCoordinate.keySet();
		for (E enumKey : enumKeys) {			
			xyDataset.addSeries(enumKey.toString(), eventsByCoordinate.get(enumKey).toArrays());
		}
		
		String xAxisLabel = "Position";
		String yAxisLabel = "Sample ID with Event";
		String title = "HeatMap for Events Per Sample";
		
		JFreeChart theChart = 
				ChartFactory.createScatterPlot(title, xAxisLabel, yAxisLabel, xyDataset, PlotOrientation.VERTICAL, true, false, false);
		XYPlot xyPlot = (XYPlot) theChart.getPlot();
		
		XYItemRenderer itemRenderer = ClusteringPlotting.getXYItemRendererHelper(1);		
		setSeriesPaintPerCluster(itemRenderer);
		xyPlot.setRenderer(itemRenderer);		
		
		xyPlot.setBackgroundPaint(ColorPastel.White.getColor());
		xyPlot.setDomainGridlinePaint(Color.white);
		xyPlot.setRangeGridlinePaint(Color.white);		
		
		LegendTitle legendTitle = theChart.getLegend();
		legendTitle.setID("Event");
		legendTitle.setItemFont(new Font("Arial", Font.BOLD, 20));
	
		//xyPlot.getRangeAxis().setRange(0, 1.02);		
		
		Font rangeAxisLabelFont = new Font("Arial", Font.BOLD, 20);
		Font rangeAxisTickFont = new Font("Arial", Font.BOLD, 20);
		xyPlot.getRangeAxis().setLabelFont(rangeAxisLabelFont);		
		xyPlot.getRangeAxis().setTickLabelFont(rangeAxisTickFont);
		xyPlot.getDomainAxis().setLabelFont(rangeAxisLabelFont);
		xyPlot.getDomainAxis().setTickLabelFont(rangeAxisTickFont);	
		GraphUtils.saveChartAsPNG(outDir + File.separator + "HeapMap.AllEvents.AllSamples", theChart, 2400, 800);


	}

	// ========================================================================
	/** Plots the event recurrence genome-wide. */
	public static void plotRecurrenceGenomeWideByEvent(EnumMapSafe<Chrom, DynamicBucketCounter> eventCount, String outDir, ClusterType eventType) {
		DoubleArrayList positionsGenomeWideAllSamples = new DoubleArrayList();
		DoubleArrayList lohCountAllSamples            = new DoubleArrayList();
		double lastPositionOnPrevChrom = 0;
		
		DoubleArrayList chromBoundaryXValue = new DoubleArrayList();
		DoubleArrayList chromBoundaryYValue = new DoubleArrayList();
		PrimitiveWrapper.WInteger lastKey = new PrimitiveWrapper.WInteger(0);		
		
		// First get the max count for the chromosome boundaries
		int maxCountAllChromosomes = 1;
		for (Chrom chrom : Chrom.Autosomes) {
			maxCountAllChromosomes = Math.max(eventCount.get(chrom).getCountMax(), maxCountAllChromosomes);
		}
		
		for (Chrom chrom : Chrom.values()) {
			if (chrom.isAutosomal()) {
				DoubleArrayList[] positionAndCount = eventCount.get(chrom).toArrayListDouble();
				for (int i = 0; i < positionAndCount[0].size(); i++) {
					positionAndCount[0].set(i, positionAndCount[0].get(i) + lastPositionOnPrevChrom);
				}
				
				positionsGenomeWideAllSamples.addAll(positionAndCount[0]);
				lohCountAllSamples.addAll           (positionAndCount[1]);
				
				boolean lastKeyExists = eventCount.get(chrom).getKeyLast(lastKey);
				lastPositionOnPrevChrom += (lastKeyExists) ? lastKey.mInt : 30000000;				
				
				double increment = (maxCountAllChromosomes) / 50.0;
				for (double d = 0; d <= maxCountAllChromosomes; d += increment) {
					chromBoundaryXValue.add(lastPositionOnPrevChrom);
					chromBoundaryYValue.add(d);
				}
			}
		}
		DefaultXYDataset posAndEventCountDataset = new DefaultXYDataset();
		posAndEventCountDataset.addSeries(eventType.name(), ArrayUtils.combineTwoDynamicArraysIntoOneStatic(positionsGenomeWideAllSamples, lohCountAllSamples));
		posAndEventCountDataset.addSeries("Boundary",       ArrayUtils.combineTwoDynamicArraysIntoOneStatic(chromBoundaryXValue, chromBoundaryYValue));
		ClusteringPlotting.plotEventSampleRecurrence(posAndEventCountDataset, outDir + File.separator + "All_Samples." + eventType.name(), getColorForEvent(eventType).getColor(), ColorPastel.Gray_80.getColor());
	}
	
	// ========================================================================
	public static void plotRecurrenceGenomeWide(EnumMapSafe<ClusterType, EnumMapSafe<Chrom, DynamicBucketCounter>> eventCounts, String outDir) {
		for (ClusterType event : ClusterType.values()) {
			plotRecurrenceGenomeWideByEvent(eventCounts.get(event), outDir, event);
		}
	}

	// ========================================================================
	public static void plotEventSampleRecurrence(XYDataset xyDataset, String outFilenameRoot, Color seriesColor, Color boundaryColor) {
		String xAxisLabel = "Position";
		String yAxisLabel = "# Samples with Event";
		String title = "Recurrence over Samples for Event: " + xyDataset.getSeriesKey(0);				
		
		JFreeChart theChart = 
				ChartFactory.createScatterPlot(title, xAxisLabel, yAxisLabel, xyDataset, PlotOrientation.VERTICAL, true, false, false);
		XYPlot xyPlot = (XYPlot) theChart.getPlot();
		
		XYItemRenderer itemRenderer = ClusteringPlotting.getXYItemRendererHelper(5);		
		//setSeriesPaintPerCluster(itemRenderer);
		itemRenderer.setSeriesPaint(0, seriesColor);
		itemRenderer.setSeriesPaint(1, boundaryColor);
		xyPlot.setRenderer(itemRenderer);		
		
		xyPlot.setBackgroundPaint(ColorPastel.White.getColor());
		xyPlot.setDomainGridlinePaint(Color.white);
		xyPlot.setRangeGridlinePaint(Color.white);		
		
		LegendTitle legendTitle = theChart.getLegend();
		legendTitle.setID("Event");
		legendTitle.setItemFont(new Font("Arial", Font.BOLD, 20));
	
		//xyPlot.getRangeAxis().setRange(0, 1.02);		
		
		Font rangeAxisLabelFont = new Font("Arial", Font.BOLD, 20);
		Font rangeAxisTickFont = new Font("Arial", Font.BOLD, 20);
		xyPlot.getRangeAxis().setLabelFont(rangeAxisLabelFont);		
		xyPlot.getRangeAxis().setTickLabelFont(rangeAxisTickFont);
		xyPlot.getDomainAxis().setLabelFont(rangeAxisLabelFont);
		xyPlot.getDomainAxis().setTickLabelFont(rangeAxisTickFont);	
		GraphUtils.saveChartAsPNG(outFilenameRoot, theChart, 2400, 800);
		
	}

	// ========================================================================
	/** Plots the VAF (variant allele frequency) of a tissue vs the genomic position. */
	public static void plotCopyNumGenomeWide(XYDataset xyDataset, String outFilenameRoot, String sampleName) {
		String xAxisLabel = "Position";
		String yAxisLabel = "Copy Number";
		String title = "Copy Number: " + sampleName;				
		
		JFreeChart theChart = 
				ChartFactory.createScatterPlot(title, xAxisLabel, yAxisLabel, xyDataset, PlotOrientation.VERTICAL, true, false, false);
		XYPlot xyPlot = (XYPlot) theChart.getPlot();
		
		XYItemRenderer itemRenderer = ClusteringPlotting.getXYItemRendererHelper(5);		
		setSeriesPaintPerCluster(itemRenderer);
		xyPlot.setRenderer(itemRenderer);		
		
		xyPlot.setBackgroundPaint(ColorPastel.White.getColor());
		xyPlot.setDomainGridlinePaint(Color.white);
		xyPlot.setRangeGridlinePaint(Color.white);		
		
		LegendTitle legendTitle = theChart.getLegend();
		legendTitle.setID("Clusters");
		legendTitle.setItemFont(new Font("Arial", Font.BOLD, 20));
	
		//xyPlot.getRangeAxis().setRange(0, 1.02);		
		
		Font rangeAxisLabelFont = new Font("Arial", Font.BOLD, 20);
		Font rangeAxisTickFont = new Font("Arial", Font.BOLD, 20);
		xyPlot.getRangeAxis().setLabelFont(rangeAxisLabelFont);		
		xyPlot.getRangeAxis().setTickLabelFont(rangeAxisTickFont);
		xyPlot.getDomainAxis().setLabelFont(rangeAxisLabelFont);
		xyPlot.getDomainAxis().setTickLabelFont(rangeAxisTickFont);	
		GraphUtils.saveChartAsPNG(outFilenameRoot, theChart, 2400, 800);
	}

	// ========================================================================
	/** Plots the VAF (variant allele frequency) of a tissue vs the genomic position. */
	public static void plotVAFGenomeWide(XYDataset xyDataset, String outFilenameRoot, String sampleName, boolean isTumor) {
		String xAxisLabel = "Position";
		String yAxisLabel = "VAF";
		String title = "VAF GenomeWide: " + sampleName + (isTumor ? " [Tumor]" : " [Normal]");
		
		JFreeChart theChart = 
				ChartFactory.createScatterPlot(title, xAxisLabel, yAxisLabel, xyDataset, PlotOrientation.VERTICAL, true, false, false);
		XYPlot xyPlot = (XYPlot) theChart.getPlot();
		
		XYItemRenderer itemRenderer = ClusteringPlotting.getXYItemRendererHelper(5);		
		setSeriesPaintPerCluster(itemRenderer);
		xyPlot.setRenderer(itemRenderer);		
		
		xyPlot.setBackgroundPaint(ColorPastel.White.getColor());
		xyPlot.setDomainGridlinePaint(Color.white);
		xyPlot.setRangeGridlinePaint(Color.white);		
		
		LegendTitle legendTitle = theChart.getLegend();
		legendTitle.setID("Clusters");
		legendTitle.setItemFont(new Font("Arial", Font.BOLD, 20));
	
		xyPlot.getRangeAxis().setRange(0, 1.02);		
		
		Font rangeAxisLabelFont = new Font("Arial", Font.BOLD, 20);
		Font rangeAxisTickFont = new Font("Arial", Font.BOLD, 20);
		xyPlot.getRangeAxis().setLabelFont(rangeAxisLabelFont);		
		xyPlot.getRangeAxis().setTickLabelFont(rangeAxisTickFont);
		xyPlot.getDomainAxis().setLabelFont(rangeAxisLabelFont);
		xyPlot.getDomainAxis().setTickLabelFont(rangeAxisTickFont);	
		GraphUtils.saveChartAsPNG(outFilenameRoot, theChart, 2400, 800);
	}

	// ========================================================================
	public static XYItemRenderer getXYItemRendererHelper(int size) {
		XYShapeRenderer xyShapeRend = new XYShapeRenderer();				
		
		//XYDotRenderer xyDotRend = new XYDotRenderer();
		Ellipse2D.Double theEllipse = new Ellipse2D.Double(0, 0, size, size);
		Rectangle2D.Double theRect = new Rectangle2D.Double(0, 0, size, size);  
		//xyDotRend.setBaseShape(theEllipse, true);
		//xyShapeRend.setBaseShape(theEllipse);
		xyShapeRend.setBaseShape(theRect);
		//xyDotRend.setDotWidth(size);
		//xyDotRend.setDotHeight(size);
		//return xyDotRend;
		return xyShapeRend;
	}

	// ========================================================================
	static void setSeriesPaintPerCluster(XYItemRenderer itemRenderer) {
		boolean allGray = false;
		boolean makeNonEventsInvisible = false;
			
		if (allGray) {
			for (ClusterType eventType : ClusterType.values()) {
				itemRenderer.setSeriesPaint(eventType.ordinal(), ColorPastel.Gray_50.getColor());
			}
			if (makeNonEventsInvisible) {
				itemRenderer.setSeriesPaint(ClusterType.HETGermline.ordinal(), ColorPastel.Gray_15.getColor());
				itemRenderer.setSeriesPaint(ClusterType.HETSomatic.ordinal(),  ColorPastel.Gray_15.getColor());
				itemRenderer.setSeriesPaint(ClusterType.Noise.ordinal(),       ColorPastel.Gray_15.getColor());
				itemRenderer.setSeriesPaint(ClusterType.Null.ordinal(),        ColorPastel.Gray_15.getColor());
			}
		} else {
			for (ClusterType eventType : ClusterType.values()) {
				itemRenderer.setSeriesPaint(eventType.ordinal(), getColorForEvent(eventType).getColor());
			}
		}
		itemRenderer.setSeriesPaint(ClusterType.Null.ordinal() + 1, ColorPastel.Black.getColor());
	}

	// ========================================================================
	public static ColorPastel getColorForEvent(ClusterType eventType) {
		switch(eventType) {
		case GainGermline: return ColorPastel.Violet;
		case GainSomatic:  return ColorPastel.Dark_Red;
		case LOH:          return ColorPastel.RGB_Blue;
		case cnLOH:        return ColorPastel.CMYK_Yellow;
		case HETGermline:  return ColorPastel.Gray_60;
		case HETSomatic:   return ColorPastel.Red_Orange;
		case Noise:        return ColorPastel.Dark_Pea_Green;
		case Null:         return ColorPastel.Gray_30;	
		default:           return null;
		}
	}
}
