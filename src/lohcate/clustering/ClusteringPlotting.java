package lohcate.clustering;

import genomeEnums.Chrom;
import genomeUtils.MaxPositionsOnChromosome;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.awt.Rectangle;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Line2D;
import java.awt.geom.Rectangle2D;
import java.awt.geom.Ellipse2D.Double;
import java.io.File;
import java.util.ArrayList;
import java.util.EnumMap;
import java.util.Iterator;
import java.util.Set;

import lohcateEnums.EventType;
import lohcateEnums.ColorPastel;
import nutils.ArrayUtils;
import nutils.ArgumentParserUtils.InputParameterBoolean;
import nutils.ArrayUtils.ParallelArrayDoubleDynamic;
import nutils.ArgumentParserUtils;
import nutils.ControlFlagBool;
import nutils.EnumMapSafe;
import nutils.GraphUtils;
import nutils.PrimitiveWrapper;
import nutils.PrimitiveWrapper.WInteger;
import nutils.counter.DynamicBucketCounter;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.LegendItem;
import org.jfree.chart.LegendItemCollection;
import org.jfree.chart.annotations.XYBoxAnnotation;
import org.jfree.chart.annotations.XYTextAnnotation;
import org.jfree.chart.axis.NumberTickUnit;
import org.jfree.chart.axis.TickUnit;
import org.jfree.chart.axis.TickUnits;
import org.jfree.chart.plot.IntervalMarker;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.ValueMarker;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.StackedXYAreaRenderer;
import org.jfree.chart.renderer.xy.StackedXYAreaRenderer2;
import org.jfree.chart.renderer.xy.StackedXYBarRenderer;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.chart.renderer.xy.XYShapeRenderer;
import org.jfree.chart.title.LegendTitle;
import org.jfree.chart.title.TextTitle;
import org.jfree.data.xy.DefaultXYDataset;
import org.jfree.data.xy.XYBarDataset;
import org.jfree.data.xy.XYDataset;
import org.jfree.ui.HorizontalAlignment;
import org.jfree.ui.Layer;

import com.carrotsearch.hppc.DoubleArrayList;

public class ClusteringPlotting {

	private static final String ChromTitleString = buildChromString();
	public static final double GenomeWidePositionDivisor = 1_000_000; 
	private static final String GenomeWideXAxisLabel = "Genomic Position (MB)";
	
	private static final int CopyNumberVAFPlotPointSize = 4;
	
	// ========================================================================
	/** Plots the VAF (variant allele frequency) of the normal tissue comapred to the tumor tissue. */
	public static void plotVAFComparison(XYDataset xyDataset, String outFilenameRoot, String sampleName) {
		String xAxisLabel = "VAF Tumor";
		String yAxisLabel = "VAF Normal";
		String title = "VAF Comparison: " + sampleName;
		
		JFreeChart theChart = 
				ChartFactory.createScatterPlot(title, xAxisLabel, yAxisLabel, xyDataset, PlotOrientation.VERTICAL, true, false, false);
		XYPlot xyPlot = (XYPlot) theChart.getPlot();
		
		XYItemRenderer itemRenderer = getXYItemRendererHelper(CopyNumberVAFPlotPointSize);		
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
		int width  = 1000;
		int height = 1000;
		GraphUtils.saveChartAsPNG(outFilenameRoot, theChart, width, height);
	}
	
	// ========================================================================
	public static<E extends Enum<E>> void plotEventsByCoordinateAcrossSamples(EnumMapSafe<E, ParallelArrayDoubleDynamic> eventsByCoordinate, MaxPositionsOnChromosome maxPosByChrom, int numSamples, String outDir) {
		DefaultXYDataset xyDataset = new DefaultXYDataset();
		
		Set<E> enumKeys = eventsByCoordinate.keySet();		
		for (E enumKey : enumKeys) {			
			double[][] eventsXY = eventsByCoordinate.get(enumKey).toArrays();
			xyDataset.addSeries(enumKey.toString(), eventsXY);
		}
		 
		ParallelArrayDoubleDynamic chromBoundaryXY = new ParallelArrayDoubleDynamic();
		double incr = numSamples / 200.0;
		for (Chrom chrom : Chrom.Autosomes) {
			long maxPos = maxPosByChrom.getMaxPositionForChrom(chrom);
			System.out.println("MaxPos:\t" + chrom + "\t" + maxPos);
			for (double yCoordinate = 0; yCoordinate <= numSamples; yCoordinate += incr) {
				chromBoundaryXY.add(maxPos, yCoordinate);
			}			
		}
		xyDataset.addSeries("Boundary", chromBoundaryXY.toArrays());
		
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
		TickUnit tickUnit =  new NumberTickUnit(1.0);
		TickUnits tickUnits = new TickUnits();
		tickUnits.add(tickUnit);		
		xyPlot.getRangeAxis().setStandardTickUnits(tickUnits);
		xyPlot.getDomainAxis().setLabelFont(rangeAxisLabelFont);
		xyPlot.getDomainAxis().setTickLabelFont(rangeAxisTickFont);	
		GraphUtils.saveChartAsPNG(outDir + File.separator + "HeapMap.AllEvents.AllSamples", theChart, 2400, 1000);


	}

	// ========================================================================
	/** Plots the event recurrence genome-wide. */
	public static void plotRecurrenceGenomeWideByEvent(EnumMapSafe<Chrom, DynamicBucketCounter> eventCount, String outDir, EventType eventType) {
		DoubleArrayList positionsGenomeWideAllSamples = new DoubleArrayList();
		DoubleArrayList eventCountAllSamples          = new DoubleArrayList();		
			
		// First get the max count for the chromosome boundaries
		int maxCountAllChromosomes = 1;
		for (Chrom chrom : Chrom.Autosomes) {
			maxCountAllChromosomes = Math.max(eventCount.get(chrom).getCountMax(), maxCountAllChromosomes);
		}
		
		for (Chrom chrom : Chrom.values()) {
			if (chrom.isAutosomal()) {				
				DoubleArrayList[] positionAndCount = eventCount.get(chrom).toArrayListDouble();
				for (int i = 0; i < positionAndCount[0].size(); i++) {
					double newPosAdj = (positionAndCount[0].get(i) + chrom.getGenomeWidePositionStart()) / GenomeWidePositionDivisor; 
					positionAndCount[0].set(i, newPosAdj);
				}
				
				positionsGenomeWideAllSamples.addAll(positionAndCount[0]);
				eventCountAllSamples.addAll         (positionAndCount[1]);				
			}
		}
		DefaultXYDataset posAndEventCountDataset = new DefaultXYDataset();
		posAndEventCountDataset.addSeries(eventType.name(), ArrayUtils.combineTwoDynamicArraysIntoOneStatic(positionsGenomeWideAllSamples, eventCountAllSamples));
		ClusteringPlotting.plotEventSampleRecurrence(posAndEventCountDataset, outDir + File.separator + "All_Samples." + eventType.name(), getColorForEvent(eventType).getColor(), maxCountAllChromosomes);
	}
	
	// ========================================================================
	public static void plotRecurrenceGenomeWide(EnumMapSafe<EventType, EnumMapSafe<Chrom, DynamicBucketCounter>> eventCounts, String outDir) {
		for (EventType event : EventType.values()) {
			plotRecurrenceGenomeWideByEvent(eventCounts.get(event), outDir, event);
		}
	}

	// ========================================================================
	public static void plotEventSampleRecurrence(XYDataset xyDataset, String outFilenameRoot, Color seriesColor, int maxYValue) {
		String xAxisLabel = GenomeWideXAxisLabel;
		String yAxisLabel = "# Samples with Event";
		String title = "Recurrence over Samples for Event: " + xyDataset.getSeriesKey(0);				
		
		JFreeChart theChart = 
				ChartFactory.createScatterPlot(title, xAxisLabel, yAxisLabel, xyDataset, PlotOrientation.VERTICAL, true, false, false);
		XYPlot xyPlot = (XYPlot) theChart.getPlot();
		//theChart.addSubtitle(getChromSubtitle());
		
		XYItemRenderer itemRenderer = ClusteringPlotting.getXYItemRendererHelper(5);		
		//setSeriesPaintPerCluster(itemRenderer);
		itemRenderer.setSeriesPaint(0, seriesColor);		
		xyPlot.setRenderer(itemRenderer);		
		
		xyPlot.setBackgroundPaint(ColorPastel.White.getColor());
		xyPlot.setDomainGridlinePaint(Color.white);
		xyPlot.setRangeGridlinePaint(Color.white);		
		
		double maxYRangeAxis = Math.max(maxYValue, 10);
		formatXYPlotBackground(xyPlot, maxYRangeAxis, false, 0, 1);
		
		
		formatLegendItems(xyPlot);
		LegendTitle legendTitle = theChart.getLegend();
		legendTitle.setID("Event");
		legendTitle.setItemFont(new Font("Arial", Font.BOLD, 32));

		//xyPlot.getRangeAxis().setRange(0, 1.02);		
		
		TickUnits tunits = new TickUnits();
		tunits.add(new NumberTickUnit(2));
		//xyPlot.getRangeAxis().setRange(0, 62);
		xyPlot.getRangeAxis().setStandardTickUnits(tunits);
		
		Font rangeAxisLabelFont = new Font("Arial", Font.BOLD, 32);
		Font rangeAxisTickFont = new Font("Arial", Font.BOLD, 20);
		xyPlot.getRangeAxis().setLabelFont(rangeAxisLabelFont);		
		xyPlot.getRangeAxis().setTickLabelFont(rangeAxisTickFont);
		xyPlot.getRangeAxis().setRange(0, maxYRangeAxis + 2);
		xyPlot.getDomainAxis().setLabelFont(rangeAxisLabelFont);
		xyPlot.getDomainAxis().setTickLabelFont(rangeAxisTickFont);	
		GraphUtils.saveChartAsPNG(outFilenameRoot, theChart, 2400, 800);
		
	}

	// ========================================================================	
	
	// ========================================================================
	public static double getMaxY(XYDataset xyDataset) {
		double maxYValue = -java.lang.Double.MAX_VALUE;
		int numSeries = xyDataset.getSeriesCount();
		for (int seriesIndex = 0; seriesIndex < numSeries; seriesIndex++) {
			int itemCount = xyDataset.getItemCount(seriesIndex);
			for (int itemIndex = 0; itemIndex < itemCount; itemIndex++) {
				double yVal = xyDataset.getYValue(seriesIndex, itemIndex);
				if (yVal > maxYValue) { maxYValue = yVal; }
			}
		}
		return maxYValue;
	}
	
	// ========================================================================
	/** Plots the VAF (variant allele frequency) of a tissue vs the genomic position. */
	public static void plotCopyNumGenomeWide(XYDataset xyDataset, String outFilenameRoot, String sampleName) {
		String xAxisLabel = GenomeWideXAxisLabel;
		String yAxisLabel = "Total Copy Number [Tumor]";
		String title = "Total Copy Number [Tumor]: " + sampleName;				
				
		JFreeChart theChart = 
				ChartFactory.createScatterPlot(title, xAxisLabel, yAxisLabel, xyDataset, PlotOrientation.VERTICAL, true, false, false);
		XYPlot xyPlot = (XYPlot) theChart.getPlot();
		//theChart.addSubtitle(getChromSubtitle());
		
		XYItemRenderer itemRenderer = ClusteringPlotting.getXYItemRendererHelper(CopyNumberVAFPlotPointSize);		
		setSeriesPaintPerCluster(itemRenderer);
		xyPlot.setRenderer(itemRenderer);		
		xyPlot.addRangeMarker(new ValueMarker(2.0));
		
		xyPlot.setBackgroundPaint(ColorPastel.White.getColor());
		xyPlot.setDomainGridlinePaint(Color.white);
		xyPlot.setRangeGridlinePaint(Color.white);		
		
		double maxCopyNumber = Math.max(4, getMaxY(xyDataset));
		formatXYPlotBackground(xyPlot, maxCopyNumber, false, 0, 0.10);
		
		formatLegendItems(xyPlot);
		LegendTitle legendTitle = theChart.getLegend();
		legendTitle.setID("Events");		
		legendTitle.setItemFont(new Font("Arial", Font.BOLD, 32));
		
		//xyPlot.getRangeAxis().setRange(0, 1.02);		
		
		Font rangeAxisLabelFont = new Font("Arial", Font.BOLD, 32);
		Font rangeAxisTickFont = new Font("Arial", Font.BOLD, 20);
		xyPlot.getRangeAxis().setLabelFont(rangeAxisLabelFont);		
		xyPlot.getRangeAxis().setTickLabelFont(rangeAxisTickFont);
		xyPlot.getRangeAxis().setRange(0, maxCopyNumber + 0.25);
		xyPlot.getDomainAxis().setLabelFont(rangeAxisLabelFont);
		xyPlot.getDomainAxis().setTickLabelFont(rangeAxisTickFont);	
		xyPlot.getDomainAxis().setRange(0, Chrom.c22.calculateGenomeWidePositionEnd() / GenomeWidePositionDivisor);
		GraphUtils.saveChartAsPNG(outFilenameRoot, theChart, 2400, 800);
	}

	// ========================================================================
	public static String buildChromString() {
		StringBuilder chromSB = new StringBuilder(2048);
		
		/* For font size 20 
		int[] numSpaces = new int[] {
				27,  //1 
				27,  //2 
				23,  //3 
				23,  //4
				20,  //5 
				19,  //6 
				17,  //7 
				16,  //8
				14,  //9 
				14, //10 
				15, //11
				13, //12
				10, //13 
				9,  //14 
				7,  //15 
				7,  //16
				7,  //17 
				5,  //18 
				4,  //19 
				4,  //20
				3,  //21 
				2
				
		};
		*/
		
		// For font size 24
		int[] numSpaces = new int[] {
				25,  //1 
				21,  //2 
				19,  //3 
				18,  //4
				18,  //5 
				15,  //6 
				15,  //7 
				14,  //8
				12,  //9 
				11, //10 
				11, //11
				11, //12
				8, //13 
				7,  //14 
				6,  //15 
				5,  //16
				5,  //17 
				5,  //18 
				3,  //19 
				3,  //20
				3,  //21 
				2
				
		};
		
		
		for (Chrom chrom : Chrom.Autosomes) {
			chromSB.append(chrom.ordinal());
			for (int i = 0; i < numSpaces[chrom.ordinal() - 1]; i++) {
				chromSB.append(" ");
			}
		}		
		return chromSB.toString();
	}

	// ========================================================================
	private static TextTitle getChromSubtitle() {
		TextTitle chromSubtitle = new TextTitle(ChromTitleString);
		chromSubtitle.setHorizontalAlignment(HorizontalAlignment.LEFT);
		chromSubtitle.setFont(new Font("Arial", Font.BOLD, 24));
		chromSubtitle.setMargin(10, 190, 0, 0);
		return chromSubtitle;
	}
	
	// ========================================================================
	private static void formatXYPlotBackground(XYPlot xyPlot, double maxYValue, boolean doWhiteBackground, double backgroundBufferYTop, double chromLabelHeightBuffer) {
		Color backgroundGreyArmP = new Color(235, 235, 235, 220);
	
		backgroundBufferYTop   += maxYValue;
		chromLabelHeightBuffer += maxYValue; 
		
		for (Chrom chrom : Chrom.Autosomes) {			
			double posStartAdj = chrom.getGenomeWidePositionStart()             / GenomeWidePositionDivisor;
			double posEndAdj   = chrom.calculateGenomeWidePositionEnd()         / GenomeWidePositionDivisor;
			double posEndPArmAdj = chrom.calculateGenomeWidePositionEndOfArmP() / GenomeWidePositionDivisor;
			double posMidAdj   = chrom.calculateGenomeWidePositionMidpoint()    / GenomeWidePositionDivisor;
			
			if (doWhiteBackground) {
				xyPlot.addAnnotation(new XYBoxAnnotation(posStartAdj, maxYValue, posEndAdj, backgroundBufferYTop, null, null, ColorPastel.White.getColor()));
			}
			xyPlot.addDomainMarker(new IntervalMarker(posStartAdj, posEndPArmAdj, backgroundGreyArmP), Layer.BACKGROUND);
			xyPlot.addDomainMarker(new ValueMarker(posStartAdj, ColorPastel.Gray_90.getColor(), new BasicStroke()));
			XYTextAnnotation chromNumber = new XYTextAnnotation("" + chrom.ordinal(), posMidAdj, chromLabelHeightBuffer);
			chromNumber.setFont(new Font("Arial", Font.BOLD, 24));
			xyPlot.addAnnotation(chromNumber);
		}
	}

	// ========================================================================
	private static void formatLegendItems(XYPlot xyPlot) {
		LegendItemCollection legendItemCollection = xyPlot.getLegendItems();
		int numItems = legendItemCollection.getItemCount();
		for (int i = 0; i < numItems; i++) {
			LegendItem lItem = legendItemCollection.get(i);
			lItem.setShape(new Rectangle(25, 25));			
		}
		xyPlot.setFixedLegendItems(legendItemCollection);
	}
	
	// ========================================================================
	/** Plots the VAF (variant allele frequency) of a tissue vs the genomic position. */
	public static void plotVAFGenomeWide(XYDataset xyDataset, String outFilenameRoot, String sampleName, boolean isTumor) {
		String xAxisLabel = GenomeWideXAxisLabel;
		String yAxisLabel = "Variant Allele Fraction (VAF)";
		String title = sampleName.substring(0, sampleName.indexOf(".")) + (isTumor ? " [Tumor]" : " [Normal]");
		
		//+ "\n" + ChromTitleString;
		
		boolean showLegend = true;
		JFreeChart theChart = 
				ChartFactory.createScatterPlot(title, xAxisLabel, yAxisLabel, xyDataset, PlotOrientation.VERTICAL, showLegend, false, false);
		XYPlot xyPlot = (XYPlot) theChart.getPlot();
		//theChart.addSubtitle(getChromSubtitle());
		theChart.getTitle().setMargin(0, 0, 10, 0);
		theChart.getTitle().setFont(new Font("Arial", Font.BOLD, 32));
		
		XYItemRenderer itemRenderer = ClusteringPlotting.getXYItemRendererHelper(CopyNumberVAFPlotPointSize);		
		setSeriesPaintPerCluster(itemRenderer);
		xyPlot.setRenderer(itemRenderer);		
		
		xyPlot.setBackgroundPaint(ColorPastel.White.getColor());
		xyPlot.setDomainGridlinePaint(Color.white);
		xyPlot.setRangeGridlinePaint(Color.white);
				
		double maxValueY = 1.0;
		formatXYPlotBackground(xyPlot, maxValueY, true, 0.04, 0.02);
		
		formatLegendItems(xyPlot);
		
		if (showLegend) {
			LegendTitle legendTitle = theChart.getLegend();
			legendTitle.setID("Clusters");
			legendTitle.setItemFont(new Font("Arial", Font.BOLD, 32));
		}
	
		xyPlot.getRangeAxis().setRange(0, 1.04);		
		xyPlot.getDomainAxis().setRange(0, Chrom.c22.calculateGenomeWidePositionEnd() / GenomeWidePositionDivisor);
		
		Font rangeAxisLabelFont = new Font("Arial", Font.BOLD, 32);
		Font rangeAxisTickFont = new Font("Arial", Font.BOLD, 20);

		//xyPlot.getDomainAxis().setTickLabelsVisible(false);
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
		
		ControlFlagBool allGray                = new ControlFlagBool(false);
		ControlFlagBool makeNonEventsInvisible = new ControlFlagBool(false);
			
		if (allGray.getValue()) {
			for (EventType eventType : EventType.values()) {
				itemRenderer.setSeriesPaint(eventType.ordinal(), ColorPastel.Gray_50.getColor());
			}
			if (makeNonEventsInvisible.getValue()) {
				itemRenderer.setSeriesPaint(EventType.HETGermline.ordinal(), ColorPastel.Gray_15.getColor());
				itemRenderer.setSeriesPaint(EventType.HETSomatic.ordinal(),  ColorPastel.Gray_15.getColor());
				itemRenderer.setSeriesPaint(EventType.Noise.ordinal(),       ColorPastel.Gray_15.getColor());
				itemRenderer.setSeriesPaint(EventType.Ignored.ordinal(),     ColorPastel.Gray_15.getColor());
			}
		} else {
			for (EventType eventType : EventType.values()) {
				itemRenderer.setSeriesPaint(eventType.ordinal(), getColorForEvent(eventType).getColor());
			}
		}
		//itemRenderer.setSeriesPaint(EventType.HETSomatic.ordinal(), getColorForEvent(EventType.HETSomatic).getColor());			
	}

	// ========================================================================
	public static ColorPastel getColorForEvent(EventType eventType) {
		System.out.println("EV:\t" + eventType);
		switch(eventType) {
		case GainGermline: return ColorPastel.Violet;
		case GainSomatic:  return ColorPastel.Dark_Red;
		case LOH:          return ColorPastel.RGB_Blue;
		case cnLOH:        return ColorPastel.Dark_Yellow;
		case HETGermline:  return ColorPastel.Gray_60;
		case HETSomatic:   return ColorPastel.Red_Orange;
		case DELHom:       return ColorPastel.Light_Cyan_Blue;
		case Noise:        return ColorPastel.Dark_Pea_Green;
		case Ignored:      return ColorPastel.Gray_30;			
		default:           return null;
		}
	}
}
