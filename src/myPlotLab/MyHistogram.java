package myPlotLab;


import java.awt.Dimension;
import java.io.File;


import javax.swing.JPanel;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.ChartRenderingInfo;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.entity.StandardEntityCollection;
import org.jfree.chart.labels.StandardCategoryToolTipGenerator;
import org.jfree.chart.plot.CategoryPlot;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.renderer.category.BarRenderer;
import org.jfree.chart.title.TextTitle;
import org.jfree.data.category.CategoryDataset;
import org.jfree.data.category.DefaultCategoryDataset;
import org.jfree.ui.ApplicationFrame;
import org.jfree.ui.RectangleEdge;
import org.jfree.ui.RefineryUtilities;

public class MyHistogram extends ApplicationFrame {
	
	/*
	 *  multiple series
	 */
	public MyHistogram(double[][] data, String[] categories, String[] series, String title, String x_lab, String y_lab, 
			int width, int hight, String output_file) { 
		super(title);
		if(categories.length!=data.length){
			System.out.println("categories.length!=data.length");
			return;
		}
		DefaultCategoryDataset dataset = new DefaultCategoryDataset(); 		
		for(int k=0;k<categories.length;k++){
			for(int i=0;i<data[k].length;i++)
				dataset.addValue(data[k][i], series[i], categories[k]);
		}
		
		JFreeChart chart = ChartFactory.createBarChart(
				title, x_lab, y_lab, 
				dataset, 
				PlotOrientation.VERTICAL, 
				true, // include legend
				true, // tooltips?
				false // URLs? 
				);	
		ChartPanel chartPanel = new ChartPanel(chart, false); 
		chartPanel.setPreferredSize(new Dimension(500, 270)); 
		setContentPane(chartPanel);	
		CategoryPlot plot = (CategoryPlot) chart.getPlot();
		BarRenderer renderer = (BarRenderer) plot.getRenderer();
		renderer.setItemMargin(0.0);
		
		 try {
			final ChartRenderingInfo chartRenderingInfo = new ChartRenderingInfo(new StandardEntityCollection());	
			ChartUtilities.saveChartAsPNG(new File(output_file), chart, width, hight, chartRenderingInfo);
			
		} catch (Exception e) {e.printStackTrace(); }
	}
	
	/*
	 *  multiple series, stacked
	 */
	public MyHistogram(double[][] data, String[] categories, String[] series, String title, String x_lab, String y_lab, 
			int width, int hight, String output_file, boolean stack) { 
		super(title);
		if(categories.length!=data.length){
			System.out.println("categories.length!=data.length");
			return;
		}
		DefaultCategoryDataset dataset = new DefaultCategoryDataset(); 		
		for(int k=0;k<categories.length;k++){
			for(int i=0;i<data[k].length;i++)
				dataset.addValue(data[k][i], series[i], categories[k]);
		}
		
		JFreeChart chart = ChartFactory.createStackedBarChart(
				title, x_lab, y_lab, 
				dataset, 
				PlotOrientation.VERTICAL, 
				true, // include legend
				true, // tooltips?
				false // URLs? 
				);	
		ChartPanel chartPanel = new ChartPanel(chart, false); 
		chartPanel.setPreferredSize(new Dimension(500, 270)); 
		setContentPane(chartPanel);	
		CategoryPlot plot = (CategoryPlot) chart.getPlot();
		BarRenderer renderer = (BarRenderer) plot.getRenderer();
		renderer.setItemMargin(0.0);
		
		 try {
			final ChartRenderingInfo chartRenderingInfo = new ChartRenderingInfo(new StandardEntityCollection());	
			ChartUtilities.saveChartAsPNG(new File(output_file), chart, width, hight, chartRenderingInfo);
		} catch (Exception e) {e.printStackTrace(); }
	}
	
	/*
	 * single series
	 */
	public MyHistogram(double[] data, String[] categories, String title, String x_lab, String y_lab, int width, int hight,
			String output_file) { 
		super(title);
		if(categories.length!=data.length){
			System.out.println("categories.length!=data.length");
			return;
		}
		DefaultCategoryDataset dataset = new DefaultCategoryDataset(); 		
		for(int k=0;k<categories.length;k++){
			dataset.addValue(data[k], x_lab, categories[k]);
		}
		JFreeChart chart = ChartFactory.createBarChart(
				title, x_lab, y_lab, 
				dataset, 
				PlotOrientation.VERTICAL, 
				false, // include legend
				true, // tooltips?
				false // URLs? 
				);	
		ChartPanel chartPanel = new ChartPanel(chart, false); 
		chartPanel.setPreferredSize(new Dimension(500, 270)); 
		setContentPane(chartPanel);		
		 try {
			final ChartRenderingInfo chartRenderingInfo = new ChartRenderingInfo(new StandardEntityCollection());	
				ChartUtilities.saveChartAsPNG(new File(output_file), chart, width, hight, chartRenderingInfo);
		}catch (Exception e) {e.printStackTrace(); }
	}
	
	/** * Starting point for the demonstration application. * * @param args ignored. */
	public static void main(String[] args) { 
		double[][] data={{1,2,5},{2,3,4},{10,9,8}};
		String[] series={"X1","X2","X3"};
		String[] category={"A","B","C"};
		MyHistogram demo = new MyHistogram(data, category, series, "Title", null,"Percentage %", 1000,400,"/Users/quan.long/areachart3.png", true);
//		MyHistogram demo = new MyHistogram(data[0], cate, "Title", "X","Y", 1000,400,"/Users/quan.long/areachart3.png");
//		demo.pack(); //RefineryUtilities.centerFrameOnScreen(demo); 
//		demo.setVisible(true);
	}

}
