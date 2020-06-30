package myPlotLab;


import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Shape;
import java.awt.geom.Area;
import java.awt.geom.Ellipse2D;
import java.io.File;
import java.util.HashMap;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.ChartRenderingInfo;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.entity.StandardEntityCollection;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.category.DefaultCategoryDataset;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.ApplicationFrame;
import org.jfree.util.ShapeList;

public class MyManhattanPlot extends ApplicationFrame{

	 public static Color[] my_faverate_colors={Color.RED,Color.BLUE,Color.ORANGE,Color.GREEN,Color.MAGENTA};	
	 public static int[] lengths_of_chrs={30427671,19698289,23459830,18585056,26975502}; //TODO
	 
	 public MyManhattanPlot(final String title, String x_lab, String y_lab, int[] series_indicator, 
	    		double[] x_values, double[] y_values, int width, int hight, String output_file) {
		super(title);
	    
		myMathLib.ArrayFuncs.remove_infinites(y_values);
	    final XYDataset dataset = createDataset(series_indicator, x_values, y_values);
	    final JFreeChart chart = createChart(dataset, title,x_lab,y_lab);
	    final ChartPanel chartPanel = new ChartPanel(chart);
	    chartPanel.setPreferredSize(new java.awt.Dimension(1000, 270));
	    setContentPane(chartPanel);
	        
	    try {
			final ChartRenderingInfo chartRenderingInfo = new ChartRenderingInfo(new StandardEntityCollection());	
			ChartUtilities.saveChartAsPNG(new File(output_file), chart, width, hight, chartRenderingInfo);
		}catch (Exception e) {e.printStackTrace(); }
	}
	
	 private XYDataset createDataset(int[] series_indicator, double[] x_values, double[] y_values) {
		if(x_values.length!=y_values.length || x_values.length!=series_indicator.length){
	    		System.out.println("x_values.length!=y_values.length || x_values.length!=series_indicator.length");
	    		return null;
	   	}
	  	final XYSeriesCollection dataset = new XYSeriesCollection();
	   	HashMap<Integer, XYSeries> the_series= new HashMap<>();
	   	for(int i=0;i<x_values.length;i++){
	   		if(the_series.containsKey(series_indicator[i])){
	  			XYSeries series=the_series.get(series_indicator[i]);
	   			series.add(x_values[i],y_values[i]);
	  		}else{
	  			XYSeries series = new XYSeries("S"+series_indicator[i]);
	   			series.add(x_values[i],y_values[i]);
	   			the_series.put(series_indicator[i], series);
	   			dataset.addSeries(series);
	 		}
	    }         
	   	return dataset;        
	 }
	 
	 
	 public MyManhattanPlot(final String title, String x_lab, String y_lab,  
	    		double[][] x_values, double[][] y_values, int width, int hight, String output_file) {
	    super(title);
	        
	    final XYDataset dataset = createDataset(x_values, y_values);
	    final JFreeChart chart = createChart(dataset, title,x_lab,y_lab);
	    final ChartPanel chartPanel = new ChartPanel(chart);
	    chartPanel.setPreferredSize(new java.awt.Dimension(1000, 270));
	    setContentPane(chartPanel);
	        
	    try {
			final ChartRenderingInfo chartRenderingInfo = new ChartRenderingInfo(new StandardEntityCollection());	
			ChartUtilities.saveChartAsPNG(new File(output_file), chart, width, hight, chartRenderingInfo);
		}catch (Exception e) {e.printStackTrace(); }
	}
	 private XYDataset createDataset(double[][] x_values, double[][] y_values) {
		 	if(x_values.length!=y_values.length ){
		    		System.out.println("x_values.length!=y_values.length");
		    		return null;
		   	}
		 	double current_x=0;
		  	final XYSeriesCollection dataset = new XYSeriesCollection();
		   	for(int i=0;i<x_values.length;i++){
		   		if(x_values[i].length!=y_values[i].length ){
		    		System.out.println("x_values[i].length!=y_values[i].length");		    		
		   		}
		   		XYSeries series = new XYSeries("S"+i);
		   		for(int j=0;j<x_values[i].length;j++){
		   			series.add(current_x+x_values[i][j], y_values[i][j]);		   			
		   		}dataset.addSeries(series);
		   		current_x+=x_values[i][x_values[i].length-1];
		    }         
		   	return dataset;        
		 }
	 
	 private JFreeChart createChart(final XYDataset dataset, String title, String x_lab, String y_lab) {
	        // create the chart...
	        final JFreeChart chart = ChartFactory.createScatterPlot(
	        	title,      // chart title
	        	x_lab,                      // x axis label
	        	y_lab,                      // y axis label
	            dataset,                  // data
	            PlotOrientation.VERTICAL,
	            false,                     // include legend
	            true,                     // tooltips
	            false                     // urls
	        );

	        // NOW DO SOME OPTIONAL CUSTOMISATION OF THE CHART...
	        chart.setBackgroundPaint(Color.white);

//	        final StandardLegend legend = (StandardLegend) chart.getLegend();
	  //      legend.setDisplaySeriesShapes(true);
	        
	        // get a reference to the plot for further customisation...
	        final XYPlot plot = chart.getXYPlot();
	        plot.setBackgroundPaint(Color.white);
	    //    plot.setAxisOffset(new Spacer(Spacer.ABSOLUTE, 5.0, 5.0, 5.0, 5.0));
	        plot.setDomainGridlinePaint(Color.white);
	        plot.setRangeGridlinePaint(Color.lightGray);
	        
	        final XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();
	        for(int k=0;k<dataset.getSeriesCount();k++){
	        	renderer.setSeriesLinesVisible(k, false);
	        	Shape the_shape=(new Ellipse2D.Double(5.0,5.0,5.0,5.0));
	        	renderer.setSeriesShape(k, the_shape);	        	
	        	renderer.setSeriesPaint(k, my_faverate_colors[k%my_faverate_colors.length] );
	        }
//	        renderer.setSeriesShapesFilled(0, true);
	        plot.setRenderer(renderer);

	        // change the auto tick unit selection to integer units only...
	        final NumberAxis rangeAxis = (NumberAxis) plot.getRangeAxis();
	        rangeAxis.setStandardTickUnits(NumberAxis.createIntegerTickUnits());
	        rangeAxis.setLabelFont(new Font("Arial", Font.BOLD, 18));
	        // OPTIONAL CUSTOMISATION COMPLETED.
	        return chart;
	        
	    }
	 
	/** * Starting point for the demonstration application. * * @param args ignored. */
	public static void main(String[] args) { 
		String file="/Users/quan.long/areachart7.png";
    	double[][] x_values={{1,2,3},{1,2},{1,2,3,4},{1,2},{1,2,3,4,5}};
    	double[][] y_values={{1,1,1},{2,2},{3,2,3,4},{4,2},{1,2,3,4,5}};
    	int[] indicators={0,0,0,1,1,1,2,2,3,3};
        final MyManhattanPlot demo = new MyManhattanPlot(null, "my_x_lab", "my_y_lab", x_values, y_values, 2000, 400, file);
	}

}
