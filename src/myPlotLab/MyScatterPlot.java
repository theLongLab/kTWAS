package myPlotLab;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.geom.Area;
import java.awt.geom.Ellipse2D;
import java.io.File;

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

public class MyScatterPlot extends ApplicationFrame{

	
	 public MyScatterPlot(final String title, String x_lab, String y_lab, String[] legend, 
	    		double[] x_values, double[][] y_values, int width, int hight, String output_file) {
	       super(title);
	        
	       final XYDataset dataset = createDataset(legend, x_values, y_values);
	       final JFreeChart chart = createChart(dataset, title,x_lab,y_lab);
	       final ChartPanel chartPanel = new ChartPanel(chart);
	       chartPanel.setPreferredSize(new java.awt.Dimension(500, 270));
	       setContentPane(chartPanel);
	        
	       try {
				final ChartRenderingInfo chartRenderingInfo = new ChartRenderingInfo(new StandardEntityCollection());	
				ChartUtilities.saveChartAsPNG(new File(output_file), chart, width, hight, chartRenderingInfo);
		} catch (Exception e) {e.printStackTrace(); }
	}
	
	 private XYDataset createDataset(String[] legend, double[] x_values, double[][] y_values) {
	    	if(legend.length!=y_values.length){
	    		System.out.println("legend.length!=y_values.length");
	    		return null;
	    	}
	    	final XYSeriesCollection dataset = new XYSeriesCollection(); 
	    	for(int i=0;i<legend.length;i++){
	    		 final XYSeries series = new XYSeries(legend[i]);
	    		 for(int k=0;k<x_values.length;k++){
	    			 series.add(x_values[k],y_values[i][k]);
	    		 }
	    		 dataset.addSeries(series);
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
	            true,                     // include legend
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
	        	renderer.setSeriesShape(k, (new Ellipse2D.Double(5.0,5.0,5.0,5.0)));
	        }
//	        renderer.setSeriesShapesFilled(0, true);
	        plot.setRenderer(renderer);

	        // change the auto tick unit selection to integer units only...
	        final NumberAxis rangeAxis = (NumberAxis) plot.getRangeAxis();
	        rangeAxis.setStandardTickUnits(NumberAxis.createIntegerTickUnits());
	        // OPTIONAL CUSTOMISATION COMPLETED.
	                
	        return chart;
	        
	    }
	 
	/** * Starting point for the demonstration application. * * @param args ignored. */
	public static void main(String[] args) { 
		String file="/Users/quan.long/areachart4.png";
    	String[] legend={"Z","Y","X"};
    	double[] x_values=new double[10];
    	double[][] y_values=new double[3][10];
    	for(int k=0;k<x_values.length;k++){
    		x_values[k]=k+1;
    		y_values[0][k]=(double)(k-10)*5.0;
    		y_values[1][k]=(double)k*k/5.0;
    		y_values[2][k]=10.0*Math.sqrt((double)k);
    	}
    	y_values[1][2]=Double.NaN;
        final MyScatterPlot demo = new MyScatterPlot(null, "my_x_lab", "my_y_lab", legend, x_values, y_values, 600, 400, file);
	}

}
