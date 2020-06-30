package myPlotLab;

import java.awt.Color;
import java.io.File;
import java.util.Arrays;

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
import org.jfree.data.xy.XYDataItem;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.ApplicationFrame;
import org.jfree.ui.RefineryUtilities;
//import org.jfree.ui.Spacer;

/**
 * A simple demonstration application showing how to create a line chart using data from an
 * {@link XYDataset}.
 *
 */
public class MyLineChart extends ApplicationFrame {

    /**
     * Creates a new demo.
     *
     * @param title  the frame title.
     */
    public MyLineChart(final String title, String x_lab, String y_lab, String[] legend, 
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
			General_fucntions.output_source(output_file+".source.csv", legend, x_values, y_values);
		} catch (Exception e) {e.printStackTrace(); }
        
    }
    
    public MyLineChart(final String title, String x_lab, String y_lab, String[] legend, 
    		double[][] xy_values, int width, int hight, String output_file, String nornalized) {
        super(title);
        if(nornalized.equals("median")){
        	for(int i=1;i<xy_values.length;i++){
        		double[] array_to_sort=xy_values[i].clone();
        		Arrays.sort(array_to_sort);
        		double median=array_to_sort[array_to_sort.length/2];
        		for(int k=0;k<xy_values[i].length;k++){
        			xy_values[i][k]/=median;
        		}
        	}
        }else if(nornalized.equals("mean")){
        	for(int i=1;i<xy_values.length;i++){
        		double sum=0;
        		for(int k=0;k<xy_values[i].length;k++){
        			sum+=xy_values[i][k];
        		}
        		double mean=sum/xy_values[i].length;
        		for(int k=0;k<xy_values[i].length;k++){
        			xy_values[i][k]/=mean;
        		}
        	}
        }
        final XYDataset dataset = createDataset(legend, xy_values);
        final JFreeChart chart = createChart(dataset, title,x_lab,y_lab);
        final ChartPanel chartPanel = new ChartPanel(chart);
        chartPanel.setPreferredSize(new java.awt.Dimension(500, 270));
        setContentPane(chartPanel);
        
        try {
			final ChartRenderingInfo chartRenderingInfo = new ChartRenderingInfo(new StandardEntityCollection());	
			ChartUtilities.saveChartAsPNG(new File(output_file), chart, width, hight, chartRenderingInfo);
		} catch (Exception e) {e.printStackTrace(); }
        return;
    }
    
    /**
     * Creates a sample dataset.
     * 
     * @return a sample dataset.
     */
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
    
    private XYDataset createDataset(String[] legend, double[][] xy_values) {
    	if(legend.length!=xy_values.length-1){
    		System.out.println("legend.length!=xy_values.length-1");
    		return null;
    	}
    	final XYSeriesCollection dataset = new XYSeriesCollection();
    	for(int i=0;i<legend.length;i++){
    		final XYSeries series = new XYSeries(legend[i]);
    		for(int k=0;k<xy_values[i+1].length;k++){
    			series.add(xy_values[0][k],xy_values[i+1][k]);
    		}
    		dataset.addSeries(series);
    	}         
    	return dataset;        
    }
    
    /**
     * Creates a chart.
     * 
     * @param dataset  the data for the chart.
     * 
     * @return a chart.
     */
    private JFreeChart createChart(final XYDataset dataset, String title, String x_lab, String y_lab) {
        
        // create the chart...
        final JFreeChart chart = ChartFactory.createXYLineChart(
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

//        final StandardLegend legend = (StandardLegend) chart.getLegend();
  //      legend.setDisplaySeriesShapes(true);
        
        // get a reference to the plot for further customisation...
        final XYPlot plot = chart.getXYPlot();
        plot.setBackgroundPaint(Color.white);
    //    plot.setAxisOffset(new Spacer(Spacer.ABSOLUTE, 5.0, 5.0, 5.0, 5.0));
        plot.setDomainGridlinePaint(Color.white);
        plot.setRangeGridlinePaint(Color.lightGray);
        
        final XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();
//        renderer.setSeriesLinesVisible(0, false);
//        renderer.setSeriesShapesVisible(1, false);
//        plot.setRenderer(renderer);

        // change the auto tick unit selection to integer units only...
        final NumberAxis rangeAxis = (NumberAxis) plot.getRangeAxis();
        rangeAxis.setStandardTickUnits(NumberAxis.createIntegerTickUnits());
        // OPTIONAL CUSTOMISATION COMPLETED.
                
        return chart;
        
    }

    // ****************************************************************************
    // * JFREECHART DEVELOPER GUIDE                                               *
    // * The JFreeChart Developer Guide, written by David Gilbert, is available   *
    // * to purchase from Object Refinery Limited:                                *
    // *                                                                          *
    // * http://www.object-refinery.com/jfreechart/guide.html                     *
    // *                                                                          *
    // * Sales are used to provide funding for the JFreeChart project - please    * 
    // * support us so that we can continue developing free software.             *
    // ****************************************************************************
    
    /**
     * Starting point for the demonstration application.
     *
     * @param args  ignored.
     */
    public static void main(final String[] args) {
    	String file="/Users/quanlong/areachart2.png";
    	String[] legend={"Z","Y","X"};
    	double[] x_values=new double[10000];
    	double[][] y_values=new double[3][10000];
    	for(int k=0;k<x_values.length;k++){
    		x_values[k]=k+1;
    		y_values[0][k]=(double)k-10;
    		y_values[1][k]=(double)k*k/10000.0;
    		y_values[2][k]=100.0*Math.sqrt((double)k);
    	}
        final MyLineChart demo = new MyLineChart(null, "my_x_lab", "my_y_lab", legend, x_values, y_values, 600, 400, file);
//       demo.pack();
//        RefineryUtilities.centerFrameOnScreen(demo);
//        demo.setVisible(true);

    }

}
	

	