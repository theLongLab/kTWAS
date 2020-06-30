package myPlotLab;

import java.awt.Color;
import java.awt.Font;
import java.awt.image.BufferedImage;
import java.io.File;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.ChartRenderingInfo;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.annotations.XYDataImageAnnotation;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.entity.StandardEntityCollection;
import org.jfree.chart.plot.Marker;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.ValueMarker;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.GrayPaintScale;
//import org.jfree.chart.renderer.LookupPaintScale;
import org.jfree.chart.renderer.PaintScale;
import org.jfree.data.general.DefaultHeatMapDataset;
import org.jfree.data.general.HeatMapDataset;
import org.jfree.data.general.HeatMapUtilities;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.ApplicationFrame;
import org.jfree.ui.Layer;
import org.jfree.ui.RectangleAnchor;
import org.jfree.ui.TextAnchor;

public class MyHeatMap extends ApplicationFrame{

	private HeatMapDataset dataset;
	
	
	public MyHeatMap(double[][] data, String title, String x_lab, String y_lab,
			int max_x, int max_y, String output_file){
		super(title);
		final JFreeChart chart = createChart(new XYSeriesCollection(), data, title, x_lab, y_lab,
				 max_x, max_y);
		final ChartPanel chartPanel = new ChartPanel(chart);
        chartPanel.setPreferredSize(new java.awt.Dimension(500, 270));
        setContentPane(chartPanel);        
        try {
			final ChartRenderingInfo chartRenderingInfo = new ChartRenderingInfo(new StandardEntityCollection());	
			ChartUtilities.saveChartAsPNG(new File(output_file), chart, 2000, 400, chartRenderingInfo);
		} catch (Exception e) {e.printStackTrace(); }
	}
	
	/*
	 * create chart
	 */
	private JFreeChart createChart(XYDataset dataset, double[][] data, String title, String x_lab, String y_lab,
			int max_x, int max_y) {

        JFreeChart chart = ChartFactory.createScatterPlot(title, x_lab,y_lab, dataset, PlotOrientation.VERTICAL,
        		true, false, false);
        this.dataset = createMapDataset(data);
        PaintScale ps = new GrayPaintScale(0.8, 1, 255);
//        PaintScale ps =new LookupPaintScale(0.0, 100.0, Color.RED);
//        for(int k=0;k<40;k++){
//        	((LookupPaintScale)ps).add((Number)(k), (new Color(k, 0, 0)));
//        }
//        ((LookupPaintScale)ps).add(0, Color.blue);
//        ((LookupPaintScale)ps).add((Number)0.1, Color.cyan);
//        ((LookupPaintScale)ps).add((Number)0.2, Color.red);
//        ((LookupPaintScale)ps).add((Number)0.3, Color.cyan);
        BufferedImage image = HeatMapUtilities.createHeatMapImage(this.dataset,ps);
       
        XYDataImageAnnotation ann = new XYDataImageAnnotation(image, 0,
                0, max_x, max_y, true);        
        XYPlot plot = (XYPlot) chart.getPlot();
        plot.setBackgroundPaint(Color.white);
        plot.setDomainGridlinePaint(Color.white);
        plot.setDomainPannable(true);
        plot.setRangePannable(true);
        plot.getRenderer().addAnnotation(ann, Layer.BACKGROUND);
        
        double[] ends={0,30427671,19698289,23459830,18585056};
        double current_bk=0;
        for(int bk=0;bk<ends.length;bk++){
        	current_bk+=ends[bk]/1000000;
        	Marker currentEnd = new ValueMarker(current_bk);
        	Color color=MyManhattanPlot.my_faverate_colors[bk % MyManhattanPlot.my_faverate_colors.length];
        	currentEnd.setPaint(color);
            currentEnd.setLabelFont(new Font("Courier", Font.BOLD, 10));
            currentEnd.setLabel("Chr"+(bk+1));
            currentEnd.setLabelPaint(color);
            currentEnd.setLabelAnchor(RectangleAnchor.TOP_RIGHT);
            currentEnd.setLabelTextAnchor(TextAnchor.TOP_LEFT);
            plot.addDomainMarker(currentEnd);
        }
//        Marker currentEnd = new ValueMarker(28);
//        currentEnd.setPaint(Color.red);
//        currentEnd.setLabelFont(new Font("Courier", Font.BOLD, 10));
//        currentEnd.setLabel("Chr2");
//        currentEnd.setLabelAnchor(RectangleAnchor.TOP_RIGHT);
//        currentEnd.setLabelTextAnchor(TextAnchor.TOP_LEFT);
//        plot.addDomainMarker(currentEnd);
        
        NumberAxis xAxis = (NumberAxis) plot.getDomainAxis();
        xAxis.setStandardTickUnits(NumberAxis.createIntegerTickUnits());
        xAxis.setLowerMargin(0.0);
        xAxis.setUpperMargin(0.0);
        NumberAxis yAxis = (NumberAxis) plot.getRangeAxis();
        yAxis.setStandardTickUnits(NumberAxis.createIntegerTickUnits());
        yAxis.setLowerMargin(0.0);
        yAxis.setUpperMargin(0.0);
        return chart;
    }
	
    private static HeatMapDataset createMapDataset(double[][] data) {
	    DefaultHeatMapDataset d = new DefaultHeatMapDataset(data.length, data[0].length, 1,
	    		data.length, 1, data[0].length);
	    for (int i = 0; i < data.length; i++) {
	        for (int j = 0; j < data[0].length; j++) {
	            d.setZValue(i, j, 1.0-data[i][j]);
	        }
	    }
	    return d;
	}
	
//    public static void main(String[] args) {
//		String output_file="/Users/quan.long/areachart6.png";
//		double data[][] =new double[2000][100];
//		for (int i = 0; i < data.length; i++) {
//            for (int j = 0; j < data[0].length; j++) {
//                data[i][j]=Math.abs(Math.sin(Math.sqrt(i * j) / 10));
//            }
//        }
//		MyHeatMap map=new MyHeatMap(data, "title", "x_lab","y_lab", 120, 100, output_file);
//
//	}

}
