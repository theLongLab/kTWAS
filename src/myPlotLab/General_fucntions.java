package myPlotLab;

import java.io.BufferedWriter;
import java.io.FileWriter;

public class General_fucntions {

	public static void output_source(String output_file, String[] legend, 
    		double[] x_values, double[][] y_values){
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(output_file));
			if(y_values.length!=legend.length){
				System.out.println("y_values.length!=legend.length. Return");return;
			}if(y_values[0].length!=x_values.length){
				System.out.println("y_values[0].length!=x_values.length. Return");return;
			}
			bw.write("X_values");
			for(int i=0;i<legend.length;i++){
				bw.write(","+legend[i]);
			}bw.write("\n");
			for(int k=0;k<x_values.length;k++){
				bw.write(x_values[k]+"");
				for(int i=0;i<legend.length;i++){
					bw.write(","+y_values[i][k]);
				}bw.write("\n");
			}bw.close();			
		}catch(Exception e){e.printStackTrace();}
		
	}

}
