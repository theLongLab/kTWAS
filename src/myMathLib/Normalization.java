package myMathLib;

import java.util.ArrayList;
import java.util.Collections;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.stat.StatUtils;

public class Normalization {
	
	public static void normalize_colomns(double[][] data){
		double sums[]=new double[data[0].length];
		for(int k=0;k<data.length;k++){
			for(int i=0;i<data[0].length;i++){
				sums[i]+=data[k][i];
			}
		}
		for(int k=0;k<data.length;k++){
			for(int i=0;i<data[0].length;i++){
				data[k][i]/=sums[i];
			}
		}
	}
	
	public static double median(ArrayList<Double> data){
		Collections.sort(data);
		int index=data.size()/2;
		if(data.size()%2==1){
			return data.get(index);
		}else{
			return (data.get(index)+data.get(index-1))/2.0;
		}
	}
	
	public static double[][] byte2double(byte[][] input){
		if(input==null)return null;
		double[][] out=new double[input.length][];
		for(int i=0;i<input.length;i++){
			out[i]=new double[input[i].length];
			for(int j=0;j<input[i].length;j++){
				out[i][j]=input[i][j];
			}
		}
		return out;
	}
	
	public static double[] byte2double(byte[] input){
		if(input==null)return null;
		double[] out=new double[input.length];
		for(int i=0;i<input.length;i++){			
			out[i]=input[i];			
		}
		return out;
	}
	
	public static void normalize_variance(double[] data){
		double sd=Math.sqrt(StatUtils.populationVariance(data));
		for(int k=0;k<data.length;k++)data[k]=data[k]/sd;
	}
	
	/*
	 *  Given a vector of data (usually simulated phenoytpe using genomes), 
	 *  add random component so that the genetic component is the value pre-specified.  
	 * 
	 *  \sigma_g/(\sigma_g + \sigma_e) = h^2 (=genetic_component);
	 */
	public static void adding_random_component(double[] data, double genetic_component) {
		double g_variance=StatFuncs.popVar_NaN(data);
		double e_variance= g_variance*(1-genetic_component)/genetic_component;
		NormalDistribution normal=new NormalDistribution(0, Math.sqrt(e_variance));
		for(int k=0;k<data.length;k++) {
			data[k]=data[k]+normal.sample();
		}
	}
}
