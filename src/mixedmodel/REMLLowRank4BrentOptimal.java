package mixedmodel;

import org.apache.commons.math3.analysis.UnivariateFunction;

public class REMLLowRank4BrentOptimal implements UnivariateFunction {
	
	public static double sigma_g2_REMLL(double[][] X_before_transformation, double[] transformed_Y, double[][] transformed_X, 
			double[] transformed2_Y, double[][] transformed2_X, double[] S, double delta, double[] beta){
		double n=X_before_transformation[0].length; // sample size
		double nk=transformed_Y.length; // smaller "sample size" after transformation.
		int num_coef=transformed_X.length;
		double sum1=0;
		for(int i=0;i<nk;i++){
			double diff_i=transformed_Y[i];
			for(int k=0;k<num_coef;k++){
				diff_i=diff_i-beta[k]*transformed_X[k][i];
			}
			sum1=sum1+(diff_i*diff_i)/(S[i]+delta);
		}
		double sum2=0;
		for(int i=0;i<n;i++){
			double diff_i=transformed2_Y[i];
			for(int k=0;k<num_coef;k++){
				diff_i=diff_i-beta[k]*transformed2_X[k][i];
			}
			sum1=sum1+(diff_i*diff_i);
		}
		return ((sum1+sum2/delta)/(n-num_coef));
	}
	
	public double value(double delta){
		return -1;
	}
	
}
