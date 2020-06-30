package myMathLib;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import org.apache.commons.math3.distribution.FDistribution;
import org.apache.commons.math3.distribution.KolmogorovSmirnovDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.distribution.TDistribution;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;

public class StatFuncs {
	
	/*
	 *  Kolmogorov Smirnov test
	 */
	public static void ks_test(double[] x, double[] y){
		double[] xx=x.clone();
		double[] yy=y.clone();
		Arrays.sort(xx);
		Arrays.sort(yy);
	}
	
	/*
	 * normality_test using Kolmogorov�Smirnov test. 
	 */
	public static double normality_test_ks(double[] x){
		double[] data=x.clone(); // don't modify the original data x!
		Arrays.sort(data);
		double mean=StatUtils.mean(data);
		double var=StatUtils.variance(data,mean);
		NormalDistribution normal=new NormalDistribution(mean, Math.sqrt(var));
		double n=data.length;
		System.out.println("n="+n);
		double Dn=0;
		for(int i=0;i<n;i++){
			double Fx=normal.cumulativeProbability(data[i]);
			double Fnx=(i+1)/n;			
			if(Math.abs(Fnx-Fx)>Dn)Dn=Math.abs(Fnx-Fx);
		}
		//System.out.println("Dn done");
		KolmogorovSmirnovDistribution ks_dist=new KolmogorovSmirnovDistribution((int)n);		
		double pvalue=1-ks_dist.cdf(Dn);
		return pvalue;
	}
	
	public static boolean variable_or_not(double[] data){
		for(int k=1;k<data.length;k++){
			if(Double.compare(data[0],data[k])!=0)return true;
		}return false;
	}
	
	public static double[] multi_reg_pvalues(double[] b, double[] sde, int n){
		TDistribution t_dis=new TDistribution(n-b.length);
		double[] t=new double[b.length];
		for(int k=0;k<b.length;k++)t[k]=Math.abs(b[k]/sde[k]);
		double[] pvalues=new double[b.length];
		for(int k=0;k<b.length;k++)pvalues[k]=2*(1-t_dis.cumulativeProbability(t[k]));
		return pvalues;
	}
	
	/*
	 *  Calculate the F-test for linear regression
	 *  n=sample size; 
	 *  p=number of parameters
	 */
	public static double total_reg_pvalue_F(double[] Y, double[] residuals, int p){
		int n =Y.length;
		if(n!=residuals.length){
			System.out.println("F_Test Error: n!=residuals.length");
			return Double.NaN;
		}
		double Y_mean=mean_no_NaN(Y);
		double SSM=0, SSE=0;
		double[] Y_hat=new double[n];
		for(int k=0;k<n;k++){
			Y_hat[k]=Y[k]-residuals[k];
			SSM=SSM+(Y_hat[k]-Y_mean)*(Y_hat[k]-Y_mean);
			SSE=SSE+residuals[k]*residuals[k];
		}
		int dfm=p-1, dfe=n-p;
		double F=(SSM/dfm)/(SSE/dfe);
		FDistribution f_dis=new FDistribution(dfm, dfe);
		double pvalue=1-f_dis.cumulativeProbability(F);
//		System.out.println("SSM="+SSM);
//		System.out.println("SSE="+SSE);
		System.out.println("F="+F);
		return pvalue;
	}
	
	/*
	 *  Calculate the F-test for linear regression
	 *  n=sample size; 
	 *  p=number of parameters
	 */
	public static double total_reg_pvalue_F_nointercept(double[] Y, double[] residuals, int p){
		int n =Y.length;
		if(n!=residuals.length){
			System.out.println("F_Test Error: n!=residuals.length");
			return Double.NaN;
		}
		double Y_mean=mean_no_NaN(Y);
		double SSM=0, SSE=0;
		double[] Y_hat=new double[n];
		for(int k=0;k<n;k++){
			Y_hat[k]=Y[k]-residuals[k];
			SSM=SSM+(Y_hat[k])*(Y_hat[k]);
			SSE=SSE+residuals[k]*residuals[k];
		}
		int dfm=p-1, dfe=n-p+1;
		double F=(SSM/dfm)/(SSE/dfe);
		FDistribution f_dis=new FDistribution(dfm, dfe);
		double pvalue=1-f_dis.cumulativeProbability(F);
//		System.out.println("SSM="+SSM);
//		System.out.println("SSE="+SSE);
		System.out.println("F="+F);
		return pvalue;
	}
	
	public static double correlationPearsons(double[] v1, double[] v2){
		PearsonsCorrelation calculator=new PearsonsCorrelation();
		return calculator.correlation(v1, v2);
	}
	
	public static double correlationSpearmans(double[] v1, double[] v2){
		SpearmansCorrelation calculator=new SpearmansCorrelation();
		return calculator.correlation(v1, v2);
	}
	
	/*
	 * calculate the likelihood of null regression model with intercept only 
	 */
	public static double null_logL(double[] y){
		double mean=StatUtils.mean(y);
		double var=StatUtils.populationVariance(y, mean);
		NormalDistribution normal=new NormalDistribution(0, Math.sqrt(var));
		double log_likelihood=0;
		for(int i=0;i<y.length;i++){
			double prob=normal.density(y[i]-mean);
			log_likelihood=log_likelihood+Math.log(prob);
		}
//		double another=0;
//		double var2=StatUtils.variance(y, mean);
//		for(int i=0;i<y.length;i++){
//			another+=y[i]*y[i]/var2;
//		}
//		another+=y.length*Math.log(2*Math.PI*var2);
//		another=-0.5*another;
//		System.out.println(another+","+log_likelihood);
		return log_likelihood;
	}
	
	/*
	 * calculate the likelihood of null regression model with transformed intercept only 
	 */
	public static double null_logL(double[] y, double[] intercept){
		double mean_y=StatUtils.mean(y);
		double mean_intercept=StatUtils.mean(intercept);		
		double beta0=mean_y/mean_intercept;
		double[] new_y=myMathLib.ArrayFuncs.linearCombination(y, intercept, 1, -beta0); 
		double var=StatUtils.populationVariance(new_y);
		NormalDistribution normal=new NormalDistribution(0, Math.sqrt(var));
		double log_likelihood=0;
		for(int i=0;i<y.length;i++){
			double prob=normal.density(new_y[i]);
			log_likelihood=log_likelihood+Math.log(prob);
		}
		return log_likelihood;
	}
	
	/*
	 * Mean with NaN skipped
	 */
	public static double mean_NaN(double[] data){
		int total_num=0;
		double sum=0;
		for(int k=0;k<data.length;k++){
			if(!Double.isNaN(data[k])){
				sum+=data[k];
				total_num++;
			}
		}return sum/total_num;
	}
	
	public static double mean_no_NaN(double[] data){
		double sum=0;
		for(int k=0;k<data.length;k++){
			if(!Double.isNaN(data[k])){
				sum+=data[k];
			}else return Double.NaN;
		}return sum/data.length;
	}
	
	public static double mean(double[] data){
		double sum=0;
		for(int k=0;k<data.length;k++){			
			sum+=data[k];			
		}return sum/data.length;
	}
	
	public static double popVar_NaN(double[] data, double mean){
		if(Double.isNaN(mean))return Double.NaN;
		int total_num=0;
		double var_sum=0;
		for(int k=0;k<data.length;k++){
			if(!Double.isNaN(data[k])){
				var_sum+=(data[k]-mean)*(data[k]-mean);
				total_num++;
			}
		}return var_sum/total_num;
	}
	
	public static double popVar_NaN(double[] data){
		double mean=mean_NaN(data);
		return popVar_NaN(data,mean);
	}
	
	public static double[] standardize(double[] data){
		double[] out=new double[data.length];
		double mean=mean_NaN(data);
		double var=popVar_NaN(data,mean);
		for(int k=0;k<data.length;k++)out[k]=(data[k]-mean)/var;
		return out;
	}
	
	public static double var_NaN(double[] data, double mean){
		if(Double.isNaN(mean))return Double.NaN;
		int total_num=0;
		double var_sum=0;
		for(int k=0;k<data.length;k++){
			if(!Double.isNaN(data[k])){
				var_sum+=(data[k]-mean)*(data[k]-mean);
				total_num++;
			}
		}return var_sum/(total_num-1);
	}
	
	public static double var(double[] data, double mean){
		if(Double.isNaN(mean))return Double.NaN;
		double var_sum=0;
		for(int k=0;k<data.length;k++){			
			var_sum+=(data[k]-mean)*(data[k]-mean);			
		}return var_sum/(data.length-1);
	}
	
	public static double var_NaN(double[] data){
		double mean=mean_NaN(data);
		return var_NaN(data,mean);
	}
	
	public static double pvalue_spearman(double r, double n, TDistribution t_dist){
		double t=r*Math.sqrt((n-2)/(1-r*r));
		double p=1-t_dist.cumulativeProbability(Math.abs(t));
		return p;
	}
	
	public static double pvalue_spearman(double r, double n){
		TDistribution t_dist=new TDistribution(n-2);
		double t=r*Math.sqrt((n-2)/(1-r*r));
		double p=1-t_dist.cumulativeProbability(Math.abs(t));
		return p;
	}
	
	public static double pvalue_Ztranform(double r, double n){
		NormalDistribution n_dist=new NormalDistribution(0, 1.0/Math.sqrt(n-3));
		double z=0.5*Math.log((1+r)/(1-r));
		double p=1-n_dist.cumulativeProbability(Math.abs(z));
		return p;
	}
	
	public static double pvalue_Ztranform(double r, NormalDistribution n_dist){
		double z=0.5*Math.log((1+r)/(1-r));
		double p=1-n_dist.cumulativeProbability(Math.abs(z));
		return p;
	}
	
	/*
	 * ready for Spearman correlation
	 */
	public static void replace_ranks(double[] data){
		HashMap<Double, ArrayList<Integer>> indexes=new HashMap<Double, ArrayList<Integer>>();
		for(int i=0;i<data.length;i++){
			if(indexes.containsKey(data[i])){
				indexes.get(data[i]).add(i);
			}else{
				ArrayList<Integer> new_list=new ArrayList<Integer>();
				new_list.add(i);
				indexes.put(data[i], new_list);
			}			
		}
		double[] data2=data.clone();
		Arrays.sort(data2);
		int step=1;
		for(int i=0;i<data2.length;i=i+step){
			ArrayList<Integer> the_list=indexes.get(data2[i]);
			step=the_list.size();
			double the_rank=i+(step-1.0)/2.0;
			for(int k=0;k<step;k++)
				data[the_list.get(k)] = the_rank+1;
		}
	}
	
	/*
	 * dim(X)=#terms x sample_size
	 */
	public static double distance_correlation(double[][] X, double[][] Y){
		int n=X[0].length;
		if(Y[0].length!=n){System.out.println("sample size not equal");}
		int dx=X.length;
		int dy=Y.length;
		double[][] a=new double[n][n], b=new double[n][n];
		for(int j=0;j<n;j++){
			a[j][j]=0; b[j][j]=0;
			for(int k=j+1;k<n;k++){
				for(int tt=0;tt<dx;tt++){
					a[j][k]=a[j][k]+(X[tt][j]-X[tt][k])*(X[tt][j]-X[tt][k]);
				}a[j][k]=Math.sqrt(a[j][k]);
				a[k][j]=a[j][k];
				for(int tt=0;tt<dy;tt++){
					b[j][k]=b[j][k]+(Y[tt][j]-Y[tt][k])*(Y[tt][j]-Y[tt][k]);
				}b[j][k]=Math.sqrt(b[j][k]);
				b[k][j]=b[j][k];
			}
		}
		double[][] A=new double[n][n], B=new double[n][n];
		double[] a_col=new double[n], a_row=new double[n];
		double[] b_col=new double[n], b_row=new double[n];
		for(int j=0;j<n;j++){
			for(int k=0;k<n;k++){
				a_row[j]+=a[j][k]; a_col[k]+=a[j][k];
				b_row[j]+=b[j][k]; b_col[k]+=b[j][k];
			}
		}
		for(int j=0;j<n;j++){
			a_row[j]=a_row[j]/n;  b_row[j]=b_row[j]/n;
		}for(int k=0;k<n;k++){
			a_col[k]=a_col[k]/n;  b_col[k]=b_col[k]/n;
		}
		double a_grand=mean(a_row), b_grand=mean(b_row);
		if(Double.compare(a_grand, mean(a_col))!=0)System.out.print("WRONG");
		if(Double.compare(b_grand, mean(b_col))!=0)System.out.print("WRONG");
		for(int j=0;j<n;j++){
			for(int k=0;k<n;k++){
				A[j][k]=a[j][k]-a_row[j]-a_col[k]+a_grand;
				B[j][k]=b[j][k]-b_row[j]-b_col[k]+b_grand;
			}
		}
		double dCov2=0, dVarX2=0, dVarY2=0;
		for(int j=0;j<n;j++){
			for(int k=0;k<n;k++){
				dCov2+=A[j][k]*B[j][k];
				dVarX2+=A[j][k]*A[j][k];
				dVarY2+=B[j][k]*B[j][k];
			}
		}
		dCov2=dCov2/(n*n);  
		dVarX2=dVarX2/(n*n);
		dVarY2=dVarY2/(n*n);
		double dVarX=Math.sqrt(dVarX2);
		double dVarY=Math.sqrt(dVarY2);
		System.out.print(Math.sqrt(dCov2));
		return Math.sqrt(dCov2)/Math.sqrt(dVarX*dVarY);
	}
	
	public static void main(String[] args){
//		double[][] X={{1,1},{1,1},{1,1}};
//		double[][] Y={{1,4},{3,3},{5,7}};
//		
//		//distance_correlation(X, Y);
//		double r=0.24, n=100;
//		System.out.println(pvalue_Ztranform(r, n));
//		System.out.println(pvalue_spearman(r, n));
		
		double[] Y={1,2,3,4,5}; double[][] Xs={{3,5,1},{2,2,-1},{4,0,0},{2,2,0},{1,3,0}};
		OLSMultipleLinearRegression reg1=new OLSMultipleLinearRegression();
		reg1.setNoIntercept(true);
		reg1.newSampleData(Y, Xs);	
		double r2=reg1.calculateAdjustedRSquared();
		double[] residuals=reg1.estimateResiduals();
		
		double P=myMathLib.StatFuncs.total_reg_pvalue_F_nointercept(Y, residuals, Xs[0].length+1);
		double[] beta=reg1.estimateRegressionParameters();
		double[] ps=multi_reg_pvalues(beta, reg1.estimateRegressionParametersStandardErrors(), Y.length);
		System.out.println();	
		for(int k=0;k<residuals.length;k++){
			System.out.println(residuals[k]);
		}System.out.println();
		for(int k=0;k<ps.length;k++){
			System.out.println(ps[k]);
		}System.out.println();	
		for(int k=0;k<beta.length;k++){
			System.out.println(beta[k]);
		}
		System.out.println("\n"+r2+"/"+P);
		
	}
}











