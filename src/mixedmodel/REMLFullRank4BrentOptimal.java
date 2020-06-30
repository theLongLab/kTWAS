package mixedmodel;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;

public class REMLFullRank4BrentOptimal implements UnivariateFunction{

	final double[][] transformed_X; // number_sites x sample_size
	final double[] transformed_Y;
	final double[] S;
	final double logXtX;
	
	public REMLFullRank4BrentOptimal(double[][] transformed_X, double[] transformed_Y, double[] S){
		this.transformed_X=transformed_X; 
		this.transformed_Y=transformed_Y; 
		this.S=S;
		Array2DRowRealMatrix XT=(new Array2DRowRealMatrix(transformed_X));
		Array2DRowRealMatrix X=(Array2DRowRealMatrix)XT.transpose();
		this.logXtX=Math.log((new LUDecomposition((X.preMultiply(XT)))).getDeterminant());
		
	}
	
	public static double[] calculate_beta(double[] transformed_Y, double[][] transformed_X, double[] S, double delta){
		int n=transformed_Y.length;
		int num_coef=transformed_X.length;
		double[][] first_term=new double[num_coef][num_coef];
		for(int k1=0;k1<num_coef;k1++){
			for(int k2=0;k2<num_coef;k2++){
				for(int i=0;i<n;i++){
					first_term[k1][k2]+=transformed_X[k1][i]*transformed_X[k2][i]/(S[i]+delta);
				}
			}			
		}
		double[] next_term=new double[num_coef];
		for(int k=0;k<num_coef;k++){
			for(int i=0;i<n;i++){
				next_term[k]+=(transformed_X[k][i]*transformed_Y[i]/(S[i]+delta));
			}
		}
//		RealMatrix matrix=MatrixUtils.blockInverse((RealMatrix)(new Array2DRowRealMatrix(first_term)),(int)n-1);	
		LUDecomposition first_term_desomposition=new LUDecomposition(new Array2DRowRealMatrix(first_term));
		RealMatrix first_term_inverse= first_term_desomposition.getSolver().getInverse() ;
		double[] beta= first_term_inverse.operate(next_term);
		return beta;
	}
	
	public static double sigma_g2_LL(double[] transformed_Y, double[][] transformed_X, double[] S, double delta, double[] beta){
		int n=transformed_Y.length;
		int num_coef=transformed_X.length;
		double sum=0;
		for(int i=0;i<n;i++){
			double diff_i=transformed_Y[i];
			for(int k=0;k<num_coef;k++){
				diff_i=diff_i-beta[k]*transformed_X[k][i];
			}
			sum=sum+(diff_i*diff_i)/(S[i]+delta);
		}
		return (sum/n);
	}
	
	public static double sigma_g2_REMLL(double[] transformed_Y, double[][] transformed_X, double[] S, double delta, double[] beta){
		int n=transformed_Y.length;
		int num_coef=transformed_X.length;
		double sum=0;
		for(int i=0;i<n;i++){
			double diff_i=transformed_Y[i];
			for(int k=0;k<num_coef;k++){
				diff_i=diff_i-beta[k]*transformed_X[k][i];
			}
			sum=sum+(diff_i*diff_i)/(S[i]+delta);
		}
		return (sum/(n-num_coef));
	}
	
//	public static double fisher_info_nobeta(double[] transformed_Y, double[] S, double h2){
//		double n=S.length;
//		double fisher_info=0;
//		double delta=1/h2-1;
//		double delta_d1=-1/(h2*h2);
//		double delta_d2=2/(h2*h2*h2);
//		
//		double first=0, second=0, third=0, noy_first=0, noy_second=0;
//		for(int i=0;i<transformed_Y.length;i++){
//			noy_first=noy_first+1.0/(S[i]+delta);
//			noy_second=noy_second+1.0/((S[i]+delta)*(S[i]+delta));
//			first=first+(transformed_Y[i])*(transformed_Y[i])/(S[i]+delta);
//			second=second+(transformed_Y[i])*(transformed_Y[i])/((S[i]+delta)*(S[i]+delta));
//			third=third+(transformed_Y[i])*(transformed_Y[i])/((S[i]+delta)*(S[i]+delta)*(S[i]+delta));
//		}
//		double LL_d1= -0.5*(n*(-second)/first+noy_first); 
//		double LL_d2= -0.5*(n*(2*third*first)-second*second)/(first*first)-noy_second;
//		fisher_info=delta_d2*LL_d1+LL_d2*delta_d1*delta_d1;
//		return fisher_info;
//	}
	
	public static double LL(double[] transformed_Y, double[][] transformed_X, double[] S, double delta){
		double[] beta=calculate_beta(transformed_Y, transformed_X, S, delta);
		double n=transformed_Y.length;
//		int num_coef=transformed_X.length;
		double log_S_delta=0;
		for(int i=0;i<n;i++)log_S_delta+=Math.log(S[i]+delta);
		double sigma_g2=sigma_g2_LL(transformed_Y, transformed_X, S, delta, beta);
		return (-0.5*(n*Math.log(2*Math.PI)+log_S_delta+n+n*Math.log(sigma_g2)));
	}
	
	public static double REMLL(double[] transformed_Y, double[][] transformed_X, double[] S, double delta, double logXtX){
		int n=transformed_Y.length;
		int num_coef=transformed_X.length;
		// calculate beta, the reason why not invoke the method is that the quadratic form is retained for other calculations
		double[][] first_term=new double[num_coef][num_coef];
		for(int k1=0;k1<num_coef;k1++){
			for(int k2=0;k2<num_coef;k2++){
				for(int i=0;i<n;i++){
					first_term[k1][k2]+=transformed_X[k1][i]*transformed_X[k2][i]/(S[i]+delta);
				}
			}			
		}
		double[] next_term=new double[num_coef];
		for(int k=0;k<num_coef;k++){
			for(int i=0;i<n;i++){
				next_term[k]+=(transformed_X[k][i]*transformed_Y[i]/(S[i]+delta));
			}
		}
		LUDecomposition first_term_desomposition=new LUDecomposition(new Array2DRowRealMatrix(first_term));
		RealMatrix first_term_inverse= first_term_desomposition.getSolver().getInverse() ;
		double[] beta= first_term_inverse.operate(next_term);
		
		double log_S_delta=0;
		for(int i=0;i<n;i++)log_S_delta+=Math.log(S[i]+delta);
		
		double sigma_g2=sigma_g2_REMLL(transformed_Y, transformed_X, S, delta, beta);
		
		
		
		double determinant_first_term=first_term_desomposition.getDeterminant();
		double REMLL= -0.5*((n-num_coef)*Math.log(2*Math.PI)+log_S_delta+(n-num_coef)+
				(n-num_coef)*Math.log(sigma_g2))+Math.log(determinant_first_term)-(logXtX);
		return REMLL;
	}
	
	@Override
	public double value(double delta) {		
		return REMLL(transformed_Y, transformed_X, S, delta, logXtX);
	}
}
