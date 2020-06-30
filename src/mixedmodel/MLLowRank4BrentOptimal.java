package mixedmodel;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;

public class MLLowRank4BrentOptimal implements UnivariateFunction{

	public final double[][] transformed_X; // number_sites x smaller sample_size due to low rank
	public final double[][] transformed2_X; // number_sites x sample_size 
	public final double[] transformed_Y; //smaller sample_size due to low rank
	public final double[] transformed2_Y;// sample_size
	public final double[] S;
	//final double logXtX; //this is for REML that has not been implemented yet.
	public final Array2DRowRealMatrix additional_first_term;
	public final double[] additional_second_term;

	 
	public MLLowRank4BrentOptimal(double[][] X_before_transformation, double[] Y_before_transformation, 
			double[] transformed_Y, double[] transformed2_Y,
			double[] S, RealMatrix U1t, RealMatrix IminusU1U1t){
		this.transformed_X=FaSTLMM.transform_X_by_Ut(U1t, X_before_transformation); 
		this.transformed2_X=FaSTLMM.transform_X_by_Ut(IminusU1U1t, X_before_transformation); 
		this.transformed_Y=transformed_Y.clone(); 
		this.transformed2_Y=transformed2_Y.clone(); 
		this.S=S.clone();
		
		int n=X_before_transformation[0].length;//n=sample size
		int num_coef=transformed_X.length;
		double[][] additional_first_term_array=new double[num_coef][num_coef];
		for(int n1=0;n1<num_coef;n1++){
			for(int n2=0;n2<num_coef;n2++){
				for(int i=0;i<n;i++){
					additional_first_term_array[n1][n2]+=(this.transformed2_X[n1][i]*this.transformed2_X[n2][i]);
				}
			}			
		}
		this.additional_first_term=new Array2DRowRealMatrix(additional_first_term_array);
		this.additional_second_term=new double[num_coef];
		for(int k=0;k<num_coef;k++){
			for(int i=0;i<n;i++){
				additional_second_term[k]+=(transformed2_X[k][i]*this.transformed2_Y[i]);
			}
		}
	}
	
	public static double[] calculate_beta(double[] transformed_Y, double[][] transformed_X, 
			Array2DRowRealMatrix additional_first_term, double[] additional_second_term, double[] S, double delta){
		int n=transformed_Y.length; // when rank is not full, this n is not the original sample size.
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
		RealMatrix first_term_sum= (new Array2DRowRealMatrix(first_term)).add(additional_first_term.scalarMultiply(1.0/delta));
		LUDecomposition first_term_desomposition=new LUDecomposition(first_term_sum);
		RealMatrix first_term_inverse= first_term_desomposition.getSolver().getInverse() ;
		double[] beta= first_term_inverse.operate(myMathLib.ArrayFuncs.linearCombination(next_term, additional_second_term, 1, 1.0/delta));
		return beta;
	}
	
	public static double sigma_g2_LL(double[] transformed_Y, double[][] transformed_X, 
			double[] transformed2_Y, double[][] transformed2_X, double[] S, double delta, double[] beta){
		double n=transformed2_Y.length; // sample size
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
			sum2=sum2+(diff_i*diff_i);
		}
		return ((sum1+sum2/delta)/n);
	}
	
	public static double LL(double[] transformed_Y, double[][] transformed_X, double[] transformed2_Y, double[][] transformed2_X,
			Array2DRowRealMatrix additional_first_term, double[] additional_second_term, double[] S, double delta){
		double[] beta=calculate_beta(transformed_Y, transformed_X, additional_first_term, additional_second_term, S, delta);
		double n=transformed2_Y.length; // sample size
		double nk=transformed_Y.length; // smaller "sample size" after transformation.
//		int num_coef=transformed_X.length;
		double log_S_delta=0;
		for(int i=0;i<nk;i++)log_S_delta+=Math.log(S[i]+delta);
		double sigma_g2=sigma_g2_LL(transformed_Y, transformed_X, transformed2_Y, transformed2_X, S, delta, beta);
		return (-0.5*(n*Math.log(2*Math.PI)+log_S_delta+(n-nk)*Math.log(delta)+
				n+n*Math.log(sigma_g2)));
	}

	
	@Override
	public double value(double delta) {		
		return LL(transformed_Y, transformed_X, transformed2_Y, transformed2_X,
				 additional_first_term, additional_second_term, S, delta);
	}
}
