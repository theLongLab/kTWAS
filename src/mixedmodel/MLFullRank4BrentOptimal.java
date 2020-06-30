package mixedmodel;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;

public class MLFullRank4BrentOptimal implements UnivariateFunction {

    final double[][] transformed_X; // number_sites x sample_size
    final double[] transformed_Y;
    final double[] S;
//	final double logXtX;

    public MLFullRank4BrentOptimal(double[][] transformed_X, double[] transformed_Y, double[] S) {
        this.transformed_X = transformed_X;
        this.transformed_Y = transformed_Y;
        this.S = S;
    }

    public static double[] calculate_beta(double[] transformed_Y, double[][] transformed_X, double[] S, double delta) {
        int n = transformed_Y.length;
        int num_coef = transformed_X.length;
        double[][] first_term = new double[num_coef][num_coef];
        for (int k1 = 0; k1 < num_coef; k1++) {
            for (int k2 = 0; k2 < num_coef; k2++) {
                for (int i = 0; i < n; i++) {
                    first_term[k1][k2] += (transformed_X[k1][i] * transformed_X[k2][i] / (S[i] + delta));
                }
            }
        }
        double[] next_term = new double[num_coef];
        for (int k = 0; k < num_coef; k++) {
            for (int i = 0; i < n; i++) {
                next_term[k] += (transformed_X[k][i] * transformed_Y[i] / (S[i] + delta));
            }
        }
//		RealMatrix matrix=MatrixUtils.blockInverse((RealMatrix)(new Array2DRowRealMatrix(first_term)),(int)n-1);	
        LUDecomposition first_term_decomposition = new LUDecomposition(new Array2DRowRealMatrix(first_term));
        RealMatrix first_term_inverse = null;
        try {
            first_term_inverse = first_term_decomposition.getSolver().getInverse();
        } catch (Exception e) {
            return null;
        }
        double[] beta = first_term_inverse.operate(next_term);
        return beta;
    }

    public static double sigma_g2_LL(double[] transformed_Y, double[][] transformed_X, double[] S, double delta, double[] beta) {
        int n = transformed_Y.length;
        int num_coef = transformed_X.length;
        double sum = 0;
        for (int i = 0; i < n; i++) {
            double diff_i = transformed_Y[i];
            for (int k = 0; k < num_coef; k++) {
                diff_i = diff_i - beta[k] * transformed_X[k][i];
            }
            sum = sum + (diff_i * diff_i) / (S[i] + delta);
        }
        return (sum / n);
    }

    public static double LL(double[] transformed_Y, double[][] transformed_X, double[] S, double delta) {
        double[] beta = calculate_beta(transformed_Y, transformed_X, S, delta);
        if (beta == null) { // singular matrix!
            return -Double.MAX_VALUE;
        }
        double n = transformed_Y.length;
//		int num_coef=transformed_X.length;
        double log_S_delta = 0;
        for (int i = 0; i < n; i++) log_S_delta += Math.log(S[i] + delta);
        double sigma_g2 = sigma_g2_LL(transformed_Y, transformed_X, S, delta, beta);
        return (-0.5 * (n * Math.log(2 * Math.PI) + log_S_delta + n + n * Math.log(sigma_g2)));
    }

    @Override
    public double value(double delta) {
        return LL(transformed_Y, transformed_X, S, delta);
    }


}
