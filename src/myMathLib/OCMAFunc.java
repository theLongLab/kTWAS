package myMathLib;

import java.io.IOException;

public class OCMAFunc {
    public double[] eigenvalues;
    public double[][] eigenvectors_asrows;

    public OCMAFunc(double[][] inp) throws IOException, InterruptedException {
        for (int i = 0; i < inp.length; i++) {
            for (int j = i + 1; j < inp.length; j++) {
                inp[i][j] = inp[j][i];
            }
        }

        MatrixFunctions calc = new MatrixFunctions();

        calc.eigen_decomp_symm(inp);

        this.eigenvalues = calc.eigen_values;
        this.eigenvectors_asrows = calc.eigen_vectors;

        calc = null;
    }

    public void modify4emma(int num) {
        double[] updated_values = new double[num];
        double[][] updated_vectors = new double[num][this.eigenvectors_asrows[0].length];
        for (int i = 0; i < num; i++) {
            updated_values[i] = this.eigenvalues[i] - 1;
            updated_vectors[i] = this.eigenvectors_asrows[i].clone();
        }
        this.eigenvalues = updated_values;
        this.eigenvectors_asrows = updated_vectors;
    }

}
