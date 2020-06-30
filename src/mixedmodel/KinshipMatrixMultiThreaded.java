package mixedmodel;

public class KinshipMatrixMultiThreaded {
    double[][] kinship;
    double[][] total_num_var_used;

    KinshipMatrixMultiThreaded(int sample_size) {
        kinship = new double[sample_size][sample_size];
        total_num_var_used = new double[sample_size][sample_size];
    }

    public void add_to_kinship(int i, int j, double val) {
        synchronized (this) {
            kinship[i][j] += val;
        }
    }

    public void increment_tot_num_var_used(int i, int j) {
        synchronized (this) {
            total_num_var_used[i][j]++;
        }
    }
}
