package mixedmodel;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

//import mixed_model.RootSearchEmmadeltaMLdLLwoZ;
//import mixed_model.RootSearchEmmadeltaREMLdLLwoZ;
import myMathLib.OCMAFunc;
import myMathLib.Test;
import flanagan.math.Matrix;
import flanagan.roots.RealRoot;

public class EMMA_OCMA {
    // constants in the calculation
    public static final double ngrids = 100;
    public static final double llim = -10;
    public static final double ulim = 10;
    public static final double esp = 1e-10;
    public Phenotype phenotype; // all samples with phenotype will have genotype in hdf5: the constructor will guarantee that.
    public int sample_size;
    // fields for MLE:
    public double mle_ML;
    public double mle_delta;
    public double mle_vg;
    public double mle_ve;
    public boolean mle_calculated = false;
    double[][] kinship_matrix; // this is corresponding to the phenotype, but generated from genotype matched kinship file.
    VariantsDouble genotype;
    int[] indexes_with_phenotype_in_genotype;
    int[] mafc;
    OCMAFunc reml_eig_L = null;
    OCMAFunc ml_eig_R = null;
    // fields for REMLE:
    double remle_REML;
    double remle_delta;
    double remle_vg;
    double remle_ve;
    double heritability;
    boolean remle_calculated = false;

    public EMMA_OCMA() {

    }

    /*
     * One has to make sure that the original kinship_matrix is corresponding to the sample_ids in genotype
     *
     *  Find the sample_ids in both genotype and ori_phenotype. generate a clean phenotype dataset for further analysis.
     */
    public EMMA_OCMA(Phenotype ori_phenotype, VariantsDouble genotype, double[][] ori_kinship_matrix) {
        this.genotype = genotype;
//		String[] id_phenotype=ori_phenotype.sample_ids;
        // find the ids
        ArrayList<Integer> useful_phe_ids = new ArrayList<>();
        HashMap<Integer, Integer> corresponding_geno_ids = new HashMap<>();
        HashMap<String, Integer> all_geno_ids = new HashMap<>();
        for (int geno_id_index = 0; geno_id_index < genotype.sample_ids.length; geno_id_index++) {
            all_geno_ids.put(genotype.sample_ids[geno_id_index], geno_id_index);
        }
        for (int phe_id_index = 0; phe_id_index < ori_phenotype.sample_ids.length; phe_id_index++) {
            if (all_geno_ids.containsKey(ori_phenotype.sample_ids[phe_id_index])) {
                int geno_id_index = all_geno_ids.get(ori_phenotype.sample_ids[phe_id_index]);
                useful_phe_ids.add(phe_id_index);
                corresponding_geno_ids.put(phe_id_index, geno_id_index);
            }
        }
        // assign the fields
        this.sample_size = useful_phe_ids.size();
        double[] the_values = new double[sample_size];
        String[] the_ids = new String[sample_size];
        this.indexes_with_phenotype_in_genotype = new int[sample_size];
        for (int k = 0; k < sample_size; k++) {
            int phe_id_index = useful_phe_ids.get(k);
            the_ids[k] = ori_phenotype.sample_ids[phe_id_index];
            the_values[k] = ori_phenotype.values[phe_id_index];
            indexes_with_phenotype_in_genotype[k] = corresponding_geno_ids.get(phe_id_index);
        }
        this.phenotype = new Phenotype(ori_phenotype.phe_id, the_ids, the_values);
        // extract kinship
        this.kinship_matrix = new double[sample_size][sample_size];
        for (int i = 0; i < sample_size; i++) {
            int gen_id_index1 = indexes_with_phenotype_in_genotype[i];
            for (int j = 0; j < sample_size; j++) {
                int gen_id_index2 = indexes_with_phenotype_in_genotype[j];
                kinship_matrix[i][j] = ori_kinship_matrix[gen_id_index1][gen_id_index2];
            }
        }//Test.write2file(the_ids, "/Users/quan.long/look/ids.new");
    }

    /*
     * One has to make sure that the original kinship_matrix is corresponding to the sample_ids in genotype
     *
     *  Find the sample_ids in both genotype and ori_phenotype. generate a clean phenotype dataset for further analysis.
     */
    public EMMA_OCMA(Phenotype ori_phenotype, VariantsDouble genotype, double[][] ori_kinship_matrix, boolean kinship_sample_selected) {
        this.genotype = genotype;
//		String[] id_phenotype=ori_phenotype.sample_ids;
        // find the ids
        ArrayList<Integer> useful_phe_ids = new ArrayList<>();
        HashMap<Integer, Integer> corresponding_geno_ids = new HashMap<>();
        HashMap<String, Integer> all_geno_ids = new HashMap<>();
        for (int geno_id_index = 0; geno_id_index < genotype.sample_ids.length; geno_id_index++) {
            all_geno_ids.put(genotype.sample_ids[geno_id_index], geno_id_index);
        }
        for (int phe_id_index = 0; phe_id_index < ori_phenotype.sample_ids.length; phe_id_index++) {
            if (all_geno_ids.containsKey(ori_phenotype.sample_ids[phe_id_index])) {
                int geno_id_index = all_geno_ids.get(ori_phenotype.sample_ids[phe_id_index]);
                useful_phe_ids.add(phe_id_index);
                corresponding_geno_ids.put(phe_id_index, geno_id_index);
            }
        }
        // assign the fields
        this.sample_size = useful_phe_ids.size();
        double[] the_values = new double[sample_size];
        String[] the_ids = new String[sample_size];
        this.indexes_with_phenotype_in_genotype = new int[sample_size];
        for (int k = 0; k < sample_size; k++) {
            int phe_id_index = useful_phe_ids.get(k);
            the_ids[k] = ori_phenotype.sample_ids[phe_id_index];
            the_values[k] = ori_phenotype.values[phe_id_index];
            indexes_with_phenotype_in_genotype[k] = corresponding_geno_ids.get(phe_id_index);
        }
        this.phenotype = new Phenotype(ori_phenotype.phe_id, the_ids, the_values);
        // NO extracting kinship since samples are already selected
        this.kinship_matrix = ori_kinship_matrix;
//		for(int i=0;i<sample_size;i++){
//			int gen_id_index1=indexes_with_phenotype_in_genotype[i];
//			for(int j=0;j<sample_size;j++){
//				int gen_id_index2=indexes_with_phenotype_in_genotype[j];
//				kinship_matrix[i][j]=ori_kinship_matrix[gen_id_index1][gen_id_index2];
//			}
//		}//Test.write2file(the_ids, "/Users/quan.long/look/ids.new");
    }

    public static OCMAFunc eigen_R_wo_Z(double[][] X, double[][] K) throws IOException, InterruptedException {
        int q = X[0].length;
        int sample_size = K.length;
        Matrix XX = new Matrix(X);
        Matrix Xt = XX.transpose();
        Matrix S = Matrix.times(Matrix.times(XX, (Matrix.times(Xt, XX)).inverse()), Xt);
        double[][] diag = new double[sample_size][sample_size];
        double[][] Kplus1 = new double[sample_size][sample_size];
        for (int i = 0; i < sample_size; i++) {
            diag[i][i] = 1;
        }
        for (int i = 0; i < sample_size; i++) {
            for (int j = 0; j < sample_size; j++) {
                Kplus1[i][j] = K[i][j] + diag[i][j];
            }
        }
        S = S.times(-1).plus(diag);

        Matrix final_array = Matrix.times(Matrix.times(S, Kplus1), S);
        for (int i = 0; i < sample_size; i++) {
            for (int j = i + 1; j < sample_size; j++)
                final_array.setElement(i, j, final_array.getElement(j, i));
        }
        OCMAFunc result = new OCMAFunc(final_array.getArrayCopy());
        diag = null;
        Kplus1 = null;
        S = null;
        Xt = null;
        XX = null;
        result.modify4emma(sample_size - q);
        return result;
    }

    public static double delta_REML_LL_wo_Z(double logdelta, double[] lambda, double[] etas) {
        int nq = etas.length;
        double delta = Math.exp(logdelta);
        double sum1 = 0, sum2 = 0;
        for (int i = 0; i < etas.length; i++) {
            sum1 += (etas[i] * etas[i] / (lambda[i] + delta));
            sum2 += (Math.log(lambda[i] + delta));
        }
        if (sum1 > 0)
            return 0.5 * (nq * (Math.log(nq / (2 * Math.PI)) - 1 - Math.log(sum1)) - sum2);
        else
            return 0.5 * (nq * (Math.log(nq / (2 * Math.PI)) - 1 - Math.log(sum1)) - sum2);
    }

    /*
     * The rest of this class is for EMMA/REML calculations
     */

    public static VarComp_Result REMLE_static(double[] y, double[][] X, double[][] kinship_matrix) {
        VarComp_Result result = new VarComp_Result();
        try {
            int n = y.length;
            int t = kinship_matrix.length;
            int q = X[0].length; // number of fixed effets;
            Matrix XX = new Matrix(X);
            Matrix Xt = XX.transpose();
            if (Matrix.times(Xt, XX).determinant() == 0) {
                System.out.println("X is singular");
                return null;
            }
            //optlogdelta <- vector(length=0)
            ArrayList<Double> optlogdelta = new ArrayList<>();
            //optLL <- vector(length=0)
            ArrayList<Double> optLL = new ArrayList<>();
            double[] etas = new double[n - q];

            OCMAFunc ml_eig_R = eigen_R_wo_Z(X, kinship_matrix);

            // etas <- crossprod(eig.R$vectors,y) has been implemented as below:

            for (int i = 0; i < etas.length; i++) {
                for (int j = 0; j < n; j++) etas[i] += ml_eig_R.eigenvectors_asrows[i][j] * y[j];
            }
            //logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
            //delta <- exp(logdelta)
            double[] logdelta = new double[(int) EMMA.ngrids + 1];
            double[] delta = new double[(int) EMMA.ngrids + 1];
            for (int i = 0; i <= EMMA.ngrids; i++) {
                logdelta[i] = i / EMMA.ngrids * (EMMA.ulim - EMMA.llim) + EMMA.llim;
                delta[i] = Math.exp(logdelta[i]);
            }
            //m <- length(logdelta)
            int m = logdelta.length; // the grids number
            //Lambdas <- matrix(eig.R$values,n-q,m) + matrix(delta,n-q,m,byrow=TRUE)
            double[][] lambdas = new double[n - q][m];
            for (int i = 0; i < n - q; i++) {
                for (int j = 0; j < m; j++) lambdas[i][j] = ml_eig_R.eigenvalues[i] + delta[j];
            }
//				//Xis <- matrix(eig.L$values,n,m) + matrix(delta,n,m,byrow=TRUE)
//				double[][] Xis= new double[n][m];
//				for(int i=0;i<n;i++){
//					for(int j=0;j<m;j++)Xis[i][j]=eig_L.eigenvalues[i]+delta[j];
//				}
            //Etasq <- matrix(etas*etas,n-q,m)
            double[][] etasq = new double[n - q][m];
            for (int i = 0; i < n - q; i++) {
                for (int j = 0; j < m; j++) etasq[i][j] = etas[i] * etas[i];
            }
//				LL <- 0.5*((n-q)*(log((n-q)/(2*pi))-1-log(colSums(Etasq/Lambdas)))-colSums(log(Lambdas)))
            double[] LL = new double[m];
            for (int i = 0; i < m; i++) {
                double sum1 = 0, sum2 = 0;
                for (int j = 0; j < n - q; j++) {
                    sum1 += (etasq[j][i] / lambdas[j][i]);
                    sum2 += (Math.log(lambdas[j][i]));
                }
                LL[i] = 0.5 * ((n - q) * (Math.log((n - q) / (2 * Math.PI)) - 1 - Math.log(sum1)) - sum2);
            }
            // dLL <- 0.5*delta*((n-q)*colSums(Etasq/(Lambdas*Lambdas))/colSums(Etasq/Lambdas)-colSums(1/Lambdas))
            double[] dLL = new double[m];
            for (int i = 0; i < m; i++) {
                double sum1 = 0, sum2 = 0, sum3 = 0;
                for (int j = 0; j < n - q; j++) {
                    sum1 += (etasq[j][i] / (lambdas[j][i] * lambdas[j][i]));
                    sum2 += (etasq[j][i] / lambdas[j][i]);
                    sum3 += (1 / lambdas[j][i]);
                }
                dLL[i] = 0.5 * delta[i] * ((n - q) * sum1 / sum2 - sum3);
            }
            if (dLL[0] < EMMA.esp) {
                optlogdelta.add((double) EMMA.llim); //optlogdelta <- append(optlogdelta, llim)
                optLL.add(EMMA.delta_REML_LL_wo_Z(EMMA.llim, ml_eig_R.eigenvalues, etas));//optLL <- append(optLL, emma.delta.ML.LL.wo.Z(llim,eig.R$values,etas,eig.L$values))
            }
            if (dLL[m - 1] > 0 - EMMA.esp) {
                optlogdelta.add((double) EMMA.ulim);
                optLL.add(delta_REML_LL_wo_Z(EMMA.ulim, ml_eig_R.eigenvalues, etas));
            }
            //for( i in 1:(m-1) )
            for (int i = 0; i < m - 1; i++) {
                if ((dLL[i] * dLL[i + 1] < 0 - EMMA.esp * EMMA.esp) && (dLL[i] > 0) && (dLL[i + 1] < 0)) {
                    RealRoot rlrt = new RealRoot();
                    RootSearchEmmadeltaREMLdLLwoZ func = new RootSearchEmmadeltaREMLdLLwoZ();
                    func.set_delta_REML_dLL_wo_Z(ml_eig_R.eigenvalues, etas);
                    double root = rlrt.brent(func, logdelta[i], logdelta[i + 1]);
                    //optlogdelta <- append(optlogdelta, r$root)
                    optlogdelta.add(root);
                    //optLL <- append(optLL, emma.delta.ML.LL.wo.Z(r$root,eig.R$values, etas, eig.L$values))
                    optLL.add(delta_REML_LL_wo_Z(root, ml_eig_R.eigenvalues, etas));
                }
            }

            int max_index = 0;
//			while(Double.isNaN(optLL.get(max_index)))max_index++;
            for (int i = 1; i < optlogdelta.size(); i++) {
                if (optLL.get(i) > optLL.get(max_index)) max_index = i;
            }
            double maxdelta = Math.exp(optlogdelta.get(max_index)); //maxdelta <- exp(optlogdelta[which.max(optLL)])
            double maxLL = optLL.get(max_index);//maxLL <- max(optLL)
            double maxva = 0;

            for (int i = 0; i < etas.length; i++) maxva += etas[i] * etas[i] / (ml_eig_R.eigenvalues[i] + maxdelta);
            maxva = maxva / (n - q);

            double maxve = maxva * maxdelta;//maxve <- maxva*maxdelta
            result = new VarComp_Result();
            result.ml = maxLL;
            result.delta = maxdelta;
            result.sigma_e = maxve;
            result.sigma_g = maxva;
//			System.out.println("Sample size: "+this.sample_size);
//			System.out.println("remle_REML="+maxLL);
//			System.out.println("remle_delta="+maxdelta);
//			System.out.println("remle_ve="+maxve);
//			System.out.println("remle_vg="+maxva);
//			System.out.println("remle_pseudo-heritability="+maxva/(maxva+maxve));
        } catch (Exception e) {
            e.printStackTrace();
        }
        return result;
    }

    public void clear() {
        this.ml_eig_R = null;
        this.reml_eig_L = null;
        this.mle_ML = Double.NaN;
        this.mle_ve = Double.NaN;
        this.mle_vg = Double.NaN;

        this.heritability = Double.NaN;
        this.remle_delta = Double.NaN;
        this.remle_REML = Double.NaN;
        this.remle_ve = Double.NaN;
        this.remle_vg = Double.NaN;
        this.kinship_matrix = null;
        this.remle_calculated = false;
        this.mle_calculated = false;
    }

    public void REMLE_null() {
        int[] indexes = {}, chrs = {};
        this.REMLE(this.phenotype.values, this.cbind(chrs, indexes), null);
    }

    public OCMAFunc eigen_L(double[][] Z) throws IOException, InterruptedException {
        double[][] K = this.kinship_matrix;
        if (Z == null) {
            return eigen_L_wo_Z();
        } else {
            System.out.println("non trivial Z part of java-EMMA has not been implemented yet");
            return null;
        }
    }

    /*
     * translated from EMMA/emma.eigen.L.wo.Z(K)
     */
    public OCMAFunc eigen_L_wo_Z() throws IOException, InterruptedException {
        double[][] K = this.kinship_matrix;
        return new OCMAFunc(K);
    }

    /*
     *  translated from EMMA/emma.eigen.R <- function(Z,K,X,complete=TRUE)
     *  X is an n x q array for fixed effects
     */
    public OCMAFunc eigen_R(double[][] Z, double[][] X) throws IOException, InterruptedException {
        double[][] K = this.kinship_matrix;
        if (Z == null) {
            return eigen_R_wo_Z(X);
        } else {
            System.out.println("non trivial Z part of java-EMMA has not been implemented yet");
            return null;
        }
    }


    /*
     * translated from EMMA/emma.delta.REML.LL.wo.Z <- function(logdelta, lambda, etas)
     */

    /*
     *  translated from EMMA/emma.eigen.R.wo.Z <- function(K, X)
     *  X is an n x q array for fixed effects
     */
    public OCMAFunc eigen_R_wo_Z(double[][] X) throws IOException, InterruptedException {
        double[][] K = this.kinship_matrix;
        int q = X[0].length;
        Matrix XX = new Matrix(X);
        Matrix Xt = XX.transpose();
        Matrix S = Matrix.times(Matrix.times(XX, (Matrix.times(Xt, XX)).inverse()), Xt);
        double[][] diag = new double[this.sample_size][this.sample_size];
        double[][] Kplus1 = new double[this.sample_size][this.sample_size];
        for (int i = 0; i < this.sample_size; i++) {
            diag[i][i] = 1;
        }
        for (int i = 0; i < this.sample_size; i++) {
            for (int j = 0; j < this.sample_size; j++) {
                Kplus1[i][j] = K[i][j] + diag[i][j];
            }
        }
        S = S.times(-1).plus(diag);

        Matrix final_array = Matrix.times(Matrix.times(S, Kplus1), S);
        for (int i = 0; i < this.sample_size; i++) {
            for (int j = i + 1; j < this.sample_size; j++)
                final_array.setElement(i, j, final_array.getElement(j, i));
        }
        OCMAFunc result = new OCMAFunc(final_array.getArrayCopy());
        diag = null;
        Kplus1 = null;
        S = null;
        Xt = null;
        XX = null;
        result.modify4emma(this.sample_size - q);
        return result;
    }

    /*
     * translated from EMMA/emma.delta.ML.LL.wo.Z <- function(logdelta, lambda, etas, xi)
     */
    public double delta_ML_LL_wo_Z(double logdelta, double[] lambda, double[] etas, double[] xi) {

        int n = xi.length;
        if (n != this.sample_size) System.out.println("length of xi wrong");
        double delta = Math.exp(logdelta);
        double sum1 = 0, sum2 = 0;
        for (int i = 0; i < etas.length; i++) {
            sum1 += (etas[i] * etas[i] / (lambda[i] + delta));
        }
        for (int i = 0; i < xi.length; i++) {
            sum2 += (Math.log(xi[i] + delta));
        }
        return 0.5 * (n * (Math.log(n / (2 * Math.PI)) - 1 - Math.log(sum1)) - sum2);
    }

    /*
     *  translated from EMMA/emma.delta.ML.dLL.wo.Z <- function(logdelta, lambda, etas, xi)
     */
    public double delta_ML_dLL_wo_Z(double logdelta, double[] lambda, double[] etas, double[] xi) {
        double n = xi.length;
        if (n != this.sample_size) System.out.println("length of xi wrong");
        double delta = Math.exp(logdelta);
        //etasq <- etas*etas
        //ldelta <- lambda+delta
        double sum1 = 0, sum2 = 0, sum3 = 0;
        for (int i = 0; i < etas.length; i++) {
            sum1 += etas[i] * etas[i] / ((lambda[i] + delta) * (lambda[i] + delta));
            sum2 += etas[i] * etas[i] / (lambda[i] + delta);
        }
        for (int i = 0; i < xi.length; i++) {
            sum3 += 1.0 / (xi[i] + delta);
        }
        return 0.5 * (n * sum1 / sum2 - sum3);
    }

    /*
     * translated from EMMA/emma.MLE <- function(y, X, K, Z=NULL, ngrids=100, llim=-10, ulim=10,
     *	esp=1e-10, eig.L = NULL, eig.R = NULL)
     *
     *	y is the phenotype, X is the fixed effect, K is the kinship and Z is the incidence matrix;
     */
    public void MLE(double[] y, double[][] X, double[][] Z) throws IOException, InterruptedException {
        int n = y.length;
        int t = this.kinship_matrix.length;
        int q = X[0].length; // number of fixed effets;
        Matrix XX = new Matrix(X);
        Matrix Xt = XX.transpose();
        if (Matrix.times(Xt, XX).determinant() == 0) {
            System.out.println("X is singular");
            return;
        }
        //optlogdelta <- vector(length=0)
        ArrayList<Double> optlogdelta = new ArrayList<>();
        //optLL <- vector(length=0)
        ArrayList<Double> optLL = new ArrayList<>();
        double[] etas = new double[n - q];
        if (Z == null) {
            if (this.reml_eig_L == null) {
                this.reml_eig_L = this.eigen_L_wo_Z();
            }
            if (this.ml_eig_R == null) {
                this.ml_eig_R = eigen_R_wo_Z(X);
            }
            // etas <- crossprod(eig.R$vectors,y) has been implemented as below:

            for (int i = 0; i < etas.length; i++) {
                for (int j = 0; j < this.sample_size; j++) etas[i] += this.ml_eig_R.eigenvectors_asrows[i][j] * y[j];
            }
            //logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
            //delta <- exp(logdelta)
            double[] logdelta = new double[(int) this.ngrids + 1];
            double[] delta = new double[(int) this.ngrids + 1];
            for (int i = 0; i <= this.ngrids; i++) {
                logdelta[i] = i / this.ngrids * (this.ulim - this.llim) + llim;
                delta[i] = Math.exp(logdelta[i]);
            }
            //m <- length(logdelta)
            int m = logdelta.length; // the grids number


            //Lambdas <- matrix(eig.R$values,n-q,m) + matrix(delta,n-q,m,byrow=TRUE)
            double[][] lambdas = new double[n - q][m];
            for (int i = 0; i < n - q; i++) {
                for (int j = 0; j < m; j++) lambdas[i][j] = ml_eig_R.eigenvalues[i] + delta[j];
            }
            //Xis <- matrix(eig.L$values,n,m) + matrix(delta,n,m,byrow=TRUE)
            double[][] Xis = new double[n][m];
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < m; j++) Xis[i][j] = reml_eig_L.eigenvalues[i] + delta[j];
            }
            //Etasq <- matrix(etas*etas,n-q,m)
            double[][] etasq = new double[n - q][m];
            for (int i = 0; i < n - q; i++) {
                for (int j = 0; j < m; j++) etasq[i][j] = etas[i] * etas[i];
            }
            //LL <- 0.5*(n*(log(n/(2*pi))-1-log(colSums(Etasq/Lambdas)))-colSums(log(Xis)))
            double[] LL = new double[m];
            for (int i = 0; i < m; i++) {
                double sum1 = 0, sum2 = 0;
                for (int j = 0; j < n - q; j++) {
                    sum1 += (etasq[j][i] / lambdas[j][i]);
                }
                for (int j = 0; j < n; j++) {
                    sum2 += (Math.log(Xis[j][i]));
                }
                LL[i] = 0.5 * (n * (Math.log(n / (2 * Math.PI)) - 1 - Math.log(sum1)) - sum2);
            }
            //dLL <- 0.5*delta*(n*colSums(Etasq/(Lambdas*Lambdas))/colSums(Etasq/Lambdas)-colSums(1/Xis))
            double[] dLL = new double[m];
            for (int i = 0; i < m; i++) {
                double sum1 = 0, sum2 = 0, sum3 = 0;
                for (int j = 0; j < n - q; j++) {
                    sum1 += (etasq[j][i] / (lambdas[j][i] * lambdas[j][i]));
                    sum2 += (etasq[j][i] / lambdas[j][i]);
                }
                for (int j = 0; j < n; j++) {
                    sum3 += (1 / Xis[j][i]);
                }
                dLL[i] = 0.5 * delta[i] * (n * sum1 / sum2 - sum3);
            }

            if (dLL[0] < this.esp) {
                optlogdelta.add((double) llim); //optlogdelta <- append(optlogdelta, llim)
                optLL.add(this.delta_ML_LL_wo_Z(this.llim, ml_eig_R.eigenvalues, etas, reml_eig_L.eigenvalues));//optLL <- append(optLL, emma.delta.ML.LL.wo.Z(llim,eig.R$values,etas,eig.L$values))
            }
            if (dLL[m - 1] > 0 - this.esp) {
                //optlogdelta <- append(optlogdelta, ulim)
                optlogdelta.add((double) this.ulim);
                //optLL <- append(optLL, emma.delta.ML.LL.wo.Z(ulim,eig.R$values,etas,eig.L$values))
                optLL.add(this.delta_ML_LL_wo_Z(this.ulim, ml_eig_R.eigenvalues, etas, reml_eig_L.eigenvalues));
            }

            //for( i in 1:(m-1) )
            for (int i = 0; i < m - 1; i++) {
                if ((dLL[i] * dLL[i + 1] < 0 - this.esp * this.esp) && (dLL[i] > 0) && (dLL[i + 1] < 0)) {
                    //r <- uniroot(emma.delta.ML.dLL.wo.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.R$values, etas=etas, xi=eig.L$values)
                    RealRoot rlrt = new RealRoot();
                    RootSearchEmmadeltaMLdLLwoZ func = new RootSearchEmmadeltaMLdLLwoZ();
                    func.set_delta_ML_dLL_wo_Z(ml_eig_R.eigenvalues, etas, reml_eig_L.eigenvalues);
                    double root = rlrt.brent(func, logdelta[i], logdelta[i + 1]);
                    //optlogdelta <- append(optlogdelta, r$root)
                    optlogdelta.add(root);
                    //optLL <- append(optLL, emma.delta.ML.LL.wo.Z(r$root,eig.R$values, etas, eig.L$values))
                    optLL.add(this.delta_ML_LL_wo_Z(root, ml_eig_R.eigenvalues, etas, reml_eig_L.eigenvalues));
                }
            }
        } else {
            System.out.println("non trivial Z part of java-EMMA has not been implemented yet");
            return;
        }
        int max_index = 0;
        for (int i = 1; i < optlogdelta.size(); i++) {
            if (optLL.get(i) > optLL.get(max_index)) max_index = i;
        }
        double maxdelta = Math.exp(optlogdelta.get(max_index)); //maxdelta <- exp(optlogdelta[which.max(optLL)])
        double maxLL = optLL.get(max_index);//maxLL <- max(optLL)
        double maxva = 0;
        if (Z == null) {
            //maxva <- sum(etas*etas/(eig.R$values+maxdelta))/n
            for (int i = 0; i < etas.length; i++) maxva += etas[i] * etas[i] / (ml_eig_R.eigenvalues[i] + maxdelta);
            maxva = maxva / n;
        } else {
            System.out.println("non trivial Z part of java-EMMA has not been implemented yet");
            //maxva <- (sum(etas.1*etas.1/(eig.R$values+maxdelta))+etas.2.sq/maxdelta)/n
        }
        double maxve = maxva * maxdelta;//maxve <- maxva*maxdelta

        //return (list(ML=maxLL,delta=maxdelta,ve=maxve,vg=maxva))
        this.mle_ML = maxLL;
        this.mle_delta = maxdelta;
        this.mle_ve = maxve;
        this.mle_vg = maxva;
        System.out.println("this.mle_ML=" + maxLL);
        System.out.println("this.mle_delta=" + maxdelta);
        System.out.println("this.mle_ve=" + maxve);
        System.out.println("this.mle_vg=" + maxva);
        return;
    }

    /*
     * translated from EMMA/emma.MLE <- function(y, X, K, Z=NULL, ngrids=100, llim=-10, ulim=10,
     *	esp=1e-10, eig.L = NULL, eig.R = NULL)
     *
     *	y is the phenotype, X is the fixed effect, K is the kinship and Z is the incidence matrix;
     */
    public void REMLE(double[] y, double[][] X, double[][] Z) {
//		Test.write2file(y, "/Users/quan.long/look/y.new");
//		Test.write2file(X, "/Users/quan.long/look/X.new");
//		Test.write2file(global_kinship_matrix, "/Users/quan.long/look/k.new");
        try {
            int n = y.length;
            int t = this.kinship_matrix.length;
            int q = X[0].length; // number of fixed effets;
            Matrix XX = new Matrix(X);
            Matrix Xt = XX.transpose();
            if (Matrix.times(Xt, XX).determinant() == 0) {
                System.out.println("X is singular");
                return;
            }
            //optlogdelta <- vector(length=0)
            ArrayList<Double> optlogdelta = new ArrayList<Double>();
            //optLL <- vector(length=0)
            ArrayList<Double> optLL = new ArrayList<Double>();
            double[] etas = new double[n - q];
            if (Z == null) {
                if (this.ml_eig_R == null) {
                    this.ml_eig_R = eigen_R_wo_Z(X);
                }
                // etas <- crossprod(eig.R$vectors,y) has been implemented as below:

                for (int i = 0; i < etas.length; i++) {
                    for (int j = 0; j < this.sample_size; j++)
                        etas[i] += this.ml_eig_R.eigenvectors_asrows[i][j] * y[j];
                }
                //logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
                //delta <- exp(logdelta)
                double[] logdelta = new double[(int) this.ngrids + 1];
                double[] delta = new double[(int) this.ngrids + 1];
                for (int i = 0; i <= this.ngrids; i++) {
                    logdelta[i] = i / this.ngrids * (this.ulim - this.llim) + llim;
                    delta[i] = Math.exp(logdelta[i]);
                }
                //m <- length(logdelta)
                int m = logdelta.length; // the grids number
                //Lambdas <- matrix(eig.R$values,n-q,m) + matrix(delta,n-q,m,byrow=TRUE)
                double[][] lambdas = new double[n - q][m];
                for (int i = 0; i < n - q; i++) {
                    for (int j = 0; j < m; j++) lambdas[i][j] = ml_eig_R.eigenvalues[i] + delta[j];
                }
//				//Xis <- matrix(eig.L$values,n,m) + matrix(delta,n,m,byrow=TRUE)
//				double[][] Xis= new double[n][m];
//				for(int i=0;i<n;i++){
//					for(int j=0;j<m;j++)Xis[i][j]=eig_L.eigenvalues[i]+delta[j];
//				}
                //Etasq <- matrix(etas*etas,n-q,m)
                double[][] etasq = new double[n - q][m];
                for (int i = 0; i < n - q; i++) {
                    for (int j = 0; j < m; j++) etasq[i][j] = etas[i] * etas[i];
                }
//				LL <- 0.5*((n-q)*(log((n-q)/(2*pi))-1-log(colSums(Etasq/Lambdas)))-colSums(log(Lambdas)))
                double[] LL = new double[m];
                for (int i = 0; i < m; i++) {
                    double sum1 = 0, sum2 = 0;
                    for (int j = 0; j < n - q; j++) {
                        sum1 += (etasq[j][i] / lambdas[j][i]);
                        sum2 += (Math.log(lambdas[j][i]));
                    }
                    LL[i] = 0.5 * ((n - q) * (Math.log((n - q) / (2 * Math.PI)) - 1 - Math.log(sum1)) - sum2);
                }
                // dLL <- 0.5*delta*((n-q)*colSums(Etasq/(Lambdas*Lambdas))/colSums(Etasq/Lambdas)-colSums(1/Lambdas))
                double[] dLL = new double[m];
                for (int i = 0; i < m; i++) {
                    double sum1 = 0, sum2 = 0, sum3 = 0;
                    for (int j = 0; j < n - q; j++) {
                        sum1 += (etasq[j][i] / (lambdas[j][i] * lambdas[j][i]));
                        sum2 += (etasq[j][i] / lambdas[j][i]);
                        sum3 += (1 / lambdas[j][i]);
                    }
                    dLL[i] = 0.5 * delta[i] * ((n - q) * sum1 / sum2 - sum3);
                }
                if (dLL[0] < this.esp) {
                    optlogdelta.add((double) llim); //optlogdelta <- append(optlogdelta, llim)
                    optLL.add(this.delta_REML_LL_wo_Z(this.llim, ml_eig_R.eigenvalues, etas));//optLL <- append(optLL, emma.delta.ML.LL.wo.Z(llim,eig.R$values,etas,eig.L$values))
                }
                if (dLL[m - 1] > 0 - this.esp) {
                    optlogdelta.add((double) this.ulim);
                    optLL.add(delta_REML_LL_wo_Z(this.ulim, ml_eig_R.eigenvalues, etas));
                }
                //for( i in 1:(m-1) )
                for (int i = 0; i < m - 1; i++) {
                    if ((dLL[i] * dLL[i + 1] < 0 - this.esp * this.esp) && (dLL[i] > 0) && (dLL[i + 1] < 0)) {
                        RealRoot rlrt = new RealRoot();
                        RootSearchEmmadeltaREMLdLLwoZ func = new RootSearchEmmadeltaREMLdLLwoZ();
                        func.set_delta_REML_dLL_wo_Z(ml_eig_R.eigenvalues, etas);
                        double root = rlrt.brent(func, logdelta[i], logdelta[i + 1]);
                        //optlogdelta <- append(optlogdelta, r$root)
                        optlogdelta.add(root);
                        //optLL <- append(optLL, emma.delta.ML.LL.wo.Z(r$root,eig.R$values, etas, eig.L$values))
                        optLL.add(delta_REML_LL_wo_Z(root, ml_eig_R.eigenvalues, etas));
                    }
                }
            } else {
                System.out.println("non trivial Z part of java-EMMA has not been implemented yet");
                return;
            }
            int max_index = 0;
//			while(Double.isNaN(optLL.get(max_index)))max_index++;
            for (int i = 1; i < optlogdelta.size(); i++) {
                if (optLL.get(i) > optLL.get(max_index)) max_index = i;
            }
            double maxdelta = Math.exp(optlogdelta.get(max_index)); //maxdelta <- exp(optlogdelta[which.max(optLL)])
            double maxLL = optLL.get(max_index);//maxLL <- max(optLL)
            double maxva = 0;
            if (Z == null) {
                for (int i = 0; i < etas.length; i++) maxva += etas[i] * etas[i] / (ml_eig_R.eigenvalues[i] + maxdelta);
                maxva = maxva / (n - q);
            }

            double maxve = maxva * maxdelta;//maxve <- maxva*maxdelta
            this.remle_REML = maxLL;
            this.remle_delta = maxdelta;
            this.remle_ve = maxve;
            this.remle_vg = maxva;
            this.heritability = maxva / (maxva + maxve);
            if (!Double.isNaN(this.remle_REML)) {
                this.remle_calculated = true;
            }
//			System.out.println("Sample size: "+this.sample_size);
//			System.out.println("remle_REML="+maxLL);
//			System.out.println("remle_delta="+maxdelta);
//			System.out.println("remle_ve="+maxve);
//			System.out.println("remle_vg="+maxva);
//			System.out.println("remle_pseudo-heritability="+maxva/(maxva+maxve));
        } catch (Exception e) {
            e.printStackTrace();
        }
        return;
    }

    public void single_ML_LRT(double[] y, double[][] X, double[][] Z) throws IOException, InterruptedException {
        this.MLE(y, X, Z);
        double M1 = this.mle_ML;
        int[] indexes0 = {};
        int[] chr0 = {};
        this.ml_eig_R = null;//TODO!!
        this.MLE(y, cbind(chr0, indexes0), Z);
        double M0 = this.mle_ML;
        double stat = 2 * (M1 - M0);
        System.out.println(Test.chi2pr(stat, 1));
    }

    /*
     * compose an X array by the indexes chosen and a column of 1 ahead:
     * It is the same to R commend cbind(1, emmadat$xs[i1,],...,emmadat$xs[ik,]). e.g. X = cbind(1,emmadat$xs[1,],emmadat$xs[3,])
     */
    public double[][] cbind(int[] chrs, int[] indexes) {
        double[][] result_X = new double[this.sample_size][indexes.length + 1];
        for (int i = 0; i < this.sample_size; i++) {
            result_X[i][0] = 1;
        }
        for (int j = 0; j < indexes.length; j++) {
            double[] data_of_this_site = this.genotype.load_one_variant_by_index(chrs[j], indexes[j]);
            for (int i = 0; i < this.sample_size; i++) {
                result_X[i][j + 1] = data_of_this_site[this.indexes_with_phenotype_in_genotype[i]];
            }
        }
        return result_X;
    }

    /*
     * return the array = diag(1/sqrt(eig_L.eigenvalues[i]+optdelta))*[eig_L.vectors]
     */
    public double[][] reml_decompositioned_array() throws IOException, InterruptedException {
        if (this.reml_eig_L == null) {
            this.reml_eig_L = this.eigen_L_wo_Z();
        }
        double[][] result = new double[this.sample_size][this.sample_size];
        double[] diag = new double[this.sample_size];
        for (int i = 0; i < this.sample_size; i++) {
            //		if(this.reml_eig_L.eigenvalues[i]>0)
            diag[i] = 1.0 / Math.sqrt(this.reml_eig_L.eigenvalues[i] + this.remle_delta);
        }
        for (int i = 0; i < this.sample_size; i++) {
            for (int j = 0; j < this.sample_size; j++) {
                result[i][j] = this.reml_eig_L.eigenvectors_asrows[i][j] * diag[i]; // TODO: Check! [i][j] or [j][i]!
            }
        }
//		Test.write2file(result, "/Users/quan.long/look/de.new");
        return result;
    }


}

