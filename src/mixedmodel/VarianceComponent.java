package mixedmodel;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;
import org.apache.commons.math3.optimization.GoalType;
import org.apache.commons.math3.optimization.univariate.BrentOptimizer;
import org.apache.commons.math3.optimization.univariate.UnivariatePointValuePair;
import org.jfree.chart.plot.MultiplePiePlot;

public class VarianceComponent {

    public KinshipMatrix kinship1;
    public KinshipMatrix kinship2;
    public MultiPhenotype phenotypes;


    public static void match_kinships_and_phenotype(String genotype4kinship1_file, String genotype4kinship2_file,
                                                    String phenotype, String kinship1_file_new, String kinship2_file_new, String phenotype_new,
                                                    String method1, String method2) {
        VariantsDouble data1 = new VariantsDouble(genotype4kinship1_file);
        VariantsDouble data2 = new VariantsDouble(genotype4kinship2_file);
        MultiPhenotype pheno = new MultiPhenotype(phenotype);
        System.out.println("Finished reading files.");
        for (int k = 0; k < pheno.num_of_pheno; k++) {
            if (pheno.phenotypes[k].sample_ids.length != pheno.phenotypes[0].sample_ids.length)
                System.out.println("Phenotype sample size is not the same");
        }
        Phenotype the_pheno = pheno.phenotypes[0];
        HashMap<String, Integer> pheno_ids = new HashMap<>();
        for (int i = 0; i < the_pheno.sample_ids.length; i++) {
            pheno_ids.put(the_pheno.sample_ids[i], i);
        }
        ArrayList<String> common_ids = new ArrayList<>();
        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(phenotype_new));
            bw.write("ID");
            for (int k = 0; k < pheno.num_of_pheno; k++) {
                bw.write("\t" + pheno.phenotypes[k].phe_id);
            }
            bw.write("\n");
            for (int i = 0; i < the_pheno.sample_ids.length; i++) {
                if (data1.sample_id2index.containsKey(the_pheno.sample_ids[i]) &&
                        data2.sample_id2index.containsKey(the_pheno.sample_ids[i])) {
                    common_ids.add(the_pheno.sample_ids[i]);
                    bw.write(the_pheno.sample_ids[i]);
                    for (int k = 0; k < pheno.num_of_pheno; k++) {
                        bw.write("\t" + pheno.phenotypes[k].values[i]);
                    }
                    bw.write("\n");
                } else {
                    System.out.println("Phenotype sample ID=" + the_pheno.sample_ids[i] + " not in kinship ids");
                }
            }
            bw.close();
            System.out.println("Finished writing pheno file.");

            System.out.println("Started generating kinship files.");
            if (method1.equals("RRM")) {
                data1.calculate_WG_RRM_kinship(kinship1_file_new, common_ids);
            } else if (method1.equals("IBS")) {
                data1.calculate_WG_IBS_kinship(kinship1_file_new + ".raw", common_ids);
                KinshipMatrix.re_scale_kinship_matrix(kinship1_file_new + ".raw", kinship1_file_new);
            } else {
                System.out.println(method1 + " is not a correct type of kinship supported by JAWAMix5");
            }

            if (method2.equals("RRM")) {
                data2.calculate_WG_RRM_kinship(kinship2_file_new, common_ids);
            } else if (method2.equals("IBS")) {
                data2.calculate_WG_IBS_kinship(kinship2_file_new + ".raw", common_ids);
                KinshipMatrix.re_scale_kinship_matrix(kinship2_file_new + ".raw", kinship2_file_new);
            } else {
                System.out.println(method1 + " is not a correct type of kinship supported by JAWAMix5");
            }
        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    /*
     * It is assumed that the dimension of kinship and all phenotypes are the same and in the same order
     * In general, constructor will guarantee that
     *
     * FastLMM used
     */
    public static void VO_sinle_kinship_fastlmm(KinshipMatrix kinship1, KinshipMatrix kinship2, MultiPhenotype phenotypes,
                                                String output) {
        try {
            System.out.println("Started solving Eigen values");
            int sample_size = kinship1.ids.length;

            EigenDecomposition eigen1 = new EigenDecomposition(new Array2DRowRealMatrix(kinship1.kinship_matrix.getData()), 0);
            RealMatrix global_Ut1 = eigen1.getVT();
            double[] S1 = eigen1.getRealEigenvalues();

            EigenDecomposition eigen2 = new EigenDecomposition(new Array2DRowRealMatrix(kinship2.kinship_matrix.getData()), 0);
            RealMatrix global_Ut2 = eigen2.getVT();
            double[] S2 = eigen2.getRealEigenvalues();

            BufferedWriter bw = new BufferedWriter(new FileWriter(output));
            bw.write("Gene,Cell-Env,Genotype\n");
            System.out.println("Started estimating VO");
            final double[] all_one = new double[sample_size];
            Arrays.fill(all_one, 1);
            for (int pheno_index = 0; pheno_index < phenotypes.num_of_pheno; pheno_index++) {
                double[] normalized_Y = myMathLib.StatFuncs.standardize(phenotypes.phenotypes[pheno_index].values);

                double[] Y_transformed1 = global_Ut1.operate(normalized_Y);
                double[] intercept1 = global_Ut1.operate(all_one);
                VarComp_Result vo1 = fullRankSolver_null_ML(Y_transformed1, intercept1, S1);

                double[] Y_transformed2 = global_Ut2.operate(normalized_Y);
                double[] intercept2 = global_Ut2.operate(all_one);
                VarComp_Result vo2 = fullRankSolver_null_ML(Y_transformed2, intercept2, S2);

                bw.write(phenotypes.phenotypes[pheno_index].phe_id + "," + vo1.sigma_g / (vo1.sigma_e + vo1.sigma_g) + "," +
                        vo2.sigma_g / (vo2.sigma_e + vo2.sigma_g) + "\n");
                System.out.println(phenotypes.phenotypes[pheno_index].phe_id + "," + vo1.sigma_g / (vo1.sigma_e + vo1.sigma_g) + "," +
                        vo2.sigma_g / (vo2.sigma_e + vo2.sigma_g));
            }
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public static void heritability_fullrank(KinshipMatrix global_kinship, MultiPhenotype phenotypes,
                                             String output) {
        try {
            System.out.println("Started solving Eigen values");
            EigenDecomposition eigen = new EigenDecomposition(new Array2DRowRealMatrix(global_kinship.kinship_matrix.getData()), 0);
            RealMatrix global_Ut = eigen.getVT();
            double[] S = eigen.getRealEigenvalues();
            int sample_size = global_kinship.ids.length;

            BufferedWriter bw = new BufferedWriter(new FileWriter(output));
            bw.write("Phenotype_ID,Pseudo-Heritability\n");
            System.out.println("Started estimating Variance Component.");
            final double[] all_one = new double[sample_size];
            Arrays.fill(all_one, 1);
            for (int pheno_index = 0; pheno_index < phenotypes.num_of_pheno; pheno_index++) {
                double[] normalized_Y = myMathLib.StatFuncs.standardize(phenotypes.phenotypes[pheno_index].values);
                double[] Y_transformed = global_Ut.operate(normalized_Y);
                double[] intercept = global_Ut.operate(all_one);
                VarComp_Result vo2 = fullRankSolver_null_ML(Y_transformed, intercept, S);
                double h2 = vo2.sigma_g / (vo2.sigma_e + vo2.sigma_g);
                double fisher_info = FaSTLMM_FullRank.fisher_info_nobeta(Y_transformed, S, h2);
                bw.write(phenotypes.phenotypes[pheno_index].phe_id + "," + h2 + "," + (1.0 / fisher_info) +
                        "," + (h2 - 1.96 / Math.sqrt(Math.abs(fisher_info * sample_size))) +
                        "," + (h2 + 1.96 / Math.sqrt(Math.abs(fisher_info * sample_size))) + "\n");
            }
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    /*
     * terms_for_kinship must in format of .tsv (i.e., phenotype format)
     * the ids of kinship terms must be the same to the ones in phenotypes.
     */
    public static void heritability_lowrank(String terms_for_kinship, MultiPhenotype phenotypes, String output) {
        try {
            MultiPhenotype terms = new MultiPhenotype(terms_for_kinship);
            int sample_size = phenotypes.sample_size_all_pheno;
            for (int k = 0; k < sample_size; k++) {
                if (!terms.sample_ids_all_pheno[k].equals(phenotypes.sample_ids_all_pheno[k])) {
                    System.out.println("Sample IDs of kinship and phenotypes not matching: " +
                            terms.sample_ids_all_pheno[k] + "/" + phenotypes.sample_ids_all_pheno[k]);
                    System.out.println("Exit!");
                    return;
                }
            }
            double[][] data_matching_phenotype = new double[terms.num_of_pheno][sample_size];
            for (int term_index = 0; term_index < terms.num_of_pheno; term_index++) {
                for (int k = 0; k < sample_size; k++) {
                    data_matching_phenotype[term_index][k] = terms.phenotypes[term_index].values[k];
                }
            }
            FaSTLMM.normalize_terms_withSqrtNumVarDivided(data_matching_phenotype);

            System.out.println("Started solving Singular Value Decomposition");

            SingularValueDecomposition local_svd = new SingularValueDecomposition(
                    (new Array2DRowRealMatrix(data_matching_phenotype)).transpose());
            double[] S = local_svd.getSingularValues();
            for (int k = 0; k < S.length; k++) S[k] = (S[k] * S[k]);
            RealMatrix local_U1t = local_svd.getUT();

            BufferedWriter bw = new BufferedWriter(new FileWriter(output));
            bw.write("Phenotype_ID,Pseudo-Heritability\n");
            System.out.println("Started estimating Variance Component.");
            final double[] all_one = new double[sample_size];
            Arrays.fill(all_one, 1);
            for (int pheno_index = 0; pheno_index < phenotypes.num_of_pheno; pheno_index++) {
                double[] Y_new_original = myMathLib.StatFuncs.standardize(phenotypes.phenotypes[pheno_index].values);
                double[] transformed_Y = local_U1t.operate(Y_new_original);
                RealMatrix local_IminusU1U1t = local_svd.getU().multiply(local_U1t).scalarMultiply(-1);
                for (int i = 0; i < sample_size; i++)
                    local_IminusU1U1t.setEntry(i, i, 1 + local_IminusU1U1t.getEntry(i, i));
                double[] Y_transformed2 = local_IminusU1U1t.operate(Y_new_original);
                double[][] X_new_original = new double[1][];
                X_new_original[0] = all_one;
                //double[][] transformed_X=new double[1][];
                //transformed_X[0]=local_U1t.operate(all_one);
                VarComp_Result result = FaSTLMM_LowRank.lowRankSolver_FixedEffects(X_new_original, Y_new_original,
                        transformed_Y, Y_transformed2, S, -Double.MAX_VALUE, local_U1t, local_IminusU1U1t);
                bw.write(phenotypes.phenotypes[pheno_index].phe_id + "," + result.sigma_g / (result.sigma_e + result.sigma_g) + "\n");

            }
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    /*
     * It is assumed that the dimension of kinship and all phenotypes are the same and in the same order
     * In general, constructor will guarantee that
     *
     * FastLMM used
     */
    public static void VO_sinle_kinship_emma(KinshipMatrix kinship1, KinshipMatrix kinship2, MultiPhenotype phenotypes,
                                             VariantsDouble genotype, String output) {
        try {
            System.out.println("Started solving Eigen values");
            int sample_size = kinship1.ids.length;

            BufferedWriter bw = new BufferedWriter(new FileWriter(output));
            bw.write("Gene,Cell-Env,Genotype\n");
            System.out.println("Started estimating VO");
            final double[] all_one = new double[sample_size];
            Arrays.fill(all_one, 1);
            for (int pheno_index = 0; pheno_index < phenotypes.num_of_pheno; pheno_index++) {
                EMMA emma1 = new EMMA(phenotypes.phenotypes[pheno_index], genotype, kinship1.kinship_matrix.getData());
                emma1.REMLE_null();
                EMMA emma2 = new EMMA(phenotypes.phenotypes[pheno_index], genotype, kinship2.kinship_matrix.getData());
                emma2.REMLE_null();

                bw.write(phenotypes.phenotypes[pheno_index].phe_id + "," + emma1.remle_vg / (emma1.remle_ve + emma1.remle_vg) + "," +
                        emma2.remle_vg / (emma2.remle_ve + emma2.remle_vg) + "\n");
                System.out.println(phenotypes.phenotypes[pheno_index].phe_id + "," + emma1.remle_vg / (emma1.remle_ve + emma1.remle_vg) + "," +
                        emma2.remle_vg / (emma2.remle_ve + emma2.remle_vg) + "\n");
            }
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    /*
     * It is assumed that the dimension of kinship and all phenotypes are the same and in the same order
     * In general constructor will guarantee that
     */
    public static void VO_double_kinship(KinshipMatrix kinship1, KinshipMatrix kinship2, MultiPhenotype phenotypes, String output) {
        try {
            System.out.println("Started solving Eigen values");
            int sample_size = kinship1.ids.length;

            EigenDecomposition eigen1 = new EigenDecomposition(new Array2DRowRealMatrix(kinship1.kinship_matrix.getData()), 0);
            RealMatrix global_Ut1 = eigen1.getVT();
            double[] S1 = eigen1.getRealEigenvalues();

            EigenDecomposition eigen2 = new EigenDecomposition(new Array2DRowRealMatrix(kinship2.kinship_matrix.getData()), 0);
            RealMatrix global_Ut2 = eigen2.getVT();
            double[] S2 = eigen2.getRealEigenvalues();

            BufferedWriter bw = new BufferedWriter(new FileWriter(output));
            bw.write("Gene,Cell-Env,Genotype,alpha\n");
            System.out.println("Started estimating VO");
            final double[] all_one = new double[sample_size];
            Arrays.fill(all_one, 1);
            for (int pheno_index = 0; pheno_index < phenotypes.num_of_pheno; pheno_index++) {
                double[] normalized_Y = myMathLib.StatFuncs.standardize(phenotypes.phenotypes[pheno_index].values);

                double[] Y_transformed1 = global_Ut1.operate(normalized_Y);
                double[] intercept1 = global_Ut1.operate(all_one);
                VarComp_Result vo1 = fullRankSolver_null_ML(Y_transformed1, intercept1, S1);

                double[] Y_transformed2 = global_Ut2.operate(normalized_Y);
                double[] intercept2 = global_Ut2.operate(all_one);
                VarComp_Result vo2 = fullRankSolver_null_ML(Y_transformed2, intercept2, S2);

                double best_ml = -Double.MAX_VALUE;
                double best_alpha = -1;
                for (double alpha = 0.1; alpha <= 0.9; alpha = alpha + 0.1) {
                    Array2DRowRealMatrix new_matrix = new Array2DRowRealMatrix(kinship1.kinship_matrix.scalarMultiply(alpha).
                            add(kinship2.kinship_matrix.scalarMultiply(1 - alpha)).getData());
                    EigenDecomposition eigen_mix = new EigenDecomposition(new_matrix, 0);
                    RealMatrix global_Ut_mix = eigen_mix.getVT();
                    double[] S_mix = eigen_mix.getRealEigenvalues();
                    double[] Y_transformed_mix = global_Ut_mix.operate(normalized_Y);
                    double[] intercept_mix = global_Ut_mix.operate(all_one);
                    VarComp_Result vo_mix = fullRankSolver_null_ML(Y_transformed_mix, intercept_mix, S_mix);
                    if (vo_mix.ml > best_ml) {
                        best_alpha = alpha;
                        best_ml = vo_mix.ml;
                    }
                }
                double last_best = best_alpha;
                for (double alpha = last_best - 0.09; alpha <= last_best + 0.09; alpha = alpha + 0.01) {
                    Array2DRowRealMatrix new_matrix = new Array2DRowRealMatrix(kinship1.kinship_matrix.scalarMultiply(alpha).
                            add(kinship2.kinship_matrix.scalarMultiply(1 - alpha)).getData());
                    EigenDecomposition eigen_mix = new EigenDecomposition(new_matrix, 0);
                    RealMatrix global_Ut_mix = eigen_mix.getVT();
                    double[] S_mix = eigen_mix.getRealEigenvalues();
                    double[] Y_transformed_mix = global_Ut_mix.operate(normalized_Y);
                    double[] intercept_mix = global_Ut_mix.operate(all_one);
                    VarComp_Result vo_mix = fullRankSolver_null_ML(Y_transformed_mix, intercept_mix, S_mix);
                    if (vo_mix.ml > best_ml) {
                        best_alpha = alpha;
                        best_ml = vo_mix.ml;
                    }
                }

                bw.write(phenotypes.phenotypes[pheno_index].phe_id + "," + vo1.sigma_g / (vo1.sigma_e + vo1.sigma_g) + "," +
                        vo2.sigma_g / (vo2.sigma_e + vo2.sigma_g) + "," + best_alpha + "\n");
                System.out.println(phenotypes.phenotypes[pheno_index].phe_id + "," + vo1.sigma_g / (vo1.sigma_e + vo1.sigma_g) + "," +
                        vo2.sigma_g / (vo2.sigma_e + vo2.sigma_g) + "," + best_alpha);
            }
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }


    }

    public static VarComp_Result fullRankSolver_null_ML(double[] Y_transformed, double[] intercept, double[] S) {
        double[][] transformed_X = new double[1][];
        transformed_X[0] = intercept;
        BrentOptimizer bo = new BrentOptimizer(FaSTLMM.brent_rel, FaSTLMM.brent_abs);
        double best_delta = -1;
        double best_ml = -(Double.MAX_VALUE);
        UnivariateFunction ml_function = new MLFullRank4BrentOptimal(transformed_X, Y_transformed, S);
        for (double min = FaSTLMM.lowlimit; min < FaSTLMM.uplimit; min += FaSTLMM.grid_step) {
            double min_delta = Math.exp(min);
            double max_delta = Math.exp(min + FaSTLMM.grid_step);
            UnivariatePointValuePair result = bo.optimize(FaSTLMM.maxEval, ml_function, GoalType.MAXIMIZE, min_delta, max_delta);
            double the_ml = result.getValue();
            if (the_ml > best_ml) {
                best_ml = the_ml;
                best_delta = result.getPoint();
            }
        }
        double[] beta_0 = MLFullRank4BrentOptimal.calculate_beta(Y_transformed, transformed_X, S, best_delta);
        double sigma_g_null = MLFullRank4BrentOptimal.sigma_g2_LL(Y_transformed, transformed_X, S, best_delta, beta_0);
        double sigma_e_null = sigma_g_null * best_delta;
        VarComp_Result result = new VarComp_Result();
        result.ml = best_ml;
        result.sigma_e = sigma_e_null;
        result.sigma_g = sigma_g_null;
        return result;
    }

    public static void temp_try_normal() {
        String working = "/Volumes/Projects/DATA/GETx/december_release/";
        String data_folder = working + "JAWAMix5_format/";
        String temp_rna_matrix = data_folder + "RPKM_GeneLevel_December.Thyroid.RNA.kinship";
        String temp_dna_matrix = data_folder + "RPKM_GeneLevel_December.Thyroid.DNA.kinship";
        String ori_pheno = data_folder + "RPKM_GeneLevel_December.Thyroid.tsv.dna-m.tsv";
        String temp_pheno = data_folder + "random_normal.tsv";
        String[] ids = new MultiPhenotype(ori_pheno).phenotypes[0].sample_ids;
        int num_phe = 100;
        NormalDistribution normal = new NormalDistribution();
        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(temp_pheno));
            bw.write("ID");
            for (int k = 0; k < num_phe; k++) bw.write("\tPhe" + k);
            bw.write("\n");
            for (int i = 0; i < ids.length; i++) {
                bw.write(ids[i]);
                for (int k = 0; k < num_phe; k++) {
                    bw.write("\t" + normal.sample());
                }
                bw.write("\n");
            }
            bw.close();
            VO_sinle_kinship_fastlmm(new KinshipMatrix(temp_rna_matrix), new KinshipMatrix(temp_dna_matrix),
                    new MultiPhenotype(temp_pheno), temp_pheno + ".vo.temp.csv");
        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    public static void main(String[] args) {
        String working = "/Volumes/Projects/DATA/GETx/december_release/";
        String data_folder = working + "JAWAMix5_format/";
        String dna_hdf5 = data_folder + "GTEx_5M_191_Dec2012.vcf.csv.hdf5";
        String rna_all_hdf5 = data_folder + "RPKM_GeneLevel_December.all.csv.hdf5";

        String exp_kinship = data_folder + "RPKM_GeneLevel_December.all.csv.hdf5.RRM.kinship";
        String dna_kinship_all_rescaled = data_folder + "GTEx_5M_191_Dec2012.vcf.csv.all.kinship.rescaled.ibs";
        String rna_all_phenotype_file = data_folder + "RPKM_GeneLevel_December.all.tsv";
        String dna_kinship_RRM = data_folder + "GTEx_5M_191_Dec2012.vcf.csv.RRMkinship.RRM";

        String dna_kinship_all_RRM = data_folder + "GTEx_5M_191_Dec2012.vcf.csv.all.RRM.kinship";

        String output_file = "/Volumes/Projects/GETx/vo/expression_VO.RRM.csv";

        temp_try_normal();
//		match_kinships_and_phenotype(dna_hdf5, rna_all_hdf5, rna_all_phenotype_file, 
//				dna_kinship_all_RRM, exp_kinship+".matched", rna_all_phenotype_file+".matched");
//		
//		VO_sinle_kinship(new KinshipMatrix(exp_kinship+".matched"), new KinshipMatrix(dna_kinship_all_RRM), 
//				new MultiPhenotype(rna_all_phenotype_file+".matched"), output_file);

    }

}
