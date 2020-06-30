package mixedmodel;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

import myMathLib.StatFuncs;

import org.apache.commons.math3.exception.MathIllegalArgumentException;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularMatrixException;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.stat.correlation.SpearmansCorrelation;
import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;

public class MultiPhenotype {

    public int num_of_pheno;
    public int sample_size_all_pheno;
    //	String[] pheno_names;
    public Phenotype[] phenotypes;
    public String[] sample_ids_all_pheno;
    public HashMap<String, Integer> phe_id2index;
    public HashMap<String, Integer> sample_id2index;

    public MultiPhenotype(ArrayList<Phenotype> phenos) {
        this.num_of_pheno = phenos.size();
        HashSet<String> all_ids = new HashSet<>();
        for (int phe_index = 0; phe_index < phenos.size(); phe_index++) {
            Phenotype the_pheno = phenos.get(phe_index);
            for (int k = 0; k < the_pheno.sample_ids.length; k++) all_ids.add(the_pheno.sample_ids[k]);
        }
        this.sample_size_all_pheno = all_ids.size();
        this.sample_ids_all_pheno = new String[sample_size_all_pheno];
        this.phenotypes = new Phenotype[num_of_pheno];
        Object[] all_ids_collection = all_ids.toArray();
        for (int k = 0; k < sample_size_all_pheno; k++) {
            this.sample_ids_all_pheno[k] = (String) all_ids_collection[k];

        }
        for (int phe_index = 0; phe_index < phenos.size(); phe_index++) {
            Phenotype the_pheno = phenos.get(phe_index);
            double[] values_with_NaN = new double[this.sample_size_all_pheno];
            the_pheno.setup_id2index();
            for (int k = 0; k < sample_size_all_pheno; k++) {
                if (the_pheno.sample_id2index.containsKey(this.sample_ids_all_pheno[k])) {
                    values_with_NaN[k] = the_pheno.values[the_pheno.sample_id2index.get(this.sample_ids_all_pheno[k])];
                } else {
                    values_with_NaN[k] = Double.NaN;
                }
            }
            Phenotype new_pheno = new Phenotype(the_pheno.phe_id, this.sample_ids_all_pheno, values_with_NaN, true);
            this.phenotypes[phe_index] = new_pheno;
        }
        this.setup_id2index();

    }

    // CHECKED -PATHUM
    public MultiPhenotype(String phenotype_file) {
        try {
            BufferedReader br = new BufferedReader(new FileReader(phenotype_file));
            String line = br.readLine(); // the header.
            String[] temp = line.split("\t");
            this.num_of_pheno = temp.length - 1;
            String[] pheno_names = new String[this.num_of_pheno];
            this.phenotypes = new Phenotype[this.num_of_pheno];
            for (int i = 0; i < this.num_of_pheno; i++) {
                pheno_names[i] = temp[i + 1];
            }
            line = br.readLine();
            while (line != null) {
                this.sample_size_all_pheno++;
                line = br.readLine();
            }
            this.sample_ids_all_pheno = new String[this.sample_size_all_pheno];
            double[][] full_phenotype = new double[num_of_pheno][this.sample_size_all_pheno];
            br = new BufferedReader(new FileReader(phenotype_file));
            line = br.readLine();
            line = br.readLine();
            int sample_index = 0;
            while (line != null) {
                temp = line.split("\t");
                this.sample_ids_all_pheno[sample_index] = temp[0];
                for (int i = 0; i < this.num_of_pheno; i++) {
                    if (temp[i + 1].equals("NA") || temp[i + 1].equals("NaN") || temp[i + 1].equals("*") || temp[i + 1].equals("N") || temp[i + 1].equals("?")) {
                        full_phenotype[i][sample_index] = Double.NaN;
                    } else {
                        full_phenotype[i][sample_index] = Double.parseDouble(temp[i + 1]);
                    }
                }
                sample_index++;
                line = br.readLine();
            }
            for (int k = 0; k < this.num_of_pheno; k++) {
                this.phenotypes[k] = new Phenotype(pheno_names[k], sample_ids_all_pheno, full_phenotype[k]);
            }
            this.setup_id2index();
            System.out.println("Finished reading " + phenotype_file);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public MultiPhenotype(String phenotype_file, String sep) {
        try {
            BufferedReader br = new BufferedReader(new FileReader(phenotype_file));
            String line = br.readLine(); // the header.
            String[] temp = line.split(sep);
            this.num_of_pheno = temp.length - 1;
            String[] pheno_names = new String[this.num_of_pheno];
            this.phenotypes = new Phenotype[this.num_of_pheno];
            for (int i = 0; i < this.num_of_pheno; i++) {
                pheno_names[i] = temp[i + 1];
            }
            line = br.readLine();
            while (line != null) {
                this.sample_size_all_pheno++;
                line = br.readLine();
            }
            this.sample_ids_all_pheno = new String[this.sample_size_all_pheno];
            double[][] full_phenotype = new double[num_of_pheno][this.sample_size_all_pheno];
            br = new BufferedReader(new FileReader(phenotype_file));
            line = br.readLine();
            line = br.readLine();
            int sample_index = 0;
            while (line != null) {
                temp = line.split(sep);
                this.sample_ids_all_pheno[sample_index] = temp[0];
                for (int i = 0; i < this.num_of_pheno; i++) {
                    if (temp[i + 1].equals("NA") || temp[i + 1].equals("NaN") || temp[i + 1].equals("*") || temp[i + 1].equals("N") || temp[i + 1].equals("?")) {
                        full_phenotype[i][sample_index] = Double.NaN;
                    } else {
                        full_phenotype[i][sample_index] = Double.parseDouble(temp[i + 1]);
                    }
                }
                sample_index++;
                line = br.readLine();
            }
            for (int k = 0; k < this.num_of_pheno; k++) {
                this.phenotypes[k] = new Phenotype(pheno_names[k], sample_ids_all_pheno, full_phenotype[k]);
            }
            this.setup_id2index();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    /*
     * in this constructor, the NaN will not be removed
     *
     * if the motivation is to do GWAS for a phenotype, there is NO case that NaN should be retained,
     * however, if we use the phenotype to calculate kinship, or using multiple phenotypes in the
     * same analysis, sometimes NaN should be retained to keep the sample size the same across all
     * phenotypes
     */

    public MultiPhenotype(String phenotype_file, boolean not_remove_NaN) {
        try {
            BufferedReader br = new BufferedReader(new FileReader(phenotype_file));
            String line = br.readLine(); // the header.
            String[] temp = line.split("\t");
            this.num_of_pheno = temp.length - 1;
            String[] pheno_names = new String[this.num_of_pheno];
            this.phenotypes = new Phenotype[this.num_of_pheno];
            for (int i = 0; i < this.num_of_pheno; i++) {
                pheno_names[i] = temp[i + 1];
            }
            line = br.readLine();
            while (line != null) {
                this.sample_size_all_pheno++;
                line = br.readLine();
            }
            this.sample_ids_all_pheno = new String[this.sample_size_all_pheno];
            double[][] full_phenotype = new double[num_of_pheno][this.sample_size_all_pheno];
            br = new BufferedReader(new FileReader(phenotype_file));
            line = br.readLine();
            line = br.readLine();
            int sample_index = 0;
            while (line != null) {
                temp = line.split("\t");
                this.sample_ids_all_pheno[sample_index] = temp[0];
                for (int i = 0; i < this.num_of_pheno; i++) {
                    if (temp[i + 1].equals("NA") || temp[i + 1].equals("NaN") || temp[i + 1].equals("*") || temp[i + 1].equals("N") || temp[i + 1].equals("?")) {
                        full_phenotype[i][sample_index] = Double.NaN;
                    } else {
                        full_phenotype[i][sample_index] = Double.parseDouble(temp[i + 1]);
                    }
                }
                sample_index++;
                line = br.readLine();
            }
            for (int k = 0; k < this.num_of_pheno; k++) {
                this.phenotypes[k] = new Phenotype(pheno_names[k], sample_ids_all_pheno, full_phenotype[k], not_remove_NaN);
            }
            this.setup_id2index();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    /*
     * phenotype-ids must be a subset of kinship-ids
     */
    public static double[][] match(KinshipMatrix kinship, Phenotype phenotype) {
        int sample_size = phenotype.sample_ids.length;
        int[] id_indexes = new int[sample_size];
        double[][] matched = new double[sample_size][sample_size];
        for (int i = 0; i < sample_size; i++) {
            if (kinship.ids2indexes.containsKey(phenotype.sample_ids[i])) {
                id_indexes[i] = kinship.ids2indexes.get(phenotype.sample_ids[i]);
            } else {
                System.out.println("!kinship.ids2indexes.containsKey(phenotype.sample_ids[i])");
                return null;
            }
        }
        for (int i = 0; i < sample_size; i++) {
            for (int j = i; j < sample_size; j++) {
                matched[i][j] = kinship.kinship_matrix.getEntry(id_indexes[i], id_indexes[j]);
                matched[j][i] = matched[i][j];
            }
        }
        //return KinshipMatrix.re_scale_kinship_matrix(matched);
        return matched;
    }

    public static double[] correlation(Phenotype p1, Phenotype p2) {
        ArrayList<String> shared_ids = new ArrayList<>();
        p1.setup_id2index();
        p2.setup_id2index();
        for (int i = 0; i < p1.sample_ids.length; i++) {
            if (p2.sample_id2index.containsKey(p1.sample_ids[i]))
                shared_ids.add(p1.sample_ids[i]);
        }
        if (shared_ids.size() <= 2) return null;
        double[] v1 = new double[shared_ids.size()];
        double[] v2 = new double[shared_ids.size()];
        for (int i = 0; i < shared_ids.size(); i++) {
            v1[i] = p1.values[p1.sample_id2index.get(shared_ids.get(i))];
            v2[i] = p2.values[p2.sample_id2index.get(shared_ids.get(i))];
        }
        double[] corr = new double[5];
        PearsonsCorrelation c0 = new PearsonsCorrelation();
        corr[0] = c0.correlation(v1, v2);
        corr[1] = StatFuncs.pvalue_spearman(Math.abs(corr[0]), shared_ids.size());
        SpearmansCorrelation c1 = new SpearmansCorrelation();
        corr[2] = c1.correlation(v1, v2);
        corr[3] = StatFuncs.pvalue_spearman(Math.abs(corr[2]), shared_ids.size());
        corr[4] = shared_ids.size();
//		for(int k=0;k<shared_ids.size();k++){
//			System.out.println(shared_ids.get(k)+"\t"+v1[k]+"\t"+v2[k]);
//		}
        return corr;
    }

    public void normalize_mean_var() {
        for (int phe_index = 0; phe_index < this.num_of_pheno; phe_index++) {
            this.phenotypes[phe_index].values = myMathLib.StatFuncs.standardize(this.phenotypes[phe_index].values);
        }
    }

    public void setup_id2index() {
        this.phe_id2index = new HashMap<>();
        for (int k = 0; k < this.num_of_pheno; k++) {
            this.phe_id2index.put(this.phenotypes[k].phe_id, k);
        }
        this.sample_id2index = new HashMap<>();
        for (int k = 0; k < this.sample_ids_all_pheno.length; k++) {
            this.sample_id2index.put(this.sample_ids_all_pheno[k], k);
        }
    }

    public void output_phenotypes_of_ids(ArrayList<String> ids_tobe_selected, String output_file) {
        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(output_file));
            bw.write("Id");
            for (int phe_index = 0; phe_index < this.num_of_pheno; phe_index++) {
                bw.write("\t" + this.phenotypes[phe_index].phe_id);
            }
            bw.write("\n");
            for (int k = 0; k < ids_tobe_selected.size(); k++) {
                int id_index = this.sample_id2index.get(ids_tobe_selected.get(k));
                bw.write(ids_tobe_selected.get(k));
                for (int phe_index = 0; phe_index < this.num_of_pheno; phe_index++) {
                    bw.write("\t" + this.phenotypes[phe_index].values[id_index]);
                }
                bw.write("\n");
            }
            bw.close();

            String folder = "/xxx/yyy/";
            File the_folder = new File(folder);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    /*
     * This function can handle object with NaN removed, therefore NOT require the same sample size
     * across all phenotype.
     */
    public void write2file(String output_file) {
        for (int phe_index = 1; phe_index < this.num_of_pheno; phe_index++) {
            if (this.phenotypes[phe_index].sample_id2index == null)
                this.phenotypes[phe_index].setup_id2index();
        }
        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(output_file));
            bw.write("Id");
            for (int phe_index = 0; phe_index < this.num_of_pheno; phe_index++) {
                bw.write("\t" + this.phenotypes[phe_index].phe_id);
            }
            bw.write("\n");
            for (int k = 0; k < this.sample_size_all_pheno; k++) {
                bw.write(this.sample_ids_all_pheno[k]);
                for (int phe_index = 0; phe_index < this.num_of_pheno; phe_index++) {
                    if (this.phenotypes[phe_index].sample_id2index.containsKey(this.sample_ids_all_pheno[k])) {
                        int index_k = this.phenotypes[phe_index].sample_id2index.get(this.sample_ids_all_pheno[k]);
                        bw.write("\t" + this.phenotypes[phe_index].values[index_k]);
                    } else {
                        bw.write("\tNA");
                    }
                }
                bw.write("\n");
            }
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    /*
     * Move the variances explained by a kinship matrix from the phenotype
     * It is assumed that all the phenotypes in the
     */
    public void convert_by_matrix_FaSTLMM(KinshipMatrix kinship) {
        this.normalize_mean_var();
        for (int phe_index = 0; phe_index < this.num_of_pheno; phe_index++) {
            double[][] matched_kinship = match(kinship, this.phenotypes[phe_index]);
            RealMatrix trans = KinshipMatrix.find_transformation_matrix_fastlmm(matched_kinship, this.phenotypes[phe_index].values);
            this.phenotypes[phe_index].values = trans.operate(this.phenotypes[phe_index].values);
        }
    }

    public void convert_by_matrix_EMMA(KinshipMatrix kinship, VariantsDouble genotype) {
        for (int phe_index = 0; phe_index < this.num_of_pheno; phe_index++) {
            if (phe_index % 1000 == 0) System.out.println("Finished " + phe_index / 1000 + "K");
            EMMA emma = new EMMA(this.phenotypes[phe_index], genotype, kinship.kinship_matrix.getData());//genotype.read_in_kinship(kinship_file));
            emma.REMLE_null();
            System.out.println(emma.remle_delta + "\t" + emma.remle_ve / emma.remle_vg);
            //double[][] matched_kinship=match(kinship, this.phenotypes[phe_index]);
            double[][] decomposed_array = emma.reml_decompositioned_array();
            emma.phenotype.generate_new_Y_by_multiplying(decomposed_array);
            this.phenotypes[phe_index].values = emma.phenotype.new_Y_4emmax;
        }
    }

    public double[] adjust_cofactors(String tobe_adj, String[] adj) {
        if (!this.phe_id2index.containsKey(tobe_adj)) {
            System.out.println(tobe_adj + "NOT exists!");
            return null;
        }
        int index_tobe = this.phe_id2index.get(tobe_adj);
        int[] indexes_adj = new int[adj.length];
        for (int k = 0; k < adj.length; k++) {
            if (!this.phe_id2index.containsKey(adj[k])) {
                System.out.println(adj[k] + "NOT exists!");
                return null;
            } else {
                indexes_adj[k] = this.phe_id2index.get(adj[k]);
            }
        }
        return adjust_cofactors(index_tobe, indexes_adj);
    }

    /*
     * before running this, one has to make sure that the NaNs are removed.
     */
    public double[] adjust_cofactors(int index_tobe, int[] indexes_adj) {
        if (this.phenotypes[index_tobe].sample_id2index == null) {
            this.phenotypes[index_tobe].setup_id2index();
        }
        double[] means = new double[indexes_adj.length];
        for (int i = 0; i < indexes_adj.length; i++) {
            means[i] = this.phenotypes[indexes_adj[i]].sample_mean();
            if (this.phenotypes[indexes_adj[i]].sample_id2index == null) {
                this.phenotypes[indexes_adj[i]].setup_id2index();
            }
        }
        double[] Y = this.phenotypes[index_tobe].values.clone();
        double[][] Xs = new double[Y.length][indexes_adj.length];
        for (int sample_index = 0; sample_index < Y.length; sample_index++) {
            String id = this.phenotypes[index_tobe].sample_ids[sample_index];
            for (int i = 0; i < indexes_adj.length; i++) {
                Phenotype the_cofactor = this.phenotypes[indexes_adj[i]];
                if (the_cofactor.sample_id2index.containsKey(id)) {
                    Xs[sample_index][i] = the_cofactor.values[the_cofactor.sample_id2index.get(id)];
                } else Xs[sample_index][i] = means[i];
            }
        }

        OLSMultipleLinearRegression reg1 = new OLSMultipleLinearRegression();
        reg1.newSampleData(Y, Xs);
        double[] residuals = reg1.estimateResiduals();
        double[] betas = reg1.estimateRegressionParameters();
        String[] new_sample_ids = this.phenotypes[index_tobe].sample_ids.clone();
        String new_phe_id = this.phenotypes[index_tobe].phe_id + ".ADJ";
        for (int i = 0; i < indexes_adj.length; i++) {
            new_phe_id = new_phe_id + "." + this.phenotypes[indexes_adj[i]].phe_id;
        }
        Phenotype adjusted = new Phenotype(new_phe_id, new_sample_ids, residuals);
        this.add_adjusted_pheno(adjusted);
        return betas;
    }

    /*
     * Add an new phenotype generated by adjusting existing phenotype
     */
    public void add_adjusted_pheno(Phenotype adj) {
        this.num_of_pheno++;
        Phenotype[] the_new = new Phenotype[this.num_of_pheno];
        for (int i = 0; i < this.num_of_pheno - 1; i++) the_new[i] = this.phenotypes[i];
        the_new[this.num_of_pheno - 1] = adj;
        this.phenotypes = the_new;
        this.phe_id2index.put(adj.phe_id, this.num_of_pheno - 1);
    }

    /*
     * regress Y to Xs on the available_IDs. Before entering this function, one has to ensure the availability of all IDs.
     */
    public void regression(int index_Y, int[] index_Xs, ArrayList<String> available_IDs, int sample_size_cutoff, BufferedWriter bw) {
        try {
            double[] Y = new double[available_IDs.size()];
            double[][] Xs = new double[available_IDs.size()][index_Xs.length];
            for (int i = 0; i < available_IDs.size(); i++) {
                String key = available_IDs.get(i);
                Y[i] = this.phenotypes[index_Y].values[this.phenotypes[index_Y].sample_id2index.get(key)];
                for (int k = 0; k < index_Xs.length; k++) {
                    Xs[i][k] = this.phenotypes[index_Xs[k]].values[this.phenotypes[index_Xs[k]].sample_id2index.get(key)];
                }
            }
            OLSMultipleLinearRegression reg1 = new OLSMultipleLinearRegression();
            try {
                reg1.newSampleData(Y, Xs);
                double[] beta = reg1.estimateRegressionParameters();
                double R2 = reg1.calculateAdjustedRSquared();
                bw.write("Sample_size=" + Y.length + "\n");
                if (Y.length < sample_size_cutoff) {
                    bw.write("Sample size smaller than " + sample_size_cutoff + ". Analysis not done.");
                    return;
                }
                bw.write("Regression coefficients: ");
                for (int k = 0; k < beta.length; k++) bw.write(beta[k] + ", ");
                bw.write("\nR2=" + R2 + "\n");
            } catch (SingularMatrixException e) {
                System.out.println("SingularMatrixException. Continued to the next regression.");
            } catch (MathIllegalArgumentException e2) {
                System.out.println("MathIllegalArgumentException. Continued to the next regression.");
            }

        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    /*
     * regress Y to Xs conditional on C
     */
    public void regression(int index_Y, int[] index_Xs, int index_C, int sample_szie_cutoff, String report_file) {
        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(report_file));
            bw.write("Naive (multiple-)regression analysis:\n" + this.phenotypes[index_Y].phe_id + "~");
            for (int k = 0; k < index_Xs.length - 1; k++) bw.write(this.phenotypes[index_Xs[k]].phe_id + "+");
            bw.write(this.phenotypes[index_Xs[index_Xs.length - 1]].phe_id + "\n");
            bw.write("Conditional on " + this.phenotypes[index_C].phe_id + "\n");

            HashMap<Integer, ArrayList<String>> C_avai_IDs = new HashMap<>();
            double[] C_values = this.phenotypes[index_C].values;
            //gather common IDs
            String[] total_IDs = this.phenotypes[index_Y].sample_ids;
            ArrayList<String> available_IDs = new ArrayList<>();
            for (int i = 0; i < total_IDs.length; i++) {
                boolean avai = true;
                for (int k = 0; k < index_Xs.length; k++) {
                    if (!this.phenotypes[index_Xs[k]].sample_id2index.containsKey(total_IDs[i])) {
                        avai = false;
                        break;
                    }
                }
                if (avai) {
                    available_IDs.add(total_IDs[i]);
                    if (this.phenotypes[index_C].sample_id2index.containsKey(total_IDs[i])) { // this ID is in C_index, so lets find out which set
                        int the_C_value = (int) (C_values[this.phenotypes[index_C].sample_id2index.get(total_IDs[i])] + 0.01); //ensure 0.999999997==1
                        if (C_avai_IDs.containsKey(the_C_value)) C_avai_IDs.get(the_C_value).add(total_IDs[i]);
                        else {
                            ArrayList<String> the_set = new ArrayList<>();
                            the_set.add(total_IDs[i]);
                            C_avai_IDs.put(the_C_value, the_set);
                        }
                    }
                }
            }
            bw.write("Combined (without stratified by the co-factor):\n");
            regression(index_Y, index_Xs, available_IDs, sample_szie_cutoff, bw);
            bw.write("\n");
            for (int the_C_value : C_avai_IDs.keySet()) {
                bw.write(this.phenotypes[index_C].phe_id + "=" + the_C_value + ":\n");
                regression(index_Y, index_Xs, C_avai_IDs.get(the_C_value), sample_szie_cutoff, bw);
                bw.write("\n");
            }
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}





























