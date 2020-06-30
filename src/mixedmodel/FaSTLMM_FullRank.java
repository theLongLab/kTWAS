package mixedmodel;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.Arrays;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.optimization.ConvergenceChecker;
import org.apache.commons.math3.optimization.GoalType;
import org.apache.commons.math3.optimization.univariate.BrentOptimizer;
import org.apache.commons.math3.optimization.univariate.UnivariatePointValuePair;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.util.MathArrays;

import ch.systemsx.cisd.base.mdarray.MDDoubleArray;
import ch.systemsx.cisd.hdf5.HDF5MDDataBlock;

/*
 * FaST LMM algorithm, Lippert C, Listgarten J, ... & Heckerman D. Nat Methods 2011  
 * 
 * essentially its results are the same to EMMA, but faster (in EMMAX time) 
 * and follows the formulation of FaST-LMM paper.	 
 */

public class FaSTLMM_FullRank extends FaSTLMM{
	
	final RealMatrix global_Ut; // eigen vectors 
	final double[] S; //eigen values
	final double[] Y_transformed;
	final double[] intercept;
	double Likelihood_null; // could be REML or ML, depending on the parameter passed to constructor.
	double sigma_g_null;
	double sigma_e_null;
	
	/* inherited:	
	public final static double machep=1e-10;//No need to use 2.220446049250313E-16; // Machine epsilon = Math.ulp(1d)
	public final static double brent_rel=Math.sqrt(machep);
	public final static double brent_abs=brent_rel;
	public final int maxEval=100;	
	public final static double uplimit=10;
	public final static double lowlimit=-10;
	public final static double grid_step=0.1;	
	
	String genotype_hdf5_file;
	VariantsDouble genotype;
	int[] indexes_with_phenotype_in_genotype;
	Phenotype phenotype; // all samples with phenotype will have genotype in the object: the constructor will guarantee that by checking the data.
	double[][] global_kinship_matrix; // this is corresponding to the phenotype, but generated from genotype matched kinship file.
	int sample_size;			
	*/
	
	/* 
	 * Before entering the constructor below, one has to guarantee the kinship file is full-ranked. 
	 * Otherwise use FaSTLMM_LowRank.java
	 */
	public FaSTLMM_FullRank(String genotype_file_hdf5, String global_kinship_file, Phenotype input_phenotype, boolean ml){
		super(input_phenotype, genotype_file_hdf5, global_kinship_file);
		EigenDecomposition eigen=new EigenDecomposition(new Array2DRowRealMatrix(this.global_kinship_matrix),0); 
		this.global_Ut=eigen.getVT();
		this.S=eigen.getRealEigenvalues();
//		double[][] havelook=this.global_Ut.getData();
		this.Y_transformed=this.global_Ut.operate(this.phenotype.values);
		double[] all_one=new double[this.sample_size];
		Arrays.fill(all_one, 1);
		this.intercept=this.global_Ut.operate(all_one);
		if(ml)this.fullRankSolver_null_ML();
		else this.fullRankSolver_null_REML();
		// set up null-REML in the model that only fixed effect is the transformed intercept. 
	}

	public void fullRankSolver_null_REML(){
		double[][] transformed_X=new double[1][];
		double[] intercept_1=new double[this.sample_size];
		Arrays.fill(intercept_1, 1);
		transformed_X[0]=this.global_Ut.operate(intercept_1);
		BrentOptimizer bo=new BrentOptimizer(FaSTLMM.brent_rel, FaSTLMM.brent_abs);
		double best_delta=-1;
		double best_reml=-(Double.MAX_VALUE);
		UnivariateFunction reml_function=new REMLFullRank4BrentOptimal(transformed_X,Y_transformed,S);
		for(double min=lowlimit;min<uplimit;min+=grid_step){
			double min_delta=Math.exp(min);
			double max_delta=Math.exp(min+grid_step);
			UnivariatePointValuePair result=bo.optimize(maxEval, reml_function, GoalType.MAXIMIZE, min_delta, max_delta);
			double the_reml=result.getValue();
			if(the_reml>best_reml){
				best_reml=the_reml;
				best_delta=result.getPoint();
			}
		}
		this.Likelihood_null=best_reml;		
		double[] beta_0=REMLFullRank4BrentOptimal.calculate_beta(Y_transformed, transformed_X, S, best_delta);
		this.sigma_g_null=REMLFullRank4BrentOptimal.sigma_g2_REMLL(Y_transformed, transformed_X, S, best_delta, beta_0);
		this.sigma_e_null=this.sigma_g_null*best_delta;
	}
	
	public void fullRankSolver_null_ML(){
		double[][] transformed_X=new double[1][];
		double[] intercept_1=new double[this.sample_size];
		Arrays.fill(intercept_1, 1);
		transformed_X[0]=this.global_Ut.operate(intercept_1);
		BrentOptimizer bo=new BrentOptimizer(FaSTLMM.brent_rel, FaSTLMM.brent_abs);
		double best_delta=-1;
		double best_ml=-(Double.MAX_VALUE);
		UnivariateFunction ml_function=new MLFullRank4BrentOptimal(transformed_X,Y_transformed,S);
		for(double min=lowlimit;min<uplimit;min+=grid_step){
			double min_delta=Math.exp(min);
			double max_delta=Math.exp(min+grid_step);
			UnivariatePointValuePair result=bo.optimize(maxEval, ml_function, GoalType.MAXIMIZE, min_delta, max_delta);
			double the_ml=result.getValue();
			if(the_ml>best_ml){
				best_ml=the_ml;
				best_delta=result.getPoint();
			}
		}
		this.Likelihood_null=best_ml;		
		double[] beta_0=MLFullRank4BrentOptimal.calculate_beta(Y_transformed, transformed_X, S, best_delta);
		this.sigma_g_null=MLFullRank4BrentOptimal.sigma_g2_LL(Y_transformed, transformed_X, S, best_delta, beta_0);
		this.sigma_e_null=this.sigma_g_null*best_delta;
	}
	
	public void analyze_REML(String output_csv_file, boolean plot, double maf_plot_threshold, double pvalue_threshold){
		double total_num_sites=0;
		for(int chr=0;chr<this.genotype.num_chrs;chr++)total_num_sites+=this.genotype.num_sites[chr];
		double corrected_threshold=pvalue_threshold/total_num_sites;
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(output_csv_file));
			bw.write("#SampleSize="+this.sample_size+ "; REML_null="+this.Likelihood_null+"; Sigma_g^2="+this.sigma_g_null+
					"; Sigma_e^2="+this.sigma_e_null+
					"; Phenotype="+this.phenotype.phe_id
					+"; GenotypeFile="+this.genotype_hdf5_file+"\n");
			bw.write("#chr,location,pvalue,AdjustedR2,coefficient,sigma_g2,MAF_count\n");
			System.out.println("start FaST-LMM");

			for(int chr=0;chr<this.genotype.num_chrs;chr++){
				int snp=0;
				for (HDF5MDDataBlock<MDDoubleArray> block : this.genotype.position_fast_blocks[chr]){
					double[][] data4thisblock=block.getData().toMatrix();
					for(int var_index=0;var_index<data4thisblock.length;var_index++){						
//						double[] X_ori=new double[this.sample_size];
//						for(int sample_index=0; sample_index<this.sample_size; sample_index++){
//							X_ori[sample_index] = data4thisblock[var_index][this.indexes_with_phenotype_in_genotype[sample_index]];
//						}
						double[] X_ori=this.extract_samples_with_phenotype(data4thisblock[var_index]);
						boolean the_same=true;
						for(int sample_index=1; sample_index<this.sample_size; sample_index++){
							if(Double.compare(X_ori[sample_index], X_ori[0])!=0)the_same=false;
						}				
						if(the_same)continue;
						double[][] Xs_after=new double[2][this.sample_size];	
						Xs_after[0]=this.intercept;
						Xs_after[1]=this.global_Ut.operate(X_ori);
						VarComp_Result the_result=fullRankSolver_FixedEffects_REML(Xs_after, Y_transformed, S, Likelihood_null);
						//double[][] result=EMMAX.reg2results_emmax(emma.phenotype.new_Y_4emmax, Xs_after);	
						if(the_result.pvalue<corrected_threshold){
							the_result.R2=the_result.beta[1]*the_result.beta[1]*StatUtils.populationVariance(X_ori)
									/StatUtils.populationVariance(this.phenotype.values);
							bw.write((chr+1)+","+this.genotype.locations[chr][snp+var_index]+","+
									the_result+","+EMMAX.mafc(X_ori)+"\n");
							bw.flush();
						}
																		
					}snp+=data4thisblock.length;
				}
				System.out.println("finished Chr"+(chr+1)+".");
			}			
			bw.close();
			if(plot) EMMAX.make_plot_one_phen(output_csv_file, maf_plot_threshold);			
		}catch(Exception e){e.printStackTrace();}
	}	
	
	public void analyze_ML(String output_csv_file, boolean plot, double maf_plot_threshold, double pvalue_threshold){
		double total_num_sites=0;
		for(int chr=0;chr<this.genotype.num_chrs;chr++)total_num_sites+=this.genotype.num_sites[chr];
		double corrected_threshold=pvalue_threshold/total_num_sites;
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(output_csv_file));
			bw.write("#SampleSize="+this.sample_size+ "; ML_null="+this.Likelihood_null+"; Sigma_g^2="+this.sigma_g_null+
					"; Sigma_e^2="+this.sigma_e_null+
					"; Phenotype="+this.phenotype.phe_id
					+"; GenotypeFile="+this.genotype_hdf5_file+"\n");
			bw.write("#chr,location,pvalue,AdjustedR2,coefficient,sigma_g2,MAF_count\n");
			System.out.println("start FaST-LMM");

			for(int chr=0;chr<this.genotype.num_chrs;chr++){
				int snp=0;
				for (HDF5MDDataBlock<MDDoubleArray> block : this.genotype.position_fast_blocks[chr]){
					double[][] data4thisblock=block.getData().toMatrix();
					for(int var_index=0;var_index<data4thisblock.length;var_index++){						
//						double[] X_ori=new double[this.sample_size];
//						for(int sample_index=0; sample_index<this.sample_size; sample_index++){
//							X_ori[sample_index] = data4thisblock[var_index][this.indexes_with_phenotype_in_genotype[sample_index]];
//						}
						double[] X_ori=this.extract_samples_with_phenotype(data4thisblock[var_index]);
						boolean the_same=true;
						for(int sample_index=1; sample_index<this.sample_size; sample_index++){
							if(Double.compare(X_ori[sample_index], X_ori[0])!=0)the_same=false;
						}				
						if(the_same)continue;
						double[][] Xs_after=new double[2][this.sample_size];	
						Xs_after[0]=this.intercept;
						Xs_after[1]=this.global_Ut.operate(X_ori);
						VarComp_Result the_result=fullRankSolver_FixedEffects_ML(Xs_after, Y_transformed, S, Likelihood_null);
						//double[][] result=EMMAX.reg2results_emmax(emma.phenotype.new_Y_4emmax, Xs_after);	
						if(the_result.pvalue<corrected_threshold){
							the_result.R2=the_result.beta[1]*the_result.beta[1]*StatUtils.populationVariance(X_ori)
									/StatUtils.populationVariance(this.phenotype.values);
							bw.write((chr+1)+","+this.genotype.locations[chr][snp+var_index]+","+
									the_result+","+EMMAX.mafc(X_ori)+"\n");
							bw.flush();
						}
																		
					}snp+=data4thisblock.length;
				}
				System.out.println("finished Chr"+(chr+1)+".");
			}			
			bw.close();
			if(plot) EMMAX.make_plot_one_phen(output_csv_file, maf_plot_threshold);			
		}catch(Exception e){e.printStackTrace();}
	}	
	
	public static VarComp_Result fullRankSolver_FixedEffects_REML(double[][] transformed_X, double[] transformed_Y, 
			double[] S, double reml_null){				
		BrentOptimizer bo=new BrentOptimizer(FaSTLMM.brent_rel, FaSTLMM.brent_abs);	
		double best_delta=-1;
		double best_reml=-(Double.MAX_VALUE);
		UnivariateFunction reml_function=new REMLFullRank4BrentOptimal(transformed_X,transformed_Y,S);
		for(double min=lowlimit;min<uplimit;min+=grid_step){
			double min_delta=Math.exp(min);
			double max_delta=Math.exp(min+grid_step);
			UnivariatePointValuePair result=bo.optimize(FaSTLMM.maxEval, reml_function, GoalType.MAXIMIZE, min_delta, max_delta);
			double the_reml=result.getValue();
			if(the_reml>best_reml){
				best_reml=the_reml;
				best_delta=result.getPoint();
			}
		}
		VarComp_Result result=new VarComp_Result();
		result.ml=best_reml;
		result.beta=REMLFullRank4BrentOptimal.calculate_beta(transformed_Y, transformed_X, S, best_delta);
		result.delta=best_delta;
		result.pvalue=myMathLib.Test.chi2pr(2*(best_reml-reml_null), 1);
		result.sigma_g=REMLFullRank4BrentOptimal.sigma_g2_REMLL(transformed_Y, transformed_X, S, best_delta, result.beta);
		result.sigma_e=result.sigma_g*best_delta;
		return result;		
	}	
	
	public static VarComp_Result fullRankSolver_FixedEffects_ML(double[][] transformed_X, double[] transformed_Y, 
			double[] S, double ml_null){				
		BrentOptimizer bo=new BrentOptimizer(FaSTLMM.brent_rel, FaSTLMM.brent_abs);	
		double best_delta=-1;
		double best_ml=-(Double.MAX_VALUE);
		UnivariateFunction ml_function=new MLFullRank4BrentOptimal(transformed_X,transformed_Y,S);
		for(double min=lowlimit;min<uplimit;min+=grid_step){
			double min_delta=Math.exp(min);
			double max_delta=Math.exp(min+grid_step);
			UnivariatePointValuePair result=bo.optimize(FaSTLMM.maxEval, ml_function, GoalType.MAXIMIZE, min_delta, max_delta);
			double the_ml=result.getValue();
			if(the_ml>best_ml){
				best_ml=the_ml;
				best_delta=result.getPoint();
			}
		}
		VarComp_Result result=new VarComp_Result();
		result.ml=best_ml;
		result.beta=MLFullRank4BrentOptimal.calculate_beta(transformed_Y, transformed_X, S, best_delta);
		if(result.beta==null){
			return null;
		}
		result.delta=best_delta;
		result.pvalue=myMathLib.Test.chi2pr(2*(best_ml-ml_null), 1);
		//System.out.println("best,"+best_ml);
		result.sigma_g=MLFullRank4BrentOptimal.sigma_g2_LL(transformed_Y, transformed_X, S, best_delta, result.beta);
		result.sigma_e=result.sigma_g*best_delta;
		return result;		
	}
	
	public static double fisher_info_nobeta(double[] transformed_Y, double[] S, double h2){
		double n=S.length;
		double fisher_info=0;
		double delta=1/h2-1;
		double delta_d1=-1/(h2*h2);
		double delta_d2=2/(h2*h2*h2);
		
		double first=0, second=0, third=0, noy_first=0, noy_second=0;
		for(int i=0;i<transformed_Y.length;i++){
			noy_first=noy_first+1.0/(S[i]+delta);
			noy_second=noy_second+1.0/((S[i]+delta)*(S[i]+delta));
			first=first+(transformed_Y[i])*(transformed_Y[i])/(S[i]+delta);
			second=second+(transformed_Y[i])*(transformed_Y[i])/((S[i]+delta)*(S[i]+delta));
			third=third+(transformed_Y[i])*(transformed_Y[i])/((S[i]+delta)*(S[i]+delta)*(S[i]+delta));
		}
		double LL_d1= -0.5*(n*(-second)/first+noy_first); 
		double LL_d2= -0.5*(n*(2*third*first)-second*second)/(first*first)-noy_second;
		fisher_info=delta_d2*LL_d1+LL_d2*delta_d1*delta_d1;
		return fisher_info;
	}
	
	public static void main(String[] args) {
		String genotype_file_hdf5="/Volumes/Projects/DATA/arabidopsis_genomes/swe180_ecker171_removebads.csv.num.mafc.hdf5";
		String global_kinship_file="/Volumes/Projects/DATA/arabidopsis_genomes/kinship_files_swe180_ecker171_removebads" +
				"/kinship.swe180_ecker171_removebads.rescaled.ibs";
		String phe_file="/Volumes/Projects/DATA/FernandoFlowCyto/genome_size.tsv";
		int phe_index=0;
		FaSTLMM_FullRank data=new FaSTLMM_FullRank(genotype_file_hdf5, global_kinship_file, new MultiPhenotype(phe_file).phenotypes[phe_index], true);
		double maf_plot_threshold=0.05;
		String output_csv_file="/Volumes/Projects/Local-kinship/genome_size/result_FaSTLMM/genomesize.resultJan26.ML.csv";
		data.analyze_ML(output_csv_file, true, maf_plot_threshold, 1000);
	}
}
