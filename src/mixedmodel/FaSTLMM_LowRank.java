package mixedmodel;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.Arrays;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;
import org.apache.commons.math3.optimization.GoalType;
import org.apache.commons.math3.optimization.univariate.BrentOptimizer;
import org.apache.commons.math3.optimization.univariate.UnivariatePointValuePair;
import org.apache.commons.math3.stat.StatUtils;

import ch.systemsx.cisd.base.mdarray.MDDoubleArray;
import ch.systemsx.cisd.hdf5.HDF5MDDataBlock;

public class FaSTLMM_LowRank extends FaSTLMM{


	
	final RealMatrix global_U1t; // eigen vectors 
	final RealMatrix global_IminusU1U1t; // In - U1*U1t
	final double[] S; //Singular values
	final double[] Y_before_transformation;
	final double[] Y_transformed;
	final double[] Y_transformed2;
	final double[] all_one;
	final int transformed_sample_size; 
	double ML_null;
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
	 * Before entering the constructor below, one has to guarantee the number of terms for use in kinship is 
	 * lower than the sample size that will for sure generate a low-rank matrix in standard model. 
	 * Otherwise use FaSTLMM_FullRank.java
	 */
	public FaSTLMM_LowRank(String genotype_file_hdf5, double[][] terms_for_kinship, String phe_file, int phe_index){
		// normalization of kinship terms and substract sample will be done in the super() constructor.
		super(genotype_file_hdf5, phe_file, terms_for_kinship, phe_index);
		SingularValueDecomposition svd=new SingularValueDecomposition(
				 (new Array2DRowRealMatrix(this.terms_for_kinship)).transpose());
		this.global_U1t=svd.getUT();
		this.S=svd.getSingularValues();
		this.transformed_sample_size=this.S.length;
		for(int k=0;k<this.S.length;k++)this.S[k]=(this.S[k]*this.S[k]);
		this.Y_before_transformation=this.phenotype.values.clone();
		this.Y_transformed=this.global_U1t.operate(this.Y_before_transformation);  // after operated, the dimension became less than sample_size		
		this.all_one=new double[this.sample_size];
		Arrays.fill(this.all_one, 1);
		if(this.Y_transformed.length!=this.transformed_sample_size){
			System.out.println("WRONG: this.Y_transformed.length!=this.transformed_sample_size");
		}
		this.global_IminusU1U1t=svd.getU().multiply(this.global_U1t).scalarMultiply(-1);
		for(int i=0;i<this.sample_size;i++)this.global_IminusU1U1t.setEntry(i, i, 1+this.global_IminusU1U1t.getEntry(i, i));
		this.Y_transformed2=this.global_IminusU1U1t.operate(this.Y_before_transformation);
		// set up null-REML in the model that only fixed effect is the transformed intercept. 
		this.lowRankSolver_null();		
	}
		
	public static double[][] read_in_terms(String file_terms_for_kinship){
		double[][] terms=null;
		try{
			
		}catch(Exception e){e.printStackTrace();}
		return terms;
	}
	public FaSTLMM_LowRank(String genotype_file_hdf5, String file_terms_for_kinship, String phe_file, int phe_index){
		// normalization of kinship terms and substract sample will be done in the super() constructor.
		super(genotype_file_hdf5, phe_file, new MultiPhenotype(file_terms_for_kinship), phe_index);
		SingularValueDecomposition svd=new SingularValueDecomposition(
				(new Array2DRowRealMatrix(this.terms_for_kinship)).transpose());
		this.global_U1t=svd.getUT();
		this.S=svd.getSingularValues();
		this.transformed_sample_size=this.S.length;
		for(int k=0;k<this.S.length;k++)this.S[k]=(this.S[k]*this.S[k]);
		this.Y_before_transformation=this.phenotype.values.clone();
		this.Y_transformed=this.global_U1t.operate(this.Y_before_transformation);  // after operated, the dimension became less than sample_size		
		this.all_one=new double[this.sample_size];
		Arrays.fill(this.all_one, 1);
		if(this.Y_transformed.length!=this.transformed_sample_size){
			System.out.println("WRONG: this.Y_transformed.length!=this.transformed_sample_size");
		}
		this.global_IminusU1U1t=svd.getU().multiply(this.global_U1t).scalarMultiply(-1);
		for(int i=0;i<this.sample_size;i++)this.global_IminusU1U1t.setEntry(i, i, 1+this.global_IminusU1U1t.getEntry(i, i));
		this.Y_transformed2=this.global_IminusU1U1t.operate(this.Y_before_transformation);
		// set up null-REML in the model that only fixed effect is the transformed intercept. 
		this.lowRankSolver_null();
	}
	
	public void lowRankSolver_null(){
		double[][] X_before_transformation=new double[1][];		
		X_before_transformation[0]=this.all_one;
		BrentOptimizer bo=new BrentOptimizer(FaSTLMM.brent_rel, FaSTLMM.brent_abs);
		double best_delta=-1;
		double best_ml=-(Double.MAX_VALUE);
		MLLowRank4BrentOptimal ml_function=new MLLowRank4BrentOptimal(X_before_transformation, this.Y_before_transformation, 
				this.Y_transformed, this.Y_transformed2, this.S, this.global_U1t, this.global_IminusU1U1t);
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
		this.ML_null=best_ml;		
		double[] beta_0=MLLowRank4BrentOptimal.calculate_beta(Y_transformed, ml_function.transformed_X, 
				ml_function.additional_first_term, ml_function.additional_second_term, S, best_delta);
		this.sigma_g_null=MLLowRank4BrentOptimal.sigma_g2_LL(Y_transformed, ml_function.transformed_X, 
				Y_transformed2, ml_function.transformed2_X,	S, best_delta, beta_0);
		this.sigma_e_null=this.sigma_g_null*best_delta;
	}
	
	public void analyze(String output_csv_file, boolean plot, double maf_plot_threshold, double pvalue_threshold){
		double total_num_sites=0;
		for(int chr=0;chr<this.genotype.num_chrs;chr++)total_num_sites+=this.genotype.num_sites[chr];
		double corrected_threshold=pvalue_threshold/total_num_sites;
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(output_csv_file));
			bw.write("#SampleSize="+this.sample_size+ "; REML_null="+this.ML_null+"; Sigma_g^2="+this.sigma_g_null+
					"; Sigma_e^2="+this.sigma_e_null+
					"; Phenotype="+this.phenotype.phe_id
					+"; GenotypeFile="+this.genotype_hdf5_file+"\n");
			bw.write("#chr,location,pvalue,AdjustedR2,coefficient,sigma_g2,MAF_count\n");
			System.out.println("started Low-Ranked FaST-LMM");
			for(int chr=0;chr<this.genotype.num_chrs;chr++){
				int snp=0;
				for (HDF5MDDataBlock<MDDoubleArray> block : this.genotype.position_fast_blocks[chr]){
					double[][] data4thisblock=block.getData().toMatrix();
					for(int var_index=0;var_index<data4thisblock.length;var_index++){
						double[] X_ori=new double[this.sample_size];
						for(int sample_index=0; sample_index<this.sample_size; sample_index++){
							X_ori[sample_index] = data4thisblock[var_index][this.indexes_with_phenotype_in_genotype[sample_index]];
						}
						boolean the_same=true;
						for(int sample_index=1; sample_index<this.sample_size; sample_index++){
							if(Double.compare(X_ori[sample_index], X_ori[0])!=0)the_same=false;
						}				
						if(the_same)continue;
						double[][] Xs_before_transformation=new double[2][];	
						Xs_before_transformation[0]=this.all_one;
						Xs_before_transformation[1]=X_ori;
						VarComp_Result the_result=lowRankSolver_FixedEffects(Xs_before_transformation, this.Y_before_transformation, 
								this.Y_transformed, this.Y_transformed2, this.S, this.ML_null, this.global_U1t, this.global_IminusU1U1t);
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
	
	public static VarComp_Result lowRankSolver_FixedEffects(double[][] X_before_transformation, double[] Y_before_transformation, 
			double[] transformed_Y, double[] transformed2_Y, double[] S, double ml_null, RealMatrix U1t, RealMatrix IminusU1U1t){				
		BrentOptimizer bo=new BrentOptimizer(FaSTLMM.brent_rel, FaSTLMM.brent_abs);	
		double best_delta=-1;
		double best_ml=-(Double.MAX_VALUE);
		MLLowRank4BrentOptimal ml_function=new MLLowRank4BrentOptimal(X_before_transformation, Y_before_transformation, 
				transformed_Y, transformed2_Y, S, U1t, IminusU1U1t);
		
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
		result.beta=MLLowRank4BrentOptimal.calculate_beta(transformed_Y, ml_function.transformed_X, 
				ml_function.additional_first_term, ml_function.additional_second_term, S, best_delta);
		result.pvalue=myMathLib.Test.chi2pr(2*(best_ml-ml_null), 1);
		result.sigma_g=MLLowRank4BrentOptimal.sigma_g2_LL(transformed_Y, ml_function.transformed_X, transformed2_Y, 
				ml_function.transformed2_X, S, best_delta, result.beta);
		result.sigma_e=result.sigma_g*best_delta;
//		double[] y_minus_beta0=myMathLib.ArrayFuncs.linearCombination(transformed_Y, transformed_X[0], 1, -result.beta[0]);
//		double[] X_summed=new double[transformed_Y.length];
//		for(int k=1;k<transformed_X.length;k++){
//			X_summed=myMathLib.ArrayFuncs.linearCombination(X_summed, transformed_X[k], 1, result.beta[k]);
//		}
		//
		
		return result;		
	}	
	
	public static void main(String[] args) {
		String genotype_file_hdf5="/Volumes/Projects/DATA/arabidopsis_genomes/swe180_ecker171_removebads.csv.num.mafc.hdf5";
		String global_kinship_file="/Volumes/Projects/DATA/arabidopsis_genomes/kinship_files_swe180_ecker171_removebads" +
				"/kinship.swe180_ecker171_removebads.rescaled.ibs";
		String phe_file="/Volumes/Projects/DATA/FernandoFlowCyto/genome_size.tsv";
		int phe_index=0;
		FaSTLMM_LowRank data=new FaSTLMM_LowRank(genotype_file_hdf5, global_kinship_file, phe_file, phe_index);
		double maf_plot_threshold=0.05;
		String output_csv_file="/Volumes/Projects/Local-kinship/genome_size/result_FaSTLMM/genomesize.resultJan26.csv";
		data.analyze(output_csv_file, true, maf_plot_threshold, 1000);
	}


}
