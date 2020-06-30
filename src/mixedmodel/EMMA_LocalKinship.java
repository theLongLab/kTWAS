package mixedmodel;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.Arrays;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

public class EMMA_LocalKinship extends EMMA {
	
	final double[] Y_global_transformed;
	final double[] intercept;
	final RealMatrix transformation_matrix; // =U'(S+delta*I)^(-1/2) 
	
	public EMMA_LocalKinship(Phenotype ori_phenotype, VariantsDouble genotype, double[][] ori_kinship_matrix){
		super(ori_phenotype, genotype, ori_kinship_matrix);
		int[] indexes={}, chrs={};
		this.REMLE(this.phenotype.values, this.cbind(chrs, indexes), null);	
		this.transformation_matrix=new Array2DRowRealMatrix(this.reml_decompositioned_array());
		
		this.Y_global_transformed=this.transformation_matrix.operate(this.phenotype.values);
		double[] all_one=new double[this.sample_size];
		Arrays.fill(all_one, 1);
		this.intercept=this.transformation_matrix.operate(all_one);		
		// remove rescaled intercept: NO NEED TO DO?		
		double mean_y=(new DescriptiveStatistics(this.Y_global_transformed)).getMean();
		double mean_beta0=(new DescriptiveStatistics(this.intercept)).getMean();
		double ratio=mean_y/mean_beta0;
		for(int k=0;k<this.Y_global_transformed.length;k++){
			this.Y_global_transformed[k]-=(ratio*this.intercept[k]);
		}
	}
	
	public void analysis_tiling_win(String mlr_output_file, int tiling_win_size, boolean plot){
		try{
			//VariantsDouble genotype=new VariantsDouble(genotype_hdf5_file);
			System.out.println(phenotype.phe_id);		
			BufferedWriter bw=new BufferedWriter(new FileWriter(mlr_output_file));
			bw.write("#SampleSize="+this.sample_size+ "; REML="+this.remle_REML+"; " +
					"h0-heritability="+this.heritability +"; delta="+this.remle_delta+"\n");
			bw.write("chr, position, p-value, variance_local, variance_global, variance_e\n");	
			final double[] all_one=new double[this.sample_size];
			Arrays.fill(all_one, 1);
			int total=0, low=0;
			double ml_null=myMathLib.StatFuncs.null_logL(this.Y_global_transformed.clone());
			for(int chr=0;chr<genotype.num_chrs;chr++){
				int last_position=genotype.locations[chr][genotype.num_sites[chr]-1];
				int start=1, end=start+tiling_win_size-1;
				while(start<last_position){
					total++;
					double[][] data_matching_phenotype=EMMAX.select_subset_of_data(genotype.load_variants_in_region(chr, start, end), this);					
					//double[][] transformed_data=data_matching_phenotype;
					if(data_matching_phenotype.length<2){
						start=start+tiling_win_size/2;
						end=start+tiling_win_size-1;
						continue;
					}
					FaSTLMM.normalize_terms_withSqrtNumVarDivided(data_matching_phenotype);
					double[][] transformed_data=EMMAX.transform_data(data_matching_phenotype, kinship_matrix);
					if(data_matching_phenotype.length<sample_size){ //low rank
						System.out.println("Low-Rank chr"+(chr+1)+" win "+start);
						low++;
						start=start+tiling_win_size/2;
						end=start+tiling_win_size-1;
						continue;
					}else if(data_matching_phenotype.length>=sample_size){							
						RealMatrix local_kinship= VariantsDouble.calculate_RRM_local_kinship_RealMatrix_noNomalization(transformed_data);
						//RealMatrix local_kinship=KinshipMatrix.re_scale_kinship_matrix(local_kinship0);
						double[] Y_new_original=this.Y_global_transformed;
						EigenDecomposition local_eigen= null;
						try{local_eigen= new EigenDecomposition((Array2DRowRealMatrix)local_kinship,0);}
						catch(Exception e){}
						if(local_eigen!=null){
							
							double[][] all_one_in_array=new double[1][];
							all_one_in_array[0]=all_one;
							double[][] transformed_X=EMMAX.transform_data(all_one_in_array, local_kinship.getData());
							VarComp_Result result0=EMMA.REMLE_static(Y_global_transformed, transformed_X, local_kinship.getData());
							
							double[] S=local_eigen.getRealEigenvalues(); //TODO:Uncommented this line and below 2 more lines to get rid of build error
							RealMatrix local_kinship_Ut= local_eigen.getVT();
							double[] transformed_Y=local_kinship_Ut.operate(Y_new_original);
							
							//double ml_null=myMathLib.StatFuncs.null_logL(transformed_Y,transformed_X[0]);
							//System.out.println(ml_null+","+this.Likelihood_null);
							VarComp_Result result=FaSTLMM_FullRank.fullRankSolver_FixedEffects_ML(transformed_X, transformed_Y, 
									S, ml_null);
							if(result!=null){
								double total_variance=result.sigma_g+result.sigma_e;//*(1+this.delta_null);
								/*bw.write((chr+1)+", "+start+", "+result.pvalue+", "+
										result.sigma_g/total_variance+", "+(result.sigma_e/total_variance)*(1/(1+this.delta_null))+", "+ 
										(result.sigma_e/total_variance)*(this.remle_delta/(1+this.delta_null))+"\n");*/	//TODO:Commented out due to not knowing what delta_null is
						}
						
						//	System.out.println("Succ finsihed chr"+(chr+1)+" win "+start);
						}
					}					
					start=start+tiling_win_size/2;
					end=start+tiling_win_size-1;
					//System.out.println("finsihed win "+start);
				}System.out.println("Succ finsihed chr"+(chr+1));
			}bw.close();
			System.out.println("Total #regions="+total+"; #Low-rank regions="+low);
			if(plot)LocalKinshipAnalyzer.make_plot_one_phen(mlr_output_file);
		}catch(Exception e){e.printStackTrace();}
	}

}
