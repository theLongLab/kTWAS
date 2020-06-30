package mixedmodel;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import nam.Cross_info;
import nam.Founder;
import nam.Imputation;
import nam.RILs;

import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;

import ch.systemsx.cisd.base.mdarray.MDDoubleArray;
import ch.systemsx.cisd.hdf5.HDF5MDDataBlock;

/**
 * Stepwise regression without mixed model
 * @author quan.long
 *
 */
public class LM {

	Phenotype phenotype; // all samples with phenotype will have genotype in hdf5: the constructor will guarantee that.
	int sample_size;		

	VariantsDouble genotype;
	int[] indexes_with_phenotype_in_genotype;
	int[] mafc;

	/*
	 * One has to make sure that the original kinship_matrix is corresponding to the sample_ids in genotype
	 * 
	 *  Find the sample_ids in both genotype and ori_phenotype. generate a clean phenotype dataset for further analysis.
	 */
	public LM(Phenotype ori_phenotype, VariantsDouble genotype){
		this.genotype=genotype;
		//		String[] id_phenotype=ori_phenotype.sample_ids;
		// find the ids
		ArrayList<Integer> useful_phe_ids=new ArrayList<Integer>();
		HashMap<Integer, Integer> corresponding_geno_ids=new HashMap<Integer, Integer>();
		HashMap<String, Integer> all_geno_ids=new HashMap<String, Integer>();
		for(int geno_id_index=0;geno_id_index<genotype.sample_ids.length;geno_id_index++){
			all_geno_ids.put(genotype.sample_ids[geno_id_index], geno_id_index);
		}
		for(int phe_id_index=0;phe_id_index<ori_phenotype.sample_ids.length;phe_id_index++){			
			if(all_geno_ids.containsKey(ori_phenotype.sample_ids[phe_id_index])){
				int geno_id_index=all_geno_ids.get(ori_phenotype.sample_ids[phe_id_index]);
				useful_phe_ids.add(phe_id_index);
				corresponding_geno_ids.put(phe_id_index, geno_id_index);					
			}				
		}
		// assign the fields
		this.sample_size=useful_phe_ids.size();
		double[] the_values=new double[sample_size];
		String[] the_ids=new String[sample_size];
		this.indexes_with_phenotype_in_genotype=new int[sample_size];
		for(int k=0;k<sample_size;k++){
			int phe_id_index=useful_phe_ids.get(k);
			the_ids[k]=ori_phenotype.sample_ids[phe_id_index];
			the_values[k]=ori_phenotype.values[phe_id_index];
			indexes_with_phenotype_in_genotype[k]=corresponding_geno_ids.get(phe_id_index);
		}		
		this.phenotype=new Phenotype(ori_phenotype.phe_id, the_ids, the_values);
	}



	public static void stepwise_analysis(String genotype_hdf5_file, String phe_file, String output_folder, int round,
			double p_value_after_correction, int min_sample_size, double maf_plot_threshold){
		try{
			VariantsDouble genotype=new VariantsDouble(genotype_hdf5_file);
			MultiPhenotype phenotypeS=new MultiPhenotype(phe_file);
			for(int phe_index=0;phe_index<phenotypeS.num_of_pheno;phe_index++){
				LM lm=new LM(phenotypeS.phenotypes[phe_index], genotype);				
				if(lm.sample_size>=min_sample_size){
					System.out.println("Running Phenotype:"+ phe_index+":"+lm.phenotype.phe_id);
					String single_marker_outfile =output_folder+"/LM."+phe_index+"_"+lm.phenotype.phe_id+".top";
					long startTime = System.currentTimeMillis();
					run_lm_or_stepwise_for_one_pheno(single_marker_outfile, p_value_after_correction,lm, maf_plot_threshold, genotype_hdf5_file, round);
					System.out.println("Time Consumed: "+ (System.currentTimeMillis() - startTime)/1000+" seconds.\n");	
				}else{
					System.out.println("Not running due to too small a sample size:" +lm.sample_size);
				}
			}			
		}catch(Exception e){e.printStackTrace();}
	}	

	public static void stepwise_analysis(String genotype_hdf5_file, String phe_file, String output_folder, int round,
			double p_value_after_correction, int min_sample_size, int phe_index, double maf_plot_threshold){
		try{
			VariantsDouble genotype=new VariantsDouble(genotype_hdf5_file);
			MultiPhenotype phenotypeS=new MultiPhenotype(phe_file);
			LM lm=new LM(phenotypeS.phenotypes[phe_index], genotype);			
			if(lm.sample_size>=min_sample_size){
				System.out.println("Running Phenotype:"+ phe_index+":"+lm.phenotype.phe_id);
				String single_marker_outfile =output_folder+"/LM."+phe_index+"_"+lm.phenotype.phe_id+".top";
				long startTime = System.currentTimeMillis();
				run_lm_or_stepwise_for_one_pheno(single_marker_outfile, p_value_after_correction,lm, maf_plot_threshold, genotype_hdf5_file, round);
				System.out.println("Time Consumed: "+ (System.currentTimeMillis() - startTime)/1000+" seconds.\n");	
			}else{
				System.out.println("Not running due to too small sample size:" +lm.sample_size);
			}					
		}catch(Exception e){e.printStackTrace();}
	}

	/*
	 * if round==1, then just run LM, otherwise run stepwise regression with the candidates in the first round LM.
	 */
	public static void run_lm_or_stepwise_for_one_pheno(String single_marker_outfile, double pvalue_threshold, LM lm, 
			double maf_plot_threshold, String genotype_hdf5_file, int round){				
		double total_num_sites=0;
		for(int chr=0;chr<lm.genotype.num_chrs;chr++)total_num_sites+=lm.genotype.num_sites[chr];
		double corrected_threshold=pvalue_threshold/total_num_sites;			
		try{			
			BufferedWriter bw_s = new BufferedWriter(new FileWriter(single_marker_outfile));
			bw_s.write("#SampleSize="+lm.sample_size+"; Phenotype="+lm.phenotype.phe_id+"; GenotypeFile="+genotype_hdf5_file+"\n");
			bw_s.write("#chr,location,pvalue,AdjustedR2,coefficient,Sd_Err,MAF_count\n");
			System.out.println("start LM");
			for(int chr=0;chr<lm.genotype.num_chrs;chr++){					
				int snp=0;
				for (HDF5MDDataBlock<MDDoubleArray> block : lm.genotype.position_fast_blocks[chr]){
					double[][] data4thisblock=block.getData().toMatrix();
					for(int var_index=0;var_index<data4thisblock.length;var_index++){
						double[][] Xs=new double[lm.sample_size][1];
						double[] X_ori=new double[lm.sample_size];
						for(int sample_index=0; sample_index<lm.sample_size; sample_index++){
							Xs[sample_index][0] = data4thisblock[var_index][lm.indexes_with_phenotype_in_genotype[sample_index]];
							X_ori[sample_index] = Xs[sample_index][0];
						}
						boolean the_same=true;
						for(int sample_index=1; sample_index<lm.sample_size; sample_index++){
							if(Double.compare(Xs[sample_index][0], Xs[0][0])!=0)the_same=false;
						}				
						if(the_same)continue;
						//OLSMultipleLinearRegression reg1=new OLSMultipleLinearRegression();									
						//reg1.newSampleData(lm.phenotype.values, Xs);
						double[][] result=EMMAX.reg2results_lm(lm.phenotype.values, Xs);
						if(result!=null){
							result[3]=new double[1];
							result[3][0]=EMMAX.mafc(X_ori);
							String out_res=EMMAX.results2string(result, corrected_threshold);
							if(out_res!=null){
								bw_s.write((chr+1)+","+lm.genotype.locations[chr][snp+var_index]+","+out_res+"\n");
							}
						}						
					}snp+=data4thisblock.length;
				}
			}
			bw_s.close();
			EMMAX.make_plot_one_phen(single_marker_outfile, maf_plot_threshold);				
			if(round>=2){ 
				// set the initial variant
				AssociationResults candidates=new AssociationResults(single_marker_outfile, corrected_threshold, maf_plot_threshold);
				candidates.generating_indexes(lm.genotype);
				int[] best_markers_indexes=new int[round];
				int[] best_markers_chrs=new int[round];
				int[] best_markers_locs=new int[round];		
				OLSMultipleLinearRegression best_reg=null;
				double[] var_explained=new double[round];
				double[] pvalues_t_test=new double[round];
				int best_P_index=0;
				// find the best P-value
				double best_var=0;
				for(int var_index=0;var_index<candidates.num_of_var;var_index++){
					if(Double.compare(best_var,candidates.AdjustedR2[var_index])<0){
						best_var=candidates.AdjustedR2[var_index];
						best_P_index=var_index;
					}
				}
				HashSet<Integer> best_indexes_set =new HashSet<Integer>();
				best_markers_indexes[0]=best_P_index;
				var_explained[0]=candidates.AdjustedR2[best_P_index];
				best_indexes_set.add(best_P_index);
				best_markers_chrs[0]=candidates.chr[best_P_index]-1;
				best_markers_locs[0]=candidates.location[best_P_index];		
				
				for(int round_index=2; round_index<=round;round_index++){
					// prepare the array including existing best markers
					double[][] Xs_ori=new double[lm.sample_size][round_index];
					for(int k=0;k<round_index-1;k++){
						int chr=best_markers_chrs[k], location=best_markers_locs[k];
						double[] the_full_list=lm.genotype.load_one_variant_by_location(chr, location);
						for(int sample_index=0; sample_index<lm.sample_size; sample_index++){
							Xs_ori[sample_index][k] = the_full_list[lm.indexes_with_phenotype_in_genotype[sample_index]];
						}
					}
					// add the new marker and do regression.
					best_P_index=0;
					double best_vc=0;	
					double best_pvalue=-1;
					System.out.println("Start Step-wise. Round="+round_index);
					for(int var_index=0;var_index<candidates.num_of_var;var_index++){
						if(best_indexes_set.contains(var_index))continue;
						int chr=candidates.chr[var_index]-1;	
						double[] the_full_list=lm.genotype.load_one_variant_by_index(chr, candidates.indexes_in_genotype[var_index]);
						for(int sample_index=0; sample_index<lm.sample_size; sample_index++){
							Xs_ori[sample_index][round_index-1] = the_full_list[lm.indexes_with_phenotype_in_genotype[sample_index]];
						}
						OLSMultipleLinearRegression reg1=new OLSMultipleLinearRegression();									
						reg1.newSampleData(lm.phenotype.values, Xs_ori);
						double[][] result=EMMAX.reg2results_lm(lm.phenotype.values, Xs_ori);
						if(result!=null && Double.compare(best_vc,result[2][0])<0){
							best_P_index=var_index;
							best_vc=result[2][0];
							best_reg=reg1;
							best_pvalue=result[1][result[1].length-1];
						}
					}							
					// fill the structure that records existing selections
					var_explained[round_index-1]=best_vc;
					pvalues_t_test[round_index-1]=best_pvalue;
					best_indexes_set.add(best_P_index);
					best_markers_chrs[round_index-1]=candidates.chr[best_P_index]-1;
					best_markers_locs[round_index-1]=candidates.location[best_P_index];
					
				}
				// output the results 
				System.out.println("Start writing stepwise results.");
				BufferedWriter bw_stw=new BufferedWriter(new FileWriter(single_marker_outfile+".stepwise.csv"));
				bw_stw.write("#Stepwise regression using variants reported in"+single_marker_outfile+"\n");
				bw_stw.write("#Round,chr,location,p-value(t),AdjustedR2\n");
				for(int k=0;k<round;k++){
					bw_stw.write((k+1)+","+(best_markers_chrs[k]+1)+","+best_markers_locs[k]+","+pvalues_t_test[k]+","+var_explained[k]+"\n");
				}
				bw_stw.write("\nRegressionParameters:\n");
				double[] coef_last=best_reg.estimateRegressionParameters();
				for(int i=0;i<round;i++)bw_stw.write(coef_last[i]+",");
				bw_stw.write(coef_last[round]+"\n");
				bw_stw.close();
			}
		}catch(Exception e){e.printStackTrace();}
	}
	
	public static void lm_select(String genotype_hdf5_file, String input_pheno, int phe_index, 
			String output_file, int[] chrs, int[] locs){
		try{
			
			int num_snp=chrs.length;			
			VariantsDouble genotype=new VariantsDouble(genotype_hdf5_file);
			MultiPhenotype phenotypeS=new MultiPhenotype(input_pheno);
			LM lm=new LM(phenotypeS.phenotypes[phe_index], genotype);	
			BufferedWriter bw=new BufferedWriter(new FileWriter(output_file));
			bw.write("#SampleSize="+lm.sample_size+"; Phenotype="+lm.phenotype.phe_id+"; GenotypeFile="+genotype_hdf5_file+"\n");
			double[][] Xs_ori=new double[lm.sample_size][num_snp];
			for(int k=0;k<num_snp;k++){				
				double[] the_full_list=lm.genotype.load_one_variant_by_location(chrs[k]-1, locs[k]);
				for(int sample_index=0; sample_index<lm.sample_size; sample_index++){
					Xs_ori[sample_index][k] = the_full_list[lm.indexes_with_phenotype_in_genotype[sample_index]];
				}
			}
			double[][] result=EMMAX.reg2results_lm(lm.phenotype.values, Xs_ori);	
			bw.write("RegressionParameters:");
			for(int k=0;k<result[0].length;k++)bw.write("\t"+result[0][k]);
			bw.write("\nIndividual_P-values(t-tests):");
			for(int k=0;k<result[1].length;k++)bw.write("\t"+result[1][k]);
			bw.write("\nAdjustedRSquared:\t"+result[2][0]);
			bw.write("\nModel_Pvalue(F-statistics):\t"+result[3][0]+"\n");
			bw.close();															
		}catch(Exception e){e.printStackTrace();}		
	}
	
	public static void stepwise_reg_nam(String genotype_hdf5_file, String phenotype_file, int phe_index,
			String output_prefix, double pvalue_threshold, double maf_plot_threshold, int round){		
		MultiPhenotype phenotypeS=new MultiPhenotype(phenotype_file);
		VariantsDouble genotype=new VariantsDouble(genotype_hdf5_file);
		if(phe_index==-1){
			for(int i=0;i<phenotypeS.num_of_pheno;i++){
				phenotypeS.phenotypes[i].substract_pedigree_mean_for_NAM();
				LM lm=new LM(phenotypeS.phenotypes[i], genotype);
				LM.run_lm_or_stepwise_for_one_pheno(output_prefix+"."+phenotypeS.phenotypes[i].phe_id+".top", 
						pvalue_threshold, lm, maf_plot_threshold, genotype_hdf5_file, round);
			}
		}
		else{
			phenotypeS.phenotypes[phe_index].substract_pedigree_mean_for_NAM();
			LM lm=new LM(phenotypeS.phenotypes[phe_index], genotype);
			LM.run_lm_or_stepwise_for_one_pheno(output_prefix+"."+phenotypeS.phenotypes[phe_index].phe_id+".top", 
					pvalue_threshold, lm, maf_plot_threshold, genotype_hdf5_file, round);
		}
		
		
	}
	
}









