package mixedmodel;

import java.util.ArrayList;
import java.util.HashMap;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.util.MathArrays;
import org.apache.commons.math3.util.MathUtils;

public class FaSTLMM {
	
	public final static double machep=1e-10;//No need to use 2.220446049250313E-16; // Machine epsilon = Math.ulp(1d)
	public final static double brent_rel=Math.sqrt(machep);
	public final static double brent_abs=brent_rel;
	public final static int maxEval=100;	
	public final static double uplimit=10;
	public final static double lowlimit=-10;
	public final static double grid_step=0.1;
	
	
	String genotype_hdf5_file;
	public VariantsDouble genotype;
	int[] indexes_with_phenotype_in_genotype;
	Phenotype phenotype; // all samples with phenotype will have genotype in the object: 
						//the constructor will guarantee that by checking the data ids.
	double[][] global_kinship_matrix=null; // the order is corresponding to the phenotype, but generated from genotype matched kinship file.
	double[][] terms_for_kinship=null;
	int sample_size;		
	
//	public FaSTLMM(String genotype_hdf5_file, String phe_file, String global_kinship_file, int phe_index){
//		try{
//			VariantsDouble genotype=new VariantsDouble(genotype_hdf5_file);
//			this.genotype_hdf5_file=genotype_hdf5_file;
//			KinshipMatrix ori_kinship= new KinshipMatrix(global_kinship_file);			
//			double[][] ori_kinship_matrix= ori_kinship.kinship_matrix.getData();//genotype.read_in_kinship(global_kinship_file);
//			Phenotype ori_phenotype=(new MultiPhenotype(phe_file)).phenotypes[phe_index];
//			this.genotype=genotype;
//			ArrayList<Integer> useful_phe_ids=new ArrayList<Integer>();
//			HashMap<Integer, Integer> corresponding_geno_ids=new HashMap<Integer, Integer>();
////			HashMap<String, Integer> all_geno_ids=new HashMap<String, Integer>();
////			for(int geno_id_index=0;geno_id_index<genotype.sample_ids.length;geno_id_index++){
////				all_geno_ids.put(genotype.sample_ids[geno_id_index], geno_id_index);
////			}
//			for(int phe_id_index=0;phe_id_index<ori_phenotype.sample_ids.length;phe_id_index++){			
//				if(this.genotype.sample_id2index.containsKey(ori_phenotype.sample_ids[phe_id_index])){
//					int geno_id_index=this.genotype.sample_id2index.get(ori_phenotype.sample_ids[phe_id_index]);
//					useful_phe_ids.add(phe_id_index);
//					corresponding_geno_ids.put(phe_id_index, geno_id_index);					
//				}				
//			}
//			// assign the fields
//			this.sample_size=useful_phe_ids.size();
//			double[] the_values=new double[sample_size];
//			String[] the_ids=new String[sample_size];
//			this.indexes_with_phenotype_in_genotype=new int[sample_size];
//			for(int k=0;k<sample_size;k++){
//				int phe_id_index=useful_phe_ids.get(k);
//				the_ids[k]=ori_phenotype.sample_ids[phe_id_index];
//				the_values[k]=ori_phenotype.values[phe_id_index];
//				indexes_with_phenotype_in_genotype[k]=corresponding_geno_ids.get(phe_id_index);
//			}		
//			this.phenotype=new Phenotype(ori_phenotype.phe_id, the_ids, the_values);
//			// extract kinship
//			double[][] global_kinship_matrix_before=new double[sample_size][sample_size];
//			for(int i=0;i<sample_size;i++){
//				int gen_id_index1=indexes_with_phenotype_in_genotype[i];
//				for(int j=0;j<sample_size;j++){
//					int gen_id_index2=indexes_with_phenotype_in_genotype[j];
//					global_kinship_matrix_before[i][j]=ori_kinship_matrix[gen_id_index1][gen_id_index2];
//				}
//			}this.global_kinship_matrix=KinshipMatrix.re_scale_kinship_matrix(global_kinship_matrix_before);
//			//Test.write2file(the_ids, "/Users/quan.long/look/ids.new");	
//		}catch(Exception e){e.printStackTrace();}
//	}
	
	
	public static int sample_size(Phenotype ori_phenotype, String genotype_hdf5_file){
		VariantsDouble genotype=new VariantsDouble(genotype_hdf5_file);
		ArrayList<Integer> useful_phe_ids=new ArrayList<Integer>();
		HashMap<Integer, Integer> corresponding_geno_ids=new HashMap<Integer, Integer>();
		for(int phe_id_index=0;phe_id_index<ori_phenotype.sample_ids.length;phe_id_index++){			
			if(genotype.sample_id2index.containsKey(ori_phenotype.sample_ids[phe_id_index])){
				int geno_id_index=genotype.sample_id2index.get(ori_phenotype.sample_ids[phe_id_index]);
				useful_phe_ids.add(phe_id_index);
				corresponding_geno_ids.put(phe_id_index, geno_id_index);					
			}				
		}
		return useful_phe_ids.size();
	}
	public FaSTLMM(Phenotype ori_phenotype, String genotype_hdf5_file, String global_kinship_file){
		try{
			VariantsDouble genotype=new VariantsDouble(genotype_hdf5_file);
			this.genotype_hdf5_file=genotype_hdf5_file;
			//double[][] ori_kinship_matrix=  new KinshipMatrix(global_kinship_file).kinship_matrix_array;//genotype.read_in_kinship(global_kinship_file);
			//double[][] ori_kinship_matrix=  genotype.read_in_kinship(global_kinship_file);
			KinshipMatrix ori_kinship= new KinshipMatrix(global_kinship_file);			
			double[][] ori_kinship_matrix= ori_kinship.kinship_matrix.getData();//genotype.read_in_kinship(global_kinship_file);
						
			this.genotype=genotype;
			ArrayList<Integer> useful_phe_ids=new ArrayList<Integer>();
			HashMap<Integer, Integer> corresponding_geno_ids=new HashMap<Integer, Integer>();
//			HashMap<String, Integer> all_geno_ids=new HashMap<String, Integer>();
//			for(int geno_id_index=0;geno_id_index<genotype.sample_ids.length;geno_id_index++){
//				all_geno_ids.put(genotype.sample_ids[geno_id_index], geno_id_index);
//			}
			for(int phe_id_index=0;phe_id_index<ori_phenotype.sample_ids.length;phe_id_index++){			
				if(this.genotype.sample_id2index.containsKey(ori_phenotype.sample_ids[phe_id_index])){
					int geno_id_index=this.genotype.sample_id2index.get(ori_phenotype.sample_ids[phe_id_index]);
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
			// extract kinship
			double[][] global_kinship_matrix_before=new double[sample_size][sample_size];
			for(int i=0;i<sample_size;i++){
				int gen_id_index1=indexes_with_phenotype_in_genotype[i];
				for(int j=0;j<sample_size;j++){
					int gen_id_index2=indexes_with_phenotype_in_genotype[j];
					global_kinship_matrix_before[i][j]=ori_kinship_matrix[gen_id_index1][gen_id_index2];
				}
			}
			ori_kinship_matrix=null;
			this.global_kinship_matrix= KinshipMatrix.re_scale_kinship_matrix(global_kinship_matrix_before);
			global_kinship_matrix_before=null;
		}catch(Exception e){e.printStackTrace();}
	}
	
	/*
	 * constructor for low-rank case when the terms for kinship is in double[][] array. 
	 */
	public FaSTLMM(String genotype_hdf5_file, String phe_file, double[][] ori_terms_for_kinship, int phe_index){
		try{
			VariantsDouble genotype=new VariantsDouble(genotype_hdf5_file);
			this.genotype_hdf5_file=genotype_hdf5_file;
			MultiPhenotype phenotypeS=new MultiPhenotype(phe_file);
			Phenotype ori_phenotype=phenotypeS.phenotypes[phe_index];
			this.genotype=genotype;
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
			// extract terms_for_kinship			
			this.terms_for_kinship=this.extract_samples_with_phenotype(ori_terms_for_kinship);
			// normalize to keep the variance of the matrix K=W*W^t is 1
			// In FaST-LMM paper, the sqrt(#SNPs) has not been applied, which I THINK will cause the VC estimates wrong.
			normalize_terms_withSqrtNumVarDivided(this.terms_for_kinship);
		}catch(Exception e){e.printStackTrace();}
	}
	
	/*
	 * constructor for low-rank case when the terms for kinship is in a .tsv file that fits in phenotype format. 
	 */
	public FaSTLMM(String genotype_hdf5_file, String phe_file, MultiPhenotype ori_terms_for_kinship, int phe_index){
		try{
			VariantsDouble genotype=new VariantsDouble(genotype_hdf5_file);
			this.genotype_hdf5_file=genotype_hdf5_file;
			MultiPhenotype phenotypeS=new MultiPhenotype(phe_file);
			Phenotype ori_phenotype=phenotypeS.phenotypes[phe_index];
			this.genotype=genotype;
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
			// extract terms_for_kinship			
			this.terms_for_kinship=this.extract_samples_with_phenotype(ori_terms_for_kinship);
			// normalize to keep the variance of the matrix K=W*W^t is 1
			// In FaST-LMM paper, the sqrt(#SNPs) has not been applied, which I THINK will cause the VC estimates wrong.
			normalize_terms_withSqrtNumVarDivided(this.terms_for_kinship);
		}catch(Exception e){e.printStackTrace();}
	}
	
	// normalize to keep the variance of each row is 1 and mean of each row is 0.
	public static void normalize_terms(double[][] data){
		int sample_size=data[0].length;
		for(int k=0;k<data.length;k++){
			double mean=StatUtils.mean(data[k]);
			double var=StatUtils.variance(data[k], mean);
			for(int i=0;i<sample_size;i++){
				data[k][i]=(data[k][i]-mean)/Math.sqrt(var);
			}
		}
	}
	
	// normalize to keep the variance of the matrix K=W*W^t is 1
	public static void normalize_terms_withSqrtNumVarDivided(double[][] terms_for_kinship){
		int sample_size=terms_for_kinship[0].length;
		double num_of_terms=terms_for_kinship.length;
		for(int k=0;k<terms_for_kinship.length;k++){
			double mean=StatUtils.mean(terms_for_kinship[k]);
			double var=StatUtils.variance(terms_for_kinship[k], mean);
			for(int i=0;i<sample_size;i++){
				terms_for_kinship[k][i]=(terms_for_kinship[k][i]-mean)/Math.sqrt(var*num_of_terms);
			}
		}
	}
	
	// dividing N for all cells
	public static void sqrtNumVarDivided(double[][] terms_for_kinship){
		int sample_size=terms_for_kinship[0].length;
		double num_of_terms=terms_for_kinship.length;
		for(int k=0;k<terms_for_kinship.length;k++){
			for(int i=0;i<sample_size;i++){
				terms_for_kinship[k][i]=terms_for_kinship[k][i]/Math.sqrt(num_of_terms);
			}
		}
	}
	
	/*
	 * X_ori in the form of d x n where n is the sample size and d is the number of fixed effects. 
	 */
	public static double[][] transform_X_by_Ut(RealMatrix Ut, double[][] X_ori){
		double[][] transformed=new double[X_ori.length][];
		for(int d=0;d<X_ori.length;d++)transformed[d]=Ut.operate(X_ori[d]);
		return transformed;
	}

	/*
	 * select subset of data that have phenotype, the input usually are just grabbed from HDF5 genotypes
	 */
	public double[][] select_subset_of_data(double[][] full_data_region){
		double[][] selected_ori=new double[full_data_region.length][this.sample_size];
		boolean[] informative=new boolean[full_data_region.length];
		for(int var_index=0;var_index<full_data_region.length;var_index++){			
			for(int sample_index=0; sample_index<this.sample_size; sample_index++){
				selected_ori[var_index][sample_index] = full_data_region[var_index][this.indexes_with_phenotype_in_genotype[sample_index]];
			}			
			for(int sample_index=1; sample_index<this.sample_size; sample_index++){
				if(Double.compare(selected_ori[var_index][sample_index], selected_ori[var_index][0])!=0){
					informative[var_index]=true;
				}
			}				
		}
		int informative_vars_num=0;
		for(int var_index=0;var_index<full_data_region.length;var_index++){		
			if(informative[var_index])informative_vars_num++;
		}
		double[][] selected_informative=new double[informative_vars_num][];
		int the_index=0;
		for(int var_index=0;var_index<full_data_region.length;var_index++){		
			if(informative[var_index]){
				selected_informative[the_index]=selected_ori[var_index];
				the_index++;
			}
		}
		return selected_informative;
	}
	
	//TODO: never run
	public double[][] extract_samples_with_phenotype(double[][] ori_terms){
		int num_var=ori_terms.length;
		double[][] extrcted_array=new double[num_var][this.sample_size];
		for(int i=0;i<sample_size;i++){
			int gen_id_index1=indexes_with_phenotype_in_genotype[i];
			for(int k=0;k<num_var;k++){
				extrcted_array[k][i]=ori_terms[k][gen_id_index1];
			}
		}return extrcted_array;
	}
	
	/*
	 * The terms are stored as if they are phenotypes
	 */
	public double[][] extract_samples_with_phenotype(MultiPhenotype ori_terms){
		int num_var=ori_terms.num_of_pheno;
		double[][] extrcted_array=new double[num_var][this.sample_size];
		for(int i=0;i<sample_size;i++){
			int gen_id_index1=indexes_with_phenotype_in_genotype[i];
			String id_in_genotype=this.genotype.sample_ids[gen_id_index1];
			int correct_index=-1;
			for(int j=0;j<ori_terms.sample_size_all_pheno;j++){
				if(ori_terms.sample_ids_all_pheno[j].equals(id_in_genotype)){
					correct_index=j;break;
				}
			}
			if(correct_index!=-1){
				for(int k=0;k<num_var;k++){
					extrcted_array[k][i]=ori_terms.phenotypes[k].values[correct_index];
				}
			}else{
				System.out.println("Individual "+id_in_genotype+" with both genotype and phenotype doesn't have kinship-terms");
				for(int k=0;k<num_var;k++){
					extrcted_array[k][i]=Double.NaN;
				}
			}
		}return extrcted_array;
	}
	
	public double[] extract_samples_with_phenotype(double[] ori_term){
		double[] extrcted_array=new double[this.sample_size];
		for(int i=0;i<sample_size;i++){
			int gen_id_index1=indexes_with_phenotype_in_genotype[i];
			extrcted_array[i]=ori_term[gen_id_index1];			
		}return extrcted_array;
	}
	
	
	
	

}
