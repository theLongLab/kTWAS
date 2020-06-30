package mixedmodel;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import myMathLib.Test;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.random.CorrelatedRandomVectorGenerator;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;

public class GBLUP {

	VariantsDouble obs_geno;
	VariantsDouble pred_geno;
	Phenotype obs_pheno;
	public Phenotype pred_pheno;
	KinshipMatrix obs_matrix;
	RelationMatrix relation_obs_pred;
	double h2;
	
	public GBLUP(VariantsDouble obs_geno, VariantsDouble pred_geno,	Phenotype obs_pheno, 
			KinshipMatrix obs_matrix, RelationMatrix relation_obs_pred,	double h2){
		this.obs_geno= obs_geno;
		this.pred_geno= pred_geno;
		this.obs_pheno= obs_pheno;
		this.obs_matrix= obs_matrix;
		this.relation_obs_pred= relation_obs_pred;
		this.h2= h2;		
	}
	
	public static String[] select_ids_random(int num_select, Phenotype pheno, VariantsDouble geno_full){
		System.out.println("Started selecting observed sample ids.");
		String[] selected=new String[num_select];
		int index=0;
		HashSet<Integer> selected_table=new HashSet<Integer>();
		String[] full=pheno.sample_ids;
		for(int k=0;k<num_select;k++){
			int the_selectedindex=(int)(Test.randomNumber()*full.length);
			if(the_selectedindex==full.length)the_selectedindex--;
			while(selected_table.contains(the_selectedindex)||Double.isNaN(pheno.values[the_selectedindex]) 
					|| (!geno_full.sample_id2index.containsKey(full[the_selectedindex]))){
				the_selectedindex=(int)(Test.randomNumber()*full.length);
				if(the_selectedindex==full.length)the_selectedindex--;
			}
			selected[index++]=full[the_selectedindex];
			selected_table.add(the_selectedindex);
		}
		System.out.println("Finished selecting observed sample ids.");
		return selected;
	}
	
	/*
	 * Given a dataset with kinship and phenotype known, partition it into two parts and try whether BLUP
	 * can predict the test part well. 
	 * 
	 * Global kinship ONLY
	 */
	public static void test_batch(Phenotype pheno_full, KinshipMatrix kinship_full, String obs_ids[], double h2,
			String output_file){
		try{
			String[] test_ids=the_rest_ids(pheno_full.sample_ids, obs_ids, kinship_full.ids2indexes);
			System.out.println("obs_ids:"+obs_ids.length+"/test_ids:"+test_ids.length);
			Phenotype obs_pheno= pheno_full.extract_sample_with_id_list(obs_ids);
			Phenotype true_pred_pheno= pheno_full.extract_sample_with_id_list(test_ids);
			KinshipMatrix obs_matrix=KinshipMatrix.submatrix(kinship_full, obs_ids);
			RelationMatrix relation_obs_pred=RelationMatrix.submatrix(kinship_full, obs_ids, test_ids);
			GBLUP gblup=new GBLUP(null, null, obs_pheno, obs_matrix, relation_obs_pred, h2);
			System.out.println("finished submatrixes.");
			gblup.estimate_global();
			double[] corr=MultiPhenotype.correlation(true_pred_pheno, gblup.pred_pheno);
			System.out.println(corr[0]+"\t"+corr[1]+"\t"+corr[2]+"\t"+corr[3]);
		}catch(Exception e){e.printStackTrace();}
	}
	

	/*
	 * Given a dataset with kinship and phenotype known, partition it into two parts and try whether BLUP
	 * can predict the test part well. 
	 * 
	 * Global kinship PLUS local terms
	 * 
	 * if geno_full==null, calculate local kinships and store the results into the paths "String[] local_kinship_files", 
	 * otherwise load them from the paths. 
	 */
	public static double[] test_batch(Phenotype pheno_full, VariantsDouble geno_full,
			KinshipMatrix kinship_full, String obs_ids[], double h2, 
			int[][] local_regions, double[] var_comp, String[] local_kinship_files, String output_file, 
			int first_sample_size){
		double[] corrs=new double[3];
		try{
			String[] test_ids=the_rest_ids(pheno_full.sample_ids, obs_ids, kinship_full.ids2indexes);
			Phenotype obs_pheno= pheno_full.extract_sample_with_id_list(obs_ids);
			//DEBUG for(int k=0;k<test_ids.length;k++)System.out.println(k+": "+test_ids[k]);
			Phenotype true_pred_pheno= pheno_full.extract_sample_with_id_list(test_ids);
			//System.out.println("sub_pheno done.");
			KinshipMatrix obs_matrix=KinshipMatrix.submatrix(kinship_full, obs_ids);
			//System.out.println("sub_kinship done.");
			RelationMatrix relation_obs_pred=RelationMatrix.submatrix(kinship_full, obs_ids, test_ids);
			//System.out.println("sub_relation done.");
			GBLUP gblup=new GBLUP(null, null, obs_pheno, obs_matrix, relation_obs_pred, h2);
			System.out.println("Finished global submatrixes.");
			gblup.estimate_global();
			double[] corr=MultiPhenotype.correlation(true_pred_pheno, gblup.pred_pheno);
			System.out.println("Global: "+corr[0]+"\t"+corr[1]+"\t"+corr[2]+"\t"+corr[3]);
			
			CausalRegions local=new CausalRegions(local_regions, var_comp, local_kinship_files, 
					obs_matrix, relation_obs_pred, obs_ids, test_ids);			
			System.out.println("finished local submatrixes.");
			gblup.estimate_local_only_sum(local, var_comp);
			double[] corr2=MultiPhenotype.correlation(true_pred_pheno, gblup.pred_pheno);
			System.out.println("Local_Sum: "+corr2[0]+"\t"+corr2[1]+"\t"+corr2[2]+"\t"+corr2[3]); 
			gblup.estimate_local_reg_(local, var_comp, geno_full, first_sample_size);
			double[] corr3=MultiPhenotype.correlation(true_pred_pheno, gblup.pred_pheno);
			corrs[0]=corr[0];
			corrs[1]=corr2[0];
			corrs[2]=corr3[0];
		}catch(Exception e){e.printStackTrace();}
		return corrs;
	}
	
	/*
	 * estimate a term
	 */
	public static double[] estimate(KinshipMatrix obs_matrix, Phenotype obs_pheno, RelationMatrix relation_obs_pred, double h2){
		double[][] data=obs_matrix.kinship_matrix.getData();
		double[][] addedI=new double[data.length][data[0].length];
		//h2*G11+(1-h2)I where G11 is the kinship matrix and I is identity matrix
		for(int i=0;i<data.length;i++){
			for(int j=0;j<data[0].length;j++){
				if(i!=j)
					addedI[i][j]= data[i][j]*h2;
				else
					addedI[i][j]= data[i][j]*h2+(1.0-h2);
			}
		}
		LUDecomposition desomposition=new LUDecomposition(new Array2DRowRealMatrix(addedI));
		RealMatrix inverse=desomposition.getSolver().getInverse();
		double[] Y_trans=inverse.operate(obs_pheno.values);
		double[] Y_pred=relation_obs_pred.relation_matrix.operate(Y_trans);
		for(int k=0;k<Y_pred.length;k++)Y_pred[k]=Y_pred[k]*h2;
		return Y_pred;
	}
	
	/*
	 * use the matrix without additional terms
	 */
	public void estimate_global(){
		System.out.println("Started estimating G-BLUP.");
		double[] Y_pred=estimate(obs_matrix, obs_pheno, relation_obs_pred, h2);
		this.pred_pheno= new Phenotype(obs_pheno.phe_id, relation_obs_pred.pred_ids.clone(), Y_pred);
		System.out.println("Finished estimating G-BLUP.");
	}
	
	/*
	 * use the matrix WITH additional terms
	 */
	public void estimate_local(CausalRegions regions, double[] var_comp){
		System.out.println("Started estimating G-BLUP.");
		double[] Y_pred=estimate(obs_matrix, obs_pheno, relation_obs_pred, h2);
		double[][] Y_pred_local=new double[var_comp.length][];
		for(int r=0;r<var_comp.length;r++){
			Y_pred_local[r]=estimate(regions.reginal_kinships[r], 
					obs_pheno, regions.relational_matrixs[r], var_comp[r]);
			for(int k=0;k<Y_pred.length;k++)Y_pred[k]+=Y_pred_local[r][k];
		}
		this.pred_pheno= new Phenotype(obs_pheno.phe_id, relation_obs_pred.pred_ids.clone(), Y_pred);
		System.out.println("Finished estimating G-BLUP.");
	}
	
	/*
	 * use the matrix WITH additional terms
	 */
	public void estimate_local_no_global(CausalRegions regions, double[] var_comp){
		System.out.println("Started estimating G-BLUP.");
		double[] Y_pred=new double[this.pred_geno.sample_ids.length];//estimate(obs_matrix, obs_pheno, relation_obs_pred, h2);
		double[][] Y_pred_local=new double[var_comp.length][];
		for(int r=0;r<var_comp.length;r++){
			Y_pred_local[r]=estimate(regions.reginal_kinships[r], 
					obs_pheno, regions.relational_matrixs[r], var_comp[r]);
			for(int k=0;k<Y_pred.length;k++)Y_pred[k]+=Y_pred_local[r][k];
		}
		this.pred_pheno= new Phenotype(obs_pheno.phe_id, relation_obs_pred.pred_ids.clone(), Y_pred);
		System.out.println("Finished estimating G-BLUP.");		
	}	
	
	
	/*
	 * use the matrix WITH additional terms
	 */
	public void estimate_global_local_sum(CausalRegions regions, double[] var_comp){
		System.out.println("Started estimating G-BLUP.");
		regions.generte_combined(var_comp);
		double[] Y_pred=estimate(regions.combined_kinship, obs_pheno, regions.combined_relation, h2);
		this.pred_pheno= new Phenotype(obs_pheno.phe_id, relation_obs_pred.pred_ids.clone(), Y_pred);
		System.out.println("Finished estimating G-BLUP.");		
	}
	
	public void estimate_local_only_sum(CausalRegions regions, double[] var_comp){
		System.out.println("Started estimating G-BLUP: local_only_sum.");
		regions.generte_combined_no_global(var_comp);
		double[] Y_pred=estimate(regions.combined_kinship, obs_pheno, regions.combined_relation, h2);
		this.pred_pheno= new Phenotype(obs_pheno.phe_id, relation_obs_pred.pred_ids.clone(), Y_pred);
		System.out.println("Finished estimating G-BLUP: local_only_sum.");		
	}

	/*
	 * use the matrix WITH additional terms
	 * var_comp[0]=h2 global
	 */
	public void estimate_local_reg(CausalRegions regions, double[] var_comp){
		double[] betas=train_betas(regions, var_comp);
		System.out.println("Started estimating G-BLUP using regression.");
		if(var_comp.length+2!=betas.length){
			System.out.println("var_comp.length+2!=betas.length");
			return;
		}
		double[] Y_pred=estimate(obs_matrix, obs_pheno, relation_obs_pred, h2);
		for(int k=0;k<Y_pred.length;k++)Y_pred[k]=Y_pred[k]*betas[1]+betas[0];
		double[][] Y_pred_local=new double[var_comp.length][];
		for(int r=0;r<var_comp.length;r++){
			Y_pred_local[r]=estimate(regions.reginal_kinships[r], 
					obs_pheno, regions.relational_matrixs[r], var_comp[r]);
			for(int k=0;k<Y_pred.length;k++)Y_pred[k]+=Y_pred_local[r][k]*betas[r+1];
		}
		this.pred_pheno= new Phenotype(obs_pheno.phe_id, relation_obs_pred.pred_ids.clone(), Y_pred);
		System.out.println("Finished estimating G-BLUP.");	
	}
	
	/*
	 * TODO!!!
	 */
	public void estimate_local_noglobal_reg(CausalRegions regions, double[] var_comp){
		double[] betas=train_betas(regions, var_comp);
		System.out.println("Started estimating G-BLUP using regression.");
		if(var_comp.length+1!=betas.length){
			System.out.println("var_comp.length+1!=betas.length");
			return;
		}
		double[] Y_pred=new double[this.pred_geno.sample_size]; //estimate(obs_matrix, obs_pheno, relation_obs_pred, h2);
		for(int k=0;k<Y_pred.length;k++)Y_pred[k]=betas[0];
		double[][] Y_pred_local=new double[var_comp.length][];
		for(int r=0;r<var_comp.length;r++){
			Y_pred_local[r]=estimate(regions.reginal_kinships[r], 
					obs_pheno, regions.relational_matrixs[r], var_comp[r]);
			for(int k=0;k<Y_pred.length;k++)Y_pred[k]+=Y_pred_local[r][k]*betas[r+1];
		}
		this.pred_pheno= new Phenotype(obs_pheno.phe_id, relation_obs_pred.pred_ids.clone(), Y_pred);
		System.out.println("Finished estimating G-BLUP.");	
	}
	
	/*
	 * The model implemented by the function is STATITICAL QUESTIONABLE!  
	 */
	public double[] train_betas(CausalRegions regions, double[] var_comp){
		System.out.println("Start training betas");
		double[] Y=this.obs_pheno.values.clone();
		int sample_size=this.obs_pheno.values.length;
		double[][] Xs_blup=new double[sample_size][1+var_comp.length];
		for(int i=0;i<sample_size;i++){
			System.out.println("~~~~ Round "+i+" ~~~~~~");
			Xs_blup[i][0]=estimate_left_one(regions.global_kinship, obs_pheno, this.obs_pheno.sample_ids[i], h2);
			for(int r=0;r<regions.reginal_kinships.length;r++){
				Xs_blup[i][1+r]=estimate_left_one(regions.reginal_kinships[r], 
						obs_pheno, this.obs_pheno.sample_ids[i], var_comp[r]);
			}
		}
		OLSMultipleLinearRegression reg1=new OLSMultipleLinearRegression();	
		reg1.newSampleData(Y, Xs_blup);		
		double[] betas=reg1.estimateRegressionParameters();
		if(betas.length!=var_comp.length+2){
			System.out.println("Regression error!");
			return null;
		}
		return betas;
	}
	
//	//temp
//	public static void train_betas(Phenotype pheno_full, KinshipMatrix kinship_full, String obs_ids[], double h2,
//			String output_file){
//		try{
//			String[] test_ids=the_rest_ids(pheno_full.sample_ids, obs_ids, kinship_full.ids2indexes);
//			System.out.println("obs_ids:"+obs_ids.length+"/test_ids:"+test_ids.length);
//			Phenotype obs_pheno= pheno_full.extract_sample_with_id_list(obs_ids);
//			Phenotype true_pred_pheno= pheno_full.extract_sample_with_id_list(test_ids);
//			KinshipMatrix obs_matrix=KinshipMatrix.submatrix(kinship_full, obs_ids);
//			RelationMatrix relation_obs_pred=RelationMatrix.submatrix(kinship_full, obs_ids, test_ids);
//			GBLUP gblup=new GBLUP(null, null, obs_pheno, obs_matrix, relation_obs_pred, h2);
//			System.out.println("finished submatrixes.");
//			gblup.estimate_global();
//			double[] corr=MultiPhenotype.correlation(true_pred_pheno, gblup.pred_pheno);
//			System.out.println(corr[0]+"\t"+corr[1]+"\t"+corr[2]+"\t"+corr[3]);
//		}catch(Exception e){e.printStackTrace();}
//	}
	
	public static double[] train_betas(Phenotype pheno_full, KinshipMatrix kinship_full, 
			CausalRegions regions, double[] var_comp, double h2, String[] train1_ids){
		if(regions.reginal_kinships.length!=var_comp.length){
			System.out.println("regions.reginal_kinships.length!=var_comp.length\n null returned.");
			return null;
		}
		String[] train2_ids=the_rest_ids(pheno_full.sample_ids, train1_ids, kinship_full.ids2indexes);
		System.out.println("1st_ids:"+train1_ids.length+"/2nd_ids:"+train2_ids.length);
		System.out.println("Start calculating BLUPs");
		Phenotype first_pheno= pheno_full.extract_sample_with_id_list(train1_ids);
		KinshipMatrix first_guys_kinship=KinshipMatrix.submatrix(kinship_full, train1_ids);		
		RelationMatrix first2second=RelationMatrix.submatrix(kinship_full, train1_ids, train2_ids);
		
		double[] Y=first_pheno.values.clone();
		//int sample_size1=first_pheno.values.length;
		int sample_size2=train2_ids.length;
		double[][] Xs_blup=new double[sample_size2][1+var_comp.length];
		double[] global=estimate(first_guys_kinship, first_pheno, first2second, h2);
		for(int i=0;i<sample_size2;i++)		Xs_blup[i][0]=global[i];
		for(int r=0;r<regions.reginal_kinships.length;r++){
			KinshipMatrix first_guys_local_kinship=KinshipMatrix.submatrix(regions.reginal_kinships[r], train1_ids);		
			RelationMatrix first2second_local=RelationMatrix.submatrix(regions.reginal_kinships[r], train1_ids, train2_ids);
			double[] local=estimate(first_guys_local_kinship, first_pheno, first2second_local, var_comp[r]);
			for(int i=0;i<sample_size2;i++)	Xs_blup[i][r+1]=local[i];
		}
		
		System.out.println("Started training betas");
		OLSMultipleLinearRegression reg1=new OLSMultipleLinearRegression();	
		Phenotype true_2nd_pheno= pheno_full.extract_sample_with_id_list(train2_ids);		
		reg1.newSampleData(true_2nd_pheno.values.clone(), Xs_blup);		
		double[] betas=reg1.estimateRegressionParameters();
		if(betas.length!=var_comp.length+2){
			System.out.println("Regression error!");
			return null;
		}
		System.out.println("Finished training betas");
		return betas;
	}
	
	/*
	 * 
	 */
	
	public void estimate_local_reg_(CausalRegions regions, double[] var_comp, VariantsDouble geno_full,
			int first_sample_size){
		String[] train1_ids=select_ids_random(first_sample_size, obs_pheno, geno_full);
		
		double[] betas=train_betas(this.obs_pheno, this.obs_matrix, regions, var_comp, this.h2, train1_ids);
		System.out.println("betas:");
		for(int i=0;i<betas.length;i++)System.out.println("\t"+betas[i]);
		System.out.println("Started estimating G-BLUP using regression.");
		if(var_comp.length+2!=betas.length){
			System.out.println("var_comp.length+2!=betas.length");
			return;
		}
		double[] Y_pred=estimate(obs_matrix, obs_pheno, relation_obs_pred, h2);
		for(int k=0;k<Y_pred.length;k++)Y_pred[k]=Y_pred[k]*betas[1]+betas[0];
		double[][] Y_pred_local=new double[var_comp.length][];
		for(int r=0;r<var_comp.length;r++){
			Y_pred_local[r]=estimate(regions.reginal_kinships[r], 
					obs_pheno, regions.relational_matrixs[r], var_comp[r]);
			for(int k=0;k<Y_pred.length;k++)Y_pred[k]+=Y_pred_local[r][k]*betas[r+2];
		}
		this.pred_pheno= new Phenotype(obs_pheno.phe_id, relation_obs_pred.pred_ids.clone(), Y_pred);
		System.out.println("Finished estimating G-BLUP.");	
	}
	
	public static double estimate_left_one(KinshipMatrix kinship, Phenotype pheno, String the_id, double var_comp){
		String[] the_left_one_id=new String[1];
		the_left_one_id[0]=the_id;
		String[] most_ids=new String[kinship.ids.length-1];
		int the_one_index=kinship.ids2indexes.get(the_id);
		for(int i=0;i<the_one_index;i++){
			most_ids[i]=kinship.ids[i];
		}for(int i=the_one_index;i<kinship.ids.length-1;i++){
			most_ids[i]=kinship.ids[i+1];
		}
		Phenotype most_guys_pheno=new Phenotype(pheno.phe_id, most_ids, pheno.extract_values_with_id_list(most_ids));
		KinshipMatrix most_guys_kinship=KinshipMatrix.submatrix(kinship, most_ids);
		RelationMatrix most2one=RelationMatrix.submatrix(kinship, most_ids, the_left_one_id);
		return estimate(most_guys_kinship, most_guys_pheno, most2one, var_comp)[0];
	}
	
	public static String[] the_rest_ids(String[] full, String[] subset, 
			HashMap<String, Integer> genotype_ids){
		ArrayList<String> full2=new ArrayList<String>();
		for(int k=0;k<full.length;k++){
			if(genotype_ids.containsKey(full[k]))
				full2.add(full[k]);
		}		
		String[] the_rest=new String[full2.size()-subset.length];
		
		HashSet<String> subset_table=new HashSet<String>();
		for(int k=0;k<subset.length;k++){
			if(subset_table.contains(subset[k])){
				System.out.println("Wrong ID: Alreday in! "+subset[k]);
			}
			subset_table.add(subset[k]);
		}
		System.out.println(full.length+"/"+full2.size()+"/"+subset.length+"/"+subset_table.size());
		int index=0;
		for(int i=0;i<full2.size();i++){
			if(!subset_table.contains(full2.get(i))){
				the_rest[index]=full2.get(i);
				index++;
			}
		}return the_rest;
	}
	
	public static String[] the_rest_ids(String[] full, String left){
		String[] the_rest=new String[full.length-1];
		int index=0;
		for(int i=0;i<full.length;i++){
			if(!full[i].equals(left))the_rest[index++]=full[i];
		}return the_rest;
	}
	
	public static void main(String[] args) {
		String eli_file="/Users/quanlong/Documents/projects2/predictions/RA_data/fromEli/aTNF/rdrw-SAMs_cept_mab_ifx_ada.pheno";
		String geno_folder="/Users/quanlong/Documents/projects2/predictions/RA_data/genotypes/";
		String pheno_folder="/Users/quanlong/Documents/projects2/predictions/RA_data/phenotypes/";
		String result_folder="/Users/quanlong/Documents/projects2/predictions/RA_data/result/";
		String out_file=result_folder+"global_blup005.csv";
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(out_file));
			bw.write("ID,pred_gen_facs,pred_clin_gen_facs\n");
			String[] types={"etanercept","infliximab","adalimumab","all"};
			double[] h2={0.05, 0.05, 0.05, 0.05};
			HashSet<String> non_NAN_ids=new HashSet<String>();
			for(int k=0;k<types.length;k++){	
				System.out.println("==============="+types[k]+"=================");
				Phenotype obs_pheno=new MultiPhenotype(pheno_folder+types[k]+".tsv").phenotypes[0];
				Phenotype comparison=new MultiPhenotype(pheno_folder+"eli/"+types[k]+".eli.test.tsv").phenotypes[0]; 
				VariantsDouble obs_geno=new VariantsDouble(geno_folder+types[k]+".train.hdf5");
				VariantsDouble pred_geno=new VariantsDouble(geno_folder+types[k]+".test.hdf5");
				KinshipMatrix obs_matrix=new KinshipMatrix(geno_folder+types[k]+".train.K.RRM");
				RelationMatrix relation_obs_pred=new RelationMatrix(geno_folder+types[k]+".train.test.K.RRM");
				double[] h2s={0.0001, 0.05,0.1,0.15, 0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1};
				//for(int i=0;i<h2s.length;i++){
				System.out.println(non_NAN_ids.size());
				GBLUP gblup=new GBLUP(obs_geno, pred_geno, obs_pheno, obs_matrix, relation_obs_pred, h2[k]);
				gblup.estimate_global();
				double[] corr=MultiPhenotype.correlation(comparison, gblup.pred_pheno);
				if(k!=types.length-1){
					for(int sample=0;sample<gblup.pred_pheno.values.length;sample++){
						non_NAN_ids.add(gblup.pred_pheno.sample_ids[sample]);
						bw.write(gblup.pred_pheno.sample_ids[sample]+","+
								200*gblup.pred_pheno.values[sample]+","+200*gblup.pred_pheno.values[sample]+"\n");
					}
					System.out.println("h2="+h2[k]+":  "+corr[0]+"\t"+corr[1]+"\t"+corr[4]);
				}
				else{
					for(int sample=0;sample<gblup.pred_pheno.values.length;sample++){
						if(!non_NAN_ids.contains(gblup.pred_pheno.sample_ids[sample])){
							bw.write(gblup.pred_pheno.sample_ids[sample]+","+
								200*gblup.pred_pheno.values[sample]+","+200*gblup.pred_pheno.values[sample]+"\n");
						}
					}
					System.out.println("h2="+h2[k]+":  "+corr[0]+"\t"+corr[1]+"\t"+corr[4]);
				}
					
				//}
			}bw.close();
		}catch(Exception e){e.printStackTrace();}	
	}

	/*
	 * chr starts from 1 in the parameter.
	 */
	public static double[] multiple_regression_prediction(ArrayList<int[]> markers, Phenotype pheno, 
			VariantsDouble geno, int training_sample_size){
		String[] training_ids=select_ids_random(training_sample_size, pheno, geno);
		String[] testing_ids=the_rest_ids(pheno.sample_ids, training_ids, geno.sample_id2index);
		double[][] snp_data_training=new double[training_ids.length][markers.size()];
		double[][] snp_data_testing=new double[markers.size()][testing_ids.length];
		
		for(int i=0;i<markers.size();i++){
			int chr=markers.get(i)[0]-1;
			int location=markers.get(i)[1];
			double[] the_snp=geno.load_one_variant_by_location(chr, location);
			//int the_moving_index_training=0;
			//int the_moving_index_testing=0;
			for(int sample_index=0;sample_index<training_ids.length;sample_index++){
				if(geno.sample_id2index.containsKey(training_ids[sample_index])){
					snp_data_training[sample_index][i]=the_snp[geno.sample_id2index.get(training_ids[sample_index])];
					//the_moving_index_training++;
				}else System.out.println("Training sample ID not exists, selection procedure wrong!");
			}			
			for(int sample_index=0;sample_index<testing_ids.length;sample_index++){
				if(geno.sample_id2index.containsKey(testing_ids[sample_index])){
					snp_data_testing[i][sample_index]=the_snp[geno.sample_id2index.get(testing_ids[sample_index])];
					//the_moving_index_testing++;
				}else System.out.println("Testing sample ID not exists, selection procedure wrong!");
			}
		}
		System.out.println("Started training betas");
		OLSMultipleLinearRegression reg1=new OLSMultipleLinearRegression();	
		Phenotype training_pheno= pheno.extract_sample_with_id_list(training_ids);		
		reg1.newSampleData(training_pheno.values.clone(), snp_data_training);		
		double[] betas=reg1.estimateRegressionParameters();
		if(betas.length!=markers.size()+1){
			System.out.println("Regression error!");
		}
		System.out.println("Finished training betas");
		
		double[] predicted_pheno=new double[testing_ids.length];
		Phenotype real_testing_pheno= pheno.extract_sample_with_id_list(testing_ids);
		for(int sample_index=0;sample_index<testing_ids.length;sample_index++){
			//int the_index_in_geno=geno.sample_id2index.get(testing_ids[sample_index]);
			predicted_pheno[sample_index]=betas[0];
			for(int i=0;i<markers.size();i++){
				predicted_pheno[sample_index]+=(betas[i+1]*snp_data_testing[i][sample_index]);				
			}
		}
		double[] corr=MultiPhenotype.correlation(new Phenotype(real_testing_pheno.phe_id, testing_ids, predicted_pheno), real_testing_pheno);
		return corr;
	}
	
	/*
	 * chr starts from 1 in the parameter.
	 * The same as above, except for training IDs are passed as parameter (instead of randomly generated).
	 */
	public static double[] multiple_regression_prediction(ArrayList<int[]> markers, Phenotype pheno, 
			VariantsDouble geno, String[] training_ids){
		//String[] training_ids=select_ids_random(training_sample_size, pheno, geno);
		String[] testing_ids=the_rest_ids(pheno.sample_ids, training_ids, geno.sample_id2index);
		double[][] snp_data_training=new double[training_ids.length][markers.size()];
		double[][] snp_data_testing=new double[markers.size()][testing_ids.length];
		
		for(int i=0;i<markers.size();i++){
			int chr=markers.get(i)[0]-1;
			int location=markers.get(i)[1];
			double[] the_snp=geno.load_one_variant_by_location(chr, location);
			//int the_moving_index_training=0;
			//int the_moving_index_testing=0;
			for(int sample_index=0;sample_index<training_ids.length;sample_index++){
				if(geno.sample_id2index.containsKey(training_ids[sample_index])){
					snp_data_training[sample_index][i]=the_snp[geno.sample_id2index.get(training_ids[sample_index])];
					//the_moving_index_training++;
				}else System.out.println("Training sample ID not exists, selection procedure wrong!");
			}			
			for(int sample_index=0;sample_index<testing_ids.length;sample_index++){
				if(geno.sample_id2index.containsKey(testing_ids[sample_index])){
					snp_data_testing[i][sample_index]=the_snp[geno.sample_id2index.get(testing_ids[sample_index])];
					//the_moving_index_testing++;
				}else System.out.println("Testing sample ID not exists, selection procedure wrong!");
			}
		}
		System.out.println("Started training betas");
		OLSMultipleLinearRegression reg1=new OLSMultipleLinearRegression();	
		Phenotype training_pheno= pheno.extract_sample_with_id_list(training_ids);		
		reg1.newSampleData(training_pheno.values.clone(), snp_data_training);		
		double[] betas=reg1.estimateRegressionParameters();
		if(betas.length!=markers.size()+1){
			System.out.println("Regression error!");
		}
		System.out.println("Finished training betas");
		
		double[] predicted_pheno=new double[testing_ids.length];
		Phenotype real_testing_pheno= pheno.extract_sample_with_id_list(testing_ids);
		for(int sample_index=0;sample_index<testing_ids.length;sample_index++){
			//int the_index_in_geno=geno.sample_id2index.get(testing_ids[sample_index]);
			predicted_pheno[sample_index]=betas[0];
			for(int i=0;i<markers.size();i++){
				predicted_pheno[sample_index]+=(betas[i+1]*snp_data_testing[i][sample_index]);				
			}
		}
		double[] corr=MultiPhenotype.correlation(new Phenotype(real_testing_pheno.phe_id, testing_ids, predicted_pheno), real_testing_pheno);
		return corr;
	}
}














