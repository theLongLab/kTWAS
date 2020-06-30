package simulations;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.time.LocalDateTime;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Random;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.time.temporal.TemporalAccessor;
import java.io.File;


import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;

import mixedmodel.EMMA;
import mixedmodel.EMMAX;
import mixedmodel.KinshipMatrix;
import mixedmodel.MultiPhenotype;
import mixedmodel.Phenotype;
import mixedmodel.VariantsDouble;
import myMathLib.Normalization;
import myMathLib.StatFuncs;
import myMathLib.Test;

public class kTWAS {	
	
	static String []  SNPs_Origin; 
	static int []  SNPs_Index; 
	static double  []  SNPs_Coeff; 
	static int Num_of_EN_SNPs;
	static int Num_of_eQTL_SNPs;
	
	public static String Gene_ENSG_ID;
	public static int  Gene_ENSG_Start;
	public static int  Gene_ENSG_End;
	public static String Gene_Chromosome;
	public static int [] EN_SNPs;
	public static double [] EN_SNP_Weights;
	public static String SNP_REF="A";
	public static String SNP_ALT="T";
	public static HashMap<String, Double> SNP_EN_Weights_Dic ;
	public static HashMap<String, String> EN_SNPs_RSID_Dic ;
	public static HashMap<String, String> EN_SNPs_Ref_Dic ;
	public static HashMap<String, String> EN_SNPs_Alt_Dic ;
	static double pheno_unit=1.0; 
	
	/*
	 * This method serves for two functions:
	 * 
	 * (1) Simulate phenotype using (multiple) SNPs
	 * (2) Simulate phenotype using (multiple) methylation sites (or regions)
	 * 
	 * As long as the input format is the same, i.e., our HDF5 file converted from the CSV-file, whether it contains
	 *  SNPs or methylation sites doesn't matter.  
	 * 
	 * The String causal_variants is in the form of ChrIndex_LocIndex;ChrIndex_LocIndex;...;ChrIndex_LocIndex
	 * Please note that all the indexes start from ZERO.
	 * 
	 * It supports the following mechanisms:
	 * (1) Additive  
	 * (2) Genetic Heterogeneity  
	 * (3) Effect size: for case/control, it is penetrance; for quantitative data, it is   
	 * 
	 */
	
//	public static void sim_phenotype_mSNPs_additive(
//			String input_genotype, String causal_variants, String output_pheno_file,
//			double heritability){
//		
//		int[][] causal_indexes=parse_causal_var_string(causal_variants);
//		int num_of_causal=causal_indexes.length;
//		// assign effect sizes
//		double[] effects=new double[num_of_causal];
//		NormalDistribution norm=new NormalDistribution();
//		for(int i=0;i<effects.length;i++){
//			effects[i]=norm.sample();
//		}
//		// assign phenotype based on genotype and effect size
//		
//		VariantsDouble genotype=new VariantsDouble(input_genotype);
//		
//		double[][] phenotypes=new double[1][genotype.sample_size]; // only one trait will be simulated
//		for(int var=0;var<num_of_causal;var++){
//			double[] genotypes=genotype.load_one_variant_by_index(causal_indexes[var][0], causal_indexes[var][1]);
//			for(int indi=0;indi<genotype.sample_size;indi++){					
//				phenotypes[0][indi]+=(effects[var]*genotypes[indi]);											
//			}
//		}	
//		
//		
//		// normalize the phenotype to ensure the specified heritability
//		Normalization.adding_random_component(phenotypes[0], heritability); 
////		normalization(phenotypes[0], heritability);
//		// write the phenotype to a file
//		write_simulated_pheno_file(phenotypes, output_pheno_file,  genotype);
//		System.out.println("Finsihed simulating phenotype with additive model and writing to "+output_pheno_file);	
//	}
	
	
	public static void sim_phenotype_mSNPs_xor(
			String input_genotype, String causal_variants, String output_pheno_file,
			double heritability, double ratio) throws IOException{
		
		BufferedWriter bw = new BufferedWriter(new FileWriter( output_pheno_file+".log", false));
		VariantsDouble genotype=new VariantsDouble(input_genotype, "csv");
		//Elastic Net
// EN + eQTL		
		int[][] causal_indexes= new int [Num_of_EN_SNPs+ Num_of_eQTL_SNPs][2];		

		for (int i=0;i< SNPs_Index.length;i++) {
			if ((SNPs_Origin[i].equals("EN")) || (SNPs_Origin[i].equals("eQTL"))) {
				causal_indexes[i][0]=0;
				causal_indexes[i][1]=SNPs_Index[i];
			}
		}
		double[][] phenotypes=new double[1][genotype.sample_size]; // only one trait will be simulated
		
		double[] genotypes_1=genotype.load_snps_by_index.get(
					(Integer.toString(causal_indexes[0][0])+":"+Integer.toString(causal_indexes[0][1])));
		double[] genotypes_2=genotype.load_snps_by_index.get(
				(Integer.toString(causal_indexes[1][0])+":"+Integer.toString(causal_indexes[1][1])));
		for(int indi=0;indi<genotype.sample_size;indi++){	
			if(Double.compare(genotypes_1[indi], genotypes_2[indi])!=0){
				phenotypes[0][indi]+=1.0;			
			}
			if (Math.abs(genotypes_1[indi]- genotypes_2[indi])==2.0) {
				phenotypes[0][indi]+=1.0;		
			}
		}
		Normalization.adding_random_component(phenotypes[0], heritability); 
		write_simulated_pheno_file(phenotypes, output_pheno_file,  genotype);
        bw.close();
		System.out.println("Finsihed simulating phenotype with xor model and writing to "
		+output_pheno_file);	
	}
	
	
	public static void sim_phenotype_mSNPs_and(
			String input_genotype, String causal_variants, String output_pheno_file,
			double heritability, double ratio) throws IOException{
		
		BufferedWriter bw = new BufferedWriter(new FileWriter( output_pheno_file+".log", false));
		VariantsDouble genotype=new VariantsDouble(input_genotype, "csv");
		//Elastic Net
// EN + eQTL		
		int[][] causal_indexes= new int [Num_of_EN_SNPs+ Num_of_eQTL_SNPs][2];		

		for (int i=0;i< SNPs_Index.length;i++) {
			if ((SNPs_Origin[i].equals("EN")) || (SNPs_Origin[i].equals("eQTL"))) {
				causal_indexes[i][0]=0;
				causal_indexes[i][1]=SNPs_Index[i];
			}
		}
		double[][] phenotypes=new double[1][genotype.sample_size]; // only one trait will be simulated
		
		double[] genotypes_1=genotype.load_snps_by_index.get(
					(Integer.toString(causal_indexes[0][0])+":"+Integer.toString(causal_indexes[0][1])));
		double[] genotypes_2=genotype.load_snps_by_index.get(
				(Integer.toString(causal_indexes[1][0])+":"+Integer.toString(causal_indexes[1][1])));
		for(int indi=0;indi<genotype.sample_size;indi++){	
			if((Double.compare(genotypes_1[indi], 0.0)!=0)  &&
					(Double.compare(genotypes_2[indi], 0.0)!=0)) {
				phenotypes[0][indi]+=1.0;			
			}
		
		}
		Normalization.adding_random_component(phenotypes[0], heritability); 
		write_simulated_pheno_file(phenotypes, output_pheno_file,  genotype);
        bw.close();
		System.out.println("Finsihed simulating phenotype with and model and writing to "
		+output_pheno_file);	
	}
	
	
	public static void sim_phenotype_mSNPs_or(
			String input_genotype, String causal_variants, String output_pheno_file,
			double heritability, double ratio) throws IOException{
		
		BufferedWriter bw = new BufferedWriter(new FileWriter( output_pheno_file+".log", false));
		VariantsDouble genotype=new VariantsDouble(input_genotype, "csv");
		//Elastic Net
// EN + eQTL		
		int[][] causal_indexes= new int [Num_of_EN_SNPs+ Num_of_eQTL_SNPs][2];		

		for (int i=0;i< SNPs_Index.length;i++) {
			if ((SNPs_Origin[i].equals("EN")) || (SNPs_Origin[i].equals("eQTL"))) {
				causal_indexes[i][0]=0;
				causal_indexes[i][1]=SNPs_Index[i];
			}
		}
		double[][] phenotypes=new double[1][genotype.sample_size]; // only one trait will be simulated
		
		double[] genotypes_1=genotype.load_snps_by_index.get(
					(Integer.toString(causal_indexes[0][0])+":"+Integer.toString(causal_indexes[0][1])));
		double[] genotypes_2=genotype.load_snps_by_index.get(
				(Integer.toString(causal_indexes[1][0])+":"+Integer.toString(causal_indexes[1][1])));
		for(int indi=0;indi<genotype.sample_size;indi++){	
			if((Double.compare(genotypes_1[indi], 0.0)!=0)  ||
					(Double.compare(genotypes_2[indi], 0.0)!=0)) {
				phenotypes[0][indi]+=1.0;			
			}
		
		}
		Normalization.adding_random_component(phenotypes[0], heritability); 
		write_simulated_pheno_file(phenotypes, output_pheno_file,  genotype);
        bw.close();
		System.out.println("Finsihed simulating phenotype or additive model and writing to "
		+output_pheno_file);	
	}
	
	
	
	
	

	public static void sim_phenotype_mSNPs_additive(
			String input_genotype, String causal_variants, String output_pheno_file,
			double heritability, double ratio) throws IOException{
		
		BufferedWriter bw = new BufferedWriter(new FileWriter( output_pheno_file+".log", false));
		VariantsDouble genotype=new VariantsDouble(input_genotype, "csv");
		//Elastic Net
		int [][] en_causal_indexes = new int [Num_of_EN_SNPs][2];
		int num_of_en_causal = Num_of_EN_SNPs;
		double[] en_effects=new double[num_of_en_causal];
		int index =0;
		for (int i=0;i< SNPs_Index.length;i++) {
			if (SNPs_Origin[i].equals("EN")) {
				en_causal_indexes[index][0]=0;
				en_causal_indexes[index][1]=SNPs_Index[i];
				en_effects[index]=SNPs_Coeff[i];
				index++;
			}
		}
		
		for(int i=0;i<en_effects.length;i++){
			System.out.println("EN EFFECT:\t"+en_effects[i]);
			bw.write(  Double.toString(en_effects[i])+"\n");
		}
		bw.write("=============================================\n");
		double[][] en_phenotypes=new double[1][genotype.sample_size]; // only one trait will be simulated
		for(int var=0;var<num_of_en_causal;var++){
			double[] en_genotypes=genotype.load_snps_by_index.get(
					(Integer.toString(en_causal_indexes[var][0])+":"+Integer.toString(en_causal_indexes[var][1])));
			for(int indi=0;indi<genotype.sample_size;indi++){
				en_phenotypes[0][indi]+=(en_effects[var]*en_genotypes[indi]);		
			}
		}
		
// eQTL
		int [][] eqtl_causal_indexes = new int [Num_of_eQTL_SNPs][2];
		int num_of_eqtl_causal = Num_of_eQTL_SNPs;
		double[] eqtl_effects=new double[num_of_eqtl_causal];
		index =0;
		for (int i=0;i< SNPs_Index.length;i++) {
			if (SNPs_Origin[i].equals("eQTL")) {
				eqtl_causal_indexes[index][0]=0;
				eqtl_causal_indexes[index][1]=SNPs_Index[i];
				eqtl_effects[index]=0.1;
				index++;
			}
		}
		
		System.out.println("-----------------------------");
		for(int i=0;i<eqtl_effects.length;i++){
			System.out.println(eqtl_effects[i]);
		}
		double[][] eqtl_phenotypes=new double[1][genotype.sample_size]; // only one trait will be simulated
		for(int var=0;var<num_of_eqtl_causal;var++){
			System.out.println(Integer.toString(eqtl_causal_indexes[var][0])+":"+Integer.toString(eqtl_causal_indexes[var][1]));
			double[] eqtl_genotypes=genotype.load_snps_by_index.get(
					(Integer.toString(eqtl_causal_indexes[var][0])+":"+Integer.toString(eqtl_causal_indexes[var][1])));
			for(int indi=0;indi<genotype.sample_size;indi++){
				eqtl_phenotypes[0][indi]+=(eqtl_effects[var]*eqtl_genotypes[indi]);	
			}
		}
		
		double en_var= StatFuncs.popVar_NaN(en_phenotypes[0]);
		double eqtl_var= StatFuncs.popVar_NaN(eqtl_phenotypes[0]);
		double coeff=0.0;
		if (ratio<Double.MIN_VALUE) {
			coeff= 1.0;
			for (int i=0;i< en_effects.length;i++) {
				en_effects[i]=0.0;
			}
		}else {
			coeff = Math.sqrt(en_var/ eqtl_var)* Math.sqrt( (1-ratio)/ratio ) ;
		}
		
		for (int i=0;i< eqtl_effects.length;i++) {
			eqtl_effects[i]= eqtl_effects[i]* coeff;
		}
		
		System.out.println("==============================");
		for(int i=0;i<eqtl_effects.length;i++){
			System.out.println(eqtl_effects[i]);
		}		
		
// EN + eQTL		
		int[][] causal_indexes= new int [Num_of_EN_SNPs+ Num_of_eQTL_SNPs][2];
		int num_of_causal=causal_indexes.length; 
		double[] effects=new double[num_of_causal];
		
		int eqtl_index =0;
		int en_index =0;
		for (int i=0;i< SNPs_Index.length;i++) {
			if (SNPs_Origin[i].equals("EN")) {
				causal_indexes[i][0]=0;
				causal_indexes[i][1]=SNPs_Index[i];
				effects[i]=en_effects[en_index];
				en_index++;
			} else if (SNPs_Origin[i].equals("eQTL")) {
				causal_indexes[i][0]=0;
				causal_indexes[i][1]=SNPs_Index[i];
				effects[i]=eqtl_effects[eqtl_index];
				eqtl_index++;
			}
		}
		System.out.println("==============================");
		System.out.println("------------------------------");
		for(int i=0;i<effects.length;i++){
			System.out.println(causal_indexes[i][1]+"\t"+effects[i]);
			bw.write(  Integer.toString(causal_indexes[i][1])+"\t"+Double.toString(effects[i]) +"\n");
		}		
		
		
		double[][] phenotypes=new double[1][genotype.sample_size]; // only one trait will be simulated
		for(int var=0;var<num_of_causal;var++){
			double[] genotypes=genotype.load_snps_by_index.get(
					(Integer.toString(causal_indexes[var][0])+":"+Integer.toString(causal_indexes[var][1])));
			for(int indi=0;indi<genotype.sample_size;indi++){					
				phenotypes[0][indi]+=(effects[var]*genotypes[indi]);											
			}
		}		
		
		Normalization.adding_random_component(phenotypes[0], heritability); 
//		normalization(phenotypes[0], heritability);
		write_simulated_pheno_file(phenotypes, output_pheno_file,  genotype);
		
		
     
        bw.close();
		
		
		System.out.println("Finsihed simulating phenotype with additive model and writing to "
		+output_pheno_file);	

//		NormalDistribution en_norm=new NormalDistribution();
//		for(int i=0;i<en_effects.length;i++){
//			en_effects[i]=en_norm.sample();
//		}

		
		
//		int[][] causal_indexes=parse_causal_var_string(causal_variants);//  
//		int num_of_causal=causal_indexes.length;
		// assign effect sizes
//		double[] effects=new double[num_of_causal];
//		NormalDistribution norm=new NormalDistribution();
//		for(int i=0;i<effects.length;i++){
//			effects[i]=norm.sample();
//		}
		
		// assign phenotype based on genotype and effect size
//		VariantsDouble genotype=new VariantsDouble(input_genotype);
//		double[][] phenotypes=new double[1][genotype.sample_size]; // only one trait will be simulated
//		for(int var=0;var<num_of_causal;var++){
//			double[] genotypes=genotype.load_one_variant_by_index(causal_indexes[var][0], causal_indexes[var][1]);
//			for(int indi=0;indi<genotype.sample_size;indi++){					
//				phenotypes[0][indi]+=(effects[var]*genotypes[indi]);											
//			}
//		}	
		// normalize the phenotype to ensure the specified heritability
//		Normalization.adding_random_component(phenotypes[0], heritability); 
//		normalization(phenotypes[0], heritability);
		// write the phenotype to a file
//		write_simulated_pheno_file(phenotypes, output_pheno_file,  genotype);
//		System.out.println("Finsihed simulating phenotype with additive model and writing to "+output_pheno_file);	
	}
	
	

	public static  void sim_phenotype_mSNPs_xor(
			String input_genotype, String causal_variants, String output_pheno_file,
			double heritability){
		
		int[][] causal_indexes=parse_causal_var_string(causal_variants);
		int num_of_causal=causal_indexes.length;
		// assign effect sizes
		double[] effects=new double[num_of_causal];
		NormalDistribution norm=new NormalDistribution();
		for(int i=0;i<effects.length;i++){
			effects[i]= 1.0;
		}
		// assign phenotype based on genotype and effect size
				
		VariantsDouble genotype=new VariantsDouble(input_genotype); 
		
		double[][] phenotypes=new double[1][genotype.sample_size]; // only one trait will be simulated
		
		for(int indi=0;indi<genotype.sample_size;indi++){	
			double [] loci_1_genotype =     genotype.load_one_variant_by_index
					(causal_indexes[0][0], causal_indexes[0][1]);
			double [] loci_2_genotype =     genotype.load_one_variant_by_index
					(causal_indexes[1][0], causal_indexes[1][1]);
			if (Double.compare(loci_1_genotype[indi] , loci_2_genotype[indi])!= 0 ) {
				phenotypes[0][indi]=  pheno_unit;	
			} 
		}
		
//		
//		for(int var=0;var<num_of_causal;var++){
//			double[] genotypes=genotype.load_one_variant_by_index
//					(causal_indexes[var][0], causal_indexes[var][1]);
//			
//			for(int indi=0;indi<genotype.sample_size;indi++){					
//				phenotypes[0][indi]+=(effects[var]*genotypes[indi]);											
//			}
//		}				
		
		
		// normalize the phenotype to ensure the specified heritability
		Normalization.adding_random_component(phenotypes[0], heritability); 
		// write the phenotype to a file
		write_simulated_pheno_file(phenotypes, output_pheno_file,  genotype);
		System.out.println("Finsihed simulating phenotype with additive model and writing to "+output_pheno_file);
	}
	
	
	
	/*
	 * Test the association of (simulated) causal variants and assess their p-value by permutations
	 * EMMAX algorithm is used to control population structure 
	 * Note that it applies to both DNA=> Phenotype and Methylation=>Phenotype. 
	 */
	public static void test_causal_emmax(String genotype_hdf5_file, String kinship_file, String pheno_file, 
			String causal_variants, int permutation, String output_file){
		try{
			int[][] causal_indexes=parse_causal_var_string(causal_variants);
			int num_of_causal=causal_indexes.length;
			Phenotype phenotype=(new MultiPhenotype(pheno_file)).phenotypes[0];
			VariantsDouble genotype=new VariantsDouble(genotype_hdf5_file);
			EMMA emma=new EMMA(phenotype, genotype,  new KinshipMatrix(kinship_file).kinship_matrix.getData());	
			int[] indexes={}, chrs={};
			emma.REMLE(emma.phenotype.values, emma.cbind(chrs, indexes), null);		
			double[][] decomposed_array = emma.reml_decompositioned_array();				
			emma.phenotype.generate_new_Y_by_multiplying(decomposed_array); 
			double[] intsept=new double[emma.sample_size];
			for(int i=0;i<emma.sample_size;i++){
				for(int j=0;j<emma.sample_size;j++){
					intsept[i]+=decomposed_array[i][j];
				}
			}	
			BufferedWriter bw=new BufferedWriter(new FileWriter(output_file));
			bw.write("ChrIndex\tLocIndex\tPermutedPValue\n");
			for(int var=0;var<num_of_causal;var++){
				double[] X_ori=genotype.load_one_variant_by_index(causal_indexes[var][0], causal_indexes[var][1]);
				double[][] Xs_after=new double[emma.sample_size][2];						
				for(int i=0;i<emma.sample_size;i++){
					Xs_after[i][0]=intsept[i];
					for(int j=0;j<emma.sample_size;j++){
						Xs_after[i][1]+=(decomposed_array[i][j]*X_ori[j]);
					}
				}
				// test the variant.
				double[] Y=emma.phenotype.new_Y_4emmax.clone();		
				double real_P=regression_P(Y, Xs_after, true);
				// use permutation to assess p-value
				double[] permuted_Ps=new double[permutation];
				for(int i_permute=0;i_permute<permutation;i_permute++){
					double[] Yp=emma.phenotype.new_Y_4emmax.clone();
					Test.randomPermute(Yp);
					permuted_Ps[i_permute]=regression_P(Yp, Xs_after, true);
				}
				Arrays.sort(permuted_Ps);
				double rank=0.0;
				for(int i_permute=0;i_permute<permutation;i_permute++){
					if(permuted_Ps[i_permute]<real_P)rank=rank+1;
					else break;
				}bw.write(causal_indexes[var][0]+"\t"+causal_indexes[var][1]+"\t"+(rank/permutation)+"\n");
			}bw.close();	
		}catch(Exception e){e.printStackTrace();}
	}
	
	/*
	 * Test the association of (simulated) causal variants and assess their p-value by permutations
	 * Linear Model is used to control population structure 
	 * Note that it applies to both DNA=> Phenotype and Methylation=>Phenotype. 
	 */
	public static void test_causal_lm(String genotype_hdf5_file, String pheno_file, 
			String causal_variants, int permutation, String output_file){
		try{
			int[][] causal_indexes=parse_causal_var_string(causal_variants);
			int num_of_causal=causal_indexes.length;
			Phenotype phenotype=(new MultiPhenotype(pheno_file)).phenotypes[0];
			VariantsDouble genotype=new VariantsDouble(genotype_hdf5_file);
			BufferedWriter bw=new BufferedWriter(new FileWriter(output_file));
			bw.write("ChrIndex\tLocIndex\tPermutedPValue\n");
			for(int var=0;var<num_of_causal;var++){
				double[] X_ori=genotype.load_one_variant_by_index(causal_indexes[var][0], causal_indexes[var][1]);
				double[][] Xs_after=new double[X_ori.length][2];						
				for(int i=0;i<X_ori.length;i++){
					Xs_after[i][0]=1.0;
					Xs_after[i][1]=(X_ori[i]);					
				}
				// test the variant.
				double[] Y=phenotype.values.clone();		
				double real_P=regression_P(Y, Xs_after, true);
				// use permutation to assess p-value
				double[] permuted_Ps=new double[permutation];
				for(int i_permute=0;i_permute<permutation;i_permute++){
					double[] Yp=phenotype.values.clone();
					Test.randomPermute(Yp);
					permuted_Ps[i_permute]=regression_P(Yp, Xs_after, true);
				}
				Arrays.sort(permuted_Ps);
				double rank=0.0;
				for(int i_permute=0;i_permute<permutation;i_permute++){
					if(permuted_Ps[i_permute]<real_P)rank=rank+1;
					else break;
				}bw.write(causal_indexes[var][0]+"\t"+causal_indexes[var][1]+"\t"+(rank/permutation)+"\n");
			}bw.close();	
		}catch(Exception e){e.printStackTrace();}
	}
	
	public static String get_causal_variants(String snp_info_fil, int num_en_snps, int num_gtex_snps)
			throws IOException{
		
		Num_of_EN_SNPs= num_en_snps;
		Num_of_eQTL_SNPs= num_gtex_snps;
		
		ArrayList<Integer> en_snps_index = new ArrayList<Integer>();
		ArrayList<Double> en_snps_coeff = new ArrayList<Double>();
		ArrayList<Integer> gtex_snps_index = new ArrayList<Integer>();
		ArrayList<Double> snps_coeff = new ArrayList<Double>();
		
    	BufferedReader br = new BufferedReader(new FileReader(snp_info_fil));
        String line = br.readLine();
        while (line.startsWith("#")) {
            line = br.readLine();
        }
//        line = br.readLine();
        int index_count = 0;
        while (line != null) {
        	String[] tmp_arr = line.split("\t");
        	if (tmp_arr[1].equals("EN")) {
        		en_snps_index.add( index_count); 
        		en_snps_coeff.add( Double.parseDouble(tmp_arr[2])); 
        		
        	} else if (tmp_arr[1].equals("eQTL")) {
        		gtex_snps_index.add(index_count);
        	}
        	snps_coeff.add(Double.parseDouble(tmp_arr[2]));
            line = br.readLine();
            index_count++ ;
        }
        br.close();
        
        
        
        if ((en_snps_index.size()< num_en_snps) || (gtex_snps_index.size()< num_gtex_snps)) {
        	System.out.println("There exists not enough SNPs in EN or eQTL for the file "+ snp_info_fil);
        	System.out.println("The number of EN SNPs:\t"+ num_en_snps+"; the number of eEQTL snps:\t"+ num_gtex_snps); 
            System.exit(0);
        } 
        
//      Collections.shuffle(en_snps_index,new Random(en_snps_index.size()));//Random
        
        for (int i=0; i<en_snps_index.size(); i++) {
        	for (int j=i; j<en_snps_index.size(); j++) {
        		if (Math.abs(en_snps_coeff.get(i)) <  Math.abs(en_snps_coeff.get(j))) {
        			Double tmp_double = en_snps_coeff.get(i);
        			en_snps_coeff.set(i, en_snps_coeff.get(j));
        			en_snps_coeff.set(j, tmp_double);
        			int tmp_int = en_snps_index.get(i);
        			en_snps_index.set(i, en_snps_index.get(j));
        			en_snps_index.set(j, tmp_int);
        		}
        	}
        }
        
        
        Collections.shuffle(gtex_snps_index,new Random(gtex_snps_index.size()));
        
        int [] sel_snps_index = new int[num_en_snps+num_gtex_snps];
        String[] sel_snps_origin = new String[num_en_snps+num_gtex_snps];
        
        
        for (int i=0; i<num_en_snps; i++) {
        	sel_snps_index[i]=en_snps_index.get(i);
        	sel_snps_origin[i]="EN";    	
        }
        
        for (int i=0; i<num_gtex_snps; i++) {
        	sel_snps_index[num_en_snps+i]=gtex_snps_index.get(i);
        	sel_snps_origin[num_en_snps+i]="eQTL";
        }
        
        for (int i=0; i<sel_snps_index.length; i++) {
        	for (int j=i; j<sel_snps_index.length; j++) {
        		if (sel_snps_index[i]> sel_snps_index[j]) {
        			int tmp_int = sel_snps_index[i];
        			sel_snps_index[i]= sel_snps_index[j];
        			sel_snps_index[j]= tmp_int;
        			String tmp_str = sel_snps_origin[i];
        			sel_snps_origin[i]= sel_snps_origin[j];
        			sel_snps_origin[j]= tmp_str;
        		}
        	}
        }
        
        
        String causal_variants="";
        for (int i=0;i< sel_snps_index.length;i++) {
        	causal_variants=causal_variants+"0_"+ Integer.toString(sel_snps_index[i])+";";
        }	
        System.out.println("causal_variants:\t"+causal_variants );

        SNPs_Index = sel_snps_index.clone();
        SNPs_Origin= sel_snps_origin.clone();
        SNPs_Coeff= new double [num_en_snps+num_gtex_snps];
        for (int i=0;i< SNPs_Coeff.length;i++) {
        	SNPs_Coeff[i]= snps_coeff.get(SNPs_Index[i]);
        }
        
        for (int i=0;i< SNPs_Index.length;i++) {
        	System.out.println( SNPs_Index[i]+"\t"+SNPs_Origin[i]+"\t"+ SNPs_Coeff[i] );
        }
		return causal_variants;
	}
	
	public static double regression_P(double[] Y, double[][] Xs, boolean no_intercept){
		OLSMultipleLinearRegression reg1=new OLSMultipleLinearRegression();	
		if(no_intercept)reg1.setNoIntercept(true);
		reg1.newSampleData(Y, Xs);	
		double pvalue=myMathLib.StatFuncs.multi_reg_pvalues(reg1.estimateRegressionParameters(), 
				reg1.estimateRegressionParametersStandardErrors(), Y.length)[1];
		return pvalue; 
	}

	/*
	 * Estimate the overall R2 of multiple mSNPs that explains the variance of a methylated site (Could be a CpG site or a DMR)
	 * Only LM is provided: this is because that the function that simulates phenotype based on a given heritability only uses LM (not EMMAX).
	 */
	public static void mSNP_aggregated_R2(String genotype_hdf5_file, String pheno_file, 
			String causal_variants, String output_file){
		try{
			int[][] causal_indexes=parse_causal_var_string(causal_variants);
			int num_of_causal=causal_indexes.length;
			Phenotype phenotype=(new MultiPhenotype(pheno_file)).phenotypes[0];
			VariantsDouble genotype=new VariantsDouble(genotype_hdf5_file);
			BufferedWriter bw=new BufferedWriter(new FileWriter(output_file));
			bw.write("The Variance component explained by the input mSNPs ("+causal_variants+") is\n");
			double[] Y=phenotype.values.clone();	
			double[][] Xs=new double[Y.length][num_of_causal];		
			for(int var=0;var<num_of_causal;var++){
				double[] X_ori=genotype.load_one_variant_by_index(causal_indexes[var][0], causal_indexes[var][1]);								
				for(int i=0;i<Y.length;i++)	Xs[i][var]=X_ori[i];				
			}
			// test the regression of all the mSNPs.
			OLSMultipleLinearRegression reg1=new OLSMultipleLinearRegression();	
			reg1.newSampleData(Y, Xs);	
			bw.write(reg1.calculateAdjustedRSquared()+"\n");			
			bw.close();	
		}catch(Exception e){e.printStackTrace();}	
	}
	
	public static void write_simulated_pheno_file(double[][] phenos, String file, VariantsDouble genotype){
		try{
			BufferedWriter bw= new BufferedWriter(new FileWriter(file));
			bw.write("ID");
			for(int i=0;i<phenos.length;i++){
				bw.write("\tSim_Pheno_"+i);
			}bw.write("\n");
			for(int indi=0;indi<genotype.sample_size;indi++){
				bw.write(genotype.sample_ids[indi]);
				for(int i=0;i<phenos.length;i++){
					bw.write("\t"+phenos[i][indi]);
				}bw.write("\n");
			}bw.close();
		}catch(Exception e){e.printStackTrace();}
	}
	
	/*
	 * adding the variance of error term
	 * change will be made to the double[] original 
	 */
	
	public static void normalization(double[] original, double heritability){
		NormalDistribution norm=new NormalDistribution();
		double sigma_g2=StatUtils.populationVariance(original);
		double sigma_e=Math.sqrt(sigma_g2*(1-heritability)/heritability); // heritability=sigma_g2/(sigma_g2+sigma_e2)
		for(int indi=0;indi<original.length;indi++){
			original[indi]+=(sigma_e*norm.sample());
		}
	}
	
	/*
	 * parse the causal_variants string into an index array.
	 */
	public static int[][] parse_causal_var_string(String causal_variants){
		String[] indexes_array=causal_variants.split(";");
		int num_of_causal=indexes_array.length;
		int[][] casual_indexes=new int[num_of_causal][2];
		for(int var=0;var<num_of_causal;var++){
			casual_indexes[var][0]=Integer.parseInt(indexes_array[var].split("_")[0]);  // ChrIndex
			casual_indexes[var][1]=Integer.parseInt(indexes_array[var].split("_")[1]);	// LocIndex
		}
		return casual_indexes;
	}	
	
	
	public static void assign_en_gtex_snps(String csv_file, String predixcan_en_snps_file, 
			String gtex_eqtl_snps_file, String output_file, String gene_name) throws IOException{
		
		ArrayList<String> chr_array = new ArrayList<String>();
		ArrayList<String> pos_array = new ArrayList<String>();

		BufferedReader br_csv = new BufferedReader(new FileReader(csv_file));
		String line = "";
		
		HashMap<String, String> origin_dict= new HashMap<String, String>();
		HashMap<String, Double> coeff_dict= new HashMap<String, Double>();
		
		while ((line = br_csv.readLine()) != null) {
			if ((!line.substring(0, 1).equals("#")) && (!line.substring(0, 1).equals("C"))) {
				line= line.replace("\n", "").replace("\r", "");
				String[] tmp_arr = line.split(",");
				chr_array.add(tmp_arr[0]);
				pos_array.add(tmp_arr[1]);
				origin_dict.put(tmp_arr[0]+":"+ tmp_arr[1], "NULL");
				coeff_dict.put(tmp_arr[0]+":"+ tmp_arr[1], 0.0);
			}
		}
		br_csv.close();
		
		int num_snps_gtex=0;
		BufferedReader br_gtex = new BufferedReader(new FileReader(gtex_eqtl_snps_file));
		while ((line = br_gtex.readLine()) != null) {
			line= line.replace("\n", "").replace("\r", "");
			String[] tmp_arr = line.split("\t");
			if (tmp_arr[0].equals(gene_name)) {
				if (origin_dict.containsKey(tmp_arr[1])) {
					origin_dict.put(tmp_arr[1], "eQTL");
					num_snps_gtex++;
				}
			}
		}
		br_gtex.close();
		
		int num_snps_en=0;
		BufferedReader br_en = new BufferedReader(new FileReader(predixcan_en_snps_file));
		while ((line = br_en.readLine()) != null) {
			line= line.replace("\n", "").replace("\r", "");
			String[] tmp_arr = line.split("\t");
			if (tmp_arr[0].equals(gene_name)) {
				if (origin_dict.containsKey(tmp_arr[1]+":"+tmp_arr[2])) {
					origin_dict.put(tmp_arr[1]+":"+ tmp_arr[2], "EN");
					coeff_dict.put(tmp_arr[1]+":"+ tmp_arr[2], Double.parseDouble(tmp_arr[3])); 
					num_snps_en++;
				}
			}
		}
		br_en.close();
		
		BufferedWriter bw = new BufferedWriter(new FileWriter( output_file, false));
		String ss="#CHR:POS\tORIGIN\tCoeff";
		bw.write(ss+ "\n");
		for (int i=0;i< chr_array.size();i++) {
			ss= chr_array.get(i)+":"+pos_array.get(i)+"\t"+origin_dict.get(chr_array.get(i)+":"+pos_array.get(i))+
					"\t"+ coeff_dict.get(chr_array.get(i)+":"+pos_array.get(i));
			bw.write(ss+ "\n");
		}
		bw.close();
		
		System.out.println("Finsih assigning PrediXcan SNPs and GTEx eQTL to csv file:\t"+output_file);
		System.out.println("Num of SNPs from eQTL:\t"+num_snps_gtex+".\nNum of SNPs from EN:\t"+num_snps_en);
		System.out.println("=================");
	}	
	
	public static void convert_csv_2_plink(String csv_file, String pheno_file, 
			 String plink_output_folder, int column_index) 
					throws IOException, InterruptedException{

		BufferedReader br_csv = new BufferedReader(new FileReader(csv_file));
		String line = br_csv.readLine();//header
        String[] temp = line.split(",");
        int sample_size = temp.length - 2;
        String [] sample_ids = new String[sample_size];
        String [] phenotypes = new String[sample_size];
        HashMap<String, String> id_dict = new HashMap<String, String>();
        
        BufferedReader br_pheno = new BufferedReader(new FileReader(pheno_file));
        line = br_pheno.readLine();//header
        int line_index =0;
        while ((line = br_pheno.readLine()) != null) {
        	line= line.replace("\n", "").replace("\r", "");
			String[] tmp_arr = line.split("\t");
			id_dict.put(tmp_arr[0], tmp_arr[column_index-1]); 
			line_index++;
        }
        br_pheno.close();
        
        for (int k = 0; k < sample_size; k++) sample_ids[k] = temp[k + 2];
        
        String tfam_file= plink_output_folder+ "/plink.tfam";
        BufferedWriter bw_tfam = new BufferedWriter(new FileWriter( tfam_file, false));
        for (int i=0;i< sample_size;i++) {
        	String ss= sample_ids[i]+"\t"+sample_ids[i]+"\t0\t0\t0\t"+ id_dict.get(sample_ids[i]);
        	bw_tfam.write(ss+"\n");
        }
        bw_tfam.close();
        
        line_index=0;
        String tped_file= plink_output_folder+ "/plink.tped";
        BufferedWriter bw_tped = new BufferedWriter(new FileWriter( tped_file, false));
        String setid_file= plink_output_folder+ "/plink.setid";
        BufferedWriter bw_setid = new BufferedWriter(new FileWriter( setid_file, false));
        String en_weight_file= plink_output_folder+ "/plink.en.weights";
		BufferedWriter bw_en_weight = new BufferedWriter(new FileWriter( en_weight_file, false));
		int snp_index = 1;
		while ((line = br_csv.readLine()) != null) {
			line= line.replace("\n", "").replace("\r", "").replace(" ", "");
			String[] tmp_arr = line.split(",");
			int pos= Integer.parseInt(tmp_arr[1]);
			if ((tmp_arr[0].equals(Gene_Chromosome)) && (pos>=Gene_ENSG_Start) && (pos<=Gene_ENSG_End)) {
				String rs=  "rs_"+Integer.toString(snp_index);
				String ref = SNP_REF;
				String alt = SNP_ALT;
				if (EN_SNPs_RSID_Dic.containsKey(tmp_arr[0]+":"+tmp_arr[1])) {
					rs= EN_SNPs_RSID_Dic.get(tmp_arr[0]+":"+tmp_arr[1]); 
				}
				if (EN_SNPs_Ref_Dic.containsKey(tmp_arr[0]+":"+tmp_arr[1])) {
					ref= EN_SNPs_Ref_Dic.get(tmp_arr[0]+":"+tmp_arr[1]); 
				}
				if (EN_SNPs_Alt_Dic.containsKey(tmp_arr[0]+":"+tmp_arr[1])) {
					alt= EN_SNPs_Alt_Dic.get(tmp_arr[0]+":"+tmp_arr[1]); 
				}
				
				bw_setid.write("SET\t"+ rs+"\n");
				String ss =tmp_arr[0]+"\t"+ rs +"\t0\t"+tmp_arr[1];
				for (int i=2;i< tmp_arr.length;i++) {
					String plink_genotype= ref+"\t"+ref; 
					if (tmp_arr[i].equals("1")) {
						plink_genotype= ref+"\t"+alt; 
					} else if (tmp_arr[i].equals("2")) {
						plink_genotype= alt+"\t"+alt; 
					}
					ss = ss+"\t"+plink_genotype;
				}
				bw_tped.write(ss+"\n");
				ss= rs+"\t";
				if (SNP_EN_Weights_Dic.containsKey(tmp_arr[0]+":"+tmp_arr[1])) {
					ss= ss+Double.toString(  SNP_EN_Weights_Dic.get(tmp_arr[0]+":"+tmp_arr[1]) ) ; 
				} else {
					ss= ss+"0.0";
				}
				bw_en_weight.write(ss+"\n");
				snp_index++;
			}
		}
		br_csv.close();
		bw_en_weight.close();
		bw_setid.close();
		bw_tped.close();
	}
	
	public static void select_plink(String plink_file, String pheno_file, 
			 String plink_output_folder, int column_index) 
					throws IOException, InterruptedException{

		BufferedReader br_plink = new BufferedReader(new FileReader(plink_file));
		String line = br_plink.readLine();
		
       String[] temp = line.split("\t");
       int sample_size = (temp.length - 4)/2;
       String [] sample_ids = new String[sample_size];
       String [] phenotypes = new String[sample_size];
       HashMap<String, String> id_dict = new HashMap<String, String>();
       br_plink.close();
       
       BufferedReader br_pheno = new BufferedReader(new FileReader(pheno_file));
//       line = br_pheno.readLine();//header
       int line_index =0;
       while ((line = br_pheno.readLine()) != null) {
       	line= line.replace("\n", "").replace("\r", "");
			String[] tmp_arr = line.split("\t");
			id_dict.put(tmp_arr[0], tmp_arr[column_index-1]); 
			sample_ids[line_index]= tmp_arr[0]; 
			line_index++;
       }
       br_pheno.close();
       
//       for (int k = 0; k < sample_size; k++) sample_ids[k] = temp[k + 2];
       
       String tfam_file= plink_output_folder+ "/plink.tfam";
       BufferedWriter bw_tfam = new BufferedWriter(new FileWriter( tfam_file, false));
       for (int i=0;i< sample_size;i++) {
       	String ss= sample_ids[i]+"\t"+sample_ids[i]+"\t0\t0\t0\t"+ id_dict.get(sample_ids[i]);
       	bw_tfam.write(ss+"\n");
       }
       bw_tfam.close();
       
       line_index=0;
       String tped_file= plink_output_folder+ "/plink.tped";
       BufferedWriter bw_tped = new BufferedWriter(new FileWriter( tped_file, false));
       String setid_file= plink_output_folder+ "/plink.setid";
       BufferedWriter bw_setid = new BufferedWriter(new FileWriter( setid_file, false));
       String en_weight_file= plink_output_folder+ "/plink.en.weights";
		BufferedWriter bw_en_weight = new BufferedWriter(new FileWriter( en_weight_file, false));
		int snp_index = 1;
		BufferedReader br_plink2 = new BufferedReader(new FileReader(plink_file));
		while ((line = br_plink2.readLine()) != null) {
			line= line.replace("\n", "").replace("\r", "").replace(" ", "");
			String[] tmp_arr = line.split("\t");
			int pos= Integer.parseInt(tmp_arr[3]);
			if ((tmp_arr[0].equals(Gene_Chromosome)) && (pos>=Gene_ENSG_Start) && (pos<=Gene_ENSG_End)) {
				String rs=  "rs_"+Integer.toString(snp_index);
				String ref = SNP_REF;
				String alt = SNP_ALT;
				if (EN_SNPs_RSID_Dic.containsKey(tmp_arr[0]+":"+tmp_arr[3])) {
					rs= EN_SNPs_RSID_Dic.get(tmp_arr[0]+":"+tmp_arr[3]); 
				}
				if (EN_SNPs_Ref_Dic.containsKey(tmp_arr[0]+":"+tmp_arr[3])) {
					ref= EN_SNPs_Ref_Dic.get(tmp_arr[0]+":"+tmp_arr[3]); 
				}
				if (EN_SNPs_Alt_Dic.containsKey(tmp_arr[0]+":"+tmp_arr[3])) {
					alt= EN_SNPs_Alt_Dic.get(tmp_arr[0]+":"+tmp_arr[3]); 
				}
				
				bw_setid.write("SET\t"+ rs+"\n");
				String ss =tmp_arr[0]+"\t"+ rs +"\t0\t"+tmp_arr[3];
				HashMap<String, Integer> allele_dict= new HashMap<String, Integer>();
				for (int i=4;i< tmp_arr.length;i++) {
					ss=ss+"\t"+ tmp_arr[i];
				}
//					if (allele_dict.containsKey( tmp_arr[i] )) {
//						allele_dict.put(tmp_arr[i], allele_dict.get(tmp_arr[i])+1);
//					} else {
//						allele_dict.put(tmp_arr[i], 1);
//					}
//				}
//				String ref_allele ="";
//				int max_allele_count=0;
//				for (String key : allele_dict.keySet()) {
//					if (allele_dict.get(key)>= max_allele_count) {
//						
//					}
//					
//				}
				
//				for (int i=0;i< (tmp_arr.length-4)/2;i++) {
//					
//					String plink_genotype= ref+"\t"+ref; 
//					if (tmp_arr[4+2*i].equals("1")) {
//						plink_genotype= ref+"\t"+alt; 
//					} else if (tmp_arr[i].equals("2")) {
//						plink_genotype= alt+"\t"+alt; 
//					}
//					ss = ss+"\t"+plink_genotype;
//				}
				
				bw_tped.write(ss+"\n");
				ss= rs+"\t";
				if (SNP_EN_Weights_Dic.containsKey(tmp_arr[0]+":"+tmp_arr[3])) {
					ss= ss+Double.toString(  SNP_EN_Weights_Dic.get(tmp_arr[0]+":"+tmp_arr[3]) ) ; 
				} else {
					ss= ss+"0.0";
				}
				bw_en_weight.write(ss+"\n");
				snp_index++;
			}
		}
		br_plink2.close();
		bw_en_weight.close();
		bw_setid.close();
		bw_tped.close();
	}
	
	
	
	public static void get_en_info(String info_file) 
					throws IOException, InterruptedException{
		String line= null;
		boolean OK= false;
		SNP_EN_Weights_Dic= new HashMap<String,Double>();
		EN_SNPs_RSID_Dic= new HashMap<String,String>();
		EN_SNPs_Ref_Dic= new HashMap<String,String>();
		EN_SNPs_Alt_Dic= new HashMap<String,String>();
		BufferedReader br_info = new BufferedReader(new FileReader(info_file));
		while ((line = br_info.readLine()) != null) {
			 line= line.replace("\n", "").replace("\r", "");
			 String[] tmp_arr = line.split("\t");
			 if (tmp_arr[0].equals(Gene_ENSG_ID)) {
				 OK= true;
				 Gene_Chromosome = tmp_arr[1];
				 Gene_ENSG_Start = Integer.parseInt(tmp_arr[2]);
				 Gene_ENSG_End = Integer.parseInt(tmp_arr[3]);
				 EN_SNPs = new int [tmp_arr.length - 4];
				 EN_SNP_Weights = new double [tmp_arr.length - 4];
				 for (int i= 4; i< tmp_arr.length;i++) {
					 String[] tmp_arr2 = tmp_arr[i].split("_");
					 EN_SNPs[i-4]= Integer.parseInt(tmp_arr2[1]);
					 EN_SNP_Weights[i-4]= Double.parseDouble(tmp_arr2[4]) ; 
					 EN_SNPs_RSID_Dic.put(tmp_arr[1]+":"+tmp_arr2[1], 
							 tmp_arr2[0]);  
					 SNP_EN_Weights_Dic.put(tmp_arr[1]+":"+tmp_arr2[1], 
							 Double.parseDouble(tmp_arr2[4]) );  
					 EN_SNPs_Ref_Dic.put(tmp_arr[1]+":"+tmp_arr2[1], 
							 tmp_arr2[2] );
					 EN_SNPs_Alt_Dic.put(tmp_arr[1]+":"+tmp_arr2[1], 
							 tmp_arr2[3] );
				 }
			 }
		}
		br_info.close();
		if (!OK) {
			System.out.println("Can not find the gene:\t"+Gene_ENSG_ID+" in the "+ info_file+ "!" );
			System.exit(0);
		}
		
	}
	
	public static void make_bed(String out_folder, String plink) 
			throws IOException, InterruptedException{
		
		ProcessBuilder CMDLine = new ProcessBuilder(plink,
				"--tfile", 
				out_folder+"plink", 
                "--make-bed",
                "--noweb",
                "--out",
                out_folder+"plink"
               );
            Process CMDProcess = CMDLine.start();
            
            BufferedReader br = new BufferedReader(new InputStreamReader(CMDProcess.getInputStream()));
            String line;
            while ((line = br.readLine()) != null) {
            	line= line.replace("\n", "");
            	System.out.println(line);
            }
            CMDProcess.waitFor();
	}
	
	public static void recode(String out_folder, String plink) 
			throws IOException, InterruptedException{
		
		ProcessBuilder CMDLine = new ProcessBuilder(plink,
				"--tfile", 
				out_folder+"plink", 
                "--recode",
                "--noweb",
                "--out",
                out_folder+"plink"
               );
            Process CMDProcess = CMDLine.start();
            
            BufferedReader br = new BufferedReader(new InputStreamReader(CMDProcess.getInputStream()));
            String line;
            while ((line = br.readLine()) != null) {
            	line= line.replace("\n", "");
            	System.out.println(line);
            }
            CMDProcess.waitFor();
	}
	
	
	public static void get_result(String out_folder) 
			throws IOException, InterruptedException{
		String line ="";
		double p_value_noweigths= 1.0;
		double p_value_enweights= 1.0;
		BufferedReader br_info = new BufferedReader(new FileReader(out_folder+"kTWAS.result"));
		while ((line = br_info.readLine()) != null) {
			 line= line.replace("\n", "").replace("\r", "");
			 String[] tmp_arr = line.split("\t");
			 if (tmp_arr[0].equals("NO_weights")) {
				 p_value_noweigths = Double.parseDouble(tmp_arr[1]);
				 
			 } else if (tmp_arr[0].equals("EN_weights")) {
				 p_value_enweights = Double.parseDouble(tmp_arr[1]);
			 }
		}
		br_info.close();
		
//		System.out.println("The p-value of SKAT on GENE "+ Gene_ENSG_ID+":\t"+ 
//				Double.toString(p_value_noweigths));
//		System.out.println("-------------------------------------");
		System.out.println("The p-value of kTWAS on GENE "+ Gene_ENSG_ID+":\t"+ 
				Double.toString(p_value_enweights));
		System.out.println("=====================================");
		
		
	}
	
	public static void write_R(String out_folder , String tpye) 
			throws IOException, InterruptedException{
		String file= out_folder+ "skat.R";
        BufferedWriter bw= new BufferedWriter(new FileWriter( file, false));
        String ss= "library(SKAT)";
        bw.write(ss+"\n");
        ss= "Generate_SSD_SetID(\""+out_folder+"plink.bed\",\""+out_folder
        		+ "plink.bim\",\""+out_folder+"plink.fam\",\""+out_folder+"plink.setid\","
        		+ "\""+out_folder+"plink.ssd\",\""+out_folder+"plink.info\")";
        bw.write(ss+"\n");
        
        ss="FAM <- Read_Plink_FAM(\""+out_folder+"plink.fam\", Is.binary = FALSE)";
        if (tpye.equals("binary")) {
        	ss="FAM <- Read_Plink_FAM(\""+out_folder+"plink.fam\", Is.binary = TRUE)";
        }
        bw.write(ss+"\n");
        ss ="y <- FAM$Phenotype";
        bw.write(ss+"\n");
//        continuous|binary"
        ss="obj <- SKAT_Null_Model(y ~ 1, out_type = \"C\")";
        if (tpye.equals("binary")) {
        	ss="obj <- SKAT_Null_Model(y ~ 1, out_type = \"D\")";
        }
        bw.write(ss+"\n");
        ss= "SSD.INFO <- Open_SSD(\""+out_folder+"plink.ssd\", \""+out_folder+"plink.info\")";
        bw.write(ss+"\n");
        ss="res_no_weight <-SKAT.SSD.All(SSD.INFO, obj )";
        bw.write(ss+"\n");
//        ss="ss= paste(\"NO_weights\",toString(res_no_weight$results$P.value),  sep = \"\\t\")";
//        bw.write(ss+"\n");
//        ss="write(ss ,\""+out_folder+"kTWAS.result\", append=TRUE)";
//        bw.write(ss+"\n");
        ss = "W_en <- Read_SNP_WeightFile(\""+out_folder+"plink.en.weights\")";
        bw.write(ss+"\n");
        ss="res_en_weight <-SKAT.SSD.All(SSD.INFO, obj,   obj.SNPWeight = W_en)";
        bw.write(ss+"\n");
        ss="ss= paste(\"EN_weights\",toString(res_en_weight$results$P.value),  sep = \"\\t\")";
        bw.write(ss+"\n");
        ss="write(ss ,\""+out_folder+"kTWAS.result\", append=TRUE)";
        bw.write(ss+"\n");
        bw.close();
	}
	
	public static void run_Skat(String out_folder, String R) 
			throws IOException, InterruptedException{
		
		ProcessBuilder CMDLine = new ProcessBuilder(R,
				out_folder+"skat.R"
               );
            Process CMDProcess = CMDLine.start();
            
            BufferedReader br = new BufferedReader(new InputStreamReader(CMDProcess.getInputStream()));
            String line;
            while ((line = br.readLine()) != null) {
            	line= line.replace("\n", "");
            	System.out.println(line);
            }
            
            BufferedReader br_error = new BufferedReader(new InputStreamReader(CMDProcess.getErrorStream()));
            while ((line = br_error.readLine()) != null) {
            	line= line.replace("\n", "");
            	System.out.println(line);
            }
            
            CMDProcess.waitFor();
	}
	
	
	public static void csv_2_plink(String csv_file, String pheno_file, 
			String snp_info_file, String rs_file, String plink_output_folder) 
					throws IOException, InterruptedException{
		
		
		
		BufferedReader br_csv = new BufferedReader(new FileReader(csv_file));
		String line = br_csv.readLine();//header
        String[] temp = line.split(",");
        int sample_size = temp.length - 2;
        String [] sample_ids = new String[sample_size];
        String [] phenotypes = new String[sample_size];
        
        BufferedReader br_pheno = new BufferedReader(new FileReader(pheno_file));
        line = br_pheno.readLine();//header
        int line_index =0;
        while ((line = br_pheno.readLine()) != null) {
        	line= line.replace("\n", "").replace("\r", "");
			String[] tmp_arr = line.split("\t");
			phenotypes[line_index] = tmp_arr[1];
			line_index++;
        }
        br_pheno.close();
        
        HashMap<String, String> rs_dict= new HashMap<String, String>();
        HashMap<String, String> ref_allele_dict= new HashMap<String, String>();
        HashMap<String, String> alt_allele_dict= new HashMap<String, String>();
        BufferedReader br_rs = new BufferedReader(new FileReader(rs_file));
        line = br_rs.readLine();//header
        while ((line = br_rs.readLine()) != null) {
        	line= line.replace("\n", "").replace("\r", "");
        	String[] tmp_arr = line.split("\t");
        	String key = tmp_arr[1];
        	rs_dict.put(key, tmp_arr[4]);
        	ref_allele_dict.put(key, tmp_arr[2]);
        	alt_allele_dict.put(key, tmp_arr[3]);
        }
        br_rs.close();
        
        for (int k = 0; k < sample_size; k++) sample_ids[k] = temp[k + 2];
        
        new File(plink_output_folder+"/metaskat_en").mkdir();
        new File(plink_output_folder+"/metaskat_eqtl").mkdir();
        
        String en_tfam_file= plink_output_folder+ "/metaskat_en/plink.tfam";
        BufferedWriter bw_en_tfam = new BufferedWriter(new FileWriter( en_tfam_file, false));
        for (int i=0;i< sample_size;i++) {
        	String ss= sample_ids[i]+"\t"+sample_ids[i]+"\t0\t0\t0\t"+ phenotypes[i];
        	bw_en_tfam.write(ss+"\n");
        	ss= sample_ids[i]+"\t"+sample_ids[i]+"\t"+phenotypes[i];
        }
        bw_en_tfam.close();
        
        String eqtl_tfam_file= plink_output_folder+ "/metaskat_eqtl/plink.tfam";
        BufferedWriter bw_eqtl_tfam = new BufferedWriter(new FileWriter( eqtl_tfam_file, false));
        for (int i=0;i< sample_size;i++) {
        	String ss= sample_ids[i]+"\t"+sample_ids[i]+"\t0\t0\t0\t"+ phenotypes[i];
        	bw_eqtl_tfam.write(ss+"\n");
        	ss= sample_ids[i]+"\t"+sample_ids[i]+"\t"+phenotypes[i];
        }
        bw_eqtl_tfam.close();
        
        
        
        
        
        String tfam_file= plink_output_folder+ "/plink.tfam";
        BufferedWriter bw_tfam = new BufferedWriter(new FileWriter( tfam_file, false));
        BufferedWriter bw_predixcan_pheno = new BufferedWriter(new FileWriter( 
        		plink_output_folder+ "/plink.predixcan.pheno", false));
        bw_predixcan_pheno.write("FID\tIID\tpheno\n");
        for (int i=0;i< sample_size;i++) {
        	String ss= sample_ids[i]+"\t"+sample_ids[i]+"\t0\t0\t0\t"+ phenotypes[i];
        	bw_tfam.write(ss+"\n");
        	ss= sample_ids[i]+"\t"+sample_ids[i]+"\t"+phenotypes[i];
        	bw_predixcan_pheno.write(ss +"\n");
        }
        bw_tfam.close();
        bw_predixcan_pheno.close();
        
        line_index=0;
        String tped_file= plink_output_folder+ "/plink.tped";
        BufferedWriter bw_tped = new BufferedWriter(new FileWriter( tped_file, false));
        String setid_file= plink_output_folder+ "/plink.setid";
        BufferedWriter bw_setid = new BufferedWriter(new FileWriter( setid_file, false));
		while ((line = br_csv.readLine()) != null) {
			line= line.replace("\n", "").replace("\r", "");
			String[] tmp_arr = line.split(",");
			String ss =tmp_arr[0]+"\t"+rs_dict.get(tmp_arr[1]) +"\t0\t"+tmp_arr[1];
			bw_setid.write("SET\t"+ rs_dict.get(tmp_arr[1])+"\n");
			line_index++;
			for (int i=2;i< tmp_arr.length;i++) {
				String ref = ref_allele_dict.get( tmp_arr[1]);
				String alt= alt_allele_dict.get( tmp_arr[1]);
				String geno= ref+"\t"+ref; 
				if (tmp_arr[i].equals("1")) {
					geno= ref+"\t"+alt; 
				} else if (tmp_arr[i].equals("2")) {
					geno= alt+"\t"+alt; 
				}
				
				ss = ss+"\t"+geno;
			}
			bw_tped.write(ss+"\n");
		}
		bw_setid.close();
		bw_tped.close();
		br_csv.close();
		
		BufferedReader br_info = new BufferedReader(new FileReader(snp_info_file));
		String en_weight_file= plink_output_folder+ "/plink.en.weights";
		BufferedWriter bw_en_weight = new BufferedWriter(new FileWriter( en_weight_file, false));
		
		String eqtl_weight_file = plink_output_folder+ "/plink.eqtl.weights";
		BufferedWriter bw_eqtl_weight = new BufferedWriter(new FileWriter( eqtl_weight_file, false));
		
		HashMap<String, String> en_rs_dict= new HashMap<String, String>();
		HashMap<String, String> eqtl_rs_dict= new HashMap<String, String>();
		
		while ((line = br_info.readLine()) != null) {
			line= line.replace("\n", "").replace("\r", "");
			if (!line.substring(0, 1).equals("#")) {
				String[] tmp_arr = line.split("\t");
				String pos = tmp_arr[0].split(":")[1];
				String ss_en ="";
				String ss_eqtl ="";
				if (tmp_arr[1].equals("EN")) {
					ss_en = rs_dict.get(pos)+"\t"+tmp_arr[2];
					ss_eqtl= rs_dict.get(pos)+"\t0.0";
					en_rs_dict.put(rs_dict.get(pos), tmp_arr[2] );
					
				} else if (tmp_arr[1].equals("eQTL")) { 
					ss_en = rs_dict.get(pos)+"\t0.0";
//					ss_eqtl= rs_dict.get(pos)+"\t1.0";
					ss_eqtl= rs_dict.get(pos)+"\t" + tmp_arr[2];
					eqtl_rs_dict.put(rs_dict.get(pos), tmp_arr[2] );
					
				}	else {
					ss_en = rs_dict.get(pos)+"\t0.0";
					ss_eqtl= rs_dict.get(pos)+"\t0.0";
				}
				bw_en_weight.write(ss_en+"\n");
				bw_eqtl_weight.write(ss_eqtl+"\n");
			}
		}
		bw_en_weight.close();
		bw_eqtl_weight.close();
		br_info.close();
		
		
		line_index=0;
        String en_tped_file= plink_output_folder+ "/metaskat_en/plink.tped";
        BufferedWriter bw_en_tped = new BufferedWriter(new FileWriter( en_tped_file, false));
        String en_setid_file= plink_output_folder+ "/metaskat_en/plink.setid";
        BufferedWriter bw_en_setid = new BufferedWriter(new FileWriter( en_setid_file, false));
        
        String eqtl_tped_file= plink_output_folder+ "/metaskat_eqtl/plink.tped";
        BufferedWriter bw_eqtl_tped = new BufferedWriter(new FileWriter( eqtl_tped_file, false));
        String eqtl_setid_file= plink_output_folder+ "/metaskat_eqtl/plink.setid";
        BufferedWriter bw_eqtl_setid = new BufferedWriter(new FileWriter( eqtl_setid_file, false));
        
        BufferedReader br_csv2 = new BufferedReader(new FileReader(csv_file));
		line = br_csv2.readLine();//header
		while ((line = br_csv2.readLine()) != null) {
			line= line.replace("\n", "").replace("\r", "");
			String[] tmp_arr = line.split(",");
			String ss =tmp_arr[0]+"\t"+rs_dict.get(tmp_arr[1]) +"\t0\t"+tmp_arr[1];
			if (en_rs_dict.containsKey(   rs_dict.get(tmp_arr[1]) )) {
				bw_en_setid.write("SET\t"+ rs_dict.get(tmp_arr[1])+"\n");
			}
			if (eqtl_rs_dict.containsKey(   rs_dict.get(tmp_arr[1]) )) {
				bw_eqtl_setid.write("SET\t"+ rs_dict.get(tmp_arr[1])+"\n");
			}
			
			for (int i=2;i< tmp_arr.length;i++) {
				String ref = ref_allele_dict.get( tmp_arr[1]);
				String alt= alt_allele_dict.get( tmp_arr[1]);
				String geno= ref+"\t"+ref; 
				if (tmp_arr[i].equals("1")) {
					geno= ref+"\t"+alt; 
				} else if (tmp_arr[i].equals("2")) {
					geno= alt+"\t"+alt; 
				}
				ss = ss+"\t"+geno;
			}
			if (eqtl_rs_dict.containsKey(   rs_dict.get(tmp_arr[1]) )) {
				bw_eqtl_tped.write(ss+"\n");
			}
			if (en_rs_dict.containsKey(   rs_dict.get(tmp_arr[1]) )) {
				bw_en_tped.write(ss+"\n");
			}
			
		}
		bw_en_tped.close();
		bw_en_setid.close();
		bw_eqtl_tped.close();
		bw_eqtl_setid.close();
		br_csv2.close();
	}

	public static void main(String[] args) throws IOException, InterruptedException {
		if(args.length==0){
			System.out.println("k-SKAT. Usage:\n"
					+ "-format CSV|PLINK\n"
					+ "-input_genotype INPUT_GENOTYPE_FILE\n"
					+ "-input_phenotype INPUT_PHENOTYPE_FILE\n"
					+ "-input_phenotype_column INPUT_PHENOTYPE_COLUMN_START:2"
					+ "-en_info_path INPUT_ELASTICNET_INFORMATION_FILE\n"
					+ "-gene INPUT_ENSEMBL_GENE_ID\n"
					+ "-output_folder OUTPUT_FOLDER_PATH\n");
			System.exit(0);
		}
		String function=args[0];
		
		if(args[0].equals("pvalue-lm")){
			if(args.length==1){
				System.out.println("Conduct association test using linear model and calcualte the P-value using permutation.\n"
						+ "Usage: <genotype_file> <causal_variants_string> <phenotype_file> <#permutations> <output_file>\n"
						+ "\tNote: \n\tgenotype_file must be an HDF5 file. \n"
						+ "\tcausal_variants_string is in the format of: ChrIndex_LocIndex;ChrIndex_LocIndex;...;ChrIndex_LocIndex. \n"
						+ "\t#permutations can be an integer between 200 and 2000.");
			}else{
				String genotype_hdf5_file=args[1];
				String causal_variants=args[2];
				String pheno_file=args[3];
				int permutation=Integer.parseInt(args[4]);
				String output_file=args[5];
				test_causal_lm(genotype_hdf5_file, pheno_file, causal_variants, permutation, output_file);
			}
			
		}else if(args[0].equals("pvalue-emmax")){
			if(args.length==1){
				System.out.println("Conduct association test using linear mixed model (EMMAX algorithm) and calcualte the P-value using permutation.\n"
						+ "Usage: <genotype_file> <causal_variants_string> <phenotype_file> <#permutations> <kinship_file> <output_file>\n"
						+ "\tNote: \n\tgenotype_file must be an HDF5 file. \n"
						+ "\tcausal_variants_string is in the format of: ChrIndex_LocIndex;ChrIndex_LocIndex;...;ChrIndex_LocIndex. \n"
						+ "\t#permutations can be an integer between 200 and 2000.\n"
						+ "\tkinship_file can be generated using jawamix5.jar.\n");
			}else{
				String genotype_hdf5_file=args[1];
				String causal_variants=args[2];
				String pheno_file=args[3];
				int permutation=Integer.parseInt(args[4]);
				String kinship_file=args[5];
				String output_file=args[6];
				test_causal_emmax(genotype_hdf5_file, kinship_file, pheno_file, causal_variants, 
						permutation, output_file);
			}
			
		}else if(args[0].equals("sim-phenotype")){
			if(args.length==1){
				System.out.println("Simulate phenotype based on genotype (of DNAm or DNA).\n"
						+ "Usage: <genotype_file> <SNP_information_file>\n"
						+ "<output_phenotype_file> <output_SNP_coefficient_file> <heritability>\n"
						+ "<model> <num of SNPs form predixcan elastic net> <num of SNPs form GTEx eQTL>\n"
						+ "<the ratio between the num of SNPs form predixcan elastic net and the num of SNPs form GTEx eQTL>\n"
						+ "\tNote:\n\t genotype_file must be an csv file. \n"
						+ "\tcausal_variants_string is in the format of: ChrIndex_LocIndex;ChrIndex_LocIndex;...;ChrIndex_LocIndex. \n"
						+ "\theritability is a number between 0 and 1.\n"
						+ "model: additive, xor, and, compensation");
				
			}else{
				String input_genotype=args[1];
				String SNP_Info_File=args[2];
				String output_pheno_file=args[3];
				double heritability=Double.parseDouble(args[4]);
				String model = args[5];
				int Num_SNPs_EN= Integer.parseInt(args[6]);
				int Num_SNPs_eQTL= Integer.parseInt(args[7]);
				double Ratio_EN_eQTL=Double.parseDouble(args[8]);
				
				String causal_variants= get_causal_variants (SNP_Info_File, Num_SNPs_EN, Num_SNPs_eQTL);
				
				if (model.equals("additive")){
					sim_phenotype_mSNPs_additive(input_genotype, causal_variants, 
							output_pheno_file, heritability, Ratio_EN_eQTL);
				} else if (model.equals("xor")){
					sim_phenotype_mSNPs_xor(input_genotype, causal_variants, 
							output_pheno_file, heritability, Ratio_EN_eQTL);
				}else if (model.equals("and")){
					sim_phenotype_mSNPs_and(input_genotype, causal_variants, 
							output_pheno_file, heritability, Ratio_EN_eQTL);
				}else if (model.equals("or")){
					sim_phenotype_mSNPs_or(input_genotype, causal_variants, 
							output_pheno_file, heritability, Ratio_EN_eQTL);
				}
			}
		}else if(args[0].equals("mSNPsR2")){
			if(args.length==1){
				System.out.println("Calcualte the total methylation variance explained by multiple mSNPs.\n"
						+ "Usage: <genotype_file> <causal_variants_string> <phenotype_file> <output_file> \n"
						+ "\tNote:\n\t genotype_file must be an HDF5 file. \n"
						+ "\tcausal_variants_string is in the format of: ChrIndex_LocIndex;ChrIndex_LocIndex;...;ChrIndex_LocIndex. \n"
						+ "");
			}else{
				String genotype_hdf5_file= args[1];
				String causal_variants= args[2];
				String pheno_file= args[3];
				String output_file= args[4];
				mSNP_aggregated_R2(genotype_hdf5_file, pheno_file, causal_variants, output_file);
			}
		}else if(args[0].equals("assign_eqtl")){
			if(args.length < 5){
				System.out.println("Assign EN snps as well as GTEx snps to the snps in the csv file.\n"
						+ "Usage: <csv_file> <predixcan_elastic_net_snps_file> <gtex_eqtl_file>\n"
						+ "<output_file> <gene_name>\n");
			} else {
				String csv_file= args[1];
				String predixcan_en_snps_file =  args[2];
				String gtex_eqtl_snps_file= args[3];
				String gene_name= args[4];
				String output_file= args[5];
				assign_en_gtex_snps(csv_file, predixcan_en_snps_file, gtex_eqtl_snps_file, 
						output_file, gene_name);
			}
			
		}else if(args[0].equals("csv_2_plink")){
			if(args.length < 4){
				System.out.println("Covert csv format to plink format.\n"
						+ "Usage: <csv_file> <phenotype file>  "
						+ "<SNP_information_file> <snp_origin> <plink output file folder>\n");
			}else {
				String csv_file= args[1];
				String pheno_file =  args[2];
				String snp_info_file= args[3];
				String rs_file= args[4];
				String plink_output_folder= args[5];
				csv_2_plink(csv_file, pheno_file,  snp_info_file, 
						rs_file, plink_output_folder);
			}
			
		} else if(args[0].equals("kTWAS")){
			final DateTimeFormatter dtf = DateTimeFormatter.ofPattern("yyyy/MM/dd HH:mm:ss");
			if(args.length < 14){
				System.out.println("kTWAS. Usage:\n"
						+ "-format CSV|PLINK\n"
						+ "-input_genotype INPUT_GENOTYPE_FILE\n"
						+ "-input_phenotype INPUT_PHENOTYPE_FILE\n"
						+ "-input_phenotype_column INPUT_PHENOTYPE_COLUMN_START:2"
						+ "-input_phenotype_type PHENOTYPE_TYPE: continuous|binary"
						+ "-en_info_path INPUT_ELASTICNET_INFORMATION_FILE\n"
						+ "-gene INPUT_ENSEMBL_GENE_ID\n"
						+ "-plink PLINK_BINARY_FILE_PATH\n"
						+ "-Rscript RSCRIPT_BINARY_FILE_PATH\n"
						+ "-output_folder OUTPUT_FOLDER_PATH\n");
				System.exit(0);
			}else {
				String genotype_file =null, phenotype_file = null, format= null, 
						info_file =null, output_folder =null, column = null, gene_id = null,
						rscript= null, plink =null, phenotype_type= null;
				
				for(int k=1;k<args.length;k++){
					if(args[k].startsWith("-")){
						if(args[k].equals("-format")) format=args[k+1];
						else if(args[k].equals("-input_genotype")) genotype_file=args[k+1];
						else if(args[k].equals("-input_phenotype")) phenotype_file=args[k+1];
						else if(args[k].equals("-input_phenotype_column")) column=args[k+1];
						else if(args[k].equals("-en_info_path")) info_file=args[k+1];
						else if(args[k].equals("-gene")) gene_id=args[k+1];
						else if(args[k].equals("-output_folder")) output_folder=args[k+1];
						else if(args[k].equals("-Rscript")) rscript=args[k+1];
						else if(args[k].equals("-plink")) plink=args[k+1];
						else if(args[k].equals("-input_phenotype_type")) phenotype_type=args[k+1];
					}
				}
				
				if(genotype_file==null){
					System.out.println("Please input the genotype file!");
					System.exit(0);
				} if(format==null){
					System.out.println("Please input genotype format: csv/ plink!");
					System.exit(0);
				}else if (phenotype_file==null){
					System.out.println("Please input the phenotype file!");
					System.exit(0);
				} else if (column==null){
					System.out.println("Please input the phenotype column number!");
					System.exit(0);
				} else if (info_file==null){
					System.out.println("Please input the gene elastic net information file!");
					System.exit(0);
				} else if (output_folder==null){
					System.out.println("Please input the output folder path!");
					System.exit(0);
				} else if (rscript==null){
					System.out.println("Please input Rscript binrary file path!");
					System.exit(0);
				} else if (gene_id==null){
					System.out.println("Please input ensembl gene id!");
					System.exit(0);
				}else if (plink==null){
					System.out.println("Please input plink binrary file path!");
					System.exit(0);
				}else if (phenotype_type==null){
					System.out.println("Please input phenotype type: continuous or binary!");
					System.exit(0);
				}
				
				output_folder=output_folder+"/";
				new File(output_folder).mkdir();
				Gene_ENSG_ID = gene_id;
				get_en_info( info_file);
				if (format.equals("csv")) {
					System.out.println("Selecting and Converting CSV file to PLINK Format.\t"+
							dtf.format(LocalDateTime.now()));
					convert_csv_2_plink( genotype_file,  phenotype_file, 
							  output_folder, Integer.parseInt(column)) ;
				} else if (format.equals("plink")) {
					select_plink( genotype_file,  phenotype_file, 
							  output_folder, Integer.parseInt(column)) ;
					System.out.println("Selecting the variants within the gene region.\t"+
							dtf.format(LocalDateTime.now()));
				} else {
					System.out.println("The input genotype format should be csv or plink format!");
					System.exit(0);
				}
				
				System.out.println("Selecting and Converting CSV file to PLINK Format Finished.\t"+
						dtf.format(LocalDateTime.now()));
				System.out.println("=====================================");
				
				System.out.println("Making Bed and Recoding the PLINK files.\t"+
						dtf.format(LocalDateTime.now()));
				make_bed(output_folder, plink);
				recode(output_folder, plink);
				System.out.println("Making Bed and Recoding the PLINK files Finished.\t"+
						dtf.format(LocalDateTime.now()));
				System.out.println("=====================================");
				
				System.out.println("Running SKAT with Elastic Net Weights.\t"+
						dtf.format(LocalDateTime.now()));
				write_R(output_folder, phenotype_type);
				run_Skat(output_folder, rscript);
				System.out.println("Running SKAT with Elastic Net Weights Finished.\t"+
						dtf.format(LocalDateTime.now()));
				System.out.println("=====================================");
				get_result(output_folder);
			}
		} else{
			System.out.println("The function "+function+" doesn't exist. Typo?\n");
		}
	}
}
