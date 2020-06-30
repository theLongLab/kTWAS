package simulations;


import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.stat.StatUtils;
import org.jfree.chart.plot.ThermometerPlot;

import mixedmodel.LocalKinshipAnalyzer;
import mixedmodel.MultiPhenotype;
import mixedmodel.VariantsByte;
import mixedmodel.VariantsDouble;
import myMathLib.Test;


import cern.jet.random.Normal;
import cern.jet.random.engine.DRand;



/**
 * Simulate the phenotype by given real genotype and parameters of the focal model
 * @author quan.long
 *
 */

public class Simulator {
	
	public static String additive_exp_pow="additive_exp_pow";
	public static String additiveInfinitesimal="infinitesimal";
	public static String additiveBidirection="additiveBidirection";
	public static String liabilityThreshold="liabilityThreshold";
//	public static String heterogeneityBackground="heterogeneityBackground"; //TODO
	public static String heterogeneityEpistasis="heterogeneityEpistasis";
	public static String heterogeneityDorminant="heterogeneityDorminant";
	public static String heterogeneityBidirection="heterogeneityBidirection";
//	public static String single_main="single_main";
	// Two-way Epistatis: num_of_causal = 2:
	public static String canalization="canalization";
	public static String compensation="compensation";
	public static String signepistasis="signepistasis";
	
	public static String with100="with100";
		
	VariantsDouble genotype;
	//VariantsByte gwas_data_byte;
	
	public double heritability;
	public double regional_percentage; 
	public int[][] qualified_variants;
	public int[][] causal_variants;  // record the causal variants for further analysis
	int[][] indexes_global_background;
	int num_variants_for_global;
	double[] global_effects_values;
	
	//double[] phenotype_values;
	int chr_causal_region; 	//the chr that contains causal variants. 
				//If chr==-1, then the causal variants are distributed in multiple chrs.
	

	/*
	 * Generate variants in whole genome
	 */
	public Simulator(String input_genotype, double heritability, 
			double min_maf, double max_maf) {
		
		this.genotype=new VariantsDouble(input_genotype);
		this.heritability=heritability;
		this.qualified_variants=this.genotype.find_vars_in_maf_range(min_maf, max_maf);
	}
	
	/*
	 * Generate variants in local region(s)
	 */
	public Simulator(String input_genotype, double heritability, 
			double min_maf, double max_maf, int chr, int start_location, int end_location) {		
		this.genotype=new VariantsDouble(input_genotype);
		this.heritability=heritability;
		this.chr_causal_region=chr;
		if(chr==-1){
			this.qualified_variants=this.genotype.find_vars_in_maf_range(min_maf, max_maf);
		}else{
			this.qualified_variants=new int[this.genotype.num_chrs][];
			this.causal_variants=new int[this.genotype.num_chrs][];
			for(int chr_index=0;chr_index<this.genotype.num_chrs;chr_index++){
				this.qualified_variants[chr_index]=new int[0];
				this.causal_variants[chr_index]=new int[0];
				
			}
			this.qualified_variants[chr]=this.genotype.find_vars_in_maf_range_in_region(chr, 
				start_location, end_location, min_maf, max_maf);
		}
	}
	
	/*
	 * Generate variants in local region(s) with global background
	 * 
	 * if (chr==-1)  look at all chromosomes; 
	 * else  check only particular region.
	 *  TODO !!!!
	 */
	public Simulator(String input_genotype, double heritability, double regional_percentage,
			double min_maf, double max_maf, int region_chr, int start_location, int end_location,
			int num_variants_for_global, double global_min_maf, double global_max_maf) {
		this.genotype=new VariantsDouble(input_genotype);
		this.heritability=heritability;
		this.regional_percentage=regional_percentage;
		this.num_variants_for_global=num_variants_for_global;		
		// generate global variants		
		this.indexes_global_background=this.genotype.find_fixed_num_of_vars_in_maf_range(global_min_maf, 
				global_max_maf, num_variants_for_global);
		this.global_effects_values=new double[this.genotype.sample_size];
		for(int the_chr=0;the_chr<this.genotype.num_chrs;the_chr++){
			for(int i=0;i<this.indexes_global_background[the_chr].length;i++){
				double[] the_genotypes=this.genotype.load_one_variant_by_index(the_chr, this.indexes_global_background[the_chr][i]);
				for(int k=0;k<this.genotype.sample_size;k++){
					this.global_effects_values[k]+=the_genotypes[k];
				}
			}
		}		
		// generate local variants
		if(region_chr==-1){
			this.qualified_variants=this.genotype.find_vars_in_maf_range(min_maf, max_maf);
		}else{
			this.qualified_variants=new int[this.genotype.num_chrs][];
			for(int chr_index=0;chr_index<this.genotype.num_chrs;chr_index++){
				if(chr_index!=region_chr)this.qualified_variants[chr_index]=new int[0];
			}
			this.qualified_variants[region_chr]=this.genotype.find_vars_in_maf_range_in_region(region_chr, 
				start_location, end_location, min_maf, max_maf);
		}
	}
	
	/*
	 * given a list of filtered qualified_vars (e.g., MAF, or S/NS), randomly select some variants for simulation
	 * multiple chr version
	 */
	public void ramdom_selected(int num_var_needed){		
		int total_candidates=0;
		int chr_num=this.qualified_variants.length;
		for(int chr=0;chr<chr_num;chr++){
			total_candidates+=this.qualified_variants[chr].length;
		}
		if(total_candidates<=num_var_needed){
			System.out.println("Total number of canditate ("+total_candidates+") is not more than needed ("+num_var_needed+"). " +
					"Use all of them");
			this.causal_variants=this.qualified_variants.clone();
			return ;
		}
		int[] locs=new int[total_candidates]; 
		int[] chrs=new int[total_candidates]; 
		int index=0;
		for(int chr=0;chr<chr_num;chr++){
			for(int k=0;k<this.qualified_variants[chr].length;k++){
				locs[index]=this.qualified_variants[chr][k];
				chrs[index]=chr;
				index++;
			}
		}
		ArrayList<Integer>[] useful=new ArrayList[chr_num];
		for(int chr=0;chr<chr_num;chr++){
			useful[chr]= new ArrayList<>();
		}
		for(int i=0;i<num_var_needed;i++){
			int the_selected_index=(int)(Test.randomNumber()*total_candidates);
			if(the_selected_index==total_candidates)the_selected_index--;
			while(useful[chrs[the_selected_index]].contains(locs[the_selected_index])){
				the_selected_index=(int)(Test.randomNumber()*total_candidates);
				if(the_selected_index==total_candidates)the_selected_index--;
			}
			useful[chrs[the_selected_index]].add(locs[the_selected_index]);			
		}
		this.causal_variants=new int[chr_num][];
		for(int chr=0;chr<chr_num;chr++){
			this.causal_variants[chr]=myFileFunctions.FileFunc.arraylist2arrayInteger(useful[chr]);
		}
	}
	
	/*
	 * given a list of filtered qualified_vars (e.g., MAF, or S/NS), randomly select some variants for simulation
	 * single chr version
	 */
	public void random_selected(int chr, int num_var_needed){		
		if(this.qualified_variants[chr].length<=num_var_needed){
			System.out.println("Total number of canditate is smaller than needed.");
			this.causal_variants[chr]= this.qualified_variants[chr].clone();
		}
		int total_candidates=this.qualified_variants[chr].length;
		ArrayList<Integer> useful= new ArrayList<>();
		for(int i=0;i<num_var_needed;i++){
			int the_selected_index=(int)(Test.randomNumber()*total_candidates);
			if(the_selected_index==total_candidates)the_selected_index--;
			while(useful.contains(this.qualified_variants[chr][the_selected_index])){
				the_selected_index=(int)(Test.randomNumber()*total_candidates);
				if(the_selected_index==total_candidates)the_selected_index--;
			}
			useful.add(this.qualified_variants[chr][the_selected_index]);			
		}
		this.causal_variants[chr]= myFileFunctions.FileFunc.arraylist2arrayInteger(useful);
	}
	
	public void output_causal(String output_causal_file){
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(output_causal_file));
			bw.write("Chr(FromOne)\tIndexInData(FromZero)\tLocation\n");
			for(int chr=0;chr<this.causal_variants.length;chr++){
				if(this.causal_variants[chr].length!=0){
					for(int k=0;k<this.causal_variants[chr].length;k++){
						int the_selected_index=this.causal_variants[chr][k];
						int the_loc=this.genotype.locations[chr][the_selected_index];
						bw.write((chr+1)+"\t"+the_selected_index+"\t"+the_loc+"\n");
					}
				}
			}bw.close();
		}catch(Exception e){e.printStackTrace();}
	}
	
	/*
	 * given a list of filtered qualified_vars (e.g., MAF, or S/NS), randomly select some variants for simulation
	 * single chr version
	 * Evenly distributed so that LD between markers will be minimal. 
	 * NOT RECOMMENDED TO BE USED WITH VERY HIGH NUMBER OF VARIANTS 
	 * 
	 *  TODO
	 */
	public static int[] evenly_selected(int[] qualified_vars, int num_var_needed){		
		if(qualified_vars.length<=num_var_needed){
			System.out.println("Total number of canditate is smaller than needed.");
			return qualified_vars;
		}
		int total_candidates=qualified_vars.length;
		int[] useful=new int[num_var_needed];
		for(int i=0;i<num_var_needed;i++){
			int the_selected_index=(int)(i*(total_candidates/(num_var_needed-1)));
			if(the_selected_index==total_candidates)the_selected_index--;
			useful[i]=qualified_vars[the_selected_index];			
		}
		return useful;
	}
	
	/*
	 * Additive model: number of SNPs, and its effects
	 * exp_ratio \in (0,1) when additive_exp_pow is chosen
	 * or exp_ratio==-1 when infinitesimal is chosen
	 * P= \Sigma a^i
	 * 
	 * When num_of_causal==1, it is just single main effect. 
	 */
	
	public void additive(String model, int num_of_causal, double exp_ratio, int num_replicates, 
			boolean global_background, String output_pheno_file, String output_causal_file){			
		double[][] phenotypes=new double[num_replicates][this.genotype.sample_size];
		for(int rep=0;rep<num_replicates;rep++){
			double[] phenotype_values=new double[this.genotype.sample_size];
			double[] effects=new double[num_of_causal];
			if(model.equals(Simulator.additiveInfinitesimal) && Double.compare(exp_ratio,-1)==0){
//				Exponential exp=new Exponential(1, (new DRand()));
				NormalDistribution norm=new NormalDistribution();
				for(int i=0;i<effects.length;i++){
					effects[i]=norm.sample();
				}			
			}else if(model.equals(Simulator.additive_exp_pow) && exp_ratio>0 && exp_ratio<1){
				effects[0]=0.5; // since we are going to do normalization to the error_weight, the absolute value of the first one doesn't matter
				for(int i=1;i<effects.length;i++){
					effects[i]=effects[i-1]*exp_ratio;
				}
			}else if(model.equals(Simulator.additiveBidirection) && Double.compare(exp_ratio,-1)==0){
				for(int i=0;i<effects.length;i++){
					NormalDistribution norm=new NormalDistribution();
					if(i%2==0)effects[i]=norm.sample();
					else effects[i]=(-1)*norm.sample();
				}
			}else{
				System.out.println("User input model \""+model+"\" when exp_ratio="+exp_ratio+" is not supported.");
				return;
			}
			ramdom_selected(num_of_causal);
			int effects_index=0;
			for(int chr=0;chr<this.causal_variants.length;chr++){
				for(int loc_index=0;loc_index<this.causal_variants[chr].length;loc_index++){
					double[] genotypes=this.genotype.load_one_variant_by_index(chr, 
							this.causal_variants[chr][loc_index]);
					for(int indi=0;indi<this.genotype.sample_size;indi++){					
						phenotype_values[indi]+=(effects[effects_index]*genotypes[indi]);											
					}
					effects_index++;
				}
			}	
			if(global_background)phenotype_values=add_background(phenotype_values,this.global_effects_values, this.regional_percentage);
			normalization(phenotype_values, this.heritability);	
			phenotypes[rep]=phenotype_values.clone();
			this.output_causal(output_causal_file+"."+rep+".txt");
		}
		this.write_simulated_pheno_file(phenotypes, output_pheno_file);	
		System.out.println("Finsihed simulating phenotype with additive model and writing to "+output_pheno_file);
	}
	
	/*
	 *  concentration_of_causal instead of a fixed number in the regions.
	 */
	public void additive2(String model, double concentration_of_causal, double exp_ratio, int num_replicates, boolean global_background,
			String output_pheno_file, String output_causal_file){			
		int total_num_candidates=0;
		int chr_num=this.qualified_variants.length;
		for(int chr=0;chr<chr_num;chr++){
			total_num_candidates+=this.qualified_variants[chr].length;
		}
		int num_of_causal=(int)(total_num_candidates*concentration_of_causal);
		this.additive(model, num_of_causal, exp_ratio, num_replicates, global_background, output_pheno_file, output_causal_file);
		System.out.println("Causal: "+total_num_candidates+"/"+num_of_causal);
		
	}
	
	public void liabilityThreshold(String model, int num_of_causal, int num_replicates, 
			boolean global_background, String output_pheno_file, String output_causal_file){			
		double[][] phenotypes=new double[num_replicates][this.genotype.sample_size];
		for(int rep=0;rep<num_replicates;rep++){
			double[] phenotype_values=new double[this.genotype.sample_size];
			double[] effects=new double[num_of_causal];
			if(model.equals(Simulator.liabilityThreshold)){
				for(int i=0;i<effects.length;i++){
					effects[i]=(i+1)%2-0.5;
				}			
			}else{
				System.out.println("User input model \""+model+"\" is not supported.");
				return;
			}
			ramdom_selected(num_of_causal);
			int effects_index=0;
			for(int chr=0;chr<this.causal_variants.length;chr++){
				for(int loc_index=0;loc_index<this.causal_variants[chr].length;loc_index++){
					double[] genotypes=this.genotype.load_one_variant_by_index(chr, 
							this.causal_variants[chr][loc_index]);
					for(int indi=0;indi<this.genotype.sample_size;indi++){					
						phenotype_values[indi]+=(effects[effects_index]*genotypes[indi]);											
					}
					effects_index++;
				}
			}
			if(global_background)phenotype_values=add_background(phenotype_values,this.global_effects_values, this.regional_percentage);
			
			// threshold is 0: 
			for(int indi=0;indi<this.genotype.sample_size;indi++){					
				if(phenotype_values[indi]>0)phenotype_values[indi]=1;
				else phenotype_values[indi]=0;
			}
			normalization(phenotype_values, this.heritability);	
			phenotypes[rep]=phenotype_values.clone();
			this.output_causal(output_causal_file+"."+rep+".txt");
		}
		this.write_simulated_pheno_file(phenotypes, output_pheno_file);	
		System.out.println("Finsihed simulating phenotype with liability model and writing to "+output_pheno_file);
	}
	
	public void liabilityThreshold2(String model, double concentration_of_causal, int num_replicates, 
			boolean global_background, String output_pheno_file, String output_causal_file){			
		int total_num_candidates=0;
		int chr_num=this.qualified_variants.length;
		for(int chr=0;chr<chr_num;chr++){
			total_num_candidates+=this.qualified_variants[chr].length;
		}
		int num_of_causal=(int)(total_num_candidates*concentration_of_causal);
		this.liabilityThreshold(model, num_of_causal, num_replicates, global_background, output_pheno_file, output_causal_file);		
	}
	
	
	/*
	 * Heterogeneity model: 
	 * 
	 * model="heterogeneityDorminant" any single mutation specified are enough to affect phenotype. 
	 * model="heterogeneityRecessive" then one need to carry a number of them.
	 * 
	 */
	public void heterogeneity(String model, int total_num_of_causal, int num_replicates, 
			int num_allele_needed_to_have_efffet, boolean global_background, String output_pheno_file,
			String output_causal_file){			
		double[][] phenotypes=new double[num_replicates][this.genotype.sample_size];
		for(int rep=0;rep<num_replicates;rep++){
			double[] phenotype_values=new double[this.genotype.sample_size];
			int[] num_carrying=new int[this.genotype.sample_size];
			if(model.equals(Simulator.heterogeneityDorminant) && num_allele_needed_to_have_efffet!=1){
				System.out.println("Under Drominant model, single mutation will have effect.\n" +
						"Please specify num_needed_to_have_effect as 1");
				return;
			}
			ramdom_selected(total_num_of_causal);
			for(int chr=0;chr<this.causal_variants.length;chr++){
				for(int loc_index=0;loc_index<this.causal_variants[chr].length;loc_index++){
					double[] genotypes=this.genotype.load_one_variant_by_index(chr, 
							this.causal_variants[chr][loc_index]);
					for(int indi=0;indi<this.genotype.sample_size;indi++){					
						num_carrying[indi]+=genotypes[indi];											
					}
				}
			}		
			for(int i=0;i<this.genotype.sample_size;i++){
				if(num_carrying[i]>=num_allele_needed_to_have_efffet)phenotype_values[i]=1;
			}
			if(global_background)phenotype_values=add_background(phenotype_values,this.global_effects_values, this.regional_percentage);
			normalization(phenotype_values, this.heritability);	
			phenotypes[rep]=phenotype_values.clone();
			this.output_causal(output_causal_file+"."+rep+".txt");
		}
		this.write_simulated_pheno_file(phenotypes, output_pheno_file);	
		System.out.println("Finsihed simulating phenotype with heterogeneity model and writing to "+output_pheno_file);
	}
	
	/*
	 * instead of specifying total_num_of_causal, we specify concentration_of_causal.
	 */
	public void heterogeneity2(String model, double concentration_of_causal, int num_replicates, int num_allele_needed_to_have_efffet, 
			boolean global_background, String output_pheno_file, String output_causal_file){			
		
		int total_num_candidates=0;
		int chr_num=this.qualified_variants.length;
		for(int chr=0;chr<chr_num;chr++){
			total_num_candidates+=this.qualified_variants[chr].length;
		}
		int total_num_of_causal=(int)(total_num_candidates*concentration_of_causal);
		this.heterogeneity(model, total_num_of_causal, num_replicates, num_allele_needed_to_have_efffet,
				global_background, output_pheno_file, output_causal_file);
	}
	
	/*
	 * 
	 * Epistasis I: "canalization"; II: compensation="compensation"; III: signepistasis;
	 * only 4 causal variants will be selected, and the one with highest MAF will be the controller. 
	 * 
	 */
	public void local_interactions(String model, int num_replicates, String output_pheno_file, String output_causal_file){
		int num_causal=4;
		if(this.chr_causal_region==-1){
			System.out.println("Function local_interactions doesn't work for multuple chrs.");
			return;
		}
		NormalDistribution normal_01=new NormalDistribution();
		double[][] phenotypes=new double[num_replicates][this.genotype.sample_size];
		for(int rep=0;rep<num_replicates;rep++){
			double[] phenotype_values=new double[this.genotype.sample_size];
			this.random_selected(this.chr_causal_region, 4);		
			int controler_index=0;
			double controler_MAF=this.genotype.mafc[chr_causal_region][causal_variants[chr_causal_region][0]];
			for(int i=1;i<num_causal;i++){
				double the_other_MAF=this.genotype.mafc[chr_causal_region][causal_variants[chr_causal_region][i]];
				if(controler_MAF<the_other_MAF){
					controler_MAF=the_other_MAF;
					controler_index=i;
				}
			}
			double[] genotypes_controller=this.genotype.load_one_variant_by_index(chr_causal_region, 
							this.causal_variants[chr_causal_region][controler_index]);
			double[][] genotypes_effector=new double[num_causal-1][];
			int effector_index=0;
			for(int i=0;i<num_causal;i++){
				if(i!=controler_index){
					genotypes_effector[effector_index]=this.genotype.load_one_variant_by_index(chr_causal_region, 
							this.causal_variants[chr_causal_region][i]);
					effector_index++;
				}
			}
			if(model.equals(canalization)){
				for(int indi=0;indi<this.genotype.sample_size;indi++){	
					if(genotypes_controller[indi]==0){
						phenotype_values[indi]=0;
					}else{ //genotypes_controller[indi] !=0
						for(int i=0;i<num_causal-1;i++){ // any effector will work
							if(genotypes_effector[i][indi]>0){
								phenotype_values[indi]=1;break;
							}
						}
					}																
				}
			}else if(model.equals(signepistasis)){
				for(int indi=0;indi<this.genotype.sample_size;indi++){	
					double sign=(genotypes_controller[indi]==0)?1:-1;
					for(int i=0;i<num_causal-1;i++){ // any effector will work
						if(genotypes_effector[i][indi]>0.01){
							phenotype_values[indi]=sign;break;
						}
					}
				}					
			}else if(model.equals(compensation)){
				for(int indi=0;indi<this.genotype.sample_size;indi++){	
					int count=(genotypes_controller[indi]==0)?0:1;
					for(int i=0;i<num_causal-1;i++){ 
						if(genotypes_effector[i][indi]>0){
							count++;
						}
					}phenotype_values[indi]=count%2;
				}	
			}
			normalization(phenotype_values, this.heritability);	
			phenotypes[rep]=phenotype_values.clone();
			this.output_causal(output_causal_file+"."+rep+".txt");
		}
		this.write_simulated_pheno_file(phenotypes, output_pheno_file);	
		System.out.println("Finsihed simulating phenotype with epistasis ("+model+") model and writing to "+output_pheno_file);
	}
	
	public void write_simulated_pheno_file(double[][] phenos, String file){
		try{
			BufferedWriter bw= new BufferedWriter(new FileWriter(file));
			bw.write("ID");
			for(int i=0;i<phenos.length;i++){
				bw.write("\tSim_Pheno_"+i);
			}bw.write("\n");
			for(int indi=0;indi<this.genotype.sample_size;indi++){
				bw.write(this.genotype.sample_ids[indi]);
				for(int i=0;i<phenos.length;i++){
					bw.write("\t"+phenos[i][indi]);
				}bw.write("\n");
			}bw.close();
		}catch(Exception e){e.printStackTrace();}
	}
	
	public static double[] add_background(double[] regional, double[] global, double regional_percentage){
		double[] sum1=regional.clone();
		double[] sum2=global.clone();
		myMathLib.Normalization.normalize_variance(sum1);
		myMathLib.Normalization.normalize_variance(sum2);
		double c1=Math.sqrt(regional_percentage);
		double c2=Math.sqrt(1-regional_percentage);
		for(int k=0;k<sum1.length;k++)sum1[k]=sum1[k]*c1+sum2[k]*c2;
		return sum1;
	}
	
	/*
	 * adding the variance of error term
	 * change will be made to the double[] original 
	 */
	public static void normalization(double[] original, double heritability){
		NormalDistribution norm=new NormalDistribution();
		double sigma_g2=StatUtils.populationVariance(original);
		double sigma_e=Math.sqrt(sigma_g2*(1-heritability)/heritability);
		for(int indi=0;indi<original.length;indi++){
			original[indi]+=(sigma_e*norm.sample());
		}
	}
	
	public static void generate_random_phenotype(String output_file, String[] sample_ids){
		try{
			BufferedWriter bw =new BufferedWriter(new FileWriter(output_file));
			bw.write("ids\trandom\n");
			for(int i=0;i<sample_ids.length;i++){
				bw.write(sample_ids[i]+"\t"+Test.randomNumber()+"\n");
			}bw.close();
		}catch(Exception e){e.printStackTrace();}
	}
	
	// Try programs while developing:
	public static void main(String[] args) {
		String input_genotype="/Volumes/Projects/DATA/arabidopsis_genomes/swe180_ecker171_removebads.csv.num2.hdf5";
		double heritability=0.7, min_maf=0.1,max_maf=0.2, global_min_maf=0.2, global_max_maf=0.5, regional_percentage=0.8;
		int chr=3, start_location=200000, end_location=400000, num_variants_for_global=1000;
		Simulator sim=new Simulator(input_genotype, heritability, regional_percentage, min_maf, max_maf, chr, start_location, 
				end_location, num_variants_for_global, global_min_maf, global_max_maf);
	}
		
}
