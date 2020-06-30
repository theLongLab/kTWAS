package mixedmodel;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashMap;

import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;


public class RareAnalyzerSynthetic {
	VariantsDouble variants;
	AssociationResults potential_synthetic_association;
	Phenotype phenotype;
	double[][] kinship;
	
	int[][][] region_coordinates; //chr x num_region x 2(start, end)
	String[][] region_names; //chr x num_region 
	int[][][] syn_of_around_region; // chr x num_region x num_syn_of_the_region
	
	int distance2region;
    double rare_threshold;
	double ld_threshold;
	
	public RareAnalyzerSynthetic(String synthetic_results_file, double gwas_pvalue_threshold, String genotype_file,
			int grouping_region_win_size, int distance2region, Phenotype phenotype, String kinship_file){
		this.variants=new VariantsDouble(genotype_file);
		this.potential_synthetic_association=new AssociationResults(synthetic_results_file, gwas_pvalue_threshold);
		this.distance2region=distance2region;
		this.phenotype=phenotype;
		this.kinship= new KinshipMatrix(kinship_file).kinship_matrix.getData();//this.variants.read_in_kinship(kinship_file);
		this.region_coordinates=null;
		
	}
	
	/*
	 * the region file should be in the form: name\tchr\tstart\tend
	 */
	public RareAnalyzerSynthetic(String synthetic_results_file, double synthetic_pvalue_threshold, double synthetic_maf_threhold, 
			String genotype_file, String regions_file, int distance2region, Phenotype phenotype, String kinship_file,
			double ld_threshold, double rare_threshold){
		System.out.println(phenotype.phe_id);
		this.variants=new VariantsDouble(genotype_file);
		this.potential_synthetic_association=new AssociationResults(synthetic_results_file, synthetic_pvalue_threshold, synthetic_maf_threhold);
		this.distance2region=distance2region;
		this.phenotype=phenotype;
		this.kinship= new KinshipMatrix(kinship_file).kinship_matrix.getData();//this.variants.read_in_kinship(kinship_file);
		this.rare_threshold=rare_threshold;
		this.ld_threshold=ld_threshold;
		try{
			BufferedReader br_region=new BufferedReader(new FileReader(regions_file));
			String line=br_region.readLine();
			ArrayList<Integer> chrs=new ArrayList<Integer>();
			HashMap<Integer, Integer> num_regions_chr=new HashMap<Integer, Integer>();
			while(line!=null){
				String[] temp=line.split("\t");
				int the_chr=Integer.parseInt(temp[1]);
				if(!chrs.contains(the_chr))chrs.add(the_chr);
				myFileFunctions.FileFunc.add2_counts_hashmap(num_regions_chr, the_chr, 1);
				line=br_region.readLine();
			}
			if(chrs.size()!=this.variants.num_chrs){
				System.out.println("Number of Chrs are not consistent between genotype file and regions file.");
				return;
			}
			this.region_coordinates=new int[this.variants.num_chrs][][];
			this.region_names=new String[this.variants.num_chrs][];
			for(int k=0;k<this.variants.num_chrs;k++){
				this.region_names[k]=new String[num_regions_chr.get(chrs.get(k))];
				this.region_coordinates[k]=new int[num_regions_chr.get(chrs.get(k))][2];
			}
			br_region=new BufferedReader(new FileReader(regions_file));
			line=br_region.readLine();
			int[] indexes=new int[variants.num_chrs];
			while(line!=null){
				String[] temp=line.split("\t");
				int the_chr=Integer.parseInt(temp[1]);
				this.region_names[the_chr-1][indexes[the_chr-1]]=temp[0];
				this.region_coordinates[the_chr-1][indexes[the_chr-1]][0]=Integer.parseInt(temp[2]);
				this.region_coordinates[the_chr-1][indexes[the_chr-1]][1]=Integer.parseInt(temp[3]);
				indexes[the_chr-1]++;
				line=br_region.readLine();
			}
			System.out.println("Finished loading data for rare variants analysis.");
		}catch(Exception e){e.printStackTrace();}		
	}
	
	
	
	/*
	 * Run 4 types of analysis:
	 * (1) naive collapse, no mixed model
	 * (2) ld collapse, no mixed model
	 * (1) naive collapse, with mixed model
	 * (1) naive collapse, with mixed model
	 */
	public void rare_association_4in1(double rare_threshold, String output_folder){
		try{
			File out_folder_dir=new File(output_folder);
			if(!out_folder_dir.exists())out_folder_dir.mkdir();
			BufferedWriter bw_mix_syn=new BufferedWriter(new FileWriter(out_folder_dir+"/"+this.phenotype.phe_id+".mix.syn.csv"));
			BufferedWriter bw_mix=new BufferedWriter(new FileWriter(out_folder_dir+"/"+this.phenotype.phe_id+".mix.nosyn.csv"));
			BufferedWriter bw_syn=new BufferedWriter(new FileWriter(out_folder_dir+"/"+this.phenotype.phe_id+".nomix.syn.csv"));
			BufferedWriter bw_nothing=new BufferedWriter(new FileWriter(out_folder_dir+"/"+this.phenotype.phe_id+".nomix.nosyn.csv"));
			String header0="#chr,location,pvalue,AdjustedR2,coefficient,Sd_Err,MAF_count,Region_id\n";
			String header1="#chr,location,pvalue,AdjustedR2,coefficient,Sd_Err,MAF_count,Region_id,syn_loc\n";
			bw_mix_syn.write("\n"+header1);bw_mix.write("\n"+header0);bw_syn.write("\n"+header1);bw_nothing.write("\n"+header0);
			EMMA emma=new EMMA(this.phenotype, this.variants, this.kinship);
			emma.REMLE_null();
			double[][] decomposed_array = emma.reml_decompositioned_array();				
			emma.phenotype.generate_new_Y_by_multiplying(decomposed_array); 
			double[] intsept=new double[emma.sample_size];
			for(int i=0;i<emma.sample_size;i++){
				for(int j=0;j<emma.sample_size;j++){
					intsept[i]+=decomposed_array[i][j];
				}
			}
			int mafc_threshold=(int)(rare_threshold*emma.sample_size*2);
			for(int chr=0;chr<variants.num_chrs;chr++){
				for(int region_index=0;region_index<this.region_coordinates[chr].length;region_index++){
					int start_location=this.region_coordinates[chr][region_index][0];
					int end_location=this.region_coordinates[chr][region_index][1];
					double[][] data=EMMAX.load_data_according_phenotype(chr, start_location, end_location, emma);
					if(data==null)continue;
					double[] collapsed=collapse_naive(mafc_threshold, data);
					int mafc=EMMAX.mafc(collapsed);
					// apply EMAAX 
					double[][] coverted_X=EMMAX.apply_array(intsept, collapsed, decomposed_array);
					if(coverted_X==null)continue;
					double[][] result_emmax=EMMAX.reg2results_emmax(emma.phenotype.new_Y_4emmax, coverted_X);
					if(result_emmax!=null){
						result_emmax[3]=new double[1];
						result_emmax[3][0]=mafc;
						String out_res=EMMAX.results2string(result_emmax, 1); // All results printed 
						if(out_res!=null){
							bw_mix.write((chr+1)+","+this.region_coordinates[chr][region_index][0]+","+out_res+","+this.region_names[chr][region_index]+"\n");
						}
					}						
					// apply LM					
					double[][] Xs=new double[emma.sample_size][1];
					for(int sample_index=0; sample_index<emma.sample_size; sample_index++){
						Xs[sample_index][0] = collapsed[sample_index];
					}
					boolean the_same=true;
					for(int sample_index=1; sample_index<emma.sample_size; sample_index++){
						if(Double.compare(Xs[sample_index][0], Xs[0][0])!=0)the_same=false;
					}				
					if(the_same)continue;
					//OLSMultipleLinearRegression reg1=new OLSMultipleLinearRegression();									
					//reg1.newSampleData(emma.phenotype.values, Xs);
					double[][] result_lm=EMMAX.reg2results_lm(emma.phenotype.values, Xs);
					if(result_lm!=null){
						result_lm[3]=new double[1];
						result_lm[3][0]=mafc;
						String out_res=EMMAX.results2string(result_lm, 1);
						if(out_res!=null){
							bw_nothing.write((chr+1)+","+this.region_coordinates[chr][region_index][0]+","+out_res+","+this.region_names[chr][region_index]+"\n");
						}
					}			
					// leverage synthetic associations
//					for(int syn_index=0;syn_index<this.syn_of_around_region[chr][region_index].length;syn_index++){
//						double[] collapsed_syn=collapse_ld(mafc_threshold, data, this.ld_threshold, chr, 
//								this.syn_of_around_region[chr][region_index][syn_index]);
//						// did EMAAX and LM
//					}					
				}
			}
			
			bw_mix.close();bw_mix_syn.close();bw_nothing.close();bw_syn.close();
			EMMAX.make_plot_one_phen(out_folder_dir+"/"+this.phenotype.phe_id+".mix.nosyn.csv", 0.05);
			EMMAX.make_plot_one_phen(out_folder_dir+"/"+this.phenotype.phe_id+".nomix.nosyn.csv", 0.05);
		}catch(Exception e){e.printStackTrace();}			
	}
	
	/*
	 * data[][] = num_snp x sample_size
	 */
	public static double[] collapse_naive(int mafc_threshold, double[][] data){
		int sample_size=data[0].length;
		double[] collapsed=new double[sample_size];
		for(int snp=0;snp<data.length;snp++){
			int[] allele_and_mafc=EMMAX.mafc_with_allele(data[snp]);
			if(allele_and_mafc!=null){
				int allele=allele_and_mafc[0], count=allele_and_mafc[1];
				if(count<mafc_threshold){ // this is a rare variant
					for(int k=0;k<sample_size;k++){
						if((int)(data[snp][k]+0.5)==allele)collapsed[k]=2;
					}					
				}
			}
		}
		return collapsed;
	}
	
	public static double[] collapse_ld(int mafc_threshold, double[][] data, double ld_threshold, int syn_chr, int syn_loc){
		return null;
	}
	
	/*
	 *  initiate this.syn_of_around_region based on 
	 *  	this.region_coordinates 
	 *      this.potential_synthetic_association
	 *      this.distance2region
	 */	
	public void find_syn_around_regions(){
	
	}
	
	/*
	 * find indexes of regions that are located surrounding the potential synthetic association.
	 */
	public int[] regions(int chr, int start_loc, int end_loc){
		return null;
	}
	
	public static void main(String[] args){
		long start=System.currentTimeMillis();
		String working_folder="/Users/quan.long/Documents/Projects/GeneMappingMethodology/rare/";
		String regions_file=working_folder+"gene_regions.txt";
		String phenotype_all=working_folder+"gen_arch_traits_70_plus.useful.tsv";
		String genotype_file=working_folder+"swe180_ecker171_removebads.csv.num.mafc.hdf5";
		String kinship_file=working_folder+"kinship.rescaled.ibs";
		String synthetic_results_file_30=working_folder+"EMMAX.30_2053_germ_284_t7.top";
		String outputput_folder=working_folder+"test_rare_germ_284";
		
		double synthetic_pvalue_threshold=0.00001;
		double synthetic_maf_threhold=0.05;
		int distance2region=100000;
		double ld_threshold=0.8;
		double rare_threshold=0.02;
		
		MultiPhenotype all_phe=new MultiPhenotype(phenotype_all);
		
		RareAnalyzerSynthetic analyzer = new RareAnalyzerSynthetic(synthetic_results_file_30, synthetic_pvalue_threshold, synthetic_maf_threhold, genotype_file, 
				regions_file, distance2region, all_phe.phenotypes[30], kinship_file, ld_threshold, rare_threshold); 
		
		analyzer.rare_association_4in1(rare_threshold, outputput_folder);
		System.out.println("Timne used: "+(System.currentTimeMillis()-start)/1000+" sec");
	}
	
}
