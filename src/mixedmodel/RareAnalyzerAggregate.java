package mixedmodel;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashMap;

import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;

public class RareAnalyzerAggregate {
	VariantsDouble variants;
	Phenotype phenotype;
	double[][] kinship;
	
	Regions regions;
    double rare_threshold;
	
	
	/*
	 * the region file should be in the form: name\tchr\tstart\tend
	 */
	public RareAnalyzerAggregate(String genotype_file, String regions_file, Phenotype phenotype, String kinship_file,
			 double rare_threshold){
		System.out.println(phenotype.phe_id);
		this.variants=new VariantsDouble(genotype_file);
		this.phenotype=phenotype;
		this.kinship= new KinshipMatrix(kinship_file).kinship_matrix.getData();//this.variants.read_in_kinship(kinship_file);
		this.rare_threshold=rare_threshold;	
		this.regions=new Regions(regions_file, this.variants);				
	}
	
	
	public RareAnalyzerAggregate(String genotype_file, int[][][] regions_array, String[][] region_names, Phenotype phenotype, String kinship_file,
			 double rare_threshold){
		System.out.println(phenotype.phe_id);
		this.variants=new VariantsDouble(genotype_file);
		this.phenotype=phenotype;
		this.kinship= new KinshipMatrix(kinship_file).kinship_matrix.getData();//this.variants.read_in_kinship(kinship_file);
		this.rare_threshold=rare_threshold;	
		this.regions=new Regions(regions_array, region_names);				
	}
	/*
	 * Run naive collapse, with mixed model
	 */
	public void rare_association(double rare_threshold, String output_folder, boolean plot){
		try{
			File out_folder_dir=new File(output_folder);
			if(!out_folder_dir.exists())out_folder_dir.mkdir();
			
			String result_file=out_folder_dir+"/Aggregate."+this.phenotype.phe_id+".csv";
			BufferedWriter bw_mix=new BufferedWriter(new FileWriter(result_file));			
			String header0="#chr,location,pvalue,AdjustedR2,coefficient,Sd_Err,MAF_count,Region_id\n";			
			bw_mix.write("\n"+header0);
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
				for(int region_index=0;region_index<this.regions.region_coordinates[chr].length;region_index++){
					int start_location=this.regions.region_coordinates[chr][region_index][0];
					int end_location=this.regions.region_coordinates[chr][region_index][1];
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
							bw_mix.write((chr+1)+","+this.regions.region_coordinates[chr][region_index][0]+","+out_res+","+
									this.regions.region_names[chr][region_index]+"\n");
						}
					}							
				}System.out.println("Aggregate finished Chr"+(chr+1));
			}			
			bw_mix.close();
			if(plot)EMMAX.make_plot_one_phen(result_file, 0.05);
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
	
}
