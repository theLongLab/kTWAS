package simulations;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Random;

import mixedmodel.VariantsDouble;
import myMathLib.Normalization;
import myMathLib.Test;

public class SimHeterogeneity {
	
	/**
	 * Simulate phenotype based on various patterns of heterogeneity
	 * (1) MAF (how rare are the variants)
	 * (2) number of regions (how heterogenous) 
	 * (3) number of variants 
	 * (4) effect sizes
	 * 		(4.1) For quantitative traits:
	 * 			  We have an additive model with a genetic variance component and residue variance component;
	 * 			  Multiple variants from the same region only contribute at most twice (reflecting compound heterozygotes)
	 * 				Phenotype = \Sigma_regions ( \Sigma_in-region_var; max 2)    
	 * 		(4.2) For categorical, i.e., case/control data:
	 * 			  Simulate a quantitative trait first, and set to case if the score is higher than a liability threshold. 
	 * 
	 * In order to control how many variants contributing to the phenotype, one can tune the following parameters:
	 * (1) MAF: the lower, the smaller number of variants can contribute. 
	 * (2) In the Case/Control study, lower down the liability threshold so that one variant is sufficient to cause the disease,
	 * 		or lower down the total number of variants. 
	 * 		(Note that, in general the quantitative traits are not caused by very few genetic variants)  
	 */
	
	static double pheno_unit=1.0;
	
	int num_candidate_regions;
	int num_causal_regions;
	int[][] candidate_causal_indexes;
	int[][] candidate_regions_locations; //chr_index, start_locus_location, end_locus_location
	HashMap<Integer, Integer> region2chr;
	double genetic_component; // a number between 0 and 1, specifying how important is genetics 
	String input_candidate_genes_file;
	int num_pheno_rounds;
	int num_subject;
	double liability_threshold;
	
	VariantsDouble genotype;  // genotype in HDF5 format
	
	/*
	 * constructor 
	 */
	public SimHeterogeneity(String genotype_hdf5, String output_var_file, int num_rounds, 
			double min_maf, double max_maf, int num_causal_regions, int num_causal_var_per_region, 
			double genetic_component, double liability_threshold, String input_candidate_genes_file) {
		this.input_candidate_genes_file=input_candidate_genes_file;
		this.genetic_component=genetic_component;
		this.liability_threshold=liability_threshold;
		this.genotype=new VariantsDouble(genotype_hdf5);
		this.region2chr=new HashMap<Integer, Integer>();
		this.candidate_regions_locations=SimHeterogeneity.read_candidate_regions(input_candidate_genes_file); //chr_index, start_locus_location, end_locus_location
		this.num_candidate_regions=candidate_regions_locations.length;
		this.num_pheno_rounds=num_rounds;
		this.num_subject=this.genotype.sample_size;
		this.sim_causal_var(output_var_file, num_rounds, min_maf, max_maf, num_causal_regions, 
				num_causal_var_per_region);
	}
	
	/*
	 * Based on the causal variants, simulate quantitative phenotype using additive model (cross-region) 
	 * and compound heterozygotes model (within-region)  
	 */

	public double[][] sim_quantitative_phenotype(String causal_var_file, String simulation_report){
		double[][] pheno_matrix=new double[this.num_pheno_rounds][this.num_subject];
		try {
			BufferedReader br=new BufferedReader(new FileReader(causal_var_file));
			BufferedWriter bw_report=new BufferedWriter(new FileWriter(simulation_report));
			String line=br.readLine();
			while(line.startsWith("#"))line=br.readLine();
			int phen_round_index=0;
			while(line!=null) {
				String[] tmp=line.split("\t");
				if(tmp.length!=this.num_causal_regions) {
					System.out.println("ERROR: tmp.length ("+tmp.length+") "
							+ "!= this.num_causal_regions ("+this.num_causal_regions+")");
				}
				String[][] loci=new String[this.num_causal_regions][];
				for(int region_i=0;region_i<this.num_causal_regions;region_i++) {
					loci[region_i]=tmp[region_i].split(" ");
				}
				// scan each individual for each particular variant to add genetic contributions. 
				for(int region_i=0;region_i<this.num_causal_regions;region_i++) {
					int[] contributed_loci_count = new int[num_subject]; // if above 2, then no increase.
					Arrays.fill(contributed_loci_count, 0);
					for(int loci_i=0;loci_i<loci[region_i].length;loci_i++) {
						String[] the_locus_string=loci[region_i][loci_i].split(":");
						//int region_index= Integer.parseInt(the_locus_string[0]);
						int chr= Integer.parseInt(the_locus_string[1]);
						int loc_index=Integer.parseInt(the_locus_string[2]);
						double[] var_alleles=genotype.load_one_variant_by_index(chr, loc_index);
						int[] minor_major_allele= find_minor_major_allele(var_alleles);
						int minor_hom=minor_major_allele[0];  // could be 0 or 2
						int major_hom=minor_major_allele[1];  // could be 0 or 2
						for(int subj_i=0;subj_i<this.num_subject;subj_i++) {
							if(Double.compare(var_alleles[subj_i], major_hom)!=0) { //not major_allele hom, should contribute phenotype, if not full for the region
								if(contributed_loci_count[subj_i]==0) { // no contribution from this region yet
									if(Double.compare(var_alleles[subj_i], 1)==0) {  // a het
										pheno_matrix[phen_round_index][subj_i]+=pheno_unit;
										bw_report.write("Round"+phen_round_index+"\t"+genotype.sample_ids[subj_i]+"\t"+loci[region_i][loci_i]+"\t"+pheno_unit+"\n");
										contributed_loci_count[subj_i]=1;
									}else if(Double.compare(var_alleles[subj_i], minor_hom)==0){  // minor hom
										pheno_matrix[phen_round_index][subj_i]+=2*pheno_unit;
										bw_report.write("Round"+phen_round_index+"\t"+genotype.sample_ids[subj_i]+"\t"+loci[region_i][loci_i]+"\t"+(pheno_unit*2)+"\n");
										contributed_loci_count[subj_i]=2;
									}else {
										System.out.println("ERROR: var_alleles[subj_i] is neither 1 nor 2: "+ var_alleles[subj_i]+" (proceeded as if it is zero)");
									}
								}else if(contributed_loci_count[subj_i]==1) { // already one contribution
									pheno_matrix[phen_round_index][subj_i]+=pheno_unit;
									bw_report.write("Round"+phen_round_index+"\t"+genotype.sample_ids[subj_i]+"\t"+loci[region_i][loci_i]+"\t"+pheno_unit+"\n");
									contributed_loci_count[subj_i]=2;
								}else if(contributed_loci_count[subj_i]!=2){ // already 2; do nothing except for a error check
									System.out.println("ERROR: contributed_loci_count[subj_i]!=2");
								}
							}
						}
					}
				}
				Normalization.adding_random_component(pheno_matrix[phen_round_index], this.genetic_component); 
				phen_round_index++;
				line=br.readLine();
			}br.close();bw_report.close();
		}catch(Exception e) {e.printStackTrace();}
		return pheno_matrix;
	}
	
	/*
	 * find the major allele and minor allele from var_alleles;
	 * minor_major_allele[0] is the minor; minor_major_allele[1] is the major
	 */
	public int[] find_minor_major_allele(double[] var_alleles) {
		int[] minor_major_allele=new int[2];
		int zero_count=0, two_count=0;
		for(int i=0;i<var_alleles.length;i++) {
			if(Double.compare(var_alleles[i], 0)==0)zero_count++;
			else if(Double.compare(var_alleles[i], 2)==0)two_count++;
		}
		if(zero_count>two_count) {
			minor_major_allele[0]=2;
			minor_major_allele[1]=0;
		}else {
			minor_major_allele[0]=0;
			minor_major_allele[1]=2;
		}
		return minor_major_allele;
	}
	
	public void output_quanty_pheno_file(String out_phenotype_file, double[][] pheno_matrix) {
		try {
			BufferedWriter bw_p=new BufferedWriter(new FileWriter(out_phenotype_file));
			bw_p.write("ID");
			for(int phe_i=0;phe_i<this.num_pheno_rounds;phe_i++) {
				bw_p.write("\tPhenotype_"+phe_i);
			}bw_p.write("\n");
			for(int subj_i=0;subj_i<this.num_subject;subj_i++) {
				bw_p.write(this.genotype.sample_ids[subj_i]);
				for(int phe_i=0;phe_i<this.num_pheno_rounds;phe_i++) {
					bw_p.write("\t"+pheno_matrix[phe_i][subj_i]);
				}bw_p.write("\n");
			}bw_p.close();
		}catch(Exception e) {e.printStackTrace();}
	}
	
	public void output_binary_pheno_file(String out_phenotype_file, double[][] pheno_matrix) {
		try {
			BufferedWriter bw_p=new BufferedWriter(new FileWriter(out_phenotype_file));
			bw_p.write("ID");
			for(int phe_i=0;phe_i<this.num_pheno_rounds;phe_i++) {
				bw_p.write("\tPhenotype_"+phe_i);
			}bw_p.write("\n");
			for(int subj_i=0;subj_i<this.num_subject;subj_i++) {
				bw_p.write(this.genotype.sample_ids[subj_i]);
				for(int phe_i=0;phe_i<this.num_pheno_rounds;phe_i++) {
					bw_p.write("\t"+(pheno_matrix[phe_i][subj_i]>=this.liability_threshold?1:0));
				}bw_p.write("\n");
			}bw_p.close();
		}catch(Exception e) {e.printStackTrace();}
	}
	
	/*
	 * generate causal variants
	 * int[][] causal_regions: #regions x 3 (3 elements are: chr_index, start_locus_index, end_locus_index (all start from zero))
	 * when there is only one region, this is within-region heterogeneity. 
	 * 
	 * The output file contains multiple lines in the form:
	 * 
	 * region_index:chr_index:loc_index region_index:chr_index:loc_index\t...\tregion_index:chr_index:loc_index
	 * NOTE: loci in the same region are separated by space; but the regions are separated by "\t".
	 * 
	 * each line stands for causal variants for a future phenotype. 
	 * 
	 */
	public void sim_causal_var(String output_var_file, int num_rounds, double min_maf, double max_maf, 
			int num_causal_regions, int num_causal_var_per_region) {
		this.num_causal_regions=num_causal_regions;
		candidate_causal_indexes=new int[num_candidate_regions][]; //#region x #loci_per_region 
		// convert locations to index:
		for(int r=0;r<num_candidate_regions;r++) {
			int chr=candidate_regions_locations[r][0];
			region2chr.put(r, chr);
			int start_location=candidate_regions_locations[r][1];
			int end_location=candidate_regions_locations[r][2];
			int[] selected_indexes=genotype.find_vars_in_maf_range_in_region(chr, start_location, end_location, min_maf, max_maf);
			candidate_causal_indexes[r]=new int[selected_indexes.length];
			for(int i=0;i<selected_indexes.length;i++) {
				candidate_causal_indexes[r][i]=selected_indexes[i];
			}
		}
		try {
			BufferedWriter bw=new BufferedWriter(new FileWriter(output_var_file));
			bw.write("# Input gene region file: "+ input_candidate_genes_file+"\n");
			bw.write("#Region_index:chr_index:loc_index:location\n");
			if(num_candidate_regions<num_causal_regions) {
				System.out.println("ERROR: num_regions<num_causal_regions!");
				System.exit(0);
			}else {
				int num_regions_with_enough_loci=0;
				for(int r=0;r<num_candidate_regions;r++) {
					if(candidate_causal_indexes[r].length>=num_causal_var_per_region)
						num_regions_with_enough_loci++;
				}if(num_regions_with_enough_loci<num_causal_regions) {
					System.out.println("ERROR: num_regions_with_enough_loci ("+num_regions_with_enough_loci+")"
							+ " < num_causal_regions ("+num_causal_regions+")!");
					System.exit(0);
				}
			}
			for(int round_i=0;round_i<num_rounds;round_i++) {
				HashSet<Integer> used_regions=new HashSet<Integer>();
				for(int region_i=0;region_i<num_causal_regions;region_i++) {
					int the_region=(int)(Test.randomNumber()*num_candidate_regions);
					while(used_regions.contains(the_region) || candidate_causal_indexes[the_region].length<num_causal_var_per_region) {
						the_region=(int)(Test.randomNumber()*num_candidate_regions);
					}
					// found a region that has enough candidates and hasn't been used; now assign variants
					HashSet<Integer> used_loci=new HashSet<Integer>();
					for(int v_index=0;v_index<num_causal_var_per_region;v_index++) {
						int the_locus=(int)(Test.randomNumber()*candidate_causal_indexes[the_region].length);
						while(used_loci.contains(the_locus)) {
							the_locus=(int)(Test.randomNumber()*candidate_causal_indexes[the_region].length);
						}
						if(v_index!=num_causal_var_per_region-1) bw.write(the_region+":"+region2chr.get(the_region)+":"+candidate_causal_indexes[the_region][the_locus]+
								":"+genotype.locations[region2chr.get(the_region)][candidate_causal_indexes[the_region][the_locus]]+" ");
						// the last locus use "\t" to indicate different regions.
						else bw.write(the_region+":"+region2chr.get(the_region)+":"+candidate_causal_indexes[the_region][the_locus]+
								":"+genotype.locations[region2chr.get(the_region)][candidate_causal_indexes[the_region][the_locus]]+"\t");
						used_loci.add(the_locus);
					}
					used_regions.add(the_region); // this one has been used now
				}
				bw.write("\n");
			}
			bw.close();			
		}catch(Exception e) {e.printStackTrace();}
	}

	/*
	 * read in causal gene file. The file can contain headers starting with "#"
	 * 
	 * the rest lines are in the form 
	 * chr_index, start_locus_location, end_locus_location
	 */
	public static int[][] read_candidate_regions(String causal_genes_file){
		HashMap<String, Integer>  chr_string2index  = generate_chr_map();
		ArrayList<String[]> regions=new ArrayList<String[]>();
		try {
			BufferedReader br=new BufferedReader(new FileReader(causal_genes_file));
			String line=br.readLine();
			while(line.startsWith("#"))line=br.readLine();
			while(line!=null) {
				String[] tmp=line.split("\t");
				regions.add(tmp);
				line=br.readLine();
			}br.close();
		}catch(Exception e) {e.printStackTrace();}
		int num_regions=regions.size();
		int[][] causal_vars=new int[num_regions][3];
		for(int k=0;k<num_regions;k++) {
			String[] region_string=regions.get(k);
			causal_vars[k][0]= chr_string2index.get(region_string[0])-1; // indexes start with ZERO!
			causal_vars[k][1]= Integer.parseInt(region_string[1]);
			causal_vars[k][2]= Integer.parseInt(region_string[2]);
		}
		return causal_vars;
	}
	
	public static HashMap<String, Integer> generate_chr_map(){
		HashMap<String, Integer> chr2index=new HashMap<String, Integer>();
		for(int k=1;k<=25;k++){
			chr2index.put("chr"+k, k);
			chr2index.put("Chr"+k, k);
			chr2index.put(""+k, k);
		}
		chr2index.put("X", 23);chr2index.put("chrX", 23);chr2index.put("ChrX", 23);
		chr2index.put("Y", 24);chr2index.put("chrY", 24);chr2index.put("ChrY", 24);
		chr2index.put("MT", 25);chr2index.put("chrMT", 25);chr2index.put("ChrMT", 25);
		return chr2index;
	}

	
	public static void main(String[] args) {
		if(args.length==0){
			System.out.println("===================================================================================\n");
			System.out.println("Simulating phenotype based on level of heterogeneity" +
					"Developer: Quan LONG Feb 2019\n" +
					"Usage: java -Xmx4g -jar heterogeneity.jar [function]");
			System.out.println("Supported functions:" +
					"\n\tsimulation" +
					"\n");
			System.exit(0);
		}
		String function=args[0];
		int number_of_parameter=10;
		if(function.equals("simulation")){
			if(args.length<number_of_parameter*2+1){
				System.out.println("Simulate variants that are causal for the phenotypic changes.");
				System.out.println("Usage: \n\t"+
						"<-genotypes\tinput_genotype_file(HDF5-format)>\n\t" +
						"<-out_var\toutput_file_recording_causal_variants>\n\t" +
						"<-out_pheno\toutput_file_for_phenotypes>\n\t" +
						"<-sim_report\tsimulation_report_file>\n\t"+
						"<-in_region_file\tinput_candidate_genes_file>\n\t" +
						"<-num_pheno_rounds\tnum_of_rounds_to_simulate>\n\t" +
						"<-num_causal_regions\tnum_of_causal_regions>\n\t" +
						"<-num_causal_var_per_region\tnum_of_causal_variants_per_region>\n\t" +
						"<-min_maf\tminimal_num_MAF_of_selected_variants>\n\t" +
						"<-max_maf\tmaximal_num_MAF_of_selected_variants>\n\t" +	
						"<-genetic_component\tgenetic_component_(a_number_between_0_and_1)>\n\t"+ 	
						"[-liability_cutoff\tliability_cutoff_above_which_it_is_a_case]\n\t"+
						"[-binary\twhether_it_is_quantitative_or_binary]\n\t\n"+
						"There are "+number_of_parameter+" mandatory parameters.\nBut you have only specified "+(args.length-1)/2+" parameters"
					);
				System.exit(0);
			}else{				
				String genotype_hdf5=null; 
				String output_var_file=null; 
				String output_pheno_file=null;
				String candidate_genes_file=null;
				String simulation_report_file=null;
				int num_rounds=10;
				int num_causal_regions=2;
				int num_causal_var_per_region=2;
				double min_maf=0.005; 
				double max_maf=0.05; 
				double genetic_component=0.7; 
				boolean quantitative = true;
				double liability_cutoff=Double.NaN;				
				for(int k=1;k<args.length;k++){
					if(args[k].startsWith("-")){
						if(args[k].equals("-genotypes"))genotype_hdf5=args[k+1];
						else if(args[k].equals("-out_var"))output_var_file=args[k+1];
						else if(args[k].equals("-out_pheno"))output_pheno_file=args[k+1];
						else if(args[k].equals("-in_region_file"))candidate_genes_file=args[k+1];
						else if(args[k].equals("-sim_report"))simulation_report_file=args[k+1];
						else if(args[k].equals("-num_pheno_rounds"))num_rounds=Integer.parseInt(args[k+1]);
						else if(args[k].equals("-num_causal_regions"))num_causal_regions=Integer.parseInt(args[k+1]);
						else if(args[k].equals("-num_causal_var_per_region"))num_causal_var_per_region=Integer.parseInt(args[k+1]);
						else if(args[k].equals("-min_maf"))min_maf=Double.parseDouble(args[k+1]);
						else if(args[k].equals("-max_maf"))max_maf=Double.parseDouble(args[k+1]);
						else if(args[k].equals("-genetic_component"))genetic_component=Double.parseDouble(args[k+1]);
						else if(args[k].equals("-liability_cutoff"))liability_cutoff=Double.parseDouble(args[k+1]);
						else if(args[k].equals("-binary"))quantitative=false;
						else {
							System.out.println("Sorry, the option "+args[k]+" is not supported. A typo?");
							System.exit(0);
						}
					}
				}
				SimHeterogeneity sim_het=new SimHeterogeneity(genotype_hdf5, output_var_file, num_rounds, 
						min_maf, max_maf, num_causal_regions, num_causal_var_per_region, 
						genetic_component, liability_cutoff, candidate_genes_file);
				double[][] pheno_matrix=sim_het.sim_quantitative_phenotype(output_var_file, simulation_report_file);
				if(quantitative) {
					sim_het.output_quanty_pheno_file(output_pheno_file, pheno_matrix);
				}else {
					sim_het.output_binary_pheno_file(output_pheno_file, pheno_matrix);
				}
				
			}
		}
		else {
			System.out.println("Sorry, the function "+args[0]+" doesn't exist");
		}

	}

}
