package simulations;


import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;


import mixedmodel.VariantsDouble;
import myMathLib.Normalization;
import myMathLib.Test;
import supportfunc.PairArrays;
import supportfunc.QuickSortArr;

public class SimHiC {
	
	/*
	 * variables
	 */
	static double pheno_unit=1.0;
	VariantsDouble genotype; //genotype in HDF5 format
	HashMap<Integer, Integer> compound2chr;
	HashMap<Integer, String> compound2region1;
	HashMap<Integer, String> compound2region2;
	int num_pheno_rounds;
	int num_candidate_compounds;
	double genetic_component; // a number between 0 and 1, specifying how important is genetics 
	String input_candidate_genes_file;
	int[][][] candidate_regions_locations;
	double liability_threshold;
	int num_subject;
	int num_causal_var_per_region;
	int flank_causal_length;
	int num_causal_compounds;
	int[][][] candidate_causal_variants_indexes;
	

	public SimHiC(String genotype_hdf5, String output_var_file, String output_compounds_file, int num_rounds, 
			double min_maf, double max_maf, int num_causal_compounds, int num_causal_var_per_region, 
			double genetic_component, double liability_threshold, 
			String input_candidate_genes_file, String type, int flank_causal_length) {
		
		this.genotype=new VariantsDouble(genotype_hdf5);
		this.compound2chr=new HashMap<Integer, Integer>();
		this.compound2region1=new HashMap<Integer, String>();
		this.compound2region2=new HashMap<Integer, String>();
		this.num_pheno_rounds=num_rounds;
		this.genetic_component=genetic_component;
		this.input_candidate_genes_file=input_candidate_genes_file;
		this.candidate_regions_locations = SimHiC.read_candidate_regions_in_compound(input_candidate_genes_file, type);
		this.num_candidate_compounds = candidate_regions_locations.length;
		this.num_subject=this.genotype.sample_size;
		this.num_causal_var_per_region = num_causal_var_per_region;
		this.liability_threshold = liability_threshold;
		this.sim_causal_var
		(output_var_file, output_compounds_file, num_rounds,  
				min_maf, max_maf, num_causal_compounds,flank_causal_length);
	}
	
	public static int[][][] read_candidate_regions_in_compound(String causal_gene_file, String type){
		/*
		 * sort the input causal gene file and remove duplication, output arr contains R1_R2
		 * causal_gene_file content:  c1 \t 1;100;200 \t 1;250;350
		 * after sort, causal_gene_content_Arr [1;100;200_1;250;350,...]
		 */
		String[] causal_gene_content_Arr = QuickSortArr.content2Arr(type, causal_gene_file);// the causal_gene has already been sorted
//		QuickSortArr.quicksort(causal_gene_content_Arr, 0, causal_gene_content_Arr.length-1, type);
		/*
		 *  causal_regions_index
		 *  ----R1--R2-----
		 *  C1 |[1,100,200], [1,250,350]
		 *  C2 |[], []
		 *  remove compounds that have less than 2 regions and don't have 2 regions on the same chromosome
		 */
		int num_compound = causal_gene_content_Arr.length;
		int [][][] causal_regions_index = new int[num_compound][2][3];
		ArrayList<Integer> index_compound_with_regions_in_same_chr= new ArrayList<Integer>();
		for (int c=0; c< num_compound;c++) {
			String[] regions4compound = causal_gene_content_Arr[c].split("_");
			if(regions4compound.length != 2) {
				System.out.print("Compound"+(c+1)+"does not have 2 regions, skip this compound");
				continue;
			}else {
				String[] region1 = regions4compound[0].split(";");
				String[] region2 = regions4compound[1].split(";");
				int region1_chr = Integer.parseInt(region1[0])-1; //chr index start with 0
				int region1_start = Integer.parseInt(region1[1]);
				int region1_end = Integer.parseInt(region1[2]);
				int region2_chr = Integer.parseInt(region2[0])-1;
				int region2_start = Integer.parseInt(region2[1]);
				int region2_end = Integer.parseInt(region2[2]);
				if(region1_chr != region2_chr) {
					System.out.println("The 2 regions from this compound come from different chr, skip this compound");
					continue;
				}else {
					index_compound_with_regions_in_same_chr.add(c);
					causal_regions_index[c][0][0] = region1_chr;
					causal_regions_index[c][0][1] = region1_start;
					causal_regions_index[c][0][2] = region1_end;
					causal_regions_index[c][1][0] = region2_chr;
					causal_regions_index[c][1][1] = region2_start;
					causal_regions_index[c][1][2] = region2_end;
				}
			}
		}
		//if some compounds are skipped, copy those not been skipped compound index to a new Arr; 
		//else, return causal_regions_index
		if(index_compound_with_regions_in_same_chr.size() != num_compound) {
			int [][][] causal_regions_index2 = new int[index_compound_with_regions_in_same_chr.size()][2][3];
			for(int c_target=0; c_target<index_compound_with_regions_in_same_chr.size();c_target++) {
				causal_regions_index2[c_target]=causal_regions_index[index_compound_with_regions_in_same_chr.get(c_target)];
			}
			return causal_regions_index2;
		}else {
			return causal_regions_index;
		}
	}
	
	/*
	 * generate causal variants
	 * int[][][] causal_regions: #compounds x 2 x 3 (3 elements are : chr_index, start_location_index, end_location_index)
	 * the output file contains multiple lines in the form:

	 * NOTE: loci in the same compound are separated by space, but compounds are separated by \t
	 * e.g. 1 variant for each part of one compound
	 * compound_index:chr_index:start_region1:end_region1:locus_index:location\s
	 * compound_index:chr_index:start_region2:end_region2:locus_index:location\s 
	 * compound_index:chr_index:start_region1:end_region1:locus_index:location\s
	 * compound_index:chr_index:start_region2:end_region2:locus_index:location\t
	 * 
	 * each line stands for causal variants for a future phenotype.
	 */
	public void sim_causal_var(String output_var_file, String output_compounds_file, int num_rounds, 
			double min_maf, double max_maf, 
			int num_causal_compounds, int flank_causal_length) {
		this.num_causal_compounds = num_causal_compounds;
		this.candidate_causal_variants_indexes = new int[this.num_candidate_compounds][2][];
		//find the index of SNPs from genotype file that has the frequency belongs to min_maf, max_maf,
		//put those indexes into an array then store them in the candidate_causal_variants_indexes
		for(int c=0;c<this.num_candidate_compounds;c++) {
			//region 1 from compound c
			int chr_region1 = candidate_regions_locations[c][0][0];//chr starts from 0. 0-21
			compound2chr.put(c, chr_region1);// can easily get the chr for each compound. Note one compound should locate at the same chr
			int start_location_region1 = candidate_regions_locations[c][0][1]-flank_causal_length; //location of start, not index of start 
			int end_location_region1 = candidate_regions_locations[c][0][2]+flank_causal_length; // location of end, not index of end
			String region1 = Integer.toString(start_location_region1)+":"+Integer.toString(end_location_region1);
			compound2region1.put(c, region1);
			int[] selected_var_indexes_region1 = genotype.find_vars_in_maf_range_in_region(chr_region1, 
					start_location_region1, end_location_region1, min_maf, max_maf);
			this.candidate_causal_variants_indexes[c][0]=new int[selected_var_indexes_region1.length];
			for(int SNP_index=0; SNP_index< selected_var_indexes_region1.length;SNP_index++) {
//				System.out.println(c+ " selected_var_indexes_regions1 ["+SNP_index+"] "+selected_var_indexes_region1[SNP_index]);
				this.candidate_causal_variants_indexes[c][0][SNP_index]=selected_var_indexes_region1[SNP_index];
			}
			//region 2 from compound c
			int chr_region2 = candidate_regions_locations[c][1][0];//chr starts from 0
			int start_location_region2 = candidate_regions_locations[c][1][1]-flank_causal_length;
			int end_location_region2 = candidate_regions_locations[c][1][2]+flank_causal_length;
			String region2 = Integer.toString(start_location_region2)+":"+Integer.toString(end_location_region2);
			compound2region2.put(c, region2);
			int[] selected_var_indexes_region2 = genotype.find_vars_in_maf_range_in_region(chr_region2, 
					start_location_region2, end_location_region2, min_maf, max_maf);
			this.candidate_causal_variants_indexes[c][1]=new int[selected_var_indexes_region2.length];
			for(int SNP_index=0; SNP_index< selected_var_indexes_region2.length; SNP_index++) {
//				System.out.println(c+ " selected_var_indexes_regions2 ["+SNP_index+"] "+selected_var_indexes_region2[SNP_index]);
				this.candidate_causal_variants_indexes[c][1][SNP_index] = selected_var_indexes_region2[SNP_index];
			}
		}
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(output_var_file));
			BufferedWriter bw_compounds_report = new BufferedWriter(new FileWriter(output_compounds_file));
			bw.write("# Input gene region file: "+ input_candidate_genes_file+"\n");
			bw.write("# compound_index:chr_index:start_location:end_location:SNP_index:SNP_location compound_index:chr_index:start_location:end_location:SNP_index:SNP_location \n");
			if(this.num_candidate_compounds < num_causal_compounds) {
				System.out.print("Error: too many causal compounds"); //should never be executed, num_causal_c=1/2, num_candidate_c~8000 
			}else {
				int num_compounds_with_enough_loci=0;
				for(int c=0; c< this.num_candidate_compounds;c++) {
					if(this.candidate_causal_variants_indexes[c][0].length>=this.num_causal_var_per_region && 
							this.candidate_causal_variants_indexes[c][1].length>=this.num_causal_var_per_region) {
						num_compounds_with_enough_loci++;
					}
				}
				if(num_compounds_with_enough_loci < num_causal_compounds) {
					System.out.println("Error: num_compounds_with_enough_loci "+num_compounds_with_enough_loci
							+" < num_causal_vompound "+ num_causal_compounds);
					System.exit(0);
				}else {
					System.out.println("The num of compounds with enough loci is "+num_compounds_with_enough_loci);
				}
			 }
			/*
			 * one compound cannot be used twice in one simulation round, but different simulation round can share the same compound
			 */
			for(int round_i=0;round_i<num_rounds;round_i++) { 
				HashSet<Integer> used_compounds = new HashSet<Integer>();
				for(int compound_i=0; compound_i<num_causal_compounds;compound_i++) {
					//select a random index of the region
					int the_compound = (int) (Test.randomNumber()*this.num_candidate_compounds);//[0,num_candidate_compounds.length-1]
					//because some compounds may only have one region as causal, in addition, there may only have one SNP in that region
					//in that case, select another compound
					while(used_compounds.contains(the_compound)||
							this.candidate_causal_variants_indexes[the_compound][0].length<this.num_causal_var_per_region||
							this.candidate_causal_variants_indexes[the_compound][1].length<this.num_causal_var_per_region) {
						the_compound = (int) (Test.randomNumber()*this.num_candidate_compounds);
					}
					//found a compound that has enough candidate causal SNPs and has not been used; now assign SNPs
					HashSet<Integer> used_loci1 = new HashSet<Integer>();
					HashSet<Integer> used_loci2 = new HashSet<Integer>();
					for(int SNP=0; SNP<this.num_causal_var_per_region; SNP++) {
						int the_locus_1 = (int)(Test.randomNumber()*this.candidate_causal_variants_indexes
								[the_compound][0].length);
						int the_locus_2 = (int)(Test.randomNumber()*this.candidate_causal_variants_indexes
								[the_compound][1].length);
						/*
						 * if used_loci1.contains(the_locus_1)&&used_loci2.contains(the_locus_2), then one SNP in one region may be used twice
						 * else if used_loci1.contains(the_locus_1)||used_loci2.contains(the_locus_2), no SNP in one region will be used twice
						 */
						while(used_loci1.contains(the_locus_1)||used_loci2.contains(the_locus_2)) {  
							the_locus_1 = (int)(Test.randomNumber()*this.candidate_causal_variants_indexes
									[the_compound][0].length);
							the_locus_2 = (int)(Test.randomNumber()*this.candidate_causal_variants_indexes
									[the_compound][1].length);
						}
						if(SNP != this.num_causal_var_per_region-1) {
							bw.write(the_compound+":"+(compound2chr.get(the_compound)+1)+":"+
									compound2region1.get(the_compound)+":"
									+this.candidate_causal_variants_indexes[the_compound][0][the_locus_1]+":"
									+genotype.locations[compound2chr.get(the_compound)][this.candidate_causal_variants_indexes[the_compound][0][the_locus_1]]
									+" "+
									the_compound+":"+(compound2chr.get(the_compound)+1)+":"+compound2region2.get(the_compound)+":"
									+this.candidate_causal_variants_indexes[the_compound][1][the_locus_2]+":"
									+genotype.locations[compound2chr.get(the_compound)][this.candidate_causal_variants_indexes[the_compound][1][the_locus_2]]
									+" ");
						}else {
							bw.write(the_compound+":"+(compound2chr.get(the_compound)+1)+":"+
									compound2region1.get(the_compound)+":"
									+this.candidate_causal_variants_indexes[the_compound][0][the_locus_1]+":"
									+genotype.locations[compound2chr.get(the_compound)][this.candidate_causal_variants_indexes[the_compound][0][the_locus_1]]
									+" "+
									the_compound+":"+(compound2chr.get(the_compound)+1)+":"+compound2region2.get(the_compound)+":"
									+this.candidate_causal_variants_indexes[the_compound][1][the_locus_2]+":"
									+genotype.locations[compound2chr.get(the_compound)][this.candidate_causal_variants_indexes[the_compound][1][the_locus_2]]
									+"\t");
						}
						used_loci1.add(the_locus_1);
						used_loci2.add(the_locus_2);
					}
					used_compounds.add(the_compound);
					bw_compounds_report.write((candidate_regions_locations[the_compound][0][0]+1)+":" //chr_index starts from 0
												+candidate_regions_locations[the_compound][0][1]+":"
												+candidate_regions_locations[the_compound][0][2]+"_"
												+(candidate_regions_locations[the_compound][1][0]+1)+":"
												+candidate_regions_locations[the_compound][1][1]+":"
												+candidate_regions_locations[the_compound][1][2]+"\t");
				}
				bw.write("\n");
				bw_compounds_report.write("\n");
			}
			bw.close();
			bw_compounds_report.close();
		}catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	public double[][] sim_quantitative_phenotype(String causal_var_file, String simulation_report, 
			String pattern, String sim_pheno){
		double[][] pheno_matrix = new double[this.num_pheno_rounds][this.num_subject];
		//initially filled up with 0.0
		double[][] real_pheno_matrix = new double[this.num_pheno_rounds][this.num_subject];
		try {
			BufferedReader br = new BufferedReader(new FileReader(causal_var_file));
			BufferedWriter bw_report = new BufferedWriter(new FileWriter(simulation_report+"."+pattern));
		    String line = br.readLine();
		    while(line.startsWith("#")) line=br.readLine();
		    int pheno_round_index = 0;
		    while(line != null) {
		    	String[] causal_compounds = line.trim().split("\t");
		    	if(causal_compounds.length != this.num_causal_compounds) {
		    		System.out.print("Error: causal_compounds.length ("+causal_compounds.length+")"
		    	+" != this.num_causal_compounds ("+this.num_causal_compounds+")");
		    	}
		    	String[][] loci = new String[this.num_causal_compounds][];
		    	for(int c_i=0;c_i<causal_compounds.length;c_i++) {
		    		loci[c_i]=causal_compounds[c_i].split(" ");
		    	}
		    	//scan each individual for each particular variant to add genetic contributions based on patterns
		    	for(int c_i=0;c_i<causal_compounds.length;c_i++) {
//		    		int[] contributed_loci_count = new int[this.num_subject];//
//		    		Arrays.fill(contributed_loci_count,0);
		    		int SNPs_in_each_region = loci[c_i].length/2;
		    		if(SNPs_in_each_region !=this.num_causal_var_per_region) {
		    			System.out.println("The number of causal SNPs per region from output var file does "
		    					+ "not match the input num_causal_var_per_region\nCheck your ourpur var file again");
		    			System.exit(0);
		    		}
		    		String[] Region1_loci = new String[SNPs_in_each_region];
		    		String[] Region2_loci = new String[SNPs_in_each_region];
		    		for(int i=0;i<SNPs_in_each_region;i++) {
		    			Region1_loci[i]=loci[c_i][2*i];//2;200;400
		    			Region2_loci[i]=loci[c_i][2*i+1];//3;400;500
		    		}
		    		//get within regional phenotype for each individual 
		    		PairArrays region1_pheno_report_string = define_region_pheno(Region1_loci);
		    		String[][] region1_report_string= region1_pheno_report_string.getStingArr();
		    		double[] region1_pheno = region1_pheno_report_string.getDoubleArr(); 
		    		PairArrays region2_pheno_report_string = define_region_pheno(Region2_loci);
		    		String[][] region2_report_string = region2_pheno_report_string.getStingArr();
		    		double[] region2_pheno = region2_pheno_report_string.getDoubleArr();
		    		
		    		
		    		//get cross-region phenotype for each individual
		    		if(pattern.equals("add")) {
			    		for(int subj_i=0;subj_i<this.num_subject;subj_i++) {
			    			pheno_matrix[pheno_round_index][subj_i]=region1_pheno[subj_i]+region2_pheno[subj_i];
			    			for(int SNP_contribute=0; SNP_contribute< region1_report_string[subj_i].length;SNP_contribute++) {
			    				if(region1_report_string[subj_i][SNP_contribute]!=null) {
//				    				System.out.println("Round"+pheno_round_index+"\t"+region1_report_string[subj_i][SNP_contribute]);
				    				bw_report.write("Round"+pheno_round_index+"\t"+region1_report_string[subj_i][SNP_contribute]+"\n");
			    				}
			    			}
			    			for(int SNP_contribute=0; SNP_contribute< region2_report_string[subj_i].length;SNP_contribute++) {
			    				if(region2_report_string[subj_i][SNP_contribute]!=null) {
			    					bw_report.write("Round"+pheno_round_index+"\t"+region2_report_string[subj_i][SNP_contribute]+"\n");
			    				}
			    			}
			    		}
		    		}else if(pattern.equals("hetero")) {//or
			    		for(int subj_i=0;subj_i<this.num_subject;subj_i++) {
			    			if((Double.compare(region1_pheno[subj_i], 0)!=0) || (Double.compare(region2_pheno[subj_i], 0)!=0)) {
			    				pheno_matrix[pheno_round_index][subj_i]=pheno_unit;
			    			}
			    			bw_report.write("=====hetero model, the following SNPs from 2 regions contribtue for 1 pheno unit in this individual==========================\n");
			    			for(int SNP_contribute=0; SNP_contribute< region1_report_string[subj_i].length;SNP_contribute++) {
			    				if(region1_report_string[subj_i][SNP_contribute]!=null) {
			    					bw_report.write("Round"+pheno_round_index+"\t"+region1_report_string[subj_i][SNP_contribute]+"\n");
			    				}
			    			}
			    			for(int SNP_contribute=0; SNP_contribute< region2_report_string[subj_i].length;SNP_contribute++) {
			    				if(region2_report_string[subj_i][SNP_contribute]!=null) {
			    					bw_report.write("Round"+pheno_round_index+"\t"+region2_report_string[subj_i][SNP_contribute]+"\n");
			    				}
			    			}
			    		}
		    		}else if(pattern.equals("both")) {//and 
			    		for(int subj_i=0;subj_i<this.num_subject;subj_i++) {
			    			if((Double.compare(region1_pheno[subj_i], 0)!=0) && 
			    					(Double.compare(region2_pheno[subj_i], 0)!=0)) {
			    				pheno_matrix[pheno_round_index][subj_i]=pheno_unit;
			    			}
			    			bw_report.write("=====both model, if both regions have causality, then this compound contribtues for 1 pheno unit in this individual==========\n");
			    			for(int SNP_contribute=0; SNP_contribute< region1_report_string[subj_i].length;SNP_contribute++) {
			    				if(region1_report_string[subj_i][SNP_contribute]!=null) {
			    					bw_report.write("Round"+pheno_round_index+"\t"+region1_report_string[subj_i][SNP_contribute]+"\n");
			    				}
			    			}
			    			for(int SNP_contribute=0; SNP_contribute< region2_report_string[subj_i].length;SNP_contribute++) {
			    				if(region2_report_string[subj_i][SNP_contribute]!=null) {
			    					bw_report.write("Round"+pheno_round_index+"\t"+region2_report_string[subj_i][SNP_contribute]+"\n");
			    				}
			    			}
			    		}
		    		}else if(pattern.equals("compensate")) {
			    		for(int subj_i=0;subj_i<this.num_subject;subj_i++) {
			    			//if region1_pheno[subj_i] equals to region2_pheno[subj_i], no contribution to pheno type
			    			if(Double.compare(region1_pheno[subj_i], region2_pheno[subj_i])!=0){
			    				double diff = Math.abs(region1_pheno[subj_i]-region2_pheno[subj_i]);
			    				if(Double.compare(diff,1.0)==0) {
			    					pheno_matrix[pheno_round_index][subj_i]=pheno_unit;
			    				}else if(Double.compare(diff, 2.0)==0) {
			    					pheno_matrix[pheno_round_index][subj_i]=2*pheno_unit;
			    				}
			    			}
			    			bw_report.write("=====Compensate model, one unit pheno in one region will compensate one unite pheno in another region========================\n");
			    			for(int SNP_contribute=0; SNP_contribute< region1_report_string[subj_i].length;SNP_contribute++) {
			    				if(region1_report_string[subj_i][SNP_contribute]!=null) {
			    					bw_report.write("Round"+pheno_round_index+"\t"+region1_report_string[subj_i][SNP_contribute]+"\n");
			    				}
			    			}
			    			for(int SNP_contribute=0; SNP_contribute< region2_report_string[subj_i].length;SNP_contribute++) {
			    				if(region2_report_string[subj_i][SNP_contribute]!=null) {
			    					bw_report.write("Round"+pheno_round_index+"\t"+region2_report_string[subj_i][SNP_contribute]+"\n");
			    				}
			    			}
			    		}
		    		}else {
		    			System.out.println("No such "+pattern+" pattern");
		    		}
		    	}
		    	//Normalization the real phenotype with a random residual
		    	real_pheno_matrix[pheno_round_index] = pheno_matrix[pheno_round_index].clone();
		    	Normalization.adding_random_component(real_pheno_matrix[pheno_round_index], this.genetic_component);
		    	pheno_round_index++;
		    	line=br.readLine();
		    }
		    
		    //output the original phenotype without residual
			BufferedWriter bw_ori_pheno = new BufferedWriter(new FileWriter(sim_pheno+".org.pheno"));
	    	bw_ori_pheno.write("ID");;
			for(int phe_i=0;phe_i<this.num_pheno_rounds;phe_i++) {
				bw_ori_pheno.write("\tPhenotype_round_"+phe_i);
			}bw_ori_pheno.write("\n");
			for(int subj_i=0;subj_i<this.num_subject;subj_i++) {
				bw_ori_pheno.write(this.genotype.sample_ids[subj_i]);
				for(int phe_i=0;phe_i<this.num_pheno_rounds;phe_i++) {
					bw_ori_pheno.write("\t"+pheno_matrix[phe_i][subj_i]);
				}bw_ori_pheno.write("\n");
			}
			//output balached phenotype for hetero/both/compensate models
			if(!pattern.equals("add")) {
				ArrayList<Integer> balanced_pheno_index = balance_case_control_pheno_index(pheno_matrix, simulation_report);
				double[][] balanced_pheno_matrix = new double[balanced_pheno_index.size()][];
				if(balanced_pheno_index.size()!=0) {//if control/case belongs to [0.9,1.1]
					for(int b_phe_i=0; b_phe_i<balanced_pheno_index.size();b_phe_i++) {
						balanced_pheno_matrix[b_phe_i] = real_pheno_matrix[balanced_pheno_index.get(b_phe_i)];
					}
				}
				BufferedWriter bw_balance_pheno = new BufferedWriter(new FileWriter(sim_pheno+".balance.pheno"));
		    	bw_balance_pheno.write("ID");;
				for(int b_phe_i=0; b_phe_i<balanced_pheno_index.size();b_phe_i++) {
					bw_balance_pheno.write("\tPhenotype_round_"+balanced_pheno_index.get(b_phe_i));
				}bw_balance_pheno.write("\n");
				for(int subj_i=0;subj_i<this.num_subject;subj_i++) {
					bw_balance_pheno.write(this.genotype.sample_ids[subj_i]);
					for(int b_phe_i=0; b_phe_i<balanced_pheno_index.size();b_phe_i++) {
						bw_balance_pheno.write("\t"+balanced_pheno_matrix[b_phe_i][subj_i]);
					}bw_balance_pheno.write("\n");
				}
				bw_balance_pheno.close();
			}
			
			br.close();
			bw_report.close();
			bw_ori_pheno.close();
		}catch(Exception e) {
			e.printStackTrace();
		}
		return real_pheno_matrix;
	}

	public static ArrayList<Integer> balance_case_control_pheno_index(
			double [][] pheno_matrix,String sim_report) {
		ArrayList<Integer> balanced_pheno_index = new ArrayList<Integer>();
		int[][] case_control_count = new int[pheno_matrix.length][2]; //[0] to count 0 in the ith pheno, [1] to count not 0 in the ith p
		for(int phe_i=0;phe_i<pheno_matrix.length;phe_i++) {
			for(int subj_i=0;subj_i<pheno_matrix[phe_i].length;subj_i++) {
				if(Double.compare(pheno_matrix[phe_i][subj_i], 0)==0){
					case_control_count[phe_i][0]+=1;
				}else {
					case_control_count[phe_i][1]+=1;
				}
			}
		}
		try{
			BufferedWriter case_control_percent = new BufferedWriter(new FileWriter(sim_report+".case.control.%.txt"));
			case_control_percent.write("Phenotype round\tcontrol percentage\n");
			for(int phe_i=0;phe_i<pheno_matrix.length;phe_i++) {
				double control_D_case = case_control_count[phe_i][0]/(1.0*case_control_count[phe_i][1]);
				case_control_percent.write("Phenotype_round_"+phe_i+"\t"+control_D_case+"\n");
				if(Double.compare(control_D_case, 1.1)<=0 && Double.compare(control_D_case, 0.9)>=0) {
					balanced_pheno_index.add(phe_i);
				}
			}
			case_control_percent.close();
		}catch (Exception e) {
			// TODO: handle exception
			e.printStackTrace();
		}
		return balanced_pheno_index;
	}
	
	//this function define within region causality
	public PairArrays define_region_pheno(String[] region_causal_loci) {
		double[] region_pheno = new double[this.num_subject];
		String[][] region_report_string = new String[this.num_subject][region_causal_loci.length];
		int[] contributed_loci_count = new int[this.num_subject]; // up to 2 mutations; 
		//1 mutation leads to an increase of phenotype by 1.0; 2 mutations lead to an increase of phenotype by 2.0;
		for(int SNP_i=0;SNP_i<region_causal_loci.length;SNP_i++) {
			String[] the_locus_string=region_causal_loci[SNP_i].split(":");
			int chr_index = Integer.parseInt(the_locus_string[1])-1; //!!!!!!chr_index starts from 0. However chr in output_var starts from 1
			int loc_index = Integer.parseInt(the_locus_string[4]);
			double[] var_alleles=genotype.load_one_variant_by_index(chr_index, loc_index);
			int[] minor_major_allele = find_minor_major_allele(var_alleles);
			int minor_homo = minor_major_allele[0];
			int major_homo = minor_major_allele[1];
			for(int subj_i=0;subj_i<this.num_subject;subj_i++) {
				if(Double.compare(var_alleles[subj_i], major_homo)!=0) {
					if(contributed_loci_count[subj_i]==0) {
						if(Double.compare(var_alleles[subj_i], 1)==0) {
							//het
							region_pheno[subj_i]+=pheno_unit;
							contributed_loci_count[subj_i]=1;
							region_report_string[subj_i][SNP_i]=(genotype.sample_ids[subj_i]+"\t"+region_causal_loci[SNP_i]+"\t"+pheno_unit);
						}else if(Double.compare(var_alleles[subj_i], minor_homo)==0) {
							//minor homo
							region_pheno[subj_i]+=2*pheno_unit;
							contributed_loci_count[subj_i]=2;
							region_report_string[subj_i][SNP_i]=(genotype.sample_ids[subj_i]+"\t"+region_causal_loci[SNP_i]+"\t"+(2*pheno_unit));
						}else {
							System.out.print("Error: var_alleles[subj_i] is neither hetero or minor homo: "
									+ var_alleles[subj_i]+" (processed as if it is major homo)");
						}
					}else if(contributed_loci_count[subj_i]==1) {
						region_pheno[subj_i]+=pheno_unit;
						contributed_loci_count[subj_i]=2;
						region_report_string[subj_i][SNP_i]=(genotype.sample_ids[subj_i]+"\t"+region_causal_loci[SNP_i]+"\t"+pheno_unit);
					}else if(contributed_loci_count[subj_i]!=2) {
						System.out.println("Error: contributed_loci_couont[subj_i]!=2");
					}
				}
			}
		}
		return new PairArrays(region_report_string, region_pheno);
	}
	
	
	//find the minor allele and assign to minor_major_allele[0]
	//find the major allele and assign to minor_major_allele[1]
	public int[] find_minor_major_allele(double[] var_alleles) {
		int[] minor_major_allele = new int[2];
		int zero_count = 0;
		int two_count =0;
		for(int i=0; i<var_alleles.length;i++) {
			if(Double.compare(var_alleles[i], 0)==0) {
				zero_count++;
			}else if(Double.compare(var_alleles[i], 2)==0) {
				two_count++;
			}
		}
		// if two_count==zero_count, we assign minor allele 0, major allele 2 
		if(zero_count>two_count) {
			minor_major_allele[0]=2;
			minor_major_allele[1]=0;
		}else {
			minor_major_allele[0]=0;
			minor_major_allele[1]=2;
		}
		return minor_major_allele;
	}
	
	public void output_quanty_pheno_file(String out_phenotype_file, double[][]pheno_matrix) {
		try {
			BufferedWriter bw_q= new BufferedWriter(new FileWriter(out_phenotype_file+".all.pheno"));
			bw_q.write("ID");;
			for(int phe_i=0;phe_i<this.num_pheno_rounds;phe_i++) {
				bw_q.write("\tPhenotype_round_"+phe_i);
			}bw_q.write("\n");
			for(int subj_i=0;subj_i<this.num_subject;subj_i++) {
				bw_q.write(this.genotype.sample_ids[subj_i]);
				for(int phe_i=0;phe_i<this.num_pheno_rounds;phe_i++) {
					bw_q.write("\t"+pheno_matrix[phe_i][subj_i]);
				}bw_q.write("\n");
			}
			bw_q.close();
		}catch(Exception e) {e.printStackTrace();}
	}
	
	
	public void output_binary_pheno_file(String out_phenotype_file, double[][]pheno_matrix) {
		try {
			BufferedWriter bw_b = new BufferedWriter(new FileWriter(out_phenotype_file+".all.pheno"));
			bw_b.write("ID");
			for(int phe_i=0;phe_i<this.num_pheno_rounds;phe_i++) {
				bw_b.write("\tPhenotype_round_"+phe_i);
			}bw_b.write("\n");
			for(int subj_i=0;subj_i<this.num_subject;subj_i++) {
				bw_b.write(this.genotype.sample_ids[subj_i]);
				for(int phe_i=0;phe_i<this.num_pheno_rounds;phe_i++) {
					int disease_state=0;
					if(Double.compare(pheno_matrix[phe_i][subj_i], this.liability_threshold)>=0) {
						disease_state=1;
					}
					bw_b.write("\t"+disease_state);
				}bw_b.write("\n");
			}
			bw_b.close();
		}catch(Exception e) {e.printStackTrace();}
	}
	
	
	public static ArrayList<String[]> return_Vars(String used_var){
		ArrayList<String[]> causal_varsArrayList = new ArrayList<String[]>();//all SNPs(chr,pos) for one simualtion are in one array
		try {
			BufferedReader used_varReader = new BufferedReader(new FileReader(used_var));
			String var_lineString = used_varReader.readLine();
			while(var_lineString != null) {
				if(!var_lineString.startsWith("#")) {
					String[] var_lineStrings = var_lineString.split(" ");
					String[] vars_in_one_compound =  new String[var_lineStrings.length];
					for(int i=0; i<var_lineStrings.length;i++) {
						String[] temp_lineStrings = var_lineStrings[i].split(":");
						String causal_var = temp_lineStrings[1]+","+temp_lineStrings[5];
						vars_in_one_compound[i]=causal_var;
					}
					causal_varsArrayList.add(vars_in_one_compound);
				}
				var_lineString=used_varReader.readLine();
			}
			used_varReader.close();
		}catch (Exception e) {
			// TODO: handle exception
			e.printStackTrace();
		}
		return(causal_varsArrayList);
	}

	public static boolean local_is_real_causal(String[] causal_SNPs_in_one_compounds, String causal_locus_start, int w) {
		boolean is_real_locus = false;
		for(int j=0;j<causal_SNPs_in_one_compounds.length;j++) {
			String[] SNPj = causal_SNPs_in_one_compounds[j].split(",");
			String[] SNP_start = causal_locus_start.split(",");
			if(SNPj[0].equals(SNP_start[0])) {
				/*
				 * becomes true only when the causal SNP from arraylist causal_SNPs is in the [] region 
				 */
				if((Integer.parseInt(SNPj[1].trim())>=Integer.parseInt(SNP_start[1].trim()))&&(Integer.parseInt(SNPj[1].trim())<=(Integer.parseInt(SNP_start[1].trim())+w))) {
					is_real_locus=true;
					break;
				}
			}
		}

		return(is_real_locus);
	}

	
	public static void local_power(ArrayList<String> pfiles_list, double adjusted_pvalue, String model, String outputfolder, int creteria,String used_var, int w) {
		int pfiles_len = pfiles_list.size();
		if(pfiles_len%2==0) {
			try {
				int half_pfiles_len=pfiles_len/2;
				int successful_count = 0;
				//read pvalues from each pvalue files (note, local and emmax has p1 and p2)
				//both local regions/SNP in p1 part and p2 part have significant pvalue lead to one successful finding 
				BufferedReader[] part1_pvalue_files = new BufferedReader[half_pfiles_len];
				BufferedReader[] part2_pvalue_files = new BufferedReader[half_pfiles_len];
				BufferedWriter power_fileBW = new BufferedWriter(new FileWriter(outputfolder+model+".creteria."+creteria+".power.log"));
				BufferedWriter p1_0_index = new BufferedWriter(new FileWriter(outputfolder+model+".part1.0.index.log"));
				BufferedWriter p2_0_index = new BufferedWriter(new FileWriter(outputfolder+model+".part2.0.index.log"));
				power_fileBW.write("#\tp1_true\tp1_total\tp2_true\tp2_total\n");
				//read all causal SNPs into an arraylist
				ArrayList<String[]> causal_varsArrayList = return_Vars(used_var);
				if(causal_varsArrayList.size()!=half_pfiles_len) {
					System.out.println("Not all causal compounds have been tested");
					System.exit(0);
				}
				//note pairewised pvalue files, line1: part 1; line2: part 2
				for(int i=0; i<half_pfiles_len;i++) {
					int part1_total_succe_count=0;
					int part1_real_succe_count=0;
					int part2_total_succe_count=0;
					int part2_real_succe_count=0;
					part1_pvalue_files[i] = new BufferedReader(new FileReader(pfiles_list.get(i))); //first half is the p value files name of l1
					part2_pvalue_files[i] = new BufferedReader(new FileReader(pfiles_list.get(i+half_pfiles_len)));//second half is the pvalue files for l2
					String part1_pvalue_lineString = part1_pvalue_files[i].readLine();
					String part2_pvalue_lineString = part2_pvalue_files[i].readLine();
					
					while(part1_pvalue_lineString != null) {
						if(!part1_pvalue_lineString.startsWith("#")&&!part1_pvalue_lineString.startsWith("chr")) {
							String[] part1_pvalue_lineStrings = part1_pvalue_lineString.split(",");
							double pvalue_of_this_line = Double.parseDouble(part1_pvalue_lineStrings[2]);
							String causal_locus_start = part1_pvalue_lineStrings[0]+","+part1_pvalue_lineStrings[1];
							// if the local/SNP pvalue is less than adjusted one
							if(Double.compare(pvalue_of_this_line, adjusted_pvalue)<=0) {
								part1_total_succe_count+=1;
								//if this region contains causal SNPs, real_succe_count+=1
								if(local_is_real_causal(causal_varsArrayList.get(i), causal_locus_start,w))part1_real_succe_count+=1;
							}
						}
						part1_pvalue_lineString=part1_pvalue_files[i].readLine();
					}
					
					while(part2_pvalue_lineString != null) {
						if(!part2_pvalue_lineString.startsWith("#")&&!part2_pvalue_lineString.startsWith("chr")) {
							String[] part2_pvalue_lineStrings = part2_pvalue_lineString.split(",");
							double pvalue_of_this_line = Double.parseDouble(part2_pvalue_lineStrings[2]);
							String causal_locus_start = part2_pvalue_lineStrings[0]+","+part2_pvalue_lineStrings[1];
							// if the local/SNP pvalue is less than adjusted one
							if(Double.compare(pvalue_of_this_line, adjusted_pvalue)<=0) {
								part2_total_succe_count+=1;
								//if this region contains causal SNPs, real_succe_count+=1
								if(local_is_real_causal(causal_varsArrayList.get(i), causal_locus_start,w))part2_real_succe_count+=1;
							}
						}
						part2_pvalue_lineString=part2_pvalue_files[i].readLine();
					}
					
					if(part1_real_succe_count==0) p1_0_index.write(i+",");
					if(part2_real_succe_count==0) p2_0_index.write(i+",");
					
					power_fileBW.write(i+"\t"+part1_real_succe_count+"\t"+part1_total_succe_count+"\t"+part2_real_succe_count+"\t"+part2_total_succe_count+"\n"); 
					if(creteria==1) {
						if((part1_real_succe_count>=1)||(part2_real_succe_count>=1)) successful_count+=1;
					}else if(creteria==2) {
						if((part1_real_succe_count>=1)&&(part2_real_succe_count>=1)) successful_count+=1;
					}
					part1_pvalue_files[i].close();
					part2_pvalue_files[i].close();
				}
				//output power
				double power = Double.valueOf(successful_count) / half_pfiles_len;
				power_fileBW.write("The power of "+half_pfiles_len+" pvalue files for model "+model+" is: "+power);
				power_fileBW.close();
				p1_0_index.close();
				p2_0_index.close();
			}catch (Exception e) {
				e.printStackTrace();
			}
		}else {
			System.out.println("Not pairwised local/emmax pfiles, please check again");
		}
	}
	
	public static void emmax_power(ArrayList<String> pfiles_list, double adjusted_pvalue, String model, String outputfolder, int creteria, String used_var) {
		int pfiles_len = pfiles_list.size();
		if(pfiles_len%2==0) {
			try {
				int successful_count = 0;
				int half_pfiles_len=pfiles_len/2;
				//read pvalues from each pvalue files (note, local and emmax has p1 and p2)
				//both local regions/SNP in p1 part and p2 part have significant pvalue lead to one successful finding 
				BufferedReader[] part1_pvalue_files = new BufferedReader[half_pfiles_len];
				BufferedReader[] part2_pvalue_files = new BufferedReader[half_pfiles_len];
				BufferedWriter power_fileBW = new BufferedWriter(new FileWriter(outputfolder+model+".creteria."+creteria+".power.log"));
				power_fileBW.write("#\tp1_true\tp1_total\tp2_true\tp2_total\n");
				//read all causal var into an arraylist 
				ArrayList<String[]> causal_varsArrayList = return_Vars(used_var);
				if(causal_varsArrayList.size()!=half_pfiles_len) {
					System.out.println("Not all causal compounds have been tested");
					System.exit(0);
				}
				//note pairewised pvalue files, first part are pfiles from part 1 and the second part are pfiles from part 2
				//creteria = 1, find one causal SNP means regional successful
				//creteria = 2. find all causal SNP means regional successful
				System.out.print("Begin pair-wised pvalue check");
				for(int i=0; i<half_pfiles_len;i++) {
					int part1_total_successful_count=0;
					int part1_real_successful_count=0;
					int part2_total_successful_count=0;
					int part2_real_successful_count=0;
					part1_pvalue_files[i] = new BufferedReader(new FileReader(pfiles_list.get(i))); //first half is the p value files name of l1
					part2_pvalue_files[i] = new BufferedReader(new FileReader(pfiles_list.get(i+half_pfiles_len)));//second half is the pvalue files for l2
					String part1_pvalue_lineString = part1_pvalue_files[i].readLine();
					String part2_pvalue_lineString = part2_pvalue_files[i].readLine();
					
					while(part1_pvalue_lineString != null) {
						if(!part1_pvalue_lineString.startsWith("#")&&!part1_pvalue_lineString.startsWith("chr")) {
							String[] part1_pvalue_lineStrings = part1_pvalue_lineString.split(",");
							double pvalue_of_this_line = Double.parseDouble(part1_pvalue_lineStrings[2]);
							String pvalue_of_this_line_postion = part1_pvalue_lineStrings[0]+","+part1_pvalue_lineStrings[1];
							// if the local/SNP pvalue is less than adjusted one
							if(Double.compare(pvalue_of_this_line, adjusted_pvalue)<=0) {
								part1_total_successful_count+=1;
								//check if this postion is the causal SNP
								String[] vars_in_one_compound = causal_varsArrayList.get(i);
								for(int k=0;k< vars_in_one_compound.length;k++) {
									if(pvalue_of_this_line_postion.equals(vars_in_one_compound[k])) {
										part1_real_successful_count++;
										break;
									}
								}
							}
						}
						part1_pvalue_lineString=part1_pvalue_files[i].readLine();
					}
					
					while(part2_pvalue_lineString != null) {
						if(!part2_pvalue_lineString.startsWith("#")&&!part2_pvalue_lineString.startsWith("chr")) {
							String[] part2_pvalue_lineStrings = part2_pvalue_lineString.split(",");
							double pvalue_of_this_line = Double.parseDouble(part2_pvalue_lineStrings[2]);
							String pvalue_of_this_line_postion = part2_pvalue_lineStrings[0]+","+part2_pvalue_lineStrings[1];
							// if the local/SNP pvalue is less than adjusted one
							if(Double.compare(pvalue_of_this_line, adjusted_pvalue)<=0) {
								part2_total_successful_count+=1;
								//check if this postion is the causal SNP 
								String[] vars_in_one_compound = causal_varsArrayList.get(i);
								for(int k=0;k< vars_in_one_compound.length;k++) {
									if(pvalue_of_this_line_postion.equals(vars_in_one_compound[k])) {
										part2_real_successful_count++;
										break;
									}
								}
							}
						}
						part2_pvalue_lineString=part2_pvalue_files[i].readLine();
					}
					
					power_fileBW.write(i+"\t"+part1_real_successful_count+"\t"+part1_total_successful_count+"\t"+part2_real_successful_count+"\t"+part2_total_successful_count+"\n");
					if(creteria==1) {
						if((part1_real_successful_count>=1)||(part2_real_successful_count>=1)) {
							successful_count+=1;
						}
					}else if(creteria==2) {
						if((part1_real_successful_count>=1)&&(part2_real_successful_count>=1)) {
							successful_count+=1;
						}
					}
					part1_pvalue_files[i].close();
					part2_pvalue_files[i].close();
				}
				//output power 
				double power = Double.valueOf(successful_count) / half_pfiles_len;
				power_fileBW.write("The power of "+half_pfiles_len+" pvalue files for model "+model+" is: "+power);
				power_fileBW.close();
			}catch (Exception e) {
				e.printStackTrace();
			}
		}else {
			System.out.println("Not pairwised local/emmax pfiles, please check again");
		}
	}
	

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		if (args.length == 0) {
			System.out.println("======================================================\n"
					+ "Developer: Qing Li Nov 2019\n"
					+ "Usage: java -Xmx4g -jar SimHiC.jar [function]\n"
					+ "Supported functions:\n"
					+ "\tsimulation\n"
					+ "\tpower\n");
			System.exit(0);
		} else {
			String function = args[0];
			int number_of_parameter=15;
			if(function.equals("simulation")) {
				if(args.length<number_of_parameter*2+1) {
					System.out.println("Usage: \n"
							+ "\t<-genotype\tinput_genotyp_file(HDF5-format)>\n"
							+ "\t<-ov\toutput_file_recording_causal_variants>\n"
							+ "\t<-oc\toutput_file_recording_causal_compounds>\n"
							+ "\t<-op\toutput_file_for_phenotypes>\n"
							+ "\t<-sr\tsimulation_report_file>\n"
							+ "\t<-info\tinput_candidate_genes_file>\n"
							+ "\t<-type\ttype_of_content_in_candidate_genes(num/SNP/compound)>\n"
							+ "\t<-nr\tnum_of_rounds_to_simulate>\n"
							+ "\t<-nc\tnum_of_causal_regions>\n"
							+ "\t<-ns\tnum_of_causal_variants_per_region>\n"
							+ "\t<-f\tflank_causal_length>\n"
							+ "\t<-min\tminimal_num_MAF_of_selected_variants>\n"
							+ "\t<-max\tmaximal_num_MAF_of_selected_variants>\n"
							+ "\t<-gc\tgenetic_component_(a_number_between_0_and_1)>\n"
							+ "\t<-p\tpattern_of_genetic_interaction_between_2_regions_of_1_compound(add/hetero/both/compensate)>\n"
							+ "\t[-th\tliability_threshold_above_which_it_is_a_case]\n"
							+ "\t[-binary\twhether_it_is_quantitative_or_binary]\n"
							+ "There are "+number_of_parameter+" mandatory paraneters. \nBut you have only specified "
									+(args.length-1)/2+" parameters");
					System.exit(0);
				}else {
					String genotype_hdf5=null;
					String output_var_file=null;
					String output_compounds_file=null;
					String output_pheno_file=null;
					String simulation_report_file=null;
					String candidate_genes_file=null;
					String type=null;
					int num_rounds=10;
					int num_causal_compounds=2;
					int num_causal_var_per_region=2;
					int flank_causal_length=0;
					double min_maf=0.005;
					double max_maf=0.05;
					double genetic_component=0.7;
					String pattern=null;
					boolean quantitative=true;
					double liability_threshold=Double.NaN;
					for(int k=1;k<args.length;k++) {
						if(args[k].startsWith("-")) {
							if(args[k].startsWith("-genotype"))genotype_hdf5=args[k+1];
							else if(args[k].equals("-ov"))output_var_file=args[k+1];
							else if(args[k].equals("-oc"))output_compounds_file=args[k+1];
							else if(args[k].equals("-op"))output_pheno_file=args[k+1];
							else if(args[k].equals("-sr"))simulation_report_file=args[k+1];
							else if(args[k].equals("-info"))candidate_genes_file=args[k+1];
							else if(args[k].equals("-type"))type=args[k+1];
							else if(args[k].equals("-nr"))num_rounds=Integer.parseInt(args[k+1]);
							else if(args[k].equals("-nc"))num_causal_compounds=Integer.parseInt(args[k+1]);
							else if(args[k].equals("-ns"))num_causal_var_per_region=Integer.parseInt(args[k+1]);
							else if(args[k].equals("-f"))flank_causal_length=Integer.parseInt(args[k+1]);
							else if(args[k].equals("-min"))min_maf=Double.parseDouble(args[k+1]);
							else if(args[k].equals("-max"))max_maf=Double.parseDouble(args[k+1]);
							else if(args[k].equals("-gc"))genetic_component=Double.parseDouble(args[k+1]);
							else if(args[k].equals("-p"))pattern=args[k+1];
							else if(args[k].equals("-th"))liability_threshold=Double.parseDouble(args[k+1]);
							else if(args[k].equals("-binary"))quantitative=false;
							else {
								System.out.println("Sorry, the option "+args[k]+"is not supported. Please check again");
								System.exit(0);
							}
						}
					}					
					SimHiC sim_hic = new SimHiC(genotype_hdf5, output_var_file, output_compounds_file, num_rounds, min_maf, max_maf, num_causal_compounds, 
							num_causal_var_per_region, genetic_component, liability_threshold, candidate_genes_file, type, flank_causal_length);
					double[][] pheno_matrix = sim_hic.sim_quantitative_phenotype(output_var_file, simulation_report_file, pattern,output_pheno_file);
					if(quantitative) {
						sim_hic.output_quanty_pheno_file(output_pheno_file, pheno_matrix);
					}else {
						sim_hic.output_binary_pheno_file(output_pheno_file, pheno_matrix);
					}
				}
			}else {
				System.out.println("Sorry, the function "+args[0]+" doesn't exist");
			}
		}
		
//		String genotype_hdf5="D:/HiCSim_test/NFBC.num.ch1-22.csv.b37.dup.hdf5" ;
//		String output_var_file="D:/HiCSim_test/output_var.txt";
//		String output_compounds_file="D:/HiCSim_test/output_compounds.txt";
//		String output_pheno_file="D:/HiCSim_test/output_pheno.tsv";
//		String simulation_report_file="D:/HiCSim_test/sim_report.txt";
//		String candidate_genes_file="D:/HiCSim_test/26499245_master_pns.compound.txtsorted.txt";
//		String type="compound";
//		int num_rounds=10;
//		int num_causal_compounds=1;
//		int num_causal_var_per_region=2;
//		int flank_causal_length=1000;
//		double min_maf=0.05;
//		double max_maf=0.4;
//		double genetic_component=0.9999;
//		String pattern="compensate";
//		boolean quantitative=true;
//		double liability_threshold=Double.NaN;
//		SimHiC sim_hic = new SimHiC(genotype_hdf5, output_var_file, output_compounds_file,num_rounds, min_maf, max_maf, num_causal_compounds, 
//		num_causal_var_per_region, genetic_component, liability_threshold, candidate_genes_file, type, flank_causal_length);
//		double[][] pheno_matrix = sim_hic.sim_quantitative_phenotype(output_var_file, simulation_report_file, pattern);
//		if(quantitative) {
//			sim_hic.output_quanty_pheno_file(output_pheno_file, pheno_matrix);
//		}else {
//			sim_hic.output_binary_pheno_file(output_pheno_file, pheno_matrix);
//		}
		
//		String pfiles_path_file ="";
//		String pvalue_adjusted_file = "";
//		String model = "";
//		String outputfolder = "";
//		SimHiC.power_caculation( pfiles_path_file, pvalue_adjusted_file,model,outputfolder);
		
	}
}
