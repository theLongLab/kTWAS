package simulations;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

import org.apache.commons.math3.distribution.NormalDistribution;

import mixedmodel.KinshipMatrix;
import mixedmodel.Mutation;
import mixedmodel.VariantsDouble;
import myMathLib.Test;

/*
 * Parameters of running the simulations:
 * 
 * (1) number of cases/controls, whether matched by IBS
 * (2) penetrance, prevalence
 * (3) score-cutoff and reliability of the score  
 * (4) number of mutations
 * (5) single causal or multiple causal in a single gene
 * (6) definition of "healthy mutations", e.g., AC cutoff
 * (7) how many top genes to take in the analysis
 * (8) parents as control 
 * (9) IBS cutoff between cases/controls
 * (10) candidate genes
 * (11) multiple setting at once to reuse some computations
 * 
 * In response to the review comments for the 2nd revision:
 * (12) multiple unaffected siblings
 * (13) multiple affected siblings
 * (14.1) altering penetrance by environmental factors
 * (14.2) altering penetrance by MAF
 * (15) output the simulation level variants and prioritization scores, for use when cross-checking with other power methods
 * (16) output the score. 
 */

public class SimPEL {

	// input data 
	VariantsDouble genotype; 
	//KinshipMatrix kinship; // As IBD check has a scaling problem, it is difficult for the users to specify them correctly.
	final int ExAC_sample_size=2*60706;
	HashMap<String, String> id2ethnicity= new HashMap<>();
	HashMap<String, ArrayList<String>> ethnicity2ids= new HashMap<>();
	HashMap<String, String> single_parent= new HashMap<>();
	HashMap<String, String> both_parents= new HashMap<>();
	HashMap<String, String> single_sibling= new HashMap<>();
	HashMap<String, String> multi_siblings= new HashMap<>();
	HashMap<String, String> gender= new HashMap<>();
 	int[][] candidate_genes; // #genes x 3
 	ArrayList<String> candidate_gene_names= new ArrayList<>();
	int[][] all_genes; // #genes x 3
	ArrayList<Mutation>[] candidate_muts; // #genes x #muts
	Mutation[] mutations;
	double[][] annotation_score; // #chr x #(annotable locs)
	
	// supporting data structure during the analysis
	// boolean[] all_labels;
	double[] causal_real_scores; // size=this.heterogeneity
	int[] causal_index_in_candidate_genes;  // size=this.heterogeneity
	int[] causal_index_in_all_genes;	 // size=this.heterogeneity
	
	int[] sim_case_indexes;
	int[] sim_control_indexes;	 // for one-one case/control study
	int[] sim_p_control_indexes2;  // for trios: the sim_control_indexes above will be one parent, and sim_control_indexes2 will be the other	
	int[] sim_s_control_indexes2;  // for multi_siblings: the sim_control_indexes above will be one sibling, and sim_s_control_indexes2 will be the other	
	//int[] casual_genes;          // size=this.heterogeneity
	//int[] mutated_allgene_index_in_cases;   // size=this.case_num
	int parent_control;
	int sibling_control;
	int trios_control;
	
	int multi_sibling_control;
	int case_same_family;
	boolean multi_affacted_sibling=false;
	
	// analytic parameters: 
	double score_weight;
	double ann_score_cutoff;
	double ann_score_confidence;
	boolean genes_file_known;
	//boolean general_pop_freq_known;
	
	// properties of the design/genetic architecture
	int num_case;	// numbers of control and case are the same, as num_case
	//double case_control_IBS;
	double general_pop_freq;
	double penetrance;
	double penetrance_var;
	double penetrance_mean;
	double prevalence;
	double prevalence_mean;
	double prevalence_var;
	boolean compound_het; 
	int heterogeneity; // how many causal genes default 1.
	
	long start;
	
	public SimPEL(
			// input files to be parsed
			String gencode_file, String candidate_genes_file, String KG_ethnicity_family_file, String KG_hdf5, String ExAC_file, String annotation_file,
			// properties of the study design
			int num_case, int heterogeneity, double general_pop_freq, 
			double penetrance_mean, double penetrance_var, 
			double prevalence_mean, double prevalence_var, 
			int parent_control, int sibling_control, int trios_control, int multi_sibling_control, boolean multi_affacted_sibling,
			// analytic parameters  
			double score_weight, double ann_score_cutoff, double ann_score_confidence, boolean genes_file_known, boolean compound_het,
			double expected_causal_maf, double impact_expected_causal_maf, double impact_expected_gxe, 
			
			// tmp_folder
			String tmp_files_folder){


		this.start=System.currentTimeMillis();
		// analytic parameters	
		this.score_weight=score_weight;
		this.ann_score_cutoff=ann_score_cutoff;
		this.ann_score_confidence=ann_score_confidence;
		this.genes_file_known=genes_file_known;
		
		// properties of the design
		this.num_case=num_case;	
		this.heterogeneity=heterogeneity;
		this.general_pop_freq=general_pop_freq;
		this.penetrance_mean=penetrance_mean;
		this.penetrance_var=penetrance_var;
		this.prevalence_mean=prevalence_mean;
		if(impact_expected_causal_maf!=0) // the penetrance will be adjusted by the expected MAF. 
			this.penetrance_mean=this.penetrance_mean-(impact_expected_causal_maf)*(expected_causal_maf);
		if(impact_expected_gxe!=0)  	  // positive impact_expected_gxe will make the penetrance higher
			this.penetrance_mean=this.penetrance_mean+impact_expected_gxe;
		this.prevalence_var=prevalence_var;
		this.compound_het=compound_het;
		this.sibling_control=sibling_control;
		this.parent_control=parent_control;
		this.trios_control=trios_control;	
		this.multi_sibling_control=multi_sibling_control;
		this.multi_affacted_sibling=multi_affacted_sibling;	
		
		// parsing input files
		this.genotype=new VariantsDouble(KG_hdf5);  // 1000G genotype file
		generate_relationship_tables(KG_ethnicity_family_file);
		//this.all_labels=new boolean[this.genotype.sample_size];
		this.mutations=new Mutation[this.genotype.sample_size];
		HashMap<String, Integer> chr2index=myFileFunctions.FileFunc.generate_chr_map();
		try{
			// read in candidate genes file			
			BufferedReader br_g=new BufferedReader(new FileReader(candidate_genes_file)); // candidate genes
			String line_g = br_g.readLine();
			ArrayList<int[]> candidate_genes_array= new ArrayList<>();
			while(line_g!=null){
				String[] tmp=line_g.split("\t");
				int[] this_region=new int[3]; // chr, start, end
				this_region[0]=chr2index.get(tmp[1]);
				this_region[1]=Integer.parseInt(tmp[2]);
				this_region[2]=Integer.parseInt(tmp[3]);
				candidate_genes_array.add(this_region);
				this.candidate_gene_names.add(tmp[0]);
				line_g=br_g.readLine();
			}br_g.close();
			this.candidate_genes=myFileFunctions.FileFunc.arraylist2intarray2(candidate_genes_array);
			System.out.println("Causal genes pool loaded: "+this.candidate_genes.length+" genes");
			if(this.candidate_genes.length<this.heterogeneity){
				System.out.println("The causal genes pool is smaller than number of causal genes! Exit");
				System.exit(0);
			}
			// read in all genes file
			br_g=new BufferedReader(new FileReader(gencode_file));
		    line_g=br_g.readLine();
		    ArrayList<int[]> all_genes_array= new ArrayList<>();
			while(line_g!=null){
				String[] tmp=line_g.split("\t");
				int[] this_region=new int[3]; // chr, start, end
				this_region[0]=chr2index.get(tmp[0]);
				this_region[1]=Integer.parseInt(tmp[3]);
				this_region[2]=Integer.parseInt(tmp[4]);
				all_genes_array.add(this_region);
				line_g=br_g.readLine();
			}
			this.all_genes=myFileFunctions.FileFunc.arraylist2intarray2(all_genes_array);	
			System.out.println("All genes loaded: "+this.all_genes.length + " genes");
			// generate candidate mutations list
			this.candidate_muts=new ArrayList[this.candidate_genes.length];
			for(int k=0;k<this.candidate_genes.length;k++)
				this.candidate_muts[k]=new ArrayList<Mutation>();
			this.generate_mut_table(ExAC_file, annotation_file, (int)(general_pop_freq*this.ExAC_sample_size));
				
			// score all locs in the pool population based on the annotation tools
			create_tmp_files(tmp_files_folder, annotation_file);
			this.annotation_score=new double[genotype.num_chrs][];
			for(int chr=0;chr<genotype.num_chrs;chr++){
				this.annotation_score[chr]=new double[genotype.locations[chr].length];
				int scored=0;
				String the_chr_mcap=tmp_files_folder+"/"+(chr+1)+".anno.txt";
				BufferedReader br_mcap=new BufferedReader(new FileReader(the_chr_mcap));
				String line_mcap=br_mcap.readLine();
				int loc_index=0;
				int loc_geno=genotype.locations[chr][loc_index];
				String[] tmp=line_mcap.split("\t");
				int loc_mcap=Integer.parseInt(tmp[1]);
				while(line_mcap!=null){
					loc_mcap=Integer.parseInt(tmp[1]);
					if(loc_geno<loc_mcap){
						while(loc_geno<loc_mcap && loc_index!=genotype.locations[chr].length-1){
							loc_index++;
							loc_geno=genotype.locations[chr][loc_index];
						}
					}
					if(loc_geno>=loc_mcap){ 
						if(loc_geno==loc_mcap){ // line matches
							this.annotation_score[chr][loc_index]=Double.parseDouble(tmp[4]);
							loc_index++;
							scored++;
						}//else: loc_mcap smaller, do nothing and read the next line. 
					}else{ // loc_index!=genotype.locations[chr].length-1; the last one!
						break; 
					}
					line_mcap=br_mcap.readLine();
					if(line_mcap!=null)tmp=line_mcap.split("\t");
					while(line_mcap!=null && Integer.parseInt(tmp[1])==loc_mcap){
						line_mcap=br_mcap.readLine();
						if(line_mcap!=null)tmp=line_mcap.split("\t");
					}
				}br_mcap.close();
				//annotation_score[chr]=new double[genotype.num_sites[chr]];
				System.out.println("Preparation: Chr"+(chr+1)+" scored: "+genotype.num_sites[chr]+"/"+scored);
			}			
		}catch(Exception e){e.printStackTrace();}
	}
	
	/*
	 * separate the annotation file into multiple chromosomes
	 */
	public static void create_tmp_files(String tmp_files_folder, String annotation_file){
		final int total_num_chr=22;
		HashMap<String, Integer> chr_map=myFileFunctions.FileFunc.generate_chr_map();
		try{
			boolean everyone_there=false;
			if((new File(tmp_files_folder)).exists()){
				everyone_there=true;
				for(int chr=0;chr<total_num_chr;chr++){
					if(!(new File(tmp_files_folder+"/"+(chr+1)+".anno.txt")).exists()){
						everyone_there=false;break;
					}
				}
			}if(!everyone_there){
				System.out.println("Generating temporary annotation files.");
				(new File(tmp_files_folder)).mkdirs();
				BufferedReader br=new BufferedReader(new FileReader(annotation_file));
				HashMap<Integer, BufferedWriter> bws= new HashMap<>();
				String line=br.readLine();
				while(line!=null){
					if(!line.startsWith("#")){
						String[] tmp=line.split("\t");
						int chr=chr_map.get(tmp[0]);
						if(bws.containsKey(chr)){
							bws.get(chr).write(line+"\n");
						}else{
							BufferedWriter bw=new BufferedWriter(new FileWriter(tmp_files_folder+"/"+(chr)+".anno.txt"));
							bw.write(line+"\n");
							bws.put(chr, bw);
						}
					}
					line=br.readLine();
				}for(int chr: bws.keySet()){
					bws.get(chr).close();
				}
				System.out.println("Temporary annotation files are ready.");
			}else{
				System.out.println("Temporary annotation files exist.");
			}
		}catch(Exception e){e.printStackTrace();}
	}
	
	/*
	 * generate gender, parents, sibling, trio, ethic groups tables.
	 */
	
	public void generate_relationship_tables(String KG_ethnicity_family_file){
		HashMap<String, String> gender_coding= new HashMap<>();
		gender_coding.put("1", "/Male"); gender_coding.put("2", "/Female");
		try{
			BufferedReader br_e=new BufferedReader(new FileReader(KG_ethnicity_family_file)); // 
			String line_e=br_e.readLine();line_e=br_e.readLine();
			HashMap<String, Integer> all_ids2index=this.genotype.sample_id2index;
			while(line_e!=null){
				String[] tmp=line_e.split("\t");
				if(all_ids2index.containsKey(tmp[1])){
					this.id2ethnicity.put(tmp[1], tmp[6]+gender_coding.get(tmp[4]));
					if(this.ethnicity2ids.containsKey(tmp[6]+gender_coding.get(tmp[4]))){
						ArrayList<String> current_group=this.ethnicity2ids.get(tmp[6]+gender_coding.get(tmp[4]));
						current_group.add(tmp[1]);
					}else{
						ArrayList<String> current_group= new ArrayList<>();
						current_group.add(tmp[1]);
						this.ethnicity2ids.put(tmp[6]+gender_coding.get(tmp[4]), current_group);
					}
				}				
				this.gender.put(tmp[1], tmp[4]);
				line_e=br_e.readLine();
			}br_e.close();
			// sibling and parent
			br_e=new BufferedReader(new FileReader(KG_ethnicity_family_file)); // 
			line_e=br_e.readLine();line_e=br_e.readLine();
			while(line_e!=null){
				String[] tmp=line_e.split("\t");
				if(!tmp[2].equals("0") && !tmp[3].equals("0") && all_ids2index.containsKey(tmp[1]) && 
						all_ids2index.containsKey(tmp[2]) && all_ids2index.containsKey(tmp[3])){
					this.both_parents.put(tmp[1], tmp[2]+";"+tmp[3]);
				}
				if(this.gender.get(tmp[1]).equals("1")){ // "1"=male: so select male parent and sibling if possible
					if((!tmp[2].equals("0")) && (!tmp[1].equals(tmp[2]))
							&& all_ids2index.containsKey(tmp[1]) && all_ids2index.containsKey(tmp[2])){ // there is a father
						this.single_parent.put(tmp[1], tmp[2]);
					}
					if(!tmp[8].equals("0")){ //there is sibling(s)
						String[] sibs =tmp[8].split(", ");
						String multi_sibs=""; int multi_sibs_found=0;
						for(int k=0;k<sibs.length;k++){
							if(this.gender.containsKey(sibs[k]) && (!tmp[1].equals(sibs[k])) && this.gender.get(sibs[k]).equals("1")
									&& all_ids2index.containsKey(tmp[1]) && all_ids2index.containsKey(sibs[k])){
								this.single_sibling.put(tmp[1], sibs[k]);  // found a sibling, with gender matched.
								if(multi_sibs_found==0){ //  found a sibling, with gender matched, contribute to multiple_sibs
									multi_sibs=sibs[k];
									multi_sibs_found=1;
								}
								else if(multi_sibs_found==1){ // already found 1, plus the current one, that will be 2.
									multi_sibs=multi_sibs+";"+sibs[k];
									this.multi_siblings.put(tmp[1], multi_sibs);
									break; 
								}
							}
						}
					}
				}else if(this.gender.get(tmp[1]).equals("2")){ // "2"=female: so select female parent and sibling if possible
					if((!tmp[3].equals("0")) && (!tmp[3].equals(tmp[1]))
							&& all_ids2index.containsKey(tmp[1]) && all_ids2index.containsKey(tmp[3])){ // there is a mother
						this.single_parent.put(tmp[1], tmp[3]);
					}
					if(!tmp[8].equals("0")){ //there is sibling(s)
						String[] sibs =tmp[8].split(", ");
						String multi_sibs=""; int multi_sibs_found=0;
						for(int k=0;k<sibs.length;k++){
							if(this.gender.containsKey(sibs[k]) && (!tmp[1].equals(sibs[k])) && this.gender.get(sibs[k]).equals("2")
									&& all_ids2index.containsKey(tmp[1]) && all_ids2index.containsKey(sibs[k])){
								this.single_sibling.put(tmp[1], sibs[k]);  // found a sibling, with gender matched.
								if(multi_sibs_found==0){ //  found a sibling, with gender matched, contribute to multiple_sibs
									multi_sibs=sibs[k];
									multi_sibs_found=1;
								}
								else if(multi_sibs_found>=1){ // already found some, plus the current one.
									multi_sibs=multi_sibs+";"+sibs[k];
									this.multi_siblings.put(tmp[1], multi_sibs);
									multi_sibs_found++; 
								}
							}
						}
					}
				}
				line_e=br_e.readLine();
			}
		}catch(Exception e){e.printStackTrace();}
	}
	
	/*  
	 *  Generate the in-region candidate mutations excluding variants that 
	 *  are above the specified allele count, i.e., common in healthy population 	   
	 */
	public void generate_mut_table(String ExAC_sites, String annotation_file, int AC_cutoff){
		HashMap<String, Integer> in_region_locs= new HashMap<>();
		for(int k=0;k<this.candidate_genes.length;k++){
			for(int loc=this.candidate_genes[k][1];loc<=this.candidate_genes[k][2];loc++){
				String the_loc=this.candidate_genes[k][0]+"_"+loc;
				in_region_locs.put(the_loc,k);
			}			
		}System.out.println("Loc_table formed: "+in_region_locs.size()+" locs");
		HashSet<String> avoid_muts= new HashSet<>();
		try{
			// record avoid list in the regions.
			//BufferedWriter bw=new BufferedWriter(new FileWriter(annotation_file+".inregion.unknown.txt"));
			BufferedReader br_a=new BufferedReader(new FileReader(ExAC_sites));
			String line_a=br_a.readLine();line_a=br_a.readLine(); // skip header
			while(line_a!=null){
				String[] tmp=line_a.split("\t");
				if(in_region_locs.containsKey(tmp[0]+"_"+tmp[1]) && Integer.parseInt(tmp[4])>=AC_cutoff){
					avoid_muts.add(tmp[0]+"_"+tmp[1]+"_"+tmp[2]+"_"+tmp[3]);
				}
				line_a=br_a.readLine();
			}br_a.close();
			System.out.println("Reference panel formed: "+avoid_muts.size()+" mutations");
			br_a=new BufferedReader(new FileReader(annotation_file));
			line_a=br_a.readLine();line_a=br_a.readLine(); // skip header
			int total_muts=0;
			while(line_a!=null){
				String[] tmp=line_a.split("\t");
				if(in_region_locs.containsKey(tmp[0]+"_"+tmp[1]) && // in region 
						Double.parseDouble(tmp[4])>=ann_score_cutoff){ 
					if(!avoid_muts.contains(tmp[0]+"_"+tmp[1]+"_"+tmp[2]+"_"+tmp[3])){ // but not in known muts
						//bw.write(line_a+"\n");
						int region_index=in_region_locs.get(tmp[0]+"_"+tmp[1]);
						Mutation the_mut=new Mutation(Integer.parseInt(tmp[0]),
								Integer.parseInt(tmp[1]),tmp[2],tmp[3], 
								Double.parseDouble(tmp[4])*(1.0+2.0*(Test.randomNumber()-1)*(1-ann_score_confidence))); // alter by confidence
						this.candidate_muts[region_index].add(the_mut);
						total_muts++;
					}
				}
				line_a=br_a.readLine();
			}
			System.out.println("Candidate mutation table generated: "+total_muts+" mutations.");
			//bw.close();
			br_a.close();
		}catch(Exception e){e.printStackTrace();}
		in_region_locs.clear();avoid_muts.clear();
	}
	
	/*
	 * Simulate case/control samples based on matching options;
	 * based on the field compond_het, we simulate one or two mutations 
	 */
	
	public void simulate_samples(){ // based on the field compond_het, we simulate one or two mutations 
		this.sim_case_indexes=new int[num_case]; 
		this.sim_control_indexes=new int[num_case];	// reset case and control labels
		this.sim_p_control_indexes2=new int[num_case];
		this.sim_s_control_indexes2=new int[num_case];
		Arrays.fill(this.sim_p_control_indexes2, -1); // sim_control_indexes2 is for trio only.
		Arrays.fill(this.sim_s_control_indexes2, -1); // sim_control_indexes2 is for multi_sib only.
		this.mutations=new Mutation[this.genotype.sample_size];  // reset all mutations
		this.causal_index_in_candidate_genes=new int[this.heterogeneity];// clear the index in previous simulation
		this.causal_index_in_all_genes=new int[this.heterogeneity]; // clear the index in previous simulation
		this.causal_real_scores=new double[this.heterogeneity]; // clear the score assigned last time! 
		// randomly select one or multiple (size=this.heterogeneity) causal genes in the candidate list.  
		HashSet<Integer> gene_indexes_set= new HashSet<>();
		for(int hg=0;hg<this.heterogeneity;hg++){
			causal_index_in_candidate_genes[hg]=(int)(Test.randomNumber()*this.candidate_genes.length);
			if(causal_index_in_candidate_genes[hg]==this.candidate_genes.length)causal_index_in_candidate_genes[hg]--;
			while(this.candidate_muts[causal_index_in_candidate_genes[hg]].size()<2||   // at least two candidate mutations
					gene_indexes_set.contains(causal_index_in_candidate_genes[hg])){ // not yet selected  
				causal_index_in_candidate_genes[hg]=(int)(Test.randomNumber()*this.candidate_genes.length);
				if(causal_index_in_candidate_genes[hg]==this.candidate_genes.length)causal_index_in_candidate_genes[hg]--;
			}
			System.out.println("Candidate Gene_"+causal_index_in_candidate_genes[hg]+" is the #"+hg+" causal gene.");
			gene_indexes_set.add(causal_index_in_candidate_genes[hg]);
		}		
		// generate cases/controls:
		HashSet<Integer> used_indexes= new HashSet<>();
		// generate samples that are multi_sibling controls
		Object[] indi_w_msibs=this.multi_siblings.keySet().toArray();
		int total_msibs=indi_w_msibs.length;
		if(total_msibs<this.multi_sibling_control){
			System.out.println("No sufficient multi-sibling control families in the population: please change to a larger pool.");
		}else{				
			for(int i=this.parent_control+this.sibling_control+this.trios_control;
				i<this.parent_control+this.sibling_control+this.trios_control+this.multi_sibling_control;i++){
				int tried=0;
				int index=(int)(Test.randomNumber()*total_msibs);
				String case_id=(String)indi_w_msibs[index];
				int case_index_in_all_samples=this.genotype.sample_id2index.get(case_id);
				String[] control_ids=this.multi_siblings.get(indi_w_msibs[index]).split(";");// two siblings
				int control_index_in_all_samples1=this.genotype.sample_id2index.get(control_ids[0]);
				int control_index_in_all_samples2=this.genotype.sample_id2index.get(control_ids[1]);
				while(index==total_msibs || used_indexes.contains(case_index_in_all_samples) 
						|| used_indexes.contains(control_index_in_all_samples1)
						|| used_indexes.contains(control_index_in_all_samples2)
						|| (control_ids.length<3 && this.multi_affacted_sibling && i==this.parent_control+this.sibling_control+this.trios_control)){
					index=(int)(Test.randomNumber()*total_msibs);
					case_id=(String)indi_w_msibs[index];
					case_index_in_all_samples=this.genotype.sample_id2index.get(case_id);
					control_ids=this.multi_siblings.get(indi_w_msibs[index]).split(";");  // could be more than 3, but take two first.
					control_index_in_all_samples1=this.genotype.sample_id2index.get(control_ids[0]);
					control_index_in_all_samples2=this.genotype.sample_id2index.get(control_ids[1]);
					tried++;
					if(tried>=total_msibs*2){
						System.out.println("No sufficient siblings in the population panel: "
								+ "\nplease change to a larger panel that contains sufficient number of large families.");
						System.exit(0);
					}
				}											
				this.sim_case_indexes[i]=case_index_in_all_samples;
				this.sim_control_indexes[i]=control_index_in_all_samples1;
				this.sim_s_control_indexes2[i]=control_index_in_all_samples2;
				used_indexes.add(this.sim_case_indexes[i]);
				used_indexes.add(this.sim_control_indexes[i]);
				used_indexes.add(this.sim_s_control_indexes2[i]);
				if(this.multi_affacted_sibling && i==this.parent_control+this.sibling_control+this.trios_control 
						&& this.multi_sibling_control>=2 && control_ids.length>=3){ // at least one sibling left
					i++;
					String affected_sib_id=control_ids[2];
					int case2_index_in_all_samples=this.genotype.sample_id2index.get(affected_sib_id);
					this.sim_case_indexes[i]=case2_index_in_all_samples;
					this.sim_control_indexes[i]=control_index_in_all_samples1;
					this.sim_s_control_indexes2[i]=control_index_in_all_samples2; // sharing two controls.
				}
			}	
		}	// generate samples that are trio controls
		Object[] indi_w_trios=this.both_parents.keySet().toArray();
		int total_trios=indi_w_trios.length;
		if(total_trios<this.trios_control){
			System.out.println("No sufficient trio families in the population: please change to a larger pool.");
		}else{				
			for(int i=this.parent_control+this.sibling_control;
					i<this.parent_control+this.sibling_control+this.trios_control;i++){
				int index=(int)(Test.randomNumber()*total_trios);
				String case_id=(String)indi_w_trios[index];
				int case_index_in_all_samples=this.genotype.sample_id2index.get(case_id);
				String[] control_ids=this.both_parents.get(indi_w_trios[index]).split(";");// two parents
				int control_index_in_all_samples1=this.genotype.sample_id2index.get(control_ids[0]);
				int control_index_in_all_samples2=this.genotype.sample_id2index.get(control_ids[1]);
				while(index==total_trios || used_indexes.contains(case_index_in_all_samples) 
						|| used_indexes.contains(control_index_in_all_samples1)
						|| used_indexes.contains(control_index_in_all_samples2)){
					index=(int)(Test.randomNumber()*total_trios);
					case_id=(String)indi_w_trios[index];
					case_index_in_all_samples=this.genotype.sample_id2index.get(case_id);
					control_ids=this.both_parents.get(indi_w_trios[index]).split(";");
					control_index_in_all_samples1=this.genotype.sample_id2index.get(control_ids[0]);
					control_index_in_all_samples2=this.genotype.sample_id2index.get(control_ids[1]);
				}											
				this.sim_case_indexes[i]=case_index_in_all_samples;
				this.sim_control_indexes[i]=control_index_in_all_samples1;
				this.sim_p_control_indexes2[i]=control_index_in_all_samples2;
				used_indexes.add(this.sim_case_indexes[i]);
				used_indexes.add(this.sim_control_indexes[i]);
				used_indexes.add(this.sim_p_control_indexes2[i]);
			}	
		}		
		// generate samples that are parental controls
		Object[] indi_w_parents=this.single_parent.keySet().toArray();
		int total_parents=indi_w_parents.length;
		if(total_parents<this.parent_control){
			System.out.println("No sufficient parent-chid pairs in the population: please change to a larger pool.");
		}else{				
			for(int i=0;i<this.parent_control;i++){
				int index=(int)(Test.randomNumber()*total_parents);
				String case_id=(String)indi_w_parents[index];
				int case_index_in_all_samples=this.genotype.sample_id2index.get(case_id);
				String control_id=this.single_parent.get(indi_w_parents[index]);
				int control_index_in_all_samples=this.genotype.sample_id2index.get(control_id);
				while(index==total_parents || used_indexes.contains(case_index_in_all_samples)
						|| used_indexes.contains(control_index_in_all_samples)){
					index=(int)(Test.randomNumber()*total_parents);
					case_id=(String)indi_w_parents[index];
					case_index_in_all_samples=this.genotype.sample_id2index.get(case_id);
					control_id=this.single_parent.get(indi_w_parents[index]);
					control_index_in_all_samples=this.genotype.sample_id2index.get(control_id);
				}		
				this.sim_case_indexes[i]=case_index_in_all_samples;
				this.sim_control_indexes[i]=control_index_in_all_samples;
				used_indexes.add(this.sim_case_indexes[i]);
				used_indexes.add(this.sim_control_indexes[i]);
			}	
		}	
		// generate samples that are sibling controls
		Object[] indi_w_siblings=this.single_sibling.keySet().toArray();
		int total_siblings=indi_w_siblings.length;
		if(total_siblings<this.sibling_control*2){
			System.out.println("No sufficient siblings in the population: please change to a larger pool.");
		}else{				
			for(int i=this.parent_control;i<this.parent_control+this.sibling_control;i++){
				int index=(int)(Test.randomNumber()*total_siblings);
				String case_id=(String)indi_w_siblings[index];
				int case_index_in_all_samples=this.genotype.sample_id2index.get(case_id);
				String control_id=this.single_sibling.get(indi_w_siblings[index]);
				int control_index_in_all_samples=this.genotype.sample_id2index.get(control_id);
				while(index==total_siblings || used_indexes.contains(case_index_in_all_samples) 
						|| used_indexes.contains(control_index_in_all_samples)){
					index=(int)(Test.randomNumber()*total_siblings);
					case_id=(String)indi_w_siblings[index];
					case_index_in_all_samples=this.genotype.sample_id2index.get(case_id);
					control_id=this.single_sibling.get(indi_w_siblings[index]);
					control_index_in_all_samples=this.genotype.sample_id2index.get(control_id);
				}											
				this.sim_case_indexes[i]=case_index_in_all_samples;
				this.sim_control_indexes[i]=control_index_in_all_samples;
				used_indexes.add(this.sim_case_indexes[i]);
				used_indexes.add(this.sim_control_indexes[i]);
			}	
		}
		
		// generate the rest random samples: matching ethic groups 
		if(genotype.sample_size<this.num_case*2){
			System.out.println("No sufficient individuals in the population: please change to a larger pool.");
		}else{
			// initiate the number of indis in each ethic group.
			HashMap<String, Integer> ethnicity_num= new HashMap<>();
			for(String it:ethnicity2ids.keySet()){
				ethnicity_num.put(it, ethnicity2ids.get(it).size());
			}
			for(int i=this.parent_control+this.sibling_control+this.trios_control+this.multi_sibling_control;i<num_case;i++){
				//select a case
				int index=(int)(Test.randomNumber()*genotype.sample_size);
				String ethnic_group=this.id2ethnicity.get(genotype.sample_ids[index]);
				while(index==genotype.sample_size || used_indexes.contains(index) || ethnicity_num.get(ethnic_group)<2){
					index=(int)(Test.randomNumber()*genotype.sample_size);
					ethnic_group=this.id2ethnicity.get(genotype.sample_ids[index]);
				}
				//select a control:
				ArrayList<String> ids_in_this_group=ethnicity2ids.get(ethnic_group);
				int index_in_group=(int)(Test.randomNumber()*ids_in_this_group.size());
				int index2=genotype.sample_id2index.get(ids_in_this_group.get(index_in_group));
				while(index_in_group==ids_in_this_group.size() || index==index2 ||
						used_indexes.contains(index2) || in_family(index, index2)!=0){
					index_in_group=(int)(Test.randomNumber()*ids_in_this_group.size());
					index2=genotype.sample_id2index.get(ids_in_this_group.get(index_in_group));
				}		
				this.sim_case_indexes[i]=index;
				this.sim_control_indexes[i]=index2;
				used_indexes.add(this.sim_case_indexes[i]);
				used_indexes.add(this.sim_control_indexes[i]);
				int ori_num=ethnicity_num.get(ethnic_group);
				ethnicity_num.put(ethnic_group,ori_num-2);					
			}			
		}
	
		// Generate mutations: 
		// (1) Cases will have mutations for sure; 
		// (2) Controls have causal mutations with the probability of [Prevalence x (1- Penetrance)]/[Penetrance x (1- Prevalence)]
		// 
		final double control_random_cut=prevalence*(1-penetrance)/(penetrance*(1-prevalence));
		for(int case_index=0;case_index<num_case;case_index++){
			int indi_index=this.sim_case_indexes[case_index];
			//assign a causal gene for this case. 	the_gene_index_in_causal \in {0,1,...,heterogeneity-1}
			// 										the_gene_index \in {0,1,...,candidate_genes.length-1}
			int the_gene_index_in_causal=(int)(Test.randomNumber()*heterogeneity);
			if(the_gene_index_in_causal==heterogeneity)the_gene_index_in_causal--;
			int the_gene_index=this.causal_index_in_candidate_genes[the_gene_index_in_causal];
			// assign a mutation for this gene & case
			if(!this.compound_het){
				int mut_index=(int)(Test.randomNumber()*this.candidate_muts[the_gene_index].size());
				if(mut_index==this.candidate_muts[the_gene_index].size())mut_index--;
				this.mutations[indi_index]=this.candidate_muts[the_gene_index].get(mut_index);
			}else{ // compound heterozygotes
				int mut_index1=(int)(Test.randomNumber()*this.candidate_muts[the_gene_index].size());
				if(mut_index1==this.candidate_muts[the_gene_index].size())mut_index1--;
				int mut_index2=(int)(Test.randomNumber()*this.candidate_muts[the_gene_index].size());
				while(mut_index2==this.candidate_muts[the_gene_index].size()||mut_index1==mut_index2){
					mut_index2=(int)(Test.randomNumber()*this.candidate_muts[the_gene_index].size());
				}
				this.mutations[indi_index]=new Mutation(this.candidate_muts[the_gene_index].get(mut_index1), 
						this.candidate_muts[the_gene_index].get(mut_index2));
			}
			if(Test.randomNumber()>=control_random_cut){ // contribute to the scores if the control doesn't have that mutation
				// this is equivalent to assigning the control a mutation. 
				if(this.sim_p_control_indexes2[case_index]==-1 && this.sim_s_control_indexes2[case_index]==-1){
					causal_real_scores[the_gene_index_in_causal]+=1+this.mutations[this.sim_case_indexes[case_index]].score*this.score_weight;													
					if(this.mutations[this.sim_case_indexes[case_index]].compound_het){
						causal_real_scores[the_gene_index_in_causal]+=this.mutations[this.sim_case_indexes[case_index]].score2*this.score_weight;
					}
				}else { // multiple family members are present. so need to check the other parent too. 
					if(Test.randomNumber()>=control_random_cut){
						causal_real_scores[the_gene_index_in_causal]+=1+this.mutations[this.sim_case_indexes[case_index]].score*this.score_weight;													
						if(this.mutations[this.sim_case_indexes[case_index]].compound_het){
							causal_real_scores[the_gene_index_in_causal]+=this.mutations[this.sim_case_indexes[case_index]].score2*this.score_weight;
						}
					}
				}
			}	
		}		
	}
	
	/*
	 * return 0 if unrelated; 1 if sibling; 2 if parent
	 */
	public int in_family(int index, int index2){
		String id=this.genotype.sample_ids[index];
		String id2=this.genotype.sample_ids[index2];
		if((this.single_parent.containsKey(id)&&this.single_parent.get(id).contains(id2))
				||(this.single_parent.containsKey(id2)&&this.single_parent.get(id2).contains(id)))return 2;
		else if((this.single_sibling.containsKey(id)&&this.single_sibling.get(id).contains(id2))
				||(this.single_sibling.containsKey(id2)&&this.single_sibling.get(id2).contains(id)))return 1;
		else return 0;		
	}
	
	/*
	 * Analysis when candidate genes are unknown 
	 * The first half of the returned array stores the standings, and the second half stores the scores.
	 */
	public double[] analysis_candidate_genes_unknown(){
		double[] standings_and_scores=new double[heterogeneity];
		//int[] int_standings=new int[this.heterogeneity];
		Arrays.fill(standings_and_scores, Double.NaN);
		//figure out the indexes of simulated casual genes in the all-genes list. 
		for(int hg=0;hg<heterogeneity;hg++){
			for(int g_i=0;g_i<this.all_genes.length;g_i++){
				if(this.candidate_genes[causal_index_in_candidate_genes[hg]][0]==this.all_genes[g_i][0] &&
						this.candidate_genes[causal_index_in_candidate_genes[hg]][1]==this.all_genes[g_i][1] &&
						this.candidate_genes[causal_index_in_candidate_genes[hg]][2]==this.all_genes[g_i][2]){
					this.causal_index_in_all_genes[hg]=g_i;
					System.out.println("The AllGene_"+this.causal_index_in_all_genes[hg]+" is the #"+hg+" causal gene.");
				}
			}
		}		
		double[][] gene_scores=new double[this.all_genes.length][num_case];		
		double[][] second_gene_scores=null;
		if(this.compound_het)second_gene_scores=new double[this.all_genes.length][num_case];		
		try{
			for(int chr=0;chr<genotype.num_chrs;chr++){
				System.out.println("Annotating genes in chr"+(chr+1));
				int largest_loc=genotype.locations[chr][genotype.locations[chr].length-1]+1;
				int[] loc2gene_index=new int[largest_loc];
				Arrays.fill(loc2gene_index, -1);
				for(int g=0;g<this.all_genes.length;g++){
					if(chr+1==this.all_genes[g][0] && this.all_genes[g][2]<largest_loc){
						 for(int loc=this.all_genes[g][1];loc<=this.all_genes[g][2];loc++){
							 loc2gene_index[loc]=g;
						 }
					}
				}
				//System.out.println("Loc2Gene map formed for chr"+(chr+1));
				for(int i=0;i<genotype.locations[chr].length;i++){
					int loc=genotype.locations[chr][i];
					if(loc2gene_index[loc]!=-1){// found the gene the loc belongs!
						int curr_gene_index=loc2gene_index[loc];
						double[] geno_data=genotype.load_one_variant_by_index(chr, i);
						for(int k=0;k<num_case;k++){
							int case_index=this.sim_case_indexes[k];
							int control_index=this.sim_control_indexes[k];
							int control_index2=this.sim_p_control_indexes2[k];
							if(control_index2==-1)control_index2=this.sim_s_control_indexes2[k];
							if(control_index2==-1){ // normal control: no trios/multi-siblings
								if(geno_data[control_index]==0 && geno_data[case_index]>0){
									if(this.annotation_score[chr][i]>gene_scores[curr_gene_index][k]){ // bigger than the best score in this gene
										gene_scores[curr_gene_index][k]=this.annotation_score[chr][i];
									}else if(this.compound_het && this.annotation_score[chr][i]>second_gene_scores[curr_gene_index][k]){
										// smaller than the best, but bigger than the second best score in this gene
										second_gene_scores[curr_gene_index][k]=this.annotation_score[chr][i];
									}
								}
							}else{  // trios: two control parents
								if((geno_data[control_index]==0 && geno_data[control_index2]==0 && geno_data[case_index]>0)||
										(geno_data[control_index]==1 && geno_data[control_index2]==1 && geno_data[case_index]==2)){
									if(this.annotation_score[chr][i]>gene_scores[curr_gene_index][k]){ // bigger than the best score in this gene
										gene_scores[curr_gene_index][k]=this.annotation_score[chr][i];
									}else if(this.compound_het && this.annotation_score[chr][i]>second_gene_scores[curr_gene_index][k]){
										// smaller than the best, but bigger than the second best score in this gene
										second_gene_scores[curr_gene_index][k]=this.annotation_score[chr][i];
									}
								}
							}
							
						}
					}					
				}//System.out.println("Score changed: "+change_counter);
			}
			// combine and output			
			double[] the_gene_score=new double[this.all_genes.length];
			for(int g=0;g<this.all_genes.length;g++){
				for(int k=0;k<num_case;k++){
					if(gene_scores[g][k]>ann_score_cutoff){
						the_gene_score[g]+=gene_scores[g][k]*score_weight+1;
						if(this.compound_het)the_gene_score[g]+=second_gene_scores[g][k]*score_weight;
					}
				}
			}
			Arrays.sort(the_gene_score); // ascending order
			for(int hg=0;hg<heterogeneity;hg++){
				the_gene_score[causal_index_in_all_genes[hg]]=0; // the same gene will not be scored here
				System.out.println("The_Real_Gene #"+hg+":\t"+causal_real_scores[hg]);
				standings_and_scores[hg+heterogeneity]=causal_real_scores[hg];
				standings_and_scores[hg]=Arrays.binarySearch(the_gene_score, causal_real_scores[hg]); 
				if(Double.isNaN(standings_and_scores[hg])){
					standings_and_scores[hg]= -100;	continue;
				}
				if(standings_and_scores[hg]<0)standings_and_scores[hg]=this.all_genes.length+standings_and_scores[hg]; // not found, convert the index to the positive one
				else if(standings_and_scores[hg]>=0)standings_and_scores[hg]=this.all_genes.length-standings_and_scores[hg];
				standings_and_scores[hg]=(int)(standings_and_scores[hg]+0.001); // ensure the machine error corrected
				System.out.println("Standings of the #"+hg+" gene: "+(standings_and_scores[hg]+1)+"/"+this.all_genes.length);					
			}					
		}catch(Exception e){e.printStackTrace();}
		return standings_and_scores;  
	}
	
	/*
	 * Analysis when candidate genes are unknown 
	 * The first half of the returned array stores the standings, and the second half stores the scores.
	 */
	public double[] analysis_candidate_genes_known(){
		double[] standings_and_score=new double[2*heterogeneity];
		//int[] int_standings=new int[this.heterogeneity];
		Arrays.fill(standings_and_score, Double.NaN);
		double[][] gene_scores=new double[this.candidate_genes.length][num_case];	
		double[][] second_gene_scores=null;
		if(this.compound_het)second_gene_scores=new double[this.candidate_genes.length][num_case];
		try{
			for(int chr=0;chr<genotype.num_chrs;chr++){
				System.out.println("Annotating genes in chr"+(chr+1));
				int largest_loc=genotype.locations[chr][genotype.locations[chr].length-1]+1;
				int[] loc2gene_index=new int[largest_loc];
				Arrays.fill(loc2gene_index, -1);
				for(int g=0;g<this.candidate_genes.length;g++){
					if(chr+1==this.candidate_genes[g][0] && this.candidate_genes[g][2]<largest_loc){
						 for(int loc=this.candidate_genes[g][1];loc<=this.candidate_genes[g][2];loc++){
							 loc2gene_index[loc]=g;
						 }
					}
				}
				//System.out.println("Loc2Gene map formed for chr"+(chr+1));
				for(int i=0;i<genotype.locations[chr].length;i++){
					int loc=genotype.locations[chr][i];
					if(loc2gene_index[loc]!=-1){// found the gene the loc belongs!
						int curr_gene_index=loc2gene_index[loc];
						double[] geno_data=genotype.load_one_variant_by_index(chr, i);
						for(int k=0;k<num_case;k++){
							int case_index=this.sim_case_indexes[k];
							int control_index=this.sim_control_indexes[k];
							int control_index2=this.sim_p_control_indexes2[k];		
							if(control_index2==-1)control_index2=this.sim_s_control_indexes2[k];
							if(control_index2==-1){ // normal control: no trios/multi-siblings
								if(geno_data[control_index]==0 && geno_data[case_index]>0){
									if(this.annotation_score[chr][i]>gene_scores[curr_gene_index][k]){ // bigger than the best score in this gene
										gene_scores[curr_gene_index][k]=this.annotation_score[chr][i];
									}else if(this.compound_het && this.annotation_score[chr][i]>second_gene_scores[curr_gene_index][k]){
										// smaller than the best, but bigger than the second best score in this gene
										second_gene_scores[curr_gene_index][k]=this.annotation_score[chr][i];
									}
								}
							}else{  // trios: two control parents
								if((geno_data[control_index]==0 && geno_data[control_index2]==0 && geno_data[case_index]>0)||
										(geno_data[control_index]==1 && geno_data[control_index2]==1 && geno_data[case_index]==2)){
									if(this.annotation_score[chr][i]>gene_scores[curr_gene_index][k]){ // bigger than the best score in this gene
										gene_scores[curr_gene_index][k]=this.annotation_score[chr][i];
									}else if(this.compound_het && this.annotation_score[chr][i]>second_gene_scores[curr_gene_index][k]){
										// smaller than the best, but bigger than the second best score in this gene
										second_gene_scores[curr_gene_index][k]=this.annotation_score[chr][i];
									}
								}
							}
						}
					}
				}//System.out.println("Score changed: "+change_counter);
			}
			// combine and output			
			double[] the_gene_score=new double[this.candidate_genes.length];
			for(int g=0;g<this.candidate_genes.length;g++){
				for(int k=0;k<num_case;k++){
					if(gene_scores[g][k]>ann_score_cutoff){
						the_gene_score[g]+=gene_scores[g][k]*score_weight+1;
						if(this.compound_het)the_gene_score[g]+=second_gene_scores[g][k]*score_weight;
					}
				}
			}
			Arrays.sort(the_gene_score); // ascending order
			for(int hg=0;hg<heterogeneity;hg++){
				the_gene_score[causal_index_in_candidate_genes[hg]]=0;// the same gene will not be scored here
				System.out.println("The_Real_Gene #"+hg+":\t"+causal_real_scores[hg]);
				standings_and_score[hg+heterogeneity]=causal_real_scores[hg];
				standings_and_score[hg]=Arrays.binarySearch(the_gene_score, causal_real_scores[hg]); 
				if(Double.isNaN(standings_and_score[hg])){
					standings_and_score[hg]= -100;	continue;
				}
				if(standings_and_score[hg]<0)standings_and_score[hg]=this.candidate_genes.length+standings_and_score[hg]; // not found, convert the index to the positive one
				else if(standings_and_score[hg]>=0)standings_and_score[hg]=this.candidate_genes.length-standings_and_score[hg];
				standings_and_score[hg]=(int)(standings_and_score[hg]+0.001); // ensure the machine error corrected
				System.out.println("Standings of the #"+hg+" gene: "+(standings_and_score[hg]+1)+"/"+this.candidate_genes.length);	
				
			}		
		}catch(Exception e){e.printStackTrace();}
		return standings_and_score;  // ensure the machine error corrected
	}
	
	public void test_power(int round, int num_top_gene, String report_file){
		double[][] standings=new double[round][];
		double[] power_count=new double[this.heterogeneity]; // power of finding 1,2,3,...all genes
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(report_file));	
			bw.write("#####################################################\n");
			bw.write("## Simulation-based Power Evaluation for design of Rare Diseases sequencing study:\n");	
			bw.write("## Parameters in this simulation:\n");			
			bw.write("#\t-num_cases:\t"+ num_case + " cases + "+num_case+" controls"+"\n");
			bw.write("#\t-sim_rounds:\t\t"+round+"\n");
			bw.write("#\t-rank_success:\t\t"+num_top_gene+"\n");
			bw.write("#\t-num_causal_genes:\t"+heterogeneity+"\n");
			bw.write("#\t-prevalence_mean:\t"+prevalence_mean+"\n");
			bw.write("#\t-prevalence_var:\t"+prevalence_var+"\n");
			bw.write("#\t-penetrance_mean:\t"+penetrance_mean+"\n");
			bw.write("#\t-penetrance_var:\t"+penetrance_var+"\n");
			bw.write("#\t-min_score:\t"+ann_score_cutoff+"\n");
			bw.write("#\t-conf_score:\t"+ann_score_confidence+"\n");
			bw.write("#\t-weight_score:\t"+score_weight+"\n");
			bw.write("#\t-pop_freq:\t"+general_pop_freq+"\n");
			bw.write("#\t-single_parents:\t"+parent_control+"\n");
			bw.write("#\t-both_parents:\t"+trios_control+"\n");
			bw.write("#\t-siblings:\t"+sibling_control+"\n");
			bw.write("#\t-include_all_genes:\t"+!genes_file_known+"\n");
			bw.write("#\t-compound_het:\t"+compound_het+"\n");
			bw.write("#####################################################\n");
			bw.flush();
			NormalDistribution ndist_pene=new NormalDistribution(penetrance_mean, penetrance_var);
			NormalDistribution ndist_pre=new NormalDistribution(prevalence_mean, prevalence_var);
			
			for(int r=0;r<round;r++){
				System.out.println("Round_"+r);
				bw.write("====Round_"+r+"====:\n");
				this.penetrance=ndist_pene.sample();
				if(this.penetrance>1)this.penetrance=1;
				if(this.penetrance<0)this.penetrance=0;
				this.prevalence=ndist_pre.sample();
				if(this.prevalence>1)this.prevalence=1;
				if(this.prevalence<0)this.prevalence=0;
				this.simulate_samples();				
				if(this.genes_file_known)standings[r]=this.analysis_candidate_genes_known();
				else	standings[r]=this.analysis_candidate_genes_unknown();
				bw.write("#Mutations:\n");
				bw.write("#Causal genes: ");
				for(int hg=0;hg<heterogeneity;hg++)
					bw.write(this.candidate_gene_names.get(causal_index_in_candidate_genes[hg])+"\t");
				bw.write("\n#Samples: "+"\n");
				for(int k=0;k<num_case;k++){
					bw.write("#\tCase:\t"+genotype.sample_ids[this.sim_case_indexes[k]]+"/"
								+this.id2ethnicity.get(genotype.sample_ids[this.sim_case_indexes[k]])+"/"
								+this.mutations[this.sim_case_indexes[k]].output()+"\n");
					bw.write("#\tCtrl:\t"+genotype.sample_ids[this.sim_control_indexes[k]]+"/"
							+this.id2ethnicity.get(genotype.sample_ids[this.sim_control_indexes[k]])+"\n");
					if(this.sim_p_control_indexes2[k]!=-1)
						bw.write("#\tCtrl:\t"+genotype.sample_ids[this.sim_p_control_indexes2[k]]+"/"
							+this.id2ethnicity.get(genotype.sample_ids[this.sim_p_control_indexes2[k]])+"\n");
					if(this.sim_s_control_indexes2[k]!=-1)
						bw.write("#\tCtrl:\t"+genotype.sample_ids[this.sim_s_control_indexes2[k]]+"/"
							+this.id2ethnicity.get(genotype.sample_ids[this.sim_s_control_indexes2[k]])+"\n");
//					if(this.mutations[this.sim_control_indexes[k]]!=null)
//						bw.write(this.mutations[this.sim_control_indexes[k]].output()+"\n");
//					else
//						bw.write("NA\n");
				}
				
				bw.write("Prioritization_score_in_round_"+r+":");				
				for(int hg=0;hg<heterogeneity;hg++){
					bw.write("\t"+(standings[r][hg+heterogeneity]));					
				}
				bw.write("\n");					
				bw.write("Rank_in_round_"+r+":");
				int count_this_round=0;
				for(int hg=0;hg<heterogeneity;hg++){
					bw.write("\t"+((int)(standings[r][hg]+0.0001)+1));
					if(standings[r][hg]+1<=num_top_gene){
						count_this_round++;
					}
				}if(count_this_round!=0){
					for(int k=0;k<count_this_round;k++)power_count[k]++;
				}
				bw.write("\n");				
				bw.flush();
			}
			bw.write("The powers :\n");
			bw.write("Finding at least "+(1)+" gene:\t"+(power_count[0]/round)+"\n");	
			for(int hg=1;hg<heterogeneity-1;hg++){
				bw.write("Finding at least "+(hg+1)+" genes:\t"+(power_count[hg]/round)+"\n");				
			}bw.write("Finding all "+heterogeneity+" genes:\t"+(power_count[heterogeneity-1]/round)+"\n");	
			bw.write("Time used:"+(System.currentTimeMillis()-start)/1000+" seconds.\n");
			bw.close();			
		}catch(Exception e){e.printStackTrace();}
		
		//return power;
	}
	
	
	public static void main(String[] args) {
		String folder="/Users/quanlong/Documents/projects/rare_annotate/";
		String KG_hdf5=folder+"g1k_all.hdf5";
		String MCAP=folder+"mcap_v1_0.txt";
		String candidate_genes_file1=folder+"HLA_ErbB.txt";
		//String candidate_genes_file1=folder+"gene_list_trusight_one_final.txt";
		String ethnicity_file=folder+"integrated_call_samples_v2.20130502.ALL.ped";
		String ExAC_file=folder+"ExAC.r0.3.1.sites.AC.txt";
		String tmp_working=folder+"tmp_chrs/";

		String gencode_file=folder+"gencode.v12.genenames.gtf";
		String report_file=folder+"report.known.t2.h2.ms2.txt";
		double prevalence_mean=0.01;
		double prevalence_var=0.01;
		double penetrance_mean=1;
		double penetrance_var=0.1;
		double ann_score_cutoff=0.025;
		double ann_score_confidence=0.95;
		double individual_weight=1;
		double general_pop_freq=0.001;
		boolean genes_file_known=true;
		boolean compound_het=false;
		int parent_control=0;
		int sibling_control=0;
		int trios_control=0;
		int multi_sibling_control=2;
		boolean multi_affacted_sibling=false;
		
		double expected_causal_maf=0;
		double impact_expected_causal_maf=0;
		double impact_expected_gxe=0; 
		
		int num_case=4;		
		int heterogeneity=2;
		int round=100;
		int num_top_gene=10;
		long start=System.currentTimeMillis();
		
		SimPEL ap=new SimPEL(gencode_file, candidate_genes_file1, ethnicity_file, KG_hdf5, ExAC_file, 
				MCAP, num_case, heterogeneity, general_pop_freq, penetrance_mean, penetrance_var, 
				prevalence_mean, prevalence_var, 
				parent_control,  sibling_control, trios_control, multi_sibling_control, multi_affacted_sibling,
				individual_weight, ann_score_cutoff, 
				ann_score_confidence, genes_file_known, compound_het, 
				expected_causal_maf, impact_expected_causal_maf, impact_expected_gxe, 
				tmp_working);
		ap.test_power(round, num_top_gene, report_file);
		System.out.println("Time used:"+(System.currentTimeMillis()-start)/1000);
	}

}
