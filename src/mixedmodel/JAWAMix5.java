package mixedmodel;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;

import Region_Emmax.Region_EMMAX;
import nam.Cross_info;
import nam.Founder;
import nam.Imputation;
import nam.RILs;
import simulations.SimPEL;


public class JAWAMix5 {

	private static String mit = "Copyright (C) 2012, QUAN LONG and QINGRUN ZHANG" +
			"\n\nPermission is hereby granted, free of charge, to any person obtaining a copy " +
			"\nof this software and associated documentation files (the \"Software\"), " +
			"\nto deal in the Software without restriction, including without limitation the " +
			"\nrights to use, copy, modify, merge, publish, distribute, sublicense, " +
			"\nand/or sell copies of the Software, and to permit persons to whom the " +
			"\nSoftware is furnished to do so, subject to the following conditions:" +
			"\n\nThe above copyright notice and this permission notice shall be included " +
			"\nin all copies or substantial portions of the Software." +
			"\n\nTHE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS " +
			"\nOR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, " +
			"\nFITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE " +
			"\nAUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER " +
			"\nLIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING " +
			"\nFROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS " +
			"\nIN THE SOFTWARE.";

	private static String JAWAMix5Intro = "JAWAMix5: JAva implementation of Whole-genome Association studies " +
			"using Mixed model based on HDF5 (r1.0.0).\n" +
			"Developer: Quan LONG & Qingrun ZHANG\n" +
			"Usage: java -Xmx2g -jar jawamix5.jar [function]";

	private static String supported_functions = "\n\temmax" +
			"\n\tlm" +
			"\n\temmax_stepwise" +
			"\n\tlm_stepwise" +
			"\n\temmax_select" +
			"\n\tlm_select" +
			"\n\tlocal" +
			"\n\tcompound" +
			"\n\trare" +
			"\n\tnam"+
			"\n\timport" +
			"\n\tchar2num" +
			"\n\ttped2num" +
			"\n\theritability"+
			"\n\tkinship"+
			"\n\trelationship"+
			"\n\thdf52csv"+
			"\n\tsimpel";

	public static void main(String[] args) {
		// Intro header if no function is specified. Provides a basic usage example.
		if(args.length==0){
			System.out.println("==========================================\n"+mit+"\n==========================================\n");
			System.out.println(JAWAMix5Intro);
			System.out.println("\nSupported functions:" + supported_functions);
			System.exit(0);
		}
		String function=args[0];	//Gets the argument supplied (name of the function to be run)


		//#####################################################################################
		// SIMPEL Function Initialization
		//#####################################################################################
		if (function.equals("simpel")) {
			if (args.length == 1) {
				System.out.println("Simulation-based Power Estimation for sequencing studies of Low-prevalence conditions.");
				System.out.println("Usage: \n\t" +
						"<-population_genotypes\tinput_genotype_file(HDF5-format)>\n\t" +
						"<-out\toutput_file>\n\t" +
						"<-causal_gene_pool\tcausal_genes_file>\n\t" +
						"<-population_pedigree\tsample_information_file>\n\t" +
						"<-all_genes\tgene_model_file_containing_all_genes>\n\t" +
						"<-mafs\tMAFs_in_healthy_populations (e.g.,ExAC)>\n\t" +
						"<-pathogenicity\tpathogenicity_score_file>\n\t" +
						"<-tmp_folder\ttemp_folder_for_working_files>\n\t" +
						"[-num_cases\t#patients (df=4)]\n\t" +
						"[-sim_rounds\t#simulation_rounds (df=100)]\n\t" +
						"[-rank_success\t#top_rank_genes_to_claim_success (df=10)]\n\t" +
						"[-num_causal_genes\tnumber_of_causal_genes (df=1)]\n\t" +
						"[-prevalence_mean\tpopulation_prevalence (df=0.01)]\n\t" +
						"[-penetrance_mean\tpenetrance_of_causal_allele (df=0.9)]\n\t" +
						"[-prevalence_var\tpopulation_prevalence (df=0.0)]\n\t" +
						"[-penetrance_var\tpenetrance_of_causal_allele (df=0.0)]\n\t" +
						"[-min_score\tpathogenicity_score_cutoff (df=0.025)]\n\t" +
						"[-conf_score\tconfidence_on_pathogenicity_score (df=0.95)]\n\t" +
						"[-weight_score\tthe_weight_multiplied_to_pathogenicity_score (df=1.0)]\n\t" +
						"[-pop_freq\tpopulation_frequency_cutoff (df=0.001)]\n\t" +
						"[-single_parents\t#controls_that_are_parents (df=0)]\n\t" +
						"[-both_parents\t#controls_that_are_trio_families (df=0)]\n\t" +
						"[-siblings\t#controls_that_are_siblings (df=0)]\n\t" +
						"[-multi_siblings\t#controls_that_are_multiple_siblings_in_the_same_family (df=0)]\n\t" +
						"[-multi_affacted_sibling\twhether_the_affacted_individuals_could_be_in_the_same_family]\n\t" +
						"[-include_all_genes\tinclude_all_genes_in_the_analysis]\n\t" +
						"[-compound_het\tconsider_compound_heterozygotes_only]\n\t" +
						"[-expected_causal_maf\texpected_MAF_of_the_causal_variants_in_healthy_population (df=0)]\n\t" +
						"[-impact_expected_causal_maf\tthe_impact_of_expected_MAF (df=0)]\n\t" +
						"[-impact_expected_gxe\tthe_impact_of_GxE_effects (df=0)]\n\t");
				System.exit(0);
			} else {
				String genotype_hdf5 = null;
				String annotation_file = null;
				String candidate_genes_file = null;
				String ethnicity_file = null;
				String ExAC_file = null;
				String gencode_file = null;
				String report_file = null;
				String tmp_folder = null;
				int num_case = 4;
				int round = 100;
				int num_top_gene = 10;
				int heterogeneity = 1;
				double prevalence_mean = 0.01;
				double penetrance_mean = 0.9;
				double prevalence_var = 0.0;
				double penetrance_var = 0.0;
				double ann_score_cutoff = 0.025;
				double ann_score_confidence = 0.95;//
				double individual_weight = 1.0;
				double general_pop_freq = 0.001;
				boolean genes_file_known = true;
				boolean compound_het = false;
				int parent_control = 0;
				int sibling_control = 0;
				int trio_control = 0;
				int multi_sibling_control = 0;
				boolean multi_affacted_sibling = false;

				double expected_causal_maf = 0;
				double impact_expected_causal_maf = 0;
				double impact_expected_gxe = 0;

				for (int k = 1; k < args.length; k++) {    //Setting arguments for SIMPEL function
					if (args[k].startsWith("-")) {
						if (args[k].equals("-population_genotypes")) genotype_hdf5 = args[k + 1];
						else if (args[k].equals("-out")) report_file = args[k + 1];
						else if (args[k].equals("-causal_gene_pool")) candidate_genes_file = args[k + 1];
						else if (args[k].equals("-population_pedigree")) ethnicity_file = args[k + 1];
						else if (args[k].equals("-all_genes")) gencode_file = args[k + 1];
						else if (args[k].equals("-mafs")) ExAC_file = args[k + 1];
						else if (args[k].equals("-pathogenicity")) annotation_file = args[k + 1];
						else if (args[k].equals("-tmp_folder")) tmp_folder = args[k + 1];
						else if (args[k].equals("-num_cases")) num_case = Integer.parseInt(args[k + 1]);
						else if (args[k].equals("-sim_rounds")) round = Integer.parseInt(args[k + 1]);
						else if (args[k].equals("-num_causal_genes")) heterogeneity = Integer.parseInt(args[k + 1]);
						else if (args[k].equals("-rank_success")) num_top_gene = Integer.parseInt(args[k + 1]);
						else if (args[k].equals("-prevalence_mean")) prevalence_mean = Double.parseDouble(args[k + 1]);
						else if (args[k].equals("-penetrance_mean")) penetrance_mean = Double.parseDouble(args[k + 1]);
						else if (args[k].equals("-prevalence_var")) prevalence_var = Double.parseDouble(args[k + 1]);
						else if (args[k].equals("-penetrance_var")) penetrance_var = Double.parseDouble(args[k + 1]);
						else if (args[k].equals("-min_score")) ann_score_cutoff = Double.parseDouble(args[k + 1]);
						else if (args[k].equals("-conf_score")) ann_score_confidence = Double.parseDouble(args[k + 1]);
						else if (args[k].equals("-weight_score")) individual_weight = Double.parseDouble(args[k + 1]);
						else if (args[k].equals("-pop_freq")) general_pop_freq = Double.parseDouble(args[k + 1]);
						else if (args[k].equals("-single_parents")) parent_control = Integer.parseInt(args[k + 1]);
						else if (args[k].equals("-both_parents")) trio_control = Integer.parseInt(args[k + 1]);
						else if (args[k].equals("-siblings")) sibling_control = Integer.parseInt(args[k + 1]);
						else if (args[k].equals("-multi_siblings"))
							multi_sibling_control = Integer.parseInt(args[k + 1]);
						else if (args[k].equals("-multi_affacted_sibling")) multi_affacted_sibling = true;
						else if (args[k].equals("-include_all_genes")) genes_file_known = false;
						else if (args[k].equals("-compound_het")) compound_het = true;
						else if (args[k].equals("-expected_causal_maf"))
							expected_causal_maf = Double.parseDouble(args[k + 1]);
						else if (args[k].equals("-impact_expected_causal_maf"))
							impact_expected_causal_maf = Double.parseDouble(args[k + 1]);
						else if (args[k].equals("-impact_expected_gxe"))
							impact_expected_gxe = Double.parseDouble(args[k + 1]);
						else {
							System.out.println("Sorry, the option " + args[k] + " is not supported. You could have misspelled the command.\n+" +
									"Please check and try again.");
							System.exit(0);
						}
					}
				}
				if (genotype_hdf5 == null || report_file == null || candidate_genes_file == null || ethnicity_file == null || gencode_file == null
						|| ExAC_file == null || annotation_file == null || tmp_folder == null) {    //Checking for Null entries in arguments
					System.out.println("Input or output files can't be null! JAWAMix5 has found the following problems:\n");
					if (genotype_hdf5 == null) System.out.println("genotype_file is missing!");
					if (report_file == null) System.out.println("report_file is missing!");
					if (candidate_genes_file == null) System.out.println("candidate_genes_file is missing!");
					if (ethnicity_file == null) System.out.println("ethnicity_file is missing!");
					if (gencode_file == null) System.out.println("gencode_file is missing!");
					if (ExAC_file == null) System.out.println("ExAC_file is missing!");
					if (annotation_file == null) System.out.println("annotation_file is missing!");
					if (tmp_folder == null) System.out.println("tmp_folder is missing!");
					System.exit(0);
				}
				if (sibling_control + parent_control > num_case) {
					System.out.println("Sorry, the number of patients having family controls cannot be larger than the total number of patients!");
					System.exit(0);
				}
				SimPEL ap = new SimPEL(gencode_file, candidate_genes_file, ethnicity_file, genotype_hdf5, ExAC_file,
						annotation_file, num_case, heterogeneity, general_pop_freq, penetrance_mean, penetrance_var, prevalence_mean, prevalence_var,
						parent_control, sibling_control, trio_control, multi_sibling_control, multi_affacted_sibling, individual_weight,
						ann_score_cutoff, ann_score_confidence, genes_file_known, compound_het,
						expected_causal_maf, impact_expected_causal_maf, impact_expected_gxe,
						tmp_folder);
				ap.test_power(round, num_top_gene, report_file);
			}
		}


		//#####################################################################################
		// Kinship Function Initialization
		//#####################################################################################
		else if (function.equals("kinship")) {    // Calculating the Kinship Matrix
			if (args.length == 1) {
				System.out.println("Compute IBS/RRM kinship matrix required for other analyses.");
				System.out.println("Usage: \n\t<-ig\tinput_genotype_file>\n\t" +
						"<-o\toutput_prefix>\n\t" +
						"[-w\ttiling_window_size(bp) (df=WG)]\n\t" +
						"[-m\tmethod (df=IBS)]\n\t" +
						"[-maf\tmin-MAF (df=0)]\n\t" +
						"[-scale\tmax_genotype_coding (df=2)]\n\t");
				System.exit(0);
			} else {
				String input = null, output_folder = null;
				int tiling_window_size = -1;
				double scale = 2.0, maf = 0;
				String method = "IBS";
				for (int k = 1; k < args.length; k++) {
					if (args[k].startsWith("-")) {    // Setting arguments for kinship function
						if (args[k].equals("-ig")) input = args[k + 1];
						else if (args[k].equals("-o")) output_folder = args[k + 1];
						else if (args[k].equals("-w")) tiling_window_size = Integer.parseInt(args[k + 1]);
						else if (args[k].equals("-m")) method = args[k + 1];
						else if (args[k].equals("-maf")) maf = Double.parseDouble(args[k + 1]);
						else if (args[k].equals("-scale")) scale = Double.parseDouble(args[k + 1]);
					}
				}
				if (input == null || output_folder == null) {
					System.out.println("Input or output-folder can't be null!");
				} else {
					if (tiling_window_size != -1) {
						LocalKinshipAnalyzer analyzer = new LocalKinshipAnalyzer(input, tiling_window_size, null);
						analyzer.calculating_kinship_tiling_windows(output_folder, scale);
					} else { //	tiling_window_size==-1, i.e., whole-genome
						if (method.equals("IBS")) {
							VariantsDouble calculator = new VariantsDouble(input);
							System.out.println("Calculating global IBS kinship for " + input);
							calculator.calculate_raw_ibs_kinship(output_folder + ".raw.IBS", scale, maf);
							KinshipMatrix.re_scale_kinship_matrix(output_folder + ".raw.IBS", output_folder + ".rescaled.IBS");
						} else if (method.equals("RRM")) {
							VariantsDouble calculator = new VariantsDouble(input);
							System.out.println("Calculating global RRM kinship for " + input);
							calculator.calculate_WG_RRM_kinship(output_folder + ".RRM", scale, maf);
							//VariantsDouble.re_scale_kinship_matrix(output_folder+".kinship.RRM", output_folder+".kinship.rescaled.RRM");
						} else {
							System.out.println("Method " + method + " is not supported. It can only be IBS or RRM.");
						}
					}
				}
			}


			//#####################################################################################
			// Char2Num Function Initialization
			// #####################################################################################
		}else if(function.equals("char2num")){
			if(args.length==1){
				System.out.println("Convert char-coded genotype CSV file to number-coded genotype CSV file ready for \"import\"");
				System.out.println("Usage: \n\t<-ig\tinput_genotype_file>\n\t" +
						"<-o\toutput_file>\n\t");
				System.exit(0);
			}else{
				String input=null, output=null;
				for(int k=1;k<args.length;k++){
					if(args[k].startsWith("-")){
						if(args[k].equals("-ig"))input=args[k+1];
						else if(args[k].equals("-o"))output=args[k+1];
			}
				}if(input==null||output==null){
					System.out.println("Input or output can't be null!");
				}else{
					VariantsDouble.char2num(input, output);
				}
			}


			//#####################################################################################
			// tped2num Function Initialization
			//#####################################################################################
		}else if(function.equals("tped2num")){
			if(args.length==1){
				System.out.println("Convert char-coded genotype PLINK tped file to number-coded genotype CSV file ready for \"import\"");
				System.out.println("Usage: \n\t<-tped\tinput_tped_file>\n\t" +
						"<-tfam\tinput_tfam_file>\n\t" +
						"<-o\toutput_file_prefix>\n\t");
				System.exit(0);
			}else{
				String input_tped=null, input_tfam=null,  output=null;
				for(int k=1;k<args.length;k++){
					if(args[k].startsWith("-")){
						if(args[k].equals("-tped"))input_tped=args[k+1];
						else if(args[k].equals("-tfam"))input_tfam=args[k+1];
						else if(args[k].equals("-o"))output=args[k+1];
			}
				}if(input_tped==null||input_tfam==null||output==null){
					System.out.println("Input or output can't be null!");
				}else{
					VariantsDouble.tped2csv_num(input_tfam, input_tped, output+".char.csv", output+".num.csv");
				}
			}

			//#####################################################################################
			// heritability Function Initialization
			//#####################################################################################
		}else if(function.equals("heritability")){
			if(args.length==1){
				//TODO: Add description of format of phenotype file expected
				System.out.println("Estimate pseudo-heritability for a large number of phenotypes for the same samples.");
				System.out.println("Usage: \n\t<-ik\tinput_kinship_file>\n\t" +
						"<-ip\tinput_phenotype_file>\n\t"+
						"<-o\toutput_file>\n\t"+
						"[-m\tmethod<fullrank|lowrank> (df=fullrank)]\n\t");
				System.exit(0);
			}else{
				String input_kinship=null, input_pheno=null,output=null;
				String method="fullrank";
				for(int k=1;k<args.length;k++){
					if(args[k].startsWith("-")){
						if(args[k].equals("-ik"))input_kinship=args[k+1];
						else if(args[k].equals("-ip"))input_pheno=args[k+1];
						else if(args[k].equals("-o"))output=args[k+1];
						else if(args[k].equals("-m"))method=args[k+1];
			}
				}if(input_kinship==null||input_pheno==null||output==null){
					System.out.println("Input or output can't be null!");
				}else{
					if(method.equals("fullrank"))
						VarianceComponent.heritability_fullrank(new KinshipMatrix(input_kinship), new MultiPhenotype(input_pheno), output);
					else if(method.equals("lowrank"))
						VarianceComponent.heritability_lowrank(input_kinship, new MultiPhenotype(input_pheno), output);
				}
			}

			//#####################################################################################
			// import Function Initialization
			//#####################################################################################
		}else if(function.equals("import")){
			if(args.length==1){
				//TODO: Add description of format of genotype file expected
				System.out.println("Import genotype from .csv file to .hdf5 file.");
				System.out.println("Usage: \n\t<-ig\tinput_genotype_file>\n\t" +
						"<-o\tout_put_hdf5_file>\n\t" +
						"[-type\ttype<byte|double> (df=double)]\n\t" +
						"[-b\tblock_size (df=5000)]\n\t");
				System.exit(0);
			}else{
				String input=null, output=null, import_type="double";
				int block_size=5000;
				for(int k=1;k<args.length;k++){
					if(args[k].startsWith("-")){
						if(args[k].equals("-ig"))input=args[k+1];
						else if(args[k].equals("-o"))output=args[k+1];
						else if(args[k].equals("-b"))block_size=Integer.parseInt(args[k+1]);
						else if(args[k].equals("-type"))import_type=args[k+1];
					}
				}if(input==null||output==null){
					System.out.println("Input or output can't be null!");
				}else{
					if(import_type.equals("double"))
						VariantsDouble.importCSV(input, output, block_size);
					else if(import_type.equals("byte"))
						VariantsByte.importCSV(input, output, block_size);
				}
			}

			//#####################################################################################
			// emmax Function Initialization
			//#####################################################################################
		}else if(function.equals("emmax")){
			if(args.length==1){
				//TODO: Add description of format of phenotype file expected
				System.out.println("Run EMMAX for phenotype(s).");
				System.out.println("Usage: \n\t<-ig\tinput_genotype_file>\n\t" +
						"<-ip\tphenotype_file>\n\t" +
						"<-o\toutput_folder>\n\t" +
						"<-ik\tkinship_file>\n\t" +
						"[-p\t<pvalue_after_multi.correct.> (df=1000)]\n\t" +
						"[-index\t<phenotype_index> (df=ALL, start from zero)]\n\t" +
						"[-min_size\t<min_sample_size> (df=100)]\n\t" +
						"[-maf\t<maf_threshold_plot> (df=0.05)]\n\t" +
						"[-plot\t<1|0> (df=1)]\n\t" +
						"[-ri\t<region infomation in format: r0 chr start end> (require pheno_index)]\n\t" +
						"[-ocma\t<to use out of core matrix calculations>]");
				System.exit(0);
			} else {
				String input_geno = "", input_pheno = null, output_folder = null, kinship = null, region_info = null;
				;
				double p_after_corr = 1000, maf_threshold_plot = 0.05;
				int phe_index = -1, min_sample_size = 100, round = 1;
				boolean plot = true;
				boolean ocma = false;
				for (int k = 1; k < args.length; k++) {
					if (args[k].startsWith("-")) {
						if (args[k].equals("-ig")) input_geno = args[k + 1];
						else if (args[k].equals("-ip")) input_pheno = args[k + 1];
						else if (args[k].equals("-ik")) kinship = args[k + 1];
						else if (args[k].equals("-o")) output_folder = args[k + 1];
						else if (args[k].equals("-index")) phe_index = Integer.parseInt(args[k + 1]);
						else if (args[k].equals("-p")) p_after_corr = Double.parseDouble(args[k + 1]);
						else if (args[k].equals("-min_size")) min_sample_size = Integer.parseInt(args[k + 1]);
						else if (args[k].equals("-maf")) maf_threshold_plot = Double.parseDouble(args[k + 1]);
						else if (args[k].equals("-plot")) {
							if (args[k + 1].equals("0")) plot = false;
						} else if (args[k].equals("-ri")) region_info = args[k + 1];
						else if (args[k].equals("-ocma")) ocma = false;
					}
				}
				if (input_geno == null || input_pheno == null || kinship == null || output_folder == null) {
					System.out.println("Input or output can't be null!");
				} else {
					if (phe_index == -1) {
						//TODO: Add check for output folder existing and create if not
						if (ocma = true) {
							EMMAX_OCMA.emmax_analysis(input_geno, input_pheno, kinship, output_folder, round, p_after_corr, min_sample_size,
									maf_threshold_plot, plot);
						} else {
							EMMAX.emmax_analysis(input_geno, input_pheno, kinship, output_folder, round, p_after_corr, min_sample_size,
									maf_threshold_plot, plot);
						}
					} else {
						//TODO: Add check for output folder existing and create if not
						if (region_info == null) {
							if (ocma = true) {
								EMMAX_OCMA.emmax_analysis(input_geno, input_pheno, kinship, output_folder, round, p_after_corr, min_sample_size,
										phe_index, maf_threshold_plot, plot);
							} else {
								EMMAX.emmax_analysis(input_geno, input_pheno, kinship, output_folder, round, p_after_corr, min_sample_size,
										phe_index, maf_threshold_plot, plot);
							}
						} else {
							int[][][] region = Region_EMMAX.read_regions_from_file(region_info);
							EMMAX.emmax_analysis_regions(input_geno, input_pheno, kinship, output_folder, 
									p_after_corr, min_sample_size, phe_index, maf_threshold_plot, plot, region);
							}
						}
					}
				}

			//#####################################################################################
			// emmax_stepwise Function Initialization
			//#####################################################################################
		} else if (function.equals("emmax_stepwise")) {
			if (args.length == 1) {
				//TODO: Add description of format of phenotype file expected
				System.out.println("Run EMMAX-based stepwise regression for phenotype(s).");
				System.out.println("Usage: \n\t<-ig\tinput_genotype_file>\n\t" +
						"<-ip\tphenotype_file>\n\t" +
						"<-o\toutput_folder>\n\t" +
						"<-ik\tkinship_file>\n\t" +
						"[-r\t<round> (df=2)]\n\t"+
						"[-p\t<pvalue_after_multi.correct.> (df=1000)]\n\t" +
						"[-index\t<phenotype_index> (df=ALL, start from zero)]\n\t" +
						"[-min_size\t<min_sample_size> (df=100)]\n\t" +
						"[-maf\t<maf_threshold_plot> (df=0.05)]\n\t"+
						"[-plot\t<1|0> (df=1)]");
				System.exit(0);
			}else{
				String input_geno="", input_pheno=null, output_folder=null, kinship=null;
				double p_after_corr=1000, maf_threshold_plot=0.05;
				int phe_index=-1, min_sample_size=100, round=2;
				boolean plot=true;
				for(int k=1; k<args.length; k++){
					if(args[k].startsWith("-")){
						if(args[k].equals("-ig"))input_geno=args[k+1];
						else if(args[k].equals("-ip"))input_pheno=args[k+1];
						else if(args[k].equals("-ik"))kinship=args[k+1];
						else if(args[k].equals("-o"))output_folder=args[k+1];
						else if(args[k].equals("-r"))round=Integer.parseInt(args[k+1]);
						else if(args[k].equals("-index"))phe_index=Integer.parseInt(args[k+1]);
						else if(args[k].equals("-p"))p_after_corr=Double.parseDouble(args[k+1]);
						else if(args[k].equals("-min_size"))min_sample_size=Integer.parseInt(args[k+1]);
						else if(args[k].equals("-maf"))maf_threshold_plot=Double.parseDouble(args[k+1]);
						else if(args[k].equals("-plot")){if(args[k+1].equals("0"))plot=false;}
					}
				}if(input_geno==null||input_pheno==null|| kinship==null||output_folder==null){
					System.out.println("Input or output can't be null!");
				}else{
					if(phe_index==-1)
						EMMAX.emmax_analysis(input_geno, input_pheno, kinship, output_folder, round, p_after_corr, min_sample_size,
								maf_threshold_plot, plot);
					else{
						EMMAX.emmax_analysis(input_geno, input_pheno, kinship, output_folder, round, p_after_corr, min_sample_size,
								phe_index, maf_threshold_plot, plot);
					}
				}
			}


			//#####################################################################################
			// emmax_select Function Initialization
			//#####################################################################################
		}else if(function.equals("emmax_select")){
			if(args.length==1){
				//TODO: Add description of format of phenotype file expected
				System.out.println("Run EMMAX one time on selected variants.");
				System.out.println("Usage: \n\t<-ig\tinput_genotype_file>\n\t" +
						"<-chrs\tchr1-chr2...-chrn>\n\t" +
						"<-locations\tloc1-loc2...-locn>\n\t" +
						"<-ip\tphenotype_file>\n\t" +
						"<-o\toutput_folder>\n\t" +
						"<-ik\tkinship_file>\n\t"+
						"[-index\t<phenotype_index> (df=0)]\n");
				System.exit(0);
			}else{
				String input_geno="", input_pheno= null, output_file=null, kinship=null;
				int[] chrs=null, locs=null;
				int phe_index=0;
				for(int k=1; k<args.length; k++){
					if(args[k].startsWith("-")){
						if(args[k].equals("-ig"))input_geno=args[k+1];
						else if(args[k].equals("-ip"))input_pheno=args[k+1];
						else if(args[k].equals("-ik"))kinship=args[k+1];
						else if(args[k].equals("-o"))output_file=args[k+1];
						else if(args[k].equals("-index"))phe_index=Integer.parseInt(args[k+1]);
						else if(args[k].equals("-chr")){
							String[] string_chrs=args[k+1].split("-");
							chrs=new int[string_chrs.length];
							for(int kk=0; kk<chrs.length; kk++)chrs[kk]=Integer.parseInt(string_chrs[kk]);
						}else if(args[k].equals("-locations")){
							String[] string_locs=args[k+1].split("-");
							if(string_locs.length!=chrs.length){
								System.out.println("Number of locations has to be the same as the number of chromosomes!\n" +
										"JAWAMix5 will now exit.");
								System.exit(0);
							}
							locs=new int[string_locs.length];
							for(int kk=0; kk<locs.length; kk++)locs[kk]=Integer.parseInt(string_locs[kk]);
						}
					}
				}if(input_geno==null||input_pheno==null|| output_file==null || kinship==null ||chrs==null||locs==null){
					System.out.println("Input or output can't be null!");
				}else{
					EMMAX.emmax_select(input_geno, input_pheno, kinship, phe_index, output_file, chrs, locs);
				}
			}


			//#####################################################################################
			// lm_select Function Initialization
			//#####################################################################################
		}else if(function.equals("lm_select")){
			if(args.length==1){
				//TODO: Add description of format of phenotype file expected
				System.out.println("Run linear regression one time on selected variants.");
				System.out.println("Usage: \n\t<-ig>\n\t" +
						"<-chrs\tchr1-chr2...-chrn>\n\t" +
						"<-locations\tloc1-loc2...-locn>\n\t" +
						"<-ip\tphenotype_file>\n\t" +
						"<-o\toutput_folder>\n\t"+
						"[-index\t<phenotype_index> (df=0)]\n");
				System.exit(0);
			}else{
				String input_geno="", input_pheno= null, output_file=null;
				int[] chrs=null, locs=null;
				int phe_index=0;
				for(int k=1; k<args.length; k++){
					if(args[k].startsWith("-")){
						if(args[k].equals("-ig"))input_geno=args[k+1];
						else if(args[k].equals("-ip"))input_pheno=args[k+1];
						else if(args[k].equals("-o"))output_file=args[k+1];
						else if(args[k].equals("-index"))phe_index=Integer.parseInt(args[k+1]);
						else if(args[k].equals("-chr")){
							String[] string_chrs=args[k+1].split("-");
							chrs=new int[string_chrs.length];
							for(int kk=0; kk<chrs.length; kk++)chrs[kk]=Integer.parseInt(string_chrs[kk]);
						}else if(args[k].equals("-locations")){
							String[] string_locs=args[k+1].split("-");
							if(string_locs.length!=chrs.length){
								System.out.println("Number of locations has to be the same as the number of chromosomes!\n"+
										"JAWAMix5 will now exit.");
								System.exit(0);
							}
							locs=new int[string_locs.length];
							for(int kk=0; kk<locs.length; kk++)locs[kk]=Integer.parseInt(string_locs[kk]);
						}
					}
				}if(input_geno==null||input_pheno==null|| output_file==null || chrs==null||locs==null){
					System.out.println("Input or output can't be null!");
				}else{
					LM.lm_select(input_geno, input_pheno, phe_index, output_file, chrs, locs);
				}
			}


			//#####################################################################################
			// lm_stepwise Function Initialization
			//#####################################################################################
		}else if(function.equals("lm_stepwise")){
			if(args.length==1){
				//TODO: Add description of format of phenotype file expected
				System.out.println("Run stepwise regression for phenotype(s) without population structure accounted.");
				System.out.println("Usage: \n\t<-ig\tinput_genotype_file>\n\t" +
						"<-ip\tphenotype_file>\n\t" +
						"<-o\toutput_folder>\n\t" +
						"[-r\t<round> (df=2)]\n\t"+
						"[-p\t<pvalue_after_multi.correct.> (df=1000)]\n\t" +
						"[-index\t<phenotype_index> (df=ALL, start from zero)]\n\t" +
						"[-min_size\t<min_sample_size> (df=100)]\n\t" +
						"[-maf\t<maf_threshold_plot> (df=0.05)]\n\t");
				System.exit(0);
			}else{
				String input_geno="", input_pheno=null, output_folder=null;
				double p_after_corr=1000, maf_threshold_plot=0.05;
				int phe_index=-1, min_sample_size=100, round=2;
				for(int k=1; k<args.length; k++){
					if(args[k].startsWith("-")){
						if(args[k].equals("-ig"))input_geno=args[k+1];
						else if(args[k].equals("-ip"))input_pheno=args[k+1];
						else if(args[k].equals("-o"))output_folder=args[k+1];
						else if(args[k].equals("-r"))round=Integer.parseInt(args[k+1]);
						else if(args[k].equals("-index"))phe_index=Integer.parseInt(args[k+1]);
						else if(args[k].equals("-p"))p_after_corr=Double.parseDouble(args[k+1]);
						else if(args[k].equals("-min_size"))min_sample_size=Integer.parseInt(args[k+1]);
						else if(args[k].equals("-maf"))maf_threshold_plot=Double.parseDouble(args[k+1]);

					}
				}if(input_geno==null||input_pheno==null|| output_folder==null){
					System.out.println("Input or output can't be null!");
				}else{
					if(phe_index==-1){
						LM.stepwise_analysis(input_geno, input_pheno, output_folder, round, p_after_corr, min_sample_size, maf_threshold_plot);
					}else{
						LM.stepwise_analysis(input_geno, input_pheno, output_folder, round, p_after_corr, min_sample_size, phe_index, maf_threshold_plot);
					}
				}
			}


			//#####################################################################################
			// lm Function Initialization
			//#####################################################################################
		}else if(function.equals("lm")){
			if(args.length==1){
				//TODO: Add description of format of phenotype file expected
				System.out.println("Run linear regression for phenotype(s) without population structure accounted.");
				System.out.println("Usage: \n\t<-ig\tinput_genotype_file>\n\t" +
						"<-ip\tphenotype_file>\n\t" +
						"<-o\toutput_folder>\n\t" +
						"[-p\t<pvalue_after_multi.correct.> (df=1000)]\n\t" +
						"[-index\t<phenotype_index> (df=ALL, start from zero)]\n\t" +
						"[-min_size\t<min_sample_size> (df=100)]\n\t" +
						"[-maf\t<maf_threshold_plot> (df=0.05)]\n\t");System.exit(0);
			}else{
				String input_geno="", input_pheno=null, output_folder=null;
				double p_after_corr=1000, maf_threshold_plot=0.05;
				int phe_index=-1, min_sample_size=100;
				final int round=1;
				for(int k=1; k<args.length; k++){
					if(args[k].startsWith("-")){
						if(args[k].equals("-ig"))input_geno=args[k+1];
						else if(args[k].equals("-ip"))input_pheno=args[k+1];
						else if(args[k].equals("-o"))output_folder=args[k+1];
						else if(args[k].equals("-index"))phe_index=Integer.parseInt(args[k+1]);
						else if(args[k].equals("-p"))p_after_corr=Double.parseDouble(args[k+1]);
						else if(args[k].equals("-min_size"))min_sample_size=Integer.parseInt(args[k+1]);
						else if(args[k].equals("-maf"))maf_threshold_plot=Double.parseDouble(args[k+1]);
					}
				}if(input_geno==null||input_pheno==null|| output_folder==null){
					System.out.println("Input or output can't be null!");
				}else{
					if(phe_index==-1){
						LM.stepwise_analysis(input_geno, input_pheno, output_folder, round, p_after_corr, min_sample_size, maf_threshold_plot);
					}else{
						LM.stepwise_analysis(input_geno, input_pheno, output_folder, round, p_after_corr, min_sample_size, phe_index, maf_threshold_plot);
					}
				}
			}


			//#####################################################################################
			// local variance component analysis Function Initialization
			//#####################################################################################
		}else if(function.equals("local")){
			if(args.length==1){
				//TODO: Add description of format of phenotype file expected
				System.out.println("Run local variance component analysis.");
				System.out.println("Usage: \n\t<-ig\tinput_genotype_file>\n\t" +
						"<-ip\tphenotype_file>\n\t" +
						"<-o\toutput_folder>\n\t" +
						"<-w\ttiling_window_size(bp)>\n\t" +
						"<-ik_l\tlocal_kinship_files_folder[not needed if -method is \"transform\"]>\n\t" +
						"<-ik_g\tglobal_kinship_file>\n\t" +
						"[-method\t<full|grid|transform> (df=\"transform\")]\n\t" +
						"[-ir\t<input_region_file>]\n\t" +
						"[-index\t<phenotype_index> (df=ALL, start from zero)]\n\t" +
						"[-min_size\t<min_sample_size> (df=40)]\n\t" +
						"[-step\t<grid_size> (df=0.01)]\n\t"+
						"[-plot\t<1|0> (df=1)]");
				System.exit(0);
			}else{
				String input_geno=null, input_pheno=null, output_folder=null, input_region=null,
						local_kinship_files_folder=null, global_kinship_file=null;
				int win_size=-1;
				int the_phe_index=-1, min_sample_size=40;
				double step=0.01;
				boolean plot=true;
				String method="transform";
				for(int k=1; k<args.length; k++){
					if(args[k].startsWith("-")){
						if(args[k].equals("-ig"))input_geno=args[k+1];
						else if(args[k].equals("-ip"))input_pheno=args[k+1];
						else if(args[k].equals("-ik_l"))local_kinship_files_folder=args[k+1];
						else if(args[k].equals("-ik_g"))global_kinship_file=args[k+1];
						else if(args[k].equals("-o"))output_folder=args[k+1];
						else if(args[k].equals("-method"))method=args[k+1];
						else if(args[k].equals("-w"))win_size=Integer.parseInt(args[k+1]);
						else if(args[k].equals("-index"))the_phe_index=Integer.parseInt(args[k+1]);
						else if(args[k].equals("-min_size"))min_sample_size=Integer.parseInt(args[k+1]);
						else if(args[k].equals("-step"))step=Double.parseDouble(args[k+1]);
						else if(args[k].equals("-plot")){if(args[k+1].equals("0"))plot=false;} else if(args[k].equals("-ir"))input_region=args[k+1];
					}
				}if(input_geno==null||input_pheno==null|| output_folder==null
						||global_kinship_file==null || (win_size==-1 && input_region==null)){
					System.out.println("Window size or input or output can't be null!");

				}else{
					LocalKinshipAnalyzer local_k=new LocalKinshipAnalyzer(input_geno, win_size, null);
					MultiPhenotype phenotypeS=new MultiPhenotype(input_pheno);
					if(the_phe_index==-1){
						System.out.println("Running all phenotypes? It is suggested to specify a phenotype index. " +
								"Otherwise it may be slow." +
								"\nLet us have a try!");
						for(int phe_index=0; phe_index<phenotypeS.num_of_pheno; phe_index++){
							if(input_region==null){
								String out_phe_file = output_folder+"Local_VO."+phe_index+"."+phenotypeS.phenotypes[phe_index].phe_id+
										".w"+win_size+".csv";
								local_k.local_VO_wins(phenotypeS.phenotypes[phe_index], input_geno, global_kinship_file,
										local_kinship_files_folder, out_phe_file, step, plot, method, min_sample_size);
							} else{
								String out_phe_file=output_folder+"Local_VO."+phe_index+"."+phenotypeS.phenotypes[phe_index].phe_id+".r"+".csv";
								local_k.local_VO_regions(phenotypeS.phenotypes[phe_index], input_geno, global_kinship_file, input_region,
										out_phe_file, plot, min_sample_size);
							}
						}
					}else{
						String out_phe_file=output_folder+"Local_VO."+the_phe_index+"."+phenotypeS.phenotypes[the_phe_index].phe_id+".w"+win_size+".csv";
						if(input_region==null){
							local_k.local_VO_wins(phenotypeS.phenotypes[the_phe_index], input_geno, global_kinship_file, local_kinship_files_folder,
									out_phe_file, step, plot, method, min_sample_size);
						}else{
							local_k.local_VO_regions(phenotypeS.phenotypes[the_phe_index], input_geno, global_kinship_file, input_region,
									out_phe_file, plot, min_sample_size);
						}
					}
				}
			}


			//#####################################################################################
			// compound - variance component analysis Function Initialization
			//#####################################################################################
		}else if(function.equals("compound")){
			if(args.length==1){
				//TODO: Add description of format of phenotype file expected
				System.out.println("Run variance component analysis for a compound.");
				System.out.println("Usage: \n\t<-ig\tinput_genotype_file>\n\t" +
						"<-ip\tphenotype_file>\n\t" +
						"<-o\toutput_folder>\n\t" +
						"<-ik_l\tcompound_kinship_files_folder[not needed if -method is \"transform\"]>\n\t" +
						"<-ik_g\tglobal_kinship_file>\n\t" +
						"[-method\t<full|grid|transform> (df=\"transform\")]\n\t" +
						"[-ic\t<input_compound_file>]\n\t" +
						"[-index\t<phenotype_index> (df=ALL, start from zero)]\n\t" +
						"[-min_size\t<min_sample_size> (df=40)]\n\t" +
						"[-step\t<grid_size> (df=0.01)]\n\t"+
						"[-plot\t<1|0> (df=1)]");
				System.exit(0);
			}else{
				String input_geno=null, input_pheno=null, output_folder=null, input_compound=null,
						compound_kinship_files_folder=null, global_kinship_file=null;
				int the_phe_index=-1, min_sample_size=40;
				boolean plot=true;
				for(int k=1; k<args.length; k++){
					if(args[k].startsWith("-")){
						if(args[k].equals("-ig"))input_geno=args[k+1];
						else if(args[k].equals("-ip"))input_pheno=args[k+1];
						else if(args[k].equals("-ik_l"))compound_kinship_files_folder=args[k+1];
						else if(args[k].equals("-ik_g"))global_kinship_file=args[k+1];
						else if(args[k].equals("-o"))output_folder=args[k+1];
						else if(args[k].equals("-index"))the_phe_index=Integer.parseInt(args[k+1]);
						else if(args[k].equals("-min_size"))min_sample_size=Integer.parseInt(args[k+1]);
						else if(args[k].equals("-plot")){if(args[k+1].equals("0"))plot=false;} else if(args[k].equals("-ic"))input_compound=args[k+1];
					}
				}if(input_geno==null||input_pheno==null|| output_folder==null
						||global_kinship_file==null || input_compound==null){
					System.out.println("input or output can't be null!");

				}else{
					CompoundAnalyzer local_k=new CompoundAnalyzer(input_geno, input_compound);
					MultiPhenotype phenotypeS=new MultiPhenotype(input_pheno);
					if(the_phe_index==-1){
						System.out.println("Running all phenotypes? It is suggested to specify a phenotype index. " +
								"Otherwise it may be slow." +
								"\nLet us have a try!");
						for(int phe_index=0; phe_index<phenotypeS.num_of_pheno; phe_index++){
							String out_phe_file=output_folder+"Compound_VO."+phe_index+"."+phenotypeS.phenotypes[phe_index].phe_id+".r"+".csv";
							local_k.compound_VO(phenotypeS.phenotypes[phe_index], input_geno, global_kinship_file, out_phe_file, plot, min_sample_size);
						}
					}else{
						String out_phe_file=output_folder+"Compound_VO."+the_phe_index+"."+phenotypeS.phenotypes[the_phe_index].phe_id+".csv";
						local_k.compound_VO(phenotypeS.phenotypes[the_phe_index], input_geno, global_kinship_file, out_phe_file, plot, min_sample_size);
					}
				}
			}

			//#####################################################################################
			// rare variant analysis Function Initialization
			//#####################################################################################
		}else if(function.equals("rare")){
			if(args.length==1){
				//TODO: Add description of format of phenotype file expected
				System.out.println("Run rare variants analysis leveraging (or not) potential synthetic associations, " +
						"with population structure accounted (or not)");
				System.out.println("Usage:\n\t" +
						"<-ig\tinput_genotype_file>\n\t" +
						"<-ip\tphenotype_file>\n\t" +
						"<-ik\tkinship_file>\n\t" +
						"<-is\tinput_single_marker_results (preferably by emmax)>\n\t" +
						"<-ir\tregions_grouping_file>\n\t" +
						"<-o\toutput_prefix>\n\t" +
						"[-index\t<phenotype_index> (df=ALL, start from zero)]\n\t" +
						"[-ld\t<LD_threshold(D-prime)> (df=0.8)]\n\t" +
						"[-rare\t<rare_threshold> (df=0.01)]\n\t" +
						"[-syn_pvalue\t<potential syntehtic association pvalue threshold> (df=0.00001)]\n\t" +
						"[-syn_maf\t<potential syntehtic association MAF threshold> (df=0.00001)]\n\t" +
						"[-dist\t<distance2region>(df=100000)]");
				System.exit(0);
			}else{
				String input_geno=null, input_pheno=null, output_prefix=null, kinship_file=null, input_single_marker_results_file=null,
						input_region_file=null;
				double ld_threshold=0.8, rare_threshold=0.01, synthetic_pvalue_threshold=0.00001, synthetic_maf_threhold=0.1;
				int the_phe_index=-1, distance2region=100000;
				for(int k=1; k<args.length; k++){
					if(args[k].startsWith("-")){
						if(args[k].equals("-ig"))input_geno=args[k+1];
						else if(args[k].equals("-ip"))input_pheno=args[k+1];
						else if(args[k].equals("-is"))input_single_marker_results_file=args[k+1];
						else if(args[k].equals("-ir"))input_region_file=args[k+1];
						else if(args[k].equals("-ik"))kinship_file=args[k+1];
						else if(args[k].equals("-o"))output_prefix=args[k+1];
						else if(args[k].equals("-index"))the_phe_index=Integer.parseInt(args[k+1]);
						else if(args[k].equals("-dist"))distance2region=Integer.parseInt(args[k+1]);
						else if(args[k].equals("-ld"))ld_threshold=Double.parseDouble(args[k+1]);
						else if(args[k].equals("-rare"))rare_threshold=Double.parseDouble(args[k+1]);
						else if(args[k].equals("-syn_pvalue"))synthetic_pvalue_threshold=Double.parseDouble(args[k+1]);
						else if(args[k].equals("-syn_maf"))synthetic_maf_threhold=Double.parseDouble(args[k+1]);
					}
				}if(input_geno==null||input_pheno==null|| kinship_file==null||output_prefix==null || input_single_marker_results_file==null
						|| input_region_file==null){
					System.out.println("Input or output can't be null!");
				}else{
					MultiPhenotype phenotypeS=new MultiPhenotype(input_pheno);
					if(the_phe_index==-1)
						for(int phe_index=0; phe_index<phenotypeS.num_of_pheno; phe_index++){
							RareAnalyzerSynthetic rare_analyzer=new RareAnalyzerSynthetic(input_single_marker_results_file,
									synthetic_pvalue_threshold, synthetic_maf_threhold,
									input_geno, input_region_file, distance2region, phenotypeS.phenotypes[phe_index], kinship_file,
									ld_threshold, rare_threshold);
							rare_analyzer.rare_association_4in1(rare_threshold, output_prefix);
						}
					else{
						RareAnalyzerSynthetic rare_analyzer=new RareAnalyzerSynthetic(input_single_marker_results_file,
								synthetic_pvalue_threshold, synthetic_maf_threhold,
								input_geno, input_region_file, distance2region, phenotypeS.phenotypes[the_phe_index], kinship_file,
								ld_threshold, rare_threshold);
						rare_analyzer.rare_association_4in1(rare_threshold, output_prefix);
					}
				}
			}


			//#####################################################################################
			// rare-aggregate variants analysis Function Initialization
			//#####################################################################################
		}else if(function.equals("rare-aggregate")){
			if(args.length==1){
				//TODO: Add description of format of phenotype file expected
				System.out.println("Run rare variants analysis using aggregate test" +
						"with population structure accounted.");
				System.out.println("Usage:\n\t" +
						"<-ig\tinput_genotype_file>\n\t" +
						"<-ip\tphenotype_file>\n\t" +
						"<-ik\tkinship_file>\n\t" +
						"<-ir\tinput_region_file>\n\t" +
						"<-o\toutput_prefix>\n\t" +
						"[-index\t<phenotype_index> (df=ALL, start from zero)]\n\t" +
						"[-rare\t<rare_threshold> (df=0.01)]\n");
				System.exit(0);
			}else{
				String input_geno=null, input_pheno=null, output_prefix=null, kinship_file=null, input_region_file=null;
				double rare_threshold=0.01;
				int win_size=50000;
				int the_phe_index=-1;
				for(int k=1; k<args.length; k++){
					if(args[k].startsWith("-")){
						if(args[k].equals("-ig"))input_geno=args[k+1];
						else if(args[k].equals("-ip"))input_pheno=args[k+1];
						else if(args[k].equals("-ir"))input_region_file=args[k+1];
						else if(args[k].equals("-ik"))kinship_file=args[k+1];
						else if(args[k].equals("-o"))output_prefix=args[k+1];
						else if(args[k].equals("-index"))the_phe_index=Integer.parseInt(args[k+1]);
						else if(args[k].equals("-rare"))rare_threshold=Double.parseDouble(args[k+1]);
						else System.out.println("Option "+args[k]+" is undefined.");
					}
				}if(input_geno==null||input_pheno==null|| kinship_file==null||output_prefix==null||input_region_file==null){
					System.out.println("Input or output can't be null!");
				}else{
					MultiPhenotype phenotypeS=new MultiPhenotype(input_pheno);
					if(the_phe_index==-1)
						for(int phe_index=0; phe_index<phenotypeS.num_of_pheno; phe_index++){
							RareAnalyzerAggregate rare_analyzer=new RareAnalyzerAggregate(input_geno, input_region_file, phenotypeS.phenotypes[phe_index], kinship_file,
									rare_threshold);
//							rare_analyzer.rare_association(rare_threshold, output_prefix);
						}
					else{
						RareAnalyzerAggregate rare_analyzer=new RareAnalyzerAggregate(input_geno, input_region_file, phenotypeS.phenotypes[the_phe_index], kinship_file,
								rare_threshold);
//						rare_analyzer.rare_association(rare_threshold, output_prefix);
					}
				}
			}


			//#####################################################################################
			// NAM (Nested Association Mapping) Function Initialization
			//#####################################################################################
		}else if(function.equals("nam")){
			System.out.println("Nested Association Mapping."+"\n");
			if(args.length<10){
				System.out.println("Usage: \n\t"
						+"<-p\tlist of RIL_pedigreee files folder>\n\t"
						+ "<-c\tCross_ID>\n\t"
						+ "<-fg\tfounder whole genome genotype file>\n\t"
						+"<-ril\tRIL phenotype>\n\t"
						+ "<-o\tOutput-prefix>\n\t"
						+"[-b\tblock size in hdf5 (df=5000)]\n\t"
						+"[-index\tphenotype_index (df=ALL, start from zero)]\n\t"
						+"[-r\tround (df=2)]\n\t"
						+"[-p\tpvalue_after_multi.correction (df=1000)]\n\t"
						+"[-maf\tmaf_threshold_plot (df=0.05)]\n\t");
				System.exit(0);
			}
			String pedigreefiles_folder =null;
			String crossid =null;
			String RILpheno =null;
			String geno250k =null;
			String output_prefix=null;
			int block_size=5000, the_phe_index=-1;

			double p_after_corr=1000;
			double maf_threshold_plot=0.05;
			int round=2;
			for(int k=0; k<args.length; k++){
				if(args[k].startsWith("-")){
					if(args[k].equals("-p")) pedigreefiles_folder=args[k+1];
					else if(args[k].equals("-c"))  crossid=args[k+1];
					else if(args[k].equals("-ril")) RILpheno=args[k+1];
					else if(args[k].equals("-fg")) geno250k=args[k+1];
					else if(args[k].equals("-o"))  output_prefix=args[k+1];
					else if(args[k].equals("-b"))  block_size=Integer.parseInt(args[k+1]);
					else if(args[k].equals("-index"))the_phe_index=Integer.parseInt(args[k+1]);
					else if(args[k].equals("-p"))p_after_corr=Double.parseDouble(args[k+1]);
						//TODO: "-p" argument used twice. Better to avoid repetition?
					else if(args[k].equals("-r"))round=Integer.parseInt(args[k+1]);
					else if(args[k].equals("-maf"))maf_threshold_plot=Double.parseDouble(args[k+1]);
				}
			}if(pedigreefiles_folder==null||crossid==null||RILpheno==null||geno250k==null||output_prefix==null){
				System.out.println("Input or output paths can't be empty!");
				System.exit(0);
			}
			VariantsDouble.run_nam_imputation(pedigreefiles_folder, crossid, RILpheno, geno250k, output_prefix, block_size);
			LM.stepwise_reg_nam(output_prefix+".geno.num.csv.hdf5", output_prefix+".pheno.tsv", the_phe_index,
					output_prefix, p_after_corr, maf_threshold_plot, round);
		}else if(function.equals("hdf52csv")){
			System.out.println("Convert HDF5 genotype file to readable CSV format"+"\n");
			if(args.length==1){
				System.out.println("Usage:\n\t" +
						"<-ig\tinput_HDF5_file>\n\t" +
						"<-o\toutput_CSV_file>\n");
				System.exit(0);
			}else{
				String input_geno=null, output_csv=null;
				for(int k=1; k<args.length; k++){
					if(args[k].startsWith("-")){
						if(args[k].equals("-ig"))input_geno=args[k+1];
						else if(args[k].equals("-o"))output_csv=args[k+1];
						else System.out.println("Option "+args[k]+" is undefined.");
					}
				}if(input_geno==null||output_csv==null){
					System.out.println("Input or output can't be null!");
				}else{
					VariantsDouble genotype=new VariantsDouble(input_geno);
					genotype.output2csv(output_csv);
				}
			}
		} else{
			System.out.println("Did you make a typo? \""+function+"\" is not a supported function. Please try again.");
			System.exit(0);
		}

	}

}
