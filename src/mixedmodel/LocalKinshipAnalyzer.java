package mixedmodel;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;


import myMathLib.Test;
import myPlotLab.MyHeatMap;
import myPlotLab.MyHistogram;
import myPlotLab.MyManhattanPlot;


public class LocalKinshipAnalyzer {
	VariantsDouble variants;
	int tiling_win_size;
	int[][][] gene_coordinates; // chr x num_of_genes_per_chr x 2[start, end]
	String[][] gene_ids; 	// chr x num_of_genes_per_chr
	double[][] global_kinship;
	//	public static int[] lengths_of_chrs={30427671,19698289,23459830,18585056,26975502};
	
	/*
	 * Constructor.
	 * the "gene_file" actually can be any file specifying interested regions in the format:
	 * id\chr\tstart\itend\n
	 * 
	 * the chr has to be numbered as "1","2,"... instead of strings like "Chr1" etc.
	 */

	public LocalKinshipAnalyzer(String hdf5_file, int win_size, String gene_file) {
		try {
			this.variants = new VariantsDouble(hdf5_file);
			this.tiling_win_size = win_size;
			if (gene_file != null) {
				BufferedReader br = new BufferedReader(new FileReader(gene_file));
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	/*
	 * This is just a template of doing anything that needs to go through the tiling windows.
	 */
	public void go_through_tiling_windows() {
		try {
			//BufferedWriter bw = new BufferedWriter(new FileWriter(outfile));
			for (int chr = 0; chr < variants.num_chrs; chr++) {
				int last_position = variants.locations[chr][variants.num_sites[chr] - 1];
				int start = 1, end = tiling_win_size;
				while (start < last_position) {
					double[][] data = variants.load_variants_in_region(chr, start, end);
					for (int k = 0; k < data.length; k++) {
						//do something
					}
					start = start + tiling_win_size / 2;
					end = start + tiling_win_size - 1;
				}
			}//bw.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	/*
	 * calculate local-kinship 
	 * scale = the biggest difference between genotypes.
	 */
	// TODO: Use of tiled window kinship matrices?
	public void calculating_kinship_tiling_windows(String output_folder, double scale) {
		System.out.println("Local_kinship for win-size= " + tiling_win_size + ".");
		try {
			for (int chr = 0; chr < variants.num_chrs; chr++) {
				int last_position = variants.locations[chr][variants.num_sites[chr] - 1];
				int start = 1, end = start + tiling_win_size - 1;
				while (start < last_position) {
					double[][] data = variants.load_variants_in_region(chr, start, end);
					String the_kin_file = output_folder + "kinship." + variants.sample_size + ".chr"
							+ (chr + 1) + "." + start + "." + end;
					BufferedWriter bw = new BufferedWriter(new FileWriter(the_kin_file + ".raw.ibs"));
					double[][] kinship = new double[variants.sample_size][variants.sample_size];
					for (int i = 0; i < variants.sample_size; i++) {
						kinship[i][i] = 1;
						for (int j = i + 1; j < variants.sample_size; j++) {
							for (int k = 0; k < data.length; k++) {
								kinship[i][j] += (scale - Math.abs(data[k][i] - data[k][j]));
							}
							kinship[i][j] = kinship[i][j] / scale / data.length;
							kinship[j][i] = kinship[i][j];
						}
					}

					for (int i = 0; i < variants.sample_size - 1; i++) {
						bw.write(variants.sample_ids[i] + ",");
					}
					bw.write(variants.sample_ids[variants.sample_size - 1] + "\n");

					for (int i = 0; i < variants.sample_size; i++) {
						for (int j = 0; j < variants.sample_size - 1; j++) {
							bw.write(kinship[i][j] + ",");
						}
						bw.write(kinship[i][variants.sample_size - 1] + "\n");
					}
					bw.close();

					KinshipMatrix.re_scale_kinship_matrix(the_kin_file + ".raw.ibs", the_kin_file + ".rescaled.ibs");
					start = start + tiling_win_size / 2;
					end = start + tiling_win_size - 1;
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public void local_VO_wins(Phenotype phenotype, String genotype_hdf5_file, String global_kinship_file, String local_kinship_folder,
			String mlr_output_file, double step, boolean plot, String method, int min_sample_size){
		int the_sample_size=FaSTLMM.sample_size(phenotype, genotype_hdf5_file);
		if(the_sample_size< min_sample_size){
			System.out.println("Overlap of sample size too small: "+the_sample_size);
			return;
		}
		if(method.equals("transform")){
			FaSTLMM_LocalKinship calculator=new FaSTLMM_LocalKinship(phenotype, genotype_hdf5_file, global_kinship_file);
			calculator.analysis_tiling_win(mlr_output_file, tiling_win_size, plot);
		}//else if(method.equals("transform")){
		//	local_VO_tranform(phenotype, genotype_hdf5_file, global_kinship_file, mlr_output_file, plot);
		//}
		else if(method.equals("full")){
			local_VO_all_grids(phenotype, genotype_hdf5_file, global_kinship_file, local_kinship_folder, 
					mlr_output_file, step, plot);
		}else if(method.equals("grid")){
			local_VO_best_grids(phenotype, genotype_hdf5_file, global_kinship_file, local_kinship_folder, 
					mlr_output_file, plot);
		}
	}
	
	public void local_VO_regions(Phenotype phenotype, String genotype_hdf5_file, String global_kinship_file, String regions_file,
			String mlr_output_file, boolean plot, int min_sample_size){
		int the_sample_size=FaSTLMM.sample_size(phenotype, genotype_hdf5_file);
		if(the_sample_size< min_sample_size){
			System.out.println("Overlap of sample size too small: "+the_sample_size);
			return;
		}
		Regions the_regions=new Regions(regions_file, variants); 
		FaSTLMM_LocalKinship calculator=new FaSTLMM_LocalKinship(phenotype, genotype_hdf5_file, global_kinship_file);
		calculator.analysis_specified_region(mlr_output_file, the_regions.region_coordinates, plot);
	}
	
	/*
	 * Check all points from 0.01 to 0.99 (if step=0.01 as default). Results: (1) find the best alpha, and (2) the distribution.
	 * It has to estimate the variance components 100 times (if step =0.01, as default).
	 */
	public void local_VO_all_grids(Phenotype phenotype, String genotype_hdf5_file, String global_kinship_file, String local_kinship_folder,
			String mlr_output_file, double step, boolean plot){
		try{
			VariantsDouble genotype=new VariantsDouble(genotype_hdf5_file);
			System.out.println(phenotype.phe_id);
			
			EMMA emma_global=new EMMA(phenotype, genotype, new KinshipMatrix(global_kinship_file).kinship_matrix.getData());//  genotype.read_in_kinship(global_kinship_file));
//			emma_global.phenotype.transform();
			emma_global.REMLE_null();
			BufferedWriter bw=new BufferedWriter(new FileWriter(mlr_output_file));
			BufferedWriter bw_dist=new BufferedWriter(new FileWriter(mlr_output_file+".dist.csv"));
			bw.write("#SampleSize="+emma_global.sample_size+ "; REML="+emma_global.remle_REML+"; " +
					"h0-heritability="+emma_global.heritability
					+"; Transform="+emma_global.phenotype.transform_approach+"(P="+emma_global.phenotype.normality_pvaule+")"+"\n");
			bw_dist.write("#SampleSize="+emma_global.sample_size+ "; REML="+emma_global.remle_REML+"; " +
					"h0-heritability="+emma_global.heritability
					+"; Transform="+emma_global.phenotype.transform_approach+"(P="+emma_global.phenotype.normality_pvaule+")"+"\n");
			bw.write("chr, position, p-value, variance_local, variance_global, h1_heritabilities\n");
			bw_dist.write("chr, position");
			for(double alpha=0;alpha<1;alpha+=step){bw_dist.write(", "+alpha);}
			bw_dist.write("\n");
			for(int chr=0;chr<variants.num_chrs;chr++){
				int last_position=variants.locations[chr][variants.num_sites[chr]-1];
				int start=1, end=start+tiling_win_size-1;
				while(start<last_position){
					bw_dist.write((chr+1)+", "+start);
	//				double[][] data=variants.load_variants_in_region(chr, start, end);
					String the_kin_file=local_kinship_folder+"kinship."+variants.sample_size+".chr"
							+(chr+1)+"."+start+"."+end+".rescaled.ibs";
					EMMA emma_local=new EMMA(emma_global.phenotype, genotype, new KinshipMatrix(the_kin_file).kinship_matrix.getData());//genotype.read_in_kinship(the_kin_file));	
					EMMA emma_new=new EMMA(emma_global.phenotype, genotype,  new KinshipMatrix(the_kin_file).kinship_matrix.getData());//genotype.read_in_kinship(the_kin_file));		
					if(emma_global.sample_size!=emma_local.sample_size){
						System.out.println("emma_global.sample_size!=emma_local.sample_size");
					}
					double[] best_h={1,0,emma_global.heritability}; // p-value, h_local, h_global
					for(double alpha=0;alpha<1;alpha+=step){
						double[][] new_kinship=new double[emma_global.sample_size][emma_global.sample_size];
						for(int i=0;i<emma_global.sample_size;i++){
							for(int j=0;j<emma_global.sample_size;j++){
								new_kinship[i][j]=emma_local.kinship_matrix[i][j]*alpha+emma_global.kinship_matrix[i][j]*(1-alpha);
							}
						}
						emma_new.clear();
						emma_new.kinship_matrix=KinshipMatrix.re_scale_kinship_matrix(new_kinship);
//						check_kinship(emma_new.kinship_matrix);
//						System.out.println("\nwin"+start+":alpha;"+alpha);
						emma_new.REMLE_null();
						bw_dist.write(", "+emma_new.remle_REML+"/"+emma_new.heritability);
						if(!Double.isNaN(emma_new.remle_REML)){
							double p_value=Test.chi2pr(-2.0*(emma_global.remle_REML-emma_new.remle_REML), 1);
							if(p_value<best_h[0]){
								best_h[0]=p_value;
								best_h[1]=emma_new.heritability*alpha;
								best_h[2]=emma_new.heritability*(1-alpha);
							}
						}
					}					
					bw.write((chr+1)+", "+start+", "+best_h[0]+", "+best_h[1]+", "+best_h[2]+", "+(best_h[1]+best_h[2])+"\n");
					bw_dist.write("\n");
					System.out.println("finsihed win"+start);
//					bw.flush();
					start=start+tiling_win_size/2;
					end=start+tiling_win_size-1;
				}
			}bw.close();
			bw_dist.close();
			if(plot)make_plot_one_phen(mlr_output_file, mlr_output_file+".dist.csv");
		}catch(Exception e){e.printStackTrace();}
	}
	
	/*
	 * without Bayesian analysis so that no need to run for all grids. 
	 * It will check 0.1, 0.2, ... , 1. and select one window to run another 10 grids to find the best estimate. 
	 */
	public void local_VO_best_grids(Phenotype phenotype, String genotype_hdf5_file, String global_kinship_file, String local_kinship_folder,
			String mlr_output_file, boolean plot){
		final double first_step=0.1, second_step=0.01;
		try{
			VariantsDouble genotype=new VariantsDouble(genotype_hdf5_file);
			System.out.println(phenotype.phe_id);
			
			EMMA emma_global=new EMMA(phenotype, genotype, new KinshipMatrix(global_kinship_file).kinship_matrix.getData());//genotype.read_in_kinship(global_kinship_file));
//			emma_global.phenotype.transform();
			emma_global.REMLE_null();
			BufferedWriter bw=new BufferedWriter(new FileWriter(mlr_output_file));
			
			bw.write("#SampleSize="+emma_global.sample_size+ "; REML="+emma_global.remle_REML+"; " +
					"h0-heritability="+emma_global.heritability
					+"; Transform="+emma_global.phenotype.transform_approach+"(P="+emma_global.phenotype.normality_pvaule+")"+"\n");
			
			bw.write("chr, position, p-value, variance_local, variance_global, h1_heritabilities\n");
			for(int chr=0;chr<variants.num_chrs;chr++){
				int last_position=variants.locations[chr][variants.num_sites[chr]-1];
				int start=1, end=start+tiling_win_size-1;
				while(start<last_position){
	//				double[][] data=variants.load_variants_in_region(chr, start, end);
					String the_kin_file=local_kinship_folder+"kinship."+variants.sample_size+".chr"
							+(chr+1)+"."+start+"."+end+".rescaled.ibs";
					EMMA emma_local=new EMMA(emma_global.phenotype, genotype, new KinshipMatrix(the_kin_file).kinship_matrix.getData());//genotype.read_in_kinship(the_kin_file));	
					EMMA emma_new=new EMMA(emma_global.phenotype, genotype, new KinshipMatrix(the_kin_file).kinship_matrix.getData());//genotype.read_in_kinship(the_kin_file));		
					if(emma_global.sample_size!=emma_local.sample_size){
						System.out.println("emma_global.sample_size!=emma_local.sample_size");
					}
					double[] best_h={1,0,emma_global.heritability}; // p-value, h_local, h_global
					double best_first_alpha=0;
					// first search
					for(double alpha=0;alpha<1;alpha+=first_step){
						double[][] new_kinship=new double[emma_global.sample_size][emma_global.sample_size];
						for(int i=0;i<emma_global.sample_size;i++){
							for(int j=0;j<emma_global.sample_size;j++){
								new_kinship[i][j]=emma_local.kinship_matrix[i][j]*alpha+emma_global.kinship_matrix[i][j]*(1-alpha);
							}
						}
						emma_new.clear();
						emma_new.kinship_matrix=KinshipMatrix.re_scale_kinship_matrix(new_kinship);
						emma_new.REMLE_null();
						if(!Double.isNaN(emma_new.remle_REML)){
							double p_value=Test.chi2pr(-2.0*(emma_global.remle_REML-emma_new.remle_REML), 1);
							if(p_value<best_h[0]){
								best_first_alpha=alpha;
								best_h[0]=p_value;
								best_h[1]=emma_new.heritability*alpha;
								best_h[2]=emma_new.heritability*(1-alpha);
							}
						}
					}	
					// search again around the best hit
					double search_start=best_first_alpha-second_step*5;
					if(search_start<=0)search_start=0;
					double search_end=best_first_alpha+second_step*5;
					if(search_end>=1)search_end=1;
					for(double alpha=search_start;alpha<=search_end;alpha+=second_step){
						double[][] new_kinship=new double[emma_global.sample_size][emma_global.sample_size];
						for(int i=0;i<emma_global.sample_size;i++){
							for(int j=0;j<emma_global.sample_size;j++){
								new_kinship[i][j]=emma_local.kinship_matrix[i][j]*alpha+emma_global.kinship_matrix[i][j]*(1-alpha);
							}
						}
						emma_new.clear();
						emma_new.kinship_matrix=KinshipMatrix.re_scale_kinship_matrix(new_kinship);
						emma_new.REMLE_null();
						if(!Double.isNaN(emma_new.remle_REML)){
							double p_value=Test.chi2pr(-2.0*(emma_global.remle_REML-emma_new.remle_REML), 1);
							if(p_value<best_h[0]){
								best_h[0]=p_value;
								best_h[1]=emma_new.heritability*alpha;
								best_h[2]=emma_new.heritability*(1-alpha);
							}
						}
					}
					bw.write((chr+1)+", "+start+", "+best_h[0]+", "+best_h[1]+", "+best_h[2]+", "+(best_h[1]+best_h[2])+"\n");
					System.out.println("finsihed win"+start);
//					bw.flush();
					start=start+tiling_win_size/2;
					end=start+tiling_win_size-1;
				}
			}bw.close();
			if(plot)make_plot_one_phen(mlr_output_file);
		}catch(Exception e){e.printStackTrace();}
	}
	
	public static double[] generate_breaks(int[] lengths_of_chrs){
		double[] breaks=new double[lengths_of_chrs.length];
		breaks[0]=0;
		for(int chr=1;chr<lengths_of_chrs.length;chr++)breaks[chr]=breaks[chr-1]+lengths_of_chrs[chr-1];
		return breaks;
	}
	/*
	 * This function makes five plots immediately after calling 
	 * public void local_VO(Phenotype phenotype, String genotype_hdf5_file, String global_kinship_file, String local_kinship_folder,
	 *		String output_file, double step)
	 *	(1) Manhattan plot for p-value of likelihood ratio results.
	 *	(2) Manhattan plot for R2 of likelihood ratio results.
	 *	(3) Heatmap for Bayesian result
	 *	(4) Histogram for regional VC distribution
	 *	(5) Histogram for total variance explained by chromosome, and chr-length as the control		  
	 */
	public static void make_plot_one_phen(String mlr_file, String bayesian_file){
		
		String phe_name=mlr_file.split("/")[mlr_file.split("/").length-1];
		try{
			if((new File(mlr_file)).exists()){
				AssociationResults mlr=new AssociationResults(mlr_file, 2);
				double[] breaks=generate_breaks(mlr.length_of_chrs);//generate_breaks(lengths_of_chrs);
				double[] locations=new double[mlr.location.length];
				double[] pvalues=new double[mlr.location.length];
				double[] var=new double[mlr.location.length];
				for(int i=0;i<locations.length;i++){
					locations[i]=(mlr.location[i]+breaks[mlr.chr[i]-1])/1000000.0;		
					pvalues[i]=-Math.log10(mlr.pvalue[i]);
					var[i]=mlr.AdjustedR2[i]*100.0;
				}	
				MyManhattanPlot plot=new MyManhattanPlot(phe_name, "Locations (Mb)", "-log10(p-value)", mlr.chr,
						locations, pvalues, 2000, 400, mlr_file+".pvalue.png");
				plot=new MyManhattanPlot(phe_name, "Locations (Mb)", "Variance Explained (%)", mlr.chr,
						locations, var, 2000, 400, mlr_file+".VarR2.png");
			}else{
				System.out.println("NOFILE: "+mlr_file);
			}
			if((new File(bayesian_file)).exists()){
				RegionalDistribution dist=new RegionalDistribution(bayesian_file);
				MyHeatMap heatmap=new MyHeatMap(dist.distribution, phe_name, "Locations (Mb)", 
						"Variance (%)", 120, 100, bayesian_file+".overall.png");
				MyHistogram hist=new MyHistogram(dist.regional_contribution, dist.regional_category,phe_name, "Variance Contributed (%)", "Proportion of Regions (%)",
						700, 400, bayesian_file+".regional.png");
				hist=new MyHistogram(dist.chr_contribution, dist.chr_names, (new String[] {"Variance Component Explained", "Chromosome Length"}), phe_name, null, "Proportion", 
						700, 400, bayesian_file+".chr.png");
			}else{
				System.out.println("NOFILE: "+bayesian_file);
			}
			
		}catch(Exception e){e.printStackTrace();}
	}
	
	public static void make_plot_one_phen(String mlr_file){		
		String phe_name=mlr_file.split("/")[mlr_file.split("/").length-1];
		try{
			if((new File(mlr_file)).exists()){
				AssociationResults mlr=new AssociationResults(mlr_file, 2);
				double[] breaks=generate_breaks(mlr.length_of_chrs);//lengths_of_chrs);
				double[] locations=new double[mlr.location.length];
				double[] pvalues=new double[mlr.location.length];
				double[] var=new double[mlr.location.length];
				for(int i=0;i<locations.length;i++){
					locations[i]=(mlr.location[i]+breaks[mlr.chr[i]-1])/1000000.0;		
					pvalues[i]=-Math.log10(mlr.pvalue[i]);
					var[i]=mlr.AdjustedR2[i]*100.0;
				}	
				MyManhattanPlot plot=new MyManhattanPlot(phe_name, "Locations (Mb)", "-log10(p-value)", mlr.chr,
						locations, pvalues, 2000, 400, mlr_file+".pvalue.png");
				plot=new MyManhattanPlot(phe_name, "Locations (Mb)", "Variance Explained (%)", mlr.chr,
						locations, var, 2000, 400, mlr_file+".VarR2.png");
			}else{
				System.out.println("FILE not exists: "+mlr_file);
			}
		}catch(Exception e){e.printStackTrace();}
	}
	
	public static void make_plot_multiple_csv(String folder){		
		
		try{
			File[] files=(new File(folder).listFiles());
			for(int file_index=0;file_index<files.length;file_index++){
				if(files[file_index].toString().endsWith(".csv")){
					String mlr_file=files[file_index].toString();
					String phe_name=mlr_file.split("/")[mlr_file.split("/").length-1];
					if((new File(mlr_file)).exists()){
						AssociationResults mlr=new AssociationResults(mlr_file, 2);
						double[] breaks=generate_breaks(mlr.length_of_chrs);//lengths_of_chrs);
						double[] locations=new double[mlr.location.length];
						double[] pvalues=new double[mlr.location.length];
						double[] var=new double[mlr.location.length];
						for(int i=0;i<locations.length;i++){
							locations[i]=(mlr.location[i]+breaks[mlr.chr[i]-1])/1000000.0;		
							pvalues[i]=-Math.log10(mlr.pvalue[i]);
							var[i]=mlr.AdjustedR2[i]*100.0;
						}	
						MyManhattanPlot plot=new MyManhattanPlot(phe_name, "Locations (Mb)", "-log10(p-value)", mlr.chr,
								locations, pvalues, 2000, 400, mlr_file+".pvalue.png");
						plot=new MyManhattanPlot(phe_name, "Locations (Mb)", "Variance Explained (%)", mlr.chr,
								locations, var, 2000, 400, mlr_file+".VarR2.png");
					}else{
						System.out.println("FILE not exists: "+mlr_file);
					}
				}
			}
			
		}catch(Exception e){e.printStackTrace();}
	}
	
	
	public static void check_kinship(double[][] kinship_matrix){
		for(int i=0;i<kinship_matrix.length;i++){
			for(int j=i+1;j<kinship_matrix.length;j++){
				if(Double.compare(kinship_matrix[i][j],kinship_matrix[j][i])!=0){
					System.out.println("kinship wrong, i,j");
				}
			}
		}
	}
	
	
	/*
	 * Transform based on global kinship and estimate the effect of local kinship by EMMA. 
	 * It only estimates the variance component once, therefore very fast.
	 */
	public void local_VO_tranform(Phenotype phenotype, String genotype_hdf5_file, String global_kinship_file, 
			String mlr_output_file, boolean plot){
		DescriptiveStatistics calculator=new DescriptiveStatistics();  
		try{
			VariantsDouble genotype=new VariantsDouble(genotype_hdf5_file);
			System.out.println(phenotype.phe_id);
			
			EMMA emma_global=new EMMA(phenotype, genotype,  new KinshipMatrix(global_kinship_file).kinship_matrix.getData());//genotype.read_in_kinship(global_kinship_file));
			emma_global.REMLE_null();
			double[][] decomposed_array = emma_global.reml_decompositioned_array();	
			emma_global.phenotype.generate_new_Y_by_multiplying(decomposed_array);
			// use the transformed Y to replace the original Y for Local-analysis.
			emma_global.phenotype.values=emma_global.phenotype.new_Y_4emmax.clone();
			// remove rescaled beta0  
			double[] scaled_beta0=EMMAX.transform_beta0(decomposed_array);
			double mean_y=(new DescriptiveStatistics(emma_global.phenotype.values)).getMean();
			double mean_beta0=(new DescriptiveStatistics(scaled_beta0)).getMean();
			double ratio=mean_y/mean_beta0;
			for(int k=0;k<decomposed_array.length;k++){
				emma_global.phenotype.values[k]-=(ratio*scaled_beta0[k]);
			}
			BufferedWriter bw=new BufferedWriter(new FileWriter(mlr_output_file));
			bw.write("#SampleSize="+emma_global.sample_size+ "; REML="+emma_global.remle_REML+"; " +
					"h0-heritability="+emma_global.heritability
					+"; Transform="+emma_global.phenotype.transform_approach+"(P="+emma_global.phenotype.normality_pvaule+")"+"\n");
			bw.write("chr, position, p-value, variance_local, variance_global, h1_heritabilities\n");			
			for(int chr=0;chr<variants.num_chrs;chr++){
				int last_position=variants.locations[chr][variants.num_sites[chr]-1];
				int start=1, end=start+tiling_win_size-1;
				while(start<last_position){					
					double[][] data=variants.load_variants_in_region(chr, start, end);					
					double[][] transformed_data=EMMAX.transform_data(data, decomposed_array, emma_global);
					EMMA emma_local=new EMMA(emma_global.phenotype, genotype, 
							VariantsDouble.calculate_RRM_local_kinship(transformed_data),true);	
					if(emma_global.sample_size!=emma_local.sample_size){
						System.out.println("emma_global.sample_size!=emma_local.sample_size");
					}					
					emma_local.REMLE_null();
					double p_value=Test.chi2pr(-2.0*(emma_local.remle_REML-myMathLib.StatFuncs.null_logL(emma_local.phenotype.values)), 1);			
					bw.write((chr+1)+", "+start+", "+p_value+", "+emma_local.remle_vg/(emma_local.remle_vg+emma_local.remle_ve)
							+", "+0+", "+emma_global.heritability+"\n");
//					System.out.println("finsihed win"+start);
//					bw.flush();
					start=start+tiling_win_size/2;
					end=start+tiling_win_size-1;
				}System.out.println("Transformed-EMMAX, finished chr"+(chr+1));
			}bw.close();
			if(plot)make_plot_one_phen(mlr_output_file);
		}catch(Exception e){e.printStackTrace();}
	}
	
	
//	public void local_VO_tranform_fastlmm(Phenotype phenotype, String genotype_hdf5_file, String global_kinship_file, 
//			String mlr_output_file, boolean plot){
//		try{
//			VariantsDouble genotype=new VariantsDouble(genotype_hdf5_file);
//			System.out.println(phenotype.phe_id);
//			
//			
//			
//			EMMA emma_global=new EMMA(phenotype, genotype, genotype.read_in_kinship(global_kinship_file));
//			emma_global.REMLE_null();
//			double[][] decomposed_array = emma_global.reml_decompositioned_array();	
//			emma_global.phenotype.generate_new_Y_by_multiplying(decomposed_array);
//			// use the transformed Y to replace the original Y for Local-analysis.
//			emma_global.phenotype.values=emma_global.phenotype.new_Y_4emmax.clone();
//			
//			BufferedWriter bw=new BufferedWriter(new FileWriter(mlr_output_file));
//			bw.write("#SampleSize="+emma_global.sample_size+ "; REML="+emma_global.remle_REML+"; " +
//					"h0-heritability="+emma_global.heritability +"\n");
//			bw.write("chr, position, p-value, variance_local, variance_global, h1_heritabilities\n");			
//			for(int chr=0;chr<variants.num_chrs;chr++){
//				int last_position=variants.locations[chr][variants.num_sites[chr]-1];
//				int start=1, end=start+tiling_win_size-1;
//				while(start<last_position){					
//					double[][] data=variants.load_variants_in_region(chr, start, end);	
//					double[][] transformed_data=EMMAX.transform_data(data, decomposed_array, emma_global);
//					EMMA emma_local=new EMMA(emma_global.phenotype, genotype, 
//							VariantsDouble.calculate_RRM_local_kinship(transformed_data),true);	
//					if(emma_global.sample_size!=emma_local.sample_size){
//						System.out.println("emma_global.sample_size!=emma_local.sample_size");
//					}					
//					emma_local.REMLE_null();
//					double p_value=Test.chi2pr(-2.0*(emma_local.remle_REML-myMathLib.StatFuncs.null_logL(emma_local.phenotype.values)), 1);			
//					bw.write((chr+1)+", "+start+", "+p_value+", "+emma_local.remle_vg/(emma_local.remle_vg+emma_local.remle_ve)
//							+", "+0+", "+emma_global.heritability+"\n");
//					System.out.println("finsihed win"+start);
////					bw.flush();
//					start=start+tiling_win_size/2;
//					end=start+tiling_win_size-1;
//				}
//			}bw.close();
//			if(plot)make_plot_one_phen(mlr_output_file);
//		}catch(Exception e){e.printStackTrace();}
//	}
	
	/*
	 * OLD program trying whether the output are the same to input
	 */
//	public void go_through_windows_no_tiling(String outfile){
//		try{
//			BufferedWriter bw=new BufferedWriter(new FileWriter(outfile));
//			for(int chr=0;chr<variants.num_chrs;chr++){
//				int last_position=variants.locations[chr][variants.num_sites[chr]-1];
//				int start=1, end=tiling_win_size;
//				while(start<last_position){
//					double[][] data=variants.load_variants_in_region(chr, start, end);
//					start=end+1;end=start+tiling_win_size-1;
//					for(int k=0;k<data.length;k++){
//					for(int i=0;i<data[k].length;i++){
//						bw.write(data[k][i]+",");
//						}bw.write("\n");
//					}
//				}
//			}bw.close();
//		}catch(Exception e){e.printStackTrace();}		
//	}
	
	
	
	public static void main(String[] args) {
		long start=System.currentTimeMillis();
//		String genotype_csv_file="/Users/quan.long/Documents/Projects/1001g/GWAS/data/data180.gwas.1.csv";
//		String phe_file="/Users/quan.long/Documents/Projects/1001g/GWAS/data/phen_raw_092910.tsv";
//		String global_kinship_file_ori="/Users/quan.long/Documents/Projects/1001g/GWAS/data/data180.gwas.1.csv.kinship";
		
		String global_kinship_file="/Volumes/Projects/DATA/arabidopsis_genomes/imputed_swe259.csv.num.hdf5.kinkinship.rescaled.ibs";
		
	
		String genotype_hdf5_file="/Volumes/Projects/DATA/arabidopsis_genomes/imputed_swe259.csv.num.hdf5";
		String genome_size="/Volumes/Projects/Local-kinship/genome_size/genome_size.tsv";
		String local_kinship_folder="/Volumes/Projects/DATA/arabidopsis_genomes/kinship_files_imputed_swe259/100K/";
		String local_kinship_VO_folder="/Users/quan.long/Documents/Projects/GeneMappingMethodology/LocalGlobal/check_programs/local_kinship_VO/";
		
		String to_plot_folder_emma="/Volumes/Projects/Local-kinship/genome_size/result_Jan11_trans/";
		String to_plot_folder_fastlmm="/Volumes/Projects/Local-kinship/genome_size/try_25Feb/";
		
		
		int win_size=25000;
		double step=0.01;
		
//		MultiPhenotype phenotypeS=new MultiPhenotype(genome_size);
//		LocalKinshipAnalyzer local_k=new LocalKinshipAnalyzer(genotype_hdf5_file, win_size, null);
//		for(int i=0;i<phenotypeS.num_of_pheno;i++){
//			Phenotype phenotype=phenotypeS.phenotypes[i];
//			local_k.local_VO(phenotype, genotype_hdf5_file,  global_kinship_file, local_kinship_folder,
//					to_plot_folder_fastlmm+"Local_VO."+i+"."+phenotype.phe_id+".w"+win_size+".csv", step, true, "transform",0);
//		}
		
		
		
		
		
//		LocalKinshipAnalyzer local_k=new LocalKinshipAnalyzer(hdf5_file, win_size, null);
//		double scale=2.0;
//		local_k.calculating_kinship_tiling_windows(local_kinship_folder, scale);
//		
//		MultiPhenotype phenotypeS=new MultiPhenotype(phe_file);
//		for(int phe_index=4;phe_index<phenotypeS.num_of_pheno;phe_index++){
//			String out_phe_file=local_kinship_VO_folder+"Local_VO."+phe_index+"_"+phenotypeS.phenotypes[phe_index].phe_id+".csv";
//			local_k.local_VO(phenotypeS.phenotypes[phe_index], hdf5_file, global_kinship_file_rescaled, local_kinship_folder,
//					out_phe_file, 0.01);
//		}
//		System.out.println((System.currentTimeMillis()-start)/1000);

		if(args.length<1){
			System.out.println("Plot .csv files: Usage <folder_path>");
			System.exit(0);
		}
		make_plot_multiple_csv(args[0]);

	}
}
