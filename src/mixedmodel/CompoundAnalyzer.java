package mixedmodel;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;


import myMathLib.Test;
import myPlotLab.MyHeatMap;
import myPlotLab.MyHistogram;
import myPlotLab.MyManhattanPlot;


public class CompoundAnalyzer {
	public VariantsDouble variants;
	public int num_compounds;
	public int[][][] compound_coordinates; // num_compound x num_of_regions_per_compound x 3[chr, start, end]
	public String[] compound_ids; 	// chr x num_of_genes_per_chr
	//double[][] global_kinship;
	
	/*
	 * Constructor.
	 * the "compound_coordinates_file" actually can be any file specifying interested regions in the format:
	 * id\tchr;start;end\tchr;start;end\tchr;start;end\n
	 * 
	 * the chr has to be numbered as "1","2,"... instead of strings like "Chr1" etc.
	 */
	public CompoundAnalyzer(String hdf5_file, String compound_coordinates_file){
		ArrayList<String> ids=new ArrayList<String>();
		ArrayList<int[][]> coordinates=new ArrayList<int[][]>();
		try{
			this.variants=new VariantsDouble(hdf5_file);
			BufferedReader br=new BufferedReader(new FileReader(compound_coordinates_file));
			String line=br.readLine();
			while(line!=null){
				String[] tmp=line.split("\t");
				ids.add(tmp[0]);
				int[][] the_coordinate=new int[tmp.length-1][3];
				for(int k=0;k<tmp.length-1;k++){
					String[] locs=tmp[k+1].split(";");
					for(int i=0;i<3;i++)the_coordinate[k][i]=Integer.parseInt(locs[i]);
				}coordinates.add(the_coordinate);
				line=br.readLine();
			}
			this.num_compounds=ids.size();	
			this.compound_coordinates=new int[this.num_compounds][][];
			this.compound_ids=new String[this.num_compounds];
			for(int k=0;k<this.num_compounds;k++){
				this.compound_coordinates[k]=coordinates.get(k);
				this.compound_ids[k]=ids.get(k);
			}
		}catch(Exception e){e.printStackTrace();}
		System.out.println("Finished loading compound file: "+compound_coordinates_file);
	}	
	
	public void compound_VO(Phenotype phenotype, String genotype_hdf5_file, String global_kinship_file, 
			String mlr_output_file, boolean plot, int min_sample_size){
		int the_sample_size=FaSTLMM.sample_size(phenotype, genotype_hdf5_file);
		if(the_sample_size< min_sample_size){
			System.out.println("Overlap of sample size too small: "+the_sample_size);
			return;
		}
		FaSTLMM_Compound calculator=new FaSTLMM_Compound(phenotype, genotype_hdf5_file, global_kinship_file);
		calculator.analysis_specified_compounds(mlr_output_file, this.compound_coordinates, plot);
	}
	
	
	
	
	public static void main(String[] args) {
		long start=System.currentTimeMillis();

		
		System.out.println((System.currentTimeMillis()-start)/1000);

		if(args.length<1){
			System.out.println("Plot .csv files: Usage <folder_path>");
			System.exit(0);
		}
		

	}
}

