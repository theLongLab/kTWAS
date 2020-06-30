package Region_Emmax;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;


import mixedmodel.EMMAX;

public class Region_EMMAX {
	//regions_file contains
	//1;1000;2000
	//2;4000;5000
	public static int[][][] read_regions_from_file(String regions_file){
		int [][][] regions;
		ArrayList<int[]> regions_listArrayList = new ArrayList<int[]>();
		try {
			BufferedReader bReader = new BufferedReader(new FileReader(regions_file));
			String lineString = bReader.readLine();
			while(lineString != null) {
				if(!lineString.startsWith("#")) {
					String[] lineStrings = lineString.split("\t");
					int[] lineIntegers = new int[lineStrings.length-1];
					for(int i=1; i<lineStrings.length;i++) {  //the region file contains r1\t 1\t 400 \t 500
						lineIntegers[i-1]=Integer.parseInt(lineStrings[i]);
					}
					regions_listArrayList.add(lineIntegers);
				}
				lineString = bReader.readLine();
			}
			bReader.close();		
		}catch (Exception e) {
			e.printStackTrace();
		}
			
		int current_chr=regions_listArrayList.get(0)[0];
		ArrayList<Integer> chrList = new ArrayList<Integer>();
		chrList.add(current_chr-1);//the chr starts from 0 !!!
		int count_chr_regions=1;
		ArrayList<Integer> chr_regionsList = new ArrayList<Integer>();
		for(int i=1; i<regions_listArrayList.size();i++) {
			if(current_chr==regions_listArrayList.get(i)[0]) {
				count_chr_regions+=1;
			}else { //different chr
				chr_regionsList.add(count_chr_regions);
				current_chr=regions_listArrayList.get(i)[0];
				chrList.add(current_chr-1);//the chr starts from 0 !!!
				count_chr_regions=1;
			}
		}chr_regionsList.add(count_chr_regions);//add the number of regions for the last chr
		if(chr_regionsList.size()!=chrList.size()) {
			System.out.println("Errors in loading the regions number for each chromosome");
			System.exit(0);
		}
		int index_4regions_list = 0;
		regions= new int[chrList.size()][][];
		for(int chr_index=0; chr_index<chrList.size(); chr_index++) {
			int num_regions_in_current_chr = chr_regionsList.get(chr_index);
			regions[chr_index] = new int[num_regions_in_current_chr][3];
			for(int chr_region_index=0;chr_region_index<num_regions_in_current_chr;chr_region_index++) {
				regions[chr_index][chr_region_index][0]=chrList.get(chr_index); //chr of the this region
				regions[chr_index][chr_region_index][1]=regions_listArrayList.get(index_4regions_list)[1]; //start locus of this region
				regions[chr_index][chr_region_index][2]=regions_listArrayList.get(index_4regions_list)[2]; //end locus of this region
				index_4regions_list++;
			}
		}
		if(regions_listArrayList.size()==index_4regions_list) {
			System.out.println("regions loaded successfully");
		}else {
			System.out.println("Errors in loading regions");
			System.exit(0);
		}
		return regions;
	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		if(args.length==0) {
			System.out.println("Input files: \n"
					+ "1)genotype file\n2)phenotype file\n3)kinship\n4)output folder\n5)regions file\n6)phenotype index");
		}else {
			String genotype_hdf5_file=args[0]; 
			String phe_file=args[1];
			String kinship_file=args[2];
			String output_folder=args[3];
			String regions_file=args[4];
			int phe_index=Integer.parseInt(args[5]);
			double p_value_after_correction=1000.0; 
			int min_sample_size=40;
			double maf_plot_threshold=0.05;
			boolean plot=false;
			int[][][] regions; 
			
//			String regions_file="C:/Users/Administrator/Desktop/regions_test.txt";
			regions = read_regions_from_file(regions_file);
			EMMAX.emmax_analysis_regions(genotype_hdf5_file, phe_file, kinship_file, output_folder, 
					p_value_after_correction, min_sample_size, phe_index, maf_plot_threshold, plot, regions);
		}
	}


}
