package mixedmodel;

import java.util.ArrayList;
import java.util.Arrays;

public class SearchAndCalculateRegionalKinships {
	
	public int[][] regions;
	public double[] var_comp;
	public String[] reginal_full_kinship_files;
	/*
	 * for simualtion only, given a full data, generate the observed and predicted local kinships
	 * 
	 * calculate and write local-kinships or load wrt whether geno_full==null
	 */
	public SearchAndCalculateRegionalKinships(String local_kinship_gwas_file, int win, 
			int peak_width, double loggedp_cutoff, String local_k_files_folder, 
			VariantsDouble geno_full){
		this.select_and_calculate_regions_p_cutoff(local_kinship_gwas_file, win, peak_width, loggedp_cutoff);
		this.reginal_full_kinship_files=new String[this.regions.length];
		for(int r=0;r<this.regions.length;r++){
			this.reginal_full_kinship_files[r]=local_k_files_folder+"chr"+(regions[r][0]+1)+"."+regions[r][1]+"."+regions[r][2]+".K";
		}		
		// calculate kinships 
		if(geno_full!=null){
			for(int r=0;r<regions.length;r++){			
				KinshipMatrix the_local_full=new KinshipMatrix(geno_full, regions[r]);
				the_local_full.write2file(reginal_full_kinship_files[r]);
				
			}		
		}
		
	}
	
	/*
	 * assign regions based on priori knowledge
	 */
	public SearchAndCalculateRegionalKinships(int[][] regions, double[] var_comp, String local_k_files_folder, 
			VariantsDouble geno_full){	
		if(regions.length!=var_comp.length){
			System.out.println("regions.length!=var_comp.length");
			return;
		}
		this.regions=regions.clone();
		this.var_comp=var_comp.clone();		
		this.reginal_full_kinship_files=new String[this.regions.length];
		for(int r=0;r<this.regions.length;r++){
			this.reginal_full_kinship_files[r]=local_k_files_folder+"chr"+(regions[r][0]+1)+"."+regions[r][1]+"."+regions[r][2]+".K";
		}		
		// calculate kinships 
		if(geno_full!=null){
			for(int r=0;r<regions.length;r++){			
				KinshipMatrix the_local_full=new KinshipMatrix(geno_full, regions[r]);
				the_local_full.write2file(reginal_full_kinship_files[r]);
				
			}		
		}		
	}
	
	public void select_and_calculate_regions_p_cutoff(String local_kinship_result_file, int win, int peak_width, 
			double loggedp_cutoff){
		AssociationResults local=new AssociationResults(local_kinship_result_file, 0.9);
		ArrayList<Integer> chr_array=new ArrayList<Integer>();
		ArrayList<Double> p_array=new ArrayList<Double>();
		ArrayList<Double> var_comp_array=new ArrayList<Double>();
		ArrayList<Integer> locations_array=new ArrayList<Integer>();
		for(int index=0;index<local.num_of_var;index++){
			double logged_p=-Math.log10(local.pvalue[index]);
			if(logged_p>=loggedp_cutoff){
				int the_index_in_array=the_index_in_array(local.chr[index]-1, local.location[index], 
						chr_array, locations_array, peak_width);
				if(the_index_in_array==-1){
					p_array.add(logged_p);
					var_comp_array.add(local.AdjustedR2[index]);
					locations_array.add(local.location[index]);
					chr_array.add(local.chr[index]-1);
				}else if(p_array.get(the_index_in_array)<logged_p){
					p_array.set(the_index_in_array, logged_p);
					var_comp_array.set(the_index_in_array, local.AdjustedR2[index]);
					locations_array.set(the_index_in_array, local.location[index]);
					chr_array.set(the_index_in_array, local.chr[index]-1);
				}
				
			}
		}
		int num_regions=p_array.size();
		this.regions=new int[num_regions][3];
		this.var_comp=new double[num_regions];
		
		for(int region_index=0;region_index<num_regions;region_index++){
			regions[region_index][0]=chr_array.get(region_index);
			regions[region_index][1]=locations_array.get(region_index);
			regions[region_index][2]=regions[region_index][1]+win;
			var_comp[region_index]=var_comp_array.get(region_index);//*p_array.get(region_index);
			System.out.println(regions[region_index][0]+":"+regions[region_index][1]+":"+regions[region_index][2]+":"+
					var_comp_array.get(region_index)+"/"+var_comp[region_index]+".");
		}
		System.out.println("Finsiehd Searching Regions.");
	}
	
	public static int the_index_in_array(int chr, int loc, ArrayList<Integer> chr_array, ArrayList<Integer> loc_array, int peak_width){
		int the_index_in_array=-1;		
		for(int k=0;k<loc_array.size();k++){
			if(chr==chr_array.get(k) && Math.abs(loc-loc_array.get(k))<peak_width)
				the_index_in_array=k;
		}return the_index_in_array;
	}
	
	/*
	 * NOT USED ANY MORE!
	 */
//	public void select_and_calculate_regions_all_chr(String local_kinship_result_file, int num_chr, int win){
//		AssociationResults local=new AssociationResults(local_kinship_result_file, 0.9);			
//		this.regions=new int[num_chr][3];
//		double[] ps=new double[num_chr];
//		this.var_comp=new double[num_chr];
//		Arrays.fill(ps, 1); Arrays.fill(this.var_comp, 0); 
//		for(int chr=0;chr<num_chr;chr++){
//			regions[chr][0]=chr;
//		}
//		for(int index=0;index<local.num_of_var;index++){
//			if(local.pvalue[index]<ps[local.chr[index]-1]){
//				ps[local.chr[index]-1]=local.pvalue[index];
//				var_comp[local.chr[index]-1]=local.AdjustedR2[index]*(-Math.log(local.pvalue[index]));
//				regions[local.chr[index]-1][1]=local.location[index];
//				regions[local.chr[index]-1][2]=local.location[index]+win;
//			}
//		}
//		for(int chr=0;chr<num_chr;chr++){
//			System.out.println(regions[chr][0]+":"+regions[chr][1]+":"+regions[chr][2]+":"+var_comp[chr]+":");
//		}
//	}
}
