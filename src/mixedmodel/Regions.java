package mixedmodel;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;

public class Regions {
	
	public int[][][] region_coordinates; //chr x num_region x 2(start, end)
	public String[][] region_names; //chr x num_region 	
	
	public Regions(int[][][] region_coordinates, String[][] region_names){ 
		this.region_coordinates=region_coordinates;
		this.region_names=region_names;
	}
	
	public Regions(String regions_file, VariantsDouble variants){
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
			if(chrs.size()!=variants.num_chrs){
				System.out.println("Number of Chrs covered by the regions is not the same to genotype file.");
				//return;
			}
			this.region_coordinates=new int[variants.num_chrs][][];
			this.region_names=new String[variants.num_chrs][];
			for(int k=0;k<variants.num_chrs;k++){
				if(chrs.contains(k+1)){
					this.region_names[k]=new String[num_regions_chr.get(k+1)];
					this.region_coordinates[k]=new int[num_regions_chr.get(k+1)][2];
				}else{
					this.region_names[k]=new String[0];
					this.region_coordinates[k]=new int[0][2];
				}
				
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
			//System.out.println("Finished loading data for rare variants analysis.");
		}catch(Exception e){e.printStackTrace();}	
	}
}
