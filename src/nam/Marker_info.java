package nam;

import java.util.ArrayList;
import java.util.HashSet;

public class Marker_info {
	public int[] chr_count; //count of marker number for each chromosome
	public ArrayList<String> chr_name;
	public String[][] id;//chr, marker
	public double[][] cM;//chr, marker
	public int[][] locations;//chr, marker
	public int chr_num;
	
	public Marker_info(ArrayList<String> chrom, ArrayList<String> marker_name, ArrayList<String>centim, ArrayList<String> location){
		HashSet<String> chr_id =new HashSet<String>();
		for(int i=0; i<chrom.size(); i++){
			chr_id.add(chrom.get(i));
		}
		this.chr_count = new int[chr_id.size()];
		this.chr_num =chr_id.size();
		this.id =new String[chr_id.size()][];
		this.cM = new double[chr_id.size()][];
		this.locations =new int[chr_id.size()][];
		this.chr_name = new ArrayList<String>();
		String cchr ="1"; int chr_mark_num =0; int chr_index =0;
		this.chr_name.add("1");
		for(int i=0; i<chrom.size(); i++){
			if(chrom.get(i).equals(cchr)){
				chr_mark_num++;
			}else{		
				this.chr_count[chr_index]=chr_mark_num;
				chr_index++;
				chr_mark_num=1;
				cchr =chrom.get(i);
				chr_name.add(cchr);
			}
		}
		this.chr_count[chr_index]=chr_mark_num;
		for(int i=0; i<chr_id.size(); i++){
			this.id[i] = new String[this.chr_count[i]];
			this.cM[i]= new double[this.chr_count[i]];
			this.locations[i] = new int[this.chr_count[i]];
		}
		
		cchr ="1";chr_mark_num=0; chr_index =0;
		for(int i=0; i<chrom.size(); i++){
//			System.out.println(chrom.get(i)+"\t"+marker_name.get(i));
			if(chrom.get(i).equals(cchr)){
				this.id[chr_index][chr_mark_num] =marker_name.get(i);
				this.cM[chr_index][chr_mark_num] =Double.parseDouble(centim.get(i));
				this.locations[chr_index][chr_mark_num] =(int)(Double.parseDouble(location.get(i))*1000000);
				chr_mark_num++;
			}else{
//				System.out.println(chrom.get(i)+"\t"+marker_name.get(i));
				chr_index++;
				chr_mark_num=0;
				this.id[chr_index][chr_mark_num] =marker_name.get(i);
				this.cM[chr_index][chr_mark_num] =Double.parseDouble(centim.get(i));
				this.locations[chr_index][chr_mark_num] =(int)(Double.parseDouble(location.get(i))*1000000);
				cchr =chrom.get(i);
				chr_mark_num++;
			}
		}
	}
	
	public Marker_info(String[] chrom, String[] marker_name, String[]centim, String[] location){
		HashSet<String> chr_id =new HashSet<String>();
		for(int i=0; i<chrom.length; i++){
			chr_id.add(chrom[i]);
		}
		this.chr_count = new int[chr_id.size()];
		this.id =new String[chr_id.size()][];
		this.cM = new double[chr_id.size()][];
		this.locations =new int[chr_id.size()][];
		this.chr_name = new ArrayList<String>();
		String cchr ="1"; int chr_mark_num =0; int chr_index =0;
		this.chr_name.add("1");
		for(int i=0; i<chrom.length; i++){
			if(chrom[i].equals(cchr)){
				chr_mark_num++;
			}else{		
				this.chr_count[chr_index]=chr_mark_num;
				chr_index++;
				chr_mark_num=1;
				cchr =chrom[i];
				chr_name.add(cchr);
			}
		}
		this.chr_count[chr_index]=chr_mark_num;
		for(int i=0; i<chr_id.size(); i++){
			this.id[i] = new String[this.chr_count[i]];
			this.cM[i]= new double[this.chr_count[i]];
			this.locations[i] = new int[this.chr_count[i]];
		}
		
		cchr ="1";chr_mark_num=0; chr_index =0;
		for(int i=0; i<chrom.length; i++){
			if(chrom[i].equals(cchr)){
				this.id[chr_index][chr_mark_num] =marker_name[i];
				this.cM[chr_index][chr_mark_num] =Double.parseDouble(centim[i]);
				this.locations[chr_index][chr_mark_num] =(int)(Double.parseDouble(location[i])*1000000);
				chr_mark_num++;
			}else{
				chr_index++;
				chr_mark_num=0;
				this.id[chr_index][chr_mark_num] =marker_name[i];
				this.cM[chr_index][chr_mark_num] =Double.parseDouble(centim[i]);
				this.locations[chr_index][chr_mark_num] =(int)(Double.parseDouble(location[i])*1000000);
				cchr =chrom[i];
			}
		}
	}

}
