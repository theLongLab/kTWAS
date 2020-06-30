package nam;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

public class RILs {
//	public  int pedNum;
	public  char[][][] genotype;   // #sample x #chr x #marker, A, B, and H missing 'N', attention! original file use D stands for missing
	public  int pedSampleSize;         // sample_size_of_pedigree do not including two founders
	public  Marker_info marker_info;     //3 x #markers // 3 fields: locus_id[chr],loc_cM[chr], loc_location[chr]
	public	String[] parents;            //2 parents' ecotypeID and the third column is crossid 
	public	String[] rilID;            //#samplesOfpedigree (no founder)
	public	String[][] phenotypes;    //(#samplesOfpedigree(not including founder) x #phenotype missing "NaN"
	public	String[] pheno_name;       //Store the name of phenotype	

//	 public RILs(File pedFiles, String phe_ril, Cross_info cross, String impute){
	public RILs(File pedFiles, String phe_ril, Cross_info cross){
			try{ 
//				read pedigree file, the first 5 rows are map info, the 6 and 7 rows are founders, and then offsprings
				if(!pedFiles.exists()) System.out.println("file "+pedFiles +" does not exist, please check ");
				else{
					ArrayList<String> chrom = new ArrayList<String>();
					ArrayList<String> marker_name = new ArrayList<String>();
					ArrayList<String> centim = new ArrayList<String>();
					ArrayList<String> location = new ArrayList<String>();
					
					BufferedReader br = new BufferedReader(new FileReader(pedFiles));
					String line = br.readLine();

					int sample =0;//samplesize include founders
					while(line!=null){				
						String[] array = line.split("\t");					
						if(array[0].equals("Locus")){
							for(int i=1; i<array.length; i++){
								marker_name.add(array[i]);
							}
						}else if(array[0].equals("Chr")){
							for(int i=1; i<array.length; i++){
								chrom.add(array[i]);
							}
						}else if(array[0].equals("CR_cM")){
							for(int i=1; i<array.length; i++){
								centim.add(array[i]);
							}
						}else if(array[0].equals("Start_Mb")){
							for(int i=1; i<array.length; i++){
								location.add(array[i]);
							}
						}else if(!array[0].equals("order")){
							sample++;
						}				
						line = br.readLine();						
					}
					this.marker_info = new Marker_info(chrom, marker_name, centim, location);
					this.pedSampleSize = sample-2;//column2 for samplesize
					this.genotype= new char[sample-2][this.marker_info.chr_num][];
					this.rilID = new String[sample-2];
					this.parents = new String[3];
					for(int s=0; s<this.genotype.length; s++){
						for(int c=0; c<this.genotype[s].length; c++){
							this.genotype[s][c] =new char[this.marker_info.chr_count[c]];
						}
					}
//read the pedigree file again, put parents info and RIL genotype
					HashSet<String> missing =new HashSet<String>();
					missing.add("N"); missing.add("D"); missing.add("_"); missing.add(".");
					br = new BufferedReader(new FileReader(pedFiles));
					line = br.readLine(); 				
					int row_index=0;
					while(line!=null){	
						row_index++;
						String[] array = line.split("\t");	
						if(row_index==6){
							if(!cross.natureName_ecoID.keySet().contains(array[0])){
								System.out.println("Sorry no cross informaiton");
							}else{
								this.parents[0]=cross.natureName_ecoID.get(array[0]);
							}
						}
						if(row_index==7){
							if(!cross.natureName_ecoID.keySet().contains(array[0])){
								System.out.println("Sorry no cross informaiton");
							}else{
								this.parents[1]=cross.natureName_ecoID.get(array[0]);
							}
						}
						this.parents[2]=cross.ecoID_crossid.get(this.parents[0]+"_"+this.parents[1]);//put cross id
						if(row_index>7){
							int rid =Integer.parseInt(array[0]);	//remove potential 0 at begining, 002 to be 2				
							this.rilID[row_index-8]=this.parents[2]+"_"+rid;//add cross id to rid, 2_2 means cross 2 sample 2
							String current_chr ="1";int chr_snp_num=0; int chr_index =0;
							for(int i=0; i<chrom.size(); i++){
								String putsnp ="";
								if(!array[i+1].equals("A") && !array[i+1].equals("B")){
									if(array[i+1].equals("C") || array[i+1].equals("H") ){
										putsnp="H";
									}else if(missing.contains(array[i+1])){
										putsnp="N";
									}
								}else{
									putsnp =array[i+1];
								}
								if(chrom.get(i).equals(current_chr)){
									this.genotype[row_index-8][chr_index][chr_snp_num]=putsnp.toCharArray()[0];
									chr_snp_num++;
								}else{
									chr_index++;
									chr_snp_num=0;
									this.genotype[row_index-8][chr_index][chr_snp_num]=putsnp.toCharArray()[0];
									current_chr =chrom.get(i);
									chr_snp_num++;
								}
							}
						}									
						line = br.readLine();
					}
//					impute RIL				
//					if(impute.equals("Y")){
					for(int id=0; id<this.genotype.length; id++){
						for(int c=0; c<this.genotype[id].length; c++){
							char left ='*'; char right='*'; int noNleft =0; int noNright=0; int stat =0;
							for(int m=0; m<this.genotype[id][c].length; m++){	
								//									System.out.println(this.genotype[id][c][m]);
								if(this.genotype[id][c][m]!='N'){
									if(stat==2) {
										right =this.genotype[id][c][m]; noNright=m;
										impute(this.genotype[id][c], left, right, noNleft, noNright,this.marker_info.locations[c]);
									}
									stat=1;
								}else{
									if(stat==1){
										left =this.genotype[id][c][m-1]; noNleft =m-1;
									}
									if(m==this.genotype[id][c].length-1){
										//											System.out.println(id+"\t"+c+"\t"+m+"\t"+left);
										noNright=m; right='*';
										impute(this.genotype[id][c], left,right, noNleft, noNright, this.marker_info.locations[c]);
									}
									stat=2;
								}
							}
						}
//						}
					}
				}

				
//	read phenotype files 
				BufferedReader brphe = new BufferedReader(new FileReader(phe_ril));
				String pheline = brphe.readLine();
				String[] phename = pheline.split("\t");
				this.phenotypes= new String[this.pedSampleSize][phename.length-2];
				this.pheno_name = new String[phename.length-2];
				for(int i=2; i<phename.length; i++){
					this.pheno_name[i-2]=phename[i];
				}
				int id_index=0;
				while(pheline!=null){
					String[]array =pheline.split("\t");
					if(!array[0].equals("family")){
						if(array[0].split("RV")[0].equals(this.parents[2])){
							String id =array[1].split("RV")[0]+"_"+array[1].split("RV")[1];
							if(id.equals(this.rilID[id_index])){
								for(int i=2; i<array.length; i++){
									if(array[i].equals(".")||array[i].equals("NA")||array[i].equals("?")||array[i].equals("/")
											||array[i].equals("*")||array[i].equals("N")||array[i].equals("-")){
										this.phenotypes[id_index][i-2]="NaN";
									}else{
										this.phenotypes[id_index][i-2]=array[i];
									}
								}
							}else{
								System.out.println("Sample id doesn't match "+id+"\t"+this.rilID[id_index]);
							}
							id_index++;
						}
					}
					pheline = brphe.readLine();					
				}
			}catch(Exception e){e.printStackTrace();}
	   }
	 
	 public  void impute(char[] ingenotype, char left, char right, int start, int end, int[] loc){
		 if(right!=left){			
				if(left=='*'){
					for(int i=start; i<end; i++){
						ingenotype[i] =right; 
					}
				}else if(right=='*'){
					for(int i=start+1; i<end+1; i++){
						ingenotype[i] =left; 
					}
				}else if(right=='H' || left=='H'){
					double r=myMathLib.Test.randomNumber();
					for(int i=start+1; i<end; i++){
						double rig_loc=loc[start];
						double lef_loc=loc[end];
						if(((loc[i]-lef_loc)/(rig_loc-lef_loc))>=r){
							ingenotype[i] =right;
						}else{
							ingenotype[i] =left;
						}
					}
				}else{
					for(int i=start+1; i<end; i++){
						ingenotype[i]='H';
					}
				}
			}else{
				for(int i=start+1; i<end; i++){
					ingenotype[i]=left;
				}
			}		 
	 }
	 
}
