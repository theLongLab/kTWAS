package nam;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;

public class Founder {
	public int[][] snps; //chr, location
	public ArrayList<String> sample_id;//ecotypeID 
	public char[][][]genotype;//founders genotype, id, chr, snp
//	public String[] pheno_names; //name of phenotypes
//	public String[][] phenotype; //phenotype values, id, phenotypes, missing phenotype NaN
	
//	public Founder(String genotype_250k, String f_pheno, Cross_info cross){	
	public Founder(String genotype_250k,  Cross_info cross){	
		try{
//			 read 250k data, put founders information in sample_id, and get to know marker number for each chromosome
			 ArrayList<Integer> founder_index = new ArrayList<Integer>();
			 this.sample_id =new ArrayList<String>();
			 int chr_snp =0; String chr ="1"; 
			 ArrayList<Integer> chr_snp_num = new ArrayList<Integer>();			 
			 BufferedReader br= new BufferedReader(new FileReader(genotype_250k));
			 String line = br.readLine();  
			 while(line!=null){
				 String[] array =line.split(",");
				 if(array[0].equals("Chromosome")){
					 for(int i=2; i<array.length; i++){
						 if(cross.ecoID_natureName.keySet().contains(array[i])){
							 founder_index.add(i);
							 this.sample_id.add(array[i]);
//							 System.out.println(i);
						 }
					 }
//					 System.out.println(cross.ecoID_natureName.size()+"\t"+this.sample_id.size());
					 if(cross.ecoID_natureName.size()!=this.sample_id.size()){
						 System.out.println("some sample are missing\t"+cross.ecoID_natureName.size()+"\t"+this.sample_id.size());
					 }
				 }else {
					 if(array[0].equals(chr)){
						 chr_snp ++;
					 }else{
						 chr_snp_num.add(chr_snp);
						 chr_snp=1; chr =array[0];
					 }
				 }
				 line = br.readLine();
			 }
			 chr_snp_num.add(chr_snp);
			 int chr_num =chr_snp_num.size();
//			 for(int i=0; i<chr_num; i++){
//				 System.out.println(chr_snp_num.get(i));
//			 }
//			 System.out.println(chr_snp_num.size());
//			 read the 250k file again, put the gentype of founders in array
			 this.genotype = new char[this.sample_id.size()][chr_num][];
			 this.snps= new int[chr_num][];

			 for(int i=0; i<chr_num; i++){
				 this.snps[i] = new int[chr_snp_num.get(i)]; 
				 for(int s =0; s<this.sample_id.size(); s++){
					 this.genotype[s][i] = new char[chr_snp_num.get(i)];
					 
				 }			 
			 }
			 br= new BufferedReader(new FileReader(genotype_250k));
			 line = br.readLine(); int marker_index =0; chr ="1"; int current_chr=0;
			 while(line!=null){
				 String[] array =line.split(",");				
				 if(!array[0].equals("Chromosome")){	
					 int sample_index =0;
					 if(array[0].equals(chr)){
						 this.snps[current_chr][marker_index]=Integer.parseInt(array[1]);
						 for(int s =2; s<array.length; s++){
							 if(founder_index.contains(s)){
								 this.genotype[sample_index][current_chr][marker_index] =array[s].toCharArray()[0];
								 sample_index++;
							 }
						 }
						 marker_index++;
					 }else{
						 current_chr++; marker_index=0; chr=array[0];
						 this.snps[current_chr][marker_index]=Integer.parseInt(array[1]);
						 for(int s =2; s<array.length; s++){
							 if(founder_index.contains(s)){
								 this.genotype[sample_index][current_chr][marker_index] =array[s].toCharArray()[0];
								 sample_index++;
							 }
						 }
						 marker_index++;
					 }
				 }
				 line = br.readLine();
			 }
////	read founder phenotype file put 
//			 BufferedReader phebr = new BufferedReader(new FileReader(f_pheno));
//			 String lfphe =  phebr.readLine();
//			 String[] phenos= lfphe.split("\t");
//			 this.pheno_names = new String[phenos.length-1];
//			 this.phenotype = new String[this.sample_id.size()][phenos.length-1];
//			 for(int i=1; i<phenos.length; i++){
//				 this.pheno_names[i-1] =phenos[i];
//			 }
//			 while(lfphe!=null){
//				 String[] array = lfphe.split("\t");
//				 if(!array[0].equals("fouders")){
//					 int sample_index =this.sample_id.indexOf(cross.natureName_ecoID.get(array[0]));
//					 for(int i=1; i<array.length; i++){
//						 if(array[i].equals(".")) this.phenotype[sample_index][i-1] ="NaN";
//						 else this.phenotype[sample_index][i-1] =array[i];						 
//					 }					
//				 }
//				 lfphe =  phebr.readLine();
//			 }
		}catch(Exception e){e.printStackTrace();}		
	}

}
