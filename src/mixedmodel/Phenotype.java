package mixedmodel;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import myMathLib.StatFuncs;
import myMathLib.Test;

import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;


public class Phenotype {
	public String phe_id;
	public String[] sample_ids;
	public double[] values;
	public double[] new_Y_4emmax;
	String transform_approach="NoTransform";
	double normality_pvaule=-1;
	public HashMap<String, Integer> sample_id2index;
	boolean NaN_removed=false;
	
	public Phenotype(String infile){
		ArrayList<String> ids= new ArrayList<>();
		ArrayList<Double> values= new ArrayList<>();
		try{
			BufferedReader br =new BufferedReader(new FileReader(infile));
			String line=br.readLine();//header
			this.phe_id= line.split("\t")[1];
			line=br.readLine();
			while(line!=null){
				String[] temp = line.split("\t");
				ids.add(temp[0]);
				if(temp[1].equals("NA")||temp[1].equals("NaN")||temp[1].equals("N")||temp[1].equals("*")){
					values.add(Double.NaN);
				}else{
					values.add(Double.parseDouble(temp[1]));
				}				
				line=br.readLine();
			}
			this.sample_ids=new String[ids.size()];
			this.values=new double[ids.size()];
			for(int k=0;k<ids.size();k++){
				this.sample_ids[k]=ids.get(k);
				this.values[k]=values.get(k);
			}
		}catch(Exception e){e.printStackTrace();}
		remove_NaN();
		setup_id2index();
	}
	
	public void setup_id2index(){
		this.sample_id2index= new HashMap<>();
		for(int i=0;i<this.sample_ids.length;i++)
			this.sample_id2index.put(this.sample_ids[i], i);
	}
	
	public Phenotype(String phe_id, String[] sample_ids, double[] values){
		this.phe_id=phe_id;
		this.sample_ids=sample_ids;
		this.values=values;
		remove_NaN();
		setup_id2index();
	}
	
	/*
	 * in this constructor, the NaN will not be removed
	 * 
	 * if the motivation is to do GWAS for a phenotype, there is NO case that NaN should be retained,
	 * however, if we use the phenotype to calculate kinship, or using multiple phenotypes in the 
	 * same analysis, sometimes NaN should be retained to keep the sample size the same across all 
	 * phenotypes 
	 */
	public Phenotype(String phe_id, String[] sample_ids, double[] values, boolean not_remove_NaN){
		this.phe_id=phe_id;
		this.sample_ids=sample_ids;
		this.values=values;
		setup_id2index();
	}
	
	public void transform(){
		this.normality_pvaule=StatFuncs.normality_test_ks(this.values);
		if(normality_pvaule>0.5){
			return;
		}
		System.out.println("Phenotype may not follow normal distribution. K-S pvalue="+normality_pvaule);
		String[] potential_transforms={"\"log\"","\"exp\"","\"sqrt\"","\"sqr\""};
		double[] pvalues=new double[potential_transforms.length];
		double[][] new_values=new double[potential_transforms.length][this.values.length];
		double min=StatUtils.min(values);
		if(min<=0){
			for(int k=0;k<this.values.length;k++){
				new_values[0][k]=Math.log(values[k]-min+1);
				new_values[1][k]=Math.exp(values[k]);
				new_values[2][k]=Math.sqrt(values[k]-min);
				new_values[3][k]=(values[k]-min)*(values[k]-min);
			}
		}else{
			for(int k=0;k<this.values.length;k++){
				new_values[0][k]=Math.log(values[k]);
				new_values[1][k]=Math.exp(values[k]);
				new_values[2][k]=Math.sqrt(values[k]);
				new_values[3][k]=(values[k])*(values[k]);
			}
		}		
		for(int i=0;i<potential_transforms.length;i++){
			pvalues[i]=StatFuncs.normality_test_ks(new_values[i]);
			System.out.println("Trying "+potential_transforms[i]+" transformation. P-value="+pvalues[i]+".");
		}
		int best_trans=-1;
		double the_pvalue=normality_pvaule;
		for(int i=0;i<potential_transforms.length;i++){
			if(pvalues[i]>the_pvalue){
				best_trans=i;
				the_pvalue=pvalues[i];
			}
		}
		if(best_trans!=-1){
			System.out.println("Used "+potential_transforms[best_trans]+" transformation. New K-S P-value="+pvalues[best_trans]);
			this.transform_approach=potential_transforms[best_trans];
			this.values=new_values[best_trans];
			normality_pvaule=pvalues[best_trans];
		}
	}
	
	public void substract_pedigree_mean_for_NAM(){
		HashMap<String, Integer> counts= new HashMap<>();
		HashMap<String, Double> sum= new HashMap<>();
		for(int i=0;i<this.sample_ids.length;i++){
			String pedigree_id=this.sample_ids[i].split("_")[0];
			myFileFunctions.FileFunc.add2_counts_hashmap(counts, pedigree_id, 1);
			myFileFunctions.FileFunc.add2_sum_hashmap(sum, pedigree_id, this.values[i]);
		}for(int i=0;i<this.sample_ids.length;i++){
			String pedigree_id=this.sample_ids[i].split("_")[0];
			this.values[i]=this.values[i]-sum.get(pedigree_id)/counts.get(pedigree_id);
		}
		System.out.println("Finished subtracting pedigree mean for NAM: "+this.phe_id);
	}
	
	void remove_NaN(){
		int nan=0;
		for(int i=0;i<sample_ids.length;i++){
			if(Double.isNaN(this.values[i])){
				nan++;
			}
		}
		String[] new_ids=new String[sample_ids.length-nan];
		double[] new_values=new double[sample_ids.length-nan];
		int index=0;
		for(int i=0;i<sample_ids.length;i++){
			if(!Double.isNaN(this.values[i])){
				new_ids[index]=sample_ids[i];
				new_values[index]=values[i];
				index++;
			}
		}this.sample_ids=new_ids;
		this.values=new_values;
		this.NaN_removed=true;
		setup_id2index();
	}
	
	/*
	 * 	generate new phenotype array by applying decomposition array, 
	 * only samples indexed in sampel_pheno_index[] applied
	 */
	public void generate_new_Y_by_multiplying(double[][] decomposed_array){	
		this.new_Y_4emmax=new double[sample_ids.length];
		if(decomposed_array.length!=values.length){
			System.out.println("Decomposition array is inconsistent to phenotype array!");
			return;
		}
		for(int i=0;i<values.length;i++){
			for(int j=0;j<values.length;j++){
				this.new_Y_4emmax[i]=this.new_Y_4emmax[i]+decomposed_array[i][j]*values[j];
			}
		}
	}
	
	public double[] extract_values_with_id_list(String[] sample_ids_list){
		this.setup_id2index();
		double[] the_data=new double[sample_ids_list.length];
		for(int i=0;i<sample_ids_list.length;i++){
			if(this.sample_id2index.containsKey(sample_ids_list[i])){
				the_data[i]=this.values[this.sample_id2index.get(sample_ids_list[i])];
			}else{
				System.out.println("!this.sample_id2index.containsKey(sample_ids_list[i])");
				the_data[i]=Double.NaN;
			}
		}
		return the_data;
	}
	
	public Phenotype extract_sample_with_id_list(String[] sample_ids_list){
		this.setup_id2index();
		double[] the_data=new double[sample_ids_list.length];
		for(int i=0;i<sample_ids_list.length;i++){
			if(this.sample_id2index.containsKey(sample_ids_list[i])){
				the_data[i]=this.values[this.sample_id2index.get(sample_ids_list[i])];
			}else{
				System.out.println("!this.sample_id2index.containsKey(sample_ids_list[i]): "+sample_ids_list[i]);
				the_data[i]=Double.NaN;
			}
		}
		return new Phenotype(this.phe_id, sample_ids_list, the_data);
	}
	
	public double[] extract_samples_with_id_list(ArrayList<String> sample_ids_list){
		this.setup_id2index();
		double[] the_data=new double[sample_ids_list.size()];
		for(int i=0;i<sample_ids_list.size();i++){
			if(this.sample_id2index.containsKey(sample_ids_list.get(i))){
				the_data[i]=this.values[this.sample_id2index.get(sample_ids_list.get(i))];
			}else{
				System.out.println("!this.sample_id2index.containsKey(sample_ids_list.get(i))");
				the_data[i]=Double.NaN;
			}
		}
		return the_data;
	}
	
	public double sample_mean(){
		if(!this.NaN_removed){
			this.remove_NaN();
			this.setup_id2index();
		}return StatFuncs.mean_no_NaN(this.values);
	}
	
	public void permute(){
		Test.randomPermute(this.values);
	}
	
	public void write2file(String out_file){
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(out_file));
			bw.write("ID\t"+this.phe_id+"\n");
			for(int k=0;k<this.sample_ids.length;k++){
				bw.write(this.sample_ids[k]+"\t"+this.values[k]+"\n");
			}bw.close();
		}catch(Exception e){e.printStackTrace();}
	}
	
	/*
	 * chr start from 0.
	 */
	public void preprocess_cofactor(VariantsDouble genotype, int chr, int loc, String output_file){
		double[] full_snp=genotype.load_one_variant_by_location(chr, loc);
		double[] matched_snp=new double[this.sample_ids.length];
		double total=0, count=0;
		HashSet<String> matched= new HashSet<>();
		for(int k=0;k<genotype.sample_size;k++){
			if(this.sample_id2index.containsKey(genotype.sample_ids[k])){
				matched_snp[this.sample_id2index.get(genotype.sample_ids[k])]=full_snp[k];
				matched.add(genotype.sample_ids[k]);
				count++;
				total+=full_snp[k];
			}
		}
		double average=total/count;
		for(int i=0;i<this.sample_ids.length;i++){
			if(!matched.contains(this.sample_ids[i])){
				matched_snp[i]=average;
			}
		}
		preprocess_cofactor(matched_snp, output_file);
		
	}
	
	public void preprocess_cofactor(double[] cofactor, String output_file){		
		OLSMultipleLinearRegression reg1=new OLSMultipleLinearRegression();
		double[] Y=this.values.clone();
		double[][] Xs=new double[Y.length][1];
		for(int m=0;m<Y.length;m++)Xs[m][0]=cofactor[m];
		reg1.newSampleData(Y, Xs);
		Phenotype regressed=new Phenotype(phe_id, sample_ids, reg1.estimateResiduals());
		regressed.write2file(output_file);
	}
	
	public void jackknife(int num, String folder){
		try{
			if(num>=this.sample_ids.length-100){
				System.out.println("There will be only "+(this.sample_ids.length-num)+" samples left!");
			}
			for(int round=0;round<this.sample_ids.length;round++){
				String out_file=folder+"jackknife"+round+".tsv";
				BufferedWriter bw=new BufferedWriter(new FileWriter(out_file));
				HashSet<String> removed= new HashSet<>();
				if(num==1){
					removed.add(this.sample_ids[round]);
				}else{
					for(int i=0;i<num;i++){
						int the_index=(int)(Test.randomNumber()*this.sample_ids.length);
						String the_sample_id=this.sample_ids[the_index];
						while(removed.contains(the_sample_id)){
							the_index=(int)(Test.randomNumber()*this.sample_ids.length);
							the_sample_id=this.sample_ids[the_index];
						}
						removed.add(the_sample_id);
					}
				}
				bw.write("ID\t"+this.phe_id+".jknf"+round+"\n");
				for(int k=0;k<this.sample_ids.length;k++){
					if(!removed.contains(this.sample_ids[k]))
						bw.write(this.sample_ids[k]+"\t"+this.values[k]+"\n");
				}bw.close();
			}			
		}catch(Exception e){e.printStackTrace();}
		
	}
}
