package mixedmodel;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import org.apache.commons.math3.linear.Array2DRowFieldMatrix;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;

import ch.systemsx.cisd.base.mdarray.MDDoubleArray;
import ch.systemsx.cisd.hdf5.HDF5MDDataBlock;

public class RelationMatrix {
	
	public RealMatrix relation_matrix;
	//public double[][] kinship_matrix_array;
	public String[] obs_ids;
	public String[] pred_ids;
	public HashMap<String, Integer> ids2indexes;	
	boolean global;
	int[] local;
	
	public RelationMatrix(RealMatrix relation_matrix, String[] obs_ids, String[] pred_ids){
		this.relation_matrix=relation_matrix;
		this.obs_ids=obs_ids;
		this.pred_ids=pred_ids;
		this.assign_ids2indexes();
	}
	public RelationMatrix(String relation_matrix_file){
		try{
			BufferedReader br=new BufferedReader(new FileReader(relation_matrix_file));
			String line=br.readLine();
			this.global=Boolean.parseBoolean(line.split("=")[1]);
			line=br.readLine();
			// TODO assign this.local
			line=br.readLine();
			String[] headers_pred=line.split(",");
			this.pred_ids=new String[headers_pred.length-1];
			for(int k=0;k<pred_ids.length;k++)this.pred_ids[k]=headers_pred[k+1];
			line=br.readLine();
			ArrayList<String> obs_ids_array=new ArrayList<String>();
			ArrayList<double[]> data_array=new ArrayList<double[]>();
			while(line!=null){
				String[] tmp=line.split(",");
				obs_ids_array.add(tmp[0]);
				double[] the_line=new double[this.pred_ids.length];
				for(int k=0;k<pred_ids.length;k++)the_line[k]=Double.parseDouble(tmp[k+1]);
				data_array.add(the_line);
				line=br.readLine();
			}
			this.obs_ids=new String[obs_ids_array.size()];
			double[][] relationship=new double[pred_ids.length][obs_ids.length];
			for(int i=0;i<obs_ids.length;i++){
				this.obs_ids[i]=obs_ids_array.get(i);
				for(int j=0;j<pred_ids.length;j++){
					relationship[j][i]=data_array.get(i)[j];
				}
			}
			this.relation_matrix=new Array2DRowRealMatrix(relationship, false);
			data_array.clear();
			this.assign_ids2indexes();
		}catch(Exception e){e.printStackTrace();}	
	}
	
	public RelationMatrix(VariantsDouble data_obs, VariantsDouble data_pred){
		this.obs_ids=data_obs.sample_ids.clone();
		this.pred_ids=data_pred.sample_ids.clone();
		this.global=true;
		this.local=null;
		this.assign_ids2indexes();
		this.calculate_WG_matrix(data_obs, data_pred);
	}
	
	public RelationMatrix(VariantsDouble data_obs, VariantsDouble data_pred,
			int[] local){
		this.obs_ids=data_obs.sample_ids.clone();
		this.pred_ids=data_pred.sample_ids.clone();
		this.global=false;
		this.local=local.clone();	
		this.assign_ids2indexes();
		this.calculate_local_matrix(data_obs, data_pred, local);
	}
	
	public void calculate_local_matrix(VariantsDouble data_obs, VariantsDouble data_pred,
			int[] region){
		double[][] relationship=new double[data_pred.sample_size][data_obs.sample_size];
		double[][] total_num_var_used=new double[this.pred_ids.length][this.obs_ids.length];	
		int chr=region[0], start=region[1], end=region[2];		
		if(data_obs.num_sites[chr]!=data_pred.num_sites[chr]){
			System.out.println("Number of variants in chromosome not equal: chr"+(chr+1));
			System.out.println(data_obs.num_sites[chr]+":"+data_pred.num_sites[chr]);
			//return;
		}
		double[][] all_snp_obs=data_obs.load_variants_in_region(chr, start, end);
		double[][] all_snp_pred=data_pred.load_variants_in_region(chr, start, end);
		if(all_snp_obs.length!=all_snp_pred.length){
			System.out.println("Number of variants in region not equal");
			System.out.println(all_snp_obs.length+":"+all_snp_pred.length);
		}
		for(int var_index=0;var_index<all_snp_obs.length;var_index++){
				//if(var_index%1000==0)System.out.println("chr"+(chr+1)+":"+var_index/1000+"K.");
			double[] snp_obs=all_snp_obs[var_index];
			double[] snp_pred=all_snp_pred[var_index];
			double mean_obs=myMathLib.StatFuncs.mean_NaN(snp_obs);
			double mean_pred=myMathLib.StatFuncs.mean_NaN(snp_pred);
			if(Double.isNaN(mean_obs)||Double.isNaN(mean_pred)){
				System.out.println("All subjects are NaN at a SNP in Chr"+(chr+1));
				continue;
			}
			double sd_obs=Math.sqrt(myMathLib.StatFuncs.var_NaN(snp_obs,mean_obs));
			double sd_pred=Math.sqrt(myMathLib.StatFuncs.var_NaN(snp_pred,mean_pred));
			double[] new_snp_obs=new double[data_obs.sample_size];
			double[] new_snp_pred=new double[data_pred.sample_size];
			for(int k=0;k<data_obs.sample_size;k++)new_snp_obs[k]=(snp_obs[k]-mean_obs)/sd_obs;
			for(int k=0;k<data_pred.sample_size;k++)new_snp_pred[k]=(snp_pred[k]-mean_pred)/sd_pred;
			for(int i=0;i<data_obs.sample_size;i++){
				for(int j=0;j<data_pred.sample_size;j++){
					if((!Double.isNaN(new_snp_obs[i])) && (!Double.isNaN(new_snp_pred[j]))){
						relationship[j][i]=relationship[j][i]+(new_snp_obs[i]*new_snp_pred[j]);
						total_num_var_used[j][i]++;
					}
				}
			}
		}	
		System.out.println("Finished relation for region: Chr"+(chr+1)+":"+start+":"+end);			
		for(int i=0;i<data_obs.sample_size;i++){
			for(int j=0;j<data_pred.sample_size;j++){
				relationship[j][i]=relationship[j][i]/total_num_var_used[j][i];
			}
		}
		this.relation_matrix= new Array2DRowRealMatrix(relationship, false); 
	}
	
	public void calculate_WG_matrix(VariantsDouble data_obs, VariantsDouble data_pred){
		double[][] relationship=new double[data_pred.sample_size][data_obs.sample_size];
		double[][] total_num_var_used=new double[this.pred_ids.length][this.obs_ids.length];	
		for(int chr=0;chr<data_obs.num_chrs;chr++){				
			if(data_obs.num_sites[chr]!=data_pred.num_sites[chr]){
				System.out.println("Number of variants not equal: chr"+(chr+1));
				int var_index2=0;
				int count_shared=0;
				for(int var_index=0;var_index<data_obs.num_sites[chr];var_index++){
					if(var_index%1000==0)System.out.println("chr"+(chr+1)+":"+var_index/1000+"K.");		
					while(var_index2<data_pred.locations[chr].length &&
							data_pred.locations[chr][var_index2]<data_obs.locations[chr][var_index]){
						var_index2++;
					}
					if(var_index2>=data_pred.locations[chr].length)break;
					if(data_pred.locations[chr][var_index2]==data_obs.locations[chr][var_index]){ // found
						count_shared++;
						double[] snp_obs=data_obs.load_one_variant_by_index(chr, var_index);
						double[] snp_pred=data_pred.load_one_variant_by_index(chr, var_index2);						
						double mean_obs=myMathLib.StatFuncs.mean(snp_obs);
						double mean_pred=myMathLib.StatFuncs.mean(snp_pred);
						if(Double.isNaN(mean_obs)||Double.isNaN(mean_pred)){
							System.out.println("All subjects are NaN at a SNP in Chr"+(chr+1));
							continue;
						}
						double sd_obs=Math.sqrt(myMathLib.StatFuncs.var(snp_obs,mean_obs));
						double sd_pred=Math.sqrt(myMathLib.StatFuncs.var(snp_pred,mean_pred));
						if(Double.isNaN(sd_obs)||Double.isNaN(sd_pred)){
							//System.out.println("SD is NaN Chr"+(chr+1));
							continue;
						}
						double[] new_snp_obs=new double[data_obs.sample_size];
						double[] new_snp_pred=new double[data_pred.sample_size];
						for(int k=0;k<data_obs.sample_size;k++)new_snp_obs[k]=(snp_obs[k]-mean_obs)/sd_obs;
						for(int k=0;k<data_pred.sample_size;k++)new_snp_pred[k]=(snp_pred[k]-mean_pred)/sd_pred;
						for(int i=0;i<data_obs.sample_size;i++){
							if(!Double.isNaN(new_snp_obs[i])){
								for(int j=0;j<data_pred.sample_size;j++){
									if(!Double.isNaN(new_snp_pred[j])){
										relationship[j][i]=relationship[j][i]+(new_snp_obs[i]*new_snp_pred[j]);
										total_num_var_used[j][i]++;
									}
								}
							}
							
						}
					}					
				}	
				System.out.println("Finished Chr"+(chr+1)+": #SNPs="+count_shared);
			}else{
				for(int var_index=0;var_index<data_obs.num_sites[chr];var_index++){
					//if(var_index%1000==0)System.out.println("chr"+(chr+1)+":"+var_index/1000+"K.");
					double[] snp_obs=data_obs.load_one_variant_by_index(chr, var_index);
					double[] snp_pred=data_pred.load_one_variant_by_index(chr, var_index);
					double mean_obs=myMathLib.StatFuncs.mean(snp_obs);
					double mean_pred=myMathLib.StatFuncs.mean(snp_pred);
					if(Double.isNaN(mean_obs)||Double.isNaN(mean_pred)){
						System.out.println("All subjects are NaN at a SNP in Chr"+(chr+1));
						continue;
					}
					double sd_obs=Math.sqrt(myMathLib.StatFuncs.var(snp_obs,mean_obs));
					double sd_pred=Math.sqrt(myMathLib.StatFuncs.var(snp_pred,mean_pred));
					if(Double.isNaN(sd_obs)||Double.isNaN(sd_pred)){
						//System.out.println("All subjects are NaN at a SNP in Chr"+(chr+1));
						continue;
					}
					double[] new_snp_obs=new double[data_obs.sample_size];
					double[] new_snp_pred=new double[data_pred.sample_size];
					for(int k=0;k<data_obs.sample_size;k++)new_snp_obs[k]=(snp_obs[k]-mean_obs)/sd_obs;
					for(int k=0;k<data_pred.sample_size;k++)new_snp_pred[k]=(snp_pred[k]-mean_pred)/sd_pred;
					for(int i=0;i<data_obs.sample_size;i++){
						if(!Double.isNaN(new_snp_obs[i])){
							for(int j=0;j<data_pred.sample_size;j++){
								if(!Double.isNaN(new_snp_pred[j])){
									relationship[j][i]=relationship[j][i]+(new_snp_obs[i]*new_snp_pred[j]);
									total_num_var_used[j][i]++;
								}
							}
						}
						
					}
				}	
				System.out.println("Finished Chr"+(chr+1));
			}
			
		}	
		for(int i=0;i<data_obs.sample_size;i++){
			for(int j=0;j<data_pred.sample_size;j++){
				relationship[j][i]=relationship[j][i]/total_num_var_used[j][i];
			}
		}
		this.relation_matrix= new Array2DRowRealMatrix(relationship, false); 
	}
	
	public RelationMatrix(String data_obs_hdf5, String data_pred_hdf5, int[] local){
		this((new VariantsDouble(data_obs_hdf5)), (new VariantsDouble(data_pred_hdf5)), local);
	}
	
	public RelationMatrix(String data_obs_hdf5, String data_pred_hdf5){
		this((new VariantsDouble(data_obs_hdf5)), (new VariantsDouble(data_pred_hdf5)));
	}
	
	public void write2file(String output_file){
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(output_file));
			bw.write("#WG="+this.global+"\n");
			bw.write("#");
			if(this.local!=null){				
				bw.write("Chr"+(this.local[0]+1)+":"+this.local[1]+":"+this.local[2]+";");
			}bw.write("\n");
			double[][] data=this.relation_matrix.getData();
			bw.write("OBS_ID");
			for(int i=0;i<pred_ids.length;i++){
				bw.write(","+pred_ids[i]);
			}bw.write("\n");
			for(int k=0;k<obs_ids.length;k++){
				bw.write(obs_ids[k]);
				for(int i=0;i<pred_ids.length;i++){
					bw.write(","+data[i][k]);
				}bw.write("\n");
			}
			bw.close();
		}catch(Exception e){e.printStackTrace();}		
	}
	
	public void assign_ids2indexes(){
		this.ids2indexes = new HashMap<String, Integer>();
		for(int k=0;k<obs_ids.length;k++)this.ids2indexes.put(obs_ids[k],k);
		for(int k=0;k<pred_ids.length;k++)this.ids2indexes.put(pred_ids[k],k);
	}
	
	/*
	 * subset_ids should be a subset of all the ids.
	 */
	public static RelationMatrix submatrix(KinshipMatrix full_matrix, String[] subset_ids_obs, String[] subset_ids_pred){
		double[][] data=new double[subset_ids_pred.length][subset_ids_obs.length];
		double[][] full_data=full_matrix.kinship_matrix.getData();
		for(int i=0;i<subset_ids_obs.length;i++){
			if(!full_matrix.ids2indexes.containsKey(subset_ids_obs[i])){
				System.out.println("RelationMatrix: Can't substract kinship since "+subset_ids_obs[i]+" isn't there!");
				return null;
			}
			int index_i=full_matrix.ids2indexes.get(subset_ids_obs[i]);
			for(int j=0;j<subset_ids_pred.length;j++){
				if(!full_matrix.ids2indexes.containsKey(subset_ids_pred[j])){
					System.out.println("RelationMatrix: Can't substract kinship since "+subset_ids_pred[j]+" isn't there!");
					return null;
				}
				int index_j=full_matrix.ids2indexes.get(subset_ids_pred[j]);
				data[j][i]=full_data[index_j][index_i];
			}
		}
		return new RelationMatrix(new Array2DRowRealMatrix(data, false), subset_ids_obs, subset_ids_pred);
	}
	
	public static void main(String[] args){
//		String folder="/Users/quanlong/Documents/projects2/predictions/RA_data/genotypes/";
//		String[] types={"etanercept","infliximab","adalimumab","all"};
//		for(int k=0;k<types.length;k++){
//			System.out.println(types[k]);
//			String data_obs_hdf5=folder+types[k]+".train.hdf5";
//			String data_pred_hdf5=folder+types[k]+".test.hdf5";
//			String output_file=folder+types[k]+".train.test.K.RRM";
//			RelationMatrix matrix=new RelationMatrix(data_obs_hdf5, data_pred_hdf5);
//			matrix.write2file(output_file);
//		}
		String data_obs_hdf5="/Users/quanlong/Documents/projects2/predictions/RA_data/genotypes/all.train.hdf5";
		String data_pred_hdf5="/Users/quanlong/Documents/projects2/predictions/RA_data/test_final/geno.final.10p.hdf5";
		String output_file="/Users/quanlong/Documents/projects2/predictions/RA_data/test_final/all.final.K.4.RRM";
		RelationMatrix matrix=new RelationMatrix(data_obs_hdf5, data_pred_hdf5);
		matrix.write2file(output_file);
	}
	
	
	
}
