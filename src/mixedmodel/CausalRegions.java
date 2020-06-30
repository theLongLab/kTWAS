package mixedmodel;

import java.util.Arrays;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;

public class CausalRegions {
	public int[][] regions;// [#regions][3]: chr,start, end
	double[] var_comp;
	public KinshipMatrix[] reginal_kinships;
	public RelationMatrix[] relational_matrixs;
	
	public KinshipMatrix global_kinship;
	public RelationMatrix global_relation;
	
	public KinshipMatrix combined_kinship;
	public RelationMatrix combined_relation;
	
	
	
	public CausalRegions(int[][] regions, double[] var_comp, String[] reginal_kinship_files, 
			String[] relation_matrix_files, KinshipMatrix global_kinship, RelationMatrix global_relation,
			VariantsDouble data_obs, VariantsDouble data_pred){
		this.regions=regions.clone();
		this.var_comp=var_comp.clone();
		for(int i=0;i<regions.length;i++)regions[i].clone();
		
		this.reginal_kinships=new KinshipMatrix[this.regions.length];
		this.relational_matrixs=new RelationMatrix[this.regions.length];
		this.global_kinship=global_kinship;
		this.global_relation=global_relation;
		// calculate kinships 
		for(int r=0;r<regions.length;r++){
			this.reginal_kinships[r]=new KinshipMatrix(data_obs, regions[r]);
			this.reginal_kinships[r].write2file(reginal_kinship_files[r]);
			this.relational_matrixs[r]=new RelationMatrix(data_obs, data_pred, regions[r]);
			this.relational_matrixs[r].write2file(relation_matrix_files[r]);
		}				
	}
	
	/*
	 * for simualtion only, given a full data, generate the observed and predicted local kinships
	 * 
	 * calculate and write local-kinships or load wrt whether geno_full==null
	 */
	public CausalRegions(int[][] regions, double[] var_comp, String[] reginal_full_kinship_files, 
			 KinshipMatrix global_kinship, RelationMatrix global_relation,
			 String[] obs_ids, String[] pred_ids){
		this.regions=regions.clone();
		this.var_comp=var_comp.clone();
		for(int i=0;i<regions.length;i++)regions[i].clone();
		
		this.reginal_kinships=new KinshipMatrix[this.regions.length];
		this.relational_matrixs=new RelationMatrix[this.regions.length];
		this.global_kinship=global_kinship;
		this.global_relation=global_relation;
		// calculate kinships 
		for(int r=0;r<regions.length;r++){
			KinshipMatrix the_local_full=new KinshipMatrix(reginal_full_kinship_files[r]);			
			this.reginal_kinships[r]=KinshipMatrix.submatrix(the_local_full, obs_ids);			
			this.relational_matrixs[r]=RelationMatrix.submatrix(the_local_full, obs_ids, pred_ids);
		}				
	}
	
	
	public CausalRegions(int[][] regions, double[] var_comp, String[] reginal_kinship_files, 
			String[] relation_matrix_files, KinshipMatrix global_kinship, RelationMatrix global_relation){
		this.regions=regions.clone();
		this.var_comp=var_comp.clone();
		for(int i=0;i<regions.length;i++)regions[i].clone();
		this.reginal_kinships=new KinshipMatrix[this.regions.length];
		this.relational_matrixs=new RelationMatrix[this.regions.length];
		this.global_kinship=global_kinship;
		this.global_relation=global_relation;
		// load kinships		
		for(int r=0;r<regions.length;r++){
			this.reginal_kinships[r]=new KinshipMatrix(reginal_kinship_files[r]);
			this.relational_matrixs[r]=new RelationMatrix(relation_matrix_files[r]);
		}				
	}
	
	public void generte_combined(double[] var_comp){
		if(var_comp!=null){
			if(var_comp.length!=regions.length){
				System.out.println("var_comp.length!=regions.length");
			}
		}
		this.combined_kinship=new KinshipMatrix(new Array2DRowRealMatrix(
				this.global_kinship.kinship_matrix.getData(),true), this.global_kinship.ids);
		this.combined_relation=new RelationMatrix(new Array2DRowRealMatrix(
				this.global_relation.relation_matrix.getData(),true), this.global_relation.obs_ids, this.global_relation.pred_ids);
		for(int r=0;r<this.reginal_kinships.length;r++){
			this.combined_kinship.kinship_matrix=
					this.combined_kinship.kinship_matrix.add
					(this.reginal_kinships[r].kinship_matrix.scalarMultiply(var_comp[r]));
			this.combined_relation.relation_matrix=
					this.combined_relation.relation_matrix.add
					(this.relational_matrixs[r].relation_matrix.scalarMultiply(var_comp[r]));
		}
	}
	
	public void generte_combined_no_global(double[] var_comp){
		if(var_comp!=null){
			if(var_comp.length!=regions.length){
				System.out.println("var_comp.length!=regions.length");
			}
		}
		this.combined_kinship=new KinshipMatrix(new Array2DRowRealMatrix(
				new double[this.global_kinship.ids.length][this.global_kinship.ids.length]), 
				this.global_kinship.ids);
		this.combined_relation=new RelationMatrix(new Array2DRowRealMatrix(
				new double[this.global_relation.pred_ids.length][this.global_relation.obs_ids.length]), 
				this.global_relation.obs_ids, this.global_relation.pred_ids);
		for(int r=0;r<this.reginal_kinships.length;r++){
			this.combined_kinship.kinship_matrix=
					this.combined_kinship.kinship_matrix.add
					(this.reginal_kinships[r].kinship_matrix.scalarMultiply(var_comp[r]));
			this.combined_relation.relation_matrix=
					this.combined_relation.relation_matrix.add
					(this.relational_matrixs[r].relation_matrix.scalarMultiply(var_comp[r]));
		}
	}
	
	public static int[][] generate_regions(String local_kinship_result_file, double[] var_comp){
		int num_chr=22,  win=100000;
		AssociationResults local=new AssociationResults(local_kinship_result_file, 0.9);
		int[][] regions=new int[num_chr][3];
		double[] ps=new double[num_chr];
		//double[] r2s=new double[num_chr];
		Arrays.fill(ps, 1); Arrays.fill(var_comp, 1); 
		for(int chr=0;chr<num_chr;chr++){
			regions[chr][0]=chr;
		}
		for(int index=0;index<local.num_of_var;index++){
			if(local.pvalue[index]<ps[local.chr[index]-1]){
				ps[local.chr[index]-1]=local.pvalue[index];
				var_comp[local.chr[index]-1]=local.AdjustedR2[index];//*(-Math.log(local.pvalue[index]));
				regions[local.chr[index]-1][1]=local.location[index];
				regions[local.chr[index]-1][2]=local.location[index]+win;
			}
		}
		for(int chr=0;chr<num_chr;chr++){
			System.out.println(regions[chr][0]+":"+regions[chr][1]+":"+regions[chr][2]+":"+var_comp[chr]+":");
		}
		return regions;
	}
	
	public static void main(String[] args){
		String geno_folder="/Users/quanlong/Documents/projects2/predictions/RA_data/genotypes/";
		String local_K_folder=geno_folder+"local_RRM/";
		String[] types={"etanercept","infliximab","adalimumab","all"};
		int win=100000;
		int[][][] regions={
				//"etanercept"
				{{0,160500001,160500001+win},{10,24350001,24350001+win}},
				//"infliximab"
				{{15,26150001,26150001+win},{11,19950001,19950001+win}},
				// "adalimumab"
				{{16,25500001,25500001+win},{2,82650001,82650001+win}},
				// "all"
				{{12,95250001,95250001+win},{13,101500001,101500001+1}}
		};
		for(int k=0;k<types.length;k++){	
			System.out.println("==============="+types[k]+"=================");
			VariantsDouble obs_geno=new VariantsDouble(geno_folder+types[k]+".train.hdf5");
			VariantsDouble pred_geno=new VariantsDouble(geno_folder+types[k]+".test.hdf5");
			int[][] the_regions=regions[k];
			String[] reginal_kinship_files=new String[the_regions.length];
			String[] relation_matrix_files=new String[the_regions.length];
			for(int i=0;i<the_regions.length;i++){
				reginal_kinship_files[i]=local_K_folder+types[k]+".train."+(the_regions[i][0]+1)+"."+
						the_regions[i][1]+"."+the_regions[i][2]+".K.RRM";
				relation_matrix_files[i]=local_K_folder+types[k]+".train.test."+(the_regions[i][0]+1)+"."+
						the_regions[i][1]+"."+the_regions[i][2]+".K.RRM";
			}
			CausalRegions cal=new CausalRegions(the_regions, null, reginal_kinship_files, relation_matrix_files,  
					null, null, obs_geno, pred_geno);
			
			
			
		}
	}
}

