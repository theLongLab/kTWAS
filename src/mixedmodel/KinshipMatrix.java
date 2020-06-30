package mixedmodel;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import myMathLib.StatFuncs;
import myMathLib.Test;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.optimization.GoalType;
import org.apache.commons.math3.optimization.univariate.BrentOptimizer;
import org.apache.commons.math3.optimization.univariate.UnivariatePointValuePair;

import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import ch.systemsx.cisd.base.mdarray.MDDoubleArray;
import ch.systemsx.cisd.hdf5.HDF5MDDataBlock;
//import flanagan.math.Matrix;
import flanagan.roots.RealRoot;

public class KinshipMatrix {
	public RealMatrix kinship_matrix;
	//public double[][] kinship_matrix_array;
	public String[] ids;
	public HashMap<String, Integer> ids2indexes= new HashMap<>();

	// CHECKED - PATHUM
	public void assign_ids2indexes() {
		for (int k = 0; k < ids.length; k++) this.ids2indexes.put(ids[k], k);
	}
	
	public KinshipMatrix(VariantsDouble genotype){					
		double[][] kinship=new double[genotype.sample_size][genotype.sample_size];	
		double[][] total_num_var_used=new double[genotype.sample_size][genotype.sample_size];	
		for(int chr=0;chr<genotype.num_chrs;chr++){				
			for (HDF5MDDataBlock<MDDoubleArray> block : genotype.position_fast_blocks[chr]){
				double[][] data4thisblock=block.getData().toMatrix();
				for(int var_index=0;var_index<data4thisblock.length;var_index++){
					double mean=myMathLib.StatFuncs.mean_NaN(data4thisblock[var_index]);
					if(Double.isNaN(mean)){
						System.out.println("All subjects are NaN at a SNP in Chr"+(chr+1));
						continue;
					}
					double sd=Math.sqrt(myMathLib.StatFuncs.var_NaN(data4thisblock[var_index],mean));
					double[] new_snp=new double[genotype.sample_size];
					for(int k=0;k<genotype.sample_size;k++)new_snp[k]=(data4thisblock[var_index][k]-mean)/sd;
					for(int i=0;i<genotype.sample_size;i++){
						for(int j=i;j<genotype.sample_size;j++){
							if((!Double.isNaN(new_snp[i])) && (!Double.isNaN(new_snp[j]))){
								kinship[i][j]=kinship[i][j]+(new_snp[i]*new_snp[j]);
								total_num_var_used[i][j]++;
							}
						}
					}
				}
			}
			System.out.println("Finished Chr"+(chr+1));
		}			
		for(int i=0;i<genotype.sample_size;i++){
			for(int j=i;j<genotype.sample_size;j++){
				kinship[i][j]=kinship[i][j]/total_num_var_used[i][j];
				kinship[j][i]=kinship[i][j];
			}
		}
		this.kinship_matrix=new Array2DRowRealMatrix(kinship, false);
		this.ids=genotype.sample_ids.clone();
		this.assign_ids2indexes();
	}
	
	public KinshipMatrix(VariantsDouble genotype, double min_MAF){					
		double[][] kinship=new double[genotype.sample_size][genotype.sample_size];	
		double[][] total_num_var_used=new double[genotype.sample_size][genotype.sample_size];	
		for(int chr=0;chr<genotype.num_chrs;chr++){				
			for (HDF5MDDataBlock<MDDoubleArray> block : genotype.position_fast_blocks[chr]){
				double[][] data4thisblock=block.getData().toMatrix();
				for(int var_index=0;var_index<data4thisblock.length;var_index++){
					if(maf(data4thisblock[var_index],2)<min_MAF)continue;
					double mean=myMathLib.StatFuncs.mean_NaN(data4thisblock[var_index]);
					if(Double.isNaN(mean)){
						System.out.println("All subjects are NaN at a SNP in Chr"+(chr+1));
						continue;
					}
					double sd=Math.sqrt(myMathLib.StatFuncs.var_NaN(data4thisblock[var_index],mean));
					double[] new_snp=new double[genotype.sample_size];
					for(int k=0;k<genotype.sample_size;k++)new_snp[k]=(data4thisblock[var_index][k]-mean)/sd;
					for(int i=0;i<genotype.sample_size;i++){
						for(int j=i;j<genotype.sample_size;j++){
							if((!Double.isNaN(new_snp[i])) && (!Double.isNaN(new_snp[j]))){
								kinship[i][j]=kinship[i][j]+(new_snp[i]*new_snp[j]);
								total_num_var_used[i][j]++;
							}
						}
					}
				}
			}
			System.out.println("Finished Chr"+(chr+1));
		}			
		for(int i=0;i<genotype.sample_size;i++){
			for(int j=i;j<genotype.sample_size;j++){
				kinship[i][j]=kinship[i][j]/total_num_var_used[i][j];
				kinship[j][i]=kinship[i][j];
			}
		}
		this.kinship_matrix=new Array2DRowRealMatrix(kinship, false);
		this.ids=genotype.sample_ids.clone();
		this.assign_ids2indexes();
	}
	
	/*
	 * int chr=local[0], start=local[1], end= local[2];
	 */
	public KinshipMatrix(VariantsDouble genotype, int[] local){					
		double[][] kinship=new double[genotype.sample_size][genotype.sample_size];	
		double[][] total_num_var_used=new double[genotype.sample_size][genotype.sample_size];	
		int chr=local[0], start=local[1], end= local[2];
		double[][] data4thisblock=genotype.load_variants_in_region(chr, start, end);
		for(int var_index=0;var_index<data4thisblock.length;var_index++){
			double mean=myMathLib.StatFuncs.mean_NaN(data4thisblock[var_index]);
			if(Double.isNaN(mean)){
				System.out.println("All subjects are NaN at a SNP in Chr"+(chr+1));
				continue;
			}
			double sd=Math.sqrt(myMathLib.StatFuncs.var_NaN(data4thisblock[var_index],mean));
			double[] new_snp=new double[genotype.sample_size];
			for(int k=0;k<genotype.sample_size;k++)new_snp[k]=(data4thisblock[var_index][k]-mean)/sd;
			for(int i=0;i<genotype.sample_size;i++){
				for(int j=i;j<genotype.sample_size;j++){
					if((!Double.isNaN(new_snp[i])) && (!Double.isNaN(new_snp[j]))){
						kinship[i][j]=kinship[i][j]+(new_snp[i]*new_snp[j]);
						total_num_var_used[i][j]++;
					}
				}
			}
		}	
		System.out.println("Finished kinship for region Chr"+(chr+1)+":"+start+":"+end);
		for(int i=0;i<genotype.sample_size;i++){
			for(int j=i;j<genotype.sample_size;j++){
				kinship[i][j]=kinship[i][j]/total_num_var_used[i][j];
				kinship[j][i]=kinship[i][j];
			}
		}
		this.kinship_matrix=new Array2DRowRealMatrix(kinship, false);
		this.ids=genotype.sample_ids.clone();
		this.assign_ids2indexes();
	}
	
	/*
	 * compound[i] encode the i_th region 
	 * int chr=compound[r_index][0], start=compound[r_index][1], end= compound[r_index][2];
	 */
	public KinshipMatrix(VariantsDouble genotype, int[][] compound){					
		double[][] kinship=new double[genotype.sample_size][genotype.sample_size];	
		double[][] total_num_var_used=new double[genotype.sample_size][genotype.sample_size];	
		
		for(int r_index=0;r_index<compound.length;r_index++){
			int chr=compound[r_index][0], start=compound[r_index][1], end= compound[r_index][2];
			double[][] data4thisblock=genotype.load_variants_in_region(chr, start, end);
			for(int var_index=0;var_index<data4thisblock.length;var_index++){
				double mean=myMathLib.StatFuncs.mean_NaN(data4thisblock[var_index]);
				if(Double.isNaN(mean)){
					System.out.println("All subjects are NaN at a SNP in Chr"+(chr+1));
					continue;
				}
				double sd=Math.sqrt(myMathLib.StatFuncs.var_NaN(data4thisblock[var_index],mean));
				double[] new_snp=new double[genotype.sample_size];
				for(int k=0;k<genotype.sample_size;k++)new_snp[k]=(data4thisblock[var_index][k]-mean)/sd;
				for(int i=0;i<genotype.sample_size;i++){
					for(int j=i;j<genotype.sample_size;j++){
						if((!Double.isNaN(new_snp[i])) && (!Double.isNaN(new_snp[j]))){
							kinship[i][j]=kinship[i][j]+(new_snp[i]*new_snp[j]);
							total_num_var_used[i][j]++;
						}
					}
				}
			}	
		}
		for(int i=0;i<genotype.sample_size;i++){
			for(int j=i;j<genotype.sample_size;j++){
				kinship[i][j]=kinship[i][j]/total_num_var_used[i][j];
				kinship[j][i]=kinship[i][j];
			}
		}
		System.out.print("Finished kinship for the compound ");
		for(int r_index=0;r_index<compound.length;r_index++){
			int chr=compound[r_index][0], start=compound[r_index][1], end= compound[r_index][2];
			System.out.print("R"+r_index+":"+(chr+1)+":"+start+":"+end+"; ");
		}System.out.println();
		this.kinship_matrix=new Array2DRowRealMatrix(kinship, false);
		this.ids=genotype.sample_ids.clone();
		this.assign_ids2indexes();
	}
	
	public KinshipMatrix(RealMatrix kinship_matrix, String[] ids){
		this.kinship_matrix=kinship_matrix;
		this.ids=ids;
		this.assign_ids2indexes();
		//this.kinship_matrix_array=this.kinship_matrix.getData();
	}
	
	public KinshipMatrix(double[][] kinship_matrix, String[] ids){
		//this.kinship_matrix_array=kinship_matrix;
		this.kinship_matrix=new Array2DRowRealMatrix(kinship_matrix, false);
		this.ids=ids;
		this.assign_ids2indexes();
	}
	
	/*
	 * reads kinship from file, while assuming matched ids are specified in parameter. 
	 */
	public KinshipMatrix(String kinship_file, String[] ids){
		double[][] kinship=new double[ids.length][ids.length];
		String[] separators={","," ","\t"};
		String sep=null;
		try{			
			BufferedReader br= new BufferedReader(new FileReader(kinship_file));
			String line=br.readLine();//first line
			for(int k=0;k<separators.length;k++){
				if(line.split(separators[k]).length==ids.length){
					sep=separators[k];
				}
			}			
			int line_index=0;
			while(line!=null){
				String[] temp=line.split(sep);
				for(int k=0;k<temp.length;k++){
					kinship[line_index][k]=Double.parseDouble(temp[k]);
				}line_index++;
				line=br.readLine();
			}
			if(line_index!=ids.length){
				System.out.println("Kinship file number of line wrong!");
			}
			this.ids=ids;
			this.kinship_matrix=new Array2DRowRealMatrix(kinship, false);
			//this.kinship_matrix_array=kinship;
			this.assign_ids2indexes();
		}catch(Exception e){e.printStackTrace();}		
	}
	
	
	/*
	 * reads kinship from file, while assuming the ids are the first line of the file. 
	 */
	// CHECKED - PATHUM
	public KinshipMatrix(String kinship_file) {
		String[] separators = {",", " ", "\t"};
		String sep = null;
		try {
			BufferedReader br = new BufferedReader(new FileReader(kinship_file));
			String line = br.readLine();//header
			int sample_size = 0;
			for (int k = 0; k < separators.length; k++) {	// Identifies the separator as the separator which causes
				// the most number of splits
				int the_len = line.split(separators[k]).length;
				if (the_len > sample_size) {
					sep = separators[k];
					sample_size = the_len;
				}
			}

			double[][] kinship = new double[sample_size][sample_size];
			this.ids = line.split(sep);
			line = br.readLine();
			int line_index = 0;
			while (line != null) {
				String[] temp = line.split(sep);
				for (int k = 0; k < temp.length; k++) {
					kinship[line_index][k] = Double.parseDouble(temp[k]);
				}
				line_index++;
				line = br.readLine();
			}
			if (line_index != ids.length) {
				System.out.println("There appears to be something wrong with the number of lines in this kinship file!");
			}
			this.kinship_matrix = new Array2DRowRealMatrix(kinship, false);
			//this.kinship_matrix_array=kinship;
			this.assign_ids2indexes();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	

	/*
	 * S_new=(n-1)S/Tr(PSP), where P= I-11'/n and 1 is the vector of ones.  
	 */
	// TODO: Find reference for calculation
	public void re_scale_kinship_matrix() {
		double[][] S = this.kinship_matrix.getData();
		int n = S.length;
		double[][] S_new = new double[n][n];
		DenseDoubleMatrix2D matrix_P = new DenseDoubleMatrix2D(n, n);

		matrix_P.assign(-1.0 / n);

		for (int i = 0; i < n; i++) {
			matrix_P.setQuick(i, i, (1.0 - 1.0 / n));
		}

		DenseDoubleMatrix2D matrix_S = new DenseDoubleMatrix2D(S);
		DenseDoubleMatrix2D PS = new DenseDoubleMatrix2D(n, n);
		DenseDoubleMatrix2D PSP = new DenseDoubleMatrix2D(n, n);
		matrix_P.zMult(matrix_S, PS);
		PS.zMult(matrix_P, PSP);

		double trace = 0;
		for (int k = 0; k < n; k++) trace += PSP.get(k, k);
		System.out.println("Rescale coef= " + (n - 1) / trace);

		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				S_new[i][j] = S[i][j] * (n - 1) / trace;
			}
		}

		this.kinship_matrix = new Array2DRowRealMatrix(S_new, false);
	}
	
	/*
	 * S_new=(n-1)S/Tr(PSP), where P= I-11'/n and 1 is the vector of ones.  
	 */
	public static double[][] re_scale_kinship_matrix(double[][] S){
		//double[][] S=this.kinship_matrix.getData();
		int n=S.length;
		double[][] S_new=new double[n][n];		
		DenseDoubleMatrix2D matrix_P= new DenseDoubleMatrix2D(n,n);
		matrix_P.assign(-1.0/n);
		for(int i=0;i<n;i++){
			matrix_P.setQuick(i, i, (1.0-1.0/n));
		}
		DenseDoubleMatrix2D matrix_S= new DenseDoubleMatrix2D(S);
		DenseDoubleMatrix2D PS= new DenseDoubleMatrix2D(n,n);
		DenseDoubleMatrix2D PSP= new DenseDoubleMatrix2D(n,n);
		matrix_P.zMult(matrix_S,PS);
		PS.zMult(matrix_P, PSP);
		double trace=0;
		for(int k=0;k<n;k++)	trace+=PSP.get(k, k);
//		System.out.println("Rescale coef="+(n-1)/trace);
		for(int i=0;i<n;i++){
			for(int j=0;j<n;j++){
				S_new[i][j]=S[i][j]*(n-1)/trace;
			}
		}
		return S_new;
	}
	
	
	/*
	 * S_new=(n-1)S/Tr(PSP), where P= I-11'/n and 1 is the vector of ones.  
	 */
	public static RealMatrix re_scale_kinship_matrix(RealMatrix matrix){
		double[][] S=matrix.getData();
		int n=S.length;
		double[][] S_new=new double[n][n];		
		DenseDoubleMatrix2D matrix_P= new DenseDoubleMatrix2D(n,n);
		matrix_P.assign(-1.0/n);
		for(int i=0;i<n;i++){
			matrix_P.setQuick(i, i, (1.0-1.0/n));
		}
		DenseDoubleMatrix2D matrix_S= new DenseDoubleMatrix2D(S);
		DenseDoubleMatrix2D PS= new DenseDoubleMatrix2D(n,n);
		DenseDoubleMatrix2D PSP= new DenseDoubleMatrix2D(n,n);
		matrix_P.zMult(matrix_S,PS);
		PS.zMult(matrix_P, PSP);
		double trace=0;
		for(int k=0;k<n;k++)	trace+=PSP.get(k, k);
		//System.out.println("Rescale coef="+(n-1)/trace);
		for(int i=0;i<n;i++){
			for(int j=0;j<n;j++){
				S_new[i][j]=S[i][j]*(n-1)/trace;
			}
		}
		return (new Array2DRowRealMatrix(S_new, false));
	}
	
	/*
	 * S_new=(n-1)S/Tr(PSP), where P= I-11'/n and 1 is the vector of ones.  
	 */
	public static RealMatrix re_scale_kinship_matrix(RealMatrix matrix, boolean print_coefficient){
		double[][] S=matrix.getData();
		int n=S.length;
		double[][] S_new=new double[n][n];		
		DenseDoubleMatrix2D matrix_P= new DenseDoubleMatrix2D(n,n);
		matrix_P.assign(-1.0/n);
		for(int i=0;i<n;i++){
			matrix_P.setQuick(i, i, (1.0-1.0/n));
		}
		DenseDoubleMatrix2D matrix_S= new DenseDoubleMatrix2D(S);
		DenseDoubleMatrix2D PS= new DenseDoubleMatrix2D(n,n);
		DenseDoubleMatrix2D PSP= new DenseDoubleMatrix2D(n,n);
		matrix_P.zMult(matrix_S,PS);
		PS.zMult(matrix_P, PSP);
		double trace=0;
		for(int k=0;k<n;k++)	trace+=PSP.get(k, k);
		System.out.println("Rescale coef="+(n-1)/trace);
		for(int i=0;i<n;i++){
			for(int j=0;j<n;j++){
				S_new[i][j]=S[i][j]*(n-1)/trace;
			}
		}
		return (new Array2DRowRealMatrix(S_new, false));
	}

	public static void re_scale_kinship_matrix(String in_kinship_file, String out_kinship_file) {
		KinshipMatrix old = new KinshipMatrix(in_kinship_file);
		old.re_scale_kinship_matrix();
		old.write2file(out_kinship_file);
	}
	
	public void write2file(String file){
		try{
			int sample_size=this.ids.length;
			BufferedWriter bw=new BufferedWriter(new FileWriter(file));
			for(int k=0;k<sample_size-1;k++){
				bw.write(ids[k]+",");
			}bw.write(ids[sample_size-1]+"\n");
			double[][] kinship=this.kinship_matrix.getData();
			for(int i=0;i<sample_size;i++){
				for(int j=0;j<sample_size-1;j++){
					bw.write(kinship[i][j]+",");
				}bw.write(kinship[i][sample_size-1]+"\n");
			}
			bw.close();
		}catch(Exception e){e.printStackTrace();}	
	}
	
	/*
	 * one has to make sure the dim matches between kinship and Y before running.
	 */
	public static RealMatrix find_transformation_matrix_fastlmm(double[][] kinship, double[] Y){
		if(kinship.length!=Y.length){
			System.out.println("Wrong: dim not match!");
			return null;
		}
		int sample_size=Y.length;
		double[][] transformed_X=new double[1][];
		double[] intercept_1=new double[sample_size];
		Arrays.fill(intercept_1, 1);	
		
		EigenDecomposition eigen=new EigenDecomposition(new Array2DRowRealMatrix(kinship, false),0); 
		RealMatrix global_Ut=eigen.getVT();
		transformed_X[0]=global_Ut.operate(intercept_1);
		double[] Y_K_transformed_null=global_Ut.operate(Y);
		double[] S=eigen.getRealEigenvalues();
		
		BrentOptimizer bo=new BrentOptimizer(FaSTLMM.brent_rel, FaSTLMM.brent_abs);
		double best_delta=-1;
		double best_ml=-(Double.MAX_VALUE);
		UnivariateFunction ml_function=new MLFullRank4BrentOptimal(transformed_X,Y_K_transformed_null,S);
		for(double min=FaSTLMM.lowlimit;min<FaSTLMM.uplimit;min+=FaSTLMM.grid_step){
			double min_delta=Math.exp(min);
			double max_delta=Math.exp(min+FaSTLMM.grid_step);
			UnivariatePointValuePair result=bo.optimize(FaSTLMM.maxEval, ml_function, GoalType.MAXIMIZE, min_delta, max_delta);
			double the_ml=result.getValue();
			if(the_ml>best_ml){
				best_ml=the_ml;
				best_delta=result.getPoint();
			}
		}		
		System.out.println(best_delta);
		double[][] transformation_matrix_array= new double[sample_size][sample_size];
		double[] diag= new double[sample_size];
		for(int i=0;i<sample_size;i++){
	//		if(this.reml_eig_L.eigenvalues[i]>0)
				diag[i]=1.0/Math.sqrt(S[i]+best_delta);
		}
		for(int i=0;i<sample_size;i++){
			for(int j=0;j<sample_size;j++){
				transformation_matrix_array[i][j]=global_Ut.getEntry(i,j)*diag[i]; 
			}
		}
		RealMatrix transformation_matrix=new Array2DRowRealMatrix(transformation_matrix_array, false);		
		return transformation_matrix;
	}
	
	public void check_normality(){
		int sample_size=this.ids.length;
		double[] data=new double[sample_size*(sample_size-1)/2];
		int index=0;
		for(int i=1;i<sample_size;i++){
			for(int j=i+1;j<sample_size;j++){
				data[index]=this.kinship_matrix.getEntry(i, j);
				index++;
			}
		}
		System.out.println(StatFuncs.normality_test_ks(data));
	}
	
	/*
	 * subset_ids should be a subset of all the ids.
	 */
	public static KinshipMatrix submatrix(KinshipMatrix full_matrix, String[] subset_ids){
		double[][] data=new double[subset_ids.length][subset_ids.length];
		double[][] full_data=full_matrix.kinship_matrix.getData();
		for(int i=0;i<subset_ids.length;i++){
			if(!full_matrix.ids2indexes.containsKey(subset_ids[i])){
				System.out.println("Can't substract kinship since "+subset_ids[i]+" isn't there!");
				return null;
			}
			int index_i=full_matrix.ids2indexes.get(subset_ids[i]);
			for(int j=0;j<subset_ids.length;j++){
				if(!full_matrix.ids2indexes.containsKey(subset_ids[j])){
					System.out.println("KinshipMatrix: Can't substract kinship since "+subset_ids[j]+" isn't there!");
					return null;
				}
				int index_j=full_matrix.ids2indexes.get(subset_ids[j]);
				data[i][j]=full_data[index_i][index_j];
			}
		}
		return new KinshipMatrix(data, subset_ids);
	}
	
	public static double maf(double[] snp, double scale){
		double counts=0, non_NaN=0;
		for(int k=0;k<snp.length;k++){
			if(!Double.isNaN(snp[k])){
				counts=counts+snp[k];
				non_NaN=non_NaN+scale;
			}
		}
		double maf=counts/non_NaN;
		if(maf>=0.5)maf=1-maf;
		return maf;
	}
}


