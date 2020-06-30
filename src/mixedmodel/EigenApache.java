package mixedmodel;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;

//import cern.colt.matrix.linalg.EigenvalueDecomposition;

//import jama.EigenvalueDecomposition;
//import jama.Matrix;

public class EigenApache {
	public double[] eigenvalues;
	public double[][] eigenvectors_asrows;
	
	public EigenApache(double[][] array){
		for(int i=0;i<array.length;i++){
			for(int j=i+1;j<array.length;j++){
				array[i][j]=array[j][i];
			}
		}
		
		Array2DRowRealMatrix matrix= new Array2DRowRealMatrix(array);
		
		EigenDecomposition eigen=new EigenDecomposition(matrix, 0);
		
//		EigenvalueDecomposition eigen =new EigenvalueDecomposition(matrix);
		this.eigenvalues=eigen.getRealEigenvalues().clone();//matrix.getSortedEigenValues().clone();
		this.eigenvectors_asrows=eigen.getVT().getData();
//		int n=this.eigenvalues.length;
//		for(int i=0;i<n/2;i++){
//			double tmp=this.eigenvalues[i];
//			this.eigenvalues[i]=this.eigenvalues[n-i-1];
//			this.eigenvalues[n-i-1]=tmp;
//			double[] tmp2=this.eigenvectors_asrows[i];
//			this.eigenvectors_asrows[i]=this.eigenvectors_asrows[n-i-1];
//			this.eigenvectors_asrows[n-i-1]=tmp2;
//		}
		matrix=null; // release memory
	}
	
	/*
	 * modification according to EMMA/emma.eigen.R.wo.Z(K,X)
	 */
	public void modify4emma(int num){
		double[] updated_values = new double[num];
		double[][] updated_vectors = new double[num][this.eigenvectors_asrows[0].length];
		for(int i=0;i<num;i++){
			updated_values[i]=this.eigenvalues[i]-1;
			updated_vectors[i]=this.eigenvectors_asrows[i].clone();
		}
		this.eigenvalues=updated_values;
		this.eigenvectors_asrows=updated_vectors;
	}
}
