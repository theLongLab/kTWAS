package myMathLib;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import javax.swing.text.StyleContext.SmallAttributeSet;

import org.apache.commons.math3.stat.clustering.Cluster;
import org.apache.commons.math3.stat.clustering.EuclideanIntegerPoint;
import org.apache.commons.math3.stat.clustering.KMeansPlusPlusClusterer;

public class ClusterApacheKPP {

	public int num_clusters;
	public double[][] center_samples;
	public int[] number_of_samples_in_cluster;
	public int[] indexes_in_genomes ; // IDs of genomes that are closest to the centers
	
	/*
	 * constructor, which actually will do all the computations
	 */
	public ClusterApacheKPP(double[][] genomes, int k_small, int k_big, int numTrials, int maxIterationsPerTrial,
			double scale){		
		double[] dbi=new double[k_big-k_small+1];
		double[][][] candidate_centroids=new double[k_big-k_small+1][][];
		int[][] num_samples_in_cluster=new int[k_big-k_small+1][];
			
		KMeansPlusPlusClusterer<EuclideanIntegerPoint> analyzer=new KMeansPlusPlusClusterer(new Random(1));
		ArrayList<EuclideanIntegerPoint> data=new ArrayList<EuclideanIntegerPoint>();
		int num_site=genomes[0].length;
		//fill data
		for(int genome_index=0;genome_index<genomes.length;genome_index++){
			int[] the_genome=new int[num_site];
			for(int j=0;j<num_site;j++)the_genome[j]=(int)(genomes[genome_index][j]*scale);
			data.add(new EuclideanIntegerPoint(the_genome));
		}// try clustering			
		int k_big_OK=k_small;
		for(int K=k_small;K<=k_big;K++){
			List<Cluster<EuclideanIntegerPoint>> result= analyzer.cluster(data, K, numTrials, maxIterationsPerTrial);
			candidate_centroids[K-k_small]=closest_samples_to_centroids(result);
			if(candidate_centroids[K-k_small]==null){
				k_big_OK=K-1;break;
			}
			dbi[K-k_small]=dbi(result);			
			num_samples_in_cluster[K-k_small]=new int[K];
			for(int k=0;k<K;k++)num_samples_in_cluster[K-k_small][k]=result.get(k).getPoints().size();
//			System.out.println(K+"===\t"+dbi[K-k_small]);
		}// see which one is the best K
		if(k_big_OK<k_small){
			System.out.println("Only Smaller number of clusters than specified lowest: too little variations.");
			this.num_clusters=k_big_OK;
			return;
		}
		int best_K=k_small;
		for(int the_k=k_small+1;the_k<=k_big_OK;the_k++){
			if(dbi[the_k-k_small]<dbi[best_K-k_small]){
				best_K=the_k;
			}
		}
		this.num_clusters=best_K;
		this.indexes_in_genomes=new int[best_K];
		this.number_of_samples_in_cluster=num_samples_in_cluster[best_K-k_small];
		this.center_samples=new double[candidate_centroids[best_K-k_small].length][];
		for(int i=0;i<this.center_samples.length;i++){
			this.center_samples[i]=candidate_centroids[best_K-k_small][i].clone();
			for(int j=0;j<this.center_samples[i].length;j++){
				this.center_samples[i][j]=this.center_samples[i][j]/scale;
			}
		}
//		myMathLib.Test.outputArray(this.centers);
		// find out who is who!!
		for(int K=0;K<this.num_clusters;K++){
			this.indexes_in_genomes[K]=0;
			double min_diff=distance(genomes[0], this.center_samples[K]);
			for(int j=1;j<genomes.length;j++){
				double current_diff=distance(genomes[j], this.center_samples[K]);
				if(min_diff>current_diff){
					min_diff=current_diff;
					this.indexes_in_genomes[K]=j;
				}
			}
//			System.out.println(min_diff+"\t"+this.indexes_in_genomes[K]+"/"+this.number_of_samples_in_cluster[K]);
		}
	}
	
	/*
	 * constructor, which actually will do all the computations: clusters with members smaller than min_count will be removed. 
	 */
	public ClusterApacheKPP(double[][] genomes, int k_small, int k_big, int numTrials, int maxIterationsPerTrial,
			double scale, int min_count){	
		this(genomes, k_small, k_big, numTrials, maxIterationsPerTrial, scale);
		if(this.num_clusters<k_small)return;
		ArrayList<Integer> qualified_clusters=new ArrayList<Integer>();
		for(int k=0;k<this.num_clusters;k++){
			if(this.number_of_samples_in_cluster[k]>=min_count){
				qualified_clusters.add(k);
			}
		}
		int new_num_clusters=qualified_clusters.size();
		double[][] new_center_samples=new double[new_num_clusters][];
		int[] new_number_of_samples_in_cluster=new int[new_num_clusters];
		int[] new_indexes_in_genomes=new int[new_num_clusters];
		for(int i=0;i<new_num_clusters;i++){
			int index=qualified_clusters.get(i);
			new_center_samples[i]=this.center_samples[index].clone();
			new_number_of_samples_in_cluster[i]=this.number_of_samples_in_cluster[index];
			new_indexes_in_genomes[i]=this.indexes_in_genomes[index];
		}
		this.center_samples=new_center_samples;
		this.indexes_in_genomes=new_indexes_in_genomes;
		this.num_clusters=new_num_clusters;
		this.number_of_samples_in_cluster=new_number_of_samples_in_cluster;
	}
	
	/*
	 * Davies-Bouldin index: Davies, D. L.; Bouldin, D. W. (1979). 
	 * "A Cluster Separation Measure". IEEE Transactions on Pattern Analysis and Machine Intelligence (2): 224. 
	 * DOI:10.1109/TPAMI.1979.4766909
	 */
	public static double dbi(double[][] centroids, double[][][] points){
		double DBI=0;
		int num_clusters=centroids.length;
		double[] Si=new double[num_clusters];
		double[][] Rij=new double[num_clusters][num_clusters];
		double[] Ri=new double[num_clusters];
		for(int i=0;i<num_clusters;i++)Si[i]=distance_from_centroid(centroids[i], points[i]);		
		for(int i=0;i<num_clusters;i++){			
			for(int j=0;j<num_clusters;j++){
				if(i!=j)Rij[i][j]=(Si[i]+Si[j])/distance(centroids[i],centroids[j]);
			}
		}
		for(int i=0;i<num_clusters;i++){
			double max=0;
			for(int j=0;j<num_clusters;j++){
				if(max<Rij[i][j])max=Rij[i][j];
			}Ri[i]=max;
		}
		for(int c=0;c<num_clusters;c++)DBI+=Ri[c];
		return DBI/num_clusters;
	}
	
	/*
	 * Transform data to make use of above double dbi(double[][] centroids, double[][][] points);
	 */
	public static double dbi(List<Cluster<EuclideanIntegerPoint>> result){
		double[][] centroids=new double[result.size()][];
		double[][][] points=new double[result.size()][][];
		for(int k=0;k<result.size();k++){
			Cluster<EuclideanIntegerPoint> the_cluster=result.get(k);
			int[] the_center=the_cluster.getCenter().getPoint();
			centroids[k]=new double[the_center.length];
			for(int i=0;i<the_center.length;i++)centroids[k][i]=the_center[i];
	
			List<EuclideanIntegerPoint> int_points = the_cluster.getPoints();
			points[k]=new double[int_points.size()][];
			for(int i=0;i<int_points.size();i++){
				int[] the_point=int_points.get(i).getPoint();
				points[k][i]=new double[the_point.length];
				for(int j=0;j<the_point.length;j++)points[k][i][j]=the_point[j];
			}
		}		
		return dbi(centroids, points);
	}
	
	public static double[][] closest_samples_to_centroids(List<Cluster<EuclideanIntegerPoint>> result){
		double[][] closest_centroids=new double[result.size()][];
		for(int K=0;K<result.size();K++){
			Cluster<EuclideanIntegerPoint> the_cluster=result.get(K);
			int[] the_center=the_cluster.getCenter().getPoint();
			List<EuclideanIntegerPoint> all_the_points = the_cluster.getPoints();
			if(all_the_points.size()==0)return null;
			double[] diffs=new double[all_the_points.size()];
			for(int i=0;i<all_the_points.size();i++){
				int[] the_point=all_the_points.get(i).getPoint();
				diffs[i]=distance(the_center, the_point);
			}
			int min_diff_index=0;
			for(int i=0;i<all_the_points.size();i++){
				if(diffs[min_diff_index]>diffs[i])min_diff_index=i;
			}
			int[] closest_centroids_int=all_the_points.get(min_diff_index).getPoint();
			closest_centroids[K]=new double[closest_centroids_int.length];		
			for(int j=0;j<closest_centroids_int.length;j++)closest_centroids[K][j]=closest_centroids_int[j];			
		}		
		return closest_centroids;
	}
	
	public static double distance_from_centroid(double[] centroid, double[][] points){
		int num_points=points.length;
		double dist=0;
		for(int k=0;k<num_points;k++){
			dist+=sqr_distance(centroid, points[k]);
		}
		return Math.sqrt(dist/num_points);
	}
	
	public static double distance(double[] v1, double[] v2){
		if(v1.length!=v2.length){
			System.out.println("Error: v1.length!=v2.length");
			return -1;
		}double dis=0;
		for(int i=0;i<v1.length;i++){
			dis=dis+(v1[i]-v2[i])*(v1[i]-v2[i]);
		}
		return Math.sqrt(dis);
	}
	
	public static double distance(int[] v1, int[] v2){
		if(v1.length!=v2.length){
			System.out.println("Error: v1.length!=v2.length");
			return -1;
		}double dis=0;
		for(int i=0;i<v1.length;i++){
			dis=dis+(v1[i]-v2[i])*(v1[i]-v2[i]);
		}
		return Math.sqrt(dis);
	}
	
	public static double sqr_distance(double[] v1, double[] v2){
		if(v1.length!=v2.length){
			System.out.println("Error: v1.length!=v2.length");
			return -1;
		}double dis=0;
		for(int i=0;i<v1.length;i++){
			dis=dis+(v1[i]-v2[i])*(v1[i]-v2[i]);
		}
		return dis;
	}
	
//	public static void main(String[] args) {
//		
//		double[][] genomes={{1,1,1,1,1,1,1,1,1,1},{1,1,1,1,1,1,0,1,1,1},{1,1,1,1,1,0,1,1,1,1},{1,1,0,1,1,1,1,1,1,1},
//				{0,0,0,0,1,1,1,1,0,0},{0,0,0,0,1,1,1,1,0,1},{0,0,0,0,1,1,1,1,1,0},{1,0,0,0,1,1,1,1,0,0},
//				{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,1,0,0},{0,0,0,1,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0}};
//		String[] ids=new String[genomes.length];
//		for(int k=0;k<genomes.length;k++)ids[k]=""+(k+1);
//		
//		double[][] c=new double[3][];
//		double[][][] points=new double[3][4][];
//		c[0]=genomes[0].clone();c[1]=genomes[4].clone();c[2]=genomes[8].clone();
//		points[0][0]=genomes[0].clone();points[0][1]=genomes[1].clone();points[0][2]=genomes[2].clone();points[0][3]=genomes[3].clone();
//		points[1][0]=genomes[4].clone();points[1][1]=genomes[5].clone();points[1][2]=genomes[6].clone();points[1][3]=genomes[7].clone();
//		points[2][0]=genomes[8].clone();points[2][1]=genomes[9].clone();points[2][2]=genomes[10].clone();points[2][3]=genomes[11].clone();
//		for(int i=0;i<c.length;i++){
//			for(int k=0;k<c[i].length;k++)c[i][k]=c[i][k]*1000;
//			for(int j=0;j<points[i].length;j++){
//				for(int k=0;k<c[i].length;k++)points[i][j][k]=points[i][j][k]*1000;
//			}
//		}
//		System.out.println(dbi(c,points));
//		
//		
//		int numTrials=1000, maxIterationsPerTrial=2000;
//		ClusterApacheKPP tryit=new ClusterApacheKPP(genomes, 2, 4, numTrials, maxIterationsPerTrial,1);
//		
//		
//		
////		double[] x={0.01,0.05,0.05,0.03,0.02,  0.4,0.3, 0.3,0.5, 0.15, 0.18, 0.14, 0.19,0.22};		
////		KMeansPlusPlusClusterer<EuclideanIntegerPoint> my_cluster=new KMeansPlusPlusClusterer(new Random());		
////		ArrayList<EuclideanIntegerPoint> data=new ArrayList<EuclideanIntegerPoint>();		
////		for(int k=0;k<x.length;k++){
////			int[] the_point=new int[1];
////			the_point[0]=(int)(Math.log10(x[k])*1000.0);
////			data.add(new EuclideanIntegerPoint(the_point));
////		}			
////		List<Cluster<EuclideanIntegerPoint>> result= my_cluster.cluster(data, 2, 5, 10);
////		for(int k=0;k<result.size();k++){
////			Cluster<EuclideanIntegerPoint> the_c=result.get(k);
////			List<EuclideanIntegerPoint> points = the_c.getPoints();
////			EuclideanIntegerPoint center=the_c.getCenter();
////			System.out.println("====");
////			System.out.println(Math.pow(10, center.getPoint()[0]/1000.0)+"\n");
////			
////			for(int i=0;i<points.size();i++)
////				System.out.println(Math.pow(10, points.get(i).getPoint()[0]/1000.0));
////		}
//
//	}

}
