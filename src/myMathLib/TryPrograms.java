package myMathLib;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.ListIterator;
import java.util.Random;
import java.util.Set;
import java.util.SortedSet;

import mixedmodel.EMMAX;
import mixedmodel.KinshipMatrix;
import mixedmodel.LocalKinshipAnalyzer;
import mixedmodel.VariantsDouble;

import org.apache.commons.math3.stat.StatUtils;





public class TryPrograms {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
//		double[] data={0.03077972, -0.66081015,  0.21684005,
//				-1.41882095,  0.13929080, -0.79641312,  1.17466404,  0.76569981, -0.56422328, -0.22769212};
//		System.out.println(StatUtils.mean(data));
//		System.out.println(StatUtils.variance(data));
		
//		System.out.println(StatFuncs.pvalue_spearman(-0.17575757575,10));
		
//		String input_csv="/Volumes/Projects/Jawamix5/debug/input.csv";
//		VariantsDouble.importCSV(input_csv, input_csv+".hdf", 10);
		
//		EMMAX.make_plot_one_phen("/Volumes/Projects/Local-kinship/simulations/have2look/176/EMMAX.0_176_Leaf_roll_10.top",0);
//		LocalKinshipAnalyzer.make_plot_one_phen("/Volumes/Projects/Local-kinship/simulations/have2look/Local_VO.0.176_Leaf_roll_10.w100000.again.csv");
		
//		Random generator = new Random();
//		for(int i=0;i<10;i++)System.out.println(generator.nextDouble());
//		System.out.println();
//		for(int i=0;i<10;i++)System.out.println(Test.randomNumber());
		
		
		String kinship_file="/Volumes/Projects/DATA/arabidopsis_genomes/kinship_files_swe180_ecker171_removebads/kinship.swe180_ecker171_removebads.RRMkinship.RRM";
		String kin2="/Volumes/Projects/Collaborations/Pathway_Fang/plink.ibd.RRM";
		KinshipMatrix matrix=new KinshipMatrix(kinship_file);
		
		int sample_size=matrix.ids.length;
		double[] data=new double[sample_size*(sample_size-1)/2];
		int index=0;
		for(int i=1;i<sample_size;i++){
			for(int j=i+1;j<sample_size;j++){
				data[index]=matrix.kinship_matrix.getEntry(i, j);
				index++;
			}
		}try{
			BufferedWriter bw=new BufferedWriter(new FileWriter("/Volumes/Projects/Collaborations/Pathway_Fang/AT.ibd.RRM.single.txt"));
			for(int i=0;i<data.length;i++){
				bw.write(data[i]+"\n");
			}bw.close();
		}catch(Exception e){e.printStackTrace();}
		
		System.out.println("Normality check:");
		matrix.check_normality();
		
	}

}
