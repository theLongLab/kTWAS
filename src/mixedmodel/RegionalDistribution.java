package mixedmodel;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;

import myMathLib.Normalization;
import myPlotLab.MyManhattanPlot;

import org.apache.commons.math3.util.MathArrays;

public class RegionalDistribution {
	
	public static double[] region_cate={0,0.01,0.02,0.04,0.08,0.16,0.32,0.64};

	public int sample_size;
	public double h0_reml;
	public double[] steps;
	
	public int num_of_region;
	public int[] chr;
	public double[] location;
	public double[][] distribution;
	public double[] regional_contribution;
	public String[] regional_category;
	
	public int chr_num;
	public double[][] chr_contribution; // [chr][0]=VC contribution; [chr][1]=chr_len
	public String[] chr_names;
	
	
	/*
	 * Constructor: the input pvalue file should be a .csv file and it is fine if separated by "," or ", ";
	 * the first line is skipped, and the second line is the header.
	 */
	public RegionalDistribution(String distribution_file){
		try{
			String sep=", ";
			BufferedReader br = new BufferedReader(new FileReader(distribution_file));
			String line=br.readLine();
			if(line.startsWith("#SampleSize="))
			this.sample_size=Integer.parseInt(line.split("; ")[0].split("=")[1]);
			this.h0_reml=Double.parseDouble(line.split("; ")[1].split("=")[1]);
			line=br.readLine();
			String[] steps_string=line.split(sep);
			this.steps=new double[steps_string.length-2];
			for(int i=2;i<steps_string.length;i++)this.steps[i-2]=Double.parseDouble(steps_string[i]);
			line=br.readLine();
			while(line!=null){
				String temp[]=line.split(sep);
				if(temp.length!=this.steps.length+2){
					line=br.readLine();
					continue;
				}
				this.num_of_region++;
				line=br.readLine();
			}
			this.chr=new int[this.num_of_region];
			this.location=new double[this.num_of_region];
			this.distribution=new double[this.num_of_region][this.steps.length];
			br = new BufferedReader(new FileReader(distribution_file));
			line=br.readLine();line=br.readLine();line=br.readLine();
			int p_index=0;
			while(line!=null){
				String temp[]=line.split(sep);
				if(temp.length!=this.steps.length+2){
					line=br.readLine();
					continue;
				}
				this.chr[p_index]=Integer.parseInt(temp[0]);
				this.location[p_index]=Integer.parseInt(temp[1]);
				for(int i=2;i<temp.length;i++){
					String[] data=temp[i].split("/");
					if(data.length!=2){
						break;
					}
					double likelihood=Double.parseDouble(data[0]);
					double h1=Double.parseDouble(data[1]);
					if(Double.isNaN(h1)||Double.isNaN(likelihood)){
						this.distribution[p_index][0]++; // local region contribute 0 to the trait with likelihood h0
					}else{
						double prob=Math.exp(likelihood-this.h0_reml);
						double var=h1*this.steps[i-2];
						boolean found=false;
						for(int tt=1;tt<this.steps.length;tt++){
							if(var<this.steps[tt]){
								this.distribution[p_index][tt-1]+=prob;
								found=true;
								break;
							}
						}if(!found)this.distribution[p_index][this.steps.length-1]+=prob;
					}
				}
				double[] normalized=MathArrays.normalizeArray(this.distribution[p_index], 1.0);
				this.distribution[p_index]=normalized;
				p_index++;
			    line=br.readLine();
			}
			this.chr_num=MyManhattanPlot.lengths_of_chrs.length; //TODO
			this.chr_names=new String[this.chr_num];
			this.chr_contribution=new double[this.chr_num][2];
			for(int t=0;t<chr_num;t++){
				this.chr_names[t]="Chr"+(t+1);
				this.chr_contribution[t][1]=MyManhattanPlot.lengths_of_chrs[t];//TODO
			}
			
			double[] regional_contributions=new double[region_cate.length];
			for(int i=0;i<this.num_of_region;i++){
				for(int k=0;k<this.steps.length;k++){
					this.chr_contribution[this.chr[i]-1][0]+=(this.steps[k]*this.distribution[i][k]);
					//identify a category:
					int the_category=0;
					for(int cate_index=1;cate_index<region_cate.length;cate_index++){
						if(this.steps[k]<region_cate[cate_index]){
							the_category=cate_index-1;
							break;
						}
					}
					regional_contributions[the_category]+=this.distribution[i][k];
				}
			}
			Normalization.normalize_colomns(this.chr_contribution);
			
			this.regional_contribution=MathArrays.normalizeArray(regional_contributions, 100.0);
			this.regional_category=new String[region_cate.length];
			for(int k=0;k<region_cate.length-1;k++)this.regional_category[k]=""+(int)(region_cate[k]*100)+"~"+(int)(region_cate[k+1]*100);
			this.regional_category[region_cate.length-1]=">="+(int)(region_cate[region_cate.length-1]*100);
//			this.regional_category=new String[this.steps.length];
//			for(double k=0;k<this.steps.length;k++)this.regional_category[(int)k]=""+((int)(((k/this.steps.length)*100)*1000))/1000.0;
		}catch(Exception e){e.printStackTrace();}
	}
}




