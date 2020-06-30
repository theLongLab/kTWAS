package mixedmodel;


import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

import nam.Cross_info;
import nam.Founder;
import nam.Imputation;
import nam.RILs;

import org.apache.commons.lang.ArrayUtils;

import cern.colt.matrix.impl.DenseDoubleMatrix2D;

import ch.systemsx.cisd.base.mdarray.MDByteArray;
import ch.systemsx.cisd.base.mdarray.MDDoubleArray;
import ch.systemsx.cisd.base.mdarray.MDIntArray;
import ch.systemsx.cisd.base.mdarray.MDShortArray;
import ch.systemsx.cisd.hdf5.HDF5CompoundType;
import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.HDF5MDDataBlock;
import ch.systemsx.cisd.hdf5.IHDF5IntReader;
import ch.systemsx.cisd.hdf5.IHDF5IntWriter;
import ch.systemsx.cisd.hdf5.IHDF5Reader;
import ch.systemsx.cisd.hdf5.IHDF5SimpleReader;
import ch.systemsx.cisd.hdf5.IHDF5SimpleWriter;
import ch.systemsx.cisd.hdf5.IHDF5Writer;

public class VariantsByte {
	
	public int sample_size;
	public int num_chrs;
	public int[] num_sites;
	int num_sites_total;
	int block_size=5000;
	int[][] locations;  // chr x num_sites[chr]
//	int[][] mafc;		// chr x num_sites[chr]
	int[] num_blocks_per_chr;
	public String[] sample_ids;
	IHDF5Reader reader;
	Iterable<HDF5MDDataBlock<MDByteArray>>[] position_fast_blocks;
//	ArrayList<HDF5MDDataBlock<MDDoubleArray>>[] block_pointers;
//	int[][] location2block_index;
//	int[][] location2offset_inblock;
	String[] position_fast_paths;
	
	public int num_genes;
	public int[] gene_chrs;
	public int[] gene_starts;
	public int[] gene_ends;
	public String[] gene_names;
	public String[] gene_ENS_ids;
	HashMap<String, Integer> gene_name2index=new HashMap<String, Integer>();
	
	/*
	 * 	ROOT:
	 *  	Group:  Genotype: 
	 *  			Attr: some info: like sample-size, number-of-variants, ... ,
	 *  			Table: Individual ids: 1D array of Strings
	 *  			Table: Gene names
	 *  			Table: Gene starts
	 *  			Table: Gene ends
	 *  			Table: Gene chrs
	 *  			N Groups: for chromosomes
	 *  				Each-groups: Two Byte arrays: individual-fast; location-fast
	 *  							 Table: Variants positions: 1D array of ints
	 *  							 Table: Annotations of the variants: 1D array of HDF5CompoundType<Annotation_g1k>
	 *  					   
	 */
	
	
	/*
	 * constructor: make an hdf5 object based on a .csv file 
	 */
	public static VariantsByte importCSV(String file_csv,String file_hdf5,int block_size){
		System.out.println("Start importing data from "+file_csv+" to HDF5 format.");
		try{
			BufferedReader br= new BufferedReader(new FileReader(file_csv));
			String line=br.readLine();//header
			String[] temp=line.split(",");
			ArrayList<Integer> num_of_sites_per_chr=new ArrayList<Integer>();
			ArrayList<ArrayList<Integer>> var_locartions_list=new ArrayList<ArrayList<Integer>>();
			ArrayList<ArrayList<Integer>> var_mafc_list=new ArrayList<ArrayList<Integer>>();
			int sample_size=temp.length-2;
			String[] sample_ids=new String[sample_size];
			for(int k=0;k<sample_size;k++)sample_ids[k]=temp[k+2];
			line=br.readLine(); // the first data line
			temp=line.split(",");
			int old_chr=Integer.parseInt(temp[0]);
			int num_of_sites_this_chr=1;
			ArrayList<Integer> locations_this_chr=new ArrayList<Integer>();
			locations_this_chr.add(Integer.parseInt(temp[1]));
			ArrayList<Integer> mafc_this_chr=new ArrayList<Integer>();
			mafc_this_chr.add(mafc(temp));
			line=br.readLine();
			while(line!=null){
				temp=line.split(",");
				if(temp.length-2!=sample_size){System.out.println("Sample size not the same!");}
				int chr=Integer.parseInt(temp[0]);
				if(chr!=old_chr){
					var_locartions_list.add(locations_this_chr);
					var_mafc_list.add(mafc_this_chr);
					num_of_sites_per_chr.add(num_of_sites_this_chr);
					locations_this_chr=new ArrayList<Integer>();
					mafc_this_chr=new ArrayList<Integer>();
					num_of_sites_this_chr=0;
					old_chr=chr;
				}
				locations_this_chr.add(Integer.parseInt(temp[1]));
				mafc_this_chr.add(mafc(temp));
				num_of_sites_this_chr++;
				line=br.readLine();
			}// the last chr:
			var_locartions_list.add(locations_this_chr);
			var_mafc_list.add(mafc_this_chr);
			num_of_sites_per_chr.add(num_of_sites_this_chr);
			int num_chrs=num_of_sites_per_chr.size();
			int[] num_sites=new int[num_chrs];
			int[][] var_positions=new int[num_chrs][];
			int[][] var_mafc=new int[num_chrs][];
			for(int chr=0;chr<num_chrs;chr++){
				num_sites[chr]=num_of_sites_per_chr.get(chr);
				var_positions[chr]=new int[num_sites[chr]];
				var_mafc[chr]=new int[num_sites[chr]];
				ArrayList<Integer> locations_the_chr=var_locartions_list.get(chr);
				ArrayList<Integer> mafc_the_chr=var_mafc_list.get(chr);
				for(int k=0;k<num_sites[chr];k++){
					var_positions[chr][k]=locations_the_chr.get(k);
					var_mafc[chr][k]=mafc_the_chr.get(k);
				}
			}
			System.out.println("Finished scanning the .csv file and found "+num_chrs+" chromosomes.\nChr\t#variants");
			for(int chr=0;chr<num_chrs;chr++){
				System.out.println("Chr"+(chr+1)+":\t"+num_sites[chr]);
			}
			File hdf5_old=new File(file_hdf5);
			if(hdf5_old.exists()){hdf5_old.delete();System.out.println(file_hdf5+" alreday exist. Deleted!");}
			IHDF5Writer writer = HDF5Factory.open(file_hdf5);
			// CODE FOR constructing above structure (all groups and tables and importing data)
			writer.createGroup("/genotype");
			writer.setIntAttribute("/genotype", "sample_size", sample_size);
			writer.setIntAttribute("/genotype", "block_size", block_size);
			writer.setIntAttribute("/genotype", "num_chrs", num_chrs);
			writer.writeStringArray("/genotype/sample_ids", sample_ids);
			int[] num_blocks_per_chr=new int[num_chrs];
			for(int chr=0;chr<num_chrs;chr++){
				num_blocks_per_chr[chr]=num_sites[chr]/block_size+1;
				writer.createGroup("/genotype/chr"+(1+chr));
				writer.writeIntArray("/genotype/chr"+(1+chr)+"/var_pos", var_positions[chr]);
				writer.writeIntArray("/genotype/chr"+(1+chr)+"/var_mafc", var_mafc[chr]);
//				writer.createStringArray("/genotype/chr"+(1+chr)+"/annotation", sample_size, num_sites[chr], sample_size, block_size);
				writer.createByteMatrix("/genotype/chr"+(1+chr)+"/indi_fast",sample_size, num_sites[chr], sample_size, block_size);
				writer.createByteMatrix("/genotype/chr"+(1+chr)+"/position_fast", num_sites[chr], sample_size, block_size, sample_size);
			}writer.writeIntArray("/genotype/num_blocks_per_chr", num_blocks_per_chr);					
			// Read data again to fill matrix:
			System.out.println("Reading data again to create HDF5 file for genotypes.");
			br= new BufferedReader(new FileReader(file_csv));
			line=br.readLine();//header
			byte[][] data4thisblock=new byte[block_size][sample_size];
			line=br.readLine();//first line
			temp=line.split(",");		
			int current_chr=Integer.parseInt(temp[0]);
			int total_var_in_chr=0;
			int current_block=0;
			int num_of_variants_in_block=0;
			while(line!=null){
				temp=line.split(",");
				int chr=Integer.parseInt(temp[0]);
				if(chr!=current_chr){
					//output the final block to h5
					writer.writeByteMatrixBlockWithOffset("/genotype/chr"+(current_chr)+"/position_fast",
							data4thisblock, num_of_variants_in_block, sample_size, current_block*block_size, 0L);
					data4thisblock=new byte[block_size][sample_size];
					num_of_variants_in_block=0;
					current_block=0;
					current_chr=chr;
					total_var_in_chr=0;
				}else{//chr==current_chr
					if(num_of_variants_in_block==block_size){
						//output one block to h5
						writer.writeByteMatrixBlock("/genotype/chr"+(current_chr)+"/position_fast",
								data4thisblock, current_block,0);
						data4thisblock=new byte[block_size][sample_size];
						num_of_variants_in_block=0;
						current_block++;
					}
				}
				for(int k=0;k<sample_size;k++){
					data4thisblock[num_of_variants_in_block][k]=Byte.parseByte(temp[k+2]);
				}
				num_of_variants_in_block++;
				total_var_in_chr++;
				line=br.readLine();
			}// write the final block:
			writer.writeByteMatrixBlockWithOffset("/genotype/chr"+(current_chr)+"/position_fast",
					data4thisblock, num_of_variants_in_block, sample_size, current_block*block_size, 0L);
			//write the indexes files		
			System.out.println("Finished writing genotypes. Now the HDF5 data is ready to use.");
			writer.close();			
			
		}catch(Exception e){e.printStackTrace();}		
		return new VariantsByte(file_hdf5);
	}
	
	/*
	 * constructor: make an hdf5 object based on a .csv file and a .info file
	 */
	public static VariantsByte importCSV(String file_csv, String file_info, String file_hdf5, int block_size){
		System.out.println("Start importing data from "+file_csv+" to HDF5 format.");
//		HashMap<String, Integer> gene_name2index=new HashMap<String, Integer>();
		try{
			File hdf5_old=new File(file_hdf5);
			if(hdf5_old.exists()){hdf5_old.delete();System.out.println(file_hdf5+" alreday exist. Deleted!");}
			IHDF5Writer writer = HDF5Factory.open(file_hdf5);
			HDF5CompoundType<Annotation_g1k> type = writer.compounds().getInferredType(Annotation_g1k.class);

			BufferedReader br= new BufferedReader(new FileReader(file_csv));
			String line=br.readLine();//header
			BufferedReader br_info= new BufferedReader(new FileReader(file_info));
			String line_info=br_info.readLine(); //header
			line_info=br_info.readLine();
			ArrayList<String> gene_name_list=new ArrayList<String>();
			ArrayList<String> gene_id_list=new ArrayList<String>();
			ArrayList<Integer> gene_starts_list=new ArrayList<Integer>();
			ArrayList<Integer> gene_ends_list=new ArrayList<Integer>();
			ArrayList<Integer> gene_chrs_list=new ArrayList<Integer>();
			HashMap<String, Integer> gene_name2index_map=new HashMap<String, Integer>();
			int num_of_genes=0;
			String last_gene_name="*";
			while(line_info!=null){
				String[] info=line_info.split("\t");
				if(!info[28].equals("*")){
					String[] vt_annotation=info[28].split(":");
					String gene_name=vt_annotation[1];
					String gene_ENS_id=vt_annotation[2];
					int chr=Integer.parseInt(info[0]);
					int loc=Integer.parseInt(info[1]);
					if(gene_name2index_map.containsKey(gene_name)){
						if(!gene_name.equals(last_gene_name)){
							//System.out.println("???: !gene_name.equals(last_gene_name)");
							//System.out.println(last_gene_name+"\n"+gene_name+"\n"+info[28]);
						}
						int gene_index=gene_name2index_map.get(gene_name);
						if(loc-1<gene_starts_list.get(gene_index)){
							gene_starts_list.set(gene_index,loc-1);
						}if(loc+1>gene_ends_list.get(gene_index)){
							gene_ends_list.set(gene_index,loc+1);
						}if(chr!=gene_chrs_list.get(gene_index)){
							System.out.println("WRONG: chr!=gene_chrs_list.get(gene_index)");
						}
					}else{
						if(gene_name.equals(last_gene_name)){
							System.out.println("WRONG: gene_name.equals(last_gene_name)");
						}
						gene_name2index_map.put(gene_name, num_of_genes);
						gene_name_list.add(gene_name);
						gene_id_list.add(gene_ENS_id);
						gene_chrs_list.add(chr);
						gene_starts_list.add(loc-1);
						gene_ends_list.add(loc+1);
						num_of_genes++;
					}
					last_gene_name=gene_name;
				}
				line_info=br_info.readLine();
			}
			System.out.println("Finished extracting gene info.");
			writer.createGroup("/genotype");
			writer.writeIntArray("/genotype/gene_chrs", myFileFunctions.FileFunc.arraylist2intarray(gene_chrs_list));
			writer.writeIntArray("/genotype/gene_starts", myFileFunctions.FileFunc.arraylist2intarray(gene_starts_list));
			writer.writeIntArray("/genotype/gene_ends", myFileFunctions.FileFunc.arraylist2intarray(gene_ends_list));
			writer.writeStringArray("/genotype/gene_names", myFileFunctions.FileFunc.arraylist2stringarray(gene_name_list));
			writer.writeStringArray("/genotype/gene_ENS_ids", myFileFunctions.FileFunc.arraylist2stringarray(gene_id_list));
			System.out.println("Finished writing gene info to HDF5.");
			BufferedWriter bw_tmp=new BufferedWriter(new FileWriter("/Users/quan.long/tmp.gene.txt"));
			for(int k=0;k<num_of_genes;k++){
				bw_tmp.write(gene_chrs_list.get(k)+"\t"+gene_starts_list.get(k)+"\t"+gene_ends_list.get(k)+"\t"+gene_name_list.get(k)+"\t"+gene_id_list.get(k)+"\n");
			}bw_tmp.close();
			
			String[] temp=line.split(",");
			ArrayList<Integer> num_of_sites_per_chr=new ArrayList<Integer>();
			ArrayList<ArrayList<Integer>> var_locartions_list=new ArrayList<ArrayList<Integer>>();
			ArrayList<ArrayList<Integer>> var_mafc_list=new ArrayList<ArrayList<Integer>>();
			
			int sample_size=temp.length-2;
			String[] sample_ids=new String[sample_size];
			for(int k=0;k<sample_size;k++)sample_ids[k]=temp[k+2];
			line=br.readLine(); // the first data line
			temp=line.split(",");
			int old_chr=Integer.parseInt(temp[0]);
			int num_of_sites_this_chr=1;
			ArrayList<Integer> locations_this_chr=new ArrayList<Integer>();
			locations_this_chr.add(Integer.parseInt(temp[1]));
			ArrayList<Integer> mafc_this_chr=new ArrayList<Integer>();
			mafc_this_chr.add(mafc(temp));
			line=br.readLine();
			while(line!=null){
				temp=line.split(",");
				if(temp.length-2!=sample_size){System.out.println("Sample size not the same!");}
				int chr=Integer.parseInt(temp[0]);
				if(chr!=old_chr){
					var_locartions_list.add(locations_this_chr);
					var_mafc_list.add(mafc_this_chr);
					num_of_sites_per_chr.add(num_of_sites_this_chr);
					System.out.println("finished scaning Chr:"+old_chr);
					locations_this_chr=new ArrayList<Integer>();
					mafc_this_chr=new ArrayList<Integer>();
					num_of_sites_this_chr=0;
					old_chr=chr;
					
				}
				locations_this_chr.add(Integer.parseInt(temp[1]));
				mafc_this_chr.add(mafc(temp));
				num_of_sites_this_chr++;
				line=br.readLine();
			}// the last chr:
			System.out.println("finished scaning Chr:"+old_chr);
			var_locartions_list.add(locations_this_chr);
			var_mafc_list.add(mafc_this_chr);
			num_of_sites_per_chr.add(num_of_sites_this_chr);
			int num_chrs=num_of_sites_per_chr.size();
			int[] num_sites=new int[num_chrs];
			int[][] var_positions=new int[num_chrs][];
			int[][] var_mafc=new int[num_chrs][];
			for(int chr=0;chr<num_chrs;chr++){
				num_sites[chr]=num_of_sites_per_chr.get(chr);
				var_positions[chr]=new int[num_sites[chr]];
				var_mafc[chr]=new int[num_sites[chr]];
				ArrayList<Integer> locations_the_chr=var_locartions_list.get(chr);
				ArrayList<Integer> mafc_the_chr=var_mafc_list.get(chr);
				for(int k=0;k<num_sites[chr];k++){
					var_positions[chr][k]=locations_the_chr.get(k);
					var_mafc[chr][k]=mafc_the_chr.get(k);
				}
			}
			System.out.println("Finished scanning the .csv file and found "+num_chrs+" chromosomes.\nChr\t#variants");
			for(int chr=0;chr<num_chrs;chr++){
				System.out.println("Chr"+(chr+1)+":\t"+num_sites[chr]);
			}
			
			// CODE FOR constructing above structure (all groups and tables and importing data)
			writer.setIntAttribute("/genotype", "sample_size", sample_size);
			writer.setIntAttribute("/genotype", "block_size", block_size);
			writer.setIntAttribute("/genotype", "num_chrs", num_chrs);
			writer.writeStringArray("/genotype/sample_ids", sample_ids);
			int[] num_blocks_per_chr=new int[num_chrs];
//			HDF5CompoundType type = Annotation_g1k.class;
			for(int chr=0;chr<num_chrs;chr++){
				num_blocks_per_chr[chr]=num_sites[chr]/block_size+1;
				writer.createGroup("/genotype/chr"+(1+chr));
				writer.writeIntArray("/genotype/chr"+(1+chr)+"/var_pos", var_positions[chr]);
				writer.writeIntArray("/genotype/chr"+(1+chr)+"/var_mafc", var_mafc[chr]);				
				writer.createCompoundArray("/genotype/chr"+(1+chr)+"/annotation",  type, sample_size, block_size);
//				writer.writeCompoundArray(, arg1)
				writer.createByteMatrix("/genotype/chr"+(1+chr)+"/indi_fast",sample_size, num_sites[chr], sample_size, block_size);
				writer.createByteMatrix("/genotype/chr"+(1+chr)+"/position_fast", num_sites[chr], sample_size, block_size, sample_size);
			}writer.writeIntArray("/genotype/num_blocks_per_chr", num_blocks_per_chr);					
			// Read data again to fill matrix:
			System.out.println("Reading data again to create HDF5 file for genotypes.");
			br= new BufferedReader(new FileReader(file_csv));
			br_info= new BufferedReader(new FileReader(file_info));
			line=br.readLine();//header
			line_info=br_info.readLine();
			byte[][] data4thisblock=new byte[block_size][sample_size];
			Annotation_g1k[] info_data4thisblock=new Annotation_g1k[block_size];
			line=br.readLine();//first line
			line_info=br_info.readLine();
			temp=line.split(",");		
			int current_chr=Integer.parseInt(temp[0]);
			int total_var_in_chr=0;
			int current_block=0;
			int num_of_variants_in_block=0;
			while(line!=null){
				temp=line.split(",");
				int chr=Integer.parseInt(temp[0]);
				if(chr!=current_chr){
					//output the final block to h5
					writer.writeByteMatrixBlockWithOffset("/genotype/chr"+(current_chr)+"/position_fast",
							data4thisblock, num_of_variants_in_block, sample_size, current_block*block_size, 0L);
					writer.writeCompoundArrayBlockWithOffset("/genotype/chr"+(current_chr)+"/annotation", type, 
							create_array_with_less_elements(info_data4thisblock, num_of_variants_in_block), current_block*block_size);
					data4thisblock=new byte[block_size][sample_size];
					info_data4thisblock=new Annotation_g1k[block_size];
					num_of_variants_in_block=0;
					current_block=0;
					current_chr=chr;
					total_var_in_chr=0;
					System.out.println("finished writing Chr:"+current_chr);
				}else{//chr==current_chr
					if(num_of_variants_in_block==block_size){
						//output one block to h5
						writer.writeByteMatrixBlock("/genotype/chr"+(current_chr)+"/position_fast",
								data4thisblock, current_block,0);
						writer.writeCompoundArrayBlock("/genotype/chr"+(current_chr)+"/annotation", type, 
								info_data4thisblock, current_block);
						data4thisblock=new byte[block_size][sample_size];
						info_data4thisblock=new Annotation_g1k[block_size];
						num_of_variants_in_block=0;
						current_block++;
					}
				}
				for(int k=0;k<sample_size;k++){
					data4thisblock[num_of_variants_in_block][k]=Byte.parseByte(temp[k+2]);
				}
				info_data4thisblock[num_of_variants_in_block]=new Annotation_g1k(line_info, gene_name2index_map);
				num_of_variants_in_block++;
				total_var_in_chr++;
				line=br.readLine();
				line_info=br_info.readLine();
			}// write the final block:
			writer.writeByteMatrixBlockWithOffset("/genotype/chr"+(current_chr)+"/position_fast",
					data4thisblock, num_of_variants_in_block, sample_size, current_block*block_size, 0L);
			writer.writeCompoundArrayBlockWithOffset("/genotype/chr"+(current_chr)+"/annotation", type, 
					create_array_with_less_elements(info_data4thisblock, num_of_variants_in_block), current_block*block_size);
			System.out.println("finished writing Chr:"+current_chr);
			//write the indexes files		
			System.out.println("Finished writing genotypes. Now the HDF5 data is ready to use.");
			writer.close();			
			
		}catch(Exception e){e.printStackTrace();}		
		return new VariantsByte(file_hdf5);
	}
	
	/*
	 * maf count
	 */
	public static int mafc(String[] temp){
		HashMap<Integer, Integer> counts=new HashMap<Integer, Integer>();
		for(int k=2;k<temp.length;k++){ // from k=2, ignore the first two columns
			int the_allele=(int)(Double.parseDouble(temp[k])+0.5);
			if(counts.containsKey(the_allele)){
				int new_count=counts.get(the_allele)+1;
				counts.put(the_allele, new_count);
			}else{counts.put(the_allele, 1);}
		}
		int[] counts_array=new int[counts.size()];
		int[] allele_array=new int[counts.size()];
		int index=0;
		for(int allele:counts.keySet()){
			allele_array[index]=allele;
			counts_array[index++]=counts.get(allele);
		}if(counts_array.length==2){
			return ((counts_array[0]>counts_array[1])?counts_array[1]:counts_array[0])*2;
		}else if(counts_array.length==3){
			int zeor=0, other=0;
//			System.out.println("THREE");
			for(int n=0;n<3;n++){
				if(allele_array[n]==0){zeor+=2*counts_array[n];}
				else if(allele_array[n]==1){zeor+=counts_array[n]; other+=counts_array[n];}
				else if(allele_array[n]==2){other+=2*counts_array[n];}
				else{return -1;}
			}
			return (zeor<other)?zeor:other;
		}else{
			return -1;
		}
	}
	
	/*
	 * constructor: load an hdf5 object
	 */
	public VariantsByte(String file_hdf5){
		this.reader=HDF5Factory.openForReading(file_hdf5);	
		
		this.gene_names=this.reader.readStringArray("/genotype/gene_names");
		this.gene_ENS_ids=this.reader.readStringArray("/genotype/gene_ENS_ids");
		this.gene_name2index=myFileFunctions.FileFunc.generate_map_for_categories(gene_names);
		this.gene_chrs=this.reader.readIntArray("/genotype/gene_chrs");
		this.gene_starts=this.reader.readIntArray("/genotype/gene_starts");
		this.gene_ends=this.reader.readIntArray("/genotype/gene_ends");
		this.num_genes=this.gene_names.length;
		
		this.sample_size = this.reader.getIntAttribute("/genotype", "sample_size");
		this.block_size=this.reader.getIntAttribute("/genotype", "block_size");
		this.num_chrs = this.reader.getIntAttribute("/genotype","num_chrs");
		this.sample_ids=this.reader.readStringArray("/genotype/sample_ids");
		
		this.num_sites=new int[this.num_chrs];
		this.locations=new int[this.num_chrs][];
//		this.mafc=new int[this.num_chrs][];
		this.position_fast_blocks=new Iterable[this.num_chrs];
		this.position_fast_paths=new String[this.num_chrs];
	
		this.num_blocks_per_chr= this.reader.readIntArray("/genotype/num_blocks_per_chr");
		for(int chr=0;chr<num_chrs;chr++){
			this.locations[chr]=this.reader.readIntArray("/genotype/chr"+(1+chr)+"/var_pos");
//			this.mafc[chr]=this.reader.readIntArray("/genotype/chr"+(1+chr)+"/var_mafc");
			this.num_sites[chr]=this.locations[chr].length;
			this.position_fast_blocks[chr]=reader.getByteMDArrayNaturalBlocks("/genotype/chr"+(1+chr)+"/position_fast");		
			this.position_fast_paths[chr]="/genotype/chr"+(1+chr)+"/position_fast";
			this.num_sites_total+=this.num_sites[chr];
		}
		for(int chr=0;chr<this.num_chrs;chr++){			
//			for (HDF5MDDataBlock<MDDoubleArray> block : this.position_fast_blocks[chr]){
//				System.out.println((chr+1)+":"+ArrayUtils.toString(block.getIndex()) + " -> "+ block.getData().toString());
//				this.block_pointers[chr].add(block);				
//			}
		}		
		System.out.println("Finsehd constructing object based on "+file_hdf5);
	}
	
	/*
	 * Return the annotations in a region:
	 * chr starts from ZERO!
	 */
	public Annotation_g1k[] load_annotations_in_region(int chr, int start_location, int end_location){
		int max=this.locations[chr][this.locations[chr].length-1];
		if(start_location> max && end_location>max)	return null;
		int start_index=Arrays.binarySearch(this.locations[chr], start_location);
		if(start_index<0){// no_found, and make use of the insertion-point returned by the binarySearch:
			start_index=-(start_index+1);
		}
		int end_index=Arrays.binarySearch(this.locations[chr], end_location);
		if(end_index<0){// no_found, and make use of the insertion-point returned by the binarySearch:
			end_index=-(end_index+1)-1;
		}
		if(end_index<start_index)return null;
		HDF5CompoundType<Annotation_g1k> type = reader.compounds().getInferredType(Annotation_g1k.class);
		return this.reader.compounds().readArrayBlockWithOffset("/genotype/chr"+(chr+1)+"/annotation", type, 
				(end_index-start_index+1), start_index);	 
	}
    
	/*
	 * Return the variants in a region:
	 * chr starts from ZERO!
	 */
	public byte[][] load_variants_in_region_byte(int chr, int start_location, int end_location){
		int max=this.locations[chr][this.locations[chr].length-1];
		if(start_location> max && end_location>max)	return null;
		int start_index=Arrays.binarySearch(this.locations[chr], start_location);
		if(start_index<0){// no_found, and make use of the insertion-point returned by the binarySearch:
			start_index=-(start_index+1);
		}
		int end_index=Arrays.binarySearch(this.locations[chr], end_location);
		if(end_index<0){// no_found, and make use of the insertion-point returned by the binarySearch:
			end_index=-(end_index+1)-1;
		}
		if(end_index<start_index)return null;
		return reader.readByteMatrixBlockWithOffset(position_fast_paths[chr], 
				(end_index-start_index+1), this.sample_size, 
				start_index,0);		 
	}
	
	public double[][] load_variants_in_region(int chr, int start_location, int end_location){
		byte[][] variants_byte= load_variants_in_region_byte(chr, start_location, end_location);
		return myMathLib.Normalization.byte2double(variants_byte);
	}
	
	/*
	 * return an indicator array in which "true" stands for needed variants that passed the maf_threshold
	 */
	public static boolean[] find_variants_above_maf(byte[][] data, double allele_freq_threshold){
		int var_num=data.length;
		boolean[] indicators=new boolean[var_num];
		for(int k=0;k<var_num;k++){
			double count=0;
			for(int i=0;i<data[k].length;i++)
				count+=data[k][i];
			double maf=count/(data[k].length*2.0);
			//if(maf>0.5)maf=1-maf;
			if(maf>allele_freq_threshold)indicators[k]=true;
		}
		return indicators;
	}
	
	/*
	 * return an indicator array in which "true" stands for needed variants that passed the maf_threshold
	 */
	public static boolean[] find_variants_below_maf(byte[][] data, double allele_freq_threshold){
		int var_num=data.length;
		boolean[] indicators=new boolean[var_num];
		for(int k=0;k<var_num;k++){
			double count=0;
			for(int i=0;i<data[k].length;i++)
				count+=data[k][i];
			double maf=count/(data[k].length*2.0);
			//if(maf>0.5)maf=1-maf;
			if(maf<allele_freq_threshold)indicators[k]=true;
		}
		return indicators;
	}
	
	public double[] load_one_variant_by_location(int chr, int location){
		int index=Arrays.binarySearch(this.locations[chr], location);
		if(index<0){
			System.out.println("Error: Variant location not in the list");
			return null;
		}
		return load_one_variant_by_index(chr, index);
	}
	
	public double[] load_one_variant_by_index(int chr, int index){
		return myMathLib.Normalization.byte2double(reader.readByteMatrixBlockWithOffset(
				position_fast_paths[chr],1, this.sample_size, index, 0)[0]);
	}
	
	public byte[] load_one_variant_by_location_byte(int chr, int location){
		int index=Arrays.binarySearch(this.locations[chr], location);
		if(index<0){
			System.out.println("Error: Variant location not in the list");
			return null;
		}
		return load_one_variant_by_index_byte(chr, index);
	}
	
	public byte[] load_one_variant_by_index_byte(int chr, int index){
		return reader.readByteMatrixBlockWithOffset(position_fast_paths[chr], 
				1, this.sample_size, index, 0)[0];
	}

	
//	public double[][] readDoubleMatrixBlock_position_fast(int chr, long blockNumberX){
//		if(blockNumberX!=this.num_blocks_per_chr[chr]-1){
//			return reader.readDoubleMatrixBlock(position_fast_paths[chr], block_size, this.sample_size, 
//					blockNumberX, 0);
//		}else{
////			this.position_fast_blocks[chr].
////			return reader.readDoubleMatrixBlockWithOffset(position_fast_paths[chr], block_size, this.sample_size, 
////					block_size*(this.num_blocks_per_chr[chr]-1),0);
//			return reader.readDoubleMatrixBlockWithOffset(position_fast_paths[chr], 
//					(int)(this.num_sites[chr]-(this.num_blocks_per_chr[chr]-1)*block_size), this.sample_size, 
//					block_size*(this.num_blocks_per_chr[chr]-1),0);
//		}
//	}
	
	public void output2csv(String csv_file){
		try{
			BufferedWriter bw =new BufferedWriter(new FileWriter(csv_file));
			bw.write("CHROM,POS");
			for(int k=0;k<this.sample_size;k++)bw.write(","+this.sample_ids[k]);
			bw.write("\n");
			for(int chr=0;chr<this.num_chrs;chr++){
				int current_variant_index=0;
				for (HDF5MDDataBlock<MDByteArray> block : this.position_fast_blocks[chr]){
//					System.out.println((chr+1)+":"+ArrayUtils.toString(block.getIndex()) + " -> "+ block.getData().toString());
					byte[][] data4thisblock=block.getData().toMatrix();
					for(int i=0;i<data4thisblock.length;i++){
						bw.write((chr+1)+","+this.locations[chr][current_variant_index]);
						for(int k=0;k<this.sample_size;k++){
							bw.write(","+data4thisblock[i][k]);
						}bw.write("\n");
						current_variant_index++;
					}
				}
			}
			bw.close();
		}catch(Exception e){e.printStackTrace();}		
	}
	
	public void output2info(String info_file){
		try{
			BufferedWriter bw =new BufferedWriter(new FileWriter(info_file));
			bw.write("CHROM,POS");
//			for(int k=0;k<this.sample_size;k++)bw.write(","+this.sample_ids[k]);
			bw.write("\n");
			HDF5CompoundType<Annotation_g1k> type = reader.compounds().getInferredType(Annotation_g1k.class);
			for(int chr=0;chr<this.num_chrs;chr++){
				int current_variant_index=0;
				for (long block_num=0; block_num<this.num_blocks_per_chr[chr]; block_num++){
					Annotation_g1k[] block_anno=this.reader.compounds().readArrayBlock("/genotype/chr"+(chr+1)+"/annotation",
							type, this.block_size, block_num);
//					System.out.println((chr+1)+":"+ArrayUtils.toString(block.getIndex()) + " -> "+ block.getData().toString());
					for(int i=0;i<block_anno.length;i++){
						bw.write((chr+1)+"\t"+this.locations[chr][current_variant_index]+"\t"+block_anno[i].toString(this.gene_names));
						bw.write("\n");
						current_variant_index++;
					}
				}
			}
			bw.close();
		}catch(Exception e){e.printStackTrace();}		
	}
	/*
	 * only for testing whether load_annotations_in_region works well.
	 */
	public void output2info2(String info_file){
		try{
			BufferedWriter bw =new BufferedWriter(new FileWriter(info_file));
			bw.write("CHROM,POS");
//			for(int k=0;k<this.sample_size;k++)bw.write(","+this.sample_ids[k]);
			bw.write("\n");
			HDF5CompoundType<Annotation_g1k> type = reader.compounds().getInferredType(Annotation_g1k.class);
			for(int chr=0;chr<this.num_chrs;chr++){
				int current_variant_index=0;
				int step=50000, max=this.locations[chr][this.locations[chr].length-1];
				for (int loc=0; loc<=max+step; loc+=step){
					Annotation_g1k[] block_anno=this.load_annotations_in_region(chr, loc, loc+step-1);
//					System.out.println((chr+1)+":"+ArrayUtils.toString(block.getIndex()) + " -> "+ block.getData().toString());
					if(block_anno==null)continue;
					for(int i=0;i<block_anno.length;i++){
						bw.write((chr+1)+"\t"+this.locations[chr][current_variant_index]
								+"\t"+block_anno[i].toString(this.gene_names));
//						bw.write((chr+1)+"\t"+current_variant_index+"\t"+block_anno[i].toString(this.gene_names));
						bw.write("\n");
						current_variant_index++;
					}
				}
			}
			bw.close();
		}catch(Exception e){e.printStackTrace();}		
	}
	
	public double[][] read_in_kinship(String kinship_file){
		double[][] kinship=new double[this.sample_size][this.sample_size];
		String[] separators={","," ","\t"};
		String sep=null;
		try{			
			BufferedReader br= new BufferedReader(new FileReader(kinship_file));
			String line=br.readLine();//header
			for(int k=0;k<separators.length;k++){
				if(line.split(separators[k]).length==this.sample_size){
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
			if(line_index!=this.sample_size){
				System.out.println("Kinship file number of line wrong!");
			}
		}catch(Exception e){e.printStackTrace();}
		return kinship;
	}
	
	
	public void calculate_raw_kinship(String output_file, double scale){
		try{			
			double[][] kinship=new double[this.sample_size][this.sample_size];			
			for(int chr=0;chr<this.num_chrs;chr++){
				for (HDF5MDDataBlock<MDByteArray> block : this.position_fast_blocks[chr]){
					byte[][] data4thisblock=block.getData().toMatrix();
					for(int var_index=0;var_index<data4thisblock.length;var_index++){
						for(int i=0;i<sample_size;i++){
							for(int j=i+1;j<sample_size;j++){
								kinship[i][j]=kinship[i][j]+(scale-Math.abs(data4thisblock[var_index][i]-data4thisblock[var_index][j]));
							}
						}
					}
				}
			}
			for(int i=0;i<sample_size;i++){
				kinship[i][i]=1;
				for(int j=i+1;j<sample_size;j++){
					kinship[i][j]=kinship[i][j]/(num_sites_total*scale);
					kinship[j][i]=kinship[i][j];
				}
			}
			BufferedWriter bw= new BufferedWriter(new FileWriter(output_file));
			for(int i=0;i<sample_size;i++){
				for(int j=0;j<sample_size-1;j++){
					bw.write(kinship[i][j]+",");
				}bw.write(kinship[i][sample_size-1]+"\n");
			}
			bw.close();
			System.out.println("Global kinship (before rescale) has been written to "+output_file);
		}catch(Exception e){e.printStackTrace();}		
	}
	
	
	/*
	 * S_new=(n-1)S/Tr(PSP), where P= I-11'/n and 1 is the vector of ones.  
	 */
	public static double[][] re_scale_kinship_matrix(double[][] S){
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
		for(int i=0;i<n;i++){
			for(int j=0;j<n;j++){
				S_new[i][j]=S[i][j]*(n-1)/trace;
			}
		}
		return S_new;
	}
	
	public static void re_scale_kinship_matrix(String in_kinship_file, String out_kinship_file){
		String[] separators={","," ","\t"};
		try{
			String sep=null;
			int sample_size=-1; 
			if(new File(in_kinship_file).exists()){
				BufferedReader br = new BufferedReader(new FileReader(in_kinship_file));
				String line=br.readLine();
				for(int k=0;k<separators.length;k++){
					if(line.split(separators[k]).length>sample_size){
						sample_size=line.split(separators[k]).length;
						sep=separators[k];
					}
				}//System.out.println("Sep=\""+sep+"\"");
//				System.out.println("Sample size=\""+sample_size+"\"");
				double[][] old_kinship=new double[sample_size][sample_size];
				int line_index=0;
				while(line!=null){
					String[] temp=line.split(sep);
					for(int i=0;i<temp.length;i++){
						old_kinship[line_index][i]=Double.parseDouble(temp[i]);
					}
					line=br.readLine();
					line_index++;
				}
				double[][] new_kinship=re_scale_kinship_matrix(old_kinship);
				BufferedWriter bw= new BufferedWriter(new FileWriter(out_kinship_file));
				for(int i=0;i<sample_size;i++){
					for(int j=0;j<sample_size-1;j++){
						bw.write(new_kinship[i][j]+",");
					}bw.write(new_kinship[i][sample_size-1]+"\n");
				}
				bw.close();
			}
		}catch(Exception e){e.printStackTrace();}		
	}
	
	/*
	 * aggregate all SNPs into one "variant" using indicators (generated by functional, rare or LD)
	 */
	public byte[] aggregate_genotype(byte[][] genotype, boolean[] indicators){
		if(genotype==null)return null;
		if(indicators.length!=genotype.length){
			System.out.println("indicators.length!=genotype.length");
			return null;
		}
		byte[] result=new byte[this.sample_size];
		for(int k=0;k<genotype.length;k++){
			if(genotype[k].length!=this.sample_size){
				System.out.println("Sample size wrong at a site!");
				continue;
			}
			for(int i=0;i<this.sample_size;i++){
				if(genotype[k][i]!=0 && indicators[k])result[i]=1;
			}
		}
		return result;
	}
	
	public byte[] aggregate_genotype_functional_rare(int chr, int start_location, int end_location, double rare){
		byte[][] genotype=this.load_variants_in_region_byte(chr, start_location, end_location);
		Annotation_g1k[] annotations=this.load_annotations_in_region(chr, start_location, end_location);
		boolean[] indicators=functional_rare_indicator(genotype, annotations, rare);
		return aggregate_genotype(genotype, indicators);
	}
	
	public byte[] aggregate_genotype_rare(int chr, int start_location, int end_location, double rare){
		byte[][] genotype=this.load_variants_in_region_byte(chr, start_location, end_location);
		Annotation_g1k[] annotations=this.load_annotations_in_region(chr, start_location, end_location);
		boolean[] indicators=rare_indicator(genotype, annotations, rare);
		return aggregate_genotype(genotype, indicators);
	}
	
	public static boolean[] functional_rare_indicator(byte[][] genotype, Annotation_g1k[] annotations, double rare){
		if(annotations==null)return null;
		if(annotations.length!=genotype.length){
			System.out.println("WRONG: annotations.length!=genotype.length");
			return null;
		}
		boolean[] indicator=new boolean[annotations.length];
		for(int k=0;k<annotations.length;k++){
			if(annotations[k].AF<rare && annotations[k].functional())indicator[k]=true;
		}
		return indicator;
	}
	
	public static boolean[] rare_indicator(byte[][] genotype, Annotation_g1k[] annotations, double rare){
		if(annotations==null)return null;
		if(annotations.length!=genotype.length){
			System.out.println("WRONG: annotations.length!=genotype.length");
			return null;
		}
		boolean[] indicator=new boolean[annotations.length];
		for(int k=0;k<annotations.length;k++){
			if(annotations[k].AF<rare)indicator[k]=true;
		}
		return indicator;
	}
	
	public static void run_nam_imputation(String pedigreefiles_folder, String crossid, String RILpheno, 
			String geno250k, String output_prefix, int block_size){
		try{
			File[] pedfiles=(new File(pedigreefiles_folder)).listFiles();
			RILs[] ril_info = new RILs[pedfiles.length];
			Cross_info cross = new Cross_info(crossid);
			for(int i=0; i<pedfiles.length; i++){
				ril_info[i] = new RILs(pedfiles[i], RILpheno, cross);
			}
			Founder founders =new Founder(geno250k, cross);				
			Imputation.imputefull(founders,  ril_info, output_prefix+".pheno.tsv", output_prefix+".geno.csv");
			System.out.println("Imputation finished:\n"+
					"Genotype file: "+ output_prefix+".geno.csv\n"+"Phenotype file: "+output_prefix+".pheno.tsv");
			VariantsDouble.char2num(output_prefix+".geno.csv", output_prefix+".geno.num.csv");
			VariantsDouble.importCSV(output_prefix+".geno.num.csv", output_prefix+".geno.num.csv.hdf5", block_size);
			
		}catch(Exception e){e.printStackTrace();}		
	}
	
	public static void char2num(String input, String output){
		//TODO
		System.out.println("Convert genotype CSV file coded with ACGT/RWSYKM into CSV file coded with 0,1,2 (0=minor hom; 1=het; 2=major hom)\n" +
				"Lines containing more than 3 alleles will be omitted.");
		HashSet<String> homs=new HashSet<String>();
		HashSet<String> hets=new HashSet<String>();
		homs.add("A");homs.add("C");homs.add("G");homs.add("T");
		hets.add("R");hets.add("W");hets.add("S");hets.add("Y");hets.add("K");hets.add("M");
		try{
			BufferedReader br=new BufferedReader(new FileReader(input));
			String line=br.readLine();
			BufferedWriter bw=new BufferedWriter(new FileWriter(output));
			bw.write(line+"\n");//header
			int sample_size= line.split(",").length-2;
			System.out.println("Sample_size="+sample_size);
			line=br.readLine();
			int un=0, one=0, three=0, four=0;
			while(line!=null){
				HashMap<String, Integer> allele_counts =new HashMap<String, Integer>();
				String[] temp=line.split(",");
				if(temp.length-2!=sample_size){ // length wrong
					System.out.println("Line length != sample_size");
					line=br.readLine();
					continue;
				}
				boolean undefined_char=false;
				for(int i=2;i<temp.length;i++){					
					if(homs.contains(temp[i])){
						myFileFunctions.FileFunc.add2_counts_hashmap(allele_counts, temp[i], 2);
					}else if(hets.contains(temp[i])){
						String[] alleles=one2two(temp[i]).split(" ");
						for(int j=0;j<2;j++)myFileFunctions.FileFunc.add2_counts_hashmap(allele_counts, alleles[j], 1);
					}else{
						System.out.println("Char not ACGT/RWSYKM: "+temp[i]);
						undefined_char=true;
						break;
					}
				}if(undefined_char || allele_counts.size()!=2){ 
					line=br.readLine();
					if(undefined_char)un++;
					else if(allele_counts.size()==1)one++;
					else if(allele_counts.size()==3)three++;
					else if(allele_counts.size()==4)four++;
					continue;
				}
				// output this line
				String[] alleles=new String[2];
				int index=0;
				for(String the_allele:allele_counts.keySet()){
					alleles[index++]=the_allele;
				}String major=alleles[0], minor=alleles[1];
				if(allele_counts.get(major)<allele_counts.get(minor)){
					String buff=major; major=minor;minor=buff;
				}bw.write(temp[0]+","+temp[1]);
				for(int i=2;i<temp.length;i++){		
					if(temp[i].equals(major))bw.write(","+0);
					else if(temp[i].equals(minor))bw.write(","+2);
					else{
						String separated=one2two(temp[i]);
						if(separated.equals(major+" "+minor)||separated.equals(minor+" "+major)){
							bw.write(","+1);
						}else{ 
							System.out.println("Error: het not consistent to homs! "+temp[i]+":"+major+"/"+minor+"Exit.");
							return;
						}
					}
				}
				bw.write("\n");				
				line=br.readLine();
			}bw.close();
			System.out.println("Finished coverting: \n"+"Underfiend chars="+un+"\n"+
			"One Allele="+one+"\n"+"Three Allele="+three+"\n"+"Fou Alleler="+four+"\n");
		}catch(Exception e){e.printStackTrace();}
	}
	
	public static String one2two(String snp){
		String result="N N";
		if(snp.equals("A")||snp.equals("C")||snp.equals("G")||snp.equals("T"))
			result=snp+" "+snp;
		else{
	//		System.out.println("HET!");
			if(snp.equals("M")) result= "A C";
			else if(snp.equals("R")) result="A G";
			else if(snp.equals("W")) result="A T";
			else if(snp.equals("S")) result="C G";
			else if(snp.equals("Y")) result="C T";
			else if(snp.equals("K")) result="G T";
//			else System.out.println(snp+": Not a correct SNP!");	
		}return result;	
	}
	/*
	 * Return the variants in a region:
	 * chr starts from ZERO!
	 */
//	public double[][] load_variants_in_region_old(int chr, int start_location, int end_location){
//		int start_index=Arrays.binarySearch(this.locations[chr], start_location);
//		if(start_index<0){// no_found, and make use of the insertion-point returned by the binarySearch:
//			start_index=-(start_index+1);
//		}
//		int end_index=Arrays.binarySearch(this.locations[chr], end_location);
//		if(end_index<0){// no_found, and make use of the insertion-point returned by the binarySearch:
//			end_index=-(end_index+1)-1;
//		}
//		int start_block=this.location2block_index[chr][start_index];
//		int start_offset=this.location2offset_inblock[chr][start_index];
//		int end_block=this.location2block_index[chr][end_index];
//		int end_offset=this.location2offset_inblock[chr][end_index];
//		if(start_block>end_block){
//			System.out.println("WRONG: start_block>end_block.");
//			return null;
//		}		
//		int total_num_vars=(block_size-start_offset)+(end_block-start_block-1)*block_size+(1+end_offset);
//		double[][] the_data=new double[total_num_vars][];
//		if(start_block==end_block){
//			double[][] buffer=reader.readDoubleMatrixBlock(position_fast_paths[chr], block_size, this.sample_size, 
//					start_block, 0);
//			if(end_offset-start_offset+1!=total_num_vars)System.out.println("Wrong: single block but numbers not consistent.");
//			for(int i=start_offset;i<=end_offset;i++){
//				the_data[i-start_offset]=buffer[i];
//			}
//		}else{
//			int index_in_the_data=0;
//			double[][] buffer=this.readDoubleMatrixBlock_position_fast(chr, start_block);
//			for(int i=start_offset;i<block_size;i++){
//				the_data[index_in_the_data++]=buffer[i];
//			}
//			for(int block_index=start_block+1;block_index<=end_block-1;block_index++){
//				buffer=this.readDoubleMatrixBlock_position_fast(chr, block_index);
//				for(int i=0;i<block_size;i++){
//					the_data[index_in_the_data++]=buffer[i];
//				}
//			}			
//			buffer=this.readDoubleMatrixBlock_position_fast(chr, end_block);
//			for(int i=0;i<=end_offset;i++){
//				the_data[index_in_the_data++]=buffer[i];
//			}
//		}
//		return the_data;
//	}
	
	
	/*
	 * OLD FUNDTION THAT NO LONGER USEFUL
	 */
//	public static double[][] create_array_with_less_elements(double[][] data4thisblock, int num_of_variants_in_block){
//		double[][] sub_array=new double[num_of_variants_in_block][data4thisblock[0].length];
//		for(int k=0;k<num_of_variants_in_block;k++){
//			for(int i=0;i<data4thisblock[0].length;i++){
//				sub_array[k][i]=data4thisblock[k][i];
//			}
//		}
//		return sub_array;
//	}
//	
	public static Annotation_g1k[] create_array_with_less_elements(Annotation_g1k[] data4thisblock, int num_of_variants_in_block){
		Annotation_g1k[] sub_array= new Annotation_g1k[num_of_variants_in_block];
		for(int k=0;k<num_of_variants_in_block;k++){
			sub_array[k]=data4thisblock[k];
		}
		return sub_array;
	}
}


