package myFileFunctions;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

public class FileFunc {
	public static void add2_counts_hashmap(HashMap<String, Integer> counts, String item, int value){
		if(counts.containsKey(item)){
			int new_count=counts.get(item)+value;
			counts.put(item, new_count);
		}else{
			counts.put(item, value);
		}
	}
	
	public static void add2_counts_hashmap(HashMap<Integer, Integer> counts, Integer item, int value){
		if(counts.containsKey(item)){
			int new_count=counts.get(item)+value;
			counts.put(item, new_count);
		}else{
			counts.put(item, value);
		}
	}
	
	public static HashMap<String, Integer> generate_chr_map(){
		HashMap<String, Integer> chr2index= new HashMap<>();
		for(int k=1;k<=25;k++){
			chr2index.put("chr"+k, k);
			chr2index.put("Chr"+k, k);
			chr2index.put(""+k, k);
		}
		chr2index.put("X", 23);chr2index.put("chrX", 23);chr2index.put("ChrX", 23);
		chr2index.put("Y", 24);chr2index.put("chrY", 24);chr2index.put("ChrY", 24);
		chr2index.put("MT", 25);chr2index.put("chrMT", 25);chr2index.put("ChrMT", 25);
		return chr2index;
	}
	public static double digital(double data, int num_dig){
		double t=1;
		for(int k=0;k<num_dig;k++)t=t*10;
		return ((int)(data*t))/t;
	}
	
	public static void add2_sum_hashmap(HashMap<String, Double> counts, String item, double value){
		if(counts.containsKey(item)){
			double new_count=counts.get(item)+value;
			counts.put(item, new_count);
		}else{
			counts.put(item, value);
		}
	}
	
	public static String[] mysplit(String line, char sep){
		char[] array=line.toCharArray();
		int number=0;
		for(int i=0;i<array.length;i++){
			if(array[i]==sep)
				number++;
		}
		String[] result=new String[number+1];
		int[] lens=new int[number+1];
		int index=0, len=0;
		for(int i=0;i<array.length;i++){
			if(array[i]!=sep)
				len++;
			else{
				lens[index]=len;
				len=0;index++;
			}				
		}lens[number]=len;
		int current_token=0, current_loc=0;
		for(int i=0;i<number+1;i++){
			char[] token=new char[lens[i]];
			for(int j=0;j<lens[i];j++){
				token[j]=array[current_loc+j];
			}
			current_loc=current_loc+lens[i]+1;
			result[i]=new String(token);
		}
		return result;
	}
	
	public static HashMap<Character, Integer> allele2index(){
		HashMap<Character, Integer> allele2index= new HashMap<>();
		allele2index.put('A', 0);allele2index.put('a', 0);
		allele2index.put('C', 1);allele2index.put('c', 1);
		allele2index.put('G', 2);allele2index.put('g', 2);
		allele2index.put('T', 3);allele2index.put('t', 3);
		return allele2index;
	}
	
	/*
	 * sort by counts[] and other[] also changed accordingly. 
	 */
	public static void sort_with(int[] counts, char[] other){
		for(int i=0;i<counts.length;i++){
			for(int j=i+1;j<counts.length;j++){
				if(counts[i]>counts[j]){
					int temp=counts[i];
					counts[i]=counts[j];
					counts[j]=temp;
					char xx = other[i];
					other[i]=other[j];
					other[j]=xx;
				}
			}
		}
	}
	
	public static ArrayList<Integer> sorted_keys(HashMap<Integer, String> data){
		ArrayList<Integer> keys= new ArrayList<>();
		for(int x: data.keySet()){
			keys.add(x);
		}
		Collections.sort(keys);
		return keys;
	}
	
	public static int[] arraylist2arrayInteger(ArrayList<Integer> data){
		int[] out=new int[data.size()];
		for(int k=0;k<out.length;k++)out[k]=data.get(k);
		return out;
	}
	
	public static double[] arraylist2arrayDouble(ArrayList<Double> data){
		double[] out=new double[data.size()];
		for(int k=0;k<out.length;k++)out[k]=data.get(k);
		return out;
	}
	
	public static String[] arraylist2arrayString(ArrayList<String> data){
		String[] out=new String[data.size()];
		for(int k=0;k<out.length;k++)out[k]=data.get(k);
		return out;
	}
	
	public static HashMap<String, Integer> generate_map_for_categories(String[] categories){
		HashMap<String, Integer> map= new HashMap<>();
		for(int i=0;i<categories.length;i++){
			map.put(categories[i], i);
		}return map;
	}
	
	public static double[] hashmap2array_sort(HashMap<Double, Integer> map){
		if(map==null)return null;
		if(map.size()==0)return new double[0];
		double[] array = new double[map.size()];
		int index=0;
		for(Iterator it=map.keySet().iterator();it.hasNext();){
			array[index++]=(Double)(it.next());
		}Arrays.sort(array);
		return array;
	}
	
	public static double[] hashset2array_sort(HashSet<Double> set){
		if(set==null)return null;
		if(set.size()==0)return new double[0];
		double[] array = new double[set.size()];
		int index=0;
		for(Iterator it=set.iterator();it.hasNext();){
			array[index++]=(Double)(it.next());
		}Arrays.sort(array);
		return array;
	}
	
	public static HashMap<String, String> generate_map_from_file(String file, int col_key, int col_value, String sep){
		HashMap<String, String> map = new HashMap<>();
		try{
			BufferedReader br=new BufferedReader(new FileReader(file));
			String line=br.readLine();
			while(line!=null){
				String[] tmp=line.split(sep);
				if(map.containsKey(tmp[col_key])){
					System.out.println("Warning: "+col_key+" alreday in the keyset of the hashmap.");
				}
				map.put(tmp[col_key], tmp[col_value]);
				line=br.readLine();
			}
		}catch(Exception e){e.printStackTrace();}
		return map;
	}
	
	public static void gunzip(String inFilename,  String outFilename){
		try {
		    // Open the compressed file
		    GZIPInputStream in = new GZIPInputStream(new FileInputStream(inFilename));

		    // Open the output file
		    OutputStream out = new FileOutputStream(outFilename);

		    // Transfer bytes from the compressed file to the output file
		    byte[] buf = new byte[1024];
		    int len;
		    while ((len = in.read(buf)) > 0) {
		        out.write(buf, 0, len);
		    }

		    // Close the file and stream
		    in.close();
		    out.close();
		} catch (IOException e) {
			//TODO: Add catch message
		}
	}
	
	public static void gzip( String inFilename, String outFilename ){
		try {
		    // Create the GZIP output stream
		    
		    GZIPOutputStream out = new GZIPOutputStream(new FileOutputStream(outFilename));

		    // Open the input file
		   
		    FileInputStream in = new FileInputStream(inFilename);

		    // Transfer bytes from the input file to the GZIP output stream
		    byte[] buf = new byte[1024];
		    int len;
		    while ((len = in.read(buf)) > 0) {
		        out.write(buf, 0, len);
		    }
		    in.close();

		    // Complete the GZIP file
		    out.finish();
		    out.close();
		} catch (IOException e) {
			//TODO: Add catch message
		}
	}
	
	public static int chr2int(String chr){
		if(chr.equals("X")){
			return 23;
		}if(chr.equals("Y")){
			return 24;
		}if(chr.equals("XY")){
			return 25;
		}if(chr.equals("MT")){
			return 26;
		}else{
			return Integer.parseInt(chr);
		}
	}
	
	public static String[] arraylist2stringarray(ArrayList<String> input){
		String[] output=new String[input.size()];
		for(int k=0;k<output.length;k++)output[k]=input.get(k);
		return output;
	}
	
	public static int[] arraylist2intarray(ArrayList<Integer> input){
		int[] output=new int[input.size()];
		for(int k=0;k<output.length;k++)output[k]=input.get(k);
		return output;
	}
	
	public static int[][] arraylist2intarray2(ArrayList<int[]> input){
		int[][] output=new int[input.size()][];
		for(int k=0;k<output.length;k++)output[k]=input.get(k).clone();
		return output;
	}
	
	/*sort file according to the column specified (starting from 0) and separator, and the number of header ignored
	 *  the result will be written into the same file
	 */
	public static void sort_file(String file, int column, int num_of_headers, String sep){
		//TODO
	}
}
