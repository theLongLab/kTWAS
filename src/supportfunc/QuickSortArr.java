package supportfunc;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;

/*
 * to sort chr_los
 */

public class QuickSortArr {	
	
	public static String[] content2Arr(String type, String input_file){
		ArrayList<String> content_list = new ArrayList<String>();
		try {
			BufferedReader inputbr = new BufferedReader(new FileReader(input_file));
			String input_line = inputbr.readLine();
			
			/*
			 * different type different read method
			 */
			if(type.equals("num") || type.equals("SNP")) { //1\n5\n3 OR 1,100\n1,400
				while(input_line != null) {
					if(!input_line.startsWith("#")) {
						String temp = input_line.trim();
						content_list.add(temp);
					}
					input_file=inputbr.readLine();
				}
			}else if(type.equals("compound")) {//c1\t1;s1;e1\t1;s2;e2
				while(input_line != null) {
					if(!input_line.startsWith("#")) {
				        String[] input_line_arr = input_line.split("\t");
				        String compound = String.valueOf(input_line_arr[1]) + "_" + input_line_arr[2]; // 1;s1;e1_1;s2;e2
				        if (!content_list.contains(compound)) {
				        	content_list.add(compound);
				        }
				        input_line = inputbr.readLine();
					}
				}
			}else {
				System.out.println("not valided type, please try again");
				System.exit(0);
			}
			inputbr.close();
				
		}catch(Exception e) {
			e.printStackTrace();
		}
		String[] contentArr =  new String[content_list.size()];
		contentArr = content_list.toArray(contentArr);
		return contentArr;
	}
	
	
	public static void output_sorted (String[] contentArr, String output_file, String type) {
		try {
			BufferedWriter outputbw = new BufferedWriter(new FileWriter(output_file));
			quicksort(contentArr, 0, contentArr.length-1, type);		
			if(type.equals("num") || type.equals("SNP")) {
				for (int i = 0; i < contentArr.length; i++) {
			        outputbw.write(contentArr[i] + "\n");
			      }
			}else if(type.equals("compound")){
				for (int i = 0; i < contentArr.length; i++) {
					String[] temp_arr = contentArr[i].split("_");
			        outputbw.write("C" + i + "\t" + temp_arr[0] + "\t" + temp_arr[1] + "\n");
			      } 
			}
			outputbw.close();
		}catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	public static void quicksort(String[] contentArr, int low, int high, String type){
		if(low < high) {
			int pi = patrition(contentArr, low, high, type);
			quicksort(contentArr, low, pi-1, type);
			quicksort(contentArr, pi+1, high, type);
		}
	}
	
	public static int patrition(String[] contentArr, int low, int high, String type) {
		/*
		 * set the high to be pivot
		 */
		int i = low -1;
		String pivot = contentArr[high];
		for(int j=low; j< high; j++) {
			
			if(type.equals("num")) {
				if(checkincreaseorderNum(contentArr[j],pivot)) {
					i+=1;
					swap(contentArr, i, j);
				}
			}else if(type.equals("SNP")) {
				if(checkincreaseorderSNP(contentArr[j],pivot)) {
					i+=1;
					swap(contentArr, i, j);
				}
			}else if(type.equals("compound")) {
				if(checkincreaseorderCompound(contentArr[j],pivot)) {
					i+=1;
					swap(contentArr, i, j);
				}
			}
		}
		swap(contentArr, i+1, high);
		return i+1;
	}
	
	public static void swap(String[] snps_arr, int i, int j) {
		String temp = snps_arr[i];
		snps_arr[i]=snps_arr[j];
		snps_arr[j]=temp;
	}
	
	public static boolean checkincreaseorderNum(String x, String y) {
		boolean increaseorder = false;
		if(Integer.parseInt(x) <= Integer.parseInt(y)) {
			increaseorder = true;
		}
		return increaseorder;
	}
	
	public static boolean checkincreaseorderSNP(String x, String y) {
		boolean increaseorder = false;
		String[] x_arr = x.split(",");
		String[] y_arr = y.split(",");
		if(Integer.parseInt(x_arr[0])<Integer.parseInt(y_arr[0])) {
			increaseorder=true; // x_CHR is smaler than y_CHR
		}else if (Integer.parseInt(x_arr[0])==Integer.parseInt(y_arr[0])) {
			if(Integer.parseInt(x_arr[1])<=Integer.parseInt(y_arr[1])) {// x_CHR = y_CHR, x_POS¡¡£¼=¡¡£ù_POS
				increaseorder=true;
			}
		}
		return increaseorder;
	}
	
	public static boolean checkincreaseorderCompound(String x, String y) {
		// TODO Auto-generated method stub
	    boolean increaseorder = false;
	    String[] x_arr = x.split("_");
	    String[] y_arr = y.split("_");
	    String[] x_arr0_arr = x_arr[0].split(";");
	    String[] y_arr0_arr = y_arr[0].split(";");
	    if (Integer.parseInt(x_arr0_arr[0]) < Integer.parseInt(y_arr0_arr[0])) {
	      increaseorder = true;
	    } else if (Integer.parseInt(x_arr0_arr[0]) == Integer.parseInt(y_arr0_arr[0]) && 
	      Integer.parseInt(x_arr0_arr[1]) <= Integer.parseInt(y_arr0_arr[1])) {
	      increaseorder = true;
	    } 
	    return increaseorder;
	}
	
	public static void main(String args[]) {
		if (args.length == 0) {
			System.out.println("inputs: \n"
					+ "1) original num/SNP/compound Info file \n"
					+ "2) type of data: num/SNP/compound \n"
					+ "outputs : \n"
					+ "1) sorted and unique HiC Info file \n" + 
		"enter all above files path");
		} else {
		    String input_file = args[0];
		    String type = args[1];
		    String output_file = String.valueOf(args[2]);
//			String input_file="C:/Users/Administrator/Desktop/26499245_master_pns.compound.txt";
//			String type="compound"; 
//			String output_file="C:/Users/Administrator/Desktop/26499245_master_pns.compound.sorted.txt";
		    String[] contentArr = content2Arr(type, input_file);
		    output_sorted(contentArr, output_file, type);
		}
	}
}
