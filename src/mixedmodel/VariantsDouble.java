package mixedmodel;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import ch.systemsx.cisd.hdf5.*;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.correlation.Covariance;

import nam.Cross_info;
import nam.Founder;
import nam.Imputation;
import nam.RILs;


//import cern.colt.matrix.impl.DenseDoubleMatrix2D;

import ch.systemsx.cisd.base.mdarray.MDDoubleArray;

public class VariantsDouble {

    public int sample_size;
    public int num_chrs;
    public int[] num_sites;
    public int[][] locations;  // chr x num_sites[chr]
    public int[][] mafc;        // chr x num_sites[chr]
    public String[] sample_ids;
    public HashMap<String, Integer> sample_id2index;
    boolean store_in_byte = false;
    long num_sites_total;
    int block_size = 5000;
    int[] num_blocks_per_chr;
    IHDF5Reader reader;
    Iterable<HDF5MDDataBlock<MDDoubleArray>>[] position_fast_blocks;
    //	ArrayList<HDF5MDDataBlock<MDDoubleArray>>[] block_pointers;
//	int[][] location2block_index;
//	int[][] location2offset_inblock;
    String[] position_fast_paths;
    public HashMap<String, double []> load_snps_by_index ;

    boolean transformed = false;
    double[] intercepts; // intercepts has to be an all-one matrix if not transformed.
    

    /*
     * 	ROOT:
     *  	Group:  Genotype:
     *  			Attr: some info: like sample-size, number-of-variants, ... ,
     *  			Table: Individual ids: 1D array of Strings
     *  			N Groups: for chromosomes
     *  				Each-groups: Two double arrays: individual-fast; location-fast
     *  							 Table: Variants positions: 1D array of ints
     *  							 Table: Annotations of the variants: 1D array of Compound-Object
     *
     */


    /*
     * constructor: load an hdf5 object
     */
    // CHECKED - PATHUM
    
    
    public VariantsDouble(String file_csv, String format) throws IOException {
    	if (format.equals("csv")){
    		this.load_snps_by_index = new HashMap<String, double []>();
    		BufferedReader br = new BufferedReader(new FileReader(file_csv));
            String line = br.readLine();//header
            String[] temp = line.split(",");
            this.sample_size = temp.length - 2;
            this.sample_ids = new String[sample_size];
            for (int k = 0; k < sample_size; k++) this.sample_ids[k] = temp[k + 2];
      
            int chr_index =0;
            int pos_index=0;
            String chr_flag = "";
            line = br.readLine(); // the first data line
            while (line != null) {
                temp = line.split(",");
                if (temp.length - 2 != sample_size) {
                    System.out.println("Sample size not the same!");
                }
//                System.out.println(temp[1]);
                if  ((temp[0].equals(chr_flag)) && (chr_index!=0)) {
                	chr_index++;
                	chr_flag= temp[0];
                	pos_index=0;
                } else {
                	double [] tmp_genotype = new double [sample_size];
                	for (int i=2;i< temp.length;i++) {
                		tmp_genotype[i-2]=Double. parseDouble(temp[i]);
                	}
                	String key = Integer.toString(chr_index)+":" + Integer.toString(pos_index);
//                	System.out.println(tmp_genotype);
                	this.load_snps_by_index.put(key, tmp_genotype);
                	pos_index++;
                	
                }
                line = br.readLine();
            }// the last chr:
            br.close();
            System.out.println("Finished scanning the .csv file and found " + ( chr_index+1) + " chromosomes.");
    	} else {
    		System.out.println("The program do not recognize the format:\t"+format);
    	}
    }
    
    public VariantsDouble(String file_hdf5) {
        this.reader = HDF5Factory.openForReading(file_hdf5);
        this.sample_size = this.reader.int32().getAttr("/genotype", "sample_size");
        this.block_size = this.reader.int32().getAttr("/genotype", "block_size");
        this.num_chrs = this.reader.int32().getAttr("/genotype", "num_chrs");
        this.sample_ids = this.reader.readStringArray("/genotype/sample_ids");
//		this.transformed=this.reader.getBooleanAttribute("/genotype","transformed");
//		this.intercepts=this.reader.readDoubleArray("/genotype/intercepts");

        this.num_sites = new int[this.num_chrs];
        this.locations = new int[this.num_chrs][];
        this.mafc = new int[this.num_chrs][];
        this.position_fast_blocks = new Iterable[this.num_chrs];
        this.position_fast_paths = new String[this.num_chrs];

        this.num_blocks_per_chr = this.reader.readIntArray("/genotype/num_blocks_per_chr");

        for (int chr = 0; chr < num_chrs; chr++) {
            this.locations[chr] = this.reader.readIntArray("/genotype/chr" + (1 + chr) + "/var_pos");
            this.mafc[chr] = this.reader.readIntArray("/genotype/chr" + (1 + chr) + "/var_mafc");
            this.num_sites[chr] = this.locations[chr].length;
            this.position_fast_blocks[chr] = reader.float64().getMDArrayNaturalBlocks("/genotype/chr" + (1 + chr) + "/position_fast");
            this.position_fast_paths[chr] = "/genotype/chr" + (1 + chr) + "/position_fast";
            this.num_sites_total += this.num_sites[chr];
        }

        this.sample_id2index = new HashMap<String, Integer>();

        for (int k = 0; k < this.sample_size; k++) this.sample_id2index.put(this.sample_ids[k], k);

        System.out.println("Finished reading " + file_hdf5 + ".");
    }

    /*
     * constructor: make an hdf5 object based on a .csv file
     */
    public static VariantsDouble importCSV(String file_csv, String file_hdf5, int block_size) {
        System.out.println("Start importing data from " + file_csv + " to HDF5 format.");
        try {
            BufferedReader br = new BufferedReader(new FileReader(file_csv));
            String line = br.readLine();//header
            String[] temp = line.split(",");
            ArrayList<Integer> num_of_sites_per_chr = new ArrayList<>();
            ArrayList<ArrayList<Integer>> var_locartions_list = new ArrayList<>();
            ArrayList<ArrayList<Integer>> var_mafc_list = new ArrayList<>();
            int sample_size = temp.length - 2;
            String[] sample_ids = new String[sample_size];
            for (int k = 0; k < sample_size; k++) sample_ids[k] = temp[k + 2];
            line = br.readLine(); // the first data line
            temp = line.split(",");
            int old_chr = Integer.parseInt(temp[0]);
            int num_of_sites_this_chr = 1;
            ArrayList<Integer> locations_this_chr = new ArrayList<>();
            locations_this_chr.add(Integer.parseInt(temp[1]));
            ArrayList<Integer> mafc_this_chr = new ArrayList<>();
            mafc_this_chr.add(mafc(temp));
            line = br.readLine();
            while (line != null) {
                temp = line.split(",");
                if (temp.length - 2 != sample_size) {
                    System.out.println("Sample size not the same!");
                }
                int chr = Integer.parseInt(temp[0]);
                if (chr != old_chr) {
                    var_locartions_list.add(locations_this_chr);
                    var_mafc_list.add(mafc_this_chr);
                    num_of_sites_per_chr.add(num_of_sites_this_chr);
                    locations_this_chr = new ArrayList<Integer>();
                    mafc_this_chr = new ArrayList<Integer>();
                    num_of_sites_this_chr = 0;
                    old_chr = chr;
                }
                locations_this_chr.add(Integer.parseInt(temp[1]));
                mafc_this_chr.add(mafc(temp));
                num_of_sites_this_chr++;
                line = br.readLine();
            }// the last chr:
            var_locartions_list.add(locations_this_chr);
            var_mafc_list.add(mafc_this_chr);
            num_of_sites_per_chr.add(num_of_sites_this_chr);
            int num_chrs = num_of_sites_per_chr.size();
            int[] num_sites = new int[num_chrs];
            int[][] var_positions = new int[num_chrs][];
            int[][] var_mafc = new int[num_chrs][];
            for (int chr = 0; chr < num_chrs; chr++) {
                num_sites[chr] = num_of_sites_per_chr.get(chr);
                var_positions[chr] = new int[num_sites[chr]];
                var_mafc[chr] = new int[num_sites[chr]];
                ArrayList<Integer> locations_the_chr = var_locartions_list.get(chr);
                ArrayList<Integer> mafc_the_chr = var_mafc_list.get(chr);
                for (int k = 0; k < num_sites[chr]; k++) {
                    var_positions[chr][k] = locations_the_chr.get(k);
                    var_mafc[chr][k] = mafc_the_chr.get(k);
                }
            }
            System.out.println("Finished scanning the .csv file and found " + num_chrs + " chromosomes.\nChr\t#variants");
            for (int chr = 0; chr < num_chrs; chr++) {
                System.out.println("Chr" + (chr + 1) + ":\t" + num_sites[chr]);
            }
            File hdf5_old = new File(file_hdf5);
            if (hdf5_old.exists()) {
                hdf5_old.delete();
                System.out.println(file_hdf5 + " already existed. It has been deleted!");
            }
            IHDF5Writer writer = HDF5Factory.open(file_hdf5);
            // CODE FOR constructing above structure (all groups and tables and importing data)
            writer.object().createGroup("/genotype");
            writer.int32().setAttr("/genotype", "sample_size", sample_size);
            writer.int32().setAttr("/genotype", "block_size", block_size);
            writer.int32().setAttr("/genotype", "num_chrs", num_chrs);
            writer.bool().setAttr("/genotype", "transformed", false);
            writer.writeStringArray("/genotype/sample_ids", sample_ids);
            double[] all_one = new double[sample_size];
            Arrays.fill(all_one, 1.0);
            writer.writeDoubleArray("/genotype/intercepts", all_one);
            int[] num_blocks_per_chr = new int[num_chrs];
            for (int chr = 0; chr < num_chrs; chr++) {
                num_blocks_per_chr[chr] = num_sites[chr] / block_size + 1;
                writer.object().createGroup("/genotype/chr" + (1 + chr));
                writer.writeIntArray("/genotype/chr" + (1 + chr) + "/var_pos", var_positions[chr]);
                writer.writeIntArray("/genotype/chr" + (1 + chr) + "/var_mafc", var_mafc[chr]);
                writer.float64().createMatrix("/genotype/chr" + (1 + chr) + "/indi_fast", sample_size, num_sites[chr], sample_size, block_size);
                writer.float64().createMatrix("/genotype/chr" + (1 + chr) + "/position_fast", num_sites[chr], sample_size, block_size, sample_size);
            }
            writer.writeIntArray("/genotype/num_blocks_per_chr", num_blocks_per_chr);
            // Read data again to fill matrix:
            System.out.println("Reading data again to create HDF5 file for genotypes.");
            br = new BufferedReader(new FileReader(file_csv));
            line = br.readLine();//header
            double[][] data4thisblock = new double[block_size][sample_size];
            line = br.readLine();//first line
            temp = line.split(",");
            int current_chr = Integer.parseInt(temp[0]);
            int total_var_in_chr = 0; // TODO: Find purpose of this variable
            int current_block = 0;
            int num_of_variants_in_block = 0;
            while (line != null) {
                temp = line.split(",");
                int chr = Integer.parseInt(temp[0]);
                if (chr != current_chr) {
                    //output the final block to h5
                    writer.float64().writeMatrixBlockWithOffset("/genotype/chr" + (current_chr) + "/position_fast",
                            data4thisblock, num_of_variants_in_block, sample_size, current_block * block_size, 0L);
                    data4thisblock = new double[block_size][sample_size];
                    num_of_variants_in_block = 0;
                    current_block = 0;
                    current_chr = chr;
                    total_var_in_chr = 0;
                } else {//chr==current_chr
                    if (num_of_variants_in_block == block_size) {
                        //output one block to h5
                        writer.float64().writeMatrixBlock("/genotype/chr" + (current_chr) + "/position_fast",
                                data4thisblock, current_block, 0);
                        data4thisblock = new double[block_size][sample_size];
                        num_of_variants_in_block = 0;
                        current_block++;
                    }
                }
                for (int k = 0; k < sample_size; k++) {
                    data4thisblock[num_of_variants_in_block][k] = Double.parseDouble(temp[k + 2]);
                }
                num_of_variants_in_block++;
                total_var_in_chr++;
                line = br.readLine();
            }// write the final block:
            writer.float64().writeMatrixBlockWithOffset("/genotype/chr" + (current_chr) + "/position_fast",
                    data4thisblock, num_of_variants_in_block, sample_size, current_block * block_size, 0L);
            //write the indexes files
            System.out.println("Finished writing genotypes. Now the HDF5 data is ready to use.");
            writer.close();

        } catch (Exception e) {
            e.printStackTrace();
        }
        return new VariantsDouble(file_hdf5);
    }

//	/*
//	 * maf count
//	 */
//	public static int mafc(String[] temp){
//		HashMap<Integer, Integer> counts=new HashMap<Integer, Integer>();
//		for(int k=2;k<temp.length;k++){ // from k=2, ignore the first two columns
//			int the_allele=(int)(Double.parseDouble(temp[k])+0.5);
//			if(counts.containsKey(the_allele)){
//				int new_count=counts.get(the_allele)+1;
//				counts.put(the_allele, new_count);
//			}else{counts.put(the_allele, 1);}
//		}
//		int[] counts_array=new int[counts.size()];
//		int[] allele_array=new int[counts.size()];
//		int index=0;
//		for(int allele:counts.keySet()){
//			allele_array[index]=allele;
//			counts_array[index++]=counts.get(allele);
//		}if(counts_array.length==2){  // hom only
//			return ((counts_array[0]>counts_array[1])?counts_array[1]:counts_array[0])*2;
//		}else if(counts_array.length==3){ // also hets
//			int zeor=0, other=0;
////			System.out.println("THREE");
//			for(int n=0;n<3;n++){
//				if(allele_array[n]==0){zeor+=2*counts_array[n];}
//				else if(allele_array[n]==1){zeor+=counts_array[n]; other+=counts_array[n];}
//				else if(allele_array[n]==2){other+=2*counts_array[n];}
//				else{return -1;}
//			}
//			return (zeor<other)?zeor:other;
//		}else{
//			return -1;
//		}
//	}

    public static int mafc(String[] temp) {
        int zero_count = 0, other_count = 0;
        for (int k = 2; k < temp.length; k++) { // from k=2, ignore the first two columns
            int the_allele = (int) (Double.parseDouble(temp[k]) + 0.1); // convert into integer
            if (the_allele == 0) {
                zero_count += 2;
            } else if (the_allele == 1) {
                zero_count += 1;
                other_count += 1;
            } else if (the_allele == 2) {
                other_count += 2;
            }
        }
        return (zero_count < other_count) ? zero_count : other_count;
    }

    /*
     * Following Yang et al AJHG 2011, calculating the inner product of the normalized genotype.
     * Genotype matrix will be normalized as w(i,j)=(x(i,j)-E(x(i)))/Sd(x(i)). Using the kinship
     * formula like this, the rescaling using Gower's centered matrix is not necessary
     *
     * The input data, as default in hdf5 storage will be num_var x sample_size;
     */
    public static double[][] calculate_RRM_local_kinship(double[][] data) {
        if (data == null) {
            return null;
        }
        int sample_size = data[0].length;
        int num_var = data.length;
        double[] means = new double[num_var];
        double[] sd = new double[num_var];
        for (int p = 0; p < num_var; p++) {
            means[p] = StatUtils.mean(data[p]);
            sd[p] = Math.sqrt(StatUtils.variance(data[p], means[p]));
        }
        double[][] data_W = new double[sample_size][num_var];
        for (int k = 0; k < sample_size; k++) {
            for (int p = 0; p < num_var; p++) {
                data_W[k][p] = (data[p][k] - means[p]) / sd[p];
            }
        }
        Array2DRowRealMatrix data_W_matrix = new Array2DRowRealMatrix(data_W);
        Array2DRowRealMatrix data_W_matrixT = (Array2DRowRealMatrix) (data_W_matrix.transpose());
        double[][] kinship = ((data_W_matrix.multiply(data_W_matrixT)).scalarMultiply(1.0 / (double) num_var)).getData();//  new double[sample_size][sample_size];
//		System.out.println("Local kinship (before rescale) has been calculted.");
        return kinship;
    }

    /*
     * return X'X after doing normalization.
     */
    public static RealMatrix calculate_RRM_local_kinship_RealMatrix(double[][] data) {
        if (data == null) {
            return null;
        }
        int sample_size = data[0].length;
        int num_var = data.length;
        double[] means = new double[num_var];
        double[] sd = new double[num_var];
        for (int p = 0; p < num_var; p++) {
            means[p] = StatUtils.mean(data[p]);
            sd[p] = Math.sqrt(StatUtils.variance(data[p], means[p]));
        }
        double[][] data_W = new double[sample_size][num_var];
        for (int k = 0; k < sample_size; k++) {
            for (int p = 0; p < num_var; p++) {
                data_W[k][p] = (data[p][k] - means[p]) / sd[p];
            }
        }
        Array2DRowRealMatrix data_W_matrix = new Array2DRowRealMatrix(data_W);
        Array2DRowRealMatrix data_W_matrixT = (Array2DRowRealMatrix) (data_W_matrix.transpose());
        RealMatrix kinship = ((data_W_matrix.multiply(data_W_matrixT)).scalarMultiply(1.0 / (double) num_var));//  new double[sample_size][sample_size];
//		System.out.println("Local kinship (before rescale) has been calculted.");
        return kinship;
    }

    /*
     * Just return X'X without doing normalization.
     * The reason why this is needed is that after transformation, one should
     * not do normalization again.
     */
    public static RealMatrix calculate_RRM_local_kinship_RealMatrix_noNomalization(double[][] data) {
        if (data == null) {
            return null;
        }
        Array2DRowRealMatrix data_W_matrix = new Array2DRowRealMatrix(data);
        Array2DRowRealMatrix data_W_matrixT = (Array2DRowRealMatrix) (data_W_matrix.transpose());
        RealMatrix kinship = data_W_matrixT.multiply(data_W_matrix);//  new double[sample_size][sample_size];
//		System.out.println("Local kinship (before rescale) has been calculted.");
        return kinship;
    }

    public static void run_nam_imputation(String pedigreefiles_folder, String crossid, String RILpheno,
                                          String geno250k, String output_prefix, int block_size) {
        try {
            File[] pedfiles = (new File(pedigreefiles_folder)).listFiles();
            RILs[] ril_info = new RILs[pedfiles.length];
            Cross_info cross = new Cross_info(crossid);
            for (int i = 0; i < pedfiles.length; i++) {
                ril_info[i] = new RILs(pedfiles[i], RILpheno, cross);
            }
            Founder founders = new Founder(geno250k, cross);
            Imputation.imputefull(founders, ril_info, output_prefix + ".pheno.tsv", output_prefix + ".geno.csv");
            System.out.println("Imputation finished:\n" +
                    "Genotype file: " + output_prefix + ".geno.csv\n" + "Phenotype file: " + output_prefix + ".pheno.tsv");
            VariantsDouble.char2num(output_prefix + ".geno.csv", output_prefix + ".geno.num.csv");
            VariantsDouble.importCSV(output_prefix + ".geno.num.csv", output_prefix + ".geno.num.csv.hdf5", block_size);

        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    /*
     * return an indicator array in which "true" stands for needed variants that passed the maf_threshold
     */
    public static boolean[] find_variants_above_maf(double[][] data, double allele_freq_threshold) {
        int var_num = data.length;
        boolean[] indicators = new boolean[var_num];
        for (int k = 0; k < var_num; k++) {
            double count = 0;
            for (int i = 0; i < data[k].length; i++)
                count += data[k][i];
            double maf = count / (data[k].length * 2.0);
            //if(maf>0.5)maf=1-maf;
            if (maf > allele_freq_threshold) indicators[k] = true;
        }
        return indicators;
    }

    /*
     * return an indicator array in which "true" stands for needed variants that passed the maf_threshold
     */
    public static boolean[] find_variants_below_maf(double[][] data, double allele_freq_threshold) {
        int var_num = data.length;
        boolean[] indicators = new boolean[var_num];
        for (int k = 0; k < var_num; k++) {
            double count = 0;
            for (int i = 0; i < data[k].length; i++)
                count += data[k][i];
            double maf = count / (data[k].length * 2.0);
            //if(maf>0.5)maf=1-maf;
            if (maf < allele_freq_threshold)
                indicators[k] = true;
        }
        return indicators;
    }

    public static void char2num(String input, String output) {
        //TODO
        System.out.println("Convert genotype CSV file coded with ACGT/RWSYKM into CSV file coded with 0,1,2" +
                " (0=major hom; 1=het; 2=minor hom)\n" +
                "Lines containing more than 3 alleles will be omitted.");
        HashSet<String> homs = new HashSet<>();
        HashSet<String> hets = new HashSet<>();
        homs.add("A");
        homs.add("C");
        homs.add("G");
        homs.add("T");
        hets.add("R");
        hets.add("W");
        hets.add("S");
        hets.add("Y");
        hets.add("K");
        hets.add("M");
        try {
            BufferedReader br = new BufferedReader(new FileReader(input));
            String line = br.readLine();
            BufferedWriter bw = new BufferedWriter(new FileWriter(output));
            bw.write(line + "\n");//header
            int sample_size = line.split(",").length - 2;
            System.out.println("Sample_size=" + sample_size);
            line = br.readLine();
            int un = 0, one = 0, three = 0, four = 0;
            while (line != null) {
                HashMap<String, Integer> allele_counts = new HashMap<>();
                String[] temp = line.split(",");
                if (temp.length - 2 != sample_size) { // length wrong
                    System.out.println("Line length != sample_size");
                    line = br.readLine();
                    continue;
                }
                boolean undefined_char = false;
                for (int i = 2; i < temp.length; i++) {
                    if (homs.contains(temp[i])) {
                        myFileFunctions.FileFunc.add2_counts_hashmap(allele_counts, temp[i], 2);
                    } else if (hets.contains(temp[i])) {
                        String[] alleles = one2two(temp[i]).split(" ");
                        for (int j = 0; j < 2; j++)
                            myFileFunctions.FileFunc.add2_counts_hashmap(allele_counts, alleles[j], 1);
                    } else {
                        //System.out.println("Char not ACGT/RWSYKM: "+temp[i]);
                        undefined_char = true;
                        //break;
                    }
                }
                if (undefined_char) un++;
                if (allele_counts.size() != 2) {
                    line = br.readLine();
                    if (allele_counts.size() == 1) one++;
                    else if (allele_counts.size() == 3) three++;
                    else if (allele_counts.size() == 4) four++;
                    continue;
                }
                // output this line
                String[] alleles = new String[2];
                int index = 0;
                for (String the_allele : allele_counts.keySet()) {
                    alleles[index++] = the_allele;
                }
                String major = alleles[0], minor = alleles[1];
                if (allele_counts.get(major) < allele_counts.get(minor)) {
                    String buff = major;
                    major = minor;
                    minor = buff;
                }
                bw.write(temp[0] + "," + temp[1]);
                for (int i = 2; i < temp.length; i++) {
                    if (temp[i].equals(major)) bw.write("," + 0);
                    else if (temp[i].equals(minor)) bw.write("," + 2);
                    else {
                        String separated = one2two(temp[i]);
                        if (separated.equals(major + " " + minor) || separated.equals(minor + " " + major)) {
                            bw.write("," + 1);
                        } else {
                            bw.write(",0");
                            //System.out.println("Error: het not consistent to homs! "+temp[i]+":"+major+"/"+minor+"Exit.");
                        }
                    }
                }
                bw.write("\n");
                line = br.readLine();
            }
            bw.close();
            System.out.println("Finished coverting: \n" + "Underfiend chars (not ACGT/RWSYKM)=" + un + "\n" +
                    "One Allele=" + one + "\n" + "Three Allele=" + three + "\n" + "Four Alleler=" + four + "\n");
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public static void tped2csv_char(String input_tfam, String input_tped, String output_csv) {
        try {
            BufferedReader br = new BufferedReader(new FileReader(input_tfam));
            String line = br.readLine();
            int sample_size = 0;
            BufferedWriter bw = new BufferedWriter(new FileWriter(output_csv));
            bw.write("CHR,LOC");
            while (line != null) {
                sample_size++;
                bw.write("," + line.split(" ")[1]);
                line = br.readLine();
            }
            bw.write("\n");
            br = new BufferedReader(new FileReader(input_tped));
            line = br.readLine();
            int skip = 0, written = 0;
            while (line != null) {
                String[] temp = line.split(" ");
                if ((temp.length - 4) / 2 != sample_size) {
                    System.out.println("Sample Size wrong for the line: " + (temp.length - 4) / 2 + "!=" + sample_size);
                }
                if (temp[0].equals("0") || temp[3].equals("0") || Integer.parseInt(temp[0]) >= 23) {
                    skip++;
                } else {
                    written++;
                    bw.write(temp[0] + "," + temp[3]);
                    for (int k = 0; k < sample_size; k++) {
                        bw.write("," + two2one(temp[4 + k * 2] + " " + temp[4 + k * 2 + 1]));
                    }
                    bw.write("\n");
                }
                line = br.readLine();
            }
            bw.close();
            System.out.println("Coverted tped to csv file: " + written + " written (" + skip + " skipped)");
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public static void tped2csv_num(String input_tfam, String input_tped,
                                    String csv_char, String csv_num) {
        tped2csv_char(input_tfam, input_tped, csv_char);
        char2num(csv_char, csv_num);
    }

    public static String one2two(String snp) {
        String result = "N N";
        if (snp.equals("A") || snp.equals("C") || snp.equals("G") || snp.equals("T"))
            result = snp + " " + snp;
        else {
            //		System.out.println("HET!");
            if (snp.equals("M")) result = "A C";
            else if (snp.equals("R")) result = "A G";
            else if (snp.equals("W")) result = "A T";
            else if (snp.equals("S")) result = "C G";
            else if (snp.equals("Y")) result = "C T";
            else if (snp.equals("K")) result = "G T";
//			else System.out.println(snp+": Not a correct SNP!");
        }
        return result;
    }

    public static String two2one(String snp) {
        String result = "N";
        if (snp.equals("A A") || snp.equals("C C") || snp.equals("G G") || snp.equals("T T"))
            result = snp.split(" ")[0];
        else {
            //		System.out.println("HET!");
            if (snp.equals("A C") || snp.equals("C A")) result = "M";
            else if (snp.equals("A G") || snp.equals("G A")) result = "R";
            else if (snp.equals("A T") || snp.equals("T A")) result = "W";
            else if (snp.equals("C G") || snp.equals("G C")) result = "S";
            else if (snp.equals("C T") || snp.equals("T C")) result = "Y";
            else if (snp.equals("G T") || snp.equals("T G")) result = "K";
//			else System.out.println(snp+": Not a correct SNP!");
        }
        return result;
    }

    /*
     * Some times, one wants to store the transformed data for further analysis:
     * TODO: never run this function yet!!
     */
    public void transform_all_data(RealMatrix transform_matrix, String new_file_hdf5) {
        if ((new LUDecomposition(transform_matrix)).getDeterminant() == 0) {
            System.out.println("Transform matrix is not full-ranked. Function Returned.");
            return;
        }
        try {
            File hdf5_old = new File(new_file_hdf5);
            if (hdf5_old.exists()) {
                hdf5_old.delete();
                System.out.println(new_file_hdf5 + " already existed. JAWAMix5 has deleted it...");
            }
            IHDF5Writer writer = HDF5Factory.open(new_file_hdf5);
            // CODE FOR constructing above structure (all groups and tables and importing data)
            writer.object().createGroup("/genotype");
            writer.int32().setAttr("/genotype", "sample_size", this.sample_size);
            writer.int32().setAttr("/genotype", "block_size", this.block_size);
            writer.int32().setAttr("/genotype", "num_chrs", this.num_chrs);
            writer.writeStringArray("/genotype/sample_ids", this.sample_ids);
            writer.bool().setAttr("/genotype", "transformed", true);
            writer.writeDoubleArray("/genotype/intercepts", transform_matrix.operate(this.intercepts));

            int[] num_blocks_per_chr = new int[num_chrs];
            for (int chr = 0; chr < num_chrs; chr++) {
                num_blocks_per_chr[chr] = num_sites[chr] / block_size + 1;
                writer.object().createGroup("/genotype/chr" + (1 + chr));
                writer.writeIntArray("/genotype/chr" + (1 + chr) + "/var_pos", this.locations[chr]);
                writer.writeIntArray("/genotype/chr" + (1 + chr) + "/var_mafc", this.mafc[chr]);
                writer.float64().createMatrix("/genotype/chr" + (1 + chr) + "/indi_fast", this.sample_size, this.num_sites[chr],
                        this.sample_size, this.block_size);
                writer.float64().createMatrix("/genotype/chr" + (1 + chr) + "/position_fast", this.num_sites[chr], this.sample_size,
                        this.block_size, this.sample_size);
            }
            writer.writeIntArray("/genotype/num_blocks_per_chr", num_blocks_per_chr);
            // Read data block to fill matrix:
            System.out.println("Transfering data to create HDF5 file for transformed genotypes Started.");

            for (int chr = 0; chr < this.num_chrs; chr++) {
                int current_block = 0;
                for (HDF5MDDataBlock<MDDoubleArray> block : this.position_fast_blocks[chr]) {
                    double[][] data4thisblock = block.getData().toMatrix();
                    double[][] transformed = FaSTLMM.transform_X_by_Ut(transform_matrix, data4thisblock);
                    writer.float64().writeMatrixBlockWithOffset("/genotype/chr" + (chr + 1) + "/position_fast",
                            data4thisblock, data4thisblock.length, sample_size, current_block * block_size, 0L);
                    current_block++;
                }
                System.out.println("finished Chr" + (chr + 1) + ".");
            }
            System.out.println("Finished transforming and writing genotypes. Now the HDF5 data is ready to use.");
            writer.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public int[][] find_vars_in_maf_range(double min_maf, double max_maf) {
        int[][] results = new int[this.num_chrs][];
        for (int chr = 0; chr < this.num_chrs; chr++) {
            ArrayList<Integer> the_qualified = new ArrayList<>();
            for (int i = 0; i < this.num_sites[chr]; i++) {
                double maf = ((double) this.mafc[chr][i]) / (this.sample_size * 2.0);
                if (maf >= min_maf && maf <= max_maf) the_qualified.add(i);
            }
            results[chr] = myFileFunctions.FileFunc.arraylist2arrayInteger(the_qualified);
        }
        return results;
    }

    public int[][] find_fixed_num_of_vars_in_maf_range(double min_maf, double max_maf, int num_vars_to_simulate) {
        int[] breaks = new int[this.num_chrs];
        breaks[0] = 0;
        for (int k = 1; k < this.num_chrs; k++) {
            breaks[k] = breaks[k - 1] + this.locations[k - 1].length;
        }
        if (!(breaks[this.num_chrs - 1] + this.locations[this.num_chrs - 1].length == this.num_sites_total)) {
            System.out.println("Breaks or this.num_sites_total incorrect!");
        }
        int[][] results = new int[this.num_chrs][];
        ArrayList<Integer>[] the_qualified = new ArrayList[this.num_chrs];
        for (int chr = 0; chr < this.num_chrs; chr++) {
            the_qualified[chr] = new ArrayList<>();
        }
        int simulated = 0;
        while (simulated < num_vars_to_simulate) {
            long location = (long) (myMathLib.Test.randomNumber() * this.num_sites_total);
            int the_location = -100;
            if (location == this.num_sites_total) location--;
            if (location == 0) location++;
            int the_chr = 0;
            for (int chr = 0; chr < this.num_chrs; chr++) {
                if (location < breaks[chr]) {
                    the_chr = chr - 1;
                    the_location = (int) (location - breaks[chr - 1]);
                    break;
                }
            }
            if (the_location == -100) { // not found, so the last chr
                the_chr = this.num_chrs - 1;
                the_location = (int) (location - breaks[this.num_chrs - 1]);
            }
            double maf = ((double) this.mafc[the_chr]
                    [the_location]) / (this.sample_size * 2.0);
            if (maf >= min_maf && maf <= max_maf) {
                the_qualified[the_chr].add(the_location);
//				System.out.println(simulated+": chr"+the_chr+"_"+the_location+"_"+maf);
                simulated++;
            }
        }
        for (int chr = 0; chr < this.num_chrs; chr++) {
            results[chr] = myFileFunctions.FileFunc.arraylist2arrayInteger(the_qualified[chr]);
        }
        return results;
    }

    public int[] find_vars_in_maf_range_in_region(int chr, int start_location, int end_location,
                                                  double min_maf, double max_maf) {
        int start_index = Arrays.binarySearch(this.locations[chr], start_location);
        if (start_index < 0) {// no_found, and make use of the insertion-point returned by the binarySearch:
            start_index = -(start_index + 1);
        }
        int end_index = Arrays.binarySearch(this.locations[chr], end_location);
        if (end_index < 0) {// no_found, and make use of the insertion-point returned by the binarySearch:
            end_index = -(end_index + 1) - 1;
        }
        ArrayList<Integer> the_qualified = new ArrayList<>();
        for (int i = start_index; i <= end_index; i++) {
            double maf = ((double) this.mafc[chr][i]) / (this.sample_size * 2.0);
            if (maf >= min_maf && maf <= max_maf) {
                the_qualified.add(i);
            }
        }
        return myFileFunctions.FileFunc.arraylist2arrayInteger(the_qualified);
    }

    /*
     * Return the variants in a region:
     * chr starts from ZERO!
     */
    // CHECKED - PATHUM
    public double[][] load_variants_in_region(int chr, int start_location, int end_location) {
        int start_index = Arrays.binarySearch(this.locations[chr], start_location);
        if (start_index < 0) {// start_location has not been found, thus we make use of the insertion-point returned by
            // the binarySearch:
            // The insertion point is defined as the point at which the key would be inserted into the array: the index
            // of the first element greater than the key, or a.length if all elements in the array are less than the
            // specified key.
            start_index = -(start_index + 1);
        }
        if (start_index >= this.locations[chr].length) return (new double[0][]);
        int end_index = Arrays.binarySearch(this.locations[chr], end_location);
        if (end_index < 0) {// end_location has not been found, thus we make use of the insertion-point returned by the
            // binarySearch:
            end_index = -(end_index + 1) - 1;
        }
        if (end_index < 0) return (new double[0][]);
        return reader.float64().readMatrixBlockWithOffset(position_fast_paths[chr],
                (end_index - start_index + 1), this.sample_size,
                start_index, 0);
    }

    /*
     * chr starts from 0.
     */
    double[] load_one_variant_by_location(int chr, int location) {
        int index = Arrays.binarySearch(this.locations[chr], location);
        if (index < 0) {
            System.out.println("Error: Variant location not in the list");
            return null;
        }
        return load_one_variant_by_index(chr, index);
    }

    public double[] load_one_variant_by_index(int chr, int index) {
//    	System.out.println(chr+"_"+  index);
        return reader.float64().readMatrixBlockWithOffset(position_fast_paths[chr],
                1, this.sample_size, index, 0)[0];
    }

//	/*
//	 * S_new=(n-1)S/Tr(PSP), where P= I-11'/n and 1 is the vector of ones.
//	 */
//	public static double[][] re_scale_kinship_matrix(double[][] S){
//		int n=S.length;
//		double[][] S_new=new double[n][n];
//		DenseDoubleMatrix2D matrix_P= new DenseDoubleMatrix2D(n,n);
//		matrix_P.assign(-1.0/n);
//		for(int i=0;i<n;i++){
//			matrix_P.setQuick(i, i, (1.0-1.0/n));
//		}
//		DenseDoubleMatrix2D matrix_S= new DenseDoubleMatrix2D(S);
//		DenseDoubleMatrix2D PS= new DenseDoubleMatrix2D(n,n);
//		DenseDoubleMatrix2D PSP= new DenseDoubleMatrix2D(n,n);
//		matrix_P.zMult(matrix_S,PS);
//		PS.zMult(matrix_P, PSP);
//		double trace=0;
//		for(int k=0;k<n;k++)	trace+=PSP.get(k, k);
//		System.out.println("Rescale coef="+(n-1)/trace);
//		for(int i=0;i<n;i++){
//			for(int j=0;j<n;j++){
//				S_new[i][j]=S[i][j]*(n-1)/trace;
//			}
//		}
//		return S_new;
//	}

//	public static void re_scale_kinship_matrix(String in_kinship_file, String out_kinship_file){
//		String[] separators={","," ","\t"};
//		try{
//			String sep=null;
//			int sample_size=-1;
//			if(new File(in_kinship_file).exists()){
//				BufferedReader br = new BufferedReader(new FileReader(in_kinship_file));
//				String header=br.readLine();
//				String line=br.readLine();
//				for(int k=0;k<separators.length;k++){
//					if(line.split(separators[k]).length>sample_size){
//						sample_size=line.split(separators[k]).length;
//						sep=separators[k];
//					}
//				}//System.out.println("Sep=\""+sep+"\"");
////				System.out.println("Sample size=\""+sample_size+"\"");
//				double[][] old_kinship=new double[sample_size][sample_size];
//				int line_index=0;
//				while(line!=null){
//					String[] temp=line.split(sep);
//					for(int i=0;i<temp.length;i++){
//						old_kinship[line_index][i]=Double.parseDouble(temp[i]);
//					}
//					line=br.readLine();
//					line_index++;
//				}
//
//				double[][] new_kinship=re_scale_kinship_matrix(old_kinship);
//				BufferedWriter bw= new BufferedWriter(new FileWriter(out_kinship_file));
//				bw.write(header+"\n");
//				for(int i=0;i<sample_size;i++){
//					for(int j=0;j<sample_size-1;j++){
//						bw.write(new_kinship[i][j]+",");
//					}bw.write(new_kinship[i][sample_size-1]+"\n");
//				}
//				bw.close();
//			}
//		}catch(Exception e){e.printStackTrace();}
//	}

    public double[][] readDoubleMatrixBlock_position_fast(int chr, long blockNumberX) {
        if (blockNumberX != this.num_blocks_per_chr[chr] - 1) {
            return reader.float64().readMatrixBlock(position_fast_paths[chr], block_size, this.sample_size,
                    blockNumberX, 0);
        } else {
//			this.position_fast_blocks[chr].
//			return reader.readDoubleMatrixBlockWithOffset(position_fast_paths[chr], block_size, this.sample_size,
//					block_size*(this.num_blocks_per_chr[chr]-1),0);
            return reader.float64().readMatrixBlockWithOffset(position_fast_paths[chr],
                    (int) (this.num_sites[chr] - (this.num_blocks_per_chr[chr] - 1) * block_size), this.sample_size,
                    block_size * (this.num_blocks_per_chr[chr] - 1), 0);
        }
    }

    public void output2csv(String csv_file) {
        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(csv_file));
            bw.write("Chromosome,Positions");
            for (int k = 0; k < this.sample_size; k++) bw.write("," + this.sample_ids[k]);
            bw.write("\n");
            for (int chr = 0; chr < this.num_chrs; chr++) {
                int current_variant_index = 0;
                for (HDF5MDDataBlock<MDDoubleArray> block : this.position_fast_blocks[chr]) {
//					System.out.println((chr+1)+":"+ArrayUtils.toString(block.getIndex()) + " -> "+ block.getData().toString());
                    double[][] data4thisblock = block.getData().toMatrix();
                    for (int i = 0; i < data4thisblock.length; i++) {
                        bw.write((chr + 1) + "," + this.locations[chr][current_variant_index]);
                        for (int k = 0; k < this.sample_size; k++) {
                            if (Double.compare(data4thisblock[i][k], (int) (data4thisblock[i][k] + 0.01)) == 0)
                                bw.write("," + (int) (data4thisblock[i][k] + 0.01));
                            else
                                bw.write("," + data4thisblock[i][k]);
                        }
                        bw.write("\n");
                        current_variant_index++;
                    }
                }
            }
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public double[][] read_in_kinship(String kinship_file) {
        double[][] kinship = new double[this.sample_size][this.sample_size];
        String[] separators = {",", " ", "\t"};
        String sep = null;
        try {
            BufferedReader br = new BufferedReader(new FileReader(kinship_file));
            String line = br.readLine();//header
            for (int k = 0; k < separators.length; k++) {
                if (line.split(separators[k]).length == this.sample_size) {
                    sep = separators[k];
                }
            }
            int line_index = 0;
            while (line != null) {
                String[] temp = line.split(sep);
                for (int k = 0; k < temp.length; k++) {
                    kinship[line_index][k] = Double.parseDouble(temp[k]);
                }
                line_index++;
                line = br.readLine();
            }
            if (line_index != this.sample_size) {
                System.out.println("Kinship file number of line wrong!");
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        return kinship;
    }

    /*
     * IBS for WG with NaN handled.
     */
    public void calculate_raw_ibs_kinship(String output_file, double scale, double min_MAF) {
        try {
            double[][] kinship = new double[this.sample_size][this.sample_size];
            double[][] total_num_var_used = new double[this.sample_size][this.sample_size];
            for (int chr = 0; chr < this.num_chrs; chr++) {
                for (HDF5MDDataBlock<MDDoubleArray> block : this.position_fast_blocks[chr]) {
                    //double start_time = System.nanoTime();
                    double[][] data4thisblock = block.getData().toMatrix();
                    for (int var_index = 0; var_index < data4thisblock.length; var_index++) {
                        if (KinshipMatrix.maf(data4thisblock[var_index], scale) < min_MAF) continue;
                        for (int i = 0; i < sample_size; i++) {
                            for (int j = i + 1; j < sample_size; j++) {
                                if ((!Double.isNaN(data4thisblock[var_index][i])) && (!Double.isNaN(data4thisblock[var_index][j]))) {
                                    kinship[i][j] = kinship[i][j] + (scale - Math.abs(data4thisblock[var_index][i] - data4thisblock[var_index][j]));
                                    total_num_var_used[i][j]++;
                                }
                            }
                        }
                    }
                    //double end_time = System.nanoTime();
                    //System.out.println("Time taken for 1 block in chr " + chr + ": " + ((end_time - start_time) / 1000000000));
                }
                System.out.println("Finished Chr" + (chr + 1));
            }
            for (int i = 0; i < sample_size; i++) {
                kinship[i][i] = 1;
                for (int j = i + 1; j < sample_size; j++) {
                    kinship[i][j] = kinship[i][j] / (total_num_var_used[i][j] * scale);
                    kinship[j][i] = kinship[i][j];
                }
            }
            BufferedWriter bw = new BufferedWriter(new FileWriter(output_file));
            for (int i = 0; i < sample_size - 1; i++) {
                bw.write(this.sample_ids[i] + ",");
            }
            bw.write(this.sample_ids[sample_size - 1] + "\n");
            for (int i = 0; i < sample_size; i++) {
                for (int j = 0; j < sample_size - 1; j++) {
                    bw.write(kinship[i][j] + ",");
                }
                bw.write(kinship[i][sample_size - 1] + "\n");
            }
            bw.close();
            System.out.println("Global kinship (before rescale) has been written to " + output_file);


        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void calculate_raw_ibs_kinship_multithreaded(String output_file, double scale, double min_MAF) {
        try {
            KinshipMatrixMultiThreaded Kinship_Matrix_Thread = new KinshipMatrixMultiThreaded(this.sample_size);

            ExecutorService es = Executors.newFixedThreadPool(1);

            for (int chr = 0; chr < this.num_chrs; chr++) {
                KinshipMultiCalculator kM = new KinshipMultiCalculator(chr, position_fast_blocks[chr], scale, min_MAF, this.sample_size, Kinship_Matrix_Thread);
                es.execute(kM);
            }
            es.shutdown();
            while (!es.isTerminated()) {
            }

            for (int i = 0; i < sample_size; i++) {
                Kinship_Matrix_Thread.kinship[i][i] = 1;
                for (int j = i + 1; j < sample_size; j++) {
                    Kinship_Matrix_Thread.kinship[i][j] = Kinship_Matrix_Thread.kinship[i][j] / (Kinship_Matrix_Thread.total_num_var_used[i][j] * scale);
                    Kinship_Matrix_Thread.kinship[j][i] = Kinship_Matrix_Thread.kinship[i][j];
                }
            }
            BufferedWriter bw = new BufferedWriter(new FileWriter(output_file));
            for (int i = 0; i < sample_size - 1; i++) {
                bw.write(this.sample_ids[i] + ",");
            }
            bw.write(this.sample_ids[sample_size - 1] + "\n");
            for (int i = 0; i < sample_size; i++) {
                for (int j = 0; j < sample_size - 1; j++) {
                    bw.write(Kinship_Matrix_Thread.kinship[i][j] + ",");
                }
                bw.write(Kinship_Matrix_Thread.kinship[i][sample_size - 1] + "\n");
            }
            bw.close();
            System.out.println("Global kinship (before rescale) has been written to " + output_file);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    /*
     * RRM for WG with NaN handled
     */
    public void calculate_WG_RRM_kinship(String output_file, double scale, double min_MAF) {
        try {
            double[][] kinship = new double[this.sample_size][this.sample_size];
            double[][] total_num_var_used = new double[this.sample_size][this.sample_size];
            for (int chr = 0; chr < this.num_chrs; chr++) {
                for (HDF5MDDataBlock<MDDoubleArray> block : this.position_fast_blocks[chr]) {
                    double[][] data4thisblock = block.getData().toMatrix();
                    for (int var_index = 0; var_index < data4thisblock.length; var_index++) {
                        if (KinshipMatrix.maf(data4thisblock[var_index], scale) < min_MAF) continue;
                        double mean = myMathLib.StatFuncs.mean_NaN(data4thisblock[var_index]);
                        if (Double.isNaN(mean)) {
                            System.out.println("All subjects are NaN at a SNP in Chr" + (chr + 1));
                            continue;
                        }
                        double sd = Math.sqrt(myMathLib.StatFuncs.var_NaN(data4thisblock[var_index], mean));
                        double[] new_snp = new double[this.sample_size];
                        for (int k = 0; k < this.sample_size; k++)
                            new_snp[k] = (data4thisblock[var_index][k] - mean) / sd;
                        for (int i = 0; i < sample_size; i++) {
                            for (int j = i; j < sample_size; j++) {
                                if ((!Double.isNaN(new_snp[i])) && (!Double.isNaN(new_snp[j]))) {
                                    kinship[i][j] = kinship[i][j] + (new_snp[i] * new_snp[j]);
                                    total_num_var_used[i][j]++;
                                }
                            }
                        }
                    }
                }
                System.out.println("Finished Chr" + (chr + 1));
            }
            for (int i = 0; i < sample_size; i++) {
                for (int j = i; j < sample_size; j++) {
                    kinship[i][j] = kinship[i][j] / total_num_var_used[i][j];
                    kinship[j][i] = kinship[i][j];
                }
            }
            BufferedWriter bw = new BufferedWriter(new FileWriter(output_file));
            for (int i = 0; i < sample_size - 1; i++) {
                bw.write(this.sample_ids[i] + ",");
            }
            bw.write(this.sample_ids[sample_size - 1] + "\n");
            for (int i = 0; i < sample_size; i++) {
                for (int j = 0; j < sample_size - 1; j++) {
                    bw.write(kinship[i][j] + ",");
                }
                bw.write(kinship[i][sample_size - 1] + "\n");
            }
            bw.close();
            System.out.println("Global kinship (before rescale) has been written to " + output_file);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    /*
     * RRM for WG of ids selected, with NaN handled
     */
    public void calculate_WG_RRM_kinship(String output_file, ArrayList<String> ids_selected) {
        int[] indexes = new int[ids_selected.size()];
        for (int k = 0; k < indexes.length; k++) {
            indexes[k] = this.sample_id2index.get(ids_selected.get(k));
        }
        try {
            double[][] kinship = new double[indexes.length][indexes.length];
            double[][] total_num_var_used = new double[indexes.length][indexes.length];
            for (int chr = 0; chr < this.num_chrs; chr++) {
                long num_NaN = 0;
                for (HDF5MDDataBlock<MDDoubleArray> block : this.position_fast_blocks[chr]) {
                    double[][] data4thisblock = block.getData().toMatrix();
                    for (int var_index = 0; var_index < data4thisblock.length; var_index++) {
                        double[] the_selected_data = new double[indexes.length];
                        for (int i = 0; i < indexes.length; i++) {
                            the_selected_data[i] = data4thisblock[var_index][indexes[i]];
                        }
                        double mean = myMathLib.StatFuncs.mean_NaN(the_selected_data);
                        if (Double.isNaN(mean)) {
                            System.out.println("All subjects are NaN at a SNP in Chr" + (chr + 1));
                            num_NaN++;
                            continue;
                        }
                        double sd = Math.sqrt(myMathLib.StatFuncs.var_NaN(the_selected_data, mean));
                        double[] new_snp = new double[indexes.length];
                        for (int k = 0; k < ids_selected.size(); k++) new_snp[k] = (the_selected_data[k] - mean) / sd;
                        for (int i = 0; i < ids_selected.size(); i++) {
                            for (int j = i; j < ids_selected.size(); j++) {
                                if ((!Double.isNaN(new_snp[i])) && (!Double.isNaN(new_snp[j]))) {
                                    kinship[i][j] = kinship[i][j] + (new_snp[i] * new_snp[j]);
                                    total_num_var_used[i][j]++;
                                }
                            }
                        }
                    }
                }
                System.out.println("Finished Chr" + (chr + 1) + ": " + num_NaN + " SNPs are all NaN");
            }
            for (int i = 0; i < indexes.length; i++) {
                for (int j = i; j < indexes.length; j++) {
                    kinship[i][j] = kinship[i][j] / total_num_var_used[i][j];
                    kinship[j][i] = kinship[i][j];
                }
            }
            BufferedWriter bw = new BufferedWriter(new FileWriter(output_file));
            for (int i = 0; i < ids_selected.size() - 1; i++) {
                bw.write(ids_selected.get(i) + ",");
            }
            bw.write(ids_selected.get(ids_selected.size() - 1) + "\n");
            for (int i = 0; i < ids_selected.size(); i++) {
                for (int j = 0; j < ids_selected.size() - 1; j++) {
                    bw.write(kinship[i][j] + ",");
                }
                bw.write(kinship[i][ids_selected.size() - 1] + "\n");
            }
            bw.close();
            System.out.println("Global RRM kinship has been written to " + output_file);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    /*
     * RRM for WG of ids selected, with NaN handled
     */
    public void calculate_WG_IBS_kinship(String output_file, ArrayList<String> ids_selected) {
        double scale = 2;
        int[] indexes = new int[ids_selected.size()];
        for (int k = 0; k < indexes.length; k++) {
            indexes[k] = this.sample_id2index.get(ids_selected.get(k));
        }
        try {
            double[][] kinship = new double[indexes.length][indexes.length];
            double[][] total_num_var_used = new double[indexes.length][indexes.length];
            for (int chr = 0; chr < this.num_chrs; chr++) {
                long num_NaN = 0;
                for (HDF5MDDataBlock<MDDoubleArray> block : this.position_fast_blocks[chr]) {
                    double[][] data4thisblock = block.getData().toMatrix();
                    for (int var_index = 0; var_index < data4thisblock.length; var_index++) {
                        double[] the_selected_data = new double[indexes.length];
                        int NaNs = 0;
                        for (int i = 0; i < indexes.length; i++) {
                            the_selected_data[i] = data4thisblock[var_index][indexes[i]];
                            if (Double.isNaN(data4thisblock[var_index][indexes[i]])) NaNs++;
                        }
                        if (NaNs == indexes.length) {
                            System.out.println("All subjects are NaN at a SNP in Chr" + (chr + 1));
                            num_NaN++;
                            continue;
                        }
                        for (int i = 0; i < ids_selected.size(); i++) {
                            for (int j = i; j < ids_selected.size(); j++) {
                                if ((!Double.isNaN(the_selected_data[i])) && (!Double.isNaN(the_selected_data[j]))) {
                                    kinship[i][j] = kinship[i][j] +
                                            (scale - Math.abs(the_selected_data[i] - the_selected_data[j]));
                                    total_num_var_used[i][j]++;
                                }
                            }
                        }
                    }
                }
                System.out.println("Finished Chr" + (chr + 1) + ": " + num_NaN + " SNPs are all NaN");
            }
            for (int i = 0; i < indexes.length; i++) {
                for (int j = i; j < indexes.length; j++) {
                    kinship[i][j] = kinship[i][j] / (total_num_var_used[i][j] * scale);
                    kinship[j][i] = kinship[i][j];
                }
            }
            BufferedWriter bw = new BufferedWriter(new FileWriter(output_file));
            for (int i = 0; i < ids_selected.size() - 1; i++) {
                bw.write(ids_selected.get(i) + ",");
            }
            bw.write(ids_selected.get(ids_selected.size() - 1) + "\n");
            for (int i = 0; i < ids_selected.size(); i++) {
                for (int j = 0; j < ids_selected.size() - 1; j++) {
                    bw.write(kinship[i][j] + ",");
                }
                bw.write(kinship[i][ids_selected.size() - 1] + "\n");
            }
            bw.close();
            System.out.println("Global IBS kinship (before rescale) has been written to " + output_file);
        } catch (Exception e) {
            e.printStackTrace();
        }
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

}


