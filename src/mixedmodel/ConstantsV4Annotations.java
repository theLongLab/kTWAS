package mixedmodel;

import java.util.HashMap;


public class ConstantsV4Annotations {

	public static String[] orders={"#CHROM","POS","ID","REF","ALT","QUAL","AA","AC","AF","AFR_AF","AMR_AF","AN","ASN_AF","AVGPOST",
		"CIEND","CIPOS","END","ERATE","EUR_AF","HOMLEN","HOMSEQ","LDAF","RSQ","SNPSOURCE","SVLEN","SVTYPE","THETA","VT","VA","NAC"};
	
	public static String[] acgt={"A","C","G","T","*"};
	
	public static String[] SNPsFunctionsal= {"nonsynonymous", "synonymous", "prematureStop", "removedStop", 
		"spliceOverlap","*"};
	public static String[] IndelsFunctionsal= {"deletionNFS", "insertionNFS", "deletionFS", "insertionFS", "startOverlap",
		"endOverlap", "spliceOverlap", "multiExonHit","*"};	
	public static String[] non_coding_functional={"miRNA","snRNA","snoRNA","rRNA","lincRNA","misc_RNA","UTR","PGENE",
		"TFPEAK","TFMOTIF","ENHANCER","*"};
	
	public static String[] var_type={"SNP","INDEL","SV","*"};
	
	public static HashMap<String, Byte> map_name2index = generate_map(ConstantsV4Annotations.orders);
	
	public static HashMap<String, Byte> actg2index = generate_map_nocase(ConstantsV4Annotations.acgt);	
	public static HashMap<String, Byte> snp_fucntion2index=generate_map(ConstantsV4Annotations.SNPsFunctionsal);
	public static HashMap<String, Byte> indel_fucntion2index=generate_map(ConstantsV4Annotations.IndelsFunctionsal);
	public static HashMap<String, Byte> noncoding_function2index=generate_map(ConstantsV4Annotations.non_coding_functional);
	public static HashMap<String, Byte> vt2index=generate_map(ConstantsV4Annotations.var_type);
	
	public static HashMap<String, Byte> generate_map(String[] type_list){
		HashMap<String, Byte> map=new HashMap<String, Byte>();
		for(byte k=0;k<type_list.length;k++){
			map.put(type_list[k], k);
		}
		return map;
	}

	public static HashMap<String, Byte> generate_map_nocase(String[] type_list){
		HashMap<String, Byte> map=new HashMap<String, Byte>();
		for(byte k=0;k<type_list.length;k++){
			map.put(type_list[k], k);
			map.put(type_list[k].toLowerCase(), k);
		}
		return map;
	}
}
