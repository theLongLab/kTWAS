package mixedmodel;

import java.util.HashMap;

public class Annotation_g1k {
	
	public byte REF;
	public byte ALT;
	public int QUAL;	
	public byte AA; // ancestral allele
	
	public int AC; //Alternate Allele Count
	public int AN; //Total Allele Count
	public double AF; // AC/AN
	
	public double AFR_AF;	
	public double AMR_AF;
	public double ASN_AF;
	public double EUR_AF;
	
	public byte VT;
//	public String VA="*";
	public byte NAC=-1;

	public byte AlleleNumber=-1;
	public int GeneName_index=-1;
	public boolean Strand;
	public byte function=-1;
	public byte TranscriptsAffected=-1;
	public byte TotalTranscripts=-1;
	
	
	public Annotation_g1k(){}
	
	public Annotation_g1k(String input_string, HashMap<String, Integer> gene_name2index){
		HashMap<String, Byte> map_name2index = ConstantsV4Annotations.map_name2index;
		HashMap<String, Byte> actg2index = ConstantsV4Annotations.actg2index;
		HashMap<String, Byte> noncoding_function2index = ConstantsV4Annotations.noncoding_function2index;
		HashMap<String, Byte> snp_fucntion2index = ConstantsV4Annotations.snp_fucntion2index;
		HashMap<String, Byte> indel_fucntion2index = ConstantsV4Annotations.indel_fucntion2index;
		HashMap<String, Byte> vt2index = ConstantsV4Annotations.vt2index;
		
		String[] temp=input_string.split("\t");
		
		if(actg2index.containsKey(temp[map_name2index.get("REF")])){
			this.REF=actg2index.get(temp[map_name2index.get("REF")]);
		}else {	this.REF=actg2index.get("*");}
		
		if(actg2index.containsKey(temp[map_name2index.get("ALT")])){
			this.ALT=actg2index.get(temp[map_name2index.get("ALT")]);
		}else {	this.ALT=actg2index.get("*");}		
		
		try{
			this.QUAL=Integer.parseInt(temp[map_name2index.get("QUAL")]);
		}catch(Exception e){this.QUAL=-1;}
		
		if(actg2index.containsKey(temp[map_name2index.get("AA")])){
			this.AA=actg2index.get(temp[map_name2index.get("AA")]);
		}else {	this.AA=actg2index.get("*");}	
		
		try{
			this.AC=Integer.parseInt(temp[map_name2index.get("AC")]);
		}catch(Exception e){this.AC=-1;}
		
		try{
			this.AN=Integer.parseInt(temp[map_name2index.get("AN")]);
		}catch(Exception e){this.AN=-1;}
		
		try{
			this.AF=Double.parseDouble(temp[map_name2index.get("AF")]);
		}catch(Exception e){this.AF=Double.NaN;}
		
		try{
			this.AFR_AF=Double.parseDouble(temp[map_name2index.get("AFR_AF")]);
		}catch(Exception e){this.AFR_AF=Double.NaN;}
		try{
			this.AMR_AF=Double.parseDouble(temp[map_name2index.get("AMR_AF")]);
		}catch(Exception e){this.AMR_AF=Double.NaN;}
		try{
			this.ASN_AF=Double.parseDouble(temp[map_name2index.get("ASN_AF")]);
		}catch(Exception e){this.ASN_AF=Double.NaN;}
		try{
			this.EUR_AF=Double.parseDouble(temp[map_name2index.get("EUR_AF")]);
		}catch(Exception e){this.EUR_AF=Double.NaN;}
	
		String nac=temp[map_name2index.get("NAC")].split(",")[0];
		if(noncoding_function2index.containsKey(nac)){
			this.NAC=noncoding_function2index.get(nac);
		}else {	this.NAC=noncoding_function2index.get("*");}
		
		if(vt2index.containsKey(temp[map_name2index.get("VT")])){
			this.VT=vt2index.get(temp[map_name2index.get("VT")]);
		}else {	this.VT=vt2index.get("*");}
		
		//process VA:
		String input=temp[map_name2index.get("VA")];
		String[] info=input.split(":");
		if(info.length<5){
//			System.out.println(input);
			return ;
		}		
		this.AlleleNumber=Byte.parseByte(info[0]);
		this.GeneName_index=gene_name2index.get(info[1]);
		if(info[3].equals("+"))this.Strand=true;
		else if(info[3].equals("-"))this.Strand=false;
		else{System.out.println("trand not + or -");}
		if(this.VT==0){
			this.function=snp_fucntion2index.get(info[4]);
		}else if(this.VT==1){
			this.function=indel_fucntion2index.get(info[4]);
		}else{
			this.function=-1;
		}String[] the_fraction=info[5].split("/");
		if(the_fraction.length!=2){
			System.out.println("the_fraction.length!=2: "+info[5]);
		}
		this.TranscriptsAffected=Byte.parseByte(the_fraction[0]);
		this.TotalTranscripts=Byte.parseByte(the_fraction[1]);
	}
	
	/*
	 * (non-Javadoc)
	 * @see java.lang.Object#toString()
	 * public byte REF;
	public byte ALT;
	public int QUAL;	
	public byte AA; // ancestral allele
	
	public int AC; //Alternate Allele Count
	public int AN; //Total Allele Count
	public double AF; // AC/AN
	
	public double AFR_AF;	
	public double AMR_AF;
	public double ASN_AF;
	public double EUR_AF;
	
	public byte VT;
	public byte NAC;

	public byte AlleleNumber;
	public int GeneName_index;
	public boolean Strand;
	public byte function;
	public byte TranscriptsAffected;
	public byte TotalTranscripts;
	 */
	
	public boolean functional(){
		if(this.function!=-1){
			if(ConstantsV4Annotations.var_type[this.VT].equals("SNP")){
				if(!ConstantsV4Annotations.SNPsFunctionsal[this.function].equals("synonymous")&&
					!ConstantsV4Annotations.SNPsFunctionsal[this.function].equals("*"))	
				return true;
			}else if(ConstantsV4Annotations.var_type[this.VT].equals("INDEL")||ConstantsV4Annotations.var_type[this.VT].equals("SV")){
				return true;
			}
		}
		if(!ConstantsV4Annotations.non_coding_functional[this.NAC].equals("*"))return true;
		return false;
	}
	
	public String toString(String[] gene_names){
		String function="*";
		if(this.function!=-1){
			if(this.VT==0)function=ConstantsV4Annotations.SNPsFunctionsal[this.function];
			else if(this.VT==1)function=ConstantsV4Annotations.IndelsFunctionsal[this.function];
		}
		String result= ConstantsV4Annotations.acgt[this.REF]+"\t"+
				ConstantsV4Annotations.acgt[this.ALT]+"\t"+
				this.QUAL+"\t"+
				ConstantsV4Annotations.acgt[this.AA]+"\t"+
				this.AC+"\t"+	this.AN+"\t"+
				this.AF+"\t"+	this.AFR_AF+"\t"+this.AMR_AF+"\t"+this.ASN_AF+"\t"+this.EUR_AF+"\t"+
				ConstantsV4Annotations.var_type[this.VT]+"\t"+
				ConstantsV4Annotations.non_coding_functional[this.NAC]+"\t";
		if(this.AlleleNumber>=1)
			result=result+this.AlleleNumber+":"+gene_names[this.GeneName_index]+":"+(this.Strand?"+":"-")+":"+function+":"+this.TranscriptsAffected+"/"+this.TotalTranscripts;
		else result=result+"*";
		return result;
	}
}
