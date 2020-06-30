package nam;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.HashSet;

import myMathLib.StatFuncs;

public class Imputation {
	
//	only write RIL sample's phenotype and genotype, do not include founders
	public static void imputefull(Founder found, RILs[] ril, String pheno, String geno){
		try{
			BufferedWriter bwphe = new BufferedWriter(new FileWriter(pheno));
			bwphe.write("ID");
			for(int i=0; i<ril[0].pheno_name.length; i++){
				bwphe.write("\t"+ril[0].pheno_name[i]);
			}bwphe.write("\n");
						
			for(int f=0; f<ril.length; f++){
				for(int id =0; id<ril[f].phenotypes.length; id++){
					bwphe.write(ril[f].rilID[id]);
					for(int p=0; p<ril[f].pheno_name.length; p++){
						bwphe.write("\t"+ril[f].phenotypes[id][p]);
					}bwphe.write("\n");
				}
			}
			bwphe.flush(); bwphe.close();
			
			BufferedWriter bwgen = new BufferedWriter(new FileWriter(geno));
			bwgen.write("Chromosome,Positions");
			for(int f=0; f<ril.length; f++){
				for(int id=0; id<ril[f].rilID.length; id++){
					bwgen.write(","+ril[f].rilID[id]);
				}
			}bwgen.write("\n");
			
			for(int c =0; c<found.snps.length; c++){			
				for(int m=0; m<found.snps[c].length; m++){
					bwgen.write((c+1)+","+found.snps[c][m]);				
					for(int f=0; f<ril.length; f++){
						int left=0; int right=0; 
						int m4chr =ril[f].marker_info.locations[c].length;
						for(int i=0; i<m4chr; i++){
							if(found.snps[c][m]>ril[f].marker_info.locations[c][i]){
								left =i;
							}
							if(found.snps[c][m]<=ril[f].marker_info.locations[c][i]){
								right=i; break;
							}
						}
						
//	start impute
						int fA_index =found.sample_id.indexOf(ril[f].parents[0]);
						int fB_index =found.sample_id.indexOf(ril[f].parents[1]);
						char A =found.genotype[fA_index][c][m];
						char B =found.genotype[fB_index][c][m];
						if(right==0){
							for(int id=0; id<ril[f].genotype.length; id++){
//								System.out.println(ril[f].rilID[id]+"\t"+ril[f].genotype[id][c][0]);
								bwgen.write(","+imputeRIL(ril[f].genotype[id][c][0], A, B));
							}
						}else if(left==(m4chr-1)){
							for(int id=0; id<ril[f].genotype.length; id++){
								bwgen.write(","+imputeRIL(ril[f].genotype[id][c][left], A, B));
							}
						}else{
							if((right-left)!=1){
								System.out.println("right and left index wrong ");
							}else{
								for(int id=0; id<ril[f].genotype.length; id++){
									if(ril[f].genotype[id][c][left]==ril[f].genotype[id][c][right]){
										bwgen.write(","+imputeRIL(ril[f].genotype[id][c][left], A, B));
									}else{								
										double lef_loc =ril[f].marker_info.locations[c][left];
										double rig_loc =ril[f].marker_info.locations[c][right];		
//										System.out.println("Recomebination rate "+(found.snps[c][m]-lef_loc)/(rig_loc-lef_loc));
										if(ril[f].genotype[id][c][left]=='H' || ril[f].genotype[id][c][right]=='H'){
											double r=myMathLib.Test.randomNumber();	
											
											if((found.snps[c][m]-lef_loc)/(rig_loc-lef_loc)>=r){
												bwgen.write(","+imputeRIL(ril[f].genotype[id][c][right], A, B));
											}else{
												bwgen.write(","+imputeRIL(ril[f].genotype[id][c][left], A, B));
											}
										}else if(ril[f].genotype[id][c][left]!='H' && ril[f].genotype[id][c][right]!='H'){
											double r1=myMathLib.Test.randomNumber();	
											double r2=myMathLib.Test.randomNumber();	
											if((found.snps[c][m]-lef_loc)/(rig_loc-lef_loc)<r1){
												bwgen.write(","+imputeRIL(ril[f].genotype[id][c][left], A, B));
											}else if((found.snps[c][m]-lef_loc)/(rig_loc-lef_loc)<r2){
												bwgen.write(","+imputeRIL('H', A, B));
											}else{
												bwgen.write(","+imputeRIL(ril[f].genotype[id][c][right], A, B));
											}
										}
										
									}
								}
							}
						}
					}
					bwgen.write("\n");bwgen.flush();
				}
			}bwgen.close();
		}catch(Exception e){e.printStackTrace();}	
	}
	public static char imputeRIL(char RILgeno, char founderA, char founderB){
		char result ='N';
		if(RILgeno=='A'){
			result =founderA;
		}else if(RILgeno=='B'){
			result =founderB;
		}else if(RILgeno=='H' || RILgeno=='C'){
			result =two2one(founderA, founderB);
		}else{System.out.println("unexpected "+RILgeno);}
		return result;
	}
	
	 public static char two2one(char snp1, char snp2){
		   char result = '*';
		   if(snp1==snp2){
			   result =snp1;
		   }else{
			   if(snp1=='A' || snp2=='A'){
				   if(snp1=='C' || snp2=='C'){
					   result ='M';
				   }else if(snp1=='G' || snp2=='G'){
					   result ='R';
				   }else if(snp1=='T' || snp2=='T'){
					   result ='W';
				   }
			   }else if(snp1=='C' || snp2=='C'){
				   if(snp1=='G' || snp2=='G'){
					   result ='S';
				   }else if(snp1=='T' || snp2=='T'){
					   result ='Y';
				   }
			   }else if(snp1=='G' || snp2=='G'){
				   if(snp1=='T' || snp2=='T'){
					   result ='K';
				   }
			   }
		   }
		   return result;
	   }

}
