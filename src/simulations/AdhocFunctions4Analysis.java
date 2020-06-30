package simulations;

import mixedmodel.LocalKinshipAnalyzer;
import mixedmodel.MultiPhenotype;

public class AdhocFunctions4Analysis {

	public static void analyze_local(String input_geno, String input_pheno, String output_folder, 
		String local_kinship_files_folder, String global_kinship_file, int win_size){
	int the_phe_index=0;
	double step=0.01;
	LocalKinshipAnalyzer local_k=new LocalKinshipAnalyzer(input_geno, win_size, null);
	MultiPhenotype phenotypeS=new MultiPhenotype(input_pheno);				
	String out_phe_file=output_folder+"Local_VO."+the_phe_index+"."+phenotypeS.phenotypes[the_phe_index].phe_id+".w"+win_size+".csv";
	local_k.local_VO_all_grids(phenotypeS.phenotypes[the_phe_index], input_geno, global_kinship_file, local_kinship_files_folder,
						out_phe_file, step, false);	
	
	}

}
