package mixedmodel;

import java.util.Arrays;

public class Regional_Results {
	public VariantsDouble genotype;
	public int chr;
	public int start_loc;
	public int end_loc;
	public int[] all_locs_in_region;
	public double[][] genotype_data_in_region;
	public int num_sites_in_region;
	
	public int[] selected_locs_indexes_in_region;
	public double[] variance_explained;
	public double[] pvalue;
	
	public Regional_Results(VariantsDouble genotype, int chr, int start_loc, int end_loc){
		this.genotype = genotype;
		this.chr = chr;
		this.start_loc = start_loc;
		this.end_loc = end_loc;

		int start_index = Arrays.binarySearch(this.genotype.locations[chr], start_loc);
		if (start_index < 0) {// no_found, and make use of the insertion-point returned by the binarySearch:
			start_index = -(start_index + 1);
		}
		int end_index = Arrays.binarySearch(this.genotype.locations[chr], end_loc);
		if (end_index < 0) {// no_found, and make use of the insertion-point returned by the binarySearch:
			end_index = -(end_index + 1) - 1;
		}
		this.genotype_data_in_region = this.genotype.reader.float64().readMatrixBlockWithOffset(genotype.position_fast_paths[chr],
				(end_index - start_index + 1), this.genotype.sample_size, start_index, 0);
		this.all_locs_in_region = new int[end_index - start_index + 1];
		if (this.genotype_data_in_region.length == this.all_locs_in_region.length) {
			this.num_sites_in_region = this.all_locs_in_region.length;
		} else {
			System.out.println("Incorrect DATA: this.genotype_data_in_region.length!=this.all_locs_in_region.length");
		}
		for (int i = start_index; i <= end_index; i++) {
			this.all_locs_in_region[i - start_index] = genotype.locations[chr][i];
		}

	}
}
