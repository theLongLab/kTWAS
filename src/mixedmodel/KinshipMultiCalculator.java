package mixedmodel;

import ch.systemsx.cisd.base.mdarray.MDDoubleArray;
import ch.systemsx.cisd.hdf5.HDF5MDDataBlock;

public class KinshipMultiCalculator implements Runnable {
    Iterable<HDF5MDDataBlock<MDDoubleArray>> position_fast_blocks;
    double scale;
    double min_MAF;
    int chr;
    int sample_size;
    KinshipMatrixMultiThreaded k;

    KinshipMultiCalculator(int chr, Iterable<HDF5MDDataBlock<MDDoubleArray>> position_fast_blocks, double scale, double min_MAF, int sample_size, KinshipMatrixMultiThreaded k) {
        this.chr = chr;
        this.position_fast_blocks = position_fast_blocks;
        this.min_MAF = min_MAF;
        this.scale = scale;
        this.sample_size = sample_size;
        this.k = k;
    }

    @Override
    public void run() {
        System.out.println("Processing Chromosome " + this.chr + "...");
        for (HDF5MDDataBlock<MDDoubleArray> block : this.position_fast_blocks) {
            double start_time = System.nanoTime();
            double[][] data4thisblock = block.getData().toMatrix();
            for (int var_index = 0; var_index < data4thisblock.length; var_index++) {
                if (KinshipMatrix.maf(data4thisblock[var_index], scale) < min_MAF) continue;
                for (int i = 0; i < sample_size; i++) {
                    for (int j = i + 1; j < sample_size; j++) {
                        if ((!Double.isNaN(data4thisblock[var_index][i])) && (!Double.isNaN(data4thisblock[var_index][j]))) {
                            k.add_to_kinship(i, j, (scale - Math.abs(data4thisblock[var_index][i] - data4thisblock[var_index][j])));
                            k.increment_tot_num_var_used(i, j);
                        }
                    }
                }
            }
            double end_time = System.nanoTime();
            System.out.println("Time taken for 1 block in chr " + this.chr + ": " + ((end_time - start_time) / 1000000000));
        }
        System.out.println("Completed Chromosome " + this.chr + "...");
    }
}
