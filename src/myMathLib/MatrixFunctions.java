package myMathLib;

import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.util.Arrays;

public class MatrixFunctions {
    double[] eigen_values;
    double[][] eigen_vectors;
    double[] single_values;
    double[][] left_single_values;
    double[][] right_single_values;

    public void eigen_decomp_symm(double[][] inp) throws IOException, InterruptedException {
        String curr_dir = System.getProperty("user.dir");

        BinaryFiles.write_bin_file(curr_dir + "/InputMatrix.bin", inp);

        ProcessBuilder pb = new ProcessBuilder("ocma", "eigen", "double", "disk",
                Integer.toString(inp.length), curr_dir + "/InputMatrix.bin", "eigenvals.bin", "eigenvecs.bin");

        Process pr = pb.start();

        System.out.println("Waiting for calculations...");

        pr.waitFor();

        eigen_values = BinaryFiles.read_bin_file(curr_dir + "/eigenvals.bin");
        eigen_vectors = BinaryFiles.read_bin_file(curr_dir + "/eigenvecs.bin", inp.length, inp.length);

        delete_file(curr_dir + "/InputMatrix.bin");
        delete_file(curr_dir + "/eigenvals.bin");
        delete_file(curr_dir + "/eigenvecs.bin");

    }

    public void svd(double[][] inp) throws IOException, InterruptedException {
        System.out.println("Starting timing...");

        long start_time = System.nanoTime();
        long full_start_time = start_time;

        String curr_dir = System.getProperty("user.dir");

        BinaryFiles.write_bin_file(curr_dir + "/InputMatrix.bin", inp);

        long end_time = System.nanoTime();

        System.out.println("Time taken to write input matrix: " + ((end_time - start_time) / 1000000000.0) + " seconds");

        start_time = System.nanoTime();

        ProcessBuilder pb = new ProcessBuilder("ocma", "singular", "double", "disk",
                Integer.toString(inp.length), Integer.toString(inp[0].length), curr_dir + "/InputMatrix.bin",
                "single_vals.bin", "left_single.bin", "right_single.bin");

        Process pr = pb.start();

        System.out.println("Waiting for calculations...");

        pr.waitFor();

        end_time = System.nanoTime();
        System.out.println("Time taken for actual calculations: " + ((end_time - start_time) / 1000000000.0) + " seconds");

        start_time = System.nanoTime();
        single_values = BinaryFiles.read_bin_file(curr_dir + "/single_vals.bin");
        left_single_values = BinaryFiles.read_bin_file(curr_dir + "/left_single.bin", inp.length, inp.length);
        right_single_values = BinaryFiles.read_bin_file(curr_dir + "/right_single.bin", inp[0].length, inp[0].length);
        end_time = System.nanoTime();
        System.out.println("Time taken to read output matrices: " + ((end_time - start_time) / 1000000000.0) + " seconds");

        start_time = System.nanoTime();
        delete_file(curr_dir + "/single_vals.bin");
        delete_file(curr_dir + "/left_single.bin");
        delete_file(curr_dir + "/right_single.bin");

        end_time = System.nanoTime();
        System.out.println("Time taken to delete files: " + ((end_time - start_time) / 1000000000.0) + " seconds");

        System.out.println("\nTotal time taken: " + ((end_time - full_start_time) / 1000000000.0) + " seconds");

    }

    private static void delete_file(String file_path) {
        File fl = new File(file_path);
        if (!fl.delete()) {
            System.out.println("Error deleting " + file_path);
        }
    }

    public static void print_matrix(double[][] inp) {
        for (double d[] : inp) {
            System.out.println(Arrays.toString(d));
        }
    }

    public static void print_matrix(double[] inp) {
        System.out.println(Arrays.toString(inp));
    }

    public static void save_matrix(double[][] inp, String fileName) throws IOException {
        RandomAccessFile stream = new RandomAccessFile(fileName, "rw");
        FileChannel channel = stream.getChannel();


        for (int i = 0; i < inp.length; i++) {
            StringBuilder builder = new StringBuilder();
            for (int j = 0; j < inp[i].length; j++)//for each column
            {
                builder.append(inp[i][j] + "");//append to the output string
                if (j < inp[j].length - 1)//if this is not the last row element
                    builder.append(" ");//then add comma (if you don't like commas you can use spaces)
            }
            builder.append("\n");
            byte[] strBytes = builder.toString().getBytes();
            ByteBuffer buffer = ByteBuffer.allocate(strBytes.length);
            buffer.put(strBytes);
            buffer.flip();
            channel.write(buffer);
        }

        stream.close();
        channel.close();
    }
}

