package myMathLib;

import java.io.*;
import java.nio.*;
import java.nio.channels.*;
import java.util.*;

import static java.lang.Math.ceil;

public class BinaryFiles {
    private static final int buff_size = 4096;

    public static double[] read_bin_file(String fileName) throws IOException {
        InputStream inp = new FileInputStream(fileName);
        File fl = new File(fileName);

        double[] data = new double[(int) ceil(fl.length() / 8.0)];
        int counter = 0;

        byte[] buff = new byte[buff_size];

        while (inp.read(buff) != -1) {
            for (int i = 0; counter < data.length && i < buff.length; i = i + 8) {
                data[data.length - counter - 1] = (ByteBuffer.wrap(Arrays.copyOfRange(buff, i, i + 8)).order(ByteOrder.LITTLE_ENDIAN).getDouble());
                counter++;
            }
        }

        return data;
    }

    public static double[][] read_bin_file(String fileName, int x_size, int y_size) throws IOException {
        InputStream inp = new FileInputStream(fileName);
        File fl = new File(fileName);
        long file_size = fl.length();
        int total_size = x_size * y_size;

        if ((file_size / 8) > (total_size)) {
            System.out.println("File contains more digits than specified matrix size! Halting process...");
            System.exit(4444);
        }

        double[][] data = new double[x_size][y_size];
        int x_counter = 0;
        int y_counter = 0;
        int counter = 0;

        byte[] buff = new byte[buff_size];

        while (inp.read(buff) != -1) {
            for (int i = 0; counter < total_size && i < buff.length; i = i + 8) {
                data[x_size - x_counter - 1][y_counter] = (ByteBuffer.wrap(Arrays.copyOfRange(buff, i, i + 8)).order(ByteOrder.LITTLE_ENDIAN).getDouble());
                y_counter++;
                counter++;

                if (y_counter >= y_size) {
                    x_counter++;
                    y_counter = 0;
                }
            }
        }

        return data;
    }

//


    public static void write_bin_file(String fileName, double[] data) throws IOException {
        RandomAccessFile fil = new RandomAccessFile(fileName, "rw");
        FileChannel inChannel = fil.getChannel();

        ByteBuffer buf = ByteBuffer.allocate(buff_size).order(ByteOrder.LITTLE_ENDIAN);
        int counter = 0;

        while (counter < data.length) {
            buf.clear();
            while (buf.position() < buf.capacity()) {
                if (counter >= data.length) break;
                buf.putDouble(data[counter]);
                counter++;
            }

            buf.flip();
            while (buf.hasRemaining()) {
                inChannel.write(buf);
            }
        }
        inChannel.close();
    }

    public static void write_bin_file(String fileName, double[][] data) throws IOException {
        RandomAccessFile fil = new RandomAccessFile(fileName, "rw");
        FileChannel inChannel = fil.getChannel();

        ByteBuffer buf = ByteBuffer.allocate(buff_size).order(ByteOrder.LITTLE_ENDIAN);

        int x_counter = 0;
        int y_counter = 0;
        int counter = 0;

        int x_size = data.length;
        int y_size = data[0].length;
        int total_size = x_size * y_size;

        while (counter < total_size) {
            buf.clear();
            while (buf.position() < buf.capacity()) {
                if (counter >= total_size) break;
                buf.putDouble(data[x_counter][y_counter]);
                counter++;
                y_counter++;

                if (y_counter >= y_size) {
                    x_counter++;
                    y_counter = 0;
                }
            }

            buf.flip();
            while (buf.hasRemaining()) {
                inChannel.write(buf);
            }
        }
        inChannel.close();
    }
}

