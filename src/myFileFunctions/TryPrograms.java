package myFileFunctions;

public class TryPrograms {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		String inFilename0="/Users/quan.long/temp.txt";
		String inFilename="/Users/quan.long/temp.txt.gz";
		String outFilename="/Users/quan.long/some.txt";
//		FileFunc.gzip(inFilename0, inFilename);
		FileFunc.gunzip(inFilename, outFilename);
		

	}

}
