package myMathLib;

import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.io.*;
import java.util.ArrayList;

//import flanagan.analysis.Stat;
//import flanagan.physprop.Diffusion;

public class Test {
	
	public static int[] iv=new int[33];
    public static int iy;
    public static int idum = -235654654;
    
    /* Jurg's original documentations:
     * 
	 * Random number generator with period of 10^8. Highly recommended by
	Press et al. (1992) "Numerical recipes" Cambridge University Press, New York

	Usage:
	Must be initialized with a negative integer (idum). At termination, 
	may write negative of last idum, that is, -idum, to file containing seed.

	Main program must have the following declarations:

	type
	integer=longint;  4-byte integer
	real=double;

	var
	iv:array[1..32] of integer;
	iy:integer;

	Thus, the global variables iv and iy are used.
	 * 
	 */
	
	public static void setSeed(int x){
		idum = x;
	}
	
	public static double randomNumber()
	{
	int ia=16807, im=2147483647;
	
	double am=1.0/im;
	int iq=127773, ir=2836, ntab=32;
    int ndiv = 1 + ((im-1)/(ntab));
	double eps=1.2e-7;
	double rnmx=1.0-eps;

	int j,k;
	double rr;

	if (idum<=0) {
		idum=-idum;
		if (idum<1) 
			idum=1;
		for (j=1;j<=ntab;j++)
			iv[j]=0;
		iy=0;
		for (j=ntab+8;j>=1;j--) {
			k = (idum/iq);
			idum=ia*(idum-k*iq)-ir*k;
			if (idum<0) idum=idum+im;
			if (j<=ntab) iv[j]=idum;
		}
		iy=iv[1];
	}
	k = (idum/iq);
	idum=ia*(idum-k*iq)-ir*k;
	if (idum<0) idum=idum+im;
	j = ((1 + iy)/ndiv);
	iy=iv[j];
	iv[j]=idum;
	rr=am*iy;
	double rrran1;
	if (rr<rnmx) rrran1=rr; 
	else rrran1=rnmx;
	
	return rrran1;
	}  	
	
  public static double pnorm(double xori){
//  {
//   Copyright (C) Jurg Ott 1988-2006.  Partly programmed by Joe
//   Terwilliger using continued fractions.  Based on several formulas in
//   Abramowitz and Stegun, "Handbook of mathematical functions", Dover, May 1968.
//
//   Calculates the upper tail probability of the normal
//   distribution for a given normal deviate, x.
//
//    4 Jan 2001  Last program version
//   16 Aug 2002  Turned into function to replace old pnorm
//  }


   double  rp,x,q,logQ,ln10;
   
  		ln10=Math.log(10);
  		rp=0.9189385332046727417803296;
  		x=Math.abs(xori);
  		if (Math.abs(x-0.5)<0.00001){
  			q=0.5;
  		}
  		else if (x>1){ 
  		    // { Formula 26.2.14 in Abramowitz and Stegun, 1968, page 932}
  			int numit=400;  		   
  		    double inter;  		   
  		    inter=x+numit+1;
  		    for (int i=1; i<= numit; i++) 
  		    	inter=(numit-i+1)/inter+x;
  		    logQ=-Math.log(inter)-0.5*Math.pow(x,2)-rp;
  		    if (x<20) q=Math.exp(logQ); else q=0;
  		    logQ =logQ/ln10;
  		}
  		else {
  		   //{Implements formula 26.2.12, page 932, in Abramowitz and Stegun}
  		    double ff,sum,count,x2,two;  		  
  		    two=2;
  		    ff=x;
  		    x2=Math.pow(x,2);
  		    count=1;
  		    sum = ff;
  		    do{
  		     count = count+two;
  		     ff = ff*(x2/count);
  		     sum = sum+ff;
  		    }while(!(ff<=1.0E-300)); //{until it underflows}
  		    q=0.5-sum*Math.exp(-0.5*x2-rp);
  		    logQ=Math.log(q)/ln10;
  		}
  		if (xori<0) 
  			q=1-q;
  		double pnorm=q;
  		return pnorm;
  }
	
  public static double chi2pr(double x2, int ndf){
	
/*	From Chi-squire value to probability
 * 
 * {INCLUDE file for Linkage Utility programs.
*	 Computes the error probability, chipr = 1-f, where f is the
*	 chi-square distribution function with ndf degrees of freedom
*	 evaluated at x2.       Function required: pnorm.
*	 Follows formulas 26.4.4 and 26.4.21 in Abramowitz and Stegun,
*	 "Handbook of mathematical functions", Dover, New York 1968}
*/
	  double z,p,x,sum,re,ch,chp, chipr;
	 
	 int  n1, n2;	
	 if ((x2>0.0) && (ndf>0)){	 
		 if (ndf==1){	  
			 x=Math.sqrt(x2);
			 chipr=2.0*pnorm(x);
		 }
		 else{  
			 if (ndf==2)
			 { chipr=Math.exp(-0.5*x2);} //{ndf=2}
			 else{               //{ndf>2}
				 n1=(ndf-1)/ 2;
				 n2=(ndf-2)/ 2;
				 if (n1==n2){
	     // {ndf is even and >2}
					 sum=1.0;
					 re=0.5*x2;
					 ch=1.0;
					 for (int i=1; i<= n2;i++){
						 ch=ch*re/i;
						 sum=sum+ch;
					 }
					 chipr=Math.exp(-re)*sum;
				 }
				 else{             //{ndf is odd and >1}
					 ch=Math.sqrt(x2);
					 z=0.39894228*Math.exp(-0.5*x2);
					 p=pnorm(ch);
					 if (ndf==3)
	    	 			{chipr=2.0*(p+z*ch);} //{ndf=3}
					 else{             // {ndf odd and >3}
						 chp=ch;
						 re=1.0;
						 for (int i=2; i<= n1;i++){
							 re=re+2.0;
							 chp=chp*x2/re;
							 ch=ch+chp;
						 }
						 chipr=2.0*(p+z*ch);
					 }
				 }
			 }
		 }
	} 
	else	chipr =1.0;
	return chipr;
  }

  	public static double chiSquarePValue(double[][] table2){
  		double table[][] = preCheck4ZeroRowAndColumn(table2);
  		double value = chiSquareValue (table);
  		int df = (table.length-1)*(table[0].length-1);
  		return chi2pr(value, df);
    }
  	
  	public static double chiSquarePValueLR(double[][] table){
  //		double table[][] = preCheck4ZeroRowAndColumn(table2);
  		double value = chiSquareValueLR(table);
  		int df = (table.length-1)*(table[0].length-1);
  		return chi2pr(value, df);
    }
  	
  	public static double chiSquareValueLR(double[][] table){
    	
	//	double table[][] = preCheck4ZeroRowAndColumn(table2);
		
		int length = table.length;
		//System.out.println("length:"+length);
		int width =  table[0].length;
		double sumtablerow[] = new double[length];
		double sumtablecolum[] = new double[width];
		double sum=0;
		
		for (int i=0;i<length;i++)
			for (int j=0;j<width;j++){
				if (table[i][j]==0)
					table[i][j]= 0.5;
				sumtablerow[i]=sumtablerow[i]+ table[i][j];
				sumtablecolum[j]=sumtablecolum[j]+ table[i][j];
				}	
		for(int i=0;i<length;i++)
			sum = sum+sumtablerow[i];
		double chi2value=0;	
		for (int i=0;i<length;i++){
			for (int j=0;j<width;j++){
				chi2value += 
					table[i][j]*(Math.log(table[i][j])-Math.log(sumtablerow[i])-
							Math.log(sumtablecolum[j]));
			}			
		}
		chi2value+=sum*Math.log(sum);
		return 2*chi2value;
	}
  
    public static double chiSquareValue (double[][] table2){
    	
		double table[][] = preCheck4ZeroRowAndColumn(table2);
		
		double chi2value=0;
		int length = table.length;
		int width =  table[0].length;
		double sumtablerow[] = new double[length];
		double sumtablecolum[] = new double[width];
		double sumfrequencycolum[] = new double[width];
		double sumfrequencyrow[] = new double[length];
		double sum=0;
		for (int i=0;i<length;i++){
			for (int j=0;j<width;j++){
				sum=sum + table[i][j];
			}
		}
//		System.out.println(sum);

		for (int i=0;i<length;i++){
			for (int j=0;j<width;j++){
				sumtablerow[i]=sumtablerow[i]+ table[i][j];
			}
			sumfrequencyrow[i]=sumtablerow[i]/sum;
		}
//	outputArray(sumtablerow);
//	outputArray(sumfrequencyrow);
 	
		for (int j=0;j<width;j++){
			for (int i=0;i<length;i++){
				sumtablecolum[j]=sumtablecolum[j]+ table[i][j];
			}		
			sumfrequencycolum[j]=sumtablecolum[j]/sum;
		}
//	outputArray(sumtablecolum);
// 	outputArray(sumfrequencycolum); 			

		chi2value=0;	
		for (int i=0;i<length;i++){
			for (int j=0;j<width;j++){
				chi2value += 
					(Math.pow((table[i][j]-sumfrequencyrow[i]*sumfrequencycolum[j]*sum),2.0))/(sumfrequencyrow[i]*sumfrequencycolum[j]*sum);
			}
		}		
		return chi2value;
	}
	
	
	/*
	 * Generate pvalue by comparing statistic with an array gain by permutations
	 */
	public static double computePvalueByArray(double[] statArray, double trueStat){
		int count = 0;
		for(int i = 0;i < statArray.length; i++ ){
			if (trueStat<=statArray[i]) count++;
		}
		double pValue = ((double)count)/((double)statArray.length);
		return pValue;
	}
	
	/*
	 * Generate pvalue by comparing pvalue with an array gain by permutations
	 */
	public static double computePvalueByPvalueArray(double[] pValuleArray, double truePValue){
		int count = 0;
		for(int i = 0;i < pValuleArray.length; i++ ){
			if (truePValue>=pValuleArray[i]) count++;
		}
		double pValue = ((double)count)/((double)pValuleArray.length);
		return pValue;
	}
	
	public static void randomPermute(double[] indicators){
//		java.util.Random r= new java.util.Random();
		int mm = indicators.length;
		for (int ii=0; ii< mm; ii++){
		  // {i1 will be a random element from ii .. mm}
//			double randomNumber = r.nextDouble();
			double randomNumber = Test.randomNumber();
			double yy = randomNumber*(mm-1-ii);
			int zz = (int)yy;
			int xx = (yy-zz>=0.5)?(zz+1):(zz);
		    int i1 = ii + (xx);
		   //{Now exchange element ii with element i1 in the array}
		    double temp = indicators[i1];
		    indicators[i1] = indicators[ii];
		    indicators[ii] = temp;
		}
	}
	
	public static void outputArray(double[][] a){
		int length = a.length;
		int width = a[0].length;
		
		for(int i=0;i<length;i++){
			System.out.print("\n");
			for(int j=0;j<a[i].length;j++){
				System.out.print(a[i][j]);
				System.out.print('\t');
			}
		}System.out.println();
	}
	
	public static void outputArray2int(double[][] a){
		int length = a.length;
		int width = a[0].length;
		
		for(int i=0;i<length;i++){
			System.out.print("\n");
			for(int j=0;j<width;j++){
				System.out.print((int)a[i][j]);
				System.out.print(' ');
			}
		}System.out.println(' ');
	}
	
	public static void outputArray(int[][] a){
		int length = a.length;
		int width = a[0].length;
		
		for(int i=0;i<length;i++){
			System.out.print("\n");
			for(int j=0;j<width;j++){
				System.out.print(a[i][j]);
				System.out.print(' ');
			}
		}System.out.println(' ');
	}
	
	public static void outputArray(boolean[][] a){
		int length = a.length;
		int width = a[0].length;
		
		for(int i=0;i<length;i++){
			System.out.print("\n");
			for(int j=0;j<width;j++){
				System.out.print(a[i][j]);
				System.out.print(' ');
			}
		}System.out.println(' ');
	}
	
	public static void outputArray(double[] a){
		int length = a.length;
	
		for(int i=0;i<length;i++){			
				System.out.print(a[i]);
				System.out.print('\t');
			}
		System.out.println(' ');
	}
	
	public static void outputArray4d(double[] a){
		int length = a.length;	
		for(int i=0;i<length;i++){
			int x=(int)(a[i]*10000);
			double y=((double)x)/10000.0;
			System.out.print(y);
			System.out.print(' ');
		}
		System.out.println(' ');
	}
	
	public static void outputArray(int[] a){
		int length = a.length;
	
		for(int i=0;i<length;i++){			
				System.out.print(a[i]);
				System.out.print(' ');
			}
		System.out.println(' ');
	}
	
	public static int[][] file2Array(String filePath){
		int length = 0;
		int width = 0;
		File file = new File(filePath);
		char[][] tempArray;
		try{
			FileReader fr0 = new FileReader(file);
			BufferedReader br0 = new BufferedReader(fr0);
			String line0 = br0.readLine();
			char[] cut0 = line0.toCharArray();
			width = cut0.length;
		
			while(line0!=null){
				length++;
				line0 = br0.readLine();
			}
		}catch(Exception e) {
				  e.printStackTrace();
		}
		tempArray = new char[length][width];
		try{			
			FileReader fr = new FileReader(file);
			BufferedReader br = new BufferedReader(fr);
			String line = br.readLine();
			int i=-1;
			while (line!=null){
				i++;
				tempArray[i]=line.toCharArray();
				line = br.readLine();
			}			
		}catch(Exception e) {
			  e.printStackTrace();
		}
		int[][] array = char2Int(tempArray);
		return array;
		
	}
	
	public static double[][] char2Double(char[][] ori){ 
		int length = ori.length;
		int width = ori[0].length;
		double[][] array = new double[length][width];
		for(int i=0;i<length;i++){
			for(int j=0;j<width;j++){
				if (ori[i][j] == '0')	array[i][j]=0;
				else if (ori[i][j] == '1')	array[i][j]=1;
				else if (ori[i][j] == '2')	array[i][j]=2;
				else if (ori[i][j] == '3')	array[i][j]=3;
				else if (ori[i][j] == '4')	array[i][j]=4;
				else if (ori[i][j] == '5')	array[i][j]=5;
				else if (ori[i][j] == '6')	array[i][j]=6;
				else if (ori[i][j] == '7')	array[i][j]=7;
				else if (ori[i][j] == '8')	array[i][j]=8;
				else if (ori[i][j] == '9')	array[i][j]=9;
			}
		}		
		return array;
	}
	
	public static int[][] char2Int(char[][] ori){ 
		int length = ori.length;
		int width = ori[0].length;
		int[][] array = new int[length][width];
		for(int i=0;i<length;i++){
			for(int j=0;j<width;j++){
				if (ori[i][j] == '0')	array[i][j]=0;
				else if (ori[i][j] == '1')	array[i][j]=1;
				else if (ori[i][j] == '2')	array[i][j]=2;
				else if (ori[i][j] == '3')	array[i][j]=3;
				else if (ori[i][j] == '4')	array[i][j]=4;
				else if (ori[i][j] == '5')	array[i][j]=5;
				else if (ori[i][j] == '6')	array[i][j]=6;
				else if (ori[i][j] == '7')	array[i][j]=7;
				else if (ori[i][j] == '8')	array[i][j]=8;
				else if (ori[i][j] == '9')	array[i][j]=9;
			}
		}		
		return array;
	}
  
//	public static void init(double[] xx, int length){
//		xx= new double[length];
//	}
	
	public static double[][] preCheck4ZeroRowAndColumn (double[][] table){
    	int length = table.length;
		int width =  table[0].length;
		// check row
		int[] NonZeroRows = new int[length]; 
		int newLen = 0;
		for(int i=0;i<length;i++)
			for(int j=0;j<width;j++)
				if (table[i][j]!=0){ 
					NonZeroRows[i]=1;
					newLen++;
					break;
				}
		//outputArray(NonZeroRows);//debug
		double[][] newTable = new double[newLen][width];
		int index = -1;
		for(int i=0;i<length;i++)
			if(NonZeroRows[i]==1){
				index++;
				for(int j=0;j<width;j++)
					newTable[index][j]=table[i][j];
			}
		//outputArray(newTable);//debug
		// check column
		int[] NonZeroColumns = new int[width]; 
		int newWid = 0;
		for(int j=0;j<width;j++)
			for(int i=0;i<newLen;i++)
				if (newTable[i][j]!=0){ 
					NonZeroColumns[j]=1;
					newWid++;
					break;
				}
		//outputArray(NonZeroColumns);//debug
		double[][] newnewTable = new double[newLen][newWid];
		index = -1;
		for(int i=0;i<width;i++)
			if(NonZeroColumns[i]==1){
				index++;
				for(int j=0;j<newLen;j++)
					newnewTable[j][index]=newTable[j][i];
			}
		return newnewTable;
    }
	
 /*	program tprob(input,output);
	 {Copyright (C) Jurg Ott 2003-2006}
	 {Interactive program to calculate two-sided tail probabilities of the
	 t distribution with a given number of degrees of freedom (whole
	 numbers only, no fractional degrees of freedom). Method based
	 on "Continuous Univariate Distributions-2", N.L. Johnson & S. Kotz,
	 Houghton Mifflin, Boston, 1970, page 96, equations 4.2 - 4.4}
 */
	
 public static double  tpr(double th, int ndf){
	/* {Input is theta = arctan( t/sqrt(ndf) ). Computes
	  function A(t|ndf) }
    */
	 double tpr;
	 double term,numer,denom,sum,sint,cost,cost2;
	 int  i7,n1,n2;
	 double half=0.5, pi = 3.141592653589793238462643;
	 
     sint=Math.sin(th);
     cost=Math.cos(th);
	//  sincos(th,sint,cost); // computes sinus(theta) and cosinus(theta) 
	  cost2=cost*cost;
	  if (ndf==1) tpr=th/pi;
	  else {                //    {ndf>1}
	   if (ndf==2) tpr=half*sint; // {ndf=2}
	   else{                 // {ndf>2}
	    n1=(ndf-1)/2;
	    n2=(ndf-2) / 2;
	    if ( n1==n2){
	                       //  {ndf is even and >2}
	     sum=1.0;
	     numer=-1.0;
	     denom=0.0;
	     term=sum;
	     for( i7=1;i7<=n2;i7++){
	      numer=numer+2;
	      denom=denom+2;
	      term=term*numer*cost2/denom;
	      sum=sum+term;
	     }
	     tpr=half*sint*sum;
	    } //{then begin ndf even and >2}
	    else{               //   {ndf is odd and >1}
	     sum=cost;
	     if (ndf>3){
	      numer=0.0;
	      denom=1.0;
	      term=sum;
	      for (i7=2;i7<=n1;i7++){
	       numer=numer+2;
	       denom=denom+2;
	       term=term*numer*cost2/denom;
	       sum=sum+term;
	      }
	     }// {if ndf>3}
	     tpr=(th+sum*sint)/pi;
	    }// {else begin ndf odd and >1}
	   }// {else begin ndf>2}
	  }// {else begin ndf>1}
	  return tpr;
  }//{tpr}

	public static double tprob(double tt, int ndf){	
	  double theta,rr;
	  tt=Math.abs(tt);
	  rr=tt/Math.sqrt(ndf);
	  theta=Math.atan(rr);
	  rr=1.0-2*tpr(theta,ndf);
	  return rr;
	 }
	
	/*
	 * sort a double arrary and it's index
	 */			
	public static int[] sort(double[] freq){
		int len=freq.length;
		int[] index = new int[len];
		for(int i=0;i<len;i++)  index[i]=i;
		for(int i=0;i<len;i++){
			for(int j=i+1;j<len;j++){
				if(freq[i]<freq[j]){
					double temp = freq[i];
					freq[i]=freq[j];
					freq[j]=temp;
					int tempIndex=index[i];
					index[i]=index[j];
					index[j]=tempIndex;
				}
			}
		}
		return index;
	}
	
	/*
	 * sort an int arrary and it's corresponding string array
	 */			
	public static void sort(int[] freq, String[] index){
		int len=freq.length;
		if(index.length!=len) System.out.println("Error: int array length not equals to String array.");
		for(int i=0;i<len;i++){
			for(int j=i+1;j<len;j++){
				if(freq[i]<freq[j]){
					int temp = freq[i];
					freq[i]=freq[j];
					freq[j]=temp;
					String tempIndex=index[i];
					index[i]=index[j];
					index[j]=tempIndex;
				}
			}
		}
	}
	
	/*
	 * sort an double array and it's corresponding string array
	 */			
	public static void sort(double[] freq, String[] index){
		int len=freq.length;
		if(index.length!=len) System.out.println("Error: int array length not equals to String array.");
		for(int i=0;i<len;i++){
			for(int j=i+1;j<len;j++){
				if(freq[i]<freq[j]){
					double temp = freq[i];
					freq[i]=freq[j];
					freq[j]=temp;
					String tempIndex=index[i];
					index[i]=index[j];
					index[j]=tempIndex;
				}
			}
		}
	}
	
	/*
	 * sort an double array and it's corresponding char array
	 */			
	public static void sort(double[] freq,char[][] index){
		int len=freq.length;
		if(index.length!=len) System.out.println("Error: int array length not equals to String array.");
		for(int i=0;i<len;i++){
			for(int j=i+1;j<len;j++){
				if(freq[i]<freq[j]){
					double temp = freq[i];
					freq[i]=freq[j];
					freq[j]=temp;
					char[] tempIndex=index[i];
					index[i]=index[j];
					index[j]=tempIndex;
				}
			}
		}
	}
	
	public static HashMap<String, Integer> allele2index(){
		HashMap<String,Integer> allele2index=new HashMap<String, Integer>();
		allele2index.put("A", 0);allele2index.put("C", 1);
		allele2index.put("G", 2);allele2index.put("T", 3);
		return allele2index;
	}
	
	public static void write2file(double[] data, String file){
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(file));
			for(int i=0;i<data.length;i++){
				bw.write(data[i]+"\n");
			}bw.flush();
		}catch(Exception e){e.printStackTrace();}
	}
	
	public static void write2file(String[] data, String file){
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(file));
			for(int i=0;i<data.length;i++){
				bw.write(data[i]+"\n");
			}bw.flush();
		}catch(Exception e){e.printStackTrace();}
	}
	
	public static void write2file(String[][] data, String file){
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(file));
			for(int i=0;i<data.length;i++){
				for(int j=0;j<data[i].length;j++)
					bw.write(data[i][j]+" ");
				bw.write("\n");
				bw.flush();
			}
		}catch(Exception e){e.printStackTrace();}
	}
	
	public static void write2file(double[][] data, String file){
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(file));
			for(int i=0;i<data.length;i++){
				for(int j=0;j<data[i].length;j++)
					bw.write(data[i][j]+" ");
				bw.write("\n");
				bw.flush();
			}
		}catch(Exception e){e.printStackTrace();}
	}
	
	public static void write2file_convert(double[][] data, String file){
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(file));
			for(int j=0;j<data[0].length;j++){
				for(int i=0;i<data.length;i++)
					bw.write(data[i][j]+" ");
				bw.write("\n");
				bw.flush();
			}
		}catch(Exception e){e.printStackTrace();}
	}
	
	public static void write2file(char[][] data, String file){
		try{
			BufferedWriter bw=new BufferedWriter(new FileWriter(file));
			for(int i=0;i<data.length;i++){
				for(int j=0;j<data[i].length;j++)
					bw.write(data[i][j]+" ");
				bw.write("\n");
				bw.flush();
			}
		}catch(Exception e){e.printStackTrace();}
	}
	
	public static double[] load_data_1(String file){
		int len = 0;
		double[] data=new double[len];
		try{
			BufferedReader br = new BufferedReader(new FileReader(file));
			String line = br.readLine();
			while(line!=null){
				len++;
				line = br.readLine();
			}
			data = new double[len];
			br = new BufferedReader(new FileReader(file));
			line = br.readLine();
			int index=0;
			while(line!=null){
				data[index++]=Double.parseDouble(line);
				line = br.readLine();
			}
		}catch(Exception e){e.printStackTrace();}
		return data;
	}
	
	public static String[] load_data_s1(String file){
		int len = 0;
		String[] data=new String[len];
		try{
			BufferedReader br = new BufferedReader(new FileReader(file));
			String line = br.readLine();
			while(line!=null){
				len++;
				line = br.readLine();
			}
			data = new String[len];
			br = new BufferedReader(new FileReader(file));
			line = br.readLine();
			int index=0;
			while(line!=null){
				data[index++]=line;
				line = br.readLine();
			}
		}catch(Exception e){e.printStackTrace();}
		return data;
	}
	
	public static double[][] load_data_2(String file){
		int len = 0;
		double[][] data=new double[len][];
		try{
			BufferedReader br = new BufferedReader(new FileReader(file));
			String line = br.readLine();
			while(line!=null){
				len++;
				line = br.readLine();
			}
			data = new double[len][];
			br = new BufferedReader(new FileReader(file));
			line = br.readLine();
			int index=0;
			while(line!=null){
				String temp[]=line.split(" ");
				data[index]=new double[temp.length];
				for(int i=0;i<temp.length;i++){
					data[index][i]=Double.parseDouble(temp[i]);
				}
				index++;
				line = br.readLine();
			}
		}catch(Exception e){e.printStackTrace();}
		return data;
	}
	
	public static String[][] load_data_s2(String file){
		int len = 0;
		String[][] data=new String[len][];
		try{
			BufferedReader br = new BufferedReader(new FileReader(file));
			String line = br.readLine();
			while(line!=null){
				len++;
				line = br.readLine();
			}
			data = new String[len][];
			br = new BufferedReader(new FileReader(file));
			line = br.readLine();
			int index=0;
			while(line!=null){
				String temp[]=line.split(" ");
				data[index]=new String[temp.length];
				for(int i=0;i<temp.length;i++){
					data[index][i]=(temp[i]);
				}
				index++;
				line = br.readLine();
			}
		}catch(Exception e){e.printStackTrace();}
		return data;
	}
	
	public static void normalize(double[] a){
		double total=0;
		for(int i=0;i<a.length;i++)
			total=total+a[i];
		for(int i=0;i<a.length;i++)
			a[i]=a[i]/total;
	}
	
	public static void normalize_rm_minus(double[] a){
		for(int i=0;i<a.length;i++){
			if(a[i]<-1)	a[i]=Double.NaN;
			else if(a[i]<0) a[i]=0;
		}
		double total=0;
		for(int i=0;i<a.length;i++)
			total=total+a[i];
		for(int i=0;i<a.length;i++)
			a[i]=a[i]/total;
	}
	
	public static double fisher_method_combineP(double[] pvalues){
		if(pvalues.length==0)
			return 1;
		double chi=0;
		for(int i=0;i<pvalues.length;i++){
			chi=chi-Math.log(pvalues[i]);
		}
		return chi2pr(chi, pvalues.length*2);
	}
	
	public static double fisher_method_combineP(ArrayList<Double> pvalues){
		if(pvalues.size()==0)
			return 1;
		double chi=0;
		for(int i=0;i<pvalues.size();i++){
			chi=chi-Math.log(pvalues.get(i));
		}
		return chi2pr(chi, pvalues.size()*2);
	}
	
	/*
	 * test for variance difference:
	 * <math>F = \frac{(N-p)}{(p-1)} 
	 * \frac{\sum_{j=1}^{p} n_j (z_{\cdot j}-z_{\cdot\cdot})^2} 
	 * {\sum_{j=1}^{p}\sum_{i=1}^{n_j} (z_{ij}-z_{\cdot j})^2}</math>
	 * 
	 * return the stat of all groups, plus the p_value
	 */
//	public static double[] brown_forsythe_test(double[][] data){
//		int p=data.length;
//		int N=0;
//		for(int j=0;j<p;j++)N=N+data[j].length;
//		double[][] dif_data=data.clone();
//		for(int j=0;j<p;j++)dif_data[j]=data[j].clone();
//		double[] medians=new double[p];
//		for(int j=0;j<p;j++){
//			Arrays.sort(dif_data[j]);
//			if(dif_data[j].length%2==1)
//				medians[j]=dif_data[j][(dif_data[j].length-1)/2];
//			else
//				medians[j]=(dif_data[j][(dif_data[j].length)/2-1]+dif_data[j][(dif_data[j].length)/2])/2;
//			for(int i=0;i<dif_data[j].length;i++)
//				dif_data[j][i]=Math.abs(dif_data[j][i]-medians[j]);
//		}
//		double[] means=new double[p]; // this time it is mean of differences
//		double all_mean=0;
//		double[] results=new double[p+1];
//		for(int j=0;j<p;j++){
//			for(int i=0;i<dif_data[j].length;i++){
//				means[j]=means[j]+dif_data[j][i];
//			}
//			all_mean+=means[j];
//			means[j]=means[j]/dif_data[j].length;
//		}all_mean=all_mean/N;
//		double up=0, down=0;
//		for(int j=0;j<p;j++){
//			up=up+data[j].length*(means[j]-all_mean)*(means[j]-all_mean);
//			for(int i=0;i<dif_data[j].length;i++){
//				double sqr=(dif_data[j][i]-means[j])*(dif_data[j][i]-means[j]);
//				down=down+sqr;
//			}
//		}
//		results[p]=Stat.fCompCDF((N-p)*up/((p-1)*down),p-1,N-p);
////		System.out.println(f_value+","+(p-1)+","+(N-p));
//		return results;
//	}
	
	/*
	 * simple test against ratio of variance:
	 */
//	public static double[] variance_ratio_test(double[][] data){
//		double[] result=new double[data.length+1];
//		int p=data.length;
//		int N=0;
//		for(int j=0;j<p;j++)N=N+data[j].length;
//		double[][] dif_data=data.clone();
//		for(int j=0;j<p;j++)dif_data[j]=data[j].clone();
//		double[] medians=new double[p];
//		for(int j=0;j<p;j++){
//			Arrays.sort(dif_data[j]);
//			if(dif_data[j].length%2==1)
//				medians[j]=dif_data[j][(dif_data[j].length-1)/2];
//			else
//				medians[j]=(dif_data[j][(dif_data[j].length)/2-1]+dif_data[j][(dif_data[j].length)/2])/2;
//			for(int i=0;i<dif_data[j].length;i++)
//				dif_data[j][i]=Math.abs(dif_data[j][i]-medians[j]);
//		}
//		double[] vars=new double[p]; // this time it is mean of differences
//		for(int j=0;j<p;j++){
//			for(int i=0;i<dif_data[j].length;i++){
//				vars[j]=vars[j]+(dif_data[j][i]*dif_data[j][i]);
//			}
//			vars[j]=vars[j]/dif_data[j].length;
//		}
//		for(int j=0;j<p;j++)result[j]=vars[j];
//		result[p]=-1;
//		if(p!=2){
////			System.out.println("not two groups");
//			return result;
//		}
//		int bigger=1;
//		if(vars[0]>vars[1])bigger=0; 
//		double f_value=vars[bigger]/vars[1-bigger];
////		System.out.println(f_value);
//		result[p]=Stat.fCompCDF(f_value, dif_data[bigger].length-1, dif_data[1-bigger].length-1);
//		return result;
//	}
	
	public static double mean(double[] data){
		double result=0;
		for(int k=0;k<data.length;k++){
			result=result+data[k];
		}
		return result/data.length;
	}
	
	/*
	 * Sum( (data[i]-mean())2 ) / (size()-1)
	 */
	public static double variance(double[] data){
		double mean=mean(data);
		double result=0;
		for(int k=0;k<data.length;k++){
			result=result+(data[k]-mean)*(data[k]-mean);
		}
		return result/(data.length-1);		
	}
	
	public static double[][] binomial(int N, int K){
        double[][] binomial = new double[N+1][K+1];

        // base cases
        for (int k = 1; k <= K; k++) binomial[0][k] = 0;
        for (int n = 0; n <= N; n++) binomial[n][0] = 1;

        // bottom-up dynamic programming
        for (int n = 1; n <= N; n++)
            for (int k = 1; k <= K; k++)
                binomial[n][k] = binomial[n-1][k-1] + binomial[n-1][k];

        return binomial;
	}
}




