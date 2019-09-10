package ijaux.scale;

//import ij.IJ;

/*
 * 	Version 1.1  8 Jun 2019
 *    		1.0  21 March 2014
 *  
 *  library utility
 *  
 */
public class SUtils {

	private SUtils() {
		// no initialization
	}
	
	
	/**
	 * 
	 * @param kernel
	 * @param sz
	 */
	public static void transp(double[] kernel, int sz) {
		if (kernel.length/sz!=sz) 
			return; // not a square kernel		

		for (int row=0; row< sz; row++) {
			for (int col=0; col< sz; col++) {
				final int i=row*sz+col;
				final double c=kernel[i];
				final int j=row+col*sz;
				kernel[i]=kernel[j];
				kernel[j]=c;
			}
		}
	}
	
	/**
	 * 
	 * @param kernel
	 * @param sz
	 */
	public static void transp(float[] kernel, int sz) {
		if (kernel.length/sz!=sz) 
			return; // not a square kernel		

		for (int row=0; row< sz; row++) {
			for (int col=0; col< sz; col++) {
				final int i=row*sz+col;
				final float c=kernel[i];
				final int j=row+col*sz;
				kernel[i]=kernel[j];
				kernel[j]=c;
			}
		}
	}
	
	
	/**
	 * @param kernel
	 */
	public static void flip(float[] kernel) {
		final int s=kernel.length-1;
		for (int i=0; i< kernel.length/2; i++) {
			final float c=kernel[i];
			kernel[i]=kernel[s-i];
			kernel[s-i]=c;
		}
	}
	
	/**
	 * @param kernel
	 */
	public static void flip(double[] kernel) {
		final int s=kernel.length-1;
		for (int i=0; i< kernel.length/2; i++) {
			final double c=kernel[i];
			kernel[i]=kernel[s-i];
			kernel[s-i]=c;
		}
	}
	
	
	/**
	 *  prints binomial coefficients
	 *
	public static void printBinCoefs(double [][] bincoefs) {
		for (int i=0; i<bincoefs.length; i++) {
			for (int j=0; j<bincoefs[i].length; j++) {
				IJ.log("C^ ["+i + "] _["+j+"] =" +bincoefs[i][j]);
			}
		}
	}
	*/
	/**
	 *  implements Matlab\Octave function linear space
	 * @param a
	 * @param b
	 * @param N
	 * @return
	 */
	public static double[] linspace(double a, double b, int n) {
		final double[] ret = new double[n];
		final double step =  (b-a)/(n-1);
		for(int i=0;i<n;i++){
			ret[i] = i*step +a;
		}
		return ret;
	}
	
	/**
	 *  implements Matlab\Octave function linear space
	 * @param a
	 * @param b
	 * @param N
	 * @return
	 */
	public static float[] linspace(float a, float b, int n) {
		final float[] ret = new float[n];
		final double step =  (b-a)/(n-1);
		// minimal loss of precision this way
		for(int i=0;i<n;i++){
			ret[i] =(float) (i*step +a);
		}
		return ret;
	}

	
	/** calculates the value of a polynomial assuming a_n + \sum a_{n-i} x^i
	 *  by Horners' method
	 * @param coef
	 * @param x
	 * @return
	 */
	public static double polyval(double[] coef, double x) {		
		final int n=coef.length-1;
		double ret=coef[0];
		for (int i=1; i<=n; i++) {			
			final double z=ret*x + coef[i];
			ret= z;	
		}	
		return ret;
	}
	
	/** calculates the value of a polynomial assuming a_0 + \sum a_i x^i
	 *  by Horners' method
	 * @param coef
	 * @param x
	 * @return value at x
	 */
	public static double polyval2(double[] coef, double x) {			
		final int n=coef.length-1;
		double ret=coef[n];
		for (int i=1; i<=n; i++) {			
			final double z=ret*x + coef[n-i];
			ret= z;			
		}
		return ret;
	}
	
	/**
	 * implements Matlab\Octave function cumulative sum cumsum 
	 * @param seq
	 * @return
	 */
	public static double[] cumsum(double[] seq) {
		double ss=0;
		double[] cs=new double[seq.length];
		for (int i=0; i<seq.length; i++) {
			ss+=seq[i];
			cs[i]=ss;
		}
		return cs;
	}

	/**
	 * Tensor product of two arrays
	 * 
	 * @param kernel
	 * @param a
	 * @param b
	 * @return
	 */
	public static float[] joinXY(float[][] kernel, int a, int b) {

		int wa=kernel[a].length; // cols
		int wb=kernel[b].length; // rows
	
		float[] jkernel=new float[wa*wb];

		for (int i=0; i<jkernel.length; i++) {
			jkernel[i]=1.0f;
		}
		if (a>=0) { // columns
			final float[] col=kernel[a]; // ->wa
			for (int c=0; c<wa; c++) { // col
				for (int r=0; r<wb; r++) { // row							
					final int idx=c + r*wa;
					jkernel[idx]*=col[c];
				}
			}
		}
		if (b>=0) { // rows
			final float[] row=kernel[b]; // ->wb
			for (int r=0; r<wb; r++) { // row	
				for (int c=0; c<wa; c++) { // col					
					final int idx=c + r *wa;
					jkernel[idx]*=row[r];
				}
			}
		}
		return jkernel;

	}
}
