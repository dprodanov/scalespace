
import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.gui.Roi;
import ij.gui.DialogListener;
import ij.measure.Calibration;
import ij.plugin.filter.ExtendedPlugInFilter;
import ij.plugin.filter.PlugInFilterRunner;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
 
import ijaux.scale.*;

import java.awt.*;
import java.io.File;
import java.util.*;
 

 
import static java.lang.Math.*;

 
import dsp.Conv;
 

/**
 * @version 	30 Dec 2020

 * 				
 *   Z. P\"usp\"oki, M. Storath, D. Sage, M. Unser., 
 *   Transforms and Operators for Directional Bioimage Analysis: 
 *   A Survey, Advances in Anatomy, Embryology and Cell Biology, vol. 219, 
 *   Springer International Publishing, ch. 3, 2016.
 * 
 * @author Dimiter Prodanov, IMEC , Sumit Kumar Vohra , KULeuven
 *
 *
 * @contents
 * The plugin computes the eigenvalues of the Structure matrix
 * 
 * 
 * @license This library is free software; you can redistribute it and/or
 *      modify it under the terms of the GNU Lesser General Public
 *      License as published by the Free Software Foundation; either
 *      version 2.1 of the License, or (at your option) any later version.
 *
 *      This library is distributed in the hope that it will be useful,
 *      but WITHOUT ANY WARRANTY; without even the implied warranty of
 *      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *       Lesser General Public License for more details.
 *
 *      You should have received a copy of the GNU Lesser General Public
 *      License along with this library; if not, write to the Free Software
 *      Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

 
public class StructureT_ implements ExtendedPlugInFilter, DialogListener   {
    @SuppressWarnings("unused")

	private PlugInFilterRunner pfr=null;

	final int flags=DOES_ALL+KEEP_PREVIEW+ NO_CHANGES;
	private String version="1.0";
    @SuppressWarnings("unused")

	private int nPasses=1;
	private int pass;

	public final static String SIGMA="LOG_sigma",MAX_LEN="G_MAX",FULL_OUTPUT="Full_out",LEN="G_len";

	private boolean isEnabled=true;

	private float[][] kernel=null;

	private ImagePlus image=null;
	public static boolean debug=IJ.debugMode;

	public boolean fulloutput=false;

	private boolean isFloat=false;
    @SuppressWarnings("unused")

	private boolean hasRoi=false;

    // noise scale
	public int sz= Prefs.getInt(LEN, 2);
	
	// smoothing scale
	public int sz2= Prefs.getInt(LEN, 4);
 
	public int max_sz= Prefs.getInt(MAX_LEN, 8);

	/* NEW VARIABLES*/

	/** A string key identifying this factory. */
	private final  String FILTER_KEY = "STRUCTURE";

	/** The pretty name of the target detector. */
	private final String FILTER_NAME = "Strcuture components";
	
  	
	/** It stores the settings of the Filter. */
	private Map< String, String > settings= new HashMap<>();
	
	/** It is the result stack*/
	private ImageStack imageStack;

	/**
	 * This method is to setup the PlugInFilter using image stored in ImagePlus 
	 * and arguments of filter
	 */
	@Override
	public int setup(String arg, ImagePlus imp) {
		image=imp;
		isFloat= (image.getType()==ImagePlus.GRAY32);
		hasRoi=imp.getRoi()!=null;
		cal=image.getCalibration();
		return  flags;
	}

	final int Ox=0, Oy=1, Oz=2;

	// It is used to check whether to calibirate or not
	private boolean doCalib = false;
	/*
	 * This variable is to calibrate the Image Window
	 */
	private Calibration cal=null;
	
	/*
	public void initialseimageStack(ImageStack img){
		this.imageStack = img;
	}
	*/
	
	/*
	 * (non-Javadoc)
	 * @see ij.plugin.filter.PlugInFilter#run(ij.process.ImageProcessor)
	 */
	@Override
	public void run(ImageProcessor ip) {
		int r = (sz-1)/2;
		GScaleSpace sp=new GScaleSpace(r);

		int r2 = (sz2-1)/2;
		GScaleSpace sp2=new GScaleSpace(r2);
		imageStack=new ImageStack(ip.getWidth(),ip.getHeight());

		imageStack = filter2(ip,sp,sp2, imageStack);
		//imageStack = filter(ip,sp,imageStack);

		image=new ImagePlus("Structure result hw="+(r),imageStack);
		image.show();
	}

	
	
	 

	
	/**
	 * 
	 * This method is helper function for both applyFilter and run method
	 * @param ip input image
	 * @param sp gaussian scale space
	 * @param sigma filter sigma
	 */
	private ImageStack filter(ImageProcessor ip, GScaleSpace sp,  ImageStack imageStack){

		ip.snapshot();

		if (!isFloat) 
			ip=ip.toFloat(0, null);

		pass++;
		double sigma=sp.getSigma();
		//System.out.println(settings.get(LEN)+"MG");
		//GScaleSpace sp=new GScaleSpace(sigma);
		float[] kernx= sp.gauss1D();
		//System.out.println("kernx :"+kernx.length);
		GScaleSpace.flip(kernx);	

		float[] kern_diff1=sp.diffGauss1D();
		//System.out.println("kernx1:"+kern_diff1.length);
		GScaleSpace.flip(kern_diff1);
		
 

		kernel=new float[4][];
		kernel[0]=kernx;
 
		kernel[1]=kern_diff1;

 

		int sz= sp.getSize();
		if (debug && pass==1) {
			FloatProcessor fpkern2=new FloatProcessor(sz,sz);

			float[][] disp= new float[2][];

			disp[0]=GScaleSpace.joinXY(kernel, 0, 1);
			disp[1]=GScaleSpace.joinXY(kernel, 1, 0);

			for (int i=0; i<sz*sz; i++)
				fpkern2.setf(i, disp[0][i]+ disp[1][i]);

			new ImagePlus("kernel sep",fpkern2).show();


		}

		FloatProcessor fpaux= (FloatProcessor) ip;

		Conv cnv=new Conv();

		FloatProcessor gradx=(FloatProcessor) fpaux.duplicate();
		FloatProcessor grady=(FloatProcessor) fpaux.duplicate();
		
 
		cnv.convolveFloat1D(gradx, kern_diff1, Ox);
		cnv.convolveFloat1D(gradx, kernx, Oy);

		cnv.convolveFloat1D(grady, kern_diff1, Oy);
		cnv.convolveFloat1D(grady, kernx, Ox);
 
		int width=ip.getWidth();
		int height=ip.getHeight();

 
		
		FloatProcessor pamp=new FloatProcessor(width, height); // amplitude of gradient
		FloatProcessor phase=new FloatProcessor(width, height); // phase of gradient
		
		FloatProcessor eigen1=new FloatProcessor(width, height); // eigenvalue 1
		FloatProcessor eigen2=new FloatProcessor(width, height); // eigenvalue 2
		FloatProcessor coher=new FloatProcessor(width, height); // coherence 

		for (int i=0; i<width*height; i++) {
			double gx=gradx.getf(i);
			double gy=grady.getf(i);

			/*
			 *  components of the Structure Tensor
			 */
			double amp=(gx*gx+gy*gy);
			if (amp==0) amp+=1e-6;
			
			final double trace=amp;
			final double det=gx*gx*gy*gy- gx*gy*gx*gy;
			final double disc= sqrt(abs(trace*trace-4.0*det));
			final double ee1=0.5*(trace+disc);
			final double ee2=0.5*(trace-disc);

			final double l1=max(ee1, ee2);
			final double l2=min(ee1, ee2);
		    double coh=0;
			if (l1+l2>0) 
			   coh=(l1-l2)/(l1+l2);
	 
			pamp.setf(i, (float) sqrt(amp));
			
			double phase1=sqrt(gy/amp);
				//	phase1=asin(phase1);
			phase.setf(i, (float) phase1);

			eigen1.setf(i, (float) l1);
			eigen2.setf(i, (float) l2);
			coher.setf(i, (float) coh);
		}

		if (fulloutput) {
			imageStack.addSlice(FILTER_KEY+"_X_diff_"+sigma, gradx);
			imageStack.addSlice(FILTER_KEY+"_Y_diff_"+sigma, grady);
		}

		imageStack.addSlice(FILTER_KEY+"_Amp_"+sigma, pamp);
		imageStack.addSlice(FILTER_KEY+"_Phase_"+sigma, phase);
		imageStack.addSlice(FILTER_KEY+"_Coh_"+sigma, phase);
		imageStack.addSlice(FILTER_KEY+"_E2_"+sigma, eigen2);
		imageStack.addSlice(FILTER_KEY+"_E1_"+sigma, eigen1);
 
		eigen2.resetMinAndMax();
 
		return imageStack;
	}

	/**
	 * 
	 * This method is helper function for both applyFilter and run method
	 * @param ip input image
	 * @param sp gaussian scale space
	 * @param sigma filter sigma
	 */
	private ImageStack filter2(ImageProcessor ip, GScaleSpace sp, GScaleSpace sp2, ImageStack imageStack){

		ip.snapshot();

		if (!isFloat) 
			ip=ip.toFloat(0, null);

		pass++;
		//System.out.println(settings.get(LEN)+"MG");
		//GScaleSpace sp=new GScaleSpace(sigma);
		float[] kernx= sp.gauss1D();
		//System.out.println("kernx :"+kernx.length);
		GScaleSpace.flip(kernx);	

		float[] kern_diff1=sp.diffGauss1D();
		//System.out.println("kernx1:"+kern_diff1.length);
		GScaleSpace.flip(kern_diff1);
		
		double sigma=sp2.getSigma();

		kernel=new float[4][];
		kernel[0]=kernx;
 
		kernel[1]=kern_diff1;

 

		int sz= sp.getSize();
		if (debug && pass==1) {
			FloatProcessor fpkern2=new FloatProcessor(sz,sz);

			float[][] disp= new float[2][];

			disp[0]=GScaleSpace.joinXY(kernel, 0, 1);
			disp[1]=GScaleSpace.joinXY(kernel, 1, 0);

			for (int i=0; i<sz*sz; i++)
				fpkern2.setf(i, disp[0][i]+ disp[1][i]);

			new ImagePlus("kernel sep",fpkern2).show();


		}

		FloatProcessor fpaux= (FloatProcessor) ip;

		Conv cnv=new Conv();

		FloatProcessor gradx=(FloatProcessor) fpaux.duplicate();
		FloatProcessor grady=(FloatProcessor) fpaux.duplicate();
		
 
		cnv.convolveFloat1D(gradx, kern_diff1, Ox);
		cnv.convolveFloat1D(gradx, kernx, Oy);

		cnv.convolveFloat1D(grady, kern_diff1, Oy);
		cnv.convolveFloat1D(grady, kernx, Ox);
 
		int width=ip.getWidth();
		int height=ip.getHeight();

 
		
		FloatProcessor pamp=new FloatProcessor(width, height); // amplitude of gradient
		FloatProcessor phase=new FloatProcessor(width, height); // phase of gradient
		
		FloatProcessor gx2=new FloatProcessor(width, height); // gx^2 
		FloatProcessor gy2=new FloatProcessor(width, height); // gy^2 
		FloatProcessor gxy=new FloatProcessor(width, height); // gx*gy 
		
		FloatProcessor eigen1=new FloatProcessor(width, height); // eigenvalue 1
		FloatProcessor eigen2=new FloatProcessor(width, height); // eigenvalue 2
		FloatProcessor coher=new FloatProcessor(width, height); // coherence 

		for (int i=0; i<width*height; i++) {
			double gx=gradx.getf(i);
			double gy=grady.getf(i);

			/*
			 *  Gradient analysis
			 */
 			final double x2=gx*gx;
			gx2.setf(i, (float) x2);
			final double y2=gy*gy;
			gy2.setf(i, (float) y2);
			
			final double xy=gx*gy;
			gxy.setf(i, (float) xy);
			double amp=sqrt(gx*gx+gy*gy);
			if (amp==0) amp+=1e-6;
			pamp.setf(i, (float)  (amp));
			
			double phase1= (gy/amp);
					phase1=asin(phase1);
			phase.setf(i, (float) phase1);
		}
		
		if (fulloutput) {
			imageStack.addSlice(FILTER_KEY+"_X_diff_"+sigma, gradx);
			imageStack.addSlice(FILTER_KEY+"_Y_diff_"+sigma, grady);
		}
		imageStack.addSlice(FILTER_KEY+"_Amp_"+sigma, pamp);
		imageStack.addSlice(FILTER_KEY+"_Phase_"+sigma, phase);
		
		// second smoothing step
		
		kernx= sp2.gauss1D();
		//System.out.println("kernx :"+kernx.length);
		GScaleSpace.flip(kernx);	

		kern_diff1=sp2.diffGauss1D();
		//System.out.println("kernx1:"+kern_diff1.length);
		GScaleSpace.flip(kern_diff1);
		
		//double sigma=sp.getSigma();
		
		cnv.convolveFloat1D(gx2, kern_diff1, Ox);
		cnv.convolveFloat1D(gx2, kernx, Oy);

		cnv.convolveFloat1D(gy2, kern_diff1, Oy);
		cnv.convolveFloat1D(gy2, kernx, Ox);
		
		cnv.convolveFloat1D(gxy, kern_diff1, Oy);
		cnv.convolveFloat1D(gxy, kernx, Ox);
		
		for (int i=0; i<width*height; i++) {
			double xx=gx2.getf(i);
			double yy=gy2.getf(i);
			double xy=gxy.getf(i);
 
			final double trace=xx+yy;
			final double det=xx*yy- xy*xy;
			final double disc= sqrt(abs(trace*trace-4.0*det));
			final double ee1=0.5*(trace+disc);
			final double ee2=0.5*(trace-disc);

			final double l1=max(ee1, ee2);
			final double l2=min(ee1, ee2);
		    double coh=0;
			if (l1+l2>0) 
			   coh=(l1-l2)/(l1+l2);
	  
			eigen1.setf(i, (float) l1);
			eigen2.setf(i, (float) l2);
			coher.setf(i, (float) coh);
		}

		


		imageStack.addSlice(FILTER_KEY+"_Coh_"+sigma, phase);
		imageStack.addSlice(FILTER_KEY+"_E2_"+sigma, eigen2);
		imageStack.addSlice(FILTER_KEY+"_E1_"+sigma, eigen1);
 
		eigen2.resetMinAndMax();
 
		return imageStack;
	}

	/**
	 * @param i
	 * @return
	 */
	public float[] getKernel(int i) {
		return kernel[i];
	}

	/* (non-Javadoc)
	 * @see ij.plugin.filter.ExtendedPlugInFilter#showDialog(ij.ImagePlus, java.lang.String, ij.plugin.filter.PlugInFilterRunner)
	 */
	@Override
	public int showDialog(ImagePlus imp, String command, PlugInFilterRunner pfr) {
		this.pfr = pfr;
		int r = (sz-1)/2;
		GenericDialog gd=new GenericDialog("Structure Tensor " + version);

		gd.addNumericField("half width (noise)", r, 2);
		gd.addNumericField("half width (signal)", 2*r, 1);
		gd.addCheckbox("Show kernel", debug);
		gd.addCheckbox("Full output", fulloutput);	
		if (cal!=null) {
			if (!cal.getUnit().equals("pixel"))
				gd.addCheckbox("units ( "+cal.getUnit() + " )", doCalib); 
		}		

		gd.addPreviewCheckbox(pfr);
		gd.addDialogListener(this);
		gd.setResizable(false);
		gd.showDialog();

		if (gd.wasCanceled()) {			
			return DONE;
		}

		return IJ.setupDialog(imp, flags);
	}




	// Called after modifications to the dialog. Returns true if valid input.
	/* (non-Javadoc)
	 * @see ij.gui.DialogListener#dialogItemChanged(ij.gui.GenericDialog, java.awt.AWTEvent)
	 */
	@Override
	public boolean dialogItemChanged(GenericDialog gd, AWTEvent e) {
		double r = (int)(gd.getNextNumber());
		double r2= (float) (gd.getNextNumber());
		debug = gd.getNextBoolean();
		fulloutput = gd.getNextBoolean();
		if (cal!=null)
			doCalib=gd.getNextBoolean();

		if (doCalib) {
			r= (r/cal.pixelWidth);
		}
		sz =  (2*(int)r+1);
		sz2 =  (2*(int)r2+1);
		if (gd.wasCanceled()) {

			return false;
		}
		return r>0;
		//return sigma>0;
	}


	/* (non-Javadoc)
	 * @see ij.plugin.filter.ExtendedPlugInFilter#setNPasses(int)
	 */
	@Override
	public void setNPasses (int nPasses) {
		this.nPasses = nPasses;
	}

	/* Saves the current settings of the plugin for further use
	 * 
	 *
	 * @param prefs
	 */
	public void savePreferences(Properties prefs) {
		prefs.put(LEN, Integer.toString(sz));
		// prefs.put(SIGMA, Float.toString(sigma));

	}

	//@Override
	public Map<String, String> getDefaultSettings() {

		settings.put(LEN, Integer.toString(sz));
		settings.put(MAX_LEN, Integer.toString(max_sz));
		settings.put(FULL_OUTPUT, Boolean.toString(fulloutput));

		return this.settings;
	}

	//@Override
	public boolean reset() {
		sz= Prefs.getInt(LEN, 2);
		max_sz= Prefs.getInt(MAX_LEN, 8);
		fulloutput= Prefs.getBoolean(FULL_OUTPUT, true);
		return true;
	}

	//@Override
	public boolean updateSettings(Map<String, String> settingsMap) {
		sz=Integer.parseInt(settingsMap.get(LEN));
		max_sz=Integer.parseInt(settingsMap.get(MAX_LEN));
		fulloutput= Boolean.parseBoolean(settingsMap.get(FULL_OUTPUT));
		
		return true;
	}

	//@Override
	public String getKey() {
		return this.FILTER_KEY;
	}
 

	//@Override
	public boolean isEnabled() {
		return isEnabled;
	}

	//@Override
	public void setEnabled(boolean isEnabled) {
		this.isEnabled= isEnabled;
	}
	
//	@Override
	public String getName() {
		return this.FILTER_NAME;
	}
	

	private double logKernel(double x){
		final double x2=x*x;
		return -(x)* exp(-0.5*x2)/(2.0*sqrt(PI));
	}
	
//	@Override
	public double[][] kernelData() {
		final int n=40;
		double [][] data=new double[2][n];
		data[0]=SUtils.linspace(-10.0, 10.0, n);
		for(int i=0; i<n; i++){
			data[1][i]=logKernel(data[0][i]);
		}
		return data;
	}

	public static void main (String[] args) {
		//StructureT_Filter_ filter=new StructureT_Filter_();
		//System.out.println("annotated fields");
		//System.out.println(filter.getAnotatedFileds());
 
		//new ImagePlus("convovled",output).show();
	 
		try {

			File f=new File(args[0]);

			if (f.exists() && f.isDirectory() ) {
				System.setProperty("plugins.dir", args[0]);
				new ImageJ();
			} else {
				throw new IllegalArgumentException();
			}
		}
		catch (Exception ex) {
			IJ.log("plugins.dir misspecified\n");
			ex.printStackTrace();
		}
		
	}


}
