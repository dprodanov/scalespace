import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.gui.DialogListener;
import ij.measure.Calibration;
import ij.plugin.filter.ExtendedPlugInFilter;
import ij.plugin.filter.PlugInFilterRunner;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ijaux.scale.*;

import java.awt.*;
import java.util.*;

import dsp.Conv;

import static java.lang.Math.*;

/**
 * @version 	1.2.5 11 Sept 2015
 * 
 * 				1.2 11 Aug 2015
 * 				1.1 27 Jun 2015
 * 				1.0  6 Oct 2014 
 * 				
 *   
 * 
 * @author Dimiter Prodanov
 * 		  IMEC
 *
 *
 * @contents
 * The plugin performs anisotropic LoG filtering. The principle is based on Michael Broadhead
 * http://works.bepress.com/cgi/viewcontent.cgi?article=1017&context=michael_broadhead 
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

 
public class ALoG_Filter_ implements ExtendedPlugInFilter, DialogListener {

	private PlugInFilterRunner pfr=null;

	final int flags=DOES_ALL+KEEP_PREVIEW+ NO_CHANGES;
	private String version="2.0";
	private int nPasses=1;
	private int pass;
 
	public final static String SIGMA="LOG_sigma", LEN="G_len";
 
	private static int sz= Prefs.getInt(LEN, 9);
	//private static float sigma=(float) Prefs.getDouble(SIGMA, 2.0f);

	private float[][] kernel=null;

	private ImagePlus image=null;
	public static boolean debug=true;//IJ.debugMode;

	public static boolean fulloutput=false;

	private boolean isFloat=false;
	
	private boolean hasRoi=false;
	
	/**
	 * 
	 */
	/* (non-Javadoc)
	 * @see ij.plugin.filter.PlugInFilter#setup(java.lang.String, ij.ImagePlus)
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

	private boolean doCalib = false;
	private Calibration cal=null;
 /*
  * (non-Javadoc)
  * @see ij.plugin.filter.PlugInFilter#run(ij.process.ImageProcessor)
  */
	@Override
	public void run(ImageProcessor ip) {
		ip.snapshot();
		 	
		if (!isFloat) 
			ip=ip.toFloat(0, null);
		
		pass++;
		int r = (sz-1)/2;
		GScaleSpace sp=new GScaleSpace(r);
		//GScaleSpace sp=new GScaleSpace(sigma);
		float[] kernx= sp.gauss1D();
		System.out.println("kernx :"+kernx.length);
		SUtils.flip(kernx);		
		float[] kern_diff2= sp.diff2Gauss1D();
		System.out.println("kernx2 :"+kern_diff2.length);
		SUtils.flip(kern_diff2);
		
		float[] kern_diff1=sp.diffGauss1D();
		//float[] kern_diff1=sp.diffNGauss1D( 1);
		System.out.println("kernx1: "+kern_diff1.length);
		SUtils.flip(kern_diff1);
		
		kernel=new float[4][];
		kernel[0]=kernx;
		kernel[1]=kern_diff2;
		kernel[2]=kern_diff1;
		
		float[] kernel2=sp.computeDiff2Kernel2D();
		kernel[3]=kernel2;
		SUtils.flip(kernel2);  // symmetric but this is the correct way
		
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
		long time=-System.nanoTime();	
		FloatProcessor fpaux= (FloatProcessor) ip;
	 		
		Conv cnv=new Conv();

		FloatProcessor gradx=(FloatProcessor) fpaux.duplicate();
		FloatProcessor grady=(FloatProcessor) fpaux.duplicate();
		FloatProcessor lap_xx=(FloatProcessor) fpaux.duplicate();
		FloatProcessor lap_yy=(FloatProcessor) fpaux.duplicate();
		FloatProcessor lap_xy=(FloatProcessor) fpaux.duplicate();

		cnv.convolveFloat1D(gradx, kern_diff1, Ox);
		cnv.convolveFloat1D(gradx, kernx, Oy);

		cnv.convolveFloat1D(grady, kern_diff1, Oy);
		cnv.convolveFloat1D(grady, kernx, Ox);

		cnv.convolveFloat1D(lap_xx, kern_diff2, Ox);
		cnv.convolveFloat1D(lap_xx, kernx, Oy);

		cnv.convolveFloat1D(lap_yy, kern_diff2, Oy);
		cnv.convolveFloat1D(lap_yy, kernx, Ox);

		cnv.convolveFloat1D(lap_xy, kern_diff1, Oy);
		cnv.convolveFloat1D(lap_xy, kern_diff1, Ox);
		int width=ip.getWidth();
		int height=ip.getHeight();

		FloatProcessor lap_t=new FloatProcessor(width, height); // tangential component
		FloatProcessor lap_o=new FloatProcessor(width, height); // orthogonal component
		FloatProcessor pamp=new FloatProcessor(width, height); // amplitude of gradient
		FloatProcessor phase=new FloatProcessor(width, height); // phase of gradient
		
		for (int i=0; i<width*height; i++) {
			double gx=gradx.getf(i);
			double gy=grady.getf(i);

			double gxy=lap_xy.getf(i);

			double gxx=lap_xx.getf(i);
			double gyy=lap_yy.getf(i);

			double lx=2.0f*gx*gy*gxy;

			gx*=gx;
			gy*=gy;		
			double dt=gy*gxx+gx*gyy;
			double dx=gx*gxx+gy*gyy;
		 		
			double amp=(gx+gy)+ 1e-6;
			
			if (abs(amp) > 1e-4) { 	
				float lt=(float)((dt-lx)/amp);
				float ot=(float)((dx+lx)/amp);
				
				if (abs(lt) <1e-8) lt=0;
				if (abs(ot) <1e-8) ot=0;	
				
				lap_t.setf(i, lt);
				lap_o.setf(i, ot);
			} 
			
			pamp.setf(i, (float) sqrt(amp));
			
			double phase1=sqrt(gy/amp);
				//	phase1=asin(phase1);
			phase.setf(i, (float) phase1);

		}
			
		ImageStack is=new ImageStack(width,height);
		int apos=4;
		
		if (fulloutput) {
			is.addSlice("X diff", gradx);
			is.addSlice("Y diff", grady);
			is.addSlice("XX diff", lap_xx);
			is.addSlice("YY diff", lap_yy);
			is.addSlice("XY diff", lap_xy);
			apos=5;
		}
		is.addSlice("Amp", pamp);
		is.addSlice("Phase", phase);
		is.addSlice("Lap T", lap_t);
		lap_o.resetMinAndMax();
		is.addSlice("Lap O", lap_o);
		
		time+=System.nanoTime();
		time/=1000.0f;
		System.out.println("elapsed time: " + time +" us");
		System.out.println("sigma: " + sp.getSigma() + 
						   " scale: " + sp.getScale() + 
						   " kernel size: "+ sp.getSize()
						   );
		
		String stackloc="";
		
		if (image.getStackSize()> 1) {
			stackloc=" z= "+image.getCurrentSlice();
			System.out.println("stack location "+stackloc);
		}
		
		
		image=new ImagePlus("ALoG result hw="+r+stackloc,is);
		image.show();
		image.setPosition(apos);
		image.getProcessor().resetMinAndMax();	
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
		GenericDialog gd=new GenericDialog("Anisotropic LoG " + version);
	
		gd.addNumericField("half width", r, 2);
		//gd.addNumericField("sigma", sigma, 1);
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
		
		//pixundo=imp.getProcessor().getPixelsCopy();
		if (gd.wasCanceled()) {			
			return DONE;
		}

		return IJ.setupDialog(imp, flags);
	}
	
	 
	
	
	// Called after modifications to the dialog. Returns true if valid input.
	/* (non-Javadoc)
	 * @see ij.gui.DialogListener#dialogItemChanged(ij.gui.GenericDialog, java.awt.AWTEvent)
	 */
	public boolean dialogItemChanged(GenericDialog gd, AWTEvent e) {
		double r = (int)(gd.getNextNumber());
		//sigma = (float) (gd.getNextNumber());
		debug = gd.getNextBoolean();
		fulloutput = gd.getNextBoolean();
		if (cal!=null)
			doCalib=gd.getNextBoolean();

		if (doCalib) {
			r= (r/cal.pixelWidth);
		}
		sz =  (2*(int)r+1);
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
   public static void savePreferences(Properties prefs) {
	   		prefs.put(LEN, Integer.toString(sz));
         // prefs.put(SIGMA, Float.toString(sigma));

   }
   
   /*
   public static void main(String[] args) {
		// TODO Auto-generated method stub
		int sz=8;
		int r = (sz-1)/2;
		GScaleSpace sp=new GScaleSpace(r);
		//GScaleSpace sp=new GScaleSpace(sigma);
		float[] kernx= sp.gauss1D();
		System.out.println("kernx :"+kernx.length);

		float[] kern_diff2= sp.diff2Gauss1D();
		System.out.println("kernx2 :"+kern_diff2.length);
		
		
		//float[] kern_diff1=sp.diffGauss1D();
		float[] kern_diff1=sp.diffNGauss1D( 1);
		System.out.println("kernx1: "+kern_diff1.length);
	
		float[] kern_diff3=sp.diffGauss1D();
				 
		System.out.println("kernx1: "+kern_diff3.length);
		

	}
	*/
	

}
