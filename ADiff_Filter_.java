import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.gui.DialogListener;
import ij.gui.Roi;
import ij.measure.Calibration;
import ij.plugin.filter.Convolver;
import ij.plugin.filter.ExtendedPlugInFilter;
import ij.plugin.filter.PlugInFilterRunner;
import ij.process.Blitter;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.LUT;
import ijaux.scale.*;

import java.awt.*;
import java.awt.event.*;
import java.util.*;

import dsp.Conv;

/**
 * @version 	1.0 10 Oct 2013
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

 /*
  *  naive implementation
  */
public class ADiff_Filter_ implements ExtendedPlugInFilter, DialogListener {

	private PlugInFilterRunner pfr=null;

	final int flags=DOES_ALL+KEEP_PREVIEW+NO_CHANGES;
	private String version="2.0";
	private int nPasses=1;
	private int pass;
 
	public final static String SIGMA="LOG_sigma", LEN="G_len", NI="G_ier", ORT="G_ortd";
 
	private static int sz= Prefs.getInt(LEN, 9);
	private static int ni= Prefs.getInt(NI, 10);
	//private static float sigma=(float) Prefs.getDouble(SIGMA, 2.0f);

	private float[][] kernel=null;

	private ImagePlus image=null;
	public static boolean debug=true;//IJ.debugMode;
	public static boolean ortdir = Prefs.getBoolean(ORT, false);
	public static boolean fulloutput=false;

	public boolean isFloat=false;
	
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
		return  flags;
	}

	final int Ox=0, Oy=1, Oz=2;
 /*
  * (non-Javadoc)
  * @see ij.plugin.filter.PlugInFilter#run(ij.process.ImageProcessor)
  */
	@Override
	public void run(ImageProcessor ip) {
		ip.snapshot();
		 	
		if (!isFloat) 
			ip=ip.toFloat(0, null);

		int width=ip.getWidth();
		int height=ip.getHeight();
		
		pass++;
		int r = (sz-1)/2;
		GScaleSpace sp=new GScaleSpace(r);
		//GScaleSpace sp=new GScaleSpace(sigma);
		float[] kernx= sp.gauss1D();
		SUtils.flip(kernx);		
		float[] kern_diff2= sp.diff2Gauss1D();
		SUtils.flip(kern_diff2);
		
		float[] kern_diff1=sp.diffGauss1D();
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
		
		final double step=0.5*Math.E*sp.getScale();
		
		Conv cnv=new Conv();
		final FloatProcessor fpaux= (FloatProcessor) ip;
		
		long time=-System.nanoTime();
		
		for (int k=0; k<ni; k++) {

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

			for (int i=0; i<width*height; i++) {
				double gx=gradx.getf(i);
				double gy=grady.getf(i);

				double gxy=lap_xy.getf(i);

				double gxx=lap_xx.getf(i);
				double gyy=lap_yy.getf(i);

				double lx=2.0*gx*gy*gxy;

				gx*=gx;
				gy*=gy;		
			
				double amp=(gx+gy);	
				//double amp=2.0*(gx+gy);		
				
				//directional laplacian components
				//float lt=(float)((dt-lx)/amp); 
				//float ot=(float)((dx+lx)/amp);
				if (ortdir) {
					double dx=gx*gxx+gy*gyy;
					float ot=0;
					if (dx!=0 || amp!=0) {
						ot=(float)((dx+lx)/amp*step);
					}
					float lap_tarr=ot +fpaux.getf(i);
					fpaux.setf(i, lap_tarr);
				} else {
					double dt=gy*gxx+gx*gyy;
					float lt=0;
					if (dt!=0 || amp!=0) {
						lt=(float)((dt-lx)/amp*step);
					}
					float lap_tarr=lt +fpaux.getf(i);
					// add oscillation suppression step
					
					fpaux.setf(i, lap_tarr);
				}

			} // end for
			
			//new ImagePlus("step "+k,(ImageProcessor)fpaux.duplicate()).show();
		} // end for
		

		
		time+=System.nanoTime();
		time/=1000.0f;
		System.out.println("elapsed time: " + time +" us");

		fpaux.resetMinAndMax();
		image.setProcessor(fpaux);
		image.updateAndDraw();
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
		GenericDialog gd=new GenericDialog("Anisotropic diffusion " + version);
		gd.addNumericField("half width", r, 1);
		gd.addNumericField("iterations", ni, 1);
		//gd.addNumericField("sigma", sigma, 1);
		gd.addCheckbox("orthogonal", ortdir);
		gd.addCheckbox("debug", debug);
	//	gd.addCheckbox("Full output", fulloutput);
	
		
		gd.addPreviewCheckbox(pfr);
		gd.addDialogListener(this);
		gd.showDialog();
		pixundo=imp.getProcessor().getPixelsCopy();
		if (gd.wasCanceled()) {			
			return DONE;
		}

		return IJ.setupDialog(imp, flags);
	}
	
	private Object pixundo;
	
	
	// Called after modifications to the dialog. Returns true if valid input.
	/* (non-Javadoc)
	 * @see ij.gui.DialogListener#dialogItemChanged(ij.gui.GenericDialog, java.awt.AWTEvent)
	 */
	public boolean dialogItemChanged(GenericDialog gd, AWTEvent e) {
		int r = (int)(gd.getNextNumber());
		ni = (int)(gd.getNextNumber());
		//sigma = (float) (gd.getNextNumber());
		ortdir = gd.getNextBoolean();
		debug = gd.getNextBoolean();
		//fulloutput = gd.getNextBoolean();
		
		 
		sz = 2*r+1;
		if (gd.wasCanceled()) {
			ImageProcessor proc=image.getProcessor();
			proc.setPixels(pixundo);
			proc.resetMinAndMax();
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
	   		prefs.put(NI, Integer.toString(ni));
	   		prefs.put(ORT, Boolean.toString(ortdir));
         // prefs.put(SIGMA, Float.toString(sigma));

   }

}
