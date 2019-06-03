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
import ijaux.scale.GScaleSpace;
import ijaux.scale.SUtils;

import java.awt.*;
import java.awt.event.*;
import java.util.*;

import dsp.Conv;

/**
 * @version 	1.5		date 23 Sept 2013
 *				- isotropic correction
 * 				1.0		date 23 Jul 2013 
 * 				Based on Mexican_Hat_Filter v 2.2
 * 				- common functionality is refactored in a library class
 * 				
 *   
 * 
 * @author Dimiter Prodanov
 * 		  IMEC
 *
 *
 * @contents
 * This pluign convolves an image with a Bi-Laplacian of Gaussian (BoG) filter
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

 
public class LoG2_Filter_ implements ExtendedPlugInFilter, DialogListener {

	private PlugInFilterRunner pfr=null;

	final int flags=DOES_ALL+SUPPORTS_MASKING+KEEP_PREVIEW;
	private String version="1.5";
	private int nPasses=1;
	private int pass;
	
	public final static String SIGMA="LOG_sigma", LEN="G_len", ISO="G_iso", ISSEP="G_SEP";

	private static int sz= Prefs.getInt(LEN, 9);
	//private static float sigma=(float) Prefs.getDouble(SIGMA, 2.0f);
	private float[][] kernel=null;

	private ImagePlus image=null;
	public static boolean debug=IJ.debugMode;

	public static boolean sep= Prefs.getBoolean(ISSEP, false);
	private static boolean isiso= Prefs.getBoolean(ISO, true);

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
		//GScaleSpace sp=new GScaleSpace(r);
		GScaleSpace sp=new GScaleSpace(r,3.0f);
		float[] kernx= sp.gauss1D();
		SUtils.flip(kernx);		

		float[] kern_diff_4= sp.diffNGauss1D(4);
		SUtils.flip(kern_diff_4);
		
		float[] kern_diff_2= sp.diffNGauss1D(2);
		SUtils.flip(kern_diff_2);		
		
		kernel=new float[4][];
		kernel[0]=kernx;
		kernel[1]=kern_diff_4;
		kernel[2]=kern_diff_2;
		float[] kernel2=sp.computeLapNKernel2D(2); // 2D kernel computation
		kernel[3]=kernel2;
			
		int sz= sp.getSize();
		
		if (debug)
			System.out.println("sz " +sz);
		
		float[][] disp= new float[3][];

		disp[0]=GScaleSpace.joinXY(kernel, 0, 1);
		disp[1]=GScaleSpace.joinXY(kernel, 1, 0);
		
		
		if (debug && pass==1) {
			FloatProcessor fp=new FloatProcessor(sz,sz);
			if (isiso) {
				disp[2]=GScaleSpace.joinXY(kernel, 2, 2);
				for (int i=0; i<sz*sz; i++)
					fp.setf(i, disp[0][i]+ disp[1][i] + 2*disp[2][i]  );
				 
			} else {
				for (int i=0; i<sz*sz; i++)
					fp.setf(i, disp[0][i]+ disp[1][i] );
			}
			new ImagePlus("kernel sep",fp).show();
			if (!sep) {
				FloatProcessor fp2=new FloatProcessor(sz,sz, kernel2);
				new ImagePlus("kernel 2D",fp2).show();
			}
		}
		long time=-System.nanoTime();	
		
		FloatProcessor fpaux= (FloatProcessor) ip;
	 		
		Conv cnv=new Conv();
		if (sep) {
			if (isiso) {
				FloatProcessor fpauxiso=(FloatProcessor) fpaux.duplicate();
								
				cnv.convolveSemiSep(fpaux, kernx, kern_diff_4);	
				for (int i=0; i<kern_diff_2.length; i++)
					kern_diff_2[i]*=Math.sqrt(2.0);
				cnv.convolveFloat1D(fpauxiso, kern_diff_2, 0); //Ox
				cnv.convolveFloat1D(fpauxiso, kern_diff_2, 1); //Oy
				
				fpaux.copyBits(fpauxiso, 0, 0, Blitter.ADD);
				System.out.println("separable & isotropic computation");
			} else {
				cnv.convolveSemiSep(fpaux, kernx, kern_diff_4);	
				System.out.println("separable & non-isotropic computation");
			}
		} else {	
			if (isiso) {			
				System.out.println("non-separable & isotropic computation");
			} else {
				for (int i=0; i<sz*sz; i++)
					kernel2[i]=disp[0][i]+ disp[1][i];
				System.out.println("non-separable & non-isotropic computation");
			} // end else
			cnv.convolveFloat(fpaux, kernel2, sz, sz);
		} // end else
	 
		time+=System.nanoTime();
		time/=1000.0f;
		System.out.println("elapsed time: " + time +" us");
		fpaux.resetMinAndMax();	
		
		/*if (convert) {
			float[] minmax=findMinAndMax(fpaux);
			final double d1 = -(fmin* minmax[1] - fmax* minmax[0])/(fmin - fmax);
			dr = dr/ (minmax[1] - minmax[0]);
			System.out.println("contrast adjustment y=ax+b \n " +
					" b " +d1 +" a " + dr);
				
			contrastAdjust(fpaux, dr, d1);
		}*/
		
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
		GenericDialog gd=new GenericDialog("Bi-LoG (BoG) " + version);
		gd.addNumericField("half width", r, 1);
		//gd.addNumericField("sigma", sigma, 1);
		gd.addCheckbox("Show kernel", debug);
		gd.addCheckbox("Separable", sep);
		gd.addCheckbox("isotropic", isiso);
		/*if (hasRoi)
			gd.addCheckbox("Brightness correct", true);*/
		
		gd.addPreviewCheckbox(pfr);
		gd.addDialogListener(this);
		gd.showDialog();
		pixundo=imp.getProcessor().getPixelsCopy();
		if (gd.wasCanceled()) {			
			//image.repaintWindow();
			return DONE;
		}
		/*if (!IJ.isMacro())
			staticSZ = sz;*/
		return IJ.setupDialog(imp, flags);
	}
	
	private Object pixundo;
	//private boolean convert=false;
	
	// Called after modifications to the dialog. Returns true if valid input.
	/* (non-Javadoc)
	 * @see ij.gui.DialogListener#dialogItemChanged(ij.gui.GenericDialog, java.awt.AWTEvent)
	 */
	public boolean dialogItemChanged(GenericDialog gd, AWTEvent e) {
		int r = (int)(gd.getNextNumber());
		//sigma = (float) (gd.getNextNumber());
		debug = gd.getNextBoolean();
		sep = gd.getNextBoolean();
		isiso = gd.getNextBoolean();
		//convert=gd.getNextBoolean();
		 
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
  		prefs.put(ISO, Boolean.toString(isiso));
  		prefs.put(ISSEP, Boolean.toString(sep));
        // prefs.put(SIGMA, Float.toString(sigma));

   }

}
