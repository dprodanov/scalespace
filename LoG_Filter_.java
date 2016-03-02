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

import java.awt.*;
import java.awt.event.*;
import java.io.File;
import java.util.*;

import dsp.Conv;

/**
 * @version 	1.1	14 Oct 2013
 * 				- moved contratAdjust -> Conv
 * 				- changed brightness adjustment factor to sigma^2		
 * 				1.1 	18 Jul 2013
 * 				- refactoring
 * 				1.0		05 Feb 2013 
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
 * This pluign convolves an image with a Mexican Hat (Laplacian of Gaussian, LoG) filter
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

 
public class LoG_Filter_ implements ExtendedPlugInFilter, DialogListener {

	private PlugInFilterRunner pfr=null;

	private final int flags=DOES_ALL+SUPPORTS_MASKING+KEEP_PREVIEW;
	private String version="2.1";
	private int nPasses=1;
	private int pass;
	
	public static boolean debug=IJ.debugMode;
	public final static String SIGMA="LOG_sigma", LEN="G_len", ISSEP="G_SEP", SCNORM="G_SCNORM";
 
	private static int sz= Prefs.getInt(LEN, 9);
	//private static float sigma=(float) Prefs.getDouble(SIGMA, 2.0f);
	public static boolean sep= Prefs.getBoolean(ISSEP, false);
	
	public static boolean scnorm= Prefs.getBoolean(SCNORM, false);
	
	private float[][] kernel=null;
	private ImagePlus image=null;	
	private boolean isFloat=false;	
	private boolean hasRoi=false;
	
	/*
	 * @param args - args[0] should point to the folder where the plugins are installed 
	 */
	public static void main(String[] args) {

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
		int r = sz;//(sz-1)/2;
		GScaleSpace sp=new GScaleSpace(r);
		//GScaleSpace sp=new GScaleSpace(sigma);
		float[] kernx= sp.gauss1D();
		GScaleSpace.flip(kernx);		
		float[] kern_diff= sp.diff2Gauss1D();
		GScaleSpace.flip(kern_diff);
		
		System.out.println("scnorm "+scnorm);
		if (scnorm) {
			double gamma=sp.getSigma(); 	 
			for (int i=0; i<kern_diff.length; i++) {
				kern_diff[i]=(float) (kern_diff[i]*gamma);
				kernx[i]=(float) (kernx[i]*gamma);
			}
		}
		kernel=new float[3][];
		kernel[0]=kernx;
		kernel[1]=kern_diff;

		float[] kernel2=sp.computeDiff2Kernel2D();
		if (scnorm) {
			double gamma=sp.getScale();
			for (int i=0; i<kern_diff.length; i++) {
				kernel2[i]=(float) (kernel2[i]*gamma);
			}
		}
		kernel[2]=kernel2;
		GScaleSpace.flip(kernel2);  // symmetric but this is the correct way
 	
		int sz= sp.getSize();
		if (debug ) {
			FloatProcessor fp=new FloatProcessor(sz,sz);

			float[][] disp= new float[2][];

			disp[0]=GScaleSpace.joinXY(kernel, 0, 1);
			disp[1]=GScaleSpace.joinXY(kernel, 1, 0);

			for (int i=0; i<sz*sz; i++)
				fp.setf(i, disp[0][i]+ disp[1][i]);

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
			cnv.convolveSemiSep(fpaux, kernx, kern_diff);			
		} else {		 
			cnv.convolveFloat(fpaux, kernel2, sz, sz);
		}
	 
		time+=System.nanoTime();
		time/=1000.0f;
		System.out.println("elapsed time: " + time +" us");
		fpaux.resetMinAndMax();	
		
		if (convert) {
 
			final double d1=0;
			final double dr=sp. getScale();	
			System.out.println("linear contrast adjustment y=ax+b \n " +
					" b= " +d1 +" a= " + dr);
			
			Conv.contrastAdjust(fpaux, dr, d1);
		}
		
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
		GenericDialog gd=new GenericDialog("Mex. Hat " + version);
		//gd.addMessage("half width");
		gd.addNumericField("hw", r, 1);
		//gd.addNumericField("sigma", sigma, 1);
		gd.addCheckbox("Show kernel", debug);
		gd.addCheckbox("Separable", sep);
		gd.addCheckbox("Scale normalize", scnorm);
		if (hasRoi)
			gd.addCheckbox("Brightness correct", true);
		
		gd.addPreviewCheckbox(pfr);
		gd.addDialogListener(this);
		gd.showDialog();
		
		pixundo=imp.getProcessor().getPixelsCopy();

		return IJ.setupDialog(imp, flags);
	}
	
	private Object pixundo;
	private boolean convert=false;
	
	// Called after modifications to the dialog. Returns true if valid input.
	/* (non-Javadoc)
	 * @see ij.gui.DialogListener#dialogItemChanged(ij.gui.GenericDialog, java.awt.AWTEvent)
	 */
	public boolean dialogItemChanged(GenericDialog gd, AWTEvent e) {
		sz = (int)(gd.getNextNumber());
		//int r = (int)(gd.getNextNumber());
		//sigma = (float) (gd.getNextNumber());
		debug = gd.getNextBoolean();
		sep = gd.getNextBoolean();
		scnorm = gd.getNextBoolean();
		convert=gd.getNextBoolean();
		
		//boolean preview=gd.getNextBoolean();
		
		//sz = 2*r+1;
		if (gd.wasCanceled()) {
			ImageProcessor proc=image.getProcessor();
			proc.setPixels(pixundo);
			proc.resetMinAndMax();
			return false;
		}
		return sz>0;
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
	   		prefs.put(ISSEP, Boolean.toString(sep));
	   		prefs.put(SCNORM, Boolean.toString(scnorm));
         // prefs.put(SIGMA, Float.toString(sigma));

   }

}
