import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.gui.DialogListener;
import ij.gui.Roi;
import ij.plugin.filter.Convolver;
import ij.plugin.filter.ExtendedPlugInFilter;
import ij.plugin.filter.PlugInFilterRunner;
import ij.process.Blitter;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ijaux.scale.GScaleSpace;


import java.awt.*;
import java.awt.event.*;
import java.io.File;
import java.util.*;

import dsp.Conv;

/**
 * @version 	 0.5
 * 				 
 *   
 * 
 * @author Dimiter Prodanov
 * 		  IMEC
 *
 *
 * @contents
 * This pluign convolves an image with a ridge detector.
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
public class LoGN2D_ implements ExtendedPlugInFilter, DialogListener {

	private PlugInFilterRunner pfr=null;

	final int flags=DOES_ALL+CONVERT_TO_FLOAT+SUPPORTS_MASKING+KEEP_PREVIEW;
	private String version="2.0";
	private int nPasses=1;
	private int pass=0;
	public final static String SIGMA="LOG_sigma", LEN="G_len", 
			ISSEP="G_SEP", GN="G_Xn", GM="G_Yn", NORM="G_Nsc";
	 
	private static int sz = Prefs.getInt(LEN, 9);
	private float[][] kernel=null;
	public static boolean debug=IJ.debugMode;

	private static int n = Prefs.getInt(GN, 1);

	private static boolean sep=Prefs.getBoolean(ISSEP, false);
	private static boolean scnorm=Prefs.getBoolean(NORM, false);
	private ImagePlus image=null;	
	private boolean isFloat=false;
	private boolean isRGB=false;

	
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

	@Override
	public int setup(String arg, ImagePlus imp) {
		isFloat= (imp.getType()==ImagePlus.GRAY32);
		isRGB = (imp.getType()==ImagePlus.COLOR_RGB);
		return  flags;
	}

	@Override
	public void run(ImageProcessor ip) {
		try {
			image.close();
		} catch (Exception Ex){
			
		}
		ImageProcessor ipaux=ip.duplicate();
	 	
		if (!isFloat && !isRGB) 
			ipaux=ipaux.toFloat(0, null);
		
		pass++;
		int r = (sz-1)/2;
		Conv cnv=new Conv();
		
		GScaleSpace sp=new GScaleSpace(r);
		//GScaleSpace sp=new GScaleSpace(sigma);
		double sigma=sp.getSigma();
		float[] kernx=sp.computeLapNKernel2D(n);
		if (scnorm) {
			GScaleSpace.scnorm(kernx,sigma,n);
		}
			
		kernel=new float[1][];
				
		kernel[0]=kernx;
		
		if (debug ) {
			FloatProcessor fp3=new FloatProcessor(sz,sz, kernx);
			new ImagePlus("Lap_"+n,fp3).show();	 
		}
		
		long time=-System.nanoTime();	
		
		FloatProcessor fpaux= (FloatProcessor) ipaux;
	
		cnv.convolveFloat(fpaux, kernx, sz, sz);
 
		time+=System.nanoTime();
		time/=1000.0f;
		System.out.println("elapsed time: " + time +" us");
		fpaux.resetMinAndMax();
		image=new ImagePlus("Convolved_"+n,fpaux);
		image.show();		
	}
	
	/* Saves the current settings of the plugin for further use
     * 
     *
    * @param prefs
    */
   public static void savePreferences(Properties prefs) {
	   		prefs.put(GN, Integer.toString(n));
	   		prefs.put(LEN, Integer.toString(sz));
	   		prefs.put(ISSEP, Boolean.toString(sep));
	   		prefs.put(NORM, Boolean.toString(scnorm));
         // prefs.put(SIGMA, Float.toString(sigma));

   }
	 
	
	public float[] getKernel(int i) {
		return kernel[i];
	}

	@Override
	public int showDialog(ImagePlus imp, String command, PlugInFilterRunner pfr) {
		this.pfr = pfr;
		int r = (sz-1)/2;
		GenericDialog gd=new GenericDialog("Power of Laplacian " + version);
		gd.addNumericField("half width", r, 1);
		gd.addNumericField("power 2 x", n, 1);
		gd.addCheckbox("Show kernel", debug);
		gd.addCheckbox("Scale nomalize", scnorm);
		gd.addPreviewCheckbox(pfr);
		gd.addDialogListener(this);
		gd.showDialog();
		if (gd.wasCanceled())
			return DONE;
		 
		 
		return IJ.setupDialog(imp, flags);
	}

	// Called after modifications to the dialog. Returns true if valid input.
	public boolean dialogItemChanged(GenericDialog gd, AWTEvent e) {
 
		int r = (int)(gd.getNextNumber());
		n = (int)(gd.getNextNumber());
		debug = gd.getNextBoolean();
		scnorm = gd.getNextBoolean();
		sz = 2*r+1;
		if (gd.wasCanceled())
			return false;
		return r>0;
	}
	
	
	@Override
	public void setNPasses (int nPasses) {
		this.nPasses = nPasses;
	}

}
