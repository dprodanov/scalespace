
import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.gui.GenericDialog;
import ij.gui.Roi;
import ij.measure.Calibration;
import ij.gui.DialogListener;
import ij.plugin.filter.ExtendedPlugInFilter;
import ij.plugin.filter.PlugInFilterRunner;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ijaux.datatype.Pair;
import ijaux.scale.*;

import java.awt.*;
import java.io.File;
import java.util.*;
import java.util.List;

/*
import org.jfree.chart.ChartFactory;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import activeSegmentation.IFilter;
*/
import dsp.Conv;

/**
 * @version   	 28 Feb 2019
 * 				 Gaussian jet of order k
 *   
 * 
 * @author Dimiter Prodanov, IMEC
 *
 *
 * @contents
 * This plug-in convolves an image with a Gaussian derivative of order  (n,m).
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
public class Gaussian_Jet_ implements ExtendedPlugInFilter, DialogListener {

    @SuppressWarnings("unused")

    private PlugInFilterRunner pfr=null;

	final int flags=DOES_ALL+CONVERT_TO_FLOAT+SUPPORTS_MASKING+KEEP_PREVIEW;
	private String version="2.0";
	   @SuppressWarnings("unused")

	private int nPasses=1;
	    @SuppressWarnings("unused")

	   private int pass=0;
	public final static String SIGMA="LOG_sigma", LEN="G_len" ,MAX_LEN="G_MAX", 
			ISSEP="G_SEP", GN="G_Xn", GM="G_Yn", SCNORM="G_SCNORM";

	private static int sz = Prefs.getInt(LEN, 5);
	//private  int max_sz= Prefs.getInt(MAX_LEN, 8);
	private float[][] kernel=null;
	public static boolean debug= IJ.debugMode;

	private static int nn = Prefs.getInt(GN, 1);
 
 
	private static boolean scnorm=Prefs.getBoolean(SCNORM, false);
	private ImagePlus image=null;	
	private boolean isFloat=false;
	private boolean isRGB=false;
	// window size in sigmas
	private static int wnd=3;
	private boolean isEnabled=true;
 


	/* NEW VARIABLES*/

	/** A string key identifying this factory. */
	private final  String FILTER_KEY = "GAUSSIAN";

	/** The pretty name of the target detector. */
	private final String FILTER_NAME = "Gaussian Jet";

	private final int TYPE=1;
	
	private Map< String, String > settings= new HashMap<String, String>();

	private ImageStack imageStack;

	/*
	 * @param args - args[0] should point to the folder where the plugins are installed 
	 */
	public static void main(String[] args) {

		try {

			File f=new File(args[0]);

			if (f.exists() && f.isDirectory() ) {
				System.setProperty("plugins.dir", args[0]);
				
			} else {
				throw new IllegalArgumentException();
			}
		}
		catch (Exception ex) {
			IJ.log("plugins.dir misspecified\n");
			ex.printStackTrace();
		}
		
		new ImageJ();
		
		GScaleSpace sp=new GScaleSpace(13);
		
		int sz= sp.getSize();
		int n=3;
		float[][] kernel=new float[n][];
		kernel[0]=sp.gauss1D();
		for (int i=1;i<n; i++) {
				float[] kerny=sp.diffNGauss1D(i);
				GScaleSpace.flip(kerny);
				kernel[i]=kerny;
		}
		
		showKernels(sz, n, kernel);
	 

	}

	private static void showKernels(int sz, int n, float[][] kernel) {
		ImageStack isc=new ImageStack(sz, sz);
		for (int i=0; i <n; i++) { 
			FloatProcessor fpkern2=new FloatProcessor(sz,sz);		 
			float[] disp=GScaleSpace.joinXY(kernel, i, n-i-1);
			fpkern2.setPixels(disp);
			isc.addSlice("Gx [" +(i)+"] Gy["+(n-i-1)+"]", fpkern2);

		}
		new ImagePlus("kernel sep",isc).show();
	}

	public void initialseimageStack(ImageStack img){
		this.imageStack = img;
	}
	
	@Override
	public int setup(String arg, ImagePlus imp) {
		image=imp;
		isFloat= (imp.getType()==ImagePlus.GRAY32);
		isRGB = (imp.getType()==ImagePlus.COLOR_RGB);
		cal=image.getCalibration();

		return  flags;
	}

	private boolean doCalib = false;
	
	/*
	 * This variable is to calibrate the Image Window
	 */
	private Calibration cal=null;

	@Override
	public void run(ImageProcessor ip) {
		
		int r = (sz-1)/2;
		
		// filter window size
		if (wnd<0)
			wnd=-wnd;
		if (wnd>5)
			wnd=5;
		GScaleSpace sp=new GScaleSpace(r, wnd);
 		
		imageStack=filter(ip, sp,  scnorm, nn+1);
	
		image=new ImagePlus("Convolved_"+nn, imageStack);
		image.setCalibration(cal);
		image.updateAndDraw();
		image.show();
		
	}


	
	private ImageStack filter(ImageProcessor ip, GScaleSpace sp, final boolean scnorm, int n){

		ImageProcessor ipaux=ip.duplicate();
		if (!isFloat && !isRGB) 
			ipaux=ipaux.toFloat(0, null);
		pass++;

		/*
		 * define an array of derivatives
		 * then convolve with each combination
		 *  d_x(k) d_y(n-k), i<=j
		 * */
		double sigma=sp.getSigma();
		 
		kernel=new float[n][];
		kernel[0]=sp.gauss1D();
		for (int i=1;i<n; i++) {
				float[] kerny=sp.diffNGauss1D(i);
				GScaleSpace.flip(kerny);
			if (scnorm) {
				scnorm(kerny,sigma,i);
			}
			kernel[i]=kerny;
		}
		
		if (debug && pass==1) {
			showKernels(sz, n, kernel);
			
		}
		
		ImageStack is=new ImageStack(ip.getWidth(), ip.getHeight());
		
		long time=-System.nanoTime();	
		Conv cnv=new Conv();
		
		for (int i=0; i <n; i++) { 
			FloatProcessor fpaux= (FloatProcessor) ipaux.duplicate();
			cnv.convolveSep(fpaux, kernel[n-i-1], kernel[i]);
			fpaux.setSnapshotPixels(null);
			is.addSlice("Gx [" +(i)+"] Gy["+(n-i-1)+"]", fpaux);
	

		}
		time+=System.nanoTime();
		time/=1000.0f;		
		System.out.println("elapsed time: " + time +" us");

		return is;

	}

	private void scnorm(float[] kern, double sigma, int n) {
		sigma=Math.pow(sigma, n);
		for (int i=0; i<kern.length; i++) {
			kern[i]*=sigma;
		}
	}

	/* Saves the current settings of the plugin for further use
	 * 
	 *
	 * @param prefs
	 */
	public static void savePreferences(Properties prefs) {
		prefs.put(GN, Integer.toString(nn));
		prefs.put(LEN, Integer.toString(sz));
		prefs.put(SCNORM, Boolean.toString(scnorm));

	}


	public float[] getKernel(int i) {
		return kernel[i];
	}

	@Override
	public int showDialog(ImagePlus imp, String command, PlugInFilterRunner pfr) {
		this.pfr = pfr;
		int r = (sz-1)/2;
		GenericDialog gd=new GenericDialog("Gaussin Jet " + version);
		//gd.addNumericField("span x sigma", wnd, 3);
		gd.addNumericField("half width", r, 1);
		gd.addNumericField("order ", nn, 0);
		gd.addCheckbox("Show kernel", debug);

		gd.addCheckbox("Scale nomalize", scnorm);
		
		gd.addPreviewCheckbox(pfr);
		gd.addDialogListener(this);
		gd.setResizable(false);		
		
		if (cal!=null) {
			if (!cal.getUnit().equals("pixel"))
				gd.addCheckbox("units ( "+cal.getUnit() + " )", doCalib); 
		}	
		
		gd.showDialog();
		
		if (gd.wasCanceled())
			return DONE;

		return IJ.setupDialog(imp, flags);
	}

	// Called after modifications to the dialog. Returns true if valid input.
	public boolean dialogItemChanged(GenericDialog gd, AWTEvent e) {
		//wnd = (int)(gd.getNextNumber());
		int r = (int)(gd.getNextNumber());
		nn = (int)(gd.getNextNumber());
		debug = gd.getNextBoolean();

		scnorm = gd.getNextBoolean();
		
		if (cal!=null)
			doCalib=gd.getNextBoolean();
		
		sz = 2*r+1;
		if (gd.wasCanceled())
			return false;

		return r>0;
	}


	//@Override
	public boolean reset() {
		sz= Prefs.getInt(LEN, 2);
		//max_sz= Prefs.getInt(MAX_LEN, 8);
		scnorm=Prefs.getBoolean(SCNORM, false);
		nn = Prefs.getInt(GN, 1);
		return true;
	}

	//@Override
	public Map<String, String> getDefaultSettings() {

		settings.put(LEN, Integer.toString(sz));
		//settings.put(MAX_LEN, Integer.toString(max_sz));
		settings.put(SCNORM, Boolean.toString(scnorm));
		settings.put(GN, Integer.toString(nn));

		return this.settings;
	}

	//@Override
	public boolean updateSettings(Map<String, String> settingsMap) {
		sz=Integer.parseInt(settingsMap.get(LEN));
		//max_sz=Integer.parseInt(settingsMap.get(MAX_LEN));
		scnorm=Boolean.parseBoolean(settingsMap.get(SCNORM));
		nn=Integer.parseInt(settingsMap.get(GN));

		return true;
	}


	//@Override
	public String getKey() {
		return FILTER_KEY;
	}

	//@Override
	public String getName() {
		return FILTER_NAME;
	}

	/* (non-Javadoc)
	 * @see ij.plugin.filter.ExtendedPlugInFilter#setNPasses(int)
	 */
	@Override
	public void setNPasses (int nPasses) {
		this.nPasses = nPasses;
	}
	
	
	//@Override
	public boolean isEnabled() {
		return isEnabled;
	}

	//@Override
	public void setEnabled(boolean isEnabled) {
		this.isEnabled= isEnabled;
	}

	//@Override
	public int getFilterType() {
		return this.TYPE;
	}

	//@Override
	public <T> T getFeatures() {
		return null;
	}

	//@Override
	public Set<String> getFeatureNames() {
		return null;
	}

	

}
