import java.awt.AWTEvent;
import java.awt.Polygon;
import java.awt.image.IndexColorModel;
import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.ListIterator;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Properties;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;



import ij.CompositeImage;
import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.gui.Overlay;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.gui.Wand;
import ij.measure.Calibration;
import ij.plugin.filter.Convolver;
import ij.plugin.filter.ExtendedPlugInFilter;
import ij.plugin.filter.PlugInFilter;
import ij.plugin.filter.PlugInFilterRunner;
import ij.plugin.frame.RoiManager;
import ij.process.Blitter;
import ij.process.ByteProcessor;
import ij.process.ColorProcessor;
import ij.process.FloatProcessor;
import ij.process.FloodFiller;
import ij.process.ImageConverter;
import ij.process.ImageProcessor;
import ij.process.LUT;
import ij.process.ShortProcessor;
import ijaux.scale.SortUtil;

/*
 * @version 	2.9  26 May 2019
 * 				- added biasing for bright or dark components
 * 				2.2.8 23 Oct 2013
 * 				- added calibration to the output images
 * 				2.2.7 10 Sept 2013
 * 				- constant boundary condition incorporated in neighborhood retrieval
 * 				2.2.6 12 Jun 2013
 * 				- output combined in a stack
 * 				2.2.5 17 May 2013
 * 				- code refactoring and cleanup
 * 				- bug fixing
 * 				2.2 22 Mar 2013
 * 				- polygon ordering by using Wand
 * 				- bug fixing
 * 				2.1 16 Mar 2013
 * 				- polygon ordering by list sorting
 * 				2.0 21 Feb 2013
 * 				- namespace change
 * 				- polygon ordering by orientation
 * 				- code refactoring
 * 				1.5 10 Feb 2012
 * 				- contour labeling added
 * 				1.0	20 Dec 2012
 * 				- blobs labeling added
 *   
 * 
 * @author Dimiter Prodanov
 * 		  IMEC
 *
 *
 * @contents
 * This pluign  computes the image zero crossings.  
 * The plugin assumes a 2D float image 
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
public class ZeroCrosser_  implements PlugInFilter {

	public static final String KSZ="KSZ", SEGM="SG", LBL="LBL", BIAS="BIAS";
	final int flags=DOES_32;
	private String version="1.0";

	private static int scale=255;

	static int sz=Prefs.getInt(KSZ,3);
	public static boolean debug=IJ.debugMode;

	private Calibration cal;
	private ImagePlus outimg;
	private IndexColorModel lut=null;

	ConcurrentHashMap<Integer,int[]> vcoords=new  ConcurrentHashMap<Integer,int[]>();
	private static boolean segment=Prefs.getBoolean(SEGM,false);
	private static boolean label=Prefs.getBoolean(LBL,false);
	private static boolean boolbias=Prefs.getBoolean(BIAS,false);
	private static RoiManager roiman=RoiManager.getInstance();
	
  
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
		if (imp==null) return DONE;
		cal=imp.getCalibration();
		if (arg.equals("about")){
			showAbout();
			return DONE;
		}
		if(IJ.versionLessThan("1.47")|| !showDialog(imp)) {
			return DONE;
		}
		else {
			makeRedLut();
			return IJ.setupDialog(imp, flags);
		}
	}


	@Override
	public void run(ImageProcessor ip) {

		final int width=ip.getWidth();
		final int height=ip.getHeight();
		final int size=width*height;

		ByteProcessor bp = detectZeros(ip, bias);
		int x=-1;
		int y=-1;
		for (int i=0; i< size; i++) {
			if (bp.get(i)>0) {
				bp.set(i, scale);			
				if (segment) {
					x=i%width ; // IJ  quirk
					y=i/width ; // IJ quirk
					vcoords.put(i,new int[]{x,y});
					//System.out.append("["+i+" - "+"( " +x+" "+y+")]\n");
				}			
			}
		}

		bp.setColorModel(lut);		
		outimg =new ImagePlus("zero crossings",bp);		
		outimg.setCalibration(cal.copy());
		outimg.show();
		if (label)
			IJ.setColumnHeadings("object\tamount");
		/*
		 * labeling the blobs & contours
		 */
		if (label) {
			ImageStack is=new ImageStack(width, height);
			/*
			 *  labeling the contours
			 */
			ShortProcessor map2=new ShortProcessor (width, height);
			int numberOfLabels=doLabels8(bp, map2, 0);
			IndexColorModel alut2 = makeRandomLut(0);
			map2.setColorModel(alut2);
			is.addSlice("blobs 4c "+numberOfLabels, map2);
			//ImagePlus img=new ImagePlus("contours 8c "+numberOfLabels,map2);
			IJ.write("contours\t"+numberOfLabels);
	
			/*
			 *  labeling the blobs
			 */
			ShortProcessor map=new ShortProcessor(width, height);
			int numberOfLabels2=doLabels4(bp, map, scale);
			IndexColorModel alut = makeRandomLut(0);
			map.setColorModel(alut);
			//ImagePlus img2=new ImagePlus("blobs 4c "+numberOfLabels,map);
			//img2.show();
			IJ.write("blobs\t"+numberOfLabels2);

			is.addSlice("blobs 4c "+numberOfLabels2, map);

			ImagePlus img=new ImagePlus("segmented",is);
			img.setCalibration(cal.copy());
			
			/*
			 *  segmenting the polygons
			 */
			if(segment) {
				Wand wand= new  Wand(map);
				final Overlay overlay= new Overlay();

				if (roiman==null) {
					roiman=new RoiManager();

				}
				if (!roiman.isVisible())
					roiman.setVisible(true);

				//int k=19;
				//int k=9;
				//int k=8;
				//int k=6;
				for (int k=1; k<numberOfLabels; k++) {
					Polygon p= getContourWand(map,k, vcoords, wand);
					//Polygon p= getContourVect(map2,k, vcoords);
					if (p!=null) {
						PolygonRoi roi=new PolygonRoi(p, Roi.POLYGON);
						roiman.addRoi(roi);
						overlay.add(roi);
					}
				}

				img.setOverlay(overlay);
			}
		/*	CompositeImage  cimg=new CompositeImage(img, 1);
		 
			cimg.setChannelLut(new LUT(alut2, 0,numberOfLabels), 1);
			cimg.setChannelLut(new LUT(alut, 0,numberOfLabels2 ), 2);
			cimg.show();*/
			img.show();
		}
		ip.reset();
	}

	/*
	private Polygon getContourOrient(ShortProcessor map,int label, 
			ConcurrentHashMap<Integer,int[]> coords) {
		final int width=map.getWidth();
		final int height=map.getHeight();
		int sz=width*height;
		LinkedHashMap<Integer,int[]> cmap= new LinkedHashMap<Integer,int[]> ();
		int cnt=0;
		//System.out.print("<< " +label+" \n");
		double xm=0; double ym=0;
		for (int i=0; i<sz; i++ ) {
			int key=map.get(i); 
			if (key==label) {	
				int[] c=new int[2];		
				c[0]=i % width;
				c[1]=i/ width;				
				xm+=c[0];
				ym+=c[1];							 
				cmap.put(cnt, c);						
				cnt++;				
			} // end if

		} // end for
		//System.out.print("\n "+cnt+" >> \n");
		xm = xm/((double)cnt);
		ym = ym/((double)cnt);

		double[] angles=new double [cnt];

		Polygon poly=new Polygon();
		int k=0;
		for (Entry<Integer, int[]> e:cmap.entrySet() ) {
			final int[] c=e.getValue();
			angles[k]=Math.atan2((double)c[1]-ym, (double)c[0]-xm)+ Math.PI;
			k++;		
		}
		//int[] index=quicksort_i(angles);
		int[] index=SortUtil.quicksort(angles, false);
		//System.out.println();
		for (int i=0; i<angles.length; i++) {
			//System.out.println( "("+i+" -> "+angles[index[i]]+"/ " +index[i]+"),");
			final int[] c=cmap.get(index[i]);
			poly.addPoint(c[0], c[1]);
		}


		return poly;
	}
*/
	private Polygon getContourWand(ShortProcessor map,int label, 
			ConcurrentHashMap<Integer,int[]> coords, Wand wand) {
		final int width=map.getWidth();
		final int height=map.getHeight();
		int sz=width*height;
		LinkedHashMap<Integer,ArrayList<int[]>> cmap= new LinkedHashMap<Integer,ArrayList<int[]>> (100);

		log("<< " +label+" \n");

		coords.entrySet();
		for (int i=0; i<sz; i++ ) {
			int key=map.get(i); 
			if (key==label) {
				distribute(i,   width,  cmap);

			} // end if

		} // end for
		//System.out.print("\n "+cnt+" >> \n");
		Set<Entry<Integer, ArrayList<int[]>>> entries=cmap.entrySet();
		//Iterator<Entry<Integer, ArrayList<int[]>>> iter=entries.iterator();
		int key=0;
		boolean first=false;
		for (Entry<Integer, ArrayList<int[]>> e:entries) {
			//ArrayList<int[]> value=e.getValue();
			//System.out.print(e.getKey() + ">>");
			if (!first) {
				key=e.getKey();
				break;
			}
			first=true;

			/*	for (int[] k: value) {
				System.out.print("("+k[0]+" "+ k[1]+"),");
			}
			System.out.print(" <<\n");
			 */
		}
		Polygon poly=new Polygon();


		try {
			ArrayList<int[]> aux=cmap.get(key);
			if (aux!=null) {
				int[] c=aux.get(0);
				int startX=c[0];
				int startY=c[1];
				log(key +" -> ("+startX+" "+ startY+")");
				wand.autoOutline(startX, startY, label, label+1);		
				poly.xpoints=wand.xpoints;
				poly.ypoints=wand.ypoints;
				poly.npoints=wand.npoints;
				return poly;
			}  
		} catch (NullPointerException e) {
			e.printStackTrace();
			return null;
		}

		return null;
	}

	final int maxd=2; // maximal metrical radius of a "circle"

	
	private void distribute(int idx, int width,  LinkedHashMap<Integer,ArrayList<int[]>> cmap) {
		int[] c=new int[2];		
		c[0]=idx % width;
		c[1]=idx / width;

		for (Entry<Integer, ArrayList<int[]>> e:cmap.entrySet()) {
			ArrayList<int[]> clist=e.getValue();
			final int lastind=clist.size()-1;
			final int[] ctail=clist.get(lastind);
			final int[] chead=clist.get(0);
			int dist= dist (ctail,c); 
			if (dist<maxd ) {
				clist.add(c);
				return;
			}
			dist= dist (chead,c); 
			if (dist<maxd ) {
				clist.add(0,c);
				return;
			}

		}
		ArrayList<int[]> alist=new ArrayList<int[]>();
		alist.add(c);
		cmap.put(idx, alist);
	}

	private int dist(int[] u, int[] v) {
		final int d= Math.max(Math.abs(u[0]-v[0]) , Math.abs(u[1]-v[1]));
		//final int d=  Math.abs(u[0]-v[0]) + Math.abs(u[1]-v[1]);
		return d;
	}

	@SuppressWarnings("unchecked")
	public <K> ArrayList<K> flip(ArrayList<K> arr) {
		ArrayList<K> ret= new ArrayList<K>();
		int c=arr.size()-1;
		if (c==1)
			return (ArrayList<K>) arr.clone();
		for (int i=c; i>=0; i--) {
			ret.add(arr.get(i));
		}
		return ret;
	}

	/*	private ArrayList<int[]>  join(LinkedHashMap<Integer,ArrayList<int[]>> cmap) {

		ArrayList<ArrayList<int[]>> pointlist=new ArrayList<ArrayList<int[]>> ();
		pointlist.addAll(cmap.values());*/
	private ArrayList<int[]>  join(ArrayList<ArrayList<int[]>> pointlist) {
		int lsize=pointlist.size();
		final int fsize=lsize;//*(lsize-1)/2;
		log("list size "+lsize);
		@SuppressWarnings("unchecked")
		ArrayList<int[]> ret=(ArrayList<int[]>) pointlist.get(0).clone();
		pointlist.remove(0);

		int[] chead=ret.get(0);
		int lastind=ret.size()-1;
		int[] ctail=ret.get(lastind);

		// we will append to the first list

		int cnt=0;
		int kind=0;
		boolean isBlob=false;
		while(!pointlist.isEmpty() && cnt<fsize && kind<lsize && !isBlob) {
			//	try {
			// getting the k-th element, presumably the first
			ArrayList<int[]> alist = pointlist.get(kind);
			int[] head=alist.get(0);
			lastind=alist.size()-1;
			int[] tail=alist.get(lastind);
			int d1=dist(chead, head);
			int d2=dist(chead, tail);
			int d3=dist(ctail, tail);
			int d4=dist(ctail, head);
			int d5=dist(chead,ctail);
			isBlob=d5==1;
			//System.out.println("d1 "+d1+", d2 "+d2+", d3 "+d3+", d4 "+d4 +" >");
			if (d1<=maxd || (d2<=maxd) || d3<=maxd || d4<=maxd) {

				if (d1<=maxd) { // append to the head, flipped, pop
					ret.addAll(0, flip(alist));
					log("joining d1: t->h+ch->ct");	
				} else {
					if (d2<=maxd) {  // append to the head,pop
						ret.addAll(0, alist);
						log("joining d2: h->t+ch->ct");		
					} else {
						if (d3<=maxd) { // append to the tail, flipped
							ret.addAll(flip(alist));
							log("joining d3: ch->ct+t->h");	
						} else {
							if (d4<=maxd) { // append to the head
								ret.addAll(alist);	
								log("joining d4: ch->ct+h->t");	
							} // d4
						} // d3
					} //d2
				} //d1
				chead=ret.get(0);
				lastind=ret.size()-1;
				ctail=ret.get(lastind);	
				pointlist.remove(kind);

				kind=0; // we reset index
				//cnt--;
				System.out.println();
				for (int[] cc: ret) {
					log("(" +cc[0]+"," +cc[1]+")");
				}
				System.out.println();
			} else {
				log("<"+kind+" : [ d1: "+d1+", d2: "+d2+", d3: "+d3+", d4: "+d4 +"]" +
						"chead "+ chead[0]+ "," +chead[1]+"; ctail "+ctail[0]+ "," +ctail[1]+" //" +
						"head " +head[0] + "," +head[1]+"; tail "+tail[0]+ "," +tail[1]+ 		
						" >");
				cnt++;
				kind++;
				/*	System.out.println("tryting: ");
					for (int[] cc: alist) {
						System.out.print("(" +cc[0]+"," +cc[1]+")");
					}
					System.out.println();*/


			}				
			lsize=pointlist.size();
			log("list size " +lsize+" blob "+isBlob);
			/*			if (lsize>1)
				kind=kind % lsize;*/
			//			} catch (Exception e) {
			//				// TODO Auto-generated catch block
			//				e.printStackTrace();
			//			}
		}// end while	

		/*	ArrayList<int[]> aux=null;
		if (isBlob && !pointlist.isEmpty() ) {
			aux=join(pointlist);
		}
		if (aux!=null){
			ret.addAll(aux);
		}*/

		return ret;
	}

	private void log(String s) {
		if (debug)
			System.out.print(s);
	}
	
	private Polygon getContourVect(ShortProcessor map,int label, 
			ConcurrentHashMap<Integer,int[]> coords) {
		final int width=map.getWidth();
		final int height=map.getHeight();
		int sz=width*height;

		LinkedHashMap<Integer,ArrayList<int[]>> umap = 
				new LinkedHashMap<Integer,ArrayList<int[]>> ();

		log("<< " +label+" \n");

		for (int i=0; i<sz; i++ ) {
			int key=map.get(i); 
			if (key==label) {	
				distribute(i,   width,  umap);
			} // end if

		} // end for
		//System.out.print("\n "+cnt+" >> \n");

		for (Entry<Integer, ArrayList<int[]>> e:umap.entrySet()) {
			ArrayList<int[]> value=e.getValue();
			log(e.getKey() + ">>");
			for (int[] k: value) {
				log("("+k[0]+" "+ k[1]+"),");
			}
			log(" <<\n");
		}

		Polygon poly=new Polygon();

		log("\n");

		ArrayList<ArrayList<int[]>> pointlist=new ArrayList<ArrayList<int[]>> ();
		pointlist.addAll(umap.values());
		ArrayList<int[]> roilist=join(pointlist);
		pointlist.clear();
		log( ">>");
		Iterator<int[]> iter=roilist.iterator();
		int[] first=iter.next();
		poly.addPoint(first[0], first[1]);
		log("("+first[0]+" "+ first[1]+"),");
		while (iter.hasNext() ) {
			int[] k=iter.next();
			poly.addPoint(k[0], k[1]);

			log("("+k[0]+" "+ k[1]+"),");
		}
		poly.addPoint(first[0], first[1]);
		System.out.print("\r\n <<");

		return poly;
	}

/*	public   int[] quicksort_i(double[] main) {
		int[] index = new int[main.length];
		for (int i=0; i<main.length; i++) {
			index[i]=i;
		}	 
		quicksort_i(main, index, 0, index.length - 1); 
		return index;
	}
	// quicksort a[left] to a[right]
	public static void quicksort_i(double[] a, int[] index, int left, int right) {
		if (right <= left) return;
		int i = partition_i(a, index, left, right);
		quicksort_i(a, index, left, i-1);
		quicksort_i(a, index, i+1, right);
	}

	// partition a[left] to a[right], assumes left < right
	private static int partition_i(double[] a, int[] index, 
			int left, int right) {
		int i = left - 1;
		int j = right;
		while (true) {
			while ( a[index[++i]] < a[index[right]])      // find item on left to swap
				;                               // a[right] acts as sentinel
			while ( a[index[right]] < a[index[--j]])      // find item on right to swap
				if (j == left) break;           // don't go out-of-bounds
			if (i >= j) break;                  // check if pointers cross
			exch_i(index, i, j);               // swap two elements into place
		}
		exch_i(index, i, right);               // swap with partition element
		return i;
	}

	// exchange a[i] and a[j]
	private static void exch_i( int[] index, int i, int j) {
		int b = index[i];
		index[i] = index[j];
		index[j] = b;
	}
*/
	/**
	 * @return
	 */
	private IndexColorModel makeRandomLut(int ubgcol) {
		boolean acceptColor = false;
		byte [] reds=new byte[256];
		byte [] greens=new byte[256];
		byte [] blues=new byte[256];

		int cnt=0;	    
		ubgcol= ubgcol & 0x000000FF;
		int bgcol= (ubgcol& 0x000000ff) | (ubgcol  >> 8) | (ubgcol  >> 16);

		log( "background " +Integer .toHexString(bgcol) +" ( "+ bgcol+ " )");

		while (cnt <256) {
			final int col = (int)(Math. random() * 16777216);
			acceptColor = !((col& 0x000000ff) < 64) &&
					(((col & 0x0000ff00) >> 8) < 64) &&
					(((col & 0x00ff0000) >> 16) < 200) ;
			if(acceptColor && col!= bgcol)	{	    			
				blues[cnt]= (byte) (col& 0x0000FF);
				greens[cnt]= (byte) ((col & 0x00FF00) >>8);
				reds[cnt]= (byte) ((col & 0xFF0000)>>16);	    
				cnt ++;
			}    
			blues[0]= (byte) ubgcol;
			greens[0]= (byte) ubgcol;
			reds[0]= (byte) ubgcol;
		}
		IndexColorModel alut=new IndexColorModel(8, 256,reds, greens, blues);
		return alut;
	}

	/** Code based on
	 * http://homepages.inf.ed.ac.uk/rbf/HIPR2/labeldemo.htm
	 */
	private int doLabels4(ImageProcessor bp, ShortProcessor map, int bgcol) {

		if (bp instanceof ColorProcessor)
			return -1;

		final int width=bp.getWidth();
		final int height=bp.getHeight();
		final int size=width*height;

		final int mwidth=map.getWidth();
		final int mheight=map.getHeight();

		if (width!=mwidth || height!= mheight)
			throw new IllegalArgumentException ("dimensions mismatch ");
		
		this.bgcol=bgcol;
		/*if (edgecorr) {
			for (int a=0;a<mwidth; a++) {
				bp.set(a, 0, bgcol);
			}

			for (int a=0;a<mheight; a++) {
				bp.set(0, a, bgcol);
			}
		}*/

		int [] labels  = new int[size/2];

		for (int i=0; i<labels.length; i++) {
			labels[i]=i ; // ramp
		}


		int[] nbs= new int[2];
		int[] nbslab= new int[2];

		int numberOfLabels =1;
		int labelColour=1; // background

		int result=0;

		for(int y=0; y<height; y++) {	 
			//labelColour=0;
			for(int x=0; x<width; x++){		
				final int val=bp.get(x, y);				
				if( val == bgcol ){
					result = 0;  //nothing here
				} else {

					//The 4 connected visited neighbours
					//neighborhood4(bp, nbs, x, y, width);
					//neighborhood4(map, nbslab, x, y, width);
					neighborhood4f(bp, nbs, x, y, width, height);
					neighborhood4f(map, nbslab, x, y, width, height);
					//label the point
					if( (nbs[0] == nbs[1]) && (nbs[0] == bgcol )) { 

						// all neighbours are 0 so gives this point a new label
						result = labelColour;
						labelColour++;
					} else { //one or more neighbours have already got labels

						int count = 0;
						int found = -1;
						for( int j=0; j<nbs.length; j++){
							if( nbs[ j ] != bgcol ){
								count +=1;
								found = j;
							}
						}
						if( count == 1 ) {
							// only one neighbour has a label, so assign the same label to this.
							result = nbslab[ found ];
						} else {
							// more than 1 neighbour has a label
							result = nbslab[ found ];
							// Equivalence the connected points
							for( int j=0; j<nbslab.length; j++){
								if( ( nbslab[ j ] != 0 ) && (nbslab[ j ] != result ) ){
									associate(labels, nbslab[ j ], result );
								} // end if
							} // end for
						} // end else

					} // end else
					map.set(x, y, result);
				} // end if			    
			} // end for
		} // end for
		//reduce labels ie 76=23=22=3 -> 76=3
		//done in reverse order to preserve sorting
		log(" labels " + labelColour);
		for( int i= labels.length -1; i > 0; i-- ){
			labels[ i ] = reduce(labels, i );
		}

		/*now labels will look something like 1=1 2=2 3=2 4=2 5=5.. 76=5 77=5
			      this needs to be condensed down again, so that there is no wasted
			      space eg in the above, the labels 3 and 4 are not used instead it jumps
			      to 5.
		 */
		if (labelColour>0) {
			int condensed[] = new int[ labelColour ]; // can't be more than nextlabel labels

			int count = 0;
			for (int i=0; i< labelColour; i++){
				if( i == labels[ i ] ) 
					condensed[ i ] = count++;
			}

			/*for( int i= condensed.length -1; i > 0; i-- ){
			System.out.println(" l " + i+ " "+condensed[ i ]);
		}*/
			numberOfLabels = count;

			// now run back through our preliminary results, replacing the raw label
			// with the reduced and condensed one, and do the scaling and offsets too
			for (int i=0; i< size; i++){
				int val=map.get(i);
				val = condensed[labels[val]];
				map.set(i, val);
			}
			return numberOfLabels;
		} else {
			return -1;
		}
	}

	/** Code based on
	 * http://homepages.inf.ed.ac.uk/rbf/HIPR2/labeldemo.htm
	 */
	/*
	private int doLabels8(ImageProcessor bp, ShortProcessor map, int bgcol, boolean edgecorr) {

		if (bp instanceof ColorProcessor)
			return -1;

		final int width=bp.getWidth();
		final int height=bp.getHeight();
		final int size=width*height;

		final int mwidth=map.getWidth();
		final int mheight=map.getHeight();

		if (width!=mwidth || height!= mheight)
			throw new IllegalArgumentException ("dimensions mismatch ");

		if (edgecorr) {
			for (int a=0;a<mwidth; a++) {
				bp.set(a, 0, bgcol);
			}

			for (int a=0;a<mheight; a++) {
				bp.set(0, a, bgcol);
			}
		}

		int [] labels  = new int[size/2];

		for (int i=0; i<labels.length; i++) {
			labels[i]=i ; // ramp
		}


		int[] nbs= new int[4];
		int[] nbslab= new int[4];

		int numberOfLabels =1;
		int labelColour=1; // background

		int result=0;

		for(int y=0; y<height; y++) {	 
			//labelColour=0;
			for(int x=0; x<width; x++){		
				final int val=bp.get(x, y);				
			      if( val == bgcol ){
			    	  result = 0;  //nothing here
			      } else {

			    	  //The 8-connected visited neighbours
			    	  neighborhood8(bp, nbs, x, y, width);
			    	  neighborhood8(map, nbslab, x, y, width);

			    	  //label the point
			    	  if( (nbs[0] == nbs[1]) && (nbs[1] == nbs[2])  && (nbs[2] == nbs[3]) && (nbs[0] == bgcol )) { 
			    		  // all neighbours are 0 so gives this point a new label
			    		  result = labelColour;
			    		  labelColour++;
			    	  } else { //one or more neighbours have already got labels

			    		  int count = 0;
			    		  int found = -1;
			    		  for( int j=0; j<nbs.length; j++){
			    			  if( nbs[ j ] != bgcol ){
			    				  count +=1;
			    				  found = j;
			    			  }
			    		  }
			    		  if( count == 1 ) {
			    			  // only one neighbour has a label, so assign the same label to this.
			    			  result = nbslab[ found ];
			    		  } else {
			    			  // more than 1 neighbour has a label
			    			  result = nbslab[ found ];
			    			  // Equivalence of the connected points
			    			  for( int j=0; j<nbslab.length; j++){
			    				  if( ( nbslab[ j ] != 0 ) && ( nbslab[ j ] != 1 )&& (nbslab[ j ] != result ) ){
			    					  associate(labels, nbslab[ j ], result );
			    				  } // end if
			    			  } // end for
			    		  } // end else

			    	  } // end else
			    	  map.set(x, y, result);
			      } // end if			    
			} // end for
		} // end for
		//reduce labels ie 76=23=22=3 -> 76=3
		//done in reverse order to preserve sorting
		for( int i= labels.length -1; i > 0; i-- ){
			labels[ i ] = reduce(labels, i );
		}

		now labels will look something like 1=1 2=2 3=2 4=2 5=5.. 76=5 77=5
			      this needs to be condensed down again, so that there is no wasted
			      space eg in the above, the labels 3 and 4 are not used instead it jumps
			      to 5.

		if (labelColour>0) {
		int condensed[] = new int[ labelColour ]; // can't be more than nextlabel labels

		int count = 0;
		for (int i=0; i< labelColour; i++){
			if( i == labels[ i ] ) 
				condensed[ i ] = count++;
		}
		numberOfLabels = count;

		// now run back through our preliminary results, replacing the raw label
		// with the reduced and condensed one, and do the scaling and offsets too
	    for (int i=0; i< size; i++){
	    	int val=map.get(i);
	    	val = condensed[ labels[ val ] ];	        
	    	map.set(i, val);

	    }
			return numberOfLabels;
		} else {
			return -1;
		}
	}
	 */
	private int doLabels8(ImageProcessor bp, ShortProcessor map, int bgcol) {

		if (bp instanceof ColorProcessor)
			return -1;

		final int width=bp.getWidth();
		final int height=bp.getHeight();
		final int size=width*height;

		final int mwidth=map.getWidth();
		final int mheight=map.getHeight();

		if (width!=mwidth || height!= mheight)
			throw new IllegalArgumentException ("dimensions mismatch ");
		
		this.bgcol=bgcol;
	/*	if (edgecorr) {
			for (int a=0;a<mwidth; a++) {
				bp.set(a, 0, bgcol);
			}

			for (int a=0;a<mheight; a++) {
				bp.set(0, a, bgcol);
			}
		}*/

		int [] labels  = new int[size/2];

		for (int i=0; i<labels.length; i++) {
			labels[i]=i ; // ramp
		}


		int[] nbs= new int[4];
		int[] nbslab= new int[4];

		int numberOfLabels =1;
		int labelColour=1; // background

		int result=labelColour;

		for(int y=0; y<height; y++) {	 
			//labelColour=0;
			for(int x=0; x<width; x++){		
				final int val=bp.get(x, y);				
				if( val == bgcol ){
					result = 0;  //nothing here
				} else {

					//The 8-connected visited neighbours
					//neighborhood8(bp, nbs, x, y, width);
					//neighborhood8(map, nbslab, x, y, width);
					neighborhood8f(bp, nbs, x, y, width, height);
					neighborhood8f(map, nbslab, x, y, width, height);
					
					//label the point
					if( (nbs[0] == nbs[1]) && (nbs[1] == nbs[2])  && (nbs[2] == nbs[3]) && (nbs[0] == bgcol )) { 
						// all neighbours are 0 so gives this point a new label
						result = labelColour;
						labelColour++;
					} else { //one or more neighbours have already got labels

						int count = 0;
						int found = -1;
						for( int j=0; j<nbs.length; j++){
							if( nbs[ j ] != bgcol ){
								count +=1;
								found = j;
							}
						}
						if( count == 1 ) {
							// only one neighbour has a label, so assign the same label to this.
							result = nbslab[ found ];
						} else {
							// more than 1 neighbour has a label
							result = nbslab[ found ];
							// Equivalence of the connected points
							for( int j=0; j<nbslab.length; j++){
								if( ( nbslab[ j ] != 0 ) && ( nbslab[ j ] != 0 )&& (nbslab[ j ] != result ) ){
									associate(labels, nbslab[ j ], result );
								} // end if
							} // end for
						} // end else

					} // end else
					map.set(x, y, result);
				} // end if			    
			} // end for
		} // end for
		//reduce labels ie 76=23=22=3 -> 76=3
		//done in reverse order to preserve sorting
		for( int i= labels.length -1; i > 0; i-- ){
			labels[ i ] = reduce(labels, i );
		}

		/*now labels will look something like 1=1 2=2 3=2 4=2 5=5.. 76=5 77=5
			      this needs to be condensed down again, so that there is no wasted
			      space eg in the above, the labels 3 and 4 are not used instead it jumps
			      to 5.
		 */
		if (labelColour>0) {
			int condensed[] = new int[ labelColour ]; // can't be more than nextlabel labels

			int count = 0;
			for (int i=0; i< labelColour; i++){
				if( i == labels[ i ] ) 
					condensed[ i ] = count++;
			}
			numberOfLabels = count;

			// now run back through our preliminary results, replacing the raw label
			// with the reduced and condensed one, and do the scaling and offsets too
			for (int i=0; i< size; i++){
				int val=map.get(i);
				val = condensed[ labels[ val ] ];	    
				//val =  labels[ val ] ;	
				map.set(i, val);

			}
			return numberOfLabels;
		} else {
			return -1;
		}
	}
	/**
	 * @param bp
	 * @param nbs
	 * @param x
	 * @param y
	 */
/*	private void neighborhood4(ImageProcessor bp, int[] nbs, int x, int y, int width) {
		if ( x <= 0 ) x=1;
		if ( x >= width ) x=width-1;
		if ( y <= 0 ) y=1;
		nbs[0]=bp.get(x-1,y); // west
		nbs[1]=bp.get(x,y-1); // south
	}*/
	
	/**
	 * @param bp - image
	 * @param nbs - neighborhood array
	 * @param x - coordinate
	 * @param y - coordinate
	 */
	private void neighborhood4f(ImageProcessor bp, int[] nbs, int x, int y, int width, int height) {
	 
		if (x-1<0 || x-1>=width) {
			nbs[0]=bgcol;
		} else {
			nbs[0]=bp.get(x-1,y); // west
		}
		if (y-1<0 || y-1>=height) {
			nbs[1]=bgcol;
		} else {
			nbs[1]=bp.get(x,y-1); // north
		}
		
	}
	
	int bgcol=0;
	private int bias;
	 

	/**
	 * @param bp
	 * @param nbs
	 * @param x
	 * @param y
	 */
/*	private void neighborhood8(ImageProcessor bp, int[] nbs, int x, int y, int width) {
		if ( x <= 0 ) x=1;
		if ( x >= width ) x=width-1;
		if ( y <= 0 ) y=1;
		nbs[0]=bp.get(x-1,y); // W
		nbs[1]=bp.get(x,y-1); // N
		nbs[2]=bp.get(x-1,y-1); // W
		nbs[3]=bp.get(x+1,y-1); // NW

	}*/
	
	
	/**
	 * @param bp
	 * @param nbs
	 * @param x
	 * @param y
	 */
	private void neighborhood8f(ImageProcessor bp, int[] nbs, int x, int y, int width, int height) {
		
		if (x-1<0 || x-1>=width) {
			nbs[0]=bgcol;
		} else {
			nbs[0]=bp.get(x-1,y); // west
		}
		
		if (y-1<0 || y-1>=height) {
			nbs[1]=bgcol;
		} else {
			nbs[1]=bp.get(x,y-1); // north
		}
		
		if (x-1<0 || x+1>=width || y-1<0 || y-1>=height) {
			nbs[2]=bgcol;
		} else {
			nbs[2]=bp.get(x-1,y-1); // north-west
		}
		
		if (x<0 || x+1>=width || y-1<0 || y-1>=height) {
			nbs[3]=bgcol;
		} else {
			nbs[3]=bp.get(x+1,y-1); // north-east
		}

	}
	/**
	 * Associate(equivalence) a with b.
	 *  a should be less than b to give some ordering (sorting)
	 * if b is already associated with some other value, then propagate
	 * down the list.
	 */
	private void associate(int[] labels, int a, int b ) {	    
		if( a > b ) {
			associate(labels, b, a );
			return;
		}
		if( ( a == b ) || ( labels[ b ] == a ) ) return;
		if( labels[ b ] == b ) {
			labels[ b ] = a;
		} else {
			associate(labels, labels[ b ], a );
			if (labels[ b ] > a) {             //***rbf new
				labels[ b ] = a;
			}
		}
	}
	
	/**
	 * Reduces the number of labels.
	 */
	private int reduce(int[] labels, int a ){

		if( labels[ a ] == a ){
			return a;
		} else {
			return reduce(labels, labels[ a ] );
		}
	}
	
	/**
	 * @param ip
	 * @param width
	 * @param height
	 * @return
	 */
	private ByteProcessor detectZeros(ImageProcessor ip, int bias) {

		final int width=ip.getWidth();
		final int height=ip.getHeight();
		ByteProcessor bp=new ByteProcessor(width, height);

		for (int y=0; y<height; y++) {
			for (int x=0; x<width-1; x++) {
				final float val=ip.getf(x, y);
				boolean zero=(val*ip.getf(x+1, y)<0 );
				//if (bias<=-1) zero=zero && (val<0);
				//else if (bias>=1) zero=zero && (val>=0);
				//if bias is between (-1, 1) then the old behavior remains
				if (zero) {
					int u=bp.get(x+1, y);
					bp.set(x+1, y, u+1);
				}
			}
		}

		for (int y=0; y<height-1; y++) {
			for (int x=0; x<width; x++) {	
				final float val=ip.getf(x, y);
				boolean zero=(val*ip.getf(x, y+1)<0 );
				//if (bias<=-1) zero=zero && (val<0);
				//else if (bias>=1) zero=zero && (val>=0);
				//if bias is between (-1, 1) then the old behavior remains
				if (zero) {
					int u=bp.get(x, y+1);
					bp.set(x, y+1, u+1);
				}
			}
		}
		return bp;
	}

	public ImagePlus getResult() {
		return outimg;
	}




	public boolean showDialog(ImagePlus imp) {

		GenericDialog gd=new GenericDialog("Input Parameters");

		gd.addNumericField("scale to", scale, 0);
		gd.addCheckbox("label", label);
		gd.addCheckbox("segment ROIs", segment);
		
		//gd.addCheckbox("bright bias", boolbias);
		gd.addCheckbox("debug", debug);

		gd.showDialog();
		if (gd.wasCanceled())
			return false;		


		scale=(int) gd.getNextNumber();
		label= gd.getNextBoolean();
		segment= gd.getNextBoolean();
		//boolbias= gd.getNextBoolean();
		debug = gd.getNextBoolean();

		if (boolbias) bias=1; else bias=-1;
		return true;
	}

	void showAbout() {
		IJ.showMessage("Zero Crosser "+version,
				"The plugin calculates zero crossings in 2D"
				);
	}
	/**
	 * 
	 */
	private void makeRedLut() {
		byte[] reds=new byte[256];
		byte[] greens=new byte[256];
		byte[] blues=new byte[256];

		for (int i=1; i<reds.length; i++) {
			if (i>=scale)
				reds[i]=-1;
		}

		lut=new IndexColorModel(8,256, reds, greens, blues);
	}




	/* Saves the current settings of the plugin for further use
	 * 
	 *
	 * @param prefs
	 */
	public static void savePreferences(Properties prefs) {

		prefs.put(KSZ, Integer.toString(sz));
		prefs.put(SEGM, Boolean.toString(segment));
		prefs.put(LBL, Boolean.toString(label));
	}

}
