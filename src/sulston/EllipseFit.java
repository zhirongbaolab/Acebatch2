package sulston;

import flanagan.math.Maximization;
import flanagan.math.MaximizationFunction;
import ij.ImagePlus;
import ij.Prefs;
import ij.gui.ImageCanvas;
import ij.gui.ImageWindow;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;

import java.awt.Color;
import java.awt.Dimension;
import java.io.File;
import java.text.DecimalFormat;

import org.rhwlab.dbaccess.EmbryoXML;
import org.rhwlab.utils.EUtils;

public class EllipseFit {

	String		iSeries;
	String 		iAnnots;
	String		iImages;
	int			iTime;
	int			iCenter;
	ImagePlus	iPlus;

	double		iZCenter;
	int			iXCenter;
	int			iYCenter;
	double		iMajor;
	double		iMinor;
	double		iAngle;
	double		iTheta;


	public int[] xCoordinates;
	public int[] yCoordinates;
	public int nCoordinates = 0;
	public double xCenter;
	public double  yCenter;
	public double major;
	public double minor;
	public double angle;
	public double theta;

	boolean record;

	public EllipseFit(String series, String imgDir, int testtime, int plane) {
		iTime = testtime;
		iCenter = plane;
		iSeries = series;
		iImages = imgDir;

	}


	public EllipseFit(EmbryoXML exml, int testTime, int plane) {
		iTime = testTime;
		iCenter = plane;
		iSeries = exml.iRecord[EmbryoXML.SERIES];
		iAnnots = exml.iRecord[EmbryoXML.ANNOTS];
		iImages = exml.iRecord[EmbryoXML.IMAGELOC];
		println("EllipseFit, " + iSeries + CS + iTime);
		//getImage();
		//maximize();
		//show();

	}

	int getImage() {
		int r = 0;
		String s = iImages + "/tif/";
		//println("getImage, " + s);
		String [] d = new File(s).list();
		if (d == null || d.length < 100) {
			println("Failure, EllipseFit could not open image for, " + iSeries);
			return 1;
		}
		String ss = d[0];
		ss = ss.substring(0, ss.length() - 13);
		String s2 = EUtils.makePaddedInt(iTime);
		ss += "-t" + s2;
		s2 = EUtils.makePaddedInt(iCenter, 2);
		ss += "-p" + s2 + ".tif";
		ss = s + ss;
		Prefs.open100Percent = true;
		ImagePlus iplus = new ImagePlus(ss);
		iPlus = iplus;
		if (iPlus.getProcessor() == null) {
			println("Failure, EllipseFit, bogus ImagePlus could not open image for, " + iSeries);
			return 1;

		}
		ImageProcessor imgproc = iPlus.getProcessor();
		Roi testRoi = new Roi(100, 100, 200, 200);
		Roi.setColor(Color.WHITE);
		iPlus.setRoi(testRoi);
		show();


		return 0;

	}

	void maximize() {
		//*
		Maximization max = new Maximization();
        MaximFunct funct = new MaximFunct(this);
        double [] sigs = {350, 250, 550, 350, 0};


        ImageProcessor proc = iPlus.getProcessor();
        int height = proc.getHeight();
        int width = proc.getWidth();
        sigs[0] = width/2;
        sigs[1] = height/2;
        sigs[2] = width * .8;
        sigs[3] = height * .8;

        double [] start = sigs;
        double [] step = {50, 50, 50, 50, 1};
        double ftol = 100;
        int nmax = 3000;
        max.addConstraint(0, -1, 200);
        max.addConstraint(1, -1, 100);
        max.addConstraint(2, -1, 200);
        max.addConstraint(3, -1, 100);
        max.addConstraint(0, 1, 700);
        max.addConstraint(1, 1, 500);
        max.addConstraint(2, 1, 700);
        max.addConstraint(3, 1, 500);
        max.nelderMead(funct, start, step, ftol, nmax);
        // get the minimum value
        double maximum = max.getMaximum();
        int niter = max.getNiter();
        //println("main, maximum=" + maximum + CS + niter);

        // get values of y and z at minimum
        StringBuffer sb = new StringBuffer(iSeries);
        double [] b = max.getParamValues();
        for (int i=0; i < b.length; i++) {
        	//println(i + CS + b[i]);
        	sb.append(C + fmt1(b[i]));
        }
        //println(sb.toString());
        int x = (int)b[0];
        int y = (int)b[1];
        int w = (int)b[2] - x;
        int h = (int)b[3] - y;
        //show(x,y,w,h);
		xCenter = b[0];
		yCenter = b[1];
		major = b[2];
		minor = b[3];
		angle = b[4];
		theta = 2 * Math.PI * angle / 360;
		iXCenter = (int)Math.round(xCenter);
		iYCenter = (int)Math.round(yCenter);
		iMajor = (int)Math.round(major);
		iMinor = (int)Math.round(minor);
		iAngle = angle;
		iTheta = theta;
		//*/
		ImageProcessor imgproc = iPlus.getProcessor();
		makeRoi(imgproc); // this is where the ellipse parms are used to find the coordinates
		PolygonRoi proi = new PolygonRoi(xCoordinates, yCoordinates, nCoordinates, Roi.POLYGON);
		//Roi testRoi = new Roi(100, 100, 200, 200);
		Roi.setColor(Color.WHITE);
		iPlus.setRoi(proi);


	}

	public void makeRoi(ImageProcessor ip) {
		record = true;
		int size = ip.getHeight()*3;
		size = Math.max(1500, size);
		xCoordinates = new int[size];
		yCoordinates = new int[size];
		nCoordinates = 0;
		drawEllipse(ip);
		//println("makeRoi,");
		record = false;
	}


	private double sqr(double x) {
		return x*x;
	}


	double [] testRoi(double x, double y, double maj, double min, double ang) {
		//double areaPenalty = 45;
		double areaPenalty = 45;
		xCenter = x;
		yCenter = y;
		major = maj;
		minor = min;
		angle = ang;
		theta = 2 * Math.PI * ang / 360;
		ImageProcessor imgproc = iPlus.getProcessor();
		///////imgproc.threshold(40);
		makeRoi(imgproc);
		PolygonRoi proi = new PolygonRoi(xCoordinates, yCoordinates, nCoordinates, Roi.POLYGON);
		iPlus.setRoi(proi);
		ImageStatistics imgstat = iPlus.getStatistics();
		double [] rtn = new double[4];
		rtn[0] = imgstat.area;
		rtn[1] = imgstat.area * imgstat.mean;
		rtn[2] = (imgstat.mean - areaPenalty) * imgstat.area;
		rtn[3] = imgstat.mean;

		return rtn;
	}

	/** Draws the ellipse on the specified image. */
	public void drawEllipse(ImageProcessor ip) {
		if (major==0.0 && minor==0.0)
			return;
		int xc = (int)Math.round(xCenter);
		int yc = (int)Math.round(yCenter);
		int maxY = ip.getHeight();
		int xmin, xmax;
		double sint, cost, rmajor2, rminor2, g11, g12, g22, k1, k2, k3;
		int x, xsave, ymin, ymax;
		int[] txmin = new int[maxY];
		int[] txmax = new int[maxY];
		double j1, j2, yr;

		sint = Math.sin(theta);
		cost = Math.cos(theta);
		rmajor2 = 1.0 / sqr(major/2);
		rminor2 = 1.0 / sqr(minor/2);
		g11 = rmajor2 * sqr(cost) + rminor2 * sqr(sint);
		g12 = (rmajor2 - rminor2) * sint * cost;
		g22 = rmajor2 * sqr(sint) + rminor2 * sqr(cost);
		k1 = -g12 / g11;
		k2 = (sqr(g12) - g11 * g22) / sqr(g11);
		k3 = 1.0 / g11;
		ymax = (int)Math.floor(Math.sqrt(Math.abs(k3 / k2)));
		if (ymax>maxY)
			ymax = maxY;
		if (ymax<1)
			ymax = 1;
		ymin = -ymax;
 		// Precalculation and use of symmetry speed things up
		for (int y=0; y<=ymax; y++) {
			//GetMinMax(y, aMinMax);
			j2 = Math.sqrt(k2 * sqr(y) + k3);
			j1 = k1 * y;
			txmin[y] = (int)Math.round(j1 - j2);
			txmax[y] = (int)Math.round(j1 + j2);
		}
		if (record) {
			xCoordinates[nCoordinates] = xc + txmin[ymax - 1];
			yCoordinates[nCoordinates] = yc + ymin;
			nCoordinates++;
		} else
			ip.moveTo(xc + txmin[ymax - 1], yc + ymin);
		for (int y=ymin; y<ymax; y++) {
			x = y<0?txmax[-y]:-txmin[y];
			if (record) {
				xCoordinates[nCoordinates] = xc + x;
				yCoordinates[nCoordinates] = yc + y;
				nCoordinates++;
			} else
				ip.lineTo(xc + x, yc + y);
		}
		for (int y=ymax; y>ymin; y--) {
			x = y<0?txmin[-y]:-txmax[y];
			if (record) {
				xCoordinates[nCoordinates] = xc + x;
				yCoordinates[nCoordinates] = yc + y;
				nCoordinates++;
			} else
				ip.lineTo(xc + x, yc + y);
		}
	}

	void show() {
		ImageProcessor iproc = iPlus.getProcessor();
		//makeRoi(iproc);
		//drawEllipse(iproc);
		ImageCanvas icanvas = new ImageCanvas(iPlus);
		icanvas.setMagnification(1.0);
		icanvas.zoomIn(10, 10);
		ImageWindow iwin = new ImageWindow(iPlus, icanvas);
		iwin.setSize(new Dimension(700, 500));
	}



	   class MaximFunct implements MaximizationFunction{
		   EllipseFit		iEllipseFit;
		   int				iCount;
		   public MaximFunct(EllipseFit ellipseFit) {
			   iEllipseFit = ellipseFit;

		   }

			public double function(double[] param) {
				double x = param[0];
				double y = param[1];
				double maj =  param[2];
				double min =  param[3];
				double ang =  param[4];
				//if (iCount++ % 10 == 0) println("function, a, " + x + CS + y + CS + maj + CS + min + CS + ang);
				double [] r = iEllipseFit.testRoi(x, y, maj, min, ang);
				//println("function, b, " + x + CS + y + CS + maj + CS + min + CS + ang);
				// TODO Auto-generated method stub
				return r[2];
			}

	   }

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

	private static void println(String s) {System.out.println(s);}
    private static void print(String s) {System.out.print(s);}
    private static final String CS = ", ", C = ",";
    private static final String TAB = "\t";
    private static final DecimalFormat DF0 = new DecimalFormat("####");
    private static final DecimalFormat DF1 = new DecimalFormat("####.#");
    private static final DecimalFormat DF4 = new DecimalFormat("####.####");
    private static String fmt4(double d) {return DF4.format(d);}
    private static String fmt1(double d) {return DF1.format(d);}
    private static String fmt0(double d) {return DF0.format(d);}


}
