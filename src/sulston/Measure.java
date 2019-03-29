package sulston;

import analysis.Embryo;
import flanagan.analysis.Regression;
import flanagan.analysis.Stat;
import flanagan.math.Minimisation;
import flanagan.plot.PlotGraph;
import org.rhwlab.image.ParsingLogic.ImageNameLogic;
import org.rhwlab.snight.*;
import sulston.EmbryoFit.MinimFunct;

import java.io.File;
import java.text.DecimalFormat;
import java.util.Vector;

public class Measure {

	/** vars */
	private Config configManager;
	private NucleiMgr nucManager;
    public MeasureCSV iMeasureCSV;

    private String iSeries;

    private String iAnnots;
    private double iTimeSlope;
    private double iTimeOffset;
    private double iEXCenter;
    private double iEYCenter;
    private double iEMajor;
    private double iEMinor;
    private double iEAngle;
    private double iZSlope;
    private double iZCenter;

    private int iMaxPlane;

    private int iTime;
    private int iCenter;
    private Embryo iEmbryo;

    // is the plotting ever used? Can't remove a lot of code here if not
    public boolean iShowDetails;

	public Measure(Config configManager, NucleiMgr nucManager) {
        this.configManager = configManager;
        this.nucManager = nucManager;

        // build the MeasureCSV object
        this.iMeasureCSV = this.configManager.getNucleiConfig().getMeasureCSV();

        this.iSeries = getSeries();
        this.iAnnots = getAnnots();
        this.iMeasureCSV.put(MeasureCSV.defaultAtt_v1[MeasureCSV.SERIES_v1], iSeries);


        this.iEmbryo = new Embryo(this.configManager, this.nucManager, this.iSeries);

        this.iShowDetails = false;
	}

    /**
     * Called externally (after the class has been built) to execute the program
     */
	public void run() {
        long start = System.currentTimeMillis();

        writeCSV("");

        int r = getTimeScale();
        //println("getTimeScale, " + meas.iMeasureCSV);
        if (r != 0) return;
        r = processNucz();
        //println("processNucz, " + meas.iMeasureCSV);
        if (r != 0) return;
        r = fitEllipse();
        if (r != 0) {
            //println("Measure.main, " + "fitEllipse failure");
            return;
        }

        writeCSV("");

        println(this.iMeasureCSV.toString());
        long last = System.currentTimeMillis();
        println("success, " + (last - start));
    }

    public void writeCSV(String suffix) {
        //String annots = iEXML.iRecord[EmbryoXML.ANNOTS];
        String annots = iAnnots;
//        String useDats = ManifestX.getManifestValue("UseDats");
//        if (useDats.equals("yes")) {
//            annots += "/dats/";
//            System.out.println("using dat prefix\n");
//        }
        annots += this.iSeries + "AuxInfo" + suffix + ".csv";
        //annots += "/dats/" + iSeries + "AuxInfo" + suffix + ".csv";

        this.iMeasureCSV.setFilePath(annots);
        this.iMeasureCSV.writeCSV();
    }

    private int getTimeScale() {
	    System.out.println("***Getting time scale***");

        EmbryoFit ef = new EmbryoFit(iSeries);
        ef.setEmbryo(iEmbryo);
        ef.getSeriesData();

        Minimisation min = new Minimisation();
        MinimFunct funct = ef.new MinimFunct(ef);

        double [] sigs = {.5, -20};
        double [] start = sigs;
        double ftol = 28;

        min.addConstraint(0, -1, 0);
        min.nelderMead(funct, start, ftol);

        // get the minimum value
        double minimum = min.getMinimum();

        // get values of y and z at minimum
        double [] b = min.getParamValues();
        iTimeSlope = b[0];
        iTimeOffset = b[1];
        iTime = (int)Math.round(iTimeSlope * SULSTONCELLSTAGE192 + iTimeOffset);

        iMeasureCSV.put(MeasureCSV.defaultAtt_v1[MeasureCSV.TSLOPE_v1], fmt4(iTimeSlope));
        iMeasureCSV.put(MeasureCSV.defaultAtt_v1[MeasureCSV.TINTERCEPT_v1], fmt4(iTimeOffset));
        iMeasureCSV.put(MeasureCSV.defaultAtt_v1[MeasureCSV.TIME_v1], fmt4(iTime));
        iMeasureCSV.put(MeasureCSV.defaultAtt_v1[MeasureCSV.ZPIXRES_v1], fmt4(iEmbryo.getZPixRes()));
        iMeasureCSV.put(MeasureCSV.defaultAtt_v1[MeasureCSV.AXIS_v1], iEmbryo.getAxis());

        //println(iSeries + C  + fmt4(iTimeSlope) + CS + fmt4(iTimeOffset) +  CS + iTime + CS + min.getMinimum() + CS + min.getNiter());
        int u = ef.compareFit(b[0], b[1], iShowDetails);
        if (minimum > 60000) {
            println("getTimeScale, big minimum, " + minimum);

        }
        return u;
    }

    private int processNucz() {
        double [] cnucZ = zposTest3();
        if (cnucZ == null) {
            println("Failure in processNucz, incomplete nuclei record.");
            return 1;
        }
        int len = cnucZ.length;
        double half = cnucZ[len - 1] / 2;
        int prev = 0;
        int next = 0;
        double cnext = 0;
        double cprev = 0;
        for (int i=0; i < len; i++) {
            if (cnucZ[i] > half) {
                next = i;
                prev = i - 1;
                cnext = cnucZ[next];
                cprev = cnucZ[prev];
                break;
            }
        }

        double frac = (half - cprev) / (cnext - cprev);
        double center = prev + frac;
        iZCenter = center;
        iCenter= (int)Math.round(center);
        int low = iCenter - 5;
        int high = iCenter + 5;
        if (high > iMaxPlane || low < 0) {
            println("Failure in processNucz, incomplete nuclei record.");
            return 1;
        }
        double [] n = new double[high - low + 1];
        double [] p = new double[high - low + 1];
        for (int i=low; i <= high; i++) {
            n[i - low] = cnucZ[i];
            p[i - low] = i;
            //println("processNucz, " + i + CS + n[i - low] + CS + p[i - low]);
        }
        Regression reg = new Regression(p, n);
        reg.linear();
        reg.supressPrint();
        reg.supressYYplot();
        double [] est = reg.getBestEstimates();
        double [] pv = reg.getPvalues();

        double lowbreak = -est[0]/est[1];
        double highbreak = (2 * half - est[0])/est[1];
        iZSlope = est[1];
        //println(iSeries + CS + fmt1(center) + CS + fmt1(iZSlope) + CS + fmt1(lowbreak) + CS + fmt1(highbreak));

        iMeasureCSV.put(MeasureCSV.defaultAtt_v1[MeasureCSV.ZSLOPE_v1], fmt4(iZSlope));
        iMeasureCSV.put(MeasureCSV.defaultAtt_v1[MeasureCSV.ZCENTER_v1], fmt4(iZCenter));


        double [][] pb = new double[4][len];
        for (int i=0; i < len; i++) {
            pb[0][i] = i;
            pb[1][i] = cnucZ[i];
            pb[2][i] = i;
            if (i >= lowbreak && i <= highbreak) {
                double x = est[1] * i + est[0];
                if (x > 0) pb[3][i] = x;
            } else if (i < lowbreak){
                pb[3][i] = 0;
            } else pb[3][i] = cnucZ[len - 1];
        }

        return 0;
    }

    private int fitEllipse() {
        String szc = iMeasureCSV.get("zc");
        double dzc = Double.parseDouble(szc);
        String stime = iMeasureCSV.get("time");
        //println("fitEllipse, " + szc + CS + stime);
        int time = Integer.parseInt(stime);
        Vector nuclei_record = this.iEmbryo.getNucManager().getNucleiRecord();
        Vector nuclei = (Vector)nuclei_record.get(time - 1);
        //println("fitEllipse, " + nuclei.size());
        Vector v = new Vector();
        for (int i=0; i < nuclei.size(); i++) {
            Nucleus n = (Nucleus)nuclei.get(i);
            //println("fitEllipse, " + n.x + CS + n.y + CS + n.z);
            if (n.status < 1) continue;
            double t = n.z - dzc;
            //if (Math.abs(t) > 2) continue;


            v.add(n);
        }
        //println("fitEllipse, count=" + v.size());
        double [] xdata = new double[v.size()];
        double [] ydata = new double[v.size()];
        double [] xdataSave = new double[v.size()];
        double [] ydataSave = new double[v.size()];
        double [] sdata = new double[v.size()];
        for (int i=0; i < v.size(); i++) {
            Nucleus n = (Nucleus)v.get(i);

            xdata[i] = n.x;
            ydata[i] = n.y;
            xdataSave[i] = n.x;
            ydataSave[i] = n.y;
            sdata[i] = n.size;


        }
        double [] da = iterativeFit(xdata, ydata);
        double xmean = da[0];
        double ymean = da[1];
        double ang = da[2];

        double slope = Math.tan(Math.toRadians(ang));

        if (iShowDetails) {
            double ymax = 500;
            int lastpoint = ydataSave.length;
            double [][] dataPlot = new double[4][ydataSave.length + 1];
            for (int j=0; j < ydataSave.length; j++) {
                dataPlot[0][j] = xdataSave[j];// - xmean;
                dataPlot[2][j] = xdataSave[j];// - xmean;
                dataPlot[1][j] = ymax - ydataSave[j];// - ymean;
                dataPlot[3][j] = ymax - ydataSave[j];
                int x = 0;
                int y = 0;
                int yc = 0;
                x = (int)dataPlot[0][j];
                y = (int)dataPlot[1][j];
                yc = (int)((x - xmean) * slope + ymean);
                if (yc >= 0 && yc <= ymax) {
                    dataPlot[3][j] = ymax - yc;
                    dataPlot[2][j] = x;
                }
                else {
                    dataPlot[3][j] = ymax - ymean;
                    dataPlot[2][j] = xmean;
                }
                Nucleus n = (Nucleus)v.get(j);
                if (n.status < 1) continue;
                double t = n.z - dzc;
                boolean ignore = Math.abs(t) > 2.5;
                //if (Math.abs(t) > 2) continue;
                if (ignore) {
                    dataPlot[0][j] = xmean;
                    dataPlot[2][j] = xmean;
                    dataPlot[1][j] = ymax - ymean;
                    dataPlot[3][j] = ymax - ymean;
                } else {
                    //println("fitEllipse, " + n.identity + CS + x + CS + y + CS + yc);
                }
                //println("fitEllipse, " + j + CS + fmt0(dataPlot[2][j]) + CS + fmt0(dataPlot[3][j]));
            }

            dataPlot[0][lastpoint] = 700;
            dataPlot[2][lastpoint] = xmean;
            dataPlot[1][lastpoint] = 500;
            dataPlot[3][lastpoint] = 500 - ymean;


            PlotGraph pg = new PlotGraph(dataPlot);
            int [] sline = new int[2];
            int [] spt = new int[2];
            sline[0] = 0;
            sline[1] = 1;
            spt[0] = 1;
            spt[1] = 0;
            pg.setLine(sline);
            pg.setPoint(spt);
            //pg.setYhigh(500);
            pg.plot();
        }
        Stat stats = new Stat(sdata);
        double smean = stats.mean_as_double();
        //double ymeanMod = statyMod.mean_as_double();
        double xmin = Double.MAX_VALUE;
        double xmax = -Double.MIN_VALUE;
        double ymin = Double.MAX_VALUE;
        double ymax = -Double.MIN_VALUE;
        for (int i=0; i < xdata.length; i++) {
            xmin = Math.min(xmin, xdata[i]);
            xmax = Math.max(xmax, xdata[i]);
            ymin = Math.min(ymin, ydata[i]);
            ymax = Math.max(ymax, ydata[i]);
        }

        double xcenter = (xmax + xmin)/2 + xmean;
        double ycenter = (ymax + ymin)/2 + ymean;
        //println("fitEllipse, B, center values, " + xcenter + CS + ycenter);
        double angUse = -ang;
        double xmajor = xmax - xmin + smean;
        double xminor = ymax - ymin + smean;

        iEXCenter = (int)Math.round(xcenter);
        iEYCenter = (int)Math.round(ycenter);
        iEMajor = (int)Math.round(xmajor);
        iEMinor = (int)Math.round(xminor);
        iEAngle = angUse;

        iMeasureCSV.put(MeasureCSV.defaultAtt_v1[MeasureCSV.EXCENTER_v1], fmt4(iEXCenter));
        iMeasureCSV.put(MeasureCSV.defaultAtt_v1[MeasureCSV.EYCENTER_v1], fmt4(iEYCenter));
        iMeasureCSV.put(MeasureCSV.defaultAtt_v1[MeasureCSV.EMAJOR_v1], fmt4(iEMajor));
        iMeasureCSV.put(MeasureCSV.defaultAtt_v1[MeasureCSV.EMINOR_v1], fmt4(iEMinor));
        iMeasureCSV.put(MeasureCSV.defaultAtt_v1[MeasureCSV.EANG_v1], fmt4(iEAngle));

        return 0;
    }


	private String getSeries() {
		println("getSeries, " + this.configManager.getConfigFileName());
		File f = new File(this.configManager.getConfigFileName());
		String fname = f.getName();
		println("getSeries, " + fname);
		String series = fname.substring(0, fname.length() - 4);
		return series;
	}

	private String getAnnots() {
		File f = new File(this.configManager.getConfigFileName());
		String annots = f.getParent();
//		String useDats = ManifestX.getManifestValue("UseDats");
//		if (useDats.equals("yes")){
//			System.out.println("using dat prefix\n");
//			annots = annots.substring(0, annots.length() - 4);
//		} else {
//			annots=annots+"/";
//        }
		return annots + ImageNameLogic.getDirectoryDelimiter(this.configManager.getConfigFileName());

	}

    private int estimateMaxPlane(Vector nuclei) {
    	double zz = 0;
    	for (int i=0; i < nuclei.size(); i++) {
    		Nucleus n = (Nucleus)nuclei.get(i);
    		zz = Math.max(zz, n.z);
    	}
    	int k = (int)Math.ceil(zz);
    	//println("estimateMaxPlane, " + k);
    	return k + 1;
    }

    private double[] zposTest3() {
		Vector nucRecord = iEmbryo.getNucManager().getNucleiRecord();
		if (iTime < 0 || iTime >= nucRecord.size()) {
			println("zposTest3 failure looking for time=" + iTime);
			return null;
		}
		Vector nuclei = (Vector)nucRecord.get(iTime);
		iMaxPlane = estimateMaxPlane(nuclei);
		//Vector nucZ = new Vector();
		double [] nucZ = new double[iMaxPlane];
		double [] z = new double[iMaxPlane];
		for (int i=0; i < nuclei.size(); i++) {
			Nucleus n = (Nucleus)nuclei.get(i);
			if (n.status == Nucleus.NILLI) continue;
			int k = (int)Math.floor(n.z);
			nucZ[k]++;
		}
		int kk = 0;
		for (int i=0; i < iMaxPlane; i++) {
			z[i] = i;
			if (i > 0) nucZ[i] += nucZ[i - 1];
		}
		return nucZ;
	}



    private double[] iterativeFit(double [] x, double [] y) {
		double [] da = new double[3];
		double angRet = 0;
		double xmean = zmean(x);
		double ymean = zmean(y);

		for (int i=0; i < 10; i++) {
			Regression r = new Regression(x, y);
			if (i < 10) {
				r.linear();
				r.supressYYplot();
			}
			else r.linearPlot();
			double [] c = r.getBestEstimates();
			double ang = Math.toDegrees(Math.atan(c[1]));
			angRet += ang;
			//println("iterativeFir, ang=" + ang);
			rotate(x, y, ang);
		}
		da[0] = xmean;
		da[1] = ymean;
		da[2] = angRet;
		return da;

	}

    private double zmean(double [] x) {
		Stat statx = new Stat(x);
		double xm = statx.mean_as_double();
		for (int i=0; i < x.length; i++) x[i] -= xm;
		return xm;
	}

    private void rotate(double [] x, double [] y, double ang) {
		for (int i=0; i < x.length; i++) {
			double [] da = handleRotation(x[i], y[i], ang);
			x[i] = da[0];
			y[i] = da[1];
		}
	}

	public double [] handleRotation(double x, double y, double ang) {
		ang = Math.toRadians(ang);
		double cosang = Math.cos(ang);
		double sinang = Math.sin(ang);
		double denom = cosang * cosang + sinang * sinang;
		double xpnum = x * cosang + y * sinang;
		double xp = xpnum / denom;
		double yp = (y - xp * sinang) / cosang;
		double [] da = new double[2];
		da[0] = xp;
		da[1] = yp;
		return da;

	}

	public String toString() {
		return iMeasureCSV.toString();
	}

	public static final int
	 SULSTONCELLSTAGE192 = 205
	 ;

	private static void println(String s) {System.out.println(s);}
    private static final DecimalFormat DF4 = new DecimalFormat("####.####");
    private static String fmt4(double d) {return DF4.format(d);}
}