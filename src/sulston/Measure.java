package sulston;

import analysis.Embryo;
import flanagan.analysis.Regression;
import flanagan.analysis.Stat;
import flanagan.math.Minimisation;
import flanagan.plot.PlotGraph;
import ij.ImagePlus;
import ij.process.ByteProcessor;
import org.rhwlab.manifest.ManifestX;
import org.rhwlab.snight.Config;
import org.rhwlab.snight.MeasureCSV;
import org.rhwlab.snight.Nucleus;
import org.rhwlab.utils.EUtils;
import sulston.EmbryoFit.MinimFunct;

import java.awt.*;
import java.io.File;
import java.text.DecimalFormat;
import java.util.Vector;

public class Measure {

	String				iSeries;
	String				iAceTreeConfigPath;
	String				iAnnots;
	double				iTimeSlope;
	double				iTimeOffset;
	double				iEXCenter;
	double				iEYCenter;
	double				iEMajor;
	double				iEMinor;
	double				iEAngle;
	double				iZSlope;
	double				iZCenter;

	int					iMaxPlane;

	//EmbryoXML			iEXML;
	int					iTime;
	int					iCenter;
	Embryo				iEmbryo;
	public MeasureCSV			iMeasureCSV;
	public boolean				iShowDetails;

	double				iZPixRes;

	static String		cEmbryoDBLoc = "/nfs/waterston/embryoDB";

	public Measure(String configPath) {
		iAceTreeConfigPath = configPath;
		iMeasureCSV = new MeasureCSV(iAceTreeConfigPath);
		iSeries = getSeries();
		iAnnots = getAnnots();
		iMeasureCSV.put(MeasureCSV.defaultAtt_v1[MeasureCSV.SERIES_v1], iSeries);


        Embryo embryo = new Embryo(iSeries, iAceTreeConfigPath);
        //embryo.setEditedTimePoints(190);
        iEmbryo = embryo;

	}

	String getSeries() {
		println("getSeries, " + iAceTreeConfigPath);
		File f = new File(iAceTreeConfigPath);
		String fname = f.getName();
		println("getSeries, " + fname);
		String series = fname.substring(0, fname.length() - 4);
		return series;
	}

	String getAnnots() {
		File f = new File(iAceTreeConfigPath);
		String annots = f.getParent();
		String useDats = ManifestX.getManifestValue("UseDats");
		if (useDats.equals("yes")){
			System.out.println("using dat prefix\n");
			annots = annots.substring(0, annots.length() - 4);
		}
		else{
			annots=annots+"/";
					}
		return annots;

	}

	public int getTimeScale() {
		EmbryoFit ef = new EmbryoFit(iSeries);
		ef.setEmbryo(iEmbryo);
		ef.getSeriesData();
        Minimisation min = new Minimisation();
        MinimFunct funct = ef.new MinimFunct(ef);
        double [] sigs = {.5, -20};
        double [] start = sigs;
        double ftol = 28;
        //double ftol = 1;
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
        int u = ef.compareFit(b[0], b[1], iShowDetails);;
        if (minimum > 60000) {
        	println("getTimeScale, big minimum, " + minimum);

        }
        return u;

	}

    int estimateMaxPlane(Vector nuclei) {
    	double zz = 0;
    	for (int i=0; i < nuclei.size(); i++) {
    		Nucleus n = (Nucleus)nuclei.get(i);
    		zz = Math.max(zz, n.z);
    	}
    	int k = (int)Math.ceil(zz);
    	//println("estimateMaxPlane, " + k);
    	return k + 1;
    }

	double [] zposTest3() {
		Vector nucRecord = iEmbryo.iNucleiMgr.getNucleiRecord();
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

	public int processNucz() {
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
		//for (int i=0; i < est.length; i++) {
		//	println("processNucz, " + i + CS + est[i] + CS + pv[i]);
		//}
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
		if (iShowDetails) {
			//PlotGraph pg = new PlotGraph(pb);
			//pg.setTitle(iSeries);
			//pg.setLine(3);
			//pg.plot();
		}

		return 0;
	}

	public int fitEllipse() {
		String szc = iMeasureCSV.get("zc");
		double dzc = Double.parseDouble(szc);
		String stime = iMeasureCSV.get("time");
		//println("fitEllipse, " + szc + CS + stime);
		int time = Integer.parseInt(stime);
		Vector nuclei_record = iEmbryo.iNucleiMgr.getNucleiRecord();
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
		//println("fitEllipseXXX, " + fmt0(xmean) + CS + fmt0(ymean) + CS + fmt4(ang) + CS + fmt4(slope));


		if (iShowDetails) {
			double ymax = 500;
			int lastpoint = ydataSave.length;
			double [][] dataPlot = new double[4][ydataSave.length + 1];
				for (int j=0; j < ydataSave.length; j++) {
					dataPlot[0][j] = xdataSave[j];// - xmean;
					dataPlot[2][j] = xdataSave[j];// - xmean;
					dataPlot[1][j] = ymax - ydataSave[j];// - ymean;
					//dataPlot[3][j] = ymean + (xdataSave[j] - xmean) * slope;
					//dataPlot[3][j] = (xdataSave[j] - xmean) * slope;
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
		//PlotGraph pg3 = new PlotGraph(xdata, ydata);
		//pg3.setLine(0);
		//pg3.setGraphTitle("after rotation");
		//pg3.setGraphTitle2("after rotation");
		//pg3.plot();
		//println("fitEllipse, B, " + xmin + CS + xmax + CS + (xmax - xmin));
		//println("fitEllipse, B, " + ymin + CS + ymax + CS + (ymax - ymin));
		//println("fitElliipse, B means,  " + xmean + CS + ymean);
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
		//println("fitEllipse, " + iMeasureCSV);

		// used test some image data results
		//test(time, dzc, xcenter, ycenter, xmajor, xminor, angUse);

		return 0;
	}

	void test(int time, double dzc, double xcenter, double ycenter, double xmajor, double xminor, double angUse) {
		Config config = new Config(iAceTreeConfigPath);

		String series = iSeries;
		String imgLoc = config.iTypicalImage;
		imgLoc = new File(imgLoc).getParent();
		imgLoc = new File(imgLoc).getParent();
		int timex = time;
		int zcenter = (int)Math.round(dzc);
		EllipseFitTest eftest = new EllipseFitTest(series, imgLoc, timex, zcenter);;
		ImagePlus iplus1 = eftest.getImage();
		eftest.show();
		eftest.maximize(xcenter, ycenter, xmajor, xminor, angUse);

	}

	double [] iterativeFit(double [] x, double [] y) {
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

	double zmean(double [] x) {
		Stat statx = new Stat(x);
		double xm = statx.mean_as_double();
		for (int i=0; i < x.length; i++) x[i] -= xm;
		return xm;
	}

	void rotate(double [] x, double [] y, double ang) {
		for (int i=0; i < x.length; i++) {
			double [] da = handleRotation(x[i], y[i], ang);
			x[i] = da[0];
			y[i] = da[1];
		}
	}




	int fitEllipse(double try3) {
		String szc = iMeasureCSV.get("zc");
		double dzc = Double.parseDouble(szc);
		String stime = iMeasureCSV.get("time");
		println("fitEllipse, " + szc + CS + stime);
		int time = Integer.parseInt(stime);
		Vector nuclei_record = iEmbryo.iNucleiMgr.getNucleiRecord();
		Vector nuclei = (Vector)nuclei_record.get(time - 1);
		println("fitEllipse, " + nuclei.size());
		Vector v = new Vector();
		for (int i=0; i < nuclei.size(); i++) {
			Nucleus n = (Nucleus)nuclei.get(i);
			println("fitEllipse, " + n.x + CS + n.y + CS + n.z);
			if (n.status < 1) continue;
			double t = n.z - dzc;
			if (Math.abs(t) > 2) continue;


			v.add(n);
		}
		println("fitEllipse, count=" + v.size());
		double [] xdata = new double[v.size()];
		double [] ydata = new double[v.size()];
		double [] sdata = new double[v.size()];
		for (int i=0; i < v.size(); i++) {
			Nucleus n = (Nucleus)v.get(i);

			xdata[i] = n.x;
			ydata[i] = n.y;
			sdata[i] = n.size;

		}
		double [] da = iterativeFit(xdata, ydata);
		double xmean = da[0];
		double ymean = da[1];
		double ang = da[2];

		//if (1 == 1) System.exit(0);


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
		PlotGraph pg3 = new PlotGraph(xdata, ydata);
		pg3.setLine(0);
		pg3.setGraphTitle("after rotation");
		pg3.setGraphTitle2("after rotation");
		pg3.plot();
		println("fitEllipse, B, " + xmin + CS + xmax + CS + (xmax - xmin));
		println("fitEllipse, B, " + ymin + CS + ymax + CS + (ymax - ymin));
		println("fitElliipse, B means,  " + xmean + CS + ymean);
		double xcenter = (xmax + xmin)/2 + xmean;
		double ycenter = (ymax + ymin)/2 + ymean;
		println("fitEllipse, B, center values, " + xcenter + CS + ycenter);
		double angUse = -ang;
		double xmajor = xmax - xmin + smean;
		double xminor = ymax - ymin + smean;

		//if (1 == 1) return 0;;
		Config config = new Config(iAceTreeConfigPath);

		String series = iSeries;
		String imgLoc = config.iTypicalImage;
		imgLoc = new File(imgLoc).getParent();
		imgLoc = new File(imgLoc).getParent();
		int timex = time;
		int zcenter = (int)Math.round(dzc);
		EllipseFitTest eftest = new EllipseFitTest(series, imgLoc, timex, zcenter);;
		ImagePlus iplus1 = eftest.getImage();
		eftest.show();
		eftest.maximize(xcenter, ycenter, xmajor, xminor, angUse);
		println("fitEllipse, B, initial values, " + xcenter + CS + ycenter + CS + xmajor + CS + xminor + CS + angUse);
		eftest.testRoi(xcenter, ycenter, xmajor, xminor, angUse);

		//if (1 == 1) return 0;

		// here, we build an image using the annotations to drop white filled circles at nuclei
		ByteProcessor bp = new ByteProcessor(iplus1.getWidth(), iplus1.getHeight());
		bp.setValue(255);
		for (int i=0; i < v.size(); i++) {
			Nucleus n = (Nucleus)v.get(i);
			Polygon p = EUtils.pCircle(n.x, n.y, n.size / 2);
			bp.fillPolygon(p);

		}

		ImagePlus iplus2 = new ImagePlus("", bp);;
		EllipseFitTest eftest2 = new EllipseFitTest(series, iplus2, timex, zcenter);
		eftest2.show();
		eftest2.maximize(xcenter, ycenter, xmajor, xminor, angUse);
		println("fitEllipse, B, initial values, " + xcenter + CS + ycenter + CS + xmajor + CS + xminor + CS + angUse);
		//eftest2.testRoi(xcenter, ycenter, xmajor, xminor, angUse);

		/*
		EllipseFitTest eftest2 = new EllipseFitTest(series, iplus1, timex, zcenter);
		eftest2.show();
		eftest2.maximize(xcenter, ycenter, xmajor, xminor, angUse);
		println("fitEllipse, B, initial values, " + xcenter + CS + ycenter + CS + xmajor + CS + xminor + CS + angUse);
		eftest2.testRoi(xcenter, ycenter, xmajor, xminor, angUse);
		*/



		return 0;
	}



	int fitEllipse(int other) {
		String szc = iMeasureCSV.get("zc");
		double dzc = Double.parseDouble(szc);
		String stime = iMeasureCSV.get("time");
		println("fitEllipse, " + szc + CS + stime);
		int time = Integer.parseInt(stime);
		Vector nuclei_record = iEmbryo.iNucleiMgr.getNucleiRecord();
		Vector nuclei = (Vector)nuclei_record.get(time - 1);
		println("fitEllipse, " + nuclei.size());
		Vector v = new Vector();
		for (int i=0; i < nuclei.size(); i++) {
			Nucleus n = (Nucleus)nuclei.get(i);
			println("fitEllipse, " + n.x + CS + n.y + CS + n.z);
			if (n.status < 1) continue;
			double t = n.z - dzc;
			if (Math.abs(t) > 2) continue;


			v.add(n);
		}
		println("fitEllipse, count=" + v.size());
		double [] xdata = new double[v.size()];
		double [] ydata = new double[v.size()];
		double [] sdata = new double[v.size()];
		for (int i=0; i < v.size(); i++) {
			Nucleus n = (Nucleus)v.get(i);

			xdata[i] = n.x;
			ydata[i] = n.y;
			sdata[i] = n.size;

		}
		Regression r = new Regression(xdata, ydata);
		r.supressYYplot();
		r.linearPlot();
		double [] c = r.getBestEstimates();
		double ang = Math.toDegrees(Math.atan(c[1]));
		println("fitEllipse, 1, " + c[0] + CS + c[1] + CS + ang);



		//ang = Math.toRadians(-30);
		Stat statx = new Stat(xdata);
		Stat staty = new Stat(ydata);
		Stat stats = new Stat(sdata);
		//Stat statyMod = new Stat(ydataMod);
		double xmean = statx.mean_as_double();
		double ymean = staty.mean_as_double();
		double smean = stats.mean_as_double();
		//double ymeanMod = statyMod.mean_as_double();
		println("fitEllipse, mean values, " + xmean + CS + ymean); // + CS + ymeanMod);
		double [] x2 = new double[v.size()];
		double [] y2 = new double[v.size()];
		double [] y2m = new double[v.size()];
		double [] x2r = new double[v.size()];
		double [] y2r = new double[v.size()];
		double xmin = Double.MAX_VALUE;
		double xmax = -Double.MIN_VALUE;
		double ymin = Double.MAX_VALUE;
		double ymax = -Double.MIN_VALUE;
		for (int i=0; i < x2.length; i++) {
			x2[i] = xdata[i] - xmean;
			y2[i] = ydata[i] - ymean;
			double [] da = handleRotation(x2[i], y2[i], 2*ang);
			x2r[i] = da[0];
			y2r[i] = da[1];
			//x2r[i] = x2[i];
			//y2r[i] = y2[i];

			StringBuffer sb = new StringBuffer("fitEllipse, B, " + i);
			sb.append(CS + fmt0(x2[i]));
			sb.append(CS + fmt0(y2[i]));
			sb.append(CS + fmt0(x2r[i]));
			sb.append(CS + fmt0(y2r[i]));
			xmin = Math.min(xmin, x2r[i]);
			xmax = Math.max(xmax, x2r[i]);
			ymin = Math.min(ymin, y2r[i]);
			ymax = Math.max(ymax, y2r[i]);
			//println(sb.toString());
		}
		Regression r2 = new Regression(x2, y2);
		r2.supressYYplot();
		r2.linearPlot();
		double [] c2 = r2.getBestEstimates();
		double ang2 = Math.toDegrees(Math.atan(c2[1]));
		println("fitEllipse, 2, " + c2[0] + CS + c2[1] + CS + ang2);

		Regression r3 = new Regression(x2r, y2r);
		r3.supressYYplot();
		r3.linearPlot();
		double [] c3 = r3.getBestEstimates();
		double ang3 = Math.toDegrees(Math.atan(c3[1]));
		println("fitEllipse, 3, " + c3[0] + CS + c3[1] + CS + ang3);

		//if (1 == 1) return 0;

		PlotGraph pg = new PlotGraph(x2, y2);
		pg.setLine(0);
		pg.setGraphTitle("before rotation");
		pg.setGraphTitle2("before rotation");
		pg.plot();
		PlotGraph pg3 = new PlotGraph(x2r, y2r);
		pg3.setLine(0);
		pg3.setGraphTitle("after rotation");
		pg3.setGraphTitle2("after rotation");
		pg3.plot();
		println("fitEllipse, B, " + xmin + CS + xmax + CS + (xmax - xmin));
		println("fitEllipse, B, " + ymin + CS + ymax + CS + (ymax - ymin));
		println("fitElliipse, B means,  " + xmean + CS + ymean);
		double xcenter = (xmax + xmin)/2 + xmean;
		double ycenter = (ymax + ymin)/2 + ymean;
		println("fitEllipse, B, center values, " + xcenter + CS + ycenter);

		//if (1 == 1) return 0;;
		Config config = new Config(iAceTreeConfigPath);

		String series = iSeries;
		String imgLoc = config.iTypicalImage;
		imgLoc = new File(imgLoc).getParent();
		imgLoc = new File(imgLoc).getParent();
		int timex = time;
		int zcenter = (int)Math.round(dzc);
		EllipseFitTest eftest = new EllipseFitTest(series, imgLoc, timex, zcenter);;
		ImagePlus iplus1 = eftest.getImage();
		iplus1.show();
		ByteProcessor bp = new ByteProcessor(iplus1.getWidth(), iplus1.getHeight());
		bp.setValue(255);
		for (int i=0; i < v.size(); i++) {
			Nucleus n = (Nucleus)v.get(i);
			Polygon p = EUtils.pCircle(n.x, n.y, n.size / 2);
			bp.fillPolygon(p);

		}

		double angUse = -2 * ang;
		double xmajor = xmax - xmin + smean;
		double xminor = ymax - ymin + smean;

		ImagePlus iplus2 = new ImagePlus("", bp);;
		eftest = new EllipseFitTest(series, iplus2, timex, zcenter);
		eftest.show();
		eftest.maximize(xcenter, ycenter, xmajor, xminor, angUse);
		//eftest.maximize(340, 248,700, 400, 0);
		println("fitEllipse, B, initial values, " + xcenter + CS + ycenter + CS + xmajor + CS + xminor + CS + angUse);
		eftest.testRoi(xcenter, ycenter, xmajor, xminor, angUse);

		EllipseFitTest eftest2 = new EllipseFitTest(series, iplus1, timex, zcenter);
		eftest2.show();
		eftest2.maximize(xcenter, ycenter, xmajor, xminor, angUse);
		println("fitEllipse, B, initial values, " + xcenter + CS + ycenter + CS + xmajor + CS + xminor + CS + angUse);
		eftest2.testRoi(xcenter, ycenter, xmajor, xminor, angUse);



		return 0;
	}

	public double [] handleRotation(double x, double y, double ang) {
		//double x = da[0];
		//double y = da[1];
		//double ang = iAng;
		//ang -= 1;
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


	int fitEllipse(boolean old) {
		Config config = new Config(iAceTreeConfigPath);
		String imgPath = config.iTypicalImage;
		String imgDir = imgPath.substring(0, imgPath.indexOf("/tif"));

		EllipseFit ef = new EllipseFit(iSeries, imgDir, iTime, iCenter);
		//EllipseFit ef = new EllipseFit(null, iTime, iCenter);
		//EllipseFit ef = new EllipseFit(iEXML, iTime, iCenter);
		int r = ef.getImage();

		//if (1 == 1) return 1;

		if (r != 0) return r;
		ef.maximize();
		if (iShowDetails) ef.show();
		//ef.show();
		iEXCenter = (int)Math.round(ef.iXCenter);
		iEYCenter = (int)Math.round(ef.iYCenter);
		iEMajor = (int)Math.round(ef.iMajor);
		iEMinor = (int)Math.round(ef.iMinor);
		iEAngle = ef.iAngle;

		iMeasureCSV.put(MeasureCSV.defaultAtt_v1[MeasureCSV.EXCENTER_v1], fmt4(iEXCenter));
		iMeasureCSV.put(MeasureCSV.defaultAtt_v1[MeasureCSV.EYCENTER_v1], fmt4(iEYCenter));
		iMeasureCSV.put(MeasureCSV.defaultAtt_v1[MeasureCSV.EMAJOR_v1], fmt4(iEMajor));
		iMeasureCSV.put(MeasureCSV.defaultAtt_v1[MeasureCSV.EMINOR_v1], fmt4(iEMinor));
		iMeasureCSV.put(MeasureCSV.defaultAtt_v1[MeasureCSV.EANG_v1], fmt4(iEAngle));
		//println("fitEllipse, " + iMeasureCSV);

		return 0;
	}

	public void writeCSV(String suffix) {
		//String annots = iEXML.iRecord[EmbryoXML.ANNOTS];
		String annots = iAnnots;
		String useDats = ManifestX.getManifestValue("UseDats");
		if (useDats.equals("yes")) {
			annots += "/dats/";
			System.out.println("using dat prefix\n");
		}
		annots += iSeries + "AuxInfo" + suffix + ".csv";
		//annots += "/dats/" + iSeries + "AuxInfo" + suffix + ".csv";
        iMeasureCSV.setFilePath(annots);
		iMeasureCSV.writeCSV();

	}

	public String toString() {
		return iMeasureCSV.toString();
	}

	public static final int
	 SULSTONCELLSTAGE192 = 205
	 ;


	/**
	 * @param args
	 */
	public static void main(String[] args) {
		//println("Measure.main, ");
		boolean showDetails = false;
		long start = System.currentTimeMillis();
		Measure meas = new Measure(args[0]);
		meas.writeCSV("");
		meas.iShowDetails = showDetails;
		int r = meas.getTimeScale();
        //println("getTimeScale, " + meas.iMeasureCSV);
		if (r != 0) return;
		r = meas.processNucz();
        //println("processNucz, " + meas.iMeasureCSV);
		if (r != 0) return;
		r = meas.fitEllipse();
		if (r != 0) {
			//println("Measure.main, " + "fitEllipse failure");
			return;
		}


		meas.writeCSV("");

		//meas.iMeasureCSV.checkHash();
		println(meas.iMeasureCSV.toString());
		long last = System.currentTimeMillis();
		println("success, " + (last - start));



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
