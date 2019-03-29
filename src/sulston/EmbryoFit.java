package sulston;

import analysis.Embryo;
import flanagan.math.MinimisationFunction;
import flanagan.plot.PlotGraph;
import org.rhwlab.snight.NucleiMgr;
import vembryo.Vembryo3;

import java.util.Vector;

public class EmbryoFit {

	String 			iSeries;
	NucleiMgr		iNucleiMgr;
	double []		iDt;
	double [] 		iDs;
	Embryo			iEmbryo;
	int				iEnd;

	public EmbryoFit(String series) {
		iSeries = series;
	}


	public void setEmbryo(Embryo emb) {
		iEmbryo = emb;
	}

	void getSeriesData() {
		Vembryo3 ves = new Vembryo3(null, false, true);
		iNucleiMgr = iEmbryo.getNucManager();
		Vector v = iNucleiMgr.getNucleiRecord();
		Vector v2 = ves.getNucleiRecord();

		int k = v.size();
		int end = Math.min(400, k);
		iEnd = end;

		//if (v.size() < end) return null;
		double [] t = new double[end];
		double [] dt = new double[end];
		double [] ds = new double[end];
		for (int i=0; i < end; i++) {
			Vector nuclei = (Vector)v.get(i);
			Vector snucs = (Vector)v2.get(i);
			t[i] = i;
			dt[i] = (nuclei.size());
			ds[i] = (snucs.size());
		}
		iDt = dt;
		iDs = ds;
	}

	// this is the function that evaluates the results of each iteration
	double testFit(double a, double b) {
		double e2 = 0;
		for (int i=0; i < 240; i++) {
			double td = a * i + b;
			int t = (int)Math.round(td);
			if (t < 0) continue;
			if (t >= iDt.length) break;
			//println("testFit, " + i + CS + t);
			double e = iDt[t] - iDs[i];
			e2 += e * e;
		}
		return e2;

	}


	int compareFit(double a, double b, boolean showDetails) {
		double [][] fit = new double[4][iDt.length];
		for (int i=0; i < iDt.length; i++) {
			double td = a * i + b;
			int t = (int)Math.round(td);
			if (t < 0) continue;
			if (t >= iDt.length) {
				println("compareFit, " + iSeries + " breaking out of processing at t=" + t);
				//return 1;
				break;
			}
			fit[0][i] = i;
			fit[2][i] = i;
			fit[1][i] = iDs[i];
			fit[3][i] = iDt[t];
		}
		if (showDetails) {
			PlotGraph pg = new PlotGraph(fit);
			pg.setLine(0);
			pg.plot();
		}
		return 0;
	}



    class MinimFunct implements MinimisationFunction{

    	EmbryoFit 	iEF;

    	public MinimFunct(EmbryoFit ef) {
    		iEF = ef;

    	}

		public double function(double[] param) {
			double e2 = 0;
			//println("function, a, " + iCount++ + CS + e2 + CS + param[0] + CS + param[1]);
			e2 = iEF.testFit(param[0], param[1]);
			//println("function, b, " + iCount++ + CS + e2 + CS + param[0] + CS + param[1]);
			return e2;
		}
    }


	/**
	 * @param args
	 */
	public static void main(String[] args) {
	}

	private static void println(String s) {System.out.println(s);}

}
