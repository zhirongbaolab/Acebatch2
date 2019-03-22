package sulston;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.text.DecimalFormat;
import java.util.Hashtable;
import java.util.Vector;

import org.rhwlab.acetree.AceTreeNoUI;
import org.rhwlab.analysis.Embryo;
import org.rhwlab.manifest.ManifestX;
import org.rhwlab.snight.Config;
import org.rhwlab.snight.Identity3;
import org.rhwlab.snight.MeasureCSV;
import org.rhwlab.snight.NucleiMgr;

public class SeriesNormalizer {

	String			iSeries;
	String			iAceTreeConfigPath;
	//EmbryoXML		iEXML;
	int				iEditedTP;
	double			iZPixRes;
	Normalizer		iNormalizer;
	Hashtable		iSNHash;
	Vector			iCells;
	String			iAxis;
	double			iTSlope;
	double			iTimeOffset;
	int				iMaxSulstonTime;
	Vector			iLeaves;
	String			iAnnots;
	String			iHeader;
	boolean			iMissingFiles;


	public SeriesNormalizer(String series, String acetreeConfigPath, int timepoints) {
		iSeries = series;
		iAceTreeConfigPath = acetreeConfigPath;
		iNormalizer = getNormalizer(acetreeConfigPath, timepoints);
		if (iNormalizer == null) return;
		iMissingFiles = false;
		iTSlope = iNormalizer.iTimeSlope;
		iTimeOffset = iNormalizer.iTimeOffset;
		iMaxSulstonTime = (int) Math.round((iEditedTP - iTimeOffset) / iTSlope);
		iSNHash = new Hashtable();
		iCells = new Vector();

	}

	public Hashtable getSNHash() {
		return iSNHash;
	}

	public double getTSlope() {
		return iTSlope;
	}

	public Vector getLeafVector() {
		return iLeaves;
	}

	public String getAnnots() {
		return iAnnots;
	}

	public String getHeader() {
		return iHeader;
	}

	public Normalizer getNormalizer(String acetreeConfigPath, int timepoints) {
		iEditedTP = timepoints;
		File f = new File(acetreeConfigPath);
		String annots = f.getParent();
		String annots_temp=annots;
		f = new File(annots);
		annots = f.getParent();
		iAnnots = annots;

		//String mcsvFilePath = annots + "/dats/";

		String mcsvFilePath = annots;
		String useDats = ManifestX.getManifestValue("UseDats");
		if (useDats.equals("yes")) mcsvFilePath += "/dats/";
		else{ 
		    mcsvFilePath=annots_temp;
		    mcsvFilePath += "/";}
		String mcsvFile = mcsvFilePath + iSeries + "AuxInfo.csv";
		println("Normalizer, " + mcsvFile);
		MeasureCSV mcsv = new MeasureCSV(mcsvFile);
		int g = mcsv.getGoodLinesRead();
		if (g != 2) {
			String [] args = new String[1];
			args[0] = iSeries;
			Measure.main(args);
			mcsv = new MeasureCSV(mcsvFile);
			g = mcsv.getGoodLinesRead();
			if (g != 2) return null;
		}
		mcsv = checkOutMeasureCSV(mcsv, mcsvFile);
		return new Normalizer(mcsv, iZPixRes);

	}

	MeasureCSV checkOutMeasureCSV(MeasureCSV mcsv, String mcsvFile) {
		int r = mcsv.isMeasured();
		if (r == 0) {
			String zpixres = mcsv.get(MeasureCSV.att[MeasureCSV.ZPIXRES]);
			String axis = mcsv.get(MeasureCSV.att[MeasureCSV.AXIS]);
			boolean need = zpixres.length() == 0;
			need = need && (axis.length() == 0) || axis.equals("XXX");
			if (need) {
				Embryo emb = new Embryo(iSeries, iAceTreeConfigPath);
				zpixres = fmt4(emb.getZPixRes());
				mcsv.put(MeasureCSV.att[MeasureCSV.ZPIXRES], zpixres);
				mcsv.put(MeasureCSV.att[MeasureCSV.AXIS], emb.getAxis());
				mcsv.writeCSV();
			}
			iAxis = mcsv.get(MeasureCSV.att[MeasureCSV.AXIS]);
			iZPixRes = Double.parseDouble(mcsv.get(MeasureCSV.att[MeasureCSV.ZPIXRES]));
		} else {
			String [] args = new String[1];
			args[0] = iSeries;
			Measure.main(args);
			mcsv = new MeasureCSV(mcsvFile);
		}
		return mcsv;

	}

	void normalizeSeries() {
		String annots = iAnnots; //iEXML.iRecord[EmbryoXML.ANNOTS];
		//String cdcsv = annots + "/dats/CD" + iSeries + ".csv";
		String cdcsv = annots;
		String useDats = ManifestX.getManifestValue("UseDats");
		if (useDats.equals("yes")) cdcsv += "/dats/";
		else cdcsv += "/";
		cdcsv += "CD" + iSeries + ".csv";
		iLeaves = getLeaves(annots);
		iAnnots = annots;
		println("normalizeSeries, " + cdcsv);
		try {
			FileInputStream fis = new FileInputStream(cdcsv);
			BufferedReader br = new BufferedReader(new InputStreamReader(fis));
			iHeader = br.readLine(); //toss header
			iHeader += ",zraw";
			while (br.ready()) {
				String s = br.readLine();
				String [] sa = s.split(C);
				int time = Integer.parseInt(sa[2]);
				if (time > iEditedTP) {
					//println("normalizeSeries, break, " + s);
					continue;
				}
				s = normalizePoint(sa);
				//println("normalizeSeries, " + s);
				addToHash(s, sa[1]);
			}
		} catch(IOException ioe) {
			ioe.printStackTrace();
			iMissingFiles = true;
			return;
		}
	}

	Vector getLeaves(String annots) {
		Vector v = new Vector();
		//String leafFile = annots + "/dats/LCD" + iSeries + ".csv";
		String leafFile = annots;
		String useDats = ManifestX.getManifestValue("UseDats");
		if (useDats.equals("yes")) leafFile += "/dats/";
		else leafFile += "/";
		leafFile += "LCD" + iSeries + ".csv";

		try {
			FileInputStream fis = new FileInputStream(leafFile);
			BufferedReader br = new BufferedReader(new InputStreamReader(fis));
			while (br.ready()) {
				String s = br.readLine();
				v.add(s);
			}
		} catch(IOException ioe) {
			ioe.printStackTrace();
			iMissingFiles = true;
		}
		return v;
	}

	void addToHash(String s, String cellName) {
		//println("addToHash, " + cellName + CS + s);
		Vector v = (Vector)iSNHash.get(cellName);
		if (v == null) {
			v = new Vector();
			iSNHash.put(cellName, v);
			iCells.add(cellName);
		}
		v.add(s);

	}

	public static void outputFromHash(Hashtable h, Vector c) {
		for (int i=0; i < c.size(); i++) {
			String cellName = (String)c.get(i);
			Vector v = (Vector)h.get(cellName);
			if (v == null) continue;
			for (int j=0; j < v.size(); j++) {
				String s = (String)v.get(j);
				println(s);
			}
		}
	}

	String normalizePoint(String [] saa) {
		// here we add unnormalized z onto the end of the array
		//if (saa[1].equals("MSpppp")) {
		//	println("normalizePoint, " + saa[0]);
		//}
		String [] sa = new String[saa.length + 1];
		for (int i=0; i < saa.length; i++) {
			sa[i] = saa[i];
		}
		sa[sa.length - 1] = saa[8];
		double dz = Double.parseDouble(sa[8]);

		int z = (int)Math.round(dz * iZPixRes);
		int x = Integer.parseInt(sa[9]);
		int y = Integer.parseInt(sa[10]);
		CellPosition cp = new CellPosition(x,y,z);
		cp.normalize(iNormalizer);
		cp.convertOrientation(iNormalizer, iAxis, iZPixRes);
		dz = cp.iNZ / iZPixRes;
		sa[9] = String.valueOf(cp.iNX);
		sa[10] = String.valueOf(cp.iNY);
		sa[8] = fmt1(dz);
		//if (saa[1].equals("MSpppp")) {
		//	println("normalizePoint, " + saa[0] + CS + sa[8] + CS + sa[9] + CS +sa[10] + CS + sa[saa.length]);
		//}
		return stringFromArray(sa);

	}

	public static String stringFromArray(String [] sa) {
		StringBuffer sb = new StringBuffer(sa[0]);
		for (int i=1; i < sa.length; i++) {
			sb.append(C + sa[i]);
		}
		return sb.toString();
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		println("SeriesNormalizer.main, ");


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
