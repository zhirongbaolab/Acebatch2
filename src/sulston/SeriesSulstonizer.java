package sulston;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.Hashtable;
import java.util.Vector;

import org.rhwlab.manifest.ManifestX;


public class SeriesSulstonizer {

	String				iSeries;
	String				iAceTreeConfigPath;
	int					iTimePoints;
	StandardSulston		iSS;
	SeriesNormalizer	iSN;
	Hashtable			iSSHash;
	double				iTSlope;
	Vector				iLeaves;
	int					iMeanLeafLength;
	int 				iMaxLeafLength;
	int					iMaxLeafStart;

	public SeriesSulstonizer(String [] args, StandardSulston stdSulston) {
		iSS = stdSulston;
		iAceTreeConfigPath = args[0];
		iTimePoints = Integer.parseInt(args[1]);
		iSeries = getSeries();
		iSN = new SeriesNormalizer(iSeries, iAceTreeConfigPath, iTimePoints);
		if (iSN.iNormalizer == null) return;
		iTSlope = iSN.getTSlope();
		iSN.normalizeSeries();
		iLeaves = iSN.getLeafVector();
		iSSHash = new Hashtable();
	}

	String getSeries() {
		File f = new File(iAceTreeConfigPath);
		String fname = f.getName();
		String series = fname.substring(0, fname.length() - 4);
		return series;
	}

	void sulstonize() {
		Vector cells = iSS.getCellsV();
		Hashtable snHash = iSN.getSNHash();
		Hashtable ssHash = iSS.getSSHash();
		for (int i=0; i < cells.size(); i++) {
			String cellName = (String)cells.get(i);
			Vector nv = (Vector)snHash.get(cellName);
			//if (nv == null) continue;
			Vector v = new Vector();
			iSSHash.put(cellName, v);
			Vector sv = (Vector)ssHash.get(cellName);

			// augment for zraw
			for (int j=0; j < sv.size(); j++) {
				String ss = (String)sv.get(j);
				//println("sulstonize, before, " + ss.split(C).length);
				ss += ",0";
				//println("sulstonize,  after, " + ss.split(C).length);
				sv.set(j, ss);
			}

			boolean leafCell = iLeaves.contains(cellName);
			boolean founderCell = false;
			String ts = (String)sv.get(0);
			String [] tsa = ts.split(C);
			founderCell = tsa[2].equals("1");

			if (!leafCell && !founderCell && nv != null) {
				//println("sulstonize, " + cellName);
				String t = (String)nv.get(0);
				String [] ta = t.split(C);
				founderCell = ta[2].equals("1");
			}
			if (leafCell && nv != null) processLeafCell(sv, nv, v, cellName);
			else if (founderCell) processFounderCell(sv,nv,v, cellName);
			else if (nv != null) {
				processCell(sv, nv, v);
				//println("sulstonize, " + sv.size() + C + nv.size() + C + cellName);
			}
		}

	}

	void preprocessLeafCells() {
		int sum = 0;
		int mint = 0;
		int maxt = 0;
		int maxStart = 0;
		String minCell = "";
		String maxCell = "";
		String maxStartCell = "";
		Vector v = new Vector();
		Hashtable snHash = iSN.getSNHash();
		Hashtable ssHash = iSS.getSSHash();
		for (int i=0; i < iLeaves.size(); i++) {
			String cellName = (String)iLeaves.get(i);
			Vector vs = (Vector)ssHash.get(cellName);
			Vector vn = (Vector)snHash.get(cellName);
			if (vs == null || vn == null) continue;
			//println("preprocessLeafCells, " + cellName + CS + (vs == null));
			String ss = (String)vs.get(0);
			String [] sa1 = ss.split(C);
			int ts = Integer.parseInt(sa1[2]);
			String sn = (String)vn.get(vn.size() - 1); // the last time in CD..csv
			String [] sa2 = sn.split(C);
			int tn = Integer.parseInt(sa2[2]);
			if (tn != iSN.iEditedTP) {
				continue;
			}

			int dtn = (int)Math.ceil(((double)vn.size())/iTSlope);
			int send = ts + dtn - 1; // this is the expected sulston end time for this leaf
			if (i == 0) {
				mint = send;
				maxt = send;
				minCell = cellName;
				maxCell = cellName;
				maxStart = ts;
				maxStartCell = cellName;
			} else {
				if (send < mint) {
					mint = send;
					minCell = cellName;
				} else if (send > maxt) {
					maxt = send;
					maxCell = cellName;
				}
				if (ts > maxStart) {
					maxStart = ts;
					maxStartCell = cellName;
				}
			}
			sum += send;
			//if (cellName.startsWith("C")) println("preprocessLeafCells, " + cellName + CS + ts + CS + tn);
		}
		int mean = sum / iLeaves.size();
		iMeanLeafLength = mean;
		iMaxLeafLength = maxt;
		iMaxLeafStart = maxStart;
		println("preprocessLeafCells, " + iMeanLeafLength + CS + maxStartCell + CS + iMaxLeafStart + CS + maxStartCell);
	}

	void processCell(Vector sv, Vector nv, Vector v) {
		double ratio = (double)nv.size() / sv.size();
		for (int i=0; i < sv.size(); i++) {
			String ss = (String)sv.get(i);
			int k = (int)Math.round(ratio * i);
			k = Math.min(k, nv.size() - 1);
			String s = (String)nv.get(k);
			s = adjustTime(ss, s);
			v.add(s);
			//println("processCell, " + s);
		}
	}

	void processLeafCell(Vector sv, Vector nv, Vector v, String cellName) {
		// in a leaf cell, we stretch or shrink the data
		// we have a certain number of data points based on the editing of this leaf
		// when sulstonized we will need a different number because our series
		// and sulstons seem to be developing at different rates
		// (it may also be due to the fact that we do not necessarily take
		// our data at one minute intervals)
		// at any rate, our Measure program has estimated the ratio of
		// our time to sulston time and we can use that to
		// control the stretch/shrink operation

		double ratio = iTSlope;

		String ss = (String)sv.get(0);
		String [] sa1 = ss.split(C);
		int ts = Integer.parseInt(sa1[2]);
		int dtn = (int)Math.ceil(((double)nv.size())/iTSlope);
		int send = ts + dtn - 1; // this is the expected sulston end time for this leaf

		String sn = (String)nv.get(nv.size() - 1); // the last time in CD..csv
		String [] sa2 = sn.split(C);
		int tn = Integer.parseInt(sa2[2]);
		int tlast = send;
		int sstart = ts;
		ratio = (double)nv.size() / ((double)(tlast - sstart));
		for (int i=sstart; i < tlast; i++) {
			int j = i - sstart;
			int k = (int)Math.ceil(ratio * j);
			k = Math.min(k, nv.size() - 1);
			// we will assign data from index k to sultonized index i
			String s = (String)nv.get(k); // the data string to use

			// it is a question of tinkering with the time values in the data string
			// to make it into the sulstonized string
			String sm = String.valueOf(i);
			// we split the string so we can tinker easily
			String [] nsa = s.split(C);
			nsa[2] = sm; // there is a time in position 2 of the split array
			int q = nsa[0].indexOf(":");
			String x = nsa[0].substring(0, q + 1);
			x += sm; // there is also a time following the colon in position 0
			nsa[0] = x;
			// after tinkering we build the string back from the modified split components
			s = SeriesNormalizer.stringFromArray(nsa);
			v.add(s);

		}

	}

	void processFounderCell(Vector sv, Vector nv, Vector v, String cellName) {
		//println("processFounderCell, " + cellName);
		//if (cellName.equals("EMS")) {
		//	println("processFounderCell, tracing, " + cellName);
		//}
		double ratio = iTSlope;
		if (nv == null) {
			for (int i=0; i < sv.size(); i++) {
				String ss = (String)sv.get(i);
				v.add(ss);

			}
			return;
		}
		int m = (int)Math.floor(nv.size() / ratio);
		int start = sv.size() - m;
		if (start > 1) {
			for (int i=0; i < start; i++) {
				String ss = (String)sv.get(i);
				v.add(ss);

			}
		}
		start = Math.max(1, start);
		for (int i=start; i < sv.size(); i++) {
			String ss = (String)sv.get(i);
			//println("processFounderCell, " + ss + CS + start + CS + i + CS + nv.size);
			int k = (int)Math.floor(ratio * (i - start));

			//println("processFounderCell, " + k + CS + nv.size() + CS + start + CS + i + CS + ss);

			//if (k == nv.size()) break;
			//k = Math.min(k, nv.size() - 1);
			if (k < nv.size()) {
				String s = (String)nv.get(k);
				s = adjustTime(ss, s);
				v.add(s);
				//println("processFounderCell, adding, " + i + CS + s);
			}
			//println("processFounderCell, " + s);
		}

	}

	String adjustTime(String ss, String s) {
		String [] ssa = ss.split(C);
		String [] sa = s.split(C);
		sa[2] = ssa[2];
		sa[0] = sa[1] + ":" + sa[2];
		return SeriesNormalizer.stringFromArray(sa);

	}

	private static final boolean DEBUG = false;

	void fileOutputFromHash(Hashtable h, Vector c) {
		//String outfile = iSN.getAnnots() + "/dats/SCD" + iSeries + ".csv";
		if (iSN.iMissingFiles) {
			println("SeriesSulstonizer.fileOutputFromHash, missing CD file, cannot proceed");
			return;
		}
		String outfile = iSN.getAnnots();
		String useDats = ManifestX.getManifestValue("UseDats");
		if (useDats.equals("yes")) outfile += "/dats/";
		else outfile += "/";
		outfile += "SCD" + iSeries + ".csv";

		//outfile = "tempOutfile.csv";
		PrintWriter pw = null;
		if (!DEBUG) {
		try {
			FileOutputStream fos = new FileOutputStream(outfile);
			pw = new PrintWriter(fos, true);

		} catch(IOException ioe) {
			ioe.printStackTrace();
			return;
		}
		// for testing only
		} else pw = new PrintWriter(System.out, true);



		pw.println(iSN.getHeader());
		Vector cv = iSS.iCells;
		for (int i=0; i < cv.size(); i++) {
			String cellName = (String)cv.get(i);
			Vector v = (Vector)h.get(cellName);
			Vector ssv = (Vector)iSS.iSSHash.get(cellName);
			for (int j=0; j < ssv.size(); j++) {
				String ss = (String)ssv.get(j);
				//ss += ",0"; // augment for zraw processing
				String s = null;
				if (v != null && v.size() > j) {
					s = (String)v.get(j);
					pw.println(s);
				} else {
					pw.println(ss);
				}
			}
		}
		if (!DEBUG) pw.close();
		println("fileOutputFromHash, " + outfile);
	}

	void averagedFileOutputFromHash(Hashtable h, Vector c) {
		//String outfile = iSN.getAnnots() + "/dats/SCA" + iSeries + ".csv";
		if (iSN.iMissingFiles) {
			println("SeriesSulstonizer.averagedFileOutputFromHash, missing CD file, cannot proceed");
			return;
		}
		String outfile = iSN.getAnnots();
		String useDats = ManifestX.getManifestValue("UseDats");
		if (useDats.equals("yes")) outfile += "/dats/";
		else outfile += "/";
		outfile += "SCA" + iSeries + ".csv";
		//outfile = "tempOutfile.csv";
		PrintWriter pw = null;
		if (!DEBUG) {
		try {
			FileOutputStream fos = new FileOutputStream(outfile);
			pw = new PrintWriter(fos, true);

		} catch(IOException ioe) {
			ioe.printStackTrace();
			return;
		}
		// for testing only
		} else pw = new PrintWriter(System.out, true);

	    pw.println(iSN.getHeader());
		for (int i=0; i < c.size(); i++) {
			String cellName = (String)c.get(i);
			//if (cellName.equals("ABarp")) TEST = true;
			Vector v = (Vector)h.get(cellName);
			boolean b = v == null || v.size() == 0;
			if (b) {
				v = (Vector)iSS.iSSHash.get(cellName);
			}
			String avg = averageData(v);
			pw.println(avg);
			TEST = false;
		}
		if (!DEBUG) pw.close();
		println("file written, " + outfile);

	}

	private static boolean TEST = false;

	String averageData(Vector v) {
		//TEST = true;
		String [] sa = ((String)v.get(0)).split(C);
		if (TEST) {
			println("averageData, TEST, ");
		}
		// base the averaging divisor on the number
		// of non zero x positions in the vector
		// to account for founder cells with leading zeros
		int xpos = 9;
		int divisor = 0;
		double [] dd = new double[sa.length];
		String [] saOut = new String[sa.length];
		for (int i=0; i < v.size(); i++) {
			String s = (String)v.get(i);
			if (TEST) println("averageData, " + s);
			sa = s.split(C);
			if (i == 0) {
				saOut[0] = sa[0];
				saOut[1] = sa[1];
				saOut[2] = sa[2];
			}
			for (int j=3; j < sa.length; j++) {
				double x = Double.parseDouble(sa[j]);
				dd[j] += x;
			}
			if (dd[xpos] > 0) divisor++; //looking at green weight to test for data presence
		}


		for (int i=3; i < dd.length; i++) {
			int m = 0;
			if (divisor > 0) m = (int)Math.round(dd[i]/divisor);
			saOut[i] = String.valueOf(m);
		}
		int zloc = 8;
		int zrawloc = 13;
		if (divisor > 0) {
			saOut[zloc] = fmt1(dd[zloc]/divisor);
			saOut[zrawloc] = fmt1(dd[zrawloc]/divisor);
		}

		return SeriesNormalizer.stringFromArray(saOut);
	}



	// args[0] = acetree config file path
	// args[1] = time points
	static void sulstonizeOne(String [] args) {
		StandardSulston stdSulston = new StandardSulston();
		SeriesSulstonizer ss = new SeriesSulstonizer(args, stdSulston);
		if (ss.iSN.iNormalizer == null) {
			println("\n*********sulstonizeOne, unable to process due to bad AuxInfo, " + ss.iSeries);
			return;
		}
		//ss.preprocessLeafCells();
		ss.sulstonize();
		ss.fileOutputFromHash(ss.iSSHash, ss.iSS.getCellsV());
		ss.averagedFileOutputFromHash(ss.iSSHash, ss.iSS.getCellsV());

	}


	public static void main(String[] args) {
		//driver();
		//if (1 == 1) System.exit(0);
		sulstonizeOne(args);
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
