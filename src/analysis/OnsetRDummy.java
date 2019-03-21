package analysis;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.text.DecimalFormat;
import java.util.Vector;

public class OnsetRDummy {

	String		iSeries;
	Vector		iOnsets;


	public OnsetRDummy(String series, String leavesFile) {
		iSeries = series;
		iOnsets = new Vector();
		makeDummyLeafReports(leavesFile);

	}

	void makeDummyLeafReports(String leavesFile) {
		try {
			FileInputStream fis = new FileInputStream(leavesFile);
			//FileInputStream fis = new FileInputStream("/nfs/waterston/biowolp/sulston/leaves671.txt");
			BufferedReader br = new BufferedReader(new InputStreamReader(fis));
			while (br.ready()) {
				String s = br.readLine();
				if (s.length() < 2) continue;
				if (s.startsWith("#")) continue;
				String [] sa = s.split(C);
				String leaf = sa[0];
				//println("makeDummyLeafReport, " + leaf);
				LeafReport lr = new LeafReport();
				lr.iLeaf = leaf;
				lr.iSeries = iSeries;
				lr.iOnsetTime = 0;
				iOnsets.add(lr);
			}
			br.close();
		} catch(IOException ioe) {
			ioe.printStackTrace();
		}
	}

	public Vector getLeafReportVector() {
		return iOnsets;
	}


	/**
	 * @param args
	 */
	public static void main(String[] args) {
		println("OnsetRDummy.main, ");
		//OnsetRDummy ord = new OnsetRDummy("20080501_C01B7_1_1");


	}
    private static void println(String s) {System.out.println(s);}
    private static void print(String s) {System.out.print(s);}
    private static final String CS = ", ", C = ",", Q="\"";
    private static final String TAB = "\t";
    private static final DecimalFormat DF0 = new DecimalFormat("####");
    private static final DecimalFormat DF1 = new DecimalFormat("####.#");
    private static final DecimalFormat DF2 = new DecimalFormat("####.##");
    private static final DecimalFormat DF4 = new DecimalFormat("####.####");
    private static String fmt4(double d) {return DF4.format(d);}
    private static String fmt2(double d) {return DF2.format(d);}
    private static String fmt1(double d) {return DF1.format(d);}
    private static String fmt0(double d) {return DF0.format(d);}

}
