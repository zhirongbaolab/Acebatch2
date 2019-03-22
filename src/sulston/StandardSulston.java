package sulston;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URL;
import java.text.DecimalFormat;
import java.util.Hashtable;
import java.util.Vector;

import org.rhwlab.acetree.AceTree;
import org.rhwlab.dbaccess.DBAccess;
import org.rhwlab.dbaccess.EmbryoXML;
import org.rhwlab.manifest.ManifestX;
import org.rhwlab.snight.DivisionCaller.Rule;

public class StandardSulston {

	//EmbryoXML		iEXML;
	Vector			iCells;
	Hashtable		iSSHash;

	public StandardSulston() {
		String series = "20081128_sulston";
		iCells = new Vector();
		iSSHash = new Hashtable();
		//iEXML = SeriesNormalizer.getEmbryoXML(series);
		//String annots = iEXML.iRecord[EmbryoXML.ANNOTS];
		//String filename = annots + "/dats/SCD" + series + ".csv";
		//processFile(filename);
		processResourceFile();

	}

	public Hashtable getSSHash() {
		return iSSHash;
	}

	public Vector getCellsV() {
		return iCells;
	}


	void processResourceFile() {
	    URL url = this.getClass().getResource("/org/rhwlab/sulston/SCD20081128_sulston.csv");
	    //URL url = AceTree.class.getResource("/org/rhwlab/snight/NewRules.txt");
	    InputStream istream = null;
	    try {
	        istream = url.openStream();
	        BufferedReader br = new BufferedReader(new InputStreamReader(istream));
	        String s;
	        br.readLine(); //toss the header
	        while (br.ready()) {
	            s = br.readLine();
				String [] sa = s.split(C);
				int time = Integer.parseInt(sa[2]);
				addToHash(s, sa[1]);
	            //println("readNewRules, " + s);
	        }
	        br.close();
	    } catch(Exception e) {
	        e.printStackTrace();
	    }

	}


	void processFile(String fileName) {
		try {
			FileInputStream fis = new FileInputStream(fileName);
			BufferedReader br = new BufferedReader(new InputStreamReader(fis));
			br.readLine(); //toss header
			while (br.ready()) {
				String s = br.readLine();
				//println("processFile, " + s);
				String [] sa = s.split(C);
				int time = Integer.parseInt(sa[2]);
				addToHash(s, sa[1]);
			}
		} catch(IOException ioe) {
			ioe.printStackTrace();
		}

	}

	void addToHash(String s, String cellName) {
		Vector v = (Vector)iSSHash.get(cellName);
		if (v == null) {
			v = new Vector();
			iSSHash.put(cellName, v);
			iCells.add(cellName);
		}
		v.add(s);

	}

	void outputFromHash() {
		for (int i=0; i < iCells.size(); i++) {
			String cellName = (String)iCells.get(i);
			Vector v = (Vector)iSSHash.get(cellName);
			for (int j=0; j < v.size(); j++) {
				String s = (String)v.get(j);
				println(s);
			}
		}
	}



	/**
	 * @param args
	 */
	public static void main(String[] args) {
		println("StandardSulston.main, ");
		StandardSulston ss = new StandardSulston();
		ss.outputFromHash();

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
