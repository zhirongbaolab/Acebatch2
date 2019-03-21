import java.text.DecimalFormat;

import org.rhwlab.manifest.ManifestX;
import sulston.Measure;


public class Measure1 {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// assuming one argument which is an acetree config file
		ManifestX.reportAndUpdateManifest();
		Measure.main(args);
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
