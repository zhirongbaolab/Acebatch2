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

    private static final String CS = ", ", C = ",";
    private static final String TAB = "\t";
    private static final DecimalFormat DF0 = new DecimalFormat("####");
    private static final DecimalFormat DF1 = new DecimalFormat("####.#");
    private static final DecimalFormat DF4 = new DecimalFormat("####.####");
}
