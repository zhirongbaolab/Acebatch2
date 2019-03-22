import analysis.SliceBkgComp7;
import org.rhwlab.manifest.ManifestX;


public class GreenExtractor1 {

	/**
	 * @param args
	 */
    public static void main(String[] args) {
		ManifestX.reportAndUpdateManifest();
        SliceBkgComp7.cSignalIsGreen = true;
        SliceBkgComp7.main(args);

    }
    private static void println(String s) {System.out.println(s);}
    private static final String CS = ", ", C = ",";
}
