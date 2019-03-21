import analysis.RedBkgComp7;
import org.rhwlab.manifest.ManifestX;


public class GreenExtractor1 {

	/**
	 * @param args
	 */
    public static void main(String[] args) {
		ManifestX.reportAndUpdateManifest();
        RedBkgComp7.cSignalIsGreen = true;
        RedBkgComp7.main(args);

    }
    private static void println(String s) {System.out.println(s);}
    private static final String CS = ", ", C = ",";
}
