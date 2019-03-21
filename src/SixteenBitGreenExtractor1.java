import analysis.SixteenBitRedBkgComp7;
import org.rhwlab.manifest.ManifestX;


public class SixteenBitGreenExtractor1 {

	/**
	 * @param args
	 */
    public static void main(String[] args) {
		ManifestX.reportAndUpdateManifest();
        SixteenBitRedBkgComp7.cSignalIsGreen = true;
        SixteenBitRedBkgComp7.main(args);

    }
    private static void println(String s) {System.out.println(s);}
    private static final String CS = ", ", C = ",";
}
