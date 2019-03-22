import analysis.StackBkgComp7;
import org.rhwlab.manifest.ManifestX;


public class SixteenBitGreenExtractor1 {

	/**
	 * @param args
	 */
    public static void main(String[] args) {
		ManifestX.reportAndUpdateManifest();
        StackBkgComp7.cSignalIsGreen = true;
        StackBkgComp7.main(args);

    }
    private static void println(String s) {System.out.println(s);}
    private static final String CS = ", ", C = ",";
}
