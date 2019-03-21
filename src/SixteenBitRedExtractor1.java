
import analysis.SixteenBitRedBkgComp7;
import org.rhwlab.manifest.ManifestX;

public class SixteenBitRedExtractor1 {

    /**
     * @param args
     * this version is to take a file with series lines
     */
    public static void main(String[] args) {
		ManifestX.reportAndUpdateManifest();
        SixteenBitRedBkgComp7.cSignalIsGreen = false;
        SixteenBitRedBkgComp7.main(args);

    }
    private static void println(String s) {System.out.println(s);}
    private static final String CS = ", ", C = ",";

}
