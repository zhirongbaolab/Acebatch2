import executables.Extractor;
import executables.Measure1;
import sulston.Measure;

public class Acebatch2 {

    public static void main(String[] args) {
        if (args.length < 1) {
            System.out.println("Usage requires at least 1 argument: Program name (Measure, Extractor)");
            System.exit(0);
        }

        // extract the xml file path
        String programName = args[0];

        if (programName != null && programName.length() > 0) {
            if (programName.toLowerCase().equals(MEASURE_String)) {
                Measure1.runMeasure(args);
            } else if (programName.toLowerCase().equals(EXTRACTOR_String)) {
                Extractor.runExtractor(args);
            } else {
                System.out.println("Couldn't determine program to execute from argument: " + programName + "\nMust be: Measure OR Extractor");
                System.exit(0);
            }
        }
    }

    private static String MEASURE_String = "measure";
    private static String EXTRACTOR_String = "extractor";
}
