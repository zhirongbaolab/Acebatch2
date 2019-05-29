package executables;

import java.io.File;
import java.text.DecimalFormat;
import java.util.Hashtable;

import org.rhwlab.manifest.ManifestX;
import org.rhwlab.snight.Config;
import org.rhwlab.snight.NucleiMgr;
import org.rhwlab.snight.XMLConfig;
import sulston.Measure;


public class Measure1 {

	/**
	 * @param args
	 */
	public static void runMeasure(String[] args) {
	    System.out.println("Beginning Measure program.\n\n");
		// there should be one argument, a path to an xml file
        // check the args
        if (args.length < 2) {
            System.out.println("Usage requires at least two arguments: 0) Program Name 1) .xml config file path");
            System.exit(0);
        }

        // extract the xml file path
        String filepath = args[1];

        // if the file doesn't have .xml appended, do so now
        if (!filepath.endsWith(".xml")) {
            filepath += ".xml";
        }

        // make sure it's a legit file before continuing
        if (new File(filepath).exists()) {
            buildConfigsAndManagers(filepath);
        } else {
            System.out.println("Supplied file does not exist. Please check and rerun: " + filepath);
            System.exit(0);
        }
	}

	private static void buildConfigsAndManagers(String configFilePath) {
	    if (configFilePath == null) {
	        System.out.println("The config file given to Measure1.buildConfigsAndManagers is null. Ensure that it's valid and rerun. Exiting...");
	        System.exit(0);
        }
        XMLConfig xmlConfigLoader = new XMLConfig();
        Hashtable<String, String> xmlConfigData = xmlConfigLoader.loadConfigDataFromXMLFile(configFilePath);

        // make a Config object which contain an ImageConfig and NucleiConfig object
        Config configManager = new Config(configFilePath);

        // make a NucleiMgr object from the config object
        NucleiMgr nucManager = new NucleiMgr(configManager.getNucleiConfig());

        // make sure that we were actually able to build a NucleiMgr
        if (nucManager.iGoodNucleiMgr == false) {
            System.out.println("ERROR: Couldn't create a NucleiMgr object. Ensure that the .zip listed in the XML is correct " +
                    "and that the program output displays the correct path to the zip.");
            System.exit(0);
        }

        if (configManager.getImageConfig() != null && configManager.getNucleiConfig() != null) {
            delegateToMeasureRoutine(configManager, nucManager);
        }
    }

    private static void delegateToMeasureRoutine(Config configManager, NucleiMgr nucManager) {
        Measure measure = new Measure(configManager, nucManager);
        if (measure != null) {
            measure.run();
        } else {
            System.out.println("Couldn't build Measure object. Exiting...");
            System.exit(0);
        }
    }
}
