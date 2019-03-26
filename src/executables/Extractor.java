package executables;

import analysis.SliceBkgComp7;
import analysis.StackBkgComp7;
import org.rhwlab.image.ParsingLogic.ImageNameLogic;
import org.rhwlab.snight.Config;
import org.rhwlab.snight.NucleiMgr;
import org.rhwlab.snight.XMLConfig;

import java.io.File;
import java.util.Hashtable;

import static org.rhwlab.image.management.ImageManager._16BIT_ID;
import static org.rhwlab.image.management.ImageManager._8BIT_ID;
import static org.rhwlab.image.management.ImageManager.getImageBitDepth;

/**
 * The main extractor class. Takes an .xml config file and a set of CLAs
 * and executes extraction
 *
 * Example usage: (run from within MATLAB)
 *  system(['java -cp acebatch2.jar executables.Extractor ',xmlname,' R, 400']);
 *
 * Created: 03/2019
 * Author: @bradenkatzman
 */

public class Extractor {

    /**
     * Parse the arguments and pass the information along to a loader
     *
     * @param args
     */
    public static void main(String[] args) {
        // check the args
        if (args.length < 3) {
            pln("Usage requires at least three arguments: 1) .xml config file path, 2) extraction color (R,G) 3) end time");
            System.exit(0);
        }

        // extract the xml file path
        String filepath = args[0];

        // make sure it's a legit file before continuing
        if (new File(filepath).exists()) {
            // the following are default values that may be updated by either CLAs or tags in the .xml
            int startTimePt = 1;
            double mid = 1.2;
            double large = 2.0;
            double blot = 1.2;

            // the second CLA is the color to extract
            String extractionColor = args[1];
            if (!extractionColor.equals(R) && !extractionColor.equals(G) && !extractionColor.equals(B)) {
                pln("Extraction color argument not properly specified. Must be either R, G, or B");
            }

            // the third CLA is the last time point
            int endTimePt = Integer.parseInt(args[2]);

            // check what other parameters have been explicitly provided as CLAs
            if (args.length > 3) {
                // the third CLA is the startTime
                startTimePt = Integer.parseInt(args[3]);
            }

            if (args.length > 6) {
                mid = Double.parseDouble(args[4]);
                large = Double.parseDouble(args[5]);
                blot = Double.parseDouble(args[6]);
            }

            // now that all of the args have been processed, figure out the relevant image properties and delegate to
            // the correct extractor
            loadDatasetAndDelegateExtraction(filepath, extractionColor, startTimePt, endTimePt, mid, large, blot);

        } else {
            pln("Provided .xml filepath does not exist.");
            System.exit(0);
        }
    }

    /**
     * Given passed command line arguments, build the Image and Nuclei config objects and then pass off to the parsing
     * method to determine the correct course of action
     *
     * @param configFilePath
     * @param extractionColor
     * @param startTimePt
     * @param endTimePt
     * @param mid
     * @param large
     * @param blot
     */
    public static void loadDatasetAndDelegateExtraction(String configFilePath, String extractionColor, int startTimePt, int endTimePt, double mid, double large, double blot) {
        XMLConfig xmlConfigLoader = new XMLConfig();
        Hashtable<String, String> xmlConfigData = xmlConfigLoader.loadConfigDataFromXMLFile(configFilePath);

        // make a Config object which contain an ImageConfig and NucleiConfig object
        Config configManager = new Config(configFilePath);

        // make a NucleiMgr object from the config object
        NucleiMgr nucManager = new NucleiMgr(configManager.getNucleiConfig());

        if (configManager.getImageConfig() != null && configManager.getNucleiConfig() != null) {
            determineImagePropsAndDelegateExtraction(configManager, nucManager, extractionColor, startTimePt, endTimePt, mid, large, blot);
        }
    }

    /**
     * This method is adopted from AceTree's ImageManager.bringUpImageSeries() and will determine the
     * type of image series being loaded (slice vs. stack) and its color channel properties. It will
     * make the necessary updates to the ImageConfig object which will be returned.
     *
     * More info: In AceTree, when bringing up the image series, a series of checks and updates are made
     * to the data pointed to in the .xml file. This includes determining the bit depth, setting certain
     * properties (based on conventions and scope properties) including whether to flip or split the images).
     * These same checks are all made here to make Acebatch2 compatible with all .xml files that are supported
     * by AceTree
     */
    public static void determineImagePropsAndDelegateExtraction(Config configManager, NucleiMgr nucManager,
                                                                       String extractionColor, int startTimePt, int endTimePt, double mid, double large, double blot) {
            // first thing we need to check if whether multiple image files (corresponding to different color channels) were provided in the config file
            // these two conditions are the result of the two conventions for supplying an <image> tag in the XML file. See documentation or ImageConfig.java
            if (!configManager.getImageConfig().areMultipleImageChannelsGiven()) {
                // only one file was provided --> let's see if it exists
                String imageFile = configManager.getImageConfig().getProvidedImageFileName();
                //System.out.println("Checking on file: " + imageFile);
                if(!new File(imageFile).exists()) {
                    System.out.println("The image listed in the config file does not exist on the system. Checking if it's an 8bit image that no longer exists...");

                    // it doesn't exist. It's likely an 8bit image file name that no longer exists, so let's do a check on the
                    // file type first (not completely reliable check) and if it's 8bit, we'll try and find a 16bit image. We can't
                    // use the normal getImageBitDepth method because it assumes a real file, and we know this one does not exist
                    if (ImageNameLogic.doesImageFollow8bitDeletedConvention(imageFile)) {
                        System.out.println("The image has an 8bit file naming convention -> trying to find it's 16bit corollary....");
                        String newFileNameAttempt = ImageNameLogic.reconfigureImagePathFrom8bitTo16bit(imageFile);
                        if (!newFileNameAttempt.equals(imageFile)) {
                            System.out.println("A 16bit file name was generated from the 8bit image file name in the config file. Checking if it exists...");
                            //System.out.println(newFileNameAttempt);
                            if (new File(newFileNameAttempt).exists()) {
                                System.out.println("16bit image file exists. Updating file in ImageConfig to: " + newFileNameAttempt);
                                configManager.getImageConfig().setProvidedImageFileName(newFileNameAttempt);
                                configManager.getImageConfig().setImagePrefixes();

                                // because the image series is now known to be 16bit stacks, we specify the two assumptions about them
                                // that AceTree makes (derived from the confocal microscope): flip and split the stack
                                // ** NOTE: if these assumptions are explicitly given, we don't override them
                                configManager.getImageConfig().setUseStack(1);

                                if (!configManager.getImageConfig().isSplitStackGiven()) {
                                    configManager.getImageConfig().setSplitStack(1);
                                }

                                if (!configManager.getImageConfig().isFlipStackGiven()) {
                                    configManager.getImageConfig().setFlipStack(1);
                                }

                                // set the starting time
                                configManager.getImageConfig().setStartingIndex(ImageNameLogic.extractTimeFromImageFileName(newFileNameAttempt));

                                // call the stack extractor for the color specified by the CLAs
                                StackBkgComp7 stackExtract = new StackBkgComp7(configManager, nucManager,
                                        extractionColor, startTimePt, endTimePt, mid, large, blot);
                                stackExtract.run();

                                return;
                            } else {
                                System.out.println("16bit image file name generated from 8bit image file name does not exist on the system. Can't bring up image series. Tried image name: " +
                                        newFileNameAttempt);
                                System.exit(0);
                            }
                        } else {
                            System.out.println("Attempt to generate 16bit image file name from 8bit image file name failed. Can't bring up image series");
                            System.exit(0);
                        }
                    } else {
                        System.out.println("Provided image file doesn't follow 8bit deleted convention. Can't bring up image series.");
                        System.exit(0);
                    }
                } else {
                    // if we've reached here, either the supplied file exists, or a 16bit corollary was found and we will now proceed with that
                    if (getImageBitDepth(imageFile) == _8BIT_ID) {
                        // load this image as the first in the image series
                        configManager.getImageConfig().setUseStack(0); // in case it isn't correctly set
                        configManager.getImageConfig().setFlipStack(0);
                        configManager.getImageConfig().setSplitStack(0);


                        // set the starting time
                        configManager.getImageConfig().setStartingIndex(ImageNameLogic.extractTimeFromImageFileName(configManager.getImageConfig().getProvidedImageFileName()));

                        // call the slice extractor for the color specified by the CLAs
                        SliceBkgComp7 sliceExtract = new SliceBkgComp7(configManager, nucManager,
                                extractionColor, startTimePt, endTimePt, mid, large, blot);
                        sliceExtract.run();

                        return;
                    } else if (getImageBitDepth(imageFile) == _16BIT_ID){
                        // we now want to check whether this image file follows the iSIM or diSPIM data hierarchy conventions. If so,
                        // we'll take advantage of that knowledge and look for other files in the series

                        // check if a second color channel can be found if we assume the iSIM data output hierarchy and format
                        String secondColorChannelFromiSIM = ImageNameLogic.findSecondiSIMColorChannel(imageFile);
                        if (!secondColorChannelFromiSIM.isEmpty()) {
                            //System.out.println("ImageManager found second channel stack by assuming iSIM data structure. Loading both channels...");
                            // the assumptions for the iSIM: don't flip, don't split
                            configManager.getImageConfig().setUseStack(1);
                            if (!configManager.getImageConfig().isSplitStackGiven()) {
                                configManager.getImageConfig().setSplitStack(0);
                            }

                            if (!configManager.getImageConfig().isFlipStackGiven()) {
                                configManager.getImageConfig().setFlipStack(0);
                            }

                            // we need to add this second color channel to the image config so that its prefix will be maintained
                            configManager.getImageConfig().addColorChannelImageToConfig(secondColorChannelFromiSIM);

                            // set the starting time
                            configManager.getImageConfig().setStartingIndex(ImageNameLogic.extractTimeFromImageFileName(imageFile));

                            // call the stack extractor for the color specified by the CLAs
                            StackBkgComp7 stackExtract = new StackBkgComp7(configManager, nucManager,
                                    extractionColor, startTimePt, endTimePt, mid, large, blot);
                            stackExtract.run();

                            return;
                        }

                        // check if a second color channel can be found is we assume the diSPIM data output hierarchy and format
                        String secondColorChannelFromdiSPIM = ImageNameLogic.findSecondDiSPIMColorChannel(imageFile);
                        if (!secondColorChannelFromdiSPIM.isEmpty()) {
                            //System.out.println("ImageManager found second channel stack by assuming diSPIM data structure. Loading both channels...");
                            // the assumptions for the diSIM: don't flip, don't split
                            configManager.getImageConfig().setUseStack(1);

                            if (!configManager.getImageConfig().isSplitStackGiven()) {
                                configManager.getImageConfig().setSplitStack(0);
                            }

                            if (!configManager.getImageConfig().isFlipStackGiven()) {
                                configManager.getImageConfig().setFlipStack(0);
                            }

                            // add the second color channel to the image config so that its prefix will be maintained
                            configManager.getImageConfig().addColorChannelImageToConfig(secondColorChannelFromdiSPIM);

                            // set the starting time
                            configManager.getImageConfig().setStartingIndex(ImageNameLogic.extractTimeFromImageFileName(imageFile));

                            // call the stack extractor for the color specified by the CLAs
                            StackBkgComp7 stackExtract = new StackBkgComp7(configManager, nucManager,
                                    extractionColor, startTimePt, endTimePt, mid, large, blot);
                            stackExtract.run();

                            return;
                        }

                        // check if this is a rare case of a 16bit slice that needs to be opened as if it was an 8bit image but with higher bit depth
                        if (ImageNameLogic.isSliceImage(imageFile)) {
                            // assumptions for the general 16bit slice: don't flip, don't split
                            configManager.getImageConfig().setUseStack(0);

                            // NOT SURE IF THESE ARE APPLICABLE IN THE SLICE CASE
                            if (configManager.getImageConfig().isSplitStackGiven()) {
                                configManager.getImageConfig().setSplitStack(0);
                            }

                            if (configManager.getImageConfig().isFlipStackGiven()) {
                                configManager.getImageConfig().setFlipStack(0);
                            }

                            // set the starting time
                            configManager.getImageConfig().setStartingIndex(ImageNameLogic.extractTimeFromImageFileName(imageFile));

                            // call the slice extractor for the color specified by the CLAs
                            SliceBkgComp7 sliceExtract = new SliceBkgComp7(configManager, nucManager,
                                    extractionColor, startTimePt, endTimePt, mid, large, blot);
                            sliceExtract.run();

                            return;
                        }

                        // if none of the above options produced a second image file containing the second color channel or determined that we have a 16bit slide
                        // we'll assume that the supplied image is from a confocal microscope, i.e. the assumptions are: split, flip
                        configManager.getImageConfig().setUseStack(1);

                        if (!configManager.getImageConfig().isSplitStackGiven()) {
                            configManager.getImageConfig().setSplitStack(1);
                        }

                        if (!configManager.getImageConfig().isFlipStackGiven()) {
                            configManager.getImageConfig().setFlipStack(1);
                        }


                        // set the starting time
                        configManager.getImageConfig().setStartingIndex(ImageNameLogic.extractTimeFromImageFileName(imageFile));
                        //this.currentImageTime = this.imageConfig.getStartingIndex();

                        // call the stack extractor for the color specified by the CLAs
                        StackBkgComp7 stackExtract = new StackBkgComp7(configManager, nucManager,
                                extractionColor, startTimePt, endTimePt, mid, large, blot);
                        stackExtract.run();

                        return;
                    }
                }
            } else {
                String imageFile = "";
                // assume that the stacks specified are from the confocal microscope, i.e. assume: flip, split
                configManager.getImageConfig().setUseStack(1);

                if (!configManager.getImageConfig().isSplitStackGiven()) {
                    configManager.getImageConfig().setSplitStack(1);
                }

                if (!configManager.getImageConfig().isFlipStackGiven()) {
                    configManager.getImageConfig().setFlipStack(1);
                }

                if (configManager.getImageConfig().getNumChannels() > 3) {
                    System.out.println("WARNING: More than three image channels were supplied in the .XML file. At this point," +
                            "AceTree only supports viewing 3 channels. All image file names " +
                            "will be loaded, but only the first three will be processed and displayed.");
                }
                // multiple images were provided in the config file. we need to query them slightly differently and then check if they exist
                String[] images = configManager.getImageConfig().getImageChannels();

                for (String s : images) {
                    if (!s.isEmpty()) {
                        imageFile = s;
                        break;
                    }
                }


                // set the starting time
                configManager.getImageConfig().setStartingIndex(ImageNameLogic.extractTimeFromImageFileName(imageFile));

                // call the stack extractor for the color specified by the CLAs
                StackBkgComp7 stackExtract = new StackBkgComp7(configManager, nucManager,
                        extractionColor, startTimePt, endTimePt, mid, large, blot);
                stackExtract.run();

                return;
            }

            System.out.println("executables.Extractor.determineImagePropsAndDelegateExtraction() reached code end.");
    }

    private static void pln(String s) {System.out.println(s);}
    private static final String R = "R";
    private static final String G = "G";
    private static final String B = "B";
}
