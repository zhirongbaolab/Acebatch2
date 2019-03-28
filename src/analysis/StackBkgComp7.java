
package analysis;

import ij.ImagePlus;
import ij.gui.OvalRoi;
import ij.io.Opener;
import ij.process.*;
import org.rhwlab.image.ParsingLogic.ImageNameLogic;
import org.rhwlab.snight.*;
import org.rhwlab.utils.EUtils;

import java.awt.*;
import java.io.File;
import java.io.FileInputStream;
import java.text.DecimalFormat;
import java.util.Collections;
import java.util.Hashtable;
import java.util.Vector;

public class StackBkgComp7 {

    //////////////////////////////////////
    /**
     * Revisions made 03/2019
     * @author bradenkatzman
     *
     *
     */
    //////////////////////////////////////

    /** vars */
    private Config configManager;
    private NucleiMgr nucManager;
    private String extractionColor;
    private int startTimePt;
    private int endTimePt;
    private double mid;
    private double large;
    private double blot;

    /**
     * Constructor called by executables.Extractor
     */
    public StackBkgComp7(Config configManager, NucleiMgr nucManager,
                         String extractionColor, int startTimePt, int endTimePt, double mid, double large, double blot) {
        // set all the parameters
        this.configManager = configManager;
        this.nucManager = nucManager;
        this.extractionColor = extractionColor;
        this.startTimePt = startTimePt;
        this.endTimePt = endTimePt;
        this.mid = mid;
        this.large = large;
        this.blot = blot;
    }

    public void run() {
        // set time var
        iTime = 1;

        // copy the nuclei record
        nuclei_record = this.nucManager.getNucleiRecord();

        // set plane vars
        iPlane = 15;
        iStartPlane = 1;
        iEndPlane = estimateHighestPlane();
        System.out.println("End plane estimated at: " + iEndPlane);

        extract();

        System.out.println("Starting save nuclei routine");
        saveNuclei();
    }

    /**
     *  Opens the provided image file and determines the number of slices in the stack
     *
     * @return
     */
    public int estimateHighestPlane() {
        int plane = 1;
        for (; plane <= this.configManager.getImageConfig().getPlaneEnd(); plane++) {
            try {
                ImagePlus ip = new Opener().openImage(this.configManager.getImageConfig().getProvidedImageFileName(), plane);
            } catch (Exception e){
                break;
            }
        }

        System.gc(); // clean up a bit

        return --plane;

//        int plane=1;
//        for (; plane < this.nucManager.getPlaneEnd(); plane++) {
//            String imageFile = iZipTifFilePath;
//            imageFile += "/" + iTifPrefixR;
//            imageFile += makeImageName(iStartTime, plane);
//            System.out.println("trying to make image "+imageFile+"\n");
//            try {
//                FileInputStream fis = new FileInputStream(imageFile);
//            } catch(Exception e) {
//                break;
//            }
//        }
//        return (--plane);
    }

    public void extract() {
        int missedFileCount = 0;
        int firstMissedFile = 0;
        println("beginning, " + this.configManager.getImageConfig().getProvidedImageFileName());
        println("parms, " + this.mid + CS + this.large + CS + this.blot + CS + this.startTimePt + CS + this.endTimePt);

        long startTime = System.currentTimeMillis();

        // iterate over all of the time points
        System.out.println("Starting extraction with: " + ImageNameLogic.appendTimeToSingle16BitTIFPrefix((this.configManager.getImageConfig().getImagePrefixes())[0], this.startTimePt));
        for (int time = this.startTimePt; time <= this.endTimePt; time++) {
            System.out.print(".");
            // do some reporting every so often
            if (time % 50 == 0) {
                System.out.println("Extracting " + extractionColor + " color data at time point: " + time);
            }

            long timeTime = System.currentTimeMillis();
            iResultsHash = new Hashtable();
            Vector cells = new Vector();

            Vector nuclei = (Vector)nuclei_record.elementAt(time - 1);
            Nucleus n = null;

            //initialize hash of results for each nucleus at this time
            for (int i=0; i < nuclei.size(); i++) {
                n = (Nucleus)nuclei.elementAt(i);
                if (n.status == Nucleus.NILLI) continue;
                iResultsHash.put(n.identity, new Result(n));
                cells.add(n.identity);
            }


            Collections.sort(cells);
            ImageStatistics bs = null;;
            String cell = null;

            // iteratve over all of the planes in the stack
            for (int p=iStartPlane; p <= iEndPlane; p++) {
                // prepare the name of the image to be processed
                String name = ImageNameLogic.appendTimeToSingle16BitTIFPrefix((this.configManager.getImageConfig().getImagePrefixes())[0], time);


                // if the stack contains multiple color channels, split and return the correct one. Otherwise, return an imageprocessor of the entire plane
                ImageProcessor ipData = getExtractionColorData(name, p);


                if (ipData == null) {
                    missedFileCount++;
                    if (firstMissedFile == 0) firstMissedFile = time;
                    continue;
                }

                int maxnorm=1;//Short.MAX_VALUE;
                ipData.setValue(0);
                // make copy for use in blotted calcs
                // and blot out all relevant nuclei
                ImageProcessor ipBlotCopy;
                ImageProcessor ipBlotTemplate;
                if (ipData instanceof ByteProcessor){
                    maxnorm=255;
                    ipBlotCopy = new ByteProcessor(ipData.getWidth(), ipData.getHeight());
                    ipBlotCopy.copyBits(ipData, 0, 0, ByteBlitter.COPY);
                    ipBlotCopy.setValue(0);
                    ipBlotTemplate = new ByteProcessor(ipData.getWidth(), ipData.getHeight());
                    ipBlotTemplate.setValue(255);
                }else
                {
                    ipBlotCopy = new ShortProcessor(ipData.getWidth(), ipData.getHeight());
                    ipBlotCopy.copyBits(ipData, 0, 0, ShortBlitter.COPY);
                    ipBlotCopy.setValue(0);
                    ipBlotTemplate = new ShortProcessor(ipData.getWidth(), ipData.getHeight());
                    ipBlotTemplate.setValue(maxnorm);
                }
                // ipblottemplate filled in with max value now  drawing over it in zeros
                //ipblotcopy is copy of image

                ipBlotTemplate.fill();
                ipBlotTemplate.setValue(0);
                for (int i=0; i < cells.size(); i++) {
                    cell = (String)cells.get(i);
                    Result result = (Result)iResultsHash.get(cell);
                    Nucleus nn = result.n;
                    double bdia = nucDiameter(nn, p, this.blot * nn.size);
                    if (bdia > 0) {
                        int r2 = (int)Math.round(bdia/2);
                        Polygon blot = EUtils.pCircle(nn.x, nn.y, r2);
                        ipBlotCopy.fillPolygon(blot);
                        ipBlotTemplate.fillPolygon(blot);
                    }
                }
                // all drawn in in black now in both over white (1) in template

                int count = 0;
                for (int i=0; i < cells.size(); i++) {
                    cell = (String)cells.get(i);
                    //if (cell.indexOf("MSppppp") != 0) continue; //##################
                    Result result = (Result)iResultsHash.get(cell);
                    n = result.n;
                    double dl = n.size * this.large;
                    double dia = nucDiameter(n, p, dl);
                    int r = (int)Math.round(dia/2);
                    if (r <= 0) continue;
                    // here if the large sphere has an intersection in this plane
                    //print("*");
                    int d = 2*r;
                    int xl = n.x - r;
                    int yl = n.y - r;
                    //OvalRoi oroi = new OvalRoi(xl, yl, d, d);
                    //ipData.setRoi(oroi);
                    ipData.setRoi(xl, yl, d, d);
                    ipBlotCopy.setRoi(xl, yl, d, d);
                    ipBlotTemplate.setRoi(xl, yl, d, d);
                    ImageProcessor ipcell = ipData.crop();//crop of original
                    ImageProcessor ipblot = ipBlotCopy.crop();//crop of blacked out original
                    ImageProcessor ipblottemplate = ipBlotTemplate.crop();//crop of nucleus blotted template
                    ImageProcessor iptemplate;
                    if (ipcell instanceof ByteProcessor){
                        iptemplate = new ByteProcessor(ipcell.getWidth(), ipcell.getHeight());
                        iptemplate.setValue(255);
                    }
                    else{iptemplate = new ShortProcessor(ipcell.getWidth(), ipcell.getHeight());
                        iptemplate.setValue(maxnorm);
                    }

                    iptemplate.fill();
                    iptemplate.setRoi(new OvalRoi(0, 0, d, d));
                    ImageProcessor mask = null;
                    try {
                        mask = iptemplate.getMask();
                        if(ipData instanceof ByteProcessor){
                            iptemplate.copyBits(mask, 0, 0, Blitter.AND);// big circle?
                            ipcell.copyBits(mask, 0, 0, Blitter.AND);//remove small circle
                            ipblot.copyBits(mask, 0, 0, Blitter.AND);
                            ipblottemplate.copyBits(mask, 0, 0, Blitter.AND);
                        }
                        else{
                            mask.multiply(1.0/255.0);//1/0 mask
                            iptemplate.copyBits(mask, 0, 0, Blitter.MULTIPLY);// big circle?
                            ipcell.copyBits(mask, 0, 0, Blitter.MULTIPLY);//remove small circle
                            ipblot.copyBits(mask, 0, 0, Blitter.MULTIPLY);
                            ipblottemplate.copyBits(mask, 0, 0, Blitter.MULTIPLY);
                        }
                    } catch(Exception e) {
                        println("test1 exception, " + e);
                        println("test1, " + iptemplate + CS + ipcell + CS + mask);
                    }

                    // iptemplate I think is now the large circle data
                    mask = null;
                    // consider the medium sphere
                    double dm = n.size * this.mid;
                    double mdia = nucDiameter(n, p, dm);
                    if (mdia > 0) {
                        // consider the nuclear sphere
                        double ndia = nucDiameter(n, p, n.size);
                        int nr = (int)Math.round(ndia/2);
                        if (nr > 0) {
                            int nd = 2*nr;
                            int nxl = n.x - nr;
                            int nyl = n.y - nr;
                            ipData.setRoi(nxl, nyl, nd, nd);
                            ImageProcessor ipnuc = ipData.crop();
                            ipnuc.setRoi(new OvalRoi(0, 0, nd, nd));
                            mask = ipnuc.getMask();

                            ImageProcessor ipbogus;
                            if(ipnuc instanceof ByteProcessor){
                                ipnuc.copyBits(mask, 0, 0, Blitter.AND);//use
                                mask = null;
                                ipbogus = new ByteProcessor(ipnuc.getWidth(), ipnuc.getHeight());
                                ipbogus.setValue(0);
                                ipbogus.fill();
                                bs = new ByteStatistics(ipbogus);
                                //println("test1, " + p + CS + bs.mean + CS + bs.area + CS + d);
                                bs = new ByteStatistics(ipnuc);
                            }else{
                                mask.multiply(1.0/255.0);//convert to 0/1 mask as
                                ipnuc.copyBits(mask, 0, 0, Blitter.MULTIPLY);//use multiply for short
                                mask = null;
                                ipbogus = new ShortProcessor(ipnuc.getWidth(), ipnuc.getHeight());
                                ipbogus.setValue(0);
                                ipbogus.fill();
                                bs = new ShortStatistics(ipbogus);
                                bs = new ShortStatistics(ipnuc);
                            }


                            //println("extract, test1, " + p + CS + fmt4(bs.mean) + CS + bs.area + CS + d);

                            //println("stats "+(bs.mean)+" "+bs.umean+" "+ipnuc.maxValue()+" "+(ipnuc.getMaxThreshold())+" "+ipnuc.getPixelValue(10,10)+"\n");
                            int[] hist=ipnuc.getHistogram();

                            result.nucPixSum += bs.mean * bs.area;
                            result.nucAreaSum += bs.area;
                        } // done with nuc

                        //                    println("extract, test2, " + result.nucPixSum + CS + result.nucAreaSum);

                        // back to annulus
                        int rmed = (int)Math.round(mdia / 2 );
                        Polygon medium = EUtils.pCircle(d/2, d/2, rmed);
                        ipcell.setValue(0);
                        ipcell.fillPolygon(medium);
                        ipblot.setValue(0);
                        ipblot.fillPolygon(medium);
                        ipblottemplate.setValue(0);
                        ipblottemplate.fillPolygon(medium);
                        iptemplate.setValue(0);
                        iptemplate.fillPolygon(medium);

                    } // end id mdia
                    ImageStatistics bs3;
                    ImageStatistics bs2;
                    ImageStatistics bs4;
                    if( ipblot instanceof ByteProcessor){
                        bs = new ByteStatistics(ipcell);
                        bs3 = new ByteStatistics(ipblot);
                        bs2 = new ByteStatistics(iptemplate);
                        bs4 = new ByteStatistics(ipblottemplate);
                    }else{
                        bs = new ShortStatistics(ipcell);
                        bs3 = new ShortStatistics(ipblot);
                        bs2 = new ShortStatistics(iptemplate);
                        bs4 = new ShortStatistics(ipblottemplate);
                    }


                    int [] ia = bs.histogram;
                    //println("extract, test, " + ia.length + CS + p + CS + bs.area);

                    double a1 = bs.mean * bs.area;

                    double a2= bs2.mean *bs2.area /maxnorm;

                    result.nucAnnulusPixSum += a1;
                    result.nucAnnulusAreaSum += a2;
                    double a3 = bs3.mean * bs3.area;
                    double a4 = bs4.mean * bs4.area / maxnorm;
                    result.nucBlottedAnnulusPixSum += a3;
                    result.nucBlottedAnnulusAreaSum += a4;

                } // endcells
            } // end planes
            for (int i=0; i < cells.size(); i++) {
                String cc = (String)cells.get(i);
                Result res = (Result)iResultsHash.get(cc);
                Nucleus nn = res.n;
                if (res.nucPixSum > 0) {
                    double expr = res.nucPixSum / res.nucAreaSum * 1000;
                    double expr2 = res.nucAnnulusPixSum / res.nucAnnulusAreaSum * 1000;
                    double expr3 = res.nucBlottedAnnulusPixSum / res.nucBlottedAnnulusAreaSum * 1000;
                    //println(""  + cc + C  + DF0.format(expr) + C + DF0.format(expr - expr2) + CS + DF0.format(expr - expr3));
                    StringBuffer sb = new StringBuffer(nn.identity);
                    sb.append(CS + nn.rwraw);
                    sb.append(CS + nn.rwcorr2);
                    sb.append(CS + nn.rwcorr3);
                    sb.append(CS + nn.rweight);
                    sb.append(CS + nn.rcount);
                    sb.append(CS + nn.rsum);

                    nn.rwraw = (int)Math.round(expr);
                    nn.rwcorr1 = 25000;
                    nn.rwcorr2 = (int)Math.round(expr2);
                    nn.rwcorr3 = (int)Math.round(expr3);

                    nn.rweight = nn.rwraw;
                    nn.rcount = (int)Math.round(res.nucAreaSum);
                    nn.rsum = (int)Math.round(res.nucPixSum);
                    StringBuffer sb2 = new StringBuffer(nn.identity);
                    sb2.append(CS + nn.rwraw);
                    sb2.append(CS + nn.rwcorr2);
                    sb2.append(CS + nn.rwcorr3);
                    sb2.append(CS + nn.rweight);
                    sb2.append(CS + nn.rcount);
                    sb2.append(CS + nn.rsum);
                    //println(sb.toString());
                    //println(sb2.toString());

                }
            }
            long endTime = System.currentTimeMillis();
            if (time % 10 == 0){
                println("extract, " + time + CS + (endTime - timeTime) + CS + (endTime - startTime));
            }
            //System.gc();

        } // end times
        if (firstMissedFile != 0) {
            println("extract, firstMissedFile, missedFileCount, " + firstMissedFile + CS + missedFileCount);
        }
    }

    protected ImageProcessor getExtractionColorData(String imageName, int plane) {
        if (imageName == null || plane < 1) {
            System.out.println("A null image name or an invalid plane was passed to getExtractionColorData(). Make sure these values are correct and rerun. Exiting...");
            System.exit(0);
        }

        FileInputStream fis;
        ImagePlus ip = null;

        try{
            //remove check
            ip = new Opener().openImage(imageName,plane);
        }
        catch (Exception e){
            System.out.println("image does not exist" + imageName);
        }


        int markerChannel = 2;
        if (extractionColor.equals(R)){
            //System.out.println("use red channel");
            markerChannel=2;}
        else if (extractionColor.equals(G)){
            //System.out.println("use green channel");
            markerChannel=1;
        } else {
            System.out.println("Extraction color not properly specified as R or G. Exiting...");
            System.exit(0);
        }


        if (ip != null) {
            // if the image is supposed to be split (i.e. it contains two color channels), split it into the correct cropped half
            if (this.configManager.getImageConfig().getSplitStack() == 1) {
                ip = splitImage(ip,markerChannel);
                return ip.getProcessor();
            } else {
                ImageProcessor iproc = ip.getProcessor();
                // check if the image is supposed to be flipped
                if (this.configManager.getImageConfig().getFlipStack() == 1) {
                    iproc.flipHorizontal();
                }
                return iproc;
            }

        } else {
            System.out.println("Image is null in getExtractionColorData. Exiting...");
            System.exit(0);
        }

        return null;
    }

    /**
     * Takes an image stack with two color channels and returns the half of it that contains the specified color
     * @param ip
     * @param markerChannel
     * @return
     */
    private ImagePlus splitImage(ImagePlus ip, int markerChannel) {
        if (ip == null) {
            System.out.println("Image passed to splitImage is null. Exiting...");
            System.exit(0);
        }

        ImageProcessor iproc = ip.getProcessor();

        // see if the image is supposed to be flipped
        if (this.configManager.getImageConfig().getFlipStack() == 1) {
            iproc.flipHorizontal();
        }

        if (markerChannel == 1) {
            // we'll assume the green side is the right side portion of the image (recall, the image may have been flipped)
            iproc.setRoi(new Rectangle(ip.getWidth()/2, 0, ip.getWidth()/2, ip.getHeight()));
        } else if (markerChannel == 2) {
            // we'll assume the red side is the left side portion of the image
            iproc.setRoi(new Rectangle(0, 0, ip.getWidth()/2, ip.getHeight()));
        } else {
            System.out.println("The marker channel supplied to splitImage() is incorrect or unsupported. Please specify R or G extraction. Exiting...");
            System.exit(0);
        }

        // process and return the cropped image
        ImageProcessor croppedImage = iproc.crop();
        ip.setProcessor(croppedImage);

        return ip;
    }






    int             iTime;
    int             iPlane;

    Vector          nuclei_record;

    public int             iStartPlane; // for testing only otherwise = 1
    public int             iEndPlane;

    Hashtable       iResultsHash;

    public void saveNuclei() {
        String fileName = configManager.getNucleiConfig().getZipFileName();
        File file = new File(fileName);

        System.out.println("saveNuclei: " + file);
        NucZipper nz = new NucZipper(file, nucManager, configManager);

    }

    protected class Result {
        Nucleus     n;
        double      nucPixSum;
        double      nucAreaSum;
        double      nucAnnulusPixSum;
        double      nucAnnulusAreaSum;
        double      nucBlottedAnnulusPixSum;
        double      nucBlottedAnnulusAreaSum;

        public Result(Nucleus nuc) {
            n = nuc;
        }
    }

    public double nucDiameter(Nucleus n, double imgPlane, double dx) {
        if (n == null) return -1; //covers some issues re currentCell and not tracking
        double r = -0.5;
        double cellPlane = n.z;
        double R = dx/2.; //pixels
        double y = (cellPlane - imgPlane)*configManager.getNucleiConfig().getZPixRes()/R;
        double r2 = 1 - y*y;
        if (r2 >= 0.) r = Math.sqrt(r2)*R;
        return 2*r;
    }


    protected static void println(String s) {System.out.println(s);}
    protected static final String CS = ", ";
    protected static final DecimalFormat DF0 = new DecimalFormat("####");
    protected static final DecimalFormat DF1 = new DecimalFormat("####.#");
    protected static final DecimalFormat DF4 = new DecimalFormat("####.####");

    private static final String R = "R";
    private static final String G = "G";
}
