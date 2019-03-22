package analysis;

import ij.ImagePlus;
import ij.gui.OvalRoi;
import ij.process.Blitter;
import ij.process.ByteBlitter;
import ij.process.ByteProcessor;
import ij.process.ByteStatistics;
import ij.process.ShortBlitter;
import ij.process.ShortProcessor;
import ij.process.ShortStatistics;
import ij.process.ImageStatistics;
import ij.process.ImageProcessor;
import ij.io.Opener;


import java.awt.Polygon;
import java.io.File;
import java.io.FileInputStream;
import java.text.DecimalFormat;
import java.util.Collections;
import java.util.Hashtable;
import java.util.Vector;
import java.util.zip.ZipEntry;

import org.rhwlab.image.ImageWindow;
import org.rhwlab.image.management.ImageConfig;
import org.rhwlab.snight.*;
import org.rhwlab.utils.EUtils;
import org.rhwlab.image.ZipImage;


public class SliceBkgComp7 {

    String          iConfigFile;
    Config			iConfig;
    public NucleiMgr       iNucleiMgr;
    int             iTime;
    int             iPlane;
    String          iImgPath;

    public String          iZipTifFilePath;
    public String          iTifPrefixR;
    Vector          nuclei_record;
    String          iTitle;

    int             iUseZip;


    protected int		iFileNameType;

    double          iKMedium;
    double          iKLarge;
    double          iKBlot;
    int             iStartTime;
    int             iEndTime;
    public int             iStartPlane; // for testing only otherwize = 1
    public int             iEndPlane;



    Hashtable       iResultsHash;

    public static boolean  cSignalIsGreen = false;

    public void extractSphere(int t, int x, int y, double z, double d) {
        ImageWindow.cZipTifFilePath = iZipTifFilePath;
        String imageFile = iTifPrefixR;
        double meanSum = 0;
        double areaSum = 0;
        Vector histV = new Vector();
        for (int p=iStartPlane; p <= iEndPlane; p++) {
            String name2 = makeImageName(t, p);
	    String name;
	    if(iConfig.iUseStack==1){    
		String [] sa = imageFile.split("/");
		 name =iConfig.iZipTifFilePath+ "/"+imageFile + name2;}
	    else{	
		 name =iConfig.iZipTifFilePath+ "/"+imageFile + name2;}

            //println("test1, " +name2+ name+"\n");
            ImageProcessor ipData = getRedData2(name,p);
            double dia = circleDiameter(z, p, d);
            int r = (int)Math.round(dia/2);
            if (r <= 0) continue;
            int dd = 2*r;
            int xl = x - r;
            int yl = y - r;
            ipData.setRoi(xl, yl, dd, dd);
            ImageProcessor ipcell = ipData.crop();
	    ImageStatistics bs;
	    if(ipcell instanceof ByteProcessor){
		    bs = new ByteStatistics(ipcell);}
	    else{
		bs=new ShortStatistics(ipcell);}

            println("extractSphere, " + p + CS + fmt1(bs.mean) + CS + bs.area);
            meanSum += bs.mean * bs.area;
            areaSum += bs.area;
            histV.add(bs.histogram);
        }
        double exp = meanSum / areaSum;
        println("extractSphere, " + fmt4(exp));
        reportHistogram(histV);

    }

    void reportHistogram(Vector v) {
    	int [] rh = new int[256];
    	for (int i=0; i < v.size(); i++) {
    		int [] ia = (int[])v.get(i);
    		for (int j=0; j < rh.length; j++) {
    			rh[j] += ia[j];
    		}

    	}
    	for (int j=0; j < rh.length; j++) {
    		println("" + j + C + rh[j]);
    	}
    }

    public double circleDiameter(double z, double imgPlane, double dx) {
        double r = -0.5;
        double cellPlane = z;
        double R = dx/2.; //pixels
        double y = (cellPlane - imgPlane)*iNucleiMgr.getZPixRes()/R;
        double r2 = 1 - y*y;
        if (r2 >= 0.) r = Math.sqrt(r2)*R;
        return 2*r;
    }


    public void extract() {
    	int missedFileCount = 0;
    	int firstMissedFile = 0;
        println("beginning, " + iConfigFile);
        println("parms, " + iKMedium + CS + iKLarge + CS + iKBlot + CS + iStartTime + CS + iEndTime);
        long startTime = System.currentTimeMillis();
        for (int time = iStartTime; time <= iEndTime; time++) {
            long timeTime = System.currentTimeMillis();
            iResultsHash = new Hashtable();
            Vector cells = new Vector();
            ImageWindow.cZipTifFilePath = iZipTifFilePath;
            String imageFile = iTifPrefixR;
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
            for (int p=iStartPlane; p <= iEndPlane; p++) {
                String name2 = makeImageName(time, p);
		String name;
	    if(iConfig.iUseStack==1){    
		String [] sa = imageFile.split("/");
		 name =iConfig.iZipTifFilePath+ "/"+imageFile + name2;}
	    else{	
		 name =iConfig.iZipTifFilePath+ "/"+imageFile + name2;}

       


	    //  String name =iConfig.iZipTifFilePath+"/"+ imageFile + name2;
	     println("8 bit about to call getreddata2, " +name2+" "+ name+"\n");
                ImageProcessor ipData = (ImageProcessor)getRedData2(name,p);
		 

                if (ipData == null) {
                	missedFileCount++;
                    if (firstMissedFile == 0) firstMissedFile = time;
                    //println("test1, missing file, " + name + ", try to continue anyway");
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
                    double bdia = nucDiameter(nn, p, iKBlot * nn.size);
                    if (bdia > 0) {
                        int r2 = (int)Math.round(bdia/2);
                        Polygon blot = EUtils.pCircle(nn.x, nn.y, r2);
                        ipBlotCopy.fillPolygon(blot);
                        ipBlotTemplate.fillPolygon(blot);
                    }
                }
		// all drawin in in black now in both over white (1) in template

                int count = 0;
                for (int i=0; i < cells.size(); i++) {
                    cell = (String)cells.get(i);
                    //if (cell.indexOf("MSppppp") != 0) continue; //##################
                    Result result = (Result)iResultsHash.get(cell);
                    n = result.n;
                    double dl = n.size * iKLarge;
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
                    double dm = n.size * iKMedium;
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
			    //println("histogram stats "+ipnuc.getHistogramMin()+ipnuc.getHistogramMax()+" "+hist.length+"\n");
			    //int hcount=0;
			    //int sum=0;
			    //for (int k=0; k<hist.length; k++) {
			    //		hcount = hcount+hist[k];
			    //	sum=sum+k*hist[k];
			    //if(hist[k]>0) println("bin "+k+" "+hist[k]);
			    //}
			    //println("manual histogram mean "+hcount+" "+sum+"\n");
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
                    ipcell = null;
                    ipblot = null;
                    iptemplate = null;
                    ipblottemplate = null;
                    //System.gc();

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
                    //new ij.gui.ImageWindow(new ImagePlus(name2, ipblot)); //##########
                    //new ij.gui.ImageWindow(new ImagePlus(name2, ipblottemplate)); //##########

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
                    println(sb.toString());
                    println(sb2.toString());

                }
            }
            long endTime = System.currentTimeMillis();
            if (time % 10 == 0)println("extract, " + time + CS + (endTime - timeTime) + CS + (endTime - startTime));
            //System.gc();

        } // end times
        if (firstMissedFile != 0) println("extract, firstMissedFile, missedFileCount, " + firstMissedFile + CS + missedFileCount);

    }



    public void loadFromFile(String filePath) {
        File f = new File(filePath);
        String parent = f.getParent();
        NucleiMgr nucMgr = new NucleiMgr(filePath);
        Config c = nucMgr.getConfig();
        iConfig = c;
        iConfigFile = c.iConfigFileName;
        iImgPath = c.iTypicalImage;
        iNucleiMgr = nucMgr;
        iPlane = 15;
        iTime = 1;
        nuclei_record = iNucleiMgr.getNucleiRecord();
        //println("loadFromFile, " + c);
        iZipTifFilePath = c.iZipTifFilePath;
	if(c.iUseStack!=1){
        iTifPrefixR = makeTifPrefixR(c.iTifPrefix);
	}
	else{
	    iTifPrefixR=c.iTifPrefix;
	}
        iStartTime = c.iStartingIndex;

        iUseZip = c.iUseZip;
        ImageWindow.cUseZip = iUseZip;

        iEndPlane = estimateHighestPlane();

        iStartPlane = 1;

    }


    public int estimateHighestPlane() {
        int plane=1;
        for (; plane < iNucleiMgr.getPlaneEnd(); plane++) {
            String imageFile = iZipTifFilePath;
            imageFile += "/" + iTifPrefixR;
            imageFile += makeImageName(iStartTime, plane);
	    System.out.println("trying to create image "+imageFile+"\n");
            try {
                FileInputStream fis = new FileInputStream(imageFile);
            } catch(Exception e) {
                break;
            }
        }
        return (--plane);
    }

    private String makeTifPrefixR(String tifPrefix) {
        if (cSignalIsGreen) return tifPrefix;
        String [] sa = tifPrefix.split("/");
        return sa[0] + "R/" + sa[1];
    }
    
    protected ImageProcessor getRedData2(String greenName,int plane) {
	
		
           
	


	System.out.println("getRedData 8bit: " + greenName +" "+ plane);
        if (ImageWindow.cUseZip == 2) return getRedZipData(greenName);
        FileInputStream fis;
        ImagePlus ip = null;
        //String ss = "/home/biowolp/AncesTree/temp2/images/050405-t050-p15.tif";
	if (iConfig.iUseStack==1){
	    ip = new Opener().openImage(greenName,plane);
	}else{
	    ip = new Opener().openImage(greenName);
	}
	
	//try {
	//       fis = new FileInputStream(greenName);
	//       byte [] ba = ImageWindow.readByteArray(fis);
	//       ip = ImageWindow.openTiff(new ByteArrayInputStream(ba), false);
	//       fis.close();
	if(iConfig.iUseStack==1) {
	    
	    int markerChannel;
	    if (!cSignalIsGreen){
		//System.out.println("use red channel");
		markerChannel=2;}
	    else{
		//System.out.println("use green channel");
		markerChannel=1;
	    }
	    
	    ip=ImageWindow.splitImage(ip,markerChannel);
	    //}/
	    
	    //ip = readData(fis);
	    ///     } catch(IOException ioe) {
	    //   System.out.println("ImageWindow.test3 exception ");
	    //   System.out.println(ioe);
	}
	//ip.setDisplayRange(0,Short.MAX_VALUE);
	//ImagePlus.getDefault16bitRange();
	//	ImagePlus.setDefault16bitRange(16);//set it to not rescale?
	if (ip != null) return ip.getProcessor();
	else return null;
    }
    

    protected ImageProcessor getRedZipData(String redName) {
        ZipImage zipImage = new ZipImage(redName);
        int k1 = redName.lastIndexOf("/") + 1;
        String ss = redName.substring(k1);
        int k2 = ss.indexOf(".");
        ss = ss.substring(0, k2);
        //System.out.println("using: " + ss);
        ZipEntry ze = null;
        if (zipImage != null) ze = zipImage.getZipEntry(ss + ".tif");
        //System.out.println("ZipEntry: " + ze);
        //if (cZipImage == null) cZipImage = new ZipImage(cZipTifFilePath);
        //ZipEntry ze = cZipImage.getZipEntry(s);
        ImagePlus ip = null;
        ip = zipImage.readData(ze, true);
        if (ip != null) return ip.getProcessor();
        else return null;
    }

    public String makeImageName(int time, int plane) {
	/**
        // typical name: t001-p15.tif
        // to be augmented later to something like: images/050405-t001-p15.tif
        // which specifies a path and prefix for the set
        StringBuffer name = new StringBuffer("t");
        name.append(EUtils.makePaddedInt(time));
        name.append("-p");
        String p = EUtils.makePaddedInt(plane, 2);
        name.append(p);
        switch(iUseZip) {
        case 0:
        case 1:
            name.append(".tif");
            break;
        default:
            name.append(".zip");
        }
        return(name.toString());
	fgf
	**/
	return imageNameHandler(time, plane);

    }

    public boolean checkExists(File f)
    {
    	return f.exists();
    }
    //this is a copy of teds image name handler routine
    //easier than trying to make the original in AceTree static...
  public String imageNameHandler(int time, int plane)
    {
    	StringBuffer namebuf = new StringBuffer("t");
        namebuf.append(EUtils.makePaddedInt(time));
        if(iConfig.iUseStack == 0)
        {
        	namebuf.append("-p");
        	String p = EUtils.makePaddedInt(plane, 2);
        	namebuf.append(p);
        }
        String original_name = namebuf.toString();
      	StringBuffer namebuf2 = new StringBuffer("t");
        namebuf2.append(String.valueOf(time));
        String new_name = namebuf2.toString();
        
        //System.out.println("AceTree.java 1548: " + iZipTifFilePath + C.Fileseparator + iTifPrefix + original_name + ".tif");
        //System.out.println("AceTree.java 1548: " + iZipTifFilePath + C.Fileseparator + iTifPrefix + new_name + ".tif");
           
        if(iFileNameType == 0)
        {
        	switch(1)
        	{
        		case 0:
        		default:
        			if(checkExists(new File(iZipTifFilePath + "/" +iConfig.iTifPrefix + original_name + ".tif"))) { iFileNameType = 1; break; }
        			if(checkExists(new File(iZipTifFilePath + "/" +iConfig.iTifPrefix + new_name + ".TIF"))) { iFileNameType = 8; break; }
        			if(checkExists(new File(iZipTifFilePath +"/" + iConfig.iTifPrefix + original_name + ".TIF"))) { iFileNameType = 2; break; }
        			if(checkExists(new File(iZipTifFilePath + "/" + iConfig.iTifPrefix + original_name + ".tiff"))) { iFileNameType = 3; break; }
					if(checkExists(new File(iZipTifFilePath + "/" + iConfig.iTifPrefix + original_name + ".TIFF"))) { iFileNameType = 4; break; }
					if(checkExists(new File(iZipTifFilePath +"/"  + iConfig.iTifPrefix + original_name + ".zip"))) { iFileNameType = 5; break; }
					if(checkExists(new File(iZipTifFilePath +"/"+ iConfig.iTifPrefix + original_name + ".ZIP"))) { iFileNameType = 6; break; }
					if(checkExists(new File(iZipTifFilePath +"/"+ iConfig.iTifPrefix + new_name + ".tif"))) { iFileNameType = 7; break; }
					if(checkExists(new File(iZipTifFilePath +"/" + iConfig.iTifPrefix + new_name + ".tiff"))) { iFileNameType = 9; break; }
					if(checkExists(new File(iZipTifFilePath +"/"+ iConfig.iTifPrefix + new_name + ".TIFF"))) { iFileNameType = 10; break; }
					if(checkExists(new File(iZipTifFilePath +"/" + iConfig.iTifPrefix + new_name + ".zip"))) { iFileNameType = 11; break; }
					if(checkExists(new File(iZipTifFilePath +"/" + iConfig.iTifPrefix + new_name + ".ZIP"))) { iFileNameType = 12; break; }
        	}
        }
        
        
        //System.out.println("AceTree.java 1557: " + iFileNameType);
        
        switch(iFileNameType)
        {
        	case 1:
        		return(original_name + ".tif");
        	case 8:
        		return(new_name + ".TIF");
        	case 2:
        		return(original_name + ".TIF");
        	case 3:
        		return(original_name + ".tiff");
        	case 4:
        		return(original_name + ".TIFF");
        	case 5:
        		return(original_name + ".zip");
        	case 6:
        		return(original_name + ".ZIP");
        	case 7:
        		return(new_name + ".tif");
        	case 9:
        		return(new_name + ".tiff");
        	case 10:
        		return(new_name + ".TIFF");
        	case 11:
        		return(new_name + ".zip");
        	case 12:
        		return(new_name + ".ZIP");
       		default:
       			return(null);
        }
    }
 

    public void setParameters(int start, int end, double kMedium, double kLarge, double kBlot) {
        iStartTime = start;
        iEndTime = end;
        iKMedium = kMedium;
        iKLarge = kLarge;
        iKBlot  = kBlot;
    }

    public void saveNuclei() {
    	String fileName = iConfig.iZipFileName;
    	File file = new File(fileName);
        System.out.println("saveNuclei: " + file);
        NucZipper nz = new NucZipper(file, iNucleiMgr);
        //nz = null;
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
    /**

    public ImageProcessor getRedData2(String greenName,int plane) {
        //FileInputStream fis;
        ImagePlus ip = null;
	// ip = ImageWindow.makeImage2(greenName);
        ip = ImageWindow.makeImage2(greenName,plane,iConfig.iUseStack);
        if (ip != null) return ip.getProcessor();
        else return null;
    }
    **/

    public double nucDiameter(Nucleus n, double imgPlane, double dx) {
        if (n == null) return -1; //covers some issues re currentCell and not tracking
        double r = -0.5;
        double cellPlane = n.z;
        double R = dx/2.; //pixels
        double y = (cellPlane - imgPlane)*iNucleiMgr.getZPixRes()/R;
        double r2 = 1 - y*y;
        if (r2 >= 0.) r = Math.sqrt(r2)*R;
        return 2*r;
    }




	/**
	 * @param args
	 */
	public static void main(String[] args) {
		if (args.length < 2) {
			println("usage requires at least two args");
			System.exit(0);
		}
		SliceBkgComp7 rbc = new SliceBkgComp7();
		rbc.loadFromFile(args[0]);
        int end = Integer.parseInt(args[1]);
        int start = 1;
        if (args.length > 2) start = Integer.parseInt(args[2]);
        double mid = 1.2;
        double large = 2.0;
        double blot = 1.2;
        if (args.length > 3) {
            mid = Double.parseDouble(args[3]);
            large = Double.parseDouble(args[4]);
            blot = Double.parseDouble(args[5]);
        }
        rbc.setParameters(start, end, mid, large, blot);

        rbc.extract();
        rbc.saveNuclei();

	}

    public SliceBkgComp7(ImageConfig imageConfig, NucleiConfig nucleiConfig,
                         String extractionColor, int startTimePt, int endTimePt, double mid, double large, double blot) {
        loadFromFile(args[0]);
        int end = Integer.parseInt(args[1]);
        int start = 1;
        if (args.length > 2) start = Integer.parseInt(args[2]);
        double mid = 1.2;
        double large = 2.0;
        double blot = 1.2;
        if (args.length > 3) {
            mid = Double.parseDouble(args[3]);
            large = Double.parseDouble(args[4]);
            blot = Double.parseDouble(args[5]);
        }
        rbc.setParameters(start, end, mid, large, blot);

        rbc.extract();
        rbc.saveNuclei();
    }

    public void run() {

    }

    protected static void println(String s) {System.out.println(s);}
    protected static void print(String s) {System.out.print(s);}
    protected static final String CS = ", ", C = ",", SP=" ";
    protected static final DecimalFormat DF0 = new DecimalFormat("####");
    protected static final DecimalFormat DF1 = new DecimalFormat("####.#");
    protected static final DecimalFormat DF4 = new DecimalFormat("####.####");
    protected static String fmt4(double d) {return DF4.format(d);}
    protected static String fmt1(double d) {return DF1.format(d);}
    protected static String fmt0(double d) {return DF0.format(d);}

}
