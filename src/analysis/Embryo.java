package analysis;

import org.rhwlab.snight.Config;
import org.rhwlab.snight.NucleiMgr;
import org.rhwlab.snight.Nucleus;
import org.rhwlab.tree.AncesTree;
import org.rhwlab.tree.Cell;
import org.rhwlab.tree.CellData;
import org.rhwlab.utils.EUtils;

import javax.swing.tree.TreeNode;
import java.io.*;
import java.text.DecimalFormat;
import java.util.Collections;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Vector;

public class Embryo {

    String          		iSeriesName;
    String          		iConfigFile;
    public NucleiMgr       	iNucleiMgr;
    double					iZPixRes;
    AncesTree       		iAncesTree;
    public Hashtable       	iCellsByName;
    Vector          		iSortedCellNames;
    int             		iTime;
    int             		iPlane;
    String          		iImgPath;

    String          		iZipTifFilePath;
    //String          iTifPrefixR;
    int             		iUseZip;
    int             		iStartTime;
    int             		iEndTime;
    int             		iStartPlane; // for testing only otherwize = 1
    int             		iEndPlane;
    Vector          		nuclei_record;

    Cell            		iRoot;
    String          		iConfigFilePath;
    Config          		iConfig;
    int             		iEditedTP;

    public String getName() {
        return iSeriesName;
    }

    public int getEditedTimePoints() {
        return iEditedTP;
    }

    public void setEditedTimePoints(int tp) {
        iEditedTP = tp;
    }

    public String getZipTifFilePath() {
        return iZipTifFilePath;
    }

    public Config getConfig() {
        return iConfig;
    }

    public int getEndPlane() {
        return iEndPlane;
    }

    public Vector getSortedCellNames() {
        return iSortedCellNames;
    }

    public double getZPixRes() {
    	return iZPixRes;
    }

    public String getAxis() {
    	String axis = iConfig.iAxisGiven;
    	if (axis.length() > 0) {
    		println("getAxis, " + iSeriesName + " axis specified as, "  + axis);
        	return axis;
    	}

    	try {
    		// look for 4 cell stage
    		Cell cABa = (Cell)iCellsByName.get("ABa");
    		Cell cABp = (Cell)iCellsByName.get("ABp");
    		Cell cEMS = (Cell)iCellsByName.get("EMS");
    		Cell cP2 = (Cell)iCellsByName.get("P2");
    		boolean good = cABa != null && cABp != null && cEMS != null && cP2 != null;
    		//println("getAxis, " + iSeriesName + CS + good);
    		int cABat1 = cABa.getTime();
    		int cABpt1 = cABp.getTime();
    		int cEMSt1 = cEMS.getTime();
    		int cP2t1 = cP2.getTime();
    		int t = Math.max(cABat1, cABpt1);
    		t = Math.max(t, cEMSt1);
    		t = Math.max(t, cP2t1);
    		Nucleus nABa = getNucleus(cABa, t);
    		Nucleus nABp = getNucleus(cABp, t);
    		Nucleus nEMS = getNucleus(cEMS, t);
    		Nucleus nP2 = getNucleus(cP2, t);

    		if (nABa.x < nP2.x) {
    			axis = "A";
    			if (nABp.y < nEMS.y) axis = axis + "DL";
    			else axis += "VR";

    		} else {
    			axis = "P";
    			if (nABp.y < nEMS.y) axis = axis + "DR";
    			else axis += "VL";
    		}
    	} catch(Exception e) {
    		e.printStackTrace();
    		axis = "XXX";
    	}
    	return axis;
    }

    Nucleus getNucleus(Cell c, int t) {
    	Vector v = (Vector)c.getCellData();
    	CellData cd = (CellData)(v.get(t - c.getTime()));
    	return cd.iNucleus;
    }


    public Embryo(String seriesName, String configPath) {
        iSeriesName = seriesName;
        iConfigFilePath = configPath;
        loadFromFile(configPath);
        AncesTree ances = new AncesTree(null, iNucleiMgr, iStartTime, iEndTime);
        iAncesTree = ances;
        iCellsByName = iAncesTree.getCellsByName();
        makeSortedCellNames();
        Cell root = ances.getRoot();
        iRoot = root;

    }

    public static Embryo normalizeLineage(Embryo emb, Embryo std, String lineage) {
        Vector notProcessed = new Vector();

        Cell m = (Cell)emb.iCellsByName.get(lineage);
        Enumeration e = m.preorderEnumeration();
        while (e.hasMoreElements()) {
            Cell c = (Cell)e.nextElement();
            String cname = c.getName();
            //if (cname.equals("Cappapaa")) {
            //    println("normalizeLineage, debug, " + cname);
            //}
            int oldStart = c.getTime();
            int newStart = ((Cell)c.getParent()).getEndTime() + 1;
            int oldEnd = c.getEndTime();
            int newEnd = oldEnd + newStart - oldStart;
            c.setTime(newStart);
            c.setEndTime(newEnd);
            //if (cname.length() > 1) break;
            //println("normalizeLineage, " + cname);
            Cell s = (Cell)std.iCellsByName.get(cname);
            if (s != null) c.setTime(s.getTime());

            // there are conditions under which we do not attempt to
            // map the data
            // we leave it unchanged in these cases
            //if (s == null) {
                //println("normalizeLineage, no std cell, " + cname);
                //continue;
            //}
            //if (c.getFateInt() != Cell.DIVIDED || s.getFateInt() != Cell.DIVIDED) continue;

            /*
            boolean canProcess = s != null;
            canProcess = canProcess && s.getFateInt() == Cell.DIVIDED;
            canProcess = canProcess && c.getFateInt() == Cell.DIVIDED;

            if (!canProcess) {
                notProcessed.add(c);
                StringBuffer sb = new StringBuffer("normalizeLineage, notProcessed");
                sb.append(CS + c.getName());
                sb.append(CS + c.getFate());
                if (s == null) sb.append(CS + "no std cell");
                else sb.append(CS + s + CS + s.getFate());
                println(sb.toString());
                continue;
            }
            */
            int endTime = c.getEndTime();
            if (s == null) {
            //    println("normalizeLineage, NOT, " + cname);
                continue;
            }
            if (s.getFateInt() != Cell.DIVIDED) {
            //    println("normalizeLineage, NOT, " + cname);
                continue;
            }

            endTime = s.getEndTime();

            //println("normalizeLineage, " + cname);

            //c.setTime(s.getTime());
            c.setEndTime(endTime);
            // now I want to change the CellData
            Vector vc = (Vector)c.getCellData();
            Vector vs = (Vector)s.getCellData();
            Vector vn = new Vector(); // holds new cell data
            vn.add(vc.get(0)); // first value from existing data
            if (vs.size() == vc.size()) continue; // nothing to do
            float r = ((float)(vc.size() - 1))/(vs.size() - 1);
            for (int i = 1; i < vs.size() - 1; i++) {
                float k = i * r;
                int k1 = (int)Math.floor(k);
                int k2 = (int)Math.ceil(k);
                float rx = k - k1;
                Nucleus n1 = ((CellData)vc.get(k1)).iNucleus;
                Nucleus n2 = ((CellData)vc.get(k2)).iNucleus;
                Nucleus n = interpolate(n1, n2, rx);
                //if (cname.equals("Cappapaa")) println("normalizeLineage, " + n);

                //println("test, " + i);
                //println("test, " + n1);
                //println("test, " + n);
                //println("test, " + n2);
                vn.add(new CellData(n));

            }
            CellData cdl = (CellData)vc.get(vc.size() - 1);
            Nucleus nl = cdl.iNucleus;
            Nucleus nlc = nl.copy();
            CellData cdc = new CellData(nlc);
            vn.add(cdc); // last value from current data
            c.setCellData(vn);
        }
        for (int i = notProcessed.size() - 1; i >= 0; i--) {
            Cell c = (Cell)notProcessed.get(i);
            println("normalizeLineage, removing, " + c.getName());
            c.removeFromParent();
        }

        return emb;
    }

/*
    public static Embryo normalizeLineage(Embryo emb, String lineage) {
        Embryo std = new Embryo("20060516_mir_57_wcherry");
        //Embryo std = new Embryo("102405_pha4red");
        Vector notProcessed = new Vector();

        Cell m = (Cell)emb.iCellsByName.get(lineage);
        Enumeration e = m.preorderEnumeration();
        while (e.hasMoreElements()) {
            Cell c = (Cell)e.nextElement();
            String cname = c.getName();
            //if (cname.equals("Cappapaa")) {
            //    println("normalizeLineage, debug, " + cname);
            //}
            int oldStart = c.getTime();
            int newStart = ((Cell)c.getParent()).getEndTime() + 1;
            int oldEnd = c.getEndTime();
            int newEnd = oldEnd + newStart - oldStart;
            c.setTime(newStart);
            c.setEndTime(newEnd);
            //if (cname.length() > 1) break;
            //println("normalizeLineage, " + cname);
            Cell s = (Cell)std.iCellsByName.get(cname);
            if (s != null) c.setTime(s.getTime());

            // there are conditions under which we do not attempt to
            // map the data
            // we leave it unchanged in these cases
            //if (s == null) {
                //println("normalizeLineage, no std cell, " + cname);
                //continue;
            //}
            //if (c.getFateInt() != Cell.DIVIDED || s.getFateInt() != Cell.DIVIDED) continue;

            int endTime = c.getEndTime();
            if (s == null) {
            //    println("normalizeLineage, NOT, " + cname);
                continue;
            }
            if (s.getFateInt() != Cell.DIVIDED) {
            //    println("normalizeLineage, NOT, " + cname);
                continue;
            }

            endTime = s.getEndTime();

            //println("normalizeLineage, " + cname);

            //c.setTime(s.getTime());
            c.setEndTime(endTime);
            // now I want to change the CellData
            Vector vc = (Vector)c.getCellData();
            Vector vs = (Vector)s.getCellData();
            Vector vn = new Vector(); // holds new cell data
            vn.add(vc.get(0)); // first value from existing data
            if (vs.size() == vc.size()) continue; // nothing to do
            float r = ((float)(vc.size() - 1))/(vs.size() - 1);
            for (int i = 1; i < vs.size() - 1; i++) {
                float k = i * r;
                int k1 = (int)Math.floor(k);
                int k2 = (int)Math.ceil(k);
                float rx = k - k1;
                Nucleus n1 = ((CellData)vc.get(k1)).iNucleus;
                Nucleus n2 = ((CellData)vc.get(k2)).iNucleus;
                Nucleus n = interpolate(n1, n2, rx);
                //if (cname.equals("Cappapaa")) println("normalizeLineage, " + n);

                //println("test, " + i);
                //println("test, " + n1);
                //println("test, " + n);
                //println("test, " + n2);
                vn.add(new CellData(n));

            }
            CellData cdl = (CellData)vc.get(vc.size() - 1);
            Nucleus nl = cdl.iNucleus;
            Nucleus nlc = nl.copy();
            CellData cdc = new CellData(nlc);
            vn.add(cdc); // last value from current data
            c.setCellData(vn);

        }
        for (int i = notProcessed.size() - 1; i >= 0; i--) {
            Cell c = (Cell)notProcessed.get(i);
            println("normalizeLineage, removing, " + c.getName());
            c.removeFromParent();
        }

        return emb;
    }
*/
    private static Nucleus interpolate(Nucleus n1, Nucleus n2, float rx) {
        Nucleus n = n1.copy();
        float b = 1 - rx;
        n.x = rounder(n1.x, n2.x, rx); //(int)Math.round(b * n1.x + rx * n2.x);
        n.y = rounder(n1.y, n2.y, rx);
        n.z = (b * n1.z + rx * n2.z);
        n.size = rounder(n1.size, n2.size, rx);
        n.weight = rounder(n1.weight, n2.weight, rx);
        n.rweight = rounder(n1.rweight, n2.rweight, rx);
        n.rsum = rounder(n1.rsum, n2.rsum, rx);
        n.rcount = rounder(n1.rcount, n2.rcount, rx);

        return n;
    }

    private static int rounder(int x1, int x2, float rx) {
        return (int)Math.round((1 - rx) * x1 + rx * x2);
    }

    public void processEmbryo() {
        //Enumeration e = root.depthFirstEnumeration();
        //while (e.hasMoreElements()) {
        //    Cell c = (Cell)e.nextElement();
        //    println("processEmbryo, " + c.getName() + CS + c.getTime() + CS + c.getEndTime() + CS + c.getFate());;
        //}
        Cell leaf = (Cell)iRoot.getFirstLeaf();
        //println("processEmbryo, " + c.getName() + CS + c.getTime() + CS + c.getEndTime() + CS + c.getFate());;
        TreeNode [] tn = leaf.getPath();
        for (int i=0; i < tn.length; i++) {
            Cell c = (Cell)tn[i];
            println("processEmbryo1, " + c.getName() + CS + c.getTime() + CS + c.getEndTime() + CS + c.getFate());;
        }
        leaf = (Cell)leaf.getNextLeaf();
        tn = leaf.getPath();
        for (int i=0; i < tn.length; i++) {
            Cell c = (Cell)tn[i];
            println("processEmbryo2, " + c.getName() + CS + c.getTime() + CS + c.getEndTime() + CS + c.getFate());;
        }
        println("processEmbryo3, ");
    }

    private static void processEmbryo2(Embryo embryo, String r) {
        Cell e = (Cell)embryo.iCellsByName.get(r);
        Enumeration df = e.breadthFirstEnumeration();
        while (df.hasMoreElements()) {
            Cell c = (Cell)df.nextElement();
            int num = c.getCellData().size();
            println("processEmbryo2, " + c.getName() + CS + c.getTime() + CS + c.getEndTime() + CS + c.getFate() + CS + num);;
        }


    }

    private void processCell(String cellName) {
        Cell c = (Cell)iCellsByName.get(cellName);
        Vector v = c.getCellData();
        for (int i=0; i < v.size(); i++) {
            CellData cd = (CellData)v.get(i);
            Nucleus n = cd.iNucleus;
            n.rweight = 25000;
        }
    }

    private void showNucleiData() {
        for (int i=50; i < 115; i++) {
            Nucleus n = (Nucleus)nuclei_record.get(i);
            println("showNucleiData, " + i + n.rweight);
        }
    }


    public void loadFromFile(String filePath) {
        File f = new File(filePath);
        //String parent = f.getParent();
        NucleiMgr nucMgr = new NucleiMgr(filePath);
        Config c = nucMgr.getConfig();
        iConfig = c;
        iImgPath = c.iTypicalImage;
        iNucleiMgr = nucMgr;
        iZPixRes = iNucleiMgr.getZPixRes();
        iPlane = 15;
        iTime = 1;
        nuclei_record = iNucleiMgr.getNucleiRecord();
        //println("loadFromFile, " + c);
        iZipTifFilePath = c.iZipTifFilePath;
        //iTifPrefixR = makeTifPrefixR(c.iTifPrefix);
        iStartTime = c.iStartingIndex;
        iEndTime = c.iEndingIndex;

        iEndPlane = estimateHighestPlane();

        iStartPlane = 1;

    }

    public int estimateHighestPlane() {
        int plane=1;
        for (; plane < 50; /*iNucleiMgr.getPlaneEnd()*/ plane++) {
            String imageFile = iZipTifFilePath;
            //imageFile += "/tif/";
            imageFile += "/";
            imageFile += iConfig.iTifPrefix;
            imageFile += makeImageName(iStartTime, plane);
            try {
                FileInputStream fis = new FileInputStream(imageFile);
            } catch(Exception e) {
                break;
            }
        }
        return (--plane);
    }

    private String makeImageName(int time, int plane) {
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
    }

    private void lookAtConfigFile() {
        try {
            FileInputStream fis = new FileInputStream(iConfigFilePath);
            BufferedReader br = new BufferedReader(new InputStreamReader(fis));
            while (br.ready()) {
                String s = br.readLine();
                println("lookAtConfig, " + s);
            }
        } catch(IOException ioe) {
            ioe.printStackTrace();
        }
        iConfig = Config.createConfigFromXMLFile(iConfigFilePath);
        iConfig.iZipFileName = "081305N.zip";
        iConfig.saveConfigXMLFile();
    }

    public AncesTree getAncesTree() {
        return iAncesTree;
    }

    private void makeSortedCellNames() {
        iSortedCellNames = new Vector();
        Enumeration e = iCellsByName.keys();
        int i = 0;
        while (e.hasMoreElements()) {
            iSortedCellNames.add ((String)e.nextElement());
        }
        Collections.sort(iSortedCellNames);
    }


    public static void main(String [] args) {

        println("Embryo.main, ");
        //Embryo embryo = new Embryo("nob1C", "/home/biowolp/0tmp/seriesalignment/20071015/20070210_nob1norm.xml");
        //embryo.processEmbryo2("Cap");
        String series = "20081014_ref-1_2_L2";
        //series = "102405_pha4red";
        //series = "20070801_hnd-1_F396";
        try {
        	FileInputStream fis = new FileInputStream("../acetreeAnalysis/AVRseries.txt");
        	BufferedReader br = new BufferedReader(new InputStreamReader(fis));
        	while(br.ready()) {
        		series = br.readLine();
        		if (series.length() < 2) continue;
//                Embryo embryo = Embryo.getEmbryoFromSeries(series);
                //embryo.getAxis();
//                println("Embryo.main, " + series + CS + embryo.getAxis());
        	}

        } catch(IOException ioe) {
        	ioe.printStackTrace();
        }
        //Embryo emb2 = Embryo.getEmbryoFromSeries("081305");
        //Embryo normEmb = Embryo.normalizeLineage(embryo, emb2, "C");
        //Embryo.processEmbryo2(embryo, "C");
        //println("\n\n");
        //Embryo.processEmbryo2(normEmb, "C");


    }


    private static void println(String s) {System.out.println(s);}
    private static void print(String s) {System.out.print(s);}
    private static final String CS = ", ";
    private static final DecimalFormat DF0 = new DecimalFormat("####");
    private static final DecimalFormat DF1 = new DecimalFormat("####.#");
    private static final DecimalFormat DF4 = new DecimalFormat("####.####");
    private static String fmt1(double d) {return DF1.format(d);}
    private static String fmt0(double d) {return DF1.format(d);}

}
