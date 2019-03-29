package analysis;

import ij.ImagePlus;
import ij.io.Opener;
import org.rhwlab.image.ParsingLogic.ImageNameLogic;
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

    /** vars */
    private Config configManager;
    private NucleiMgr nucManager;

    private AncesTree iAncesTree;

    private String iSeriesName;

    public Hashtable iCellsByName;
    private Vector iSortedCellNames;
    private int iTime;
    private int iPlane;


    private int iStartTime;
    private int iEndTime;
    private int iStartPlane; // for testing only otherwize = 1
    private int iEndPlane;
    private Vector nuclei_record;

    private Cell iRoot;

    public Embryo(Config configManager, NucleiMgr nucManager, String seriesName) {
        this.configManager = configManager;
        this.nucManager = nucManager;
        iSeriesName = seriesName;

        // set the plane vars
        iStartPlane = 1;
        iPlane = 15;
        iEndPlane = estimateHighestPlane();

        // set the time vars
        this.iStartTime = this.configManager.getImageConfig().getStartingIndex();
        this.iTime = 1;
        this.iEndTime = this.configManager.getImageConfig().getEndingIndex();

        nuclei_record = this.nucManager.getNucleiRecord();

        AncesTree ances = new AncesTree(null, this.nucManager, iStartTime, iEndTime);
        iAncesTree = ances;
        iCellsByName = iAncesTree.getCellsByName();
        makeSortedCellNames();
        Cell root = ances.getRoot();
        iRoot = root;
    }

    public double getZPixRes() {
        return configManager.getNucleiConfig().getZPixRes();
    }

    public String getAxis() {
    	String axis = this.configManager.iAxisGiven;
    	if (axis == null) axis = "";
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

    private Nucleus getNucleus(Cell c, int t) {
    	Vector v = c.getCellData();
    	CellData cd = (CellData)(v.get(t - c.getTime()));
    	return cd.iNucleus;
    }


    public int estimateHighestPlane() {
        if (this.configManager.getImageConfig().getUseStack() == 1) { // stack mode
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
        } else { // slice mode
            int plane = 1;
            for (; plane <= this.configManager.getImageConfig().getPlaneEnd(); plane++) {
                String imageName = ImageNameLogic.appendTimeAndPlaneTo8BittifPrefix(this.configManager.getImageConfig().getImagePrefixes()[0], iTime, plane);
                if (!new File(imageName).exists()) {
                    // make sure this isn't a situation where the first plane couldn't be found. Otherwise this next call will return 0
                    if (plane == 1) {
                        System.out.println("The highest plane couldn't not be determined. This file does not exist: " + imageName + ". Exiting...");
                        System.exit(0);
                    }

                    return --plane;
                }
            }

            return --plane;
        }
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

    public NucleiMgr getNucManager() {
        return this.nucManager;
    }

    private static void println(String s) {System.out.println(s);}
}
