package vembryo;

import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.text.DecimalFormat;
import java.util.Collections;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.List;
import java.util.ListIterator;
import java.util.Vector;

import javax.swing.ImageIcon;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JScrollPane;
import javax.swing.tree.TreeNode;

import org.rhwlab.dbaccess.DBAccess;
import org.rhwlab.dbaccess.EmbryoXML;
import org.rhwlab.manifest.ManifestX;
import org.rhwlab.snight.Nucleus;
import org.rhwlab.tree.Cell;
import org.rhwlab.tree.CellData;
import org.rhwlab.tree.VTreeImpl;

public class Vembryo4 {

    String          iSeries;
    String          iCsvPath;
    public Cell     iRoot;
    Hashtable       iCellsByName;
    int             iEndTime;
    int				iLastTime;
    Vector			nuclei_record;
    String			iError;
    String			iPrefix;

    // this is a later constructor when I decided there
    // should be a series name member
    // this is essentially a constructor using only the series name
    // maybe it can be cleaned up later
    public Vembryo4(String series, String prefix) {
    	iPrefix = prefix;
        initNucRec();
        iCsvPath = getPathForSeries(series);
        iSeries = series;
        //println("Vembryo4, " + iCsvPath);
        buildFounders();
        try {
        	buildOutTree();
        } catch(Exception e) {
        	e.printStackTrace();
        	iError = iSeries;
        	println("Vembryo4, error, " + iError);
        }
        Vector v = new Vector();
        for (int i=0; i < iEndTime; i++) {
        	v.add(nuclei_record.get(i));
        }
        nuclei_record = v;
    }

    public void initNucRec() {
		iEndTime = 600;
		iLastTime = 0;
		nuclei_record = new Vector();
		for (int i=1; i <= iEndTime; i++) {
			nuclei_record.add(i - 1, new Vector());
		}
		iEndTime = 1;

    }

    public Cell getRoot() {
    	return iRoot;
    }

    public int getEndTime() {
    	return iEndTime;
    }

    public Hashtable getCellsByName() {
        return iCellsByName;
    }

    public Cell getCell(String cellName) {
        Cell c = (Cell)iCellsByName.get(cellName);
        return c;
    }

    public void reportLeafCells(String rootName, int colorChoice) {
    	Cell root = (Cell)iCellsByName.get(rootName);
    	Enumeration e = root.preorderEnumeration();
    	while (e.hasMoreElements()) {
    		Cell c = (Cell)e.nextElement();
    		if (c.isLeaf() && c.getTime() < iEndTime) {
    			Vector cd = c.getCellData();
    			int avgExpr = getAvgExpr(cd, colorChoice);
    			//println("reportLeafCells, " + c.getName() + CS + c.getTime() + CS + iEndTime + CS + avgExpr);
    		}
    	}
    }

    private int getAvgExpr(Vector cd, int colorChoice) {
    	double avgExpr = 0;
    	for (int i=0; i < cd.size(); i++) {
    		CellData cdd = (CellData)cd.get(i);
    		avgExpr += getExpression(cdd, colorChoice);
    	}
    	return (int)Math.round(avgExpr / cd.size());
    }

    private int getExpression(CellData cccd, int colorChoice) {
        Nucleus n = cccd.iNucleus;
        int red = n.rweight; // - n.rwcorr3; // here using blot corrected value

        switch (colorChoice) {
            case 1:
            	red -= n.rwcorr1; //global
            	break;
            case 2:
            	red -= n.rwcorr2; //local
            	break;
            case 3:
            	red -= n.rwcorr3; //blot
            	break;
            case 4:
            	red -= n.rwcorr4; //cross
            	break;
            default:
            	break;
        }
        return red;

    }



    public Vector getNucleiRecord() {
    	return nuclei_record;
    }

    public String [] getPathNames(String leafName, int numCells) {
        //println("getPathNames, " + leafName);
        int nUse = numCells;
        int nFilled = 0;
        String [] sa = new String[numCells];
        Cell leaf = (Cell)iCellsByName.get(leafName);
        while (leaf == null) {
            leafName = leafName.substring(0, leafName.length() - 1);
            leaf = (Cell)iCellsByName.get(leafName);
            nUse--;
        }
        for (int i=nUse; i <= numCells; i++) {
            sa[i - 1] = leafName;
            nFilled++;
        }
        // we have now padded the leaf end of sa with a good leafName

        // now look to the root end
        TreeNode [] tna = leaf.getPath();
        for (int i = 0; i < nUse; i++) {
            int k = tna.length - i - 1;
            k = Math.max(k, 0);
            Cell cc = (Cell)tna[k];
            //println("getPathName, " + i + CS + (tna.length - i - 1) + CS + (sa.length - i - 1) + CS + cc.getName());
            sa[nUse - i - 1] = cc.getName();
        }
        return sa;
    }

    // 20080729 -- hack introduced for testing servlet access
    private String getPathForSeries(String series) {
        String path = null;
        //DBAccess.cDBLocation = "/nfs/waterston/embryoDB/";
        if (DBAccess.cDBLocation == null) {
            ManifestX.reportAndUpdateManifest();
            DBAccess.cDBLocation = ManifestX.getManifestValue("DBLocation");
        }
        try {
            EmbryoXML exml = new EmbryoXML(series);
            String annots = exml.iRecord[EmbryoXML.ANNOTS];
            path = annots + "/dats/" + iPrefix + series + ".csv";
            String endTime = exml.iRecord[EmbryoXML.EDITEDTP];
            iEndTime = Integer.parseInt(endTime);
        } catch(Exception fnfe) {
            fnfe.printStackTrace();
        }
        return path;
    }

    // 20080729 -- hack introduced for testing servlet access
    public static boolean checkPathForSeries(String series) {
    	if (series == null || series.length() == 0) return true;
        String path = null;
        DBAccess.cDBLocation = "/nfs/waterston/embryoDBnew/";
        if (DBAccess.cDBLocation == null) {
            ManifestX.reportAndUpdateManifest();
            DBAccess.cDBLocation = ManifestX.getManifestValue("DBLocation");
        }
        try {
            EmbryoXML exml = new EmbryoXML(series);
            String annots = exml.iRecord[EmbryoXML.ANNOTS];
            path = annots + "/dats/CD" + series + ".csv";
            //String endTime = exml.iRecord[EmbryoXML.EDITEDTP];
            //iEndTime = Integer.parseInt(endTime);
        } catch(Exception fnfe) {
            fnfe.printStackTrace();
        }
        int locCD = path.indexOf("CD");
        String s1 = path.substring(0, locCD);
        String s2 = path.substring(locCD);
        path = s1 + "S" + s2;
        File f = new File(path);
        return f.exists();
    }



    public void buildOutTree() throws Exception {
        //try {
            FileInputStream fis = new FileInputStream(iCsvPath);
            BufferedReader br = new BufferedReader(new InputStreamReader(fis));
            br.readLine(); //toss the header line
            int count = 0;
            while (br.ready()) {
                String s = br.readLine();
                if (s.length() < 2) continue;
                if (s.startsWith("#")) continue;
                //println("buildOutTree, " + s);
                processLine(s);
                //if (count++ > 100) break;
            }
            br.close();
    /*
    } catch(Exception e) {
        	println("buildOutTree, exception");
            e.printStackTrace();
            System.exit(0);
        }
        */
    }

    private void processLine(String s) throws Exception {
        //if (s.startsWith("Caa")) println("processLine, " + s);
        String [] sa = s.split(",");
        String cn = sa[1];
        if (cn.equals("Caaaaa")) {
            //println("processLine, " + s);
        }
        if (cn.startsWith("N")) return; //Nucs are still in there somehow

        int ts = Integer.parseInt(sa[STARTTIME]);
        iEndTime = Math.max(iEndTime, ts);
        //println("Vembryo4.processLine, " + iEndTime + CS + ts + CS + s);
        //if (iEndTime >= 304) System.exit(1);
        //if (ts > iEndTime) return;
        Cell c = (Cell)iCellsByName.get(cn);
        boolean existingCell = (c != null);
        if (existingCell) {
            int time = c.getTime();
            if (time == 0) {
                c.setTime(ts);
            }
        } else {
            c = new Cell(cn, Cell.LARGEENDTIME, ts);
        }
        Cell p = processParent(c);
        if (p == null) {
            //println("processLine error, " + c.getName());
            return;
        }
        if (!existingCell) {
            p.add(c);
            c.setEndFate(Cell.DIED);
            iCellsByName.put(cn, c);
        }
        //println("processLine, " + c.getName());
        addCellData(c, sa, ts);


        //Cell p = (Cell)iCellsByName.get(pn);
        //p.add(c);
        //int v = 0;
    }

    public void showCellsVsTime() {
    	int max = Math.min(300, nuclei_record.size());
    	for (int i=230; i <= 260; i+=1) {
    		Vector v = (Vector)nuclei_record.get(i - 1);
    		//println("showCellsVsTime, sulston time=" + i + CS + v.size());
    		//for (int j=0; j < v.size(); j++) {
    		//	Nucleus n = (Nucleus)v.get(j);
    		//	println("showCellsVsTime, B, " + n.identity);
    		//}
    		println("showCellsVsTime, " + i + CS + v.size());
    	}
    }

    private static final int
         CELLNAME = 1
        ,STARTTIME = 2
        ,RAWRED = 3
        ,GLOBAL = 4
        ,LOCAL = 5
        ,BLOT = 6
        ,CROSS = 7
        ,Z = 8
        ,X = 9
        ,Y = 10
        ,SIZE = 11
        ,GWEIGHT = 12
        ;





    private void addCellData(Cell c, String [] sa, int time) {
        //println("addCellData, " + c.getName() + CS + time);
        String name = c.getName();
        c.iEndingIndex = time;
        c.setEndTime(time);
        Vector cdv = c.getCellData();
        Nucleus n = new Nucleus();
        try {
            n.identity = name;
            n.rweight = Integer.parseInt(sa[RAWRED]);
            n.rwraw = Integer.parseInt(sa[RAWRED]);
            n.rwcorr1 = n.rwraw - Integer.parseInt(sa[GLOBAL]);
            n.rwcorr2 = n.rwraw - Integer.parseInt(sa[LOCAL]);
            n.rwcorr3 = n.rwraw - Integer.parseInt(sa[BLOT]);
            n.rwcorr4 = n.rwraw - Integer.parseInt(sa[CROSS]);
            //n.rweight = n.rwraw - n.rwcorr4; just testing
            n.z = Float.parseFloat(sa[Z]);
            n.x = Integer.parseInt(sa[X]);
            n.y = Integer.parseInt(sa[Y]);
            n.size = Integer.parseInt(sa[SIZE]);
            n.weight = Integer.parseInt(sa[GWEIGHT]);
            CellData cd = new CellData(n);
            cdv.add(cd);
        } catch(ArrayIndexOutOfBoundsException aiob) {
            //println("addCellData, " + name + CS + iCsvPath);
            //aiob.printStackTrace();
            //System.exit(0);
        }
        //println("addCellData, " + n);

        if (time > 1) {
        	Vector nv = (Vector)nuclei_record.get(time - 1);
            nv.add(n);
        }
        iLastTime = Math.max(iLastTime,time);

    }

    private String parentName(String cn) {
        if (cn.startsWith("Z")) {
            int h = 3;
        }

        int w = cn.length();
        if (w == 1) {
            if (cn.equals("E")) return "EMS";
            else if (cn.equals("C")) return "P2";
            else return "P3";
        }
        char x = cn.charAt(w - 1);
        if (Character.isLowerCase(x)) return cn.substring(0, w - 1);
        else {
            if (cn.equals("P0")) return null;
            if (cn.equals("P1")) return "P0";
        	if (cn.equals("EMS")) return "P1";
            if (cn.equals("P2")) return "P1";
            if (cn.equals("P3")) return "P2";
            if (cn.equals("P4")) return "P3";
            if (cn.equals("MS")) return "EMS";
            if (cn.equals("AB")) return "P0";
        }
        return "P4";
    }

    private Cell processParent(Cell c) {
        String pn = parentName(c.getName());
        Cell p = null;
        if (pn != null) p = (Cell)iCellsByName.get(pn);
        if (p == null) {
            //println("Vembryo4.processParent, " + iCsvPath);
            //System.exit(0);
        }
        else {
        	p.setEndTime(c.getTime() - 1);
        	p.setEndFate(Cell.DIVIDED);
        }
        //p.add(c);
        return p;

    }

    private void buildFoundersX() {
        for (int i=0; i < FOUNDERS.length; i++) {
            String cn = FOUNDERS[i];
            Cell c = new Cell(cn, Cell.LARGEENDTIME, 0);

        }

    }

    public void buildFounders() {
        iCellsByName = new Hashtable();

        String cn = "P0";
        Cell c = new Cell(cn, Cell.LARGEENDTIME, 0);
        iRoot = c;
        iCellsByName.put(cn, c);
        cn = "AB";
        c = new Cell(cn, Cell.LARGEENDTIME, 0);
        iCellsByName.put(cn, c);
        iRoot.add(c);
        cn = "P1";
        c = new Cell(cn, Cell.LARGEENDTIME, 0);
        iCellsByName.put(cn, c);
        iRoot.add(c);

        Cell cr = c; // P1
        cn = "EMS";
        c = new Cell(cn, Cell.LARGEENDTIME, 0);
        iCellsByName.put(cn, c);
        cr.add(c);
        cn = "P2";
        c = new Cell(cn, Cell.LARGEENDTIME, 0);
        iCellsByName.put(cn, c);
        cr.add(c);

        cr = (Cell)iCellsByName.get("EMS");
        cn = "MS";
        c = new Cell(cn, Cell.LARGEENDTIME, 0);
        iCellsByName.put(cn, c);
        cr.add(c);
        cn = "E";
        c = new Cell(cn, Cell.LARGEENDTIME, 0);
        iCellsByName.put(cn, c);
        cr.add(c);

        cr = (Cell)iCellsByName.get("P2");
        cn = "C";
        c = new Cell(cn, Cell.LARGEENDTIME, 0);
        iCellsByName.put(cn, c);
        cr.add(c);
        cn = "P3";
        c = new Cell(cn, Cell.LARGEENDTIME, 0);
        iCellsByName.put(cn, c);
        cr.add(c);

        cr = (Cell)iCellsByName.get("P3");
        cn = "D";
        c = new Cell(cn, Cell.LARGEENDTIME, 0);
        iCellsByName.put(cn, c);
        cr.add(c);
        cn = "P4";
        c = new Cell(cn, Cell.LARGEENDTIME, 0);
        iCellsByName.put(cn, c);
        cr.add(c);

        // try for fix to series starting at ABal
        cr = (Cell)iCellsByName.get("AB");
        cn = "ABa";
        c = new Cell(cn, Cell.LARGEENDTIME, 0);
        iCellsByName.put(cn, c);
        cr.add(c);
        cn = "ABp";
        c = new Cell(cn, Cell.LARGEENDTIME, 0);
        iCellsByName.put(cn, c);
        cr.add(c);

        cr = (Cell)iCellsByName.get("ABa");
        cn = "ABal";
        c = new Cell(cn, Cell.LARGEENDTIME, 0);
        iCellsByName.put(cn, c);
        cr.add(c);
        cn = "ABar";
        c = new Cell(cn, Cell.LARGEENDTIME, 0);
        iCellsByName.put(cn, c);
        cr.add(c);

        cr = (Cell)iCellsByName.get("ABp");
        cn = "ABpl";
        c = new Cell(cn, Cell.LARGEENDTIME, 0);
        iCellsByName.put(cn, c);
        cr.add(c);
        cn = "ABpr";
        c = new Cell(cn, Cell.LARGEENDTIME, 0);
        iCellsByName.put(cn, c);
        cr.add(c);


    }


    private static final String [] FOUNDERS =
    {"P0", "P1", "P2", "P3", "P4", "AB", "EMS", "E", "MS", "C", "D"};

    private String isFounder(String cn) {
        for (int i=0; i < FOUNDERS.length; i++) {
            if (cn.equals(FOUNDERS[i])) return FOUNDERS[i];
        }
        return null;
    }


    public Vector getNucs(int time) {
    	return (Vector)nuclei_record.get(time - 1);
    }

    public Hashtable getCellsByNameHash() {
    	return iCellsByName;
    }




    public static void main(String [] args) {
    }

    private static void println(String s) {System.out.println(s);}
    private static void print(String s) {System.out.print(s);}
    private static final String CS = ", ";
    private static final String TAB = "\t";
    private static final DecimalFormat DF0 = new DecimalFormat("####");
    private static final DecimalFormat DF1 = new DecimalFormat("####.#");
    private static final DecimalFormat DF4 = new DecimalFormat("####.####");
    private static String fmt4(double d) {return DF4.format(d);}
    private static String fmt1(double d) {return DF1.format(d);}
    private static String fmt0(double d) {return DF0.format(d);}

}
