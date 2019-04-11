package vembryo;

import org.rhwlab.snight.NucleiMgr;
import org.rhwlab.snight.Nucleus;
import org.rhwlab.tree.Cell;
import org.rhwlab.tree.CellData;

import java.io.*;
import java.net.URL;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Vector;

public class Vembryo3 {

    String          iCsvPath;
    Cell            iRoot;
    Hashtable       iCellsByName;
    Vector          iRoots;
    int             iEndTime;
    int				iLastTime;
    Vector			nuclei_record;
    NucleiMgr		iNucleiMgr;

    public Vembryo3(String csvPath, boolean fake, boolean resource) {
        iCsvPath = csvPath;
        initNucRec();
        buildFounders();
        buildOutTree(resource);
        if (fake) {
        	addFakeData();
        	addSuccessorInfo();
        }

    }

    public Vector getNucleiRecord() {
    	return nuclei_record;
    }

    public void initNucRec() {
		iEndTime = 600;
		iLastTime = 0;
		nuclei_record = new Vector();
		for (int i=1; i <= iEndTime; i++) {
			nuclei_record.add(i - 1, new Vector());
		}
		iNucleiMgr = new NucleiMgr();
		iNucleiMgr.setNucleiRecord(nuclei_record);
    }

    public void buildTreeFromResourceFile() {
        System.out.println("Building tree from resource file: SCD20081128_sulston.csv");
	    URL url = this.getClass().getResource("/sulston/SCD20081128_sulston.csv");

	    InputStream istream = null;
	    try {
	        istream = url.openStream();
	        BufferedReader br = new BufferedReader(new InputStreamReader(istream));
	        String s;
	        br.readLine(); // toss the header
	        while (br.ready()) {
	            s = br.readLine();
                if (s.length() < 2) continue;
                if (s.startsWith("#")) continue;
                processLine(s);
	        }
	        br.close();
	    } catch(Exception e) {
	        e.printStackTrace();
	    }


    }

    public void buildOutTree(boolean resource) {
    	if (resource) {
    		buildTreeFromResourceFile();
    		return;
    	}
        try {
            FileInputStream fis = new FileInputStream(iCsvPath);
            BufferedReader br = new BufferedReader(new InputStreamReader(fis));
            br.readLine(); //toss the header line
            while (br.ready()) {
                String s = br.readLine();
                if (s.length() < 2) continue;
                if (s.startsWith("#")) continue;
                processLine(s);
            }
            br.close();
        } catch(Exception e) {
            e.printStackTrace();
        }
    }

    void addFakeData() {
        Cell c = (Cell)iRoot;
        Enumeration e = c.preorderEnumeration();
        while (e.hasMoreElements()) {
        	// treat each cell like a parent wholes successor indices we will assign
            Cell cc = (Cell)e.nextElement();
            Vector vp = cc.getCellData();
            for (int i=0; i < vp.size(); i++) {
            	CellData cd = (CellData)vp.get(i);
            	Nucleus nuc = cd.iNucleus;
            	nuc.x = 350;
            	nuc.y = 250;
            	nuc.z = 15;
            	nuc.size = 20;
            }
        }

    }

    void addSuccessorInfo() {
    	//println("addSuccessorInfo, *********************************");
        Cell c = (Cell)iRoot;
        Enumeration e = c.preorderEnumeration();
        while (e.hasMoreElements()) {
        	// treat each cell like a parent wholes successor indices we will assign
            Cell cc = (Cell)e.nextElement();
            int endTime = cc.getEndTime();
            //if (endTime < 52) println("addSuccessorInfo, " + cc.getName() + CS + cc.getTime() + CS + endTime);

            Vector vp = cc.getCellData();
            int size = vp.size();
            if (size == 0) continue;

            // first pretend there are no divisions so we propagate forward
            CellData cdB = (CellData)vp.get(0);
            Nucleus nucNow = cdB.iNucleus;
            for (int i=1; i < vp.size(); i++) {
            	CellData cd = (CellData)vp.get(i);
            	Nucleus nucNext = cd.iNucleus;
            	nucNow.successor1 = nucNext.index;
            	nucNow.successor2 = Nucleus.NILLI;
            	//if (endTime < 17) println("addSuccessorInfo, " + i + CS + vp.size() + CS + nucNow);
            	nucNow = nucNext;

            }

            CellData cdp = (CellData)vp.get(vp.size() - 1);
            Nucleus p = cdp.iNucleus;
        	//if (endTime < 17) println("addSuccessorInfo, bbb, " + (vp.size() - 1) + CS + p);
            int k = cc.getChildCount();
            if (k < 2) continue;

            Cell d0 = (Cell)cc.getChildAt(0);
            Vector vd0 = d0.getCellData();
            CellData cd0 = (CellData)vd0.get(0);
            Nucleus nd0 = cd0.iNucleus;
            p.successor1 = nd0.index;

            Cell d1 = (Cell)cc.getChildAt(1);
            Vector vd1 = d1.getCellData();
            CellData cd1 = (CellData)vd1.get(0);
            Nucleus nd1 = cd1.iNucleus;
            p.successor2 = nd1.index;
        }


    }

    private void processLine(String s) throws Exception {
        //println("processLine, " + s);
        String [] sa = s.split(",");
        String cn = sa[1];
        if (cn.equals("P0")) return;
        //if (cn.equals("D")) {
        //    println("processLine, " + s);
        //}
        if (cn.startsWith("N")) return; //Nucs are still in there somehow

        int ts = Integer.parseInt(sa[STARTTIME]);
        if (ts > iEndTime) return;

        Cell c = (Cell)iCellsByName.get(cn);
        boolean newCell = (c == null) || c.getTime() == 0;
        if (c == null) c = new Cell(cn, Cell.LARGEENDTIME, ts);
        if (newCell) c.setTime(ts);

        Cell p = processParent(c);
        if (p == null) {
            println("processLine error, " + c.getName());
            return;
        }
        if (newCell) {
        	//println("processLine, " + c.getName() + CS + p.getName());
            p.add(c);
            p.setEndFate(Cell.DIVIDED);
            //p.setEndTime(time)
            c.setEndFate(Cell.DIED);
            iCellsByName.put(cn, c);
        }
        //println("processLine, " + c.getName());
        if (ts > 0) addCellData(c, sa, ts);

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
    	Vector nv = (Vector)nuclei_record.get(time - 1);
        int index = nv.size() + 1;

        String name = c.getName();
        c.iEndingIndex = time;
        c.setEndTime(time);
        Vector cdv = c.getCellData();
        int predecessor = -1;
        if (time > 1) predecessor = getPredecessor(c);
        Nucleus n = new Nucleus();
        try {
        	n.index = index;
        	n.predecessor = predecessor;
        	n.status = 1;
            n.identity = name;
            n.rweight = Integer.parseInt(sa[RAWRED]);
            n.rwraw = Integer.parseInt(sa[RAWRED]);
            n.rwcorr1 = n.rwraw - Integer.parseInt(sa[GLOBAL]);
            n.rwcorr2 = n.rwraw - Integer.parseInt(sa[LOCAL]);
            n.rwcorr3 = n.rwraw - Integer.parseInt(sa[BLOT]);
            n.rwcorr4 = n.rwraw - Integer.parseInt(sa[CROSS]);
            n.z = Float.parseFloat(sa[Z]);
            n.x = Integer.parseInt(sa[X]);
            n.y = Integer.parseInt(sa[Y]);
            n.size = Integer.parseInt(sa[SIZE]);
            n.weight = Integer.parseInt(sa[GWEIGHT]);
            CellData cd = new CellData(n);
            cdv.add(cd);
            nv.add(n);

        } catch(ArrayIndexOutOfBoundsException aiob) {
            println("addCellData, " + name + CS + iCsvPath);
            aiob.printStackTrace();
            System.exit(0);
        }

        iLastTime = Math.max(iLastTime,time);

    }

    private int getPredecessor(Cell c) {
    	Vector cdv = c.getCellData();
    	CellData cd = null;
    	int size = cdv.size();
    	if (size == 0) {
    		Cell p = (Cell)c.getParent();
    		cdv = p.getCellData();
    		size = cdv.size();
    		if (size == 0) return -1;
    	}
    	cd = (CellData)cdv.get(size - 1);
    	if (cd != null) {
    		return cd.iNucleus.index;

    	}

    	return 0;
    }

    private String parentName(String cn) {
        if (cn.startsWith("Z")) {
            int h = 3;
        }

        int w = cn.length();

        char x = 'X';
        if (w > 1) x = cn.charAt(w - 1);
        if (Character.isLowerCase(x)) return cn.substring(0, w - 1);
        else {
            if (cn.equals("AB")) return "P0";
        	if (cn.equals("P1")) return "P0";
            if (cn.equals("EMS")) return "P1";
            if (cn.equals("P2")) return "P1";
            if (cn.equals("P3")) return "P2";
            if (cn.equals("C")) return "P2";
            if (cn.equals("P4")) return "P3";
            if (cn.equals("D")) return "P3";
            if (cn.equals("MS")) return "EMS";
            if (cn.equals("E")) return "EMS";
        }
        return "P4";
    }

    private Cell processParent(Cell c) {
        String pn = parentName(c.getName());
        Cell p = (Cell)iCellsByName.get(pn);
        p.setEndTime(c.getTime() - 1);
        //p.add(c);
        return p;

    }

    private void buildFounders() {
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

    private static void println(String s) {System.out.println(s);}
    private static final String CS = ", ", C = ",";
}