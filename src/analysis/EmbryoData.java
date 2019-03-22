package analysis;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.Enumeration;
import java.util.Vector;

import org.rhwlab.snight.Config;
import org.rhwlab.snight.Nucleus;
import org.rhwlab.tree.Cell;
import org.rhwlab.tree.CellData;


public class EmbryoData {

    Embryo      iEmbryo;


    public EmbryoData(String configFile) {
        String name = (new File(configFile)).getName();
        iEmbryo = new Embryo(name, configFile);
    }

    public EmbryoData(Embryo embryo) {
        iEmbryo = embryo;
    }


    public Vector findLeavesOfInterest(String rootName, int endTime) {
        Vector leaves = new Vector();
        Cell root = (Cell)iEmbryo.iCellsByName.get(rootName);
        int k = root.getLeafCount();
        //println("rootName, " + k);
        int count = 0;
        Enumeration e = root.depthFirstEnumeration();
        while (e.hasMoreElements()) {
            Cell c = (Cell)e.nextElement();
            if (c.isLeaf()) {
                while(c.getTime() > endTime) {
                	c = (Cell)c.getParent();
                }
                if (leaves.contains(c)) continue;
                leaves.add(c);
            }
        }
        return leaves;
    }


    private static final String ROOTNAME = "P0";

    // gets a Vector of strings where the string is the image of a csv file line
    // and a file name
    private void exportData(Vector outVec, String fileName) {
    	println("exportData, " + fileName);
    	//if (1 == 1) return;
        FileOutputStream fos = null;
        PrintWriter pw = null;
        try {
            fos = new FileOutputStream(fileName);

        } catch(IOException ioe) {
            ioe.printStackTrace();
            return;
        }
        pw = new PrintWriter(fos);
        for (int i=0; i < outVec.size(); i++) {
            pw.println((String)outVec.get(i));
        }
        pw.close();
    }

    // note that this function does gets much more than red data
    // the name goes back to its earlier embodiment before we added
    // the rest of the data into the csv files
    // also type will normally be "all" so all the red data fields will be included
    //
    // note that z positions are coded: a z of 15.2 is given as an integer 152
    private Vector getRedDataV(Nucleus n, String type) {
        Vector v = new Vector();
        if (!type.equals("all")) {
            v.add(new Integer(n.getCorrectedRed(type)));
        } else {
            int ntypes = Config.REDCHOICE.length;
            for (int i=0; i < ntypes; i++) {
                v.add(new Integer(n.getCorrectedRed(Config.REDCHOICE[i])));
            }
        }
        v.add(new Integer(Math.round(10 * n.z))); //coding z as int of value * 10
        v.add(new Integer(n.x));
        v.add(new Integer(n.y));
        v.add(new Integer(n.size));
        v.add(new Integer(n.weight));
        return v;
    }

    // provides the csv representation of what is in the vector
    // that is a string with fields separated by commas
    // where each field is the string representation of the data
    private String redV2S(Vector v) {
        int zloc = 5; //where to find z in the Vector
        StringBuffer sb = new StringBuffer();
        sb.append(C + ((Integer)v.get(0)).intValue());
        for (int i=1; i < v.size(); i++) {
            if (i == zloc) {
                // this is coded z plane
                double z = (double)((Integer)v.get(i)).intValue();
                sb.append(C + z/10);
            } else {
                sb.append(C + ((Integer)v.get(i)).intValue());
            }
        }
        return sb.toString();
    }


    public void makeVector3(String rootName, int endTime, String type, String fileName, String leavesName) {
    	boolean debug = false;
        Vector outVec = new Vector();
        Vector outAvgVec = new Vector();
        Cell root = (Cell)iEmbryo.iCellsByName.get(rootName);
        Vector leaves = findLeavesOfInterest(rootName, endTime);
        //Vector rootLeaves = findLeavesOfInterest(rootName, endTime);
        //rootLeaves.add(0, root);
        String header = "cellTime,cell,time";
        if (!type.equals("all")) header += C + type;
        else {
            header += C + "none" + C + "global" + C + "local" + C + "blot" + C + "cross";
            header += C + "z" + C + "x" + C + "y" + C + "size" + C + "gweight";
        }
        int zloc = 5; //where to find z in the Vector
        outVec.add(header);
        outAvgVec.add(header);
        int leafCount = 0;
        Nucleus nlast = null;
        int tlast = 0;
        Enumeration e = root.preorderEnumeration();
        while (e.hasMoreElements()) {
            Vector va = new Vector(); //to hold all the Vectors of red data
            Cell c = (Cell)e.nextElement();
            debug = false;
            //if (c.getName().startsWith("ABplp")) debug = true;
            //println("makeVector3, " + c.getName());
            int tstart = c.getTime();
            Vector cdv = c.getCellData();
            String cellTime = "";
            String s = "";
            for (int k=0; k < cdv.size(); k++) {
                tlast = tstart + k;
                if (tlast > endTime) break;
                CellData cd = (CellData)cdv.get(k);
                Nucleus n = cd.iNucleus;
                nlast = n;
                cellTime = c.getName() + ":" + (tstart + k);
                s = cellTime + C + c.getName() + C + (tstart + k);
                Vector v = getRedDataV(n, type);
                s = s + redV2S(v);


                if (debug) {
                	int r = (Integer)v.get(3);
                	println("makeVector3, " + r + CS + s);
                }
                //s = s + C + n.z;
                outVec.add(s);
                va.add(v);
                //println(s);
                //println("" + n);

            }
            // now generate the line of cell averaged data
            int numPoints = va.size();
            if (numPoints > 0) {
                //int time = (tstart + tlast)/2;
                int time = tstart;
                cellTime = c.getName() + ":" + time;
                s = cellTime + C + c.getName() + C + time;
                Vector v = (Vector)va.get(0);
                int numVals = v.size();
                int [] avg = new int[numVals];
                for (int i=0; i < numVals; i++) {
                    avg[i] = 0;
                }
                for (int j=0; j < va.size(); j++) {
                    v = (Vector)va.get(j);
                    for (int i=0; i < numVals; i++) {
                        avg[i] += ((Integer)v.get(i)).intValue();
                    }
                }
                for (int i=0; i < numVals; i++) {
                    if (i == zloc) {
                        // this is coded z plane
                        double z = (double)(Math.round(avg[i]/numPoints));
                        s += C + z/10;
                    } else {
                        s += C + ((int)Math.round(avg[i]/numPoints));
                    }
                }
                //println(s);
                outAvgVec.add(s);
            }
        }
        exportData(outVec, fileName);
        fileName = fileName.replaceFirst("CD", "CA");
        exportData(outAvgVec, fileName);
        Vector leafNames = new Vector();
        for (int i=0; i < leaves.size(); i++) {
            String s = ((Cell)leaves.get(i)).getName();
            leafNames.add(s);
        }
        exportData(leafNames, leavesName);

        // now generate and output the root to leaves paths
        Vector [] paths = new Vector[leaves.size()];
        for (int i=0; i < leaves.size(); i++) {
            paths[i] = new Vector();
            Cell cLeaf = (Cell)leaves.get(i);
            e = cLeaf.pathFromAncestorEnumeration(root);
            while (e.hasMoreElements()) {
                Cell ca = (Cell)e.nextElement();
                paths[i].add(ca.getName());
                //println(ca.getName());
            }
        }
        //println("\n\n\n");
        boolean done = false;
        //while (!done) {
            StringBuffer sb = new StringBuffer();
            done = true;
            for (int i=0; i < paths.length; i++) {
                Vector v = paths[i];
                sb = new StringBuffer();
                while (v.size() > 0) {
                    String s = "";
                    try {
                        s = (String)v.remove(0);
                        done = false;
                    } catch(ArrayIndexOutOfBoundsException aiobe) {
                    }
                    sb.append(s);
                    if (v.size() > 0) sb.append(CS);
                }
            }
    }

    public void arraySpreadSheet(String fileName, String type) {
        println("arraySpreadSheet, " + fileName + CS + type);
        Cell root = (Cell)iEmbryo.iCellsByName.get(ROOTNAME);
        int endTime = iEmbryo.getEditedTimePoints();
        File f = new File(fileName);
        String parent = f.getParent();
        String leafFileName = "";

        if (parent == null) leafFileName = "L" + f.getName();
        else leafFileName = parent + "/L" + f.getName();

        makeVector3(root.getName(), endTime, type, fileName, leafFileName);
    }

    private void leavesSpreadSheet(String fileName, String type) {
        println("leavesSpreadSheet, " + fileName + CS + type);
        //if (1 == 1) return;
        Cell root = (Cell)iEmbryo.iCellsByName.get(ROOTNAME);
        int endTime = iEmbryo.getEditedTimePoints();
        Vector leaves = findLeavesOfInterest(root.getName(), endTime);
        Vector [] paths = new Vector[leaves.size()];
        Cell cLeaf = null;
        int jlast = 0;
        for (int i=0; i < leaves.size(); i++) {
            jlast = 0;
            paths[i] = new Vector();
            paths[i].add(iEmbryo.getName());
            cLeaf = (Cell)leaves.get(i);
            paths[i].add(cLeaf.getName());
            Enumeration e = cLeaf.pathFromAncestorEnumeration(root);
            while (e.hasMoreElements()) {
                Cell ca = (Cell)e.nextElement();
                //println(ca.getName());
                //paths[i].add(ca.getName());
                Vector cd = ca.getCellData();
                for (int j=0; j < cd.size(); j++) {
                    if (++jlast > endTime) break;
                    CellData cdo = (CellData)cd.get(j);
                    Nucleus n = cdo.iNucleus;
                    paths[i].add(String.valueOf(n.getCorrectedRed(type)));
                    //println("leavesSpreadSheet, " + n.rwraw);
                }
            }
            //println("leavesSpreadSheet, " + i + CS + jlast + CS + cLeaf.getName());
        }
        Vector outVec = new Vector();
        for (int i=0; i < paths.length; i++) {
            String s = dataLine(paths[i]);
            outVec.add(s);
        }
        exportData(outVec, fileName);

    }


    private String dataLine(Vector v) {
        StringBuffer sb = new StringBuffer();
        sb.append(v.get(0));
        for (int i=1; i < v.size(); i++) {
            sb.append(C + (String)v.get(i));
        }
        return sb.toString();
    }



    public static void showArgs(String [] args) {
        if (args == null) return;
        StringBuffer sb = new StringBuffer("showArgs: ");
        for (int i=0; i < args.length; i++) {
            sb.append(args[i] + SP);
        }
        println(sb.toString());

    }

    public static void spreadsheet1(String [] args) throws Exception {
        showArgs(args);
        String e = "";
        String p = "";
        String type = "none";
        e = args[0];
        type = args[1];
        p = args[2];
        //if (p.charAt(p.length() - 1) != '/') p += "/";

        //if (args.length > 2) type = args[2];
        //String ext = ".csv";
        Embryo embryo = new Embryo(args[0], args[3]);
        embryo.setEditedTimePoints(Integer.parseInt(args[4]));
        EmbryoData ed = new EmbryoData(embryo);
        ed.leavesSpreadSheet(p, type);

    }


    public static void spreadsheet2(String [] args) throws Exception {
        //acebatch2.org.rhwlab.analysis.EmbryoData
    	showArgs(args);
        Embryo embryo = new Embryo(args[0], args[3]);
        embryo.setEditedTimePoints(Integer.parseInt(args[4]));
        EmbryoData ed = new EmbryoData(embryo);
        ed.arraySpreadSheet(args[2], args[1]);

    }

    // needs 5 arguments
    // 0 = series name
    // 1 = "all"
    // 2 = CD..csv file name
    // 3 = acetree config file path
    // 4 = number of time points
    public static void main(String [] args) {
    	try {
    		spreadsheet2(args);
    	} catch(Exception e) {
    		e.printStackTrace();
    	}
    }



    private static void println(String s) {System.out.println(s);}
    private static void print(String s) {System.out.print(s);}
    private static final String CS = ", ", C = ",", SP=" ";
    private static final DecimalFormat DF0 = new DecimalFormat(" 0000");
    private static final DecimalFormat DF1 = new DecimalFormat("####.#");
    private static final DecimalFormat DF2 = new DecimalFormat("00.00");
    private static final DecimalFormat DF4 = new DecimalFormat("####.####");
    private static final DecimalFormat DF9 = new DecimalFormat("000000000");
    private static String fmt4(double d) {return DF4.format(d);}
    private static String fmt1(double d) {return DF1.format(d);}
    private static String fmt0(double d) {return DF0.format(d);}

}
