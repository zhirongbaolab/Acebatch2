package analysis;

import java.text.DecimalFormat;

public class LeafReport {

	public String		iLeaf;
	public boolean		iExpressing;
	public boolean		iOnset;
	public int			iPeakExp;
	public String		iOnsetCell;
	public int			iOnsetTime;
	public double		iFractionLifetime;
	public String		iFate;
	public String		iGene;
	public String		iSeries;

	public LeafReport() {
	}

	public LeafReport(String s) {
		String [] sa = s.split(C);
		iLeaf = sa[LEAF];
		iExpressing = Boolean.parseBoolean(sa[EXPRESSING]);
		iOnset = Boolean.parseBoolean(sa[ONSET]);
		iPeakExp = Integer.parseInt(sa[PEAK]);
		iOnsetCell = sa[ONSETCELL];
		iOnsetTime = Integer.parseInt(sa[ONSETTIME]);
		iFractionLifetime = Double.parseDouble(sa[FRACTIONLIFETIME]);
		iFate = sa[FATE];
		iGene = sa[GENE];
		iSeries = sa[SERIES];

	}

	public String toString() {
		StringBuffer sb = new StringBuffer(NULL + iLeaf);
		sb.append(C + iExpressing);
		sb.append(C + iOnset);
		sb.append(C + iPeakExp);
		sb.append(C + iOnsetCell);
		sb.append(C + iOnsetTime);
		sb.append(C + fmt2(iFractionLifetime));
		sb.append(C + iFate);
		sb.append(C + iGene);
		sb.append(C + iSeries);
		return sb.toString();

	}

	public String toRList() {
		StringBuffer sb = new StringBuffer("c(");
		sb.append(rString(iLeaf));
		sb.append(C + rBoolean(iExpressing));
		sb.append(C + rBoolean(iOnset));
		sb.append(C + iPeakExp);
		sb.append(C + rString(iOnsetCell));
		sb.append(C + iOnsetTime);
		sb.append(C + fmt2(iFractionLifetime));
		sb.append(C + rString(iFate));
		sb.append(C + rString(iGene));
		sb.append(C + rString(iSeries));
		sb.append(")");

		return sb.toString();
	}

	String rBoolean(boolean b) {
		if (b) return "T";
		return "F";
	}

	String rString(String x) {
		return Q + x + Q;
	}

	public static final int
		 LEAF = 0
		,EXPRESSING = 1
		,ONSET= 2
		,PEAK = 3
		,ONSETCELL = 4
		,ONSETTIME = 5
		,FRACTIONLIFETIME = 6
		,FATE = 7
		,GENE = 8
		,SERIES = 9



		;



	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		//String s = "pha-4,ABalpaaaaa,e3VL,126,ABalpaa,1,20061215_pha4I2L_11,8109";
		String s = "ABalpaaaaa,true,true,8109,ABalpaa,126,0.65,e3VL,pha-4,20061215_pha4I2L_11";
		LeafReport lr = new LeafReport(s);
		println("rList, " + lr.toRList());
		println("main, " + lr);
		lr = new LeafReport();
		println("main, " + lr);


	}
    public static void println(String s) {System.out.println(s);}
    public static void print(String s) {System.out.print(s);}
    public static final String CS = ", ", C = ",", Q = "\"";
    public static final String TAB = "\t", NULL = "";
    public static final DecimalFormat DF0 = new DecimalFormat("####");
    public static final DecimalFormat DF1 = new DecimalFormat("####.#");
    public static final DecimalFormat DF2 = new DecimalFormat("####.##");
    public static final DecimalFormat DF4 = new DecimalFormat("####.####");
    public static final DecimalFormat DFE = new DecimalFormat("0.####E0");

    public static String fmt4(double d) {return DF4.format(d);}
    public static String fmt2(double d) {return DF2.format(d);}
    public static String fmt1(double d) {return DF1.format(d);}
    public static String fmt0(double d) {return DF0.format(d);}
    public static String fmtE(double d) {return DFE.format(d);}

}
