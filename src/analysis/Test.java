package analysis;

import java.text.DecimalFormat;

import flanagan.analysis.Stat;

public class Test {

	void test() {
	}

	public static String ttest(double [] da1, double [] da2) {
		return null;
	}

	void rdata(double [] da, double mu, double sig) {
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		println("Test.main, ");
		Test t = new Test();
		t.test();
		//t.rdata();

	}
    private static void println(String s) {System.out.println(s);}
    private static void print(String s) {System.out.print(s);}
    private static final String CS = ", ", C = ",", Q = "\"";
    private static final String TAB = "\t";
    private static final DecimalFormat DF0 = new DecimalFormat("####");
    private static final DecimalFormat DF1 = new DecimalFormat("####.#");
    private static final DecimalFormat DF4 = new DecimalFormat("####.####");
    private static String fmt4(double d) {return DF4.format(d);}
    private static String fmt1(double d) {return DF1.format(d);}
    private static String fmt0(double d) {return DF0.format(d);}
    private static final DecimalFormat DFE = new DecimalFormat("0.###E0");
    private static String fmtE(double d) {return DFE.format(d);}

}
