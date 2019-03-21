package sulston;

import org.rhwlab.snight.MeasureCSV;

public class Normalizer {

	String		iSeries;
	double		iTimeSlope;
	double		iTimeOffset;
	double		iXCenter;
	double		iYCenter;
	double		iZCenter;
	double		iMajor;
	double		iMinor;
	double		iAngle;
	double		iZSlope;
	int			iTime;
	double		iZPixRes;

	public Normalizer(String [] sa) {
		iZPixRes = 11.1;
		iSeries = sa[0];
		iTimeSlope = Double.parseDouble(sa[1]);
		iTimeOffset = Double.parseDouble(sa[2]);
		iXCenter = Double.parseDouble(sa[3]);
		iYCenter = Double.parseDouble(sa[4]);
		iMajor = Double.parseDouble(sa[5]);
		iMinor = Double.parseDouble(sa[6]);
		iAngle = Double.parseDouble(sa[7]);
		iZCenter = Double.parseDouble(sa[8]);
		iZSlope = Double.parseDouble(sa[9]);
		iTime = Integer.parseInt(sa[10]);
	}

	public Normalizer(MeasureCSV measureCSV, double zPixRes) {
		iZPixRes = zPixRes;
	    iSeries = measureCSV.get(MeasureCSV.att[MeasureCSV.SERIES]);
		iTimeSlope = Double.parseDouble(measureCSV.get(MeasureCSV.att[MeasureCSV.TSLOPE]));
		iTimeOffset = Double.parseDouble(measureCSV.get(MeasureCSV.att[MeasureCSV.TINTERCEPT]));
		iXCenter = Double.parseDouble(measureCSV.get(MeasureCSV.att[MeasureCSV.EXCENTER]));
		iYCenter = Double.parseDouble(measureCSV.get(MeasureCSV.att[MeasureCSV.EYCENTER]));
		iMajor = Double.parseDouble(measureCSV.get(MeasureCSV.att[MeasureCSV.EMAJOR]));
		iMinor = Double.parseDouble(measureCSV.get(MeasureCSV.att[MeasureCSV.EMINOR]));
		iAngle = Double.parseDouble(measureCSV.get(MeasureCSV.att[MeasureCSV.EANG]));
		iZCenter = Double.parseDouble(measureCSV.get(MeasureCSV.att[MeasureCSV.ZCENTER]));
		iZSlope = Double.parseDouble(measureCSV.get(MeasureCSV.att[MeasureCSV.ZSLOPE]));
		iTime = Integer.parseInt(measureCSV.get(MeasureCSV.att[MeasureCSV.TIME]));


	}

	public String toString() {
		StringBuffer sb = new StringBuffer();
		sb.append(iSeries);
		sb.append(C + iTimeSlope);
		sb.append(C + iTimeOffset);
		sb.append(C + iXCenter);
		sb.append(C + iYCenter);
		sb.append(C + iMajor);
		sb.append(C + iMinor);
		sb.append(C + iAngle);
		sb.append(C + iZCenter);
		sb.append(C + iZSlope);
		sb.append(C + iTime);

		return sb.toString();
	}

	static final String C = ",";


	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
