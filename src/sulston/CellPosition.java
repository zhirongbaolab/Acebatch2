package sulston;

import org.rhwlab.snight.DivisionCaller;
import org.rhwlab.snight.MeasureCSV;

public class CellPosition {
	public String		iSeries;
	public String		iName;
	public int			iX;
	public int			iY;
	public int			iZ;
	public int			iNX;
	public int			iNY;
	public int			iNZ;
	public int			iDia;

	public CellPosition() {

	}

	public CellPosition(int x, int y, int z) {
		// we work with z in pixels in this class
		this();
		iX = x;
		iY = y;
		iZ = z;
	}

	public CellPosition(String [] sa) {
		iX = Integer.parseInt(sa[0]);
		iY = Integer.parseInt(sa[1]);
		iZ = Integer.parseInt(sa[2]);
		iNX = Integer.parseInt(sa[3]);
		iNY = Integer.parseInt(sa[4]);
		iNZ = Integer.parseInt(sa[5]);
		iDia = Integer.parseInt(sa[6]);
		iName = sa[7];
		iSeries = sa[8];

	}

	public void normalize(Normalizer norm) {
		double ddelx = iX - norm.iXCenter;
		double ddely = iY - norm.iYCenter;
		double ddelz = iZ - norm.iZCenter * norm.iZPixRes;
		double ang = Math.toRadians(-norm.iAngle);
		double [] u = DivisionCaller.handleRotation(ddelx, ddely, ang);
		ddelx = u[0];
		ddely = u[1];
		double dMajor = Double.parseDouble(MeasureCSV.defaultAtt[MeasureCSV.EMAJOR]);
		double dMinor = Double.parseDouble(MeasureCSV.defaultAtt[MeasureCSV.EMINOR]);
		double dSlope = Double.parseDouble(MeasureCSV.defaultAtt[MeasureCSV.ZSLOPE]);
		double dxc = Double.parseDouble(MeasureCSV.defaultAtt[MeasureCSV.EXCENTER]);
		double dyc = Double.parseDouble(MeasureCSV.defaultAtt[MeasureCSV.EYCENTER]);
		double dzc = Double.parseDouble(MeasureCSV.defaultAtt[MeasureCSV.ZCENTER]);
		ddelx *= dMajor/norm.iMajor;
		ddely *= dMinor/norm.iMinor;
		ddelz *= dSlope/norm.iZSlope;

		iNX = (int)Math.round(ddelx + dxc);
		iNY = (int)Math.round(ddely + dyc);
		iNZ = (int)Math.round(ddelz + dzc * norm.iZPixRes);

	}


	public void convertOrientation(Normalizer norm, String axis, double zpixres) {
		if (axis.equals("ADL")) return;

		else {
			int z = iNZ;
			int y = iNY;
			int x = iNX;
			//int xcenter = (int)Math.round(norm.iXCenter);
			//int ycenter = (int)Math.round(norm.iYCenter);
			//int zcenter = (int)Math.round(norm.iZCenter * zpixres);
			int xcenter = (int)Math.round(Double.parseDouble(MeasureCSV.defaultAtt[MeasureCSV.EXCENTER]));
			int ycenter = (int)Math.round(Double.parseDouble(MeasureCSV.defaultAtt[MeasureCSV.EYCENTER]));
			int zcenter = (int)Math.round(Double.parseDouble(MeasureCSV.defaultAtt[MeasureCSV.ZCENTER]) * zpixres);
			axis = axis.toUpperCase();
			if (axis.equals("AVR")) {
				y = ycenter - (y - ycenter);
				z = zcenter - (z - zcenter);
			} else if (axis.equals("PDR")) {
				x = xcenter - (x - xcenter);
				z = zcenter - (z - zcenter);
			} else if (axis.equals("PVL")) {
				x = xcenter - (x - xcenter);
				y = ycenter - (y - ycenter);
			}
			iNX = x;
			iNY = y;
			iNZ = z;
		}

	}


	public String toString() {
		StringBuffer sb = new StringBuffer();
		sb.append(iX);
		sb.append(C + iY);
		sb.append(C + iZ);
		sb.append(C + iNX);
		sb.append(C + iNY);
		sb.append(C + iNZ);
		sb.append(C + iDia);
		sb.append(C + iName);
		sb.append(C + iSeries);

		return sb.toString();
	}

	static final String C = ",";

}
