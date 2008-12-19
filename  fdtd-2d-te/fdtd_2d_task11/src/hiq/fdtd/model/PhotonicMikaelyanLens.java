package hiq.fdtd.model;

import java.util.ArrayList;
import java.util.List;

public class PhotonicMikaelyanLens {
	private int N;

	private List<Double> radiuses; // size of N
	private List<Double> lConstants; // size of N

	private double n0_mik;
	private double l_mik;

	/** length of PC lens */
	private double l;
	/** lattice constant */
	private double d;
	/** PC refractive index */
	private double n_pc;
	/** height of PC lens */
	private double b;
	private double z0;

	public double getZ0() {
		return z0;
	}

	public void setZ0(double z0) {
		this.z0 = z0;
	}

	private double mik_refrIndex(double y) {
		double res;
		res = n0_mik / Math.cosh((Math.PI * y) / (2 * l_mik));
		return res;
	}

	private double getRadii(double y) {
		double res;

		// first solution
		double temp1 = d / (2 * (n_pc - 1));
		double temp2 = n_pc - (l_mik / l) * mik_refrIndex(y);
		res = temp1 * temp2;

		// second solution

		// double temp1 = (n_pc - mik_refrIndex(y)) / (n_pc - 1) / Math.PI;
		// res = d * Math.sqrt(temp1);

		return res;
	}

	public PhotonicMikaelyanLens(double n0_mik, double l_mik, double z0,
			double l, double d, double n_pc, double b) {
		super();
		this.n0_mik = n0_mik;
		this.l_mik = l_mik;
		this.l = l;
		this.d = d;
		this.n_pc = n_pc;
		this.b = b;
		this.z0 = z0;

		generateLens();
	}

	public void generateLens() {
		this.N = (int) Math.floor(b / d);
		double y_part = (b - N * d) / 2;
		this.radiuses = new ArrayList<Double>();
		this.lConstants = new ArrayList<Double>();
		double y_mid = b / 2;
		for (int i = 0; i < N; i++) {
			double y = y_part + d / 2 + i * d;
			double radii = this.getRadii(Math.abs(y - y_mid));

			lConstants.add(d);
			radiuses.add(radii);
		}
	}

	public int getN() {
		return N;
	}

	public void setN(int n) {
		N = n;
	}

	public double getN_pc() {
		return n_pc;
	}

	public void setN_pc(double n_pc) {
		this.n_pc = n_pc;
	}

	public double getL() {
		return l;
	}

	public void setL(double l) {
		this.l = l;
	}

	public double getY0(int i) {
		double y_part = (b - N * d) / 2;
		return (y_part + d / 2 + i * d);
	}

	public double getRadii(int i) {
		return this.radiuses.get(i);
	}

	public double getLConstant(int i) {
		return this.lConstants.get(i);
	}

	public List<Double> getRadiuses() {
		return radiuses;
	}

	public void setRadiuses(List<Double> radiuses) {
		this.radiuses = radiuses;
	}

	public double getD() {
		return d;
	}

	public void setD(double d) {
		this.d = d;
	}
	
	public String toString(){
		
		String res = "";
		
		res += ("l = " + l + "; ");
		res += ("d = " + d + "; ");
		res += ("n_pc = " + n_pc + "; ");
		res += ("b = " + b + "; ");
		res += ("z0 = " + z0 + "; ");
		
		res += ("n0_mik = " + n0_mik + "; ");
		res += ("l_mik = " + l_mik + "; ");
		
		res += "RADIUSES (nm): ";
		for (Double d: this.getRadiuses()){
			res += (d * 1e9 + "; \n");
		}
		return res;
	}

}
