package hiq.fdtd.genetic;

import hiq.fdtd.fdtd.FDTD;
import hiq.fdtd.genetic.util.Helper;
import hiq.fdtd.model.PhotonicMikaelyanLens;

import java.util.List;
import java.util.ArrayList;

public class Chromo {

	private static final double MUT_POSS = 0.1; // mutation possibility

	public int M; // number of encoded features
	public int u;
	public int N; // total number of bits

	private int[] chromo = new int[N];

	/**
	 * @return the chromo
	 */
	public int[] getChromo() {
		return chromo;
	}

	/**
	 * @param chromo
	 *            the chromo to set
	 */
	public void setChromo(int[] chromo) {
		this.chromo = chromo;
	}

	public List<Integer> getEntity() {

		List<Integer> radiuses = new ArrayList<Integer>();
		for (int i = 0; i < M; i++) {
			int[] g = new int[u];
			for (int j = 0; j < u; j++) {
				g[j] = chromo[i * u + j];
			}
			radiuses.add(Helper.grayToDecimal(g));
		}
		return radiuses;
	}

	public void mutation() {
		for (int i = 0; i < N; i++) {
			double r = Math.random();
			if (r < MUT_POSS) {
				chromo[i] = (int) Math.abs(chromo[i] - 1);
			}
		}
	}

	public String toString() {
		return super.toString();
	}

	public static Chromo getChromo(PhotonicMikaelyanLens lens, FDTD fdtd) {
		List<Double> radiuses = lens.getRadiuses();
		Chromo crm = new Chromo();

		if (radiuses.size() % 2 == 0) {
			crm.M = radiuses.size() / 2;
		} else {
			crm.M = (radiuses.size() + 1) / 2;
		}
		crm.u = Helper.binaryView(fdtd.fl(lens.getD())).length;
		crm.N = crm.M * crm.u;
		int[] crom = new int[crm.N];
		for (int i = 0; i < crm.M; i++) {
			int[] f = Helper.decimalToGray(fdtd.fl(radiuses.get(i)), crm.u);
			for (int j = 0; j < crm.u; j++) {
				crom[j + i * crm.u] = f[j];
			}

		}
		System.out.println("u = " + crm.u);
		System.out.println("M = " + crm.M);
		crm.setChromo(crom);
		return crm;
	}
}
