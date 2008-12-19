package hiq.fdtd.genetic;

import hiq.fdtd.fdtd.FDTD;
import hiq.fdtd.model.PhotonicMikaelyanLens;

import java.util.ArrayList;
import java.util.List;
import java.io.*;

public class Genetic {

	private static final double EPSILON = 0.001;
	private static final int AMOUNT = 8; // how many entities in population
	private static final int AGES = 100; // how many AGES
	private int iteration = 0;
	private int correct = 0;

	private double[][] story = new double[AGES][AMOUNT];
	private Chromo[][] storyPoints = new Chromo[AGES][AMOUNT];

	private double[] fields; // for roulet selection

	PhotonicMikaelyanLens lens;
	PhotonicMikaelyanLens resLens;

	public FDTD obj;

	public static void main(String[] args) {
		new Genetic();
	}

	public Genetic() {

		obj = new FDTD();

		this.lens = new PhotonicMikaelyanLens(1.5, 3e-6, obj.tf_z0_ + 0.5e-6,
				3e-6, 0.25e-6, 1.5, obj.tf_y_ - obj.tf_y0_);

		List<Chromo> sols = generateFirstPopulation();

		while (!stopCondition(sols)) {
			List<Chromo> ancs = new ArrayList<Chromo>();

			preRouletSelection(sols);
			for (int i = 0; i < AMOUNT / 2; i++) {
				List<Chromo> parents = rouletSelection(sols);
				ancs.addAll(makeChildren(parents));
			}

			// mutation(ancs);
			sols = ancs;

			System.out.println("ITERATION: " + ++iteration);
		}

		obj.exportAll();
		try {
			exportResult();
			exportStory();
		} catch (Exception e) {

		}

	}

	private void exportStory() throws FileNotFoundException, IOException {

		PrintWriter pw = new PrintWriter("story.txt");
		for (int i = 0; i < AGES; i++) {
			for (int j = 0; j < AMOUNT; j++) {
				pw.print(story[i][j] + "{" + storyPoints[i][j] + "}" + "; ");
			}
			pw.println();
		}

		pw.flush();
		pw.close();
	}

	private boolean stopCondition(List<Chromo> ancs) {

		boolean res = iteration >= AGES;
		return res;

	}

	private void mutation(List<Chromo> ancs) {
		for (Chromo c : ancs) {
			c.mutation();
		}
	}

	private List<Chromo> makeChildren(List<Chromo> parents) {
		int N = parents.get(0).N;
		int M = parents.get(0).M;
		int u = parents.get(0).u;

		int divisionPoint = (int) Math.floor(Math.random() * (N - 1)) + 1;

		int[] a = new int[N];
		int[] b = new int[N];

		int[] p1 = parents.get(0).getChromo();
		int[] p2 = parents.get(1).getChromo();

		for (int i = 0; i < divisionPoint; i++) {
			a[i] = p1[i];
			b[i] = p2[i];
		}

		for (int i = divisionPoint; i < N; i++) {
			a[i] = p2[i];
			b[i] = p1[i];
		}

		Chromo c1 = new Chromo();
		c1.N = N;
		c1.M = M;
		c1.u = u;
		c1.setChromo(a);

		/*
		 * Chromo c2 = new Chromo(); c2.N = N; c2.M = M; c2.u = u;
		 * c2.setChromo(b);
		 */
		Chromo c2 = new Chromo();
		c2.N = N;
		c2.M = M;
		c2.u = u;
		c2.setChromo(a);
		c2.mutation();

		List<Chromo> res = new ArrayList<Chromo>();
		res.add(c1);
		res.add(c2);

		return res;

	}

	private void preRouletSelection(List<Chromo> sols) {

		double[] squares = new double[AMOUNT];
		double summ = 0;
		double[] fitn = new double[AMOUNT];

		for (int i = 0; i < AMOUNT; i++) {
			fitn[i] = fitness(sols.get(i));
			story[iteration][i] = fitn[i];
			storyPoints[iteration][i] = sols.get(i);
			summ += fitn[i];
		}

		for (int i = 0; i < AMOUNT; i++) {
			squares[i] = fitn[i] / summ;
		}

		fields = new double[AMOUNT];

		fields[0] = squares[0];

		for (int i = 1; i < AMOUNT; i++) {
			fields[i] = fields[i - 1] + squares[i];
		}

	}

	private List<Chromo> rouletSelection(List<Chromo> sols) {

		List<Chromo> res = new ArrayList<Chromo>();

		double r1 = Math.random();
		double r2 = Math.random();

		int k = 0;
		int l = 0;

		for (int i = AMOUNT - 2; i > 0; i--) {
			if (r1 > fields[i]) {
				k = i + 1;
				break;
			}
		}

		for (int i = AMOUNT - 2; i > 0; i--) {
			if (r2 > fields[i]) {
				l = i + 1;
				break;
			}
		}

		System.out.println("roulet selection: " + k + ", " + l);

		res.add(sols.get(k));
		res.add(sols.get(l));

		return res;
	}

	/**
	 * This function must be higher if entity is more appropriate.
	 * 
	 * @param c
	 * @return
	 */
	private double fitness(Chromo c) {
		double res = 0;
		res = function(c);
		return (res * res * res * res);
	}

	private double function(Chromo c) {
		List<Double> rs = new ArrayList<Double>();
		List<Integer> list = c.getEntity();

		int l = lens.getRadiuses().size();
		if (l % 2 == 0) {
			for (int i = 0; i < list.size(); i++) {
				rs.add(list.get(i) * obj.ds);
			}
			for (int i = 0; i < list.size(); i++) {
				rs.add(list.get(list.size() - i - 1) * obj.ds);
			}
		}

		if (l % 2 == 1) {
			for (int i = 0; i < list.size(); i++) {
				rs.add(list.get(i) * obj.ds);
			}
			for (int i = 0; i < list.size() - 1; i++) {
				rs.add(list.get(list.size() - i - 2) * obj.ds);
			}
		}

		resLens = new PhotonicMikaelyanLens(1.5, 3e-6, obj.tf_z0_ + 0.5e-6,
				3e-6, 0.25e-6, 1.5, obj.tf_y_ - obj.tf_y0_);

		resLens.setRadiuses(rs);

		obj.init();
		obj.initPML();

		obj.setPhotonicMikLens(resLens);

		obj.run();
		double res = obj.maxValue(obj.fl(obj.tf_z0_), obj.fl(obj.tf_z_));
		System.out.println("MAX: "
				+ obj.maxValue(obj.fl(obj.tf_z0_), obj.fl(obj.tf_z_)));
		return res;

	}

	private List<Chromo> generateFirstPopulation() {
		List<Chromo> res = new ArrayList<Chromo>();
		for (int i = 0; i < this.AMOUNT; i++) {
			res.add(Chromo.getChromo(lens, obj));
		}
		return res;
	}

	/**
	 * 
	 * @throws FileNotFoundException
	 */
	private void exportResult() throws FileNotFoundException {
		PrintWriter pw = new PrintWriter("result/result.txt");
		pw.println("AMOUNT: " + AMOUNT);
		pw.println("AGES: " + AGES);
		pw.println("AMOUNT: " + AMOUNT);

		pw.println();

		pw.println(resLens);
		pw.flush();
		pw.close();
	}
}