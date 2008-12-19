package hiq.fdtd.fdtd;

import hiq.fdtd.model.PhotonicMikaelyanLens;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.PrintWriter;
import java.util.Calendar;
import java.util.StringTokenizer;

public class FDTD {

	int j, k, n;
	int ph_n, ph_m, ph_1D_N;

	private static final String RESULT_PATH = "result/";

	private static final double PI = 3.14159265358979;
	private static final double c = 2.99792458e8;
	private static final double ep0 = 8.85418782e-12;
	private static final double mu0 = 12.5663706144e-7;

	private double Len = 1.55e-6;
	private int L = 30;

	public double ds = Len / L;
	private double dt = ds / (Math.sqrt(2) * c);

	private double ysize_ = 7e-6;
	private double zsize_ = 17e-6;

	private double y1_rel = 0.14;
	private double y2_rel = 0.86;

	private double z1_rel = 0.059;
	private double z2_rel = 0.941;

	/*
	 * private double z1_rel = 0.08; private double z2_rel = 0.92;
	 */
	public double tf_y0_ = y1_rel * ysize_;
	public double tf_y_ = y2_rel * ysize_;

	public double tf_z0_ = z1_rel * zsize_;
	public double tf_z_ = z2_rel * zsize_;

	public int ysize = fl(ysize_);
	public int zsize = fl(zsize_);

	private int tf_y0 = fl(tf_y0_);
	private int tf_y = fl(tf_y_);

	private int tf_z0 = fl(tf_z0_);
	private int tf_z = fl(tf_z_);

	private int tsize = 500;

	private int aver = 1;

	private double[][] Ex = new double[ysize][zsize];
	private double[][] Exy = new double[ysize][zsize];
	private double[][] Exz = new double[ysize][zsize];

	private double[][] Hy = new double[ysize][zsize];
	private double[][] Hz = new double[ysize][zsize];

	public double[][] IntExAver = new double[ysize][zsize];

	public double[][] ep = new double[ysize][zsize];
	private double[][] mu = new double[ysize][zsize];

	private double N_rel = 0.64;
	private int pml_N = (int) (Math.floor(N_rel * Len / ds) - 2);
	private double sigm_max = 0.5e5;

	private double[] sig_y = new double[pml_N];
	private double[] sig_ys = new double[pml_N];
	private double[] sig_z = new double[pml_N];
	private double[] sig_zs = new double[pml_N];

	private int focus;
	private int real_focus;
	private double f_s = 0;

	private double E0 = 1;
	private double H0 = 1;

	private double Q = 0 * PI;

	private double f = c / Len;
	private double w = 2 * PI * f;
	private double K = 2 * PI / Len;
	private double Ky = K * Math.sin(Q);
	private double Kz = K * Math.cos(Q);

	double C1 = dt / ds;

	int aver_time = (int) Math.floor(aver * (1 / (dt * f)) + 0.5);

	int t2 = tsize - 1;
	int del_t = (int) Math.floor(1 / (4 * f * dt) + 0.5);
	int t1 = t2 - del_t;

	public static void main(String[] args) {
		new FDTD();
	}

	public FDTD() {

		init();
		initPML();
		setObject();

		long time1 = Calendar.getInstance().getTimeInMillis();

		run();

		long time2 = Calendar.getInstance().getTimeInMillis();
		System.out.println("Computation finished. Time: " + (time2 - time1)
				/ 1000 + " sec");

		exportAll();

	}

	public FDTD(int h) {

	}

	public void exportAll() {
		exportInten();
		exportPar();
		exportEps();
	}

	/**
	 * This method is used for setting object
	 */
	private void setObject() {
		/*
		double media = 2.83;
		double z0 = tf_z0_;

		setRect2(media, z0, 8e-6);
		setRect2(media, z0 + 8e-6, 4e-6, 0.25e-6);
		setGradientRI(media, 1, z0 + 12e-6, tf_z_, 0.25e-6);

		double r1 = 93e-9;
		double r2 = 108e-9;
		double r3 = 125e-9;

		double z1 = z0 + 5e-6;
		double z2 = z1 + 3e-6;
		for (int i = 0; i < 3; i++) {
			setPhRod2(1, ysize_ / 2 + i * 0.25e-6 - ds, z1, z2, 0.25e-6, r1);
			setPhRod2(1, ysize_ / 2 - i * 0.25e-6, z1, z2, 0.25e-6, r1);
		}

		for (int i = 3; i < 6; i++) {
			setPhRod2(1, ysize_ / 2 + i * 0.25e-6, z1, z2, 0.25e-6, r2);
			setPhRod2(1, ysize_ / 2 - i * 0.25e-6, z1, z2, 0.25e-6, r2);
		}

		for (int i = 6; i < 9; i++) {
			setPhRod2(1, ysize_ / 2 + i * 0.25e-6 - ds, z1, z2, 0.25e-6, r3);
			setPhRod2(1, ysize_ / 2 - i * 0.25e-6, z1, z2, 0.25e-6, r3);
		}

		// setMikaelyanLens(media, tf_z0_ + 1e-6, 3e-6);
		// setPhMikLens2(media, tf_z0_ + 1e-6, 3e-6, media, 3e-6, 0.25e-6);
		 * 
		 */
		
		setSph(2, ysize_ / 2, zsize_ / 2, 2.5e-6);

	}

	public void init() {
		for (j = 0; j < ysize; j++) {
			for (k = 0; k < zsize; k++) {
				Ex[j][k] = 0;
				Exy[j][k] = 0;
				Exz[j][k] = 0;
				Hy[j][k] = 0;
				Hz[j][k] = 0;

				IntExAver[j][k] = 0;

				ep[j][k] = ep0;
				mu[j][k] = mu0;

			}
		}
	}

	public void initPML() {
		System.out.println("PML: " + pml_N);
		System.out.println("PML: " + pml_N * ds);
		for (int i = 0; i < pml_N; i++) {

			sig_y[i] = getSig(sigm_max, i);
			sig_z[i] = getSig(sigm_max, i);

			sig_ys[i] = sig_y[i] * (mu0 / ep0);
			sig_zs[i] = sig_z[i] * (mu0 / ep0);
		}
	}

	double getSig(double sig_max, int n) {

		double s = ((double) n) / (pml_N - 1);
		double t = sig_max * Math.pow(s, 2);
		return t;

	}

	private double qu(double x) {
		return Math.pow(x, 2);
	}

	public void run() {

		for (n = 0; n < tsize; n++) {

			connectHy();
			nstep_Hy();

			nstep_Hz();
			connectHz();

			nstep_Ex();
			connectEx();

			pmlCond();

			if (n >= (tsize - aver_time) && n <= tsize) {
				for (j = 0; j < ysize; j++) {
					for (k = 0; k < zsize; k++) {
						double intensity = qu(Ex[j][k]) + qu(Hy[j][k])
								+ qu(Hz[j][k]);
						IntExAver[j][k] += (intensity / aver_time);
					}
				}
			}

			System.out.println("Computation in progress. Step " + n + "/"
					+ tsize);

		}

		System.out.println("Computation is succesfully finished");

		int f_s_ = (int) Math.floor(f_s / ds + 0.5);
		real_focus = findMax(f_s_, tf_z - 1);

	}

	void nstep_Hy() {

		for (k = pml_N; k < zsize - pml_N; k++) {
			for (j = pml_N; j < ysize - pml_N; j++)
				Hy[j][k] -= (C1 / mu[j][k]) * (Ex[j][k + 1] - Ex[j][k]);
		}

	}

	void nstep_Hz() {

		for (k = pml_N; k < zsize - pml_N; k++) {
			for (j = pml_N; j < ysize - pml_N; j++)
				Hz[j][k] += (C1 / mu[j][k]) * (Ex[j + 1][k] - Ex[j][k]);
		}
	}

	void nstep_Ex() {
		for (k = pml_N; k < zsize - pml_N; k++) {
			for (j = pml_N; j < ysize - pml_N; j++)
				Ex[j][k] += (C1 / ep[j][k])
						* ((Hz[j][k] - Hz[j - 1][k]) - (Hy[j][k] - Hy[j][k - 1]));
		}
	}

	void connectHy() {

		for (j = tf_y0; j <= tf_y; j++) {

			Hy[j][tf_z0 - 1] += (C1 / mu[j][tf_z0 - 1])
					* planeWaveEx(j, tf_z0, n); // left boundary
			Hy[j][tf_z] -= (C1 / mu[j][tf_z]) * planeWaveEx(j, tf_z, n); // right
			// boundary

		}
	}

	void connectEx() {

		for (j = tf_y0; j <= tf_y; j++) {

			Ex[j][tf_z0] += (C1 / ep[j][tf_z0]) * planeWaveHy(j, tf_z0 - 1, n); // left
			// boundary
			Ex[j][tf_z] -= (C1 / ep[j][tf_z]) * planeWaveHy(j, tf_z, n); // right
			// boundary

		}

		for (k = tf_z0; k <= tf_z; k++) {

			Ex[tf_y][k] += (C1 / ep[tf_y][k]) * planeWaveHz(tf_y, k, n); // top
			// boundary
			Ex[tf_y0][k] -= (C1 / ep[tf_y0][k])
					* planeWaveHz((tf_y0 - 1), k, n); // bottom boundary

		}
	}

	void connectHz() {

		for (k = tf_z0; k <= tf_z; k++) {

			Hz[tf_y][k] += (C1 / mu[tf_y][k]) * planeWaveEx(tf_y, k, n); // top
			// boundary
			Hz[tf_y0 - 1][k] -= (C1 / mu[tf_y0 - 1][k])
					* planeWaveEx(tf_y0, k, n); // bottom boundary

		}

	}

	double planeWaveEx(int j_, int k_, int n_) {

		double res = E0 * Math.sin(Ky * j_ * ds + Kz * k_ * ds - w * n_ * dt);
		return res;

	}

	double planeWaveHy(int j_, int k_, int n_) {

		return Math.cos(Q)
				* (H0 / (c * mu0))
				* Math.sin(Ky * j_ * ds + Kz * (k_ + 0.5) * ds - w * (n_ + 0.5)
						* dt);

	}

	double planeWaveHz(int j_, int k_, int n_) {

		return Math.sin(Q)
				* (H0 / (c * mu0))
				* Math.sin(Ky * (j_ + 0.5) * ds + Kz * k_ * ds - w * (n_ + 0.5)
						* dt);

	}

	void pmlCond() {

		pml_left();
		pml_lt();
		pml_top();
		pml_rt();
		pml_right();
		pml_rb();
		pml_bottom();
		pml_lb();

		pml_save();

	}

	void pml_left() {
		// Hy
		for (j = pml_N; j < ysize - pml_N; j++) {
			for (k = 0; k < pml_N - 1; k++)
				Hy[j][k] = ((1 - 0.5 * (dt * sig_zs[pml_N - k - 1] / mu[j][k]))
						* Hy[j][k] - (C1 / mu[j][k] * (Exy[j][k + 1]
						+ Exz[j][k + 1] - Exy[j][k] - Exz[j][k])))
						/ (1 + 0.5 * (dt * sig_zs[pml_N - k - 1] / mu[j][k]));

			k = pml_N - 1;

			Hy[j][k] = ((1 - 0.5 * (dt * sig_zs[0] / mu[j][k])) * Hy[j][k] - (C1
					/ mu[j][k] * (Ex[j][k + 1] - Exy[j][k] - Exz[j][k])))
					/ (1 + 0.5 * (dt * sig_zs[0] / mu[j][k]));
		}

		// Hz
		for (j = pml_N; j < ysize - pml_N; j++) {
			for (k = 0; k < pml_N; k++) {

				Hz[j][k] += (C1 / mu[j][k])
						* (Exy[j + 1][k] + Exz[j + 1][k] - Exy[j][k] - Exz[j][k]);

			}
		}

		// Exy
		for (j = pml_N; j < ysize - pml_N; j++) {
			for (k = 1; k < pml_N; k++) {

				Exy[j][k] += C1 / ep[j][k] * (Hz[j][k] - Hz[j - 1][k]);

			}
		}

		// Exz
		for (j = pml_N; j < ysize - pml_N; j++) {
			for (k = 1; k < pml_N; k++) {

				Exz[j][k] = ((1 - 0.5 * (dt * sig_z[pml_N - k - 1] / ep[j][k]))
						* Exz[j][k] - (C1 / ep[j][k] * (Hy[j][k] - Hy[j][k - 1])))
						/ (1 + 0.5 * (dt * sig_z[pml_N - k - 1] / ep[j][k]));

			}
		}
	}

	void pml_lt() {

		// Hy
		for (j = ysize - pml_N; j < ysize; j++) {
			for (k = 0; k < pml_N; k++)

				Hy[j][k] = ((1 - 0.5 * (dt * sig_zs[pml_N - k - 1] / mu[j][k]))
						* Hy[j][k] - (C1 / mu[j][k] * (Exy[j][k + 1]
						+ Exz[j][k + 1] - Exy[j][k] - Exz[j][k])))
						/ (1 + 0.5 * (dt * sig_zs[pml_N - k - 1] / mu[j][k]));

		}

		// Hz
		for (j = ysize - pml_N; j < ysize - 1; j++) {
			for (k = 0; k < pml_N; k++) {

				Hz[j][k] = ((1 - 0.5 * (dt * sig_ys[pml_N + j - ysize] / mu[j][k]))
						* Hz[j][k] + (C1 / mu[j][k] * (Exy[j + 1][k]
						+ Exz[j + 1][k] - Exy[j][k] - Exz[j][k])))
						/ (1 + 0.5 * (dt * sig_ys[pml_N + j - ysize] / mu[j][k]));

			}
		}

		// Exy
		for (j = ysize - pml_N; j < ysize; j++) {
			for (k = 1; k < pml_N; k++) {

				Exy[j][k] = ((1 - 0.5 * (dt * sig_y[pml_N + j - ysize] / ep[j][k]))
						* Exy[j][k] + (C1 / ep[j][k] * (Hz[j][k] - Hz[j - 1][k])))
						/ (1 + 0.5 * (dt * sig_y[pml_N + j - ysize] / ep[j][k]));

			}
		}

		// Exz
		for (j = ysize - pml_N; j < ysize; j++) {
			for (k = 1; k < pml_N; k++) {

				Exz[j][k] = ((1 - 0.5 * (dt * sig_z[pml_N - k - 1] / ep[j][k]))
						* Exz[j][k] - (C1 / ep[j][k] * (Hy[j][k] - Hy[j][k - 1])))
						/ (1 + 0.5 * (dt * sig_z[pml_N - k - 1] / ep[j][k]));

			}
		}

	}

	void pml_top() {

		// Hy

		for (j = ysize - pml_N; j < ysize; j++) {
			for (k = pml_N; k < zsize - pml_N; k++) {

				Hy[j][k] -= (C1 / mu[j][k])
						* (Exy[j][k + 1] + Exz[j][k + 1] - Exy[j][k] - Exz[j][k]);
			}
		}

		// Hz
		for (j = ysize - pml_N; j < ysize - 1; j++) {
			for (k = pml_N; k < zsize - pml_N; k++) {

				Hz[j][k] = ((1 - 0.5 * (dt * sig_ys[pml_N + j - ysize] / mu[j][k]))
						* Hz[j][k] + (C1 / mu[j][k] * (Exy[j + 1][k]
						+ Exz[j + 1][k] - Exy[j][k] - Exz[j][k])))
						/ (1 + 0.5 * (dt * sig_ys[pml_N + j - ysize] / mu[j][k]));

			}
		}

		// Exy

		for (j = ysize - pml_N; j < ysize; j++) {
			for (k = pml_N; k < zsize - pml_N; k++) {

				Exy[j][k] = ((1 - 0.5 * (dt * sig_y[pml_N + j - ysize] / ep[j][k]))
						* Exy[j][k] + (C1 / ep[j][k] * (Hz[j][k] - Hz[j - 1][k])))
						/ (1 + 0.5 * (dt * sig_y[pml_N + j - ysize] / ep[j][k]));

			}
		}

		// Exz
		for (j = ysize - pml_N; j < ysize; j++) {
			for (k = pml_N; k < zsize - pml_N; k++) {

				Exz[j][k] -= (C1 / ep[j][k]) * (Hy[j][k] - Hy[j][k - 1]);

			}
		}
	}

	void pml_rt() {

		// Hy

		for (j = ysize - pml_N; j < ysize; j++) {
			for (k = zsize - pml_N; k < zsize - 1; k++)

				Hy[j][k] = ((1 - 0.5 * (dt * sig_zs[pml_N + k - zsize] / mu[j][k]))
						* Hy[j][k] - (C1 / mu[j][k] * (Exy[j][k + 1]
						+ Exz[j][k + 1] - Exy[j][k] - Exz[j][k])))
						/ (1 + 0.5 * (dt * sig_zs[pml_N + k - zsize] / mu[j][k]));

		}

		// Hz
		for (j = ysize - pml_N; j < ysize - 1; j++) {
			for (k = zsize - pml_N; k < zsize; k++) {

				Hz[j][k] = ((1 - 0.5 * (dt * sig_ys[pml_N + j - ysize] / mu[j][k]))
						* Hz[j][k] + (C1 / mu[j][k] * (Exy[j + 1][k]
						+ Exz[j + 1][k] - Exy[j][k] - Exz[j][k])))
						/ (1 + 0.5 * (dt * sig_ys[pml_N + j - ysize] / mu[j][k]));

			}
		}

		// Exy

		for (j = ysize - pml_N; j < ysize; j++) {
			for (k = zsize - pml_N; k < zsize; k++) {

				Exy[j][k] = ((1 - 0.5 * (dt * sig_y[pml_N + j - ysize] / ep[j][k]))
						* Exy[j][k] + (C1 / ep[j][k] * (Hz[j][k] - Hz[j - 1][k])))
						/ (1 + 0.5 * (dt * sig_y[pml_N + j - ysize] / ep[j][k]));

			}
		}

		// Exz
		for (j = ysize - pml_N; j < ysize; j++) {
			for (k = zsize - pml_N; k < zsize; k++) {

				Exz[j][k] = ((1 - 0.5 * (dt * sig_z[pml_N + k - zsize] / ep[j][k]))
						* Exz[j][k] - (C1 / ep[j][k] * (Hy[j][k] - Hy[j][k - 1])))
						/ (1 + 0.5 * (dt * sig_z[pml_N + k - zsize] / ep[j][k]));

			}
		}
	}

	void pml_right() {

		// Hy
		for (j = pml_N; j < ysize - pml_N; j++)
			for (k = zsize - pml_N; k < zsize - 1; k++) {

				Hy[j][k] = ((1 - 0.5 * (dt * sig_zs[pml_N + k - zsize] / mu[j][k]))
						* Hy[j][k] - (C1 / mu[j][k] * (Exy[j][k + 1]
						+ Exz[j][k + 1] - Exy[j][k] - Exz[j][k])))
						/ (1 + 0.5 * (dt * sig_zs[pml_N + k - zsize] / mu[j][k]));

			}

		// Hz
		for (j = pml_N; j < ysize - pml_N; j++) {
			for (k = zsize - pml_N; k < zsize; k++) {

				Hz[j][k] += (C1 / mu[j][k])
						* (Exy[j + 1][k] + Exz[j + 1][k] - Exy[j][k] - Exz[j][k]);

			}
		}

		// Exy
		for (j = pml_N; j < ysize - pml_N; j++) {
			for (k = zsize - pml_N; k < zsize; k++) {

				Exy[j][k] += (C1 / ep[j][k]) * (Hz[j][k] - Hz[j - 1][k]);

			}
		}

		// Exz
		for (j = pml_N; j < ysize - pml_N; j++) {
			for (k = zsize - pml_N; k < zsize; k++) {

				Exz[j][k] = ((1 - 0.5 * (dt * sig_z[pml_N + k - zsize] / ep[j][k]))
						* Exz[j][k] - (C1 / ep[j][k] * (Hy[j][k] - Hy[j][k - 1])))
						/ (1 + 0.5 * (dt * sig_z[pml_N + k - zsize] / ep[j][k]));

			}
		}
	}

	void pml_rb() {

		// Hy

		for (j = 0; j < pml_N; j++) {
			for (k = zsize - pml_N; k < zsize - 1; k++)

				Hy[j][k] = ((1 - 0.5 * (dt * sig_zs[pml_N + k - zsize] / mu[j][k]))
						* Hy[j][k] - (C1 / mu[j][k] * (Exy[j][k + 1]
						+ Exz[j][k + 1] - Exy[j][k] - Exz[j][k])))
						/ (1 + 0.5 * (dt * sig_zs[pml_N + k - zsize] / mu[j][k]));

		}

		// Hz
		for (j = 0; j < pml_N; j++) {
			for (k = zsize - pml_N; k < zsize; k++) {

				Hz[j][k] = ((1 - 0.5 * (dt * sig_ys[pml_N - j - 1] / mu[j][k]))
						* Hz[j][k] + (C1 / mu[j][k] * (Exy[j + 1][k]
						+ Exz[j + 1][k] - Exy[j][k] - Exz[j][k])))
						/ (1 + 0.5 * (dt * sig_ys[pml_N - j - 1] / mu[j][k]));

			}
		}

		// Exy

		for (j = 1; j < pml_N; j++) {
			for (k = zsize - pml_N; k < zsize; k++) {

				Exy[j][k] = ((1 - 0.5 * (dt * sig_y[pml_N - j - 1] / ep[j][k]))
						* Exy[j][k] + (C1 / ep[j][k] * (Hz[j][k] - Hz[j - 1][k])))
						/ (1 + 0.5 * (dt * sig_y[pml_N - j - 1] / ep[j][k]));

			}
		}

		// Exz
		for (j = 1; j < pml_N; j++) {
			for (k = zsize - pml_N; k < zsize; k++) {

				Exz[j][k] = ((1 - 0.5 * (dt * sig_z[pml_N + k - zsize] / ep[j][k]))
						* Exz[j][k] - (C1 / ep[j][k] * (Hy[j][k] - Hy[j][k - 1])))
						/ (1 + 0.5 * (dt * sig_z[pml_N + k - zsize] / ep[j][k]));

			}
		}
	}

	void pml_bottom() {

		// Hy
		for (j = 0; j < pml_N; j++) {
			for (k = pml_N; k < zsize - pml_N; k++)

				Hy[j][k] -= (C1 / mu[j][k])
						* (Exy[j][k + 1] + Exz[j][k + 1] - Exy[j][k] - Exz[j][k]);

		}

		// Hz
		for (j = 0; j < pml_N - 1; j++) {
			for (k = pml_N; k < zsize - pml_N; k++) {

				Hz[j][k] = ((1 - 0.5 * (dt * sig_ys[pml_N - j - 1] / mu[j][k]))
						* Hz[j][k] + (C1 / mu[j][k] * (Exy[j + 1][k]
						+ Exz[j + 1][k] - Exy[j][k] - Exz[j][k])))
						/ (1 + 0.5 * (dt * sig_ys[pml_N - j - 1] / mu[j][k]));

			}

		}

		j = pml_N - 1;
		for (k = pml_N; k < zsize - pml_N; k++) {

			Hz[j][k] = ((1 - 0.5 * (dt * sig_ys[pml_N - j - 1] / mu[j][k]))
					* Hz[j][k] + (C1 / mu[j][k] * (Ex[j + 1][k] - Exy[j][k] - Exz[j][k])))
					/ (1 + 0.5 * (dt * sig_ys[pml_N - j - 1] / mu[j][k]));

		}

		// Exy
		for (j = 1; j < pml_N; j++) {
			for (k = pml_N; k < zsize - pml_N; k++) {

				Exy[j][k] = ((1 - 0.5 * (dt * sig_y[pml_N - j - 1] / ep[j][k]))
						* Exy[j][k] + (C1 / ep[j][k] * (Hz[j][k] - Hz[j - 1][k])))
						/ (1 + 0.5 * (dt * sig_y[pml_N - j - 1] / ep[j][k]));

			}
		}

		// Exz
		for (j = 1; j < pml_N; j++) {
			for (k = pml_N; k < zsize - pml_N; k++) {

				Exz[j][k] -= (C1 / ep[j][k]) * (Hy[j][k] - Hy[j][k - 1]);

			}
		}
	}

	void pml_lb() {

		// Hy

		for (j = 0; j < pml_N; j++) {
			for (k = 0; k < pml_N; k++)

				Hy[j][k] = ((1 - 0.5 * (dt * sig_zs[pml_N - k - 1] / mu[j][k]))
						* Hy[j][k] - (C1 / mu[j][k] * (Exy[j][k + 1]
						+ Exz[j][k + 1] - Exy[j][k] - Exz[j][k])))
						/ (1 + 0.5 * (dt * sig_zs[pml_N - k - 1] / mu[j][k]));

		}

		// Hz
		for (j = 0; j < pml_N; j++) {
			for (k = 0; k < pml_N; k++) {

				Hz[j][k] = ((1 - 0.5 * (dt * sig_ys[pml_N - j - 1] / mu[j][k]))
						* Hz[j][k] + (C1 / mu[j][k] * (Exy[j + 1][k]
						+ Exz[j + 1][k] - Exy[j][k] - Exz[j][k])))
						/ (1 + 0.5 * (dt * sig_ys[pml_N - j - 1] / mu[j][k]));

			}
		}

		// Exy

		for (j = 1; j < pml_N; j++) {
			for (k = 1; k < pml_N; k++) {

				Exy[j][k] = ((1 - 0.5 * (dt * sig_y[pml_N - j - 1] / ep[j][k]))
						* Exy[j][k] + (C1 / ep[j][k] * (Hz[j][k] - Hz[j - 1][k])))
						/ (1 + 0.5 * (dt * sig_y[pml_N - j - 1] / ep[j][k]));

			}
		}

		// Exz
		for (j = 1; j < pml_N; j++) {
			for (k = 1; k < pml_N; k++) {

				Exz[j][k] = ((1 - 0.5 * (dt * sig_z[pml_N - k - 1] / ep[j][k]))
						* Exz[j][k] - (C1 / ep[j][k] * (Hy[j][k] - Hy[j][k - 1])))
						/ (1 + 0.5 * (dt * sig_z[pml_N - k - 1] / ep[j][k]));

			}
		}
	}

	private void pml_save() {

		for (j = 0; j < ysize; j++) {
			for (k = 0; k < pml_N; k++)
				Ex[j][k] = Exy[j][k] + Exz[j][k];

			for (k = zsize - pml_N; k < zsize; k++)
				Ex[j][k] = Exy[j][k] + Exz[j][k];
		}

		for (k = pml_N; k < zsize - pml_N; k++) {
			for (j = 0; j < pml_N; j++)
				Ex[j][k] = Exy[j][k] + Exz[j][k];

			for (j = ysize - pml_N; j < ysize; j++)
				Ex[j][k] = Exy[j][k] + Exz[j][k];
		}

	}

	private void exportInten() {

		System.out.println("Exporting");
		PrintWriter pw;
		try {
			pw = new PrintWriter(RESULT_PATH + "intenEx.dat");

			for (j = 0; j < ysize; j++) {
				pw.print(IntExAver[j][0]);
				for (k = 1; k < zsize; k++) {
					pw.print(", " + IntExAver[j][k]);
				}
				pw.println();
			}

			pw.flush();
			pw.close();
		} catch (FileNotFoundException e) {

		}
	}

	private int findMax(int k1, int k2) {

		double res = -1;
		int num = 0;
		for (int k = k1; k <= k2; k++)
			if (IntExAver[ysize / 2][k] >= res) {
				res = IntExAver[ysize / 2][k];
				num = k;
			}
		return num;
	}

	public double maxValue(int k1, int k2) {

		double res = -1;
		for (int k = k1; k <= k2; k++)
			if (IntExAver[ysize / 2][k] >= res) {
				res = IntExAver[ysize / 2][k];
			}
		return res;
	}

	private void exportPar() {

		PrintWriter pw;
		try {
			pw = new PrintWriter(RESULT_PATH + "pars.txt");
			pw.println(ysize + ", " + zsize + ", " + tf_y0 + ", " + tf_y + ", "
					+ tf_z0 + ", " + tf_z + ", " + ds + ", " + dt + ", " + L
					+ ", " + real_focus + ", " + ph_n + ", " + ph_m + ", "
					+ ph_1D_N);

			pw.flush();
			pw.close();
		} catch (FileNotFoundException e) {
		}
	}

	private void exportEps() {
		try {
			PrintWriter pw = new PrintWriter(RESULT_PATH + "refr_index.dat");
			for (j = 0; j < ysize; j++) {
				pw.print(ep[j][0]);
				for (k = 1; k < zsize; k++) {
					pw.print(", " + ep[j][k]);
				}
				pw.println();
			}
			pw.flush();
			pw.close();
		} catch (FileNotFoundException e) {
		}

	}

	void setRect2(double ref_ind, double z0, double dz) {
		setRect(ref_ind * ref_ind, tf_y0, fl(z0), tf_y - tf_y0 + 1, fl(dz));
	}

	void setRect2(double ref_ind, double z0, double dz, double r) {

		int k0 = (int) Math.floor(0.5 + z0 / ds);
		int dk = (int) Math.floor(0.5 + dz / ds);

		setRect(ref_ind * ref_ind, fl(ysize_ / 2 - r), k0, fl(2 * r), dk);
	}

	void setRect(double eps, int j0, int k0, int dj, int dk) {

		for (j = j0; j < j0 + dj; j++) {
			for (k = k0; k < k0 + dk; k++) {
				ep[j][k] = eps * ep0;
			}
		}
	}

	void setPhRod(double n, int j0, int k0, int k, int D, int R) {

		double eps = n * n;

		ph_n = (int) Math.floor((k - k0) / D + 0.5);
		int del = (int) Math.floor(((k - k0) - ph_n * D) / 2 + 0.5);
		int d1 = (int) Math.floor(D / 2);

		for (int u = 0; u < ph_n; u++) {
			setSph(eps, j0, k0 + del + d1 + u * D, R);
		}

	}

	void setPhRod2(double n, double y0, double z0, double z, double D, double R) {

		int j0 = fl(y0);
		int k0 = fl(z0);
		int k = fl(z);
		int D_ = fl(D);
		int R_ = fl(R);

		setPhRod(n, j0, k0, k, D_, R_);

	}

	void setSph(double eps, int j0, int k0, int R) {
		for (j = j0 - R; j < j0 + R + 1; j++) {
			for (k = k0 - R; k < k0 + R + 1; k++) {
				if (Math.pow((j - j0), 2) + Math.pow((k - k0), 2) < Math.pow(R,
						2)) {
					ep[j][k] = eps * ep0;
				}
			}
		}
	}

	void setSph(double eps, double y0, double z0, double R) {
		setSph(eps, fl(y0), fl(z0), fl(R));
	}

	// convert double to the int on the net
	public int fl(double d) {
		return (int) Math.floor(d / ds + 0.5);
	}

	public void setPhotonicMikLens(PhotonicMikaelyanLens lens) {

		double n = lens.getN_pc();
		double z0 = lens.getZ0();
		double L = lens.getL();

		this.setRect2(n, z0, L);

		for (int i = 0; i < lens.getN(); i++) {
			double y0 = tf_y0_ + lens.getY0(i);
			double d = lens.getLConstant(i);
			double r = lens.getRadii(i);
			setPhRod2(1, y0, z0, z0 + L, d, r);

		}
	}

	void setPhMikLens(double ref_ind, int k0, int a, double n0, int Lmik, int D) {

		int M = (int) Math.floor((tf_y - tf_y0) / (2 * D) + 0.5);
		int N = (int) Math.floor(a / D + 0.5);

		int del = (int) Math.floor(((tf_y - tf_y0) - M * D) / 2 + 0.5);

		int y_mid = (int) Math.floor(ysize / 2 + 0.5);

		for (int v = 0; v < M; v++) {

			double r_ = getPhMikRadii(ref_ind, n0, a * ds, Lmik * ds, D * ds,
					N, v);

			int r = fl(r_);
			System.out.println("real: " + r_ * 1e6 + " | getting: " + r * ds
					* 1e6);

			setPhRod(1, y_mid + v * D - 1, k0, k0 + a, D, r);
			setPhRod(1, y_mid - v * D, k0, k0 + a, D, r);
		}

	}

	private void setPhMikLens2(double ref_ind, double z0, double a, double n0,
			double Lmik, double D) {
		f_s = z0 + a;
		this.setRect2(ref_ind, z0, a);
		setPhMikLens(ref_ind, fl(z0), fl(a), n0, fl(Lmik), fl(D));

	}

	double getPhMikRadii(double ref_ind, double n0, double a, double Lmik,
			double D, int N, int l) {

		double r1 = 1 / Math.cosh(PI * D * l / (2 * Lmik));
		double r2 = D / (2 * (ref_ind - 1));
		double r3 = n0 * r1 * Lmik / a;
		double res = r2 * (ref_ind - r3);

		return res;

	}

	/**
	 * PC-lens with constant lattice period.
	 * 
	 * @param hole
	 *            refractive index of the holes
	 * @param media
	 *            refractive index of the holes
	 * @param z0
	 *            z of first hole (bottom left)
	 * @param y0
	 *            y of first hole (bottom left)
	 * @param z_n
	 *            number of holes along the z axis
	 * @param y_n
	 *            number of holes along the y axis
	 * @param r
	 *            radii oh the holes
	 * @param d
	 *            lattice constant
	 */
	private void setPhLensWithCP(double hole, double media, double z0,
			double y0, int z_n, int y_n, double r, double d) {
		setRect2(media, z0 - d / 2, z_n * d);

		for (int i = 0; i < z_n; i++) {
			for (int j = 0; j < y_n; j++) {
				setSph(qu(hole), y0 + d * j, z0 + d * i, r);
			}
		}
	}

	/**
	 * PC-lens with constant lattice period.
	 * 
	 * @param hole
	 *            refractive index of the holes
	 * @param media
	 *            refractive index of the holes
	 * @param z0
	 *            z of first hole (bottom left)
	 * @param y0
	 *            y of first hole (bottom left)
	 * @param z_n
	 *            number of holes along the z axis
	 * @param y_n
	 *            number of holes along the y axis
	 * @param r
	 *            radii oh the holes
	 * @param d
	 *            lattice constant
	 */
	private void setTringularPhLensWithCP(double hole, double media, double z0,
			double y0, int z_n, int y_n, double r, double d) {
		setRect2(media, z0 - d / 2, z_n * d);

		for (int i = 0; i < z_n; i++) {
			if (i % 2 == 0) {
				for (int j = 0; j < y_n; j++) {
					setSph(qu(hole), y0 + d * j, z0 + d * i, r);
				}
			} else {
				for (int j = 0; j < y_n - 1; j++) {
					setSph(qu(hole), y0 + d * (j + 0.5), z0 + d * i, r);
				}
			}
		}
	}

	public double maxEps() {
		double max = 0;

		for (int j = 0; j < ysize; j++) {
			for (int k = 0; k < zsize; k++) {
				if (ep[j][k] > max) {
					max = ep[j][k];
				}
			}
		}

		return max;
	}

	public double maxInt() {
		double max = 0;
		for (int j = 0; j < ysize; j++) {
			for (int k = 0; k < zsize; k++) {
				if (IntExAver[j][k] > max) {
					max = ep[j][k];
				}
			}
		}
		return max;
	}

	/**
	 * Sets Mikaelyan lens with height of tfsf-height
	 * 
	 * @param n0
	 *            refractive index in center
	 * @param z0
	 *            the left boundary
	 * @param L
	 *            length of the lens
	 */
	public void setMikaelyanLens(double n0, double z0, double L) {
		int k0 = fl(z0);
		int dk = fl(L);
		for (int k = k0; k < k0 + dk; k++) {
			for (int j = tf_y0; j < tf_y; j++) {
				double n_mik = nmik(n0, (ysize_ / 2 - j * ds), L);
				ep[j][k] *= (n_mik * n_mik);
			}
		}
	}

	public double nmik(double n0, double x, double L) {
		return (n0 / Math.cosh(PI * Math.abs(x) / 2 / L));
	}

	public void loadEpslon(File file) {
		try {
			BufferedReader br = new BufferedReader(new FileReader(file));
			for (int j = 0; j < ysize; j++) {
				String str = br.readLine();
				StringTokenizer st = new StringTokenizer(str, ",");
				for (int k = 0; k < zsize; k++) {
					ep[j][k] = Double.valueOf(st.nextToken()).doubleValue();
				}
			}
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	private void setGradientRI(double n1, double n2, double z1, double z2,
			double r) {
		int j_mid = fl(ysize_ / 2);
		int r_ = fl(r);
		int k1 = fl(z1);
		int k2 = fl(z2);

		double a = (n2 - n1) / (z2 - z1);
		double b = (n1 * z2 - z1 * n2) / (z2 - z1);

		for (int k = k1; k < k2; k++) {
			double newn = a * k * ds + b;
			for (int j = j_mid - r_; j < j_mid + r_; j++) {
				ep[j][k] = ep0 * Math.pow(newn, 2);
			}
		}
	}
	
}
