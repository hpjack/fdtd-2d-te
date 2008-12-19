package hiq.fdtd.genetic.util;

public class Helper {
	public static int[] binaryView(int a) {
		int n = 0;
		for (int i = 0; i < 1024; i++) {
			if (Math.pow(2, i) > a) {
				n = i;
				break;
			}
		}

		int[] res = new int[n];
		int b = a;
		for (int i = n - 1; i >= 0; i--) {
			int c = b - ((int) Math.pow(2, i));
			if (c >= 0) {
				res[n - i - 1] = 1;
				b = c;
			} else {
				res[n - i - 1] = 0;
			}
		}
		return res;

	}

	public static int decView(int[] a) {
		int res = 0;
		int n = a.length;
		for (int i = 0; i < n; i++) {
			res += ((int) a[i] * Math.pow(2, n - i - 1));
		}
		return res;
	}

	public static int[] binaryView(int a, int max) {
		int n = max;

		int[] res = new int[n];
		int b = a;
		for (int i = n - 1; i >= 0; i--) {
			int c = b - ((int) Math.pow(2, i));
			if (c >= 0) {
				res[n - i - 1] = 1;
				b = c;
			} else {
				res[n - i - 1] = 0;
			}
		}
		return res;
	}

	public static int[] decimalToGray(int a, int max) {
		int[] binary = Helper.binaryView(a, max);
		int size = binary.length;
		int[] gray = new int[size];
		for (int i = size - 1; i > 0; i--) {
			gray[i] = (binary[i] + binary[i - 1]) % 2;
		}
		gray[0] = binary[0];
		return gray;
	}

	public static int grayToDecimal(int [] gray){
		int size = gray.length;
		int [] binary = new int [size];
		for (int i = size - 1; i > 0; i--){
			for (int j = i; j >= 0; j--){
				binary[i] += gray[j];
			}
			binary[i] %= 2;
		}
		binary[0] = gray[0];
		return Helper.decView(binary);
	}
}
