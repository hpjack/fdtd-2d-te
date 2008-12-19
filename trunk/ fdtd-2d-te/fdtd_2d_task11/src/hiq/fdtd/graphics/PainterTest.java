package hiq.fdtd.graphics;

import hiq.fdtd.fdtd.FDTD;

public class PainterTest {
	public static void main(String[] args) {
		FDTD fdtd = new FDTD();
		ObjectPainter op = new ObjectPainter(fdtd);
	}
}
