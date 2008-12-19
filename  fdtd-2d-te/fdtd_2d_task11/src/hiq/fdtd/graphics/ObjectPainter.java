package hiq.fdtd.graphics;

import hiq.fdtd.fdtd.FDTD;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;

public class ObjectPainter extends AbstractPainter {

	public ObjectPainter(FDTD fdtd) {
		super(fdtd);
	}

	@Override
	public void paint(Graphics g) {
		super.paint(g);
		Graphics2D g2 = (Graphics2D) g;
		
		double k1 = (double)fdtd.zsize / (double)xsize;
		double k2 = (double)fdtd.ysize / (double)ysize;
		
		double eps_max = fdtd.maxEps();
		
		for (int i = 0; i < xsize; i++){
			int u = (int) Math.floor(k1 * i);
			for (int j = 0; j < ysize; j++){
				int v = (int) Math.floor(k2 * j);
				
				double eps = fdtd.ep[v][u];
				int p = 256 - (int) Math.floor(eps / eps_max * 256);
				g2.setColor(new Color(p, p, p));
				g2.drawRect(i, j, 1, 1);
			}
		}
	}
}
