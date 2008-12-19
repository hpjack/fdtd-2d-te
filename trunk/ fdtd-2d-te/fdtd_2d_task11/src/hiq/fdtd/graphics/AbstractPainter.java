package hiq.fdtd.graphics;

import hiq.fdtd.fdtd.FDTD;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.event.ComponentEvent;
import java.awt.event.ComponentListener;

import javax.swing.JFrame;

public abstract class AbstractPainter extends JFrame {
	protected FDTD fdtd;

	protected int xsize = 800;
	protected int ysize = 500;

	public AbstractPainter(FDTD fdtd) {
		this.fdtd = fdtd;
		paintGUI();
	}

	private void paintGUI() {
		this.setSize(xsize, ysize);
		this.setTitle("object");
		this.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		this.setVisible(true);

		this.addComponentListener(new ComponentListener() {
			public void componentResized(ComponentEvent evt) {
				Component c = (Component) evt.getSource();

				// Get new size
				Dimension newSize = c.getSize();
				xsize = newSize.width;
				ysize = newSize.height;

				repaint();
			}

			@Override
			public void componentHidden(ComponentEvent arg0) {
				// TODO Auto-generated method stub

			}

			@Override
			public void componentMoved(ComponentEvent arg0) {
				// TODO Auto-generated method stub

			}

			@Override
			public void componentShown(ComponentEvent arg0) {
				// TODO Auto-generated method stub

			}
		});
	}
}
