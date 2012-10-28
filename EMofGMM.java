
import java.net.*;
import java.io.*;
import java.nio.*;
import java.util.regex.*;
import java.util.*;
import java.awt.*;
import java.awt.event.*;
import java.awt.geom.*;
import java.text.*;
import javax.swing.*;
import javax.swing.tree.*;
import javax.swing.border.*;
import javax.swing.table.*;
import javax.swing.event.*;
import java.nio.channels.*;
import java.awt.image.*;
import javax.imageio.*;
import Jama.*;

class MyUtil {
	public static void drawPoint(Graphics g, Matrix point, Color c, int size) {
		g.setColor(c);
		g.fillOval((int)(point.get(0, 0) - (size / 2)), (int)(point.get(1, 0) - (size / 2)), size, size);
	}
	public static void drawPoints(Graphics g, Vector points, Color c, int size) {
		for (int i = 0; i < points.size(); i++) {
			drawPoint(g, (Matrix)points.get(i), c, size);
		}
	}
	public static void drawGaussian(Graphics g, Matrix center, Matrix covar, Color c, int size) {
		Graphics2D g2 = (Graphics2D)g;
		AffineTransform old = g2.getTransform();
		EigenvalueDecomposition e = covar.eig();
		Matrix diag = e.getD();
		diag.set(0, 0, Math.sqrt(diag.get(0, 0)));
		diag.set(1, 1, Math.sqrt(diag.get(1, 1)));
		diag.timesEquals(2.0);
		
		Matrix t = e.getV().times(diag.times(new Matrix(new double[] {.01, 0, 0, .01}, 2)));
		g2.setTransform(new AffineTransform(
			t.get(0, 0), t.get(1, 0), t.get(0, 1), t.get(1, 1), center.get(0, 0), center.get(1, 0)));
		
		g2.setColor(c);
		g2.drawOval(-100, -100, 200, 200);
		
		g2.setTransform(old);
		drawPoint(g, center, c, size);
	}
	public static Matrix makePoint(double x, double y) {
		return new Matrix(new double[] {x, y}, 2);
	}
	
	public static Matrix findCentroid(Vector points) {
		Matrix mean = makePoint(0, 0);
		for (int i = 0; i < points.size(); i++) {
			Matrix p = (Matrix)points.get(i);
			mean.plusEquals(p);
		}
		return mean.times(1.0 / points.size());
	}
	public static Matrix findCentroid(Vector points, Matrix weights) {
		Matrix mean = makePoint(0, 0);
		double sum = 0.0;
		for (int i = 0; i < points.size(); i++) {
			Matrix p = (Matrix)points.get(i);
			double weight = weights.get(i, 0);
			
			mean.plusEquals(p.times(weight));
			sum += weight;
		}
		return mean.times(1.0 / sum);
	}
	
	public static Matrix findCovar(Vector points, Matrix mean) {
		Matrix covar = new Matrix(2, 2, 0.0);
		for (int i = 0; i < points.size(); i++) {
			Matrix p = ((Matrix)points.get(i)).minus(mean);
			covar.plusEquals(p.times(p.transpose()));
		}
		return covar.times(1.0 / points.size());
	}
	public static Matrix findCovar(Vector points, Matrix weights, Matrix mean) {
		Matrix covar = new Matrix(2, 2, 0.0);
		double sum = 0.0;
		for (int i = 0; i < points.size(); i++) {
			Matrix p = ((Matrix)points.get(i)).minus(mean);
			double weight = weights.get(i, 0);
			covar.plusEquals(p.times(p.transpose()).times(weight));
			sum += weight;
		}
		return covar.times(1.0 / sum);
	}
	
	public static void normalizeRows(Matrix m) {
		for (int r = 0; r < m.getRowDimension(); r++) {
			double sum = 0.0;
			for (int c = 0; c < m.getColumnDimension(); c++) {
				sum += m.get(r, c);
			}
			for (int c = 0; c < m.getColumnDimension(); c++) {
				m.set(r, c, m.get(r, c) / sum);
			}
		}
	}
	
	public static boolean isValid(Matrix m) {
		for (int r = 0; r < m.getRowDimension(); r++) {
			for (int c = 0; c < m.getColumnDimension(); c++) {
				Double d = new Double(m.get(r, c));
				if (d.isInfinite() || d.isNaN()) {
					return false;
				}
			}
		}
		return true;
	}
}

class Gaussian {
	public Matrix mean;
	public Matrix covar;
		
	public Gaussian(Vector points) {
		mean = MyUtil.findCentroid(points);
		covar = MyUtil.findCovar(points, mean);
	}
	public Gaussian(Vector points, Matrix weights) {
		mean = MyUtil.findCentroid(points, weights);
		covar = MyUtil.findCovar(points, weights, mean);
	}
	public Gaussian(Matrix mean, double sd) {
		this.mean = mean;
		covar = new Matrix(new double[] {sd * sd, 0, 0, sd * sd}, 2);
	}
	
	public void draw(Graphics g, Color c) {
		MyUtil.drawGaussian(g, mean, covar, c, 8);
	}
	public double p(Matrix x) {
		x = x.minus(mean);
		return (1.0 / Math.sqrt(Math.pow(Math.PI * 2.0, 2.0) * covar.det())) *
			Math.exp(-0.5 * x.transpose().times(covar.inverse().times(x)).get(0, 0));
	}
	public boolean isValid() {
		try {
			covar.inverse();
		} catch (Exception e) {
			return false;
		}
		return MyUtil.isValid(mean) && MyUtil.isValid(covar);
	}
}

class EMofGMM extends JFrame {
	
	Vector points = new Vector();
	Vector gaussians = new Vector();
	
	public static void main(String[] args) throws Exception {
		new EMofGMM();
	}
	public EMofGMM() {
		setTitle("EM of GMM");
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		setSize(640, 480);
		
		getContentPane().setLayout(new BorderLayout());
		
		// make it so you can add points
		JPanel c = new JPanel() {
			public void paint(Graphics g) {
				super.paint(g);
				MyUtil.drawPoints(g, points, Color.BLACK, 4);
				for (int i = 0; i < gaussians.size(); i++)
					((Gaussian)gaussians.get(i)).draw(g, Color.RED);
			}
		};
		c.addMouseListener(new MouseAdapter() {
			public void mousePressed(MouseEvent e) {
				Matrix point = MyUtil.makePoint(e.getX(), e.getY());
				if (e.getButton() == MouseEvent.BUTTON1) {
					points.add(point);
				} else {
					gaussians.add(new Gaussian(point, 25));
				}
				repaint();
			}
		});
		c.setBackground(Color.WHITE);
		getContentPane().add(c);
		
		// make it so you can step the algorithm
		JPanel p = new JPanel();
		JButton b;
		p.add(b = new JButton("Step"));
		b.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				step();
				repaint();
			}
		});
		p.add(b = new JButton("Clear points"));
		b.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				clearPoints();
				repaint();
			}
		});
		p.add(b = new JButton("Clear means"));
		b.addActionListener(new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				clearMeans();
				repaint();
			}
		});
		getContentPane().add(p, BorderLayout.SOUTH);
		
		show();
	}
	
	public void clearPoints() {
		points.clear();
	}
	public void clearMeans() {
		gaussians.clear();
	}	
	public void step() {
		// find the extent to which each point is associated with each gaussian
		Matrix weights = new Matrix(points.size(), gaussians.size());
		for (int pi = 0; pi < points.size(); pi++) {
			Matrix p = (Matrix)points.get(pi);			
			for (int gi = 0; gi < gaussians.size(); gi++) {
				Gaussian g = (Gaussian)gaussians.get(gi);
				weights.set(pi, gi, g.p(p));
			}
		}
		MyUtil.normalizeRows(weights);
		
		// create a new set of guassians based on the point-weights
		Vector newGaussians = new Vector();
		for (int gi = 0; gi < gaussians.size(); gi++) {
			Matrix giWeights = weights.getMatrix(0, points.size() - 1, gi, gi);
			Gaussian g = new Gaussian(points, giWeights);
			if (g.isValid())
				newGaussians.add(g);
		}
		gaussians = newGaussians;
	}
}
