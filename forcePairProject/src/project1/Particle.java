package project1;

public class Particle {
	private String type;
	private double mass;
	private double epsilon;
	private double sigma;
	private double gamma;
	private double[] coord;
	
	public Particle(String ty, double m, double e, double s, double g, double[] cd) {
		type = ty;
		mass = m;
		epsilon = e;
		sigma = s;
		gamma = g;
		coord = cd;
	}
	
	public String getType() {
		return type;
	}
	public double getMass() {
		return mass;
	}
	public double getGamma() {
		return gamma;
	}
	public double getEpsilon() {
		return epsilon;
	}
	public double[] getCoord() {
		return coord;
	}
	public double getCoordAtPos(int pos) {
		return coord[pos];
	}
	public void setCoord(double[] c) {
		coord = c;
	}
	public void setCoordAtPos(double c, int pos) {
		coord[pos] = c;
	}
}
