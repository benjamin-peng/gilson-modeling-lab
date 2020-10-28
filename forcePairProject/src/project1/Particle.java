package project1;

public class Particle {
	private String type;
	private double mass;
	private double epsilon;
	private double sigma;
	private double gamma;
	private double[] coord;
	private double[] momentum;
	private double[] momentumPrime;
	private double[] force;
	
	public Particle(String ty, double m, double e, double s, double g, double[] cd) {
		type = ty;
		mass = m;
		epsilon = e;
		sigma = s;
		gamma = g;
		coord = cd;
		momentum = new double[]{0, 0, 0};
		momentumPrime = new double[]{0, 0, 0};
		force = new double[]{0, 0, 0};
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
	public double getSigma() {
		return sigma;
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
	public double getSpeed() {
		double speed = Math.sqrt(Math.pow(momentum[0]/mass, 2) + Math.pow(momentum[1]/mass, 2) + Math.pow(momentum[2]/mass, 2));
		return speed;
	}
	public void setCoord(double[] c) {
		coord = c;
	}
	public void setCoordAtPos(double c, int pos) {
		coord[pos] = c;
	}
	
	public double[] getMomentum(boolean prime) {
		if (prime) return momentumPrime;
		else return momentum;
	}
	public double getMomentumAtPos(boolean prime, int pos) {
		if (prime) return momentumPrime[pos];
		else return momentum[pos];
	}
	public void setMomentum(boolean prime, double[] m) {
		if (prime) momentumPrime = m;
		else momentum = m;
	}
	public void setMomentumAtPos(boolean prime, double m, int pos) {
		if (prime) momentumPrime[pos] = m;
		else momentum[pos] = m;
	}
	
	public double[] getForce() {
		return force;
	}
	public double getForceAtPos(int pos) {
		return force[pos];
	}
	public void setForce(double[] f) {
		force = f;
	}
	public void setForceAtPos(double f, int pos) {
		force[pos] = f;
	}
}
