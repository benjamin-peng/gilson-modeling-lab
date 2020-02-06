package project1;

import java.util.*;
public class VerletAlgorithm {
	
	private double mass, epsilon, sigma, gamma, temp, nAtoms;
	
	// constructor
	public VerletAlgorithm(double m, double ep, double sig, double gam, double t, double n) {
		mass = m;
		epsilon = ep;
		sigma = sig;
		gamma = gam;
		temp = t;
		nAtoms = n;
	}
	
	public double[][] randomCoord(int nAtoms) {
		Random rand = new Random();
		double[][] coord = new double[nAtoms][nAtoms];
		for (int i = 0; i < nAtoms; i++) {
			for (int j = 0; j < nAtoms; j++) {
				coord[i][j] = -10.0 + (10.0 - (-10.0)) * rand.nextDouble();
			}
		}
		return coord[][];
	}
	
	
	
}
