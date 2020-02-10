package project1;

import java.util.*;
import java.lang.Math;
public class VerletAlgorithm {
	
	private double mass, epsilon, sigma, gamma, temp, b, timestep;
	private int nAtoms, steps;
	private double[][] coord;
	private double[][] coordPrime; // temp holding array for the next timestep while preserving the last timestep's coords for calculation of force
	private double[][] momentum;
	private double[][] force;
	Random ran = new Random(); 
	
	// constructor
	public VerletAlgorithm() {
	}

	public VerletAlgorithm(double m, double ep, double sig, double gam, double B, double t, double time, int n, int nStep) {
		mass = m;
		epsilon = ep;
		sigma = sig;
		gamma = gam;
		temp = t;
		b = B;
		nAtoms = n;
		timestep = time;
		steps = nStep;
		setMomentum();
		setPrime();
	}
	
	public double[][] getRandomCoord() {
		double[][] coord = new double[nAtoms][nAtoms];
		Random rand = new Random();
		for (int i = 0; i < nAtoms; i++) {
			for (int j = 0; j < 3; j++) {
				//coord[i][j] = 0.0;
				coord[i][j] = -10.0 + (10.0 - (-10.0)) * rand.nextDouble();
			}
		}
		return coord;
	}
	
	public double[][] getPresetCoord() {
		double[][] coord = {{0.0, 0.0, 0.3}, {0.1, 0.2, -0.3}, {0.0, 0.0, 0.0}};
		return coord;
	}
	
	public double getLJForce(double r, double vector) {
		if (r != 0) {
			return ((48 * epsilon)/(r*r)) * ((Math.pow((sigma / r), 12)) - .5*(Math.pow((sigma / r), 6))) * vector;
		}
		else return 0.0;
	}
	
	public double[][] getForce() {
		double[][] force = new double[nAtoms][3];
		double[] distanceVector = new double[3];
		double r;
		
		for (int i = 0; i < nAtoms; i++) {
			force[i][0] = 0.0;
			force[i][1] = 0.0;
			force[i][2] = 0.0;
			
			for (int j = 0; j < nAtoms; j++) {
				distanceVector[0] = coord[i][0] - coord[j][0];
				distanceVector[1] = coord[i][1] - coord[j][1];
			    distanceVector[2] = coord[i][2] - coord[j][2];
			    
			    r = Math.sqrt(Math.pow(distanceVector[0], 2) + Math.pow(distanceVector[1], 2) + Math.pow(distanceVector[2], 2));
			    
			    force[i][0] += getLJForce(r, distanceVector[0]);
			    force[i][0] += getLJForce(r, distanceVector[1]);
			    force[i][0] += getLJForce(r, distanceVector[2]);
			}
		}
		
		return force;
	}
	
	public double[][] getMomentum() {
		double[][] momentum = new double[nAtoms][3];
		for (int i = 0; i < nAtoms; i++) {
			for (int j = 0; j < 3; j++) {
				momentum[i][j] = 0.0;
			}
		}
		return momentum;
	}
	
	public double[][] getPrime() {
		double[][] coordPrime = new double[nAtoms][3];
		for (int i = 0; i < nAtoms; i++) {
			for (int j = 0; j < 3; j++) {
				coordPrime[i][j] = 0.0;
			}
		}
		return coordPrime;
	}
	
	public void setMomentum() {
		momentum = getMomentum();
	}
	
	public void setForce() {
		force = getForce();
	}
	
	public void setPrime() {
		coordPrime = getPrime();
	}
	
	public double[][] getNewCoord() {
		for (int i = 0; i < nAtoms; i++) {
			setForce();
			for (int j = 0; j < 3; j++) {
				momentum[i][j] = momentum[i][j] + 0.5 * timestep * force[i][j];
                coord[i][j] = coord[i][j] + 0.5 * timestep * momentum[i][j] / mass;
                coordPrime[i][j] = Math.exp(-gamma * timestep) * momentum[i][j] + Math.sqrt(1 - Math.exp(-2 * gamma * timestep)) * Math.sqrt(mass * b * temp) * (ran.nextGaussian() * .5 + .5);
                coord[i][j] = coord[i][j] + .5 * timestep * coordPrime[i][j] / mass;
                momentum[i][j] = coordPrime[i][j]  + .5 * timestep * force[i][j];
			}
			
		}
		
		for (int i = 0; i < nAtoms; i++) {
			for (int j = 0; j < 3; j++) {
				setForce();
				momentum[i][j] = coordPrime[i][j]  + .5 * timestep * force[i][j];
			}
		}
		return coord;
	}
	
	
	public void setRandomCoord() {
		coord = getRandomCoord();
	}
	
	public void setPresetCoord() {
		coord = getPresetCoord();
	}
	
	//prints in .xyz file format
	public void printCoord() {
		for (int k = 0; k < steps; k++) {
			System.out.println(nAtoms);
			System.out.println(k);
			for (int i = 0; i < nAtoms; i++) {		
				System.out.print("Ar");
				for (int j = 0; j < nAtoms; j++) {
					System.out.print(" " + coord[i][j]);
					if ((j + 1) % 3 == 0) System.out.print("\n");
				}
			}
			getNewCoord();
		}
	}
	
	
	
}
