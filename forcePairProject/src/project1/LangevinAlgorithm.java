package project1;

import java.util.*;
import java.io.*; 

import java.lang.Math;
public class LangevinAlgorithm {
	
	private double baseGamma, baseTemp, b, timestep, steps;
	private int[] nAtoms;
	private double avogadro = 6.0221409e23;
	private ArrayList<Particle> particleList;
	ArrayList<ArrayList<Integer>> LJPairList;
	ArrayList<ArrayList<Integer>> covPairList;
	ArrayList<ArrayList<Integer>> enzPairList;


	//List<Double> speedList = new ArrayList<Double>();
	//List<Double> kineticEList = new ArrayList<Double>();
	Random ran = new Random(); 
	
	// constructors
	public LangevinAlgorithm() {
		//TODO: add default constructor
	}

	public LangevinAlgorithm(double[] m, double[] ep, double[] sig, double gam, double B, double t, double time, /*double thrust,*/int[] n, double nStep) {
		nAtoms = new int[3];
		nAtoms = n.clone();
		LJPairList = new ArrayList<ArrayList<Integer>>();
		covPairList = new ArrayList<ArrayList<Integer>>();
		enzPairList = new ArrayList<ArrayList<Integer>>();
		baseGamma = gam / 10;
		baseTemp = t;
		b = B;
		timestep = time;
		steps = nStep;
		//creates a particle list where each particle has its own coord
		particleList = new ArrayList<>();
		
		double[] coord = {0, 0, 0};
		double[] coord1 = {6e-9, 6e-9, 6e-9};

		
		for (int i = 0; i < n[1]; i++) {
			particleList.add(new Particle("S", m[2], ep[2], sig[2], gam, getRandomCoord()));//getRandomCoord()));
		}
		for (int i = 0; i < n[0]; i++) {
			particleList.add(new Particle("E", m[0], ep[0], sig[0], gam, getRandomCoord()));
		}
		for (int i = 0; i < n[0]; i++) {
			double[] closeCoord = particleList.get(i + n[1]).getCoord().clone();
			
			for (int j = 0; j < 3; j++) {
				closeCoord[j] += 3.464e-9;
			}
			
			particleList.add(new Particle("A", m[1], ep[1], sig[1], gam, closeCoord));//getRandomCloseCoord(particleList.get(i + n[1]).getCoord())));
		}
		

		//0 position = enzyme #, 1 position = substrate #

	}
	
	public void setPairLists() {
		for (int i = 0; i < nAtoms[1]; i++) {
			for (int j = i + 1; j < nAtoms[1]; j++) {
				ArrayList<Integer> forcePair = new ArrayList<Integer>();
				forcePair.add(i);
				forcePair.add(j);
				LJPairList.add(forcePair);
			}
		}
		
		
		//creates LJ list with substrate and active sites, double check this works
		for (int i = 0; i < nAtoms[1]; i++) {
			for (int k = nAtoms[1] + nAtoms[0] + i; k < nAtoms[1] + nAtoms[0] * 2; k++) {
				ArrayList<Integer> forcePair = new ArrayList<Integer>();
				forcePair.add(i);
				forcePair.add(k);
				LJPairList.add(forcePair);
			}
		}
		
		
		//TODO: Fix enz pair list and add cov pair list
		for (int i = nAtoms[1]; i < nAtoms[0] + nAtoms[1]; i++) {
			ArrayList<Integer> forcePair = new ArrayList<Integer>();
			forcePair.add(i);
			forcePair.add(i + nAtoms[0]);
			enzPairList.add(forcePair);
		}
		

		
		for (int i = 0; i < nAtoms[1]; i++) {
			for (int j = nAtoms[1] + nAtoms[0]; j < 2 * nAtoms[0] + nAtoms[1]; j++) {
				ArrayList<Integer> forcePair = new ArrayList<Integer>();
				forcePair.add(i);
				forcePair.add(j);
				covPairList.add(forcePair);
			}
		}
		
		System.out.println("electrostatic pair list");
		System.out.println(Arrays.toString(covPairList.toArray()));
	}
	
	public double[] getRandomCoord() {
		double[] coord = new double[3];
		Random rand = new Random();
		for (int i = 0; i < 3; i++) {
			coord[i] = -1e-7 + (1e-7 - (-1e-7)) * rand.nextDouble();
		}
		System.out.println(Arrays.toString(coord));
		return coord;
	}
	
	public double[] getRandomCloseCoord(double[] coord) {
		for (int i = 0; i < 3; i++) {
			coord[i] += 3.464e-9;
		}
		return coord;
	}
	

	public void setRandomCoord() {
		for (Particle p : particleList) {
			p.setCoord(getRandomCoord());
		}
	}

	

	public double getLJForce(double r, double vector, double epsilon, double sigma) {
		if (r != 0.0) {
			return ((48 * epsilon)/(r*r)) * ((Math.pow((sigma / r), 12)) - .5*(Math.pow((sigma / r), 6))) * vector;
		}
		else return 0.0;
	}
	
	public double getCovForce(double q1, double q2, double r) {
		return 8.988e9 * q1 * q1 / (r * r);
	}
	
	public void setForce() {

		setNewForce();
		
		for (ArrayList<Integer> pair : LJPairList) {
			
			double[] distanceVector = new double[3];
			double r;
			double eps = Math.sqrt(particleList.get(pair.get(0)).getEpsilon() * particleList.get(pair.get(1)).getEpsilon()) ;
			double sig = (particleList.get(pair.get(0)).getSigma() + particleList.get(pair.get(1)).getSigma()) / 2;
			distanceVector[0] = particleList.get(pair.get(0)).getCoordAtPos(0) - particleList.get(pair.get(1)).getCoordAtPos(0);
			distanceVector[1] = particleList.get(pair.get(0)).getCoordAtPos(1) - particleList.get(pair.get(1)).getCoordAtPos(1);
			distanceVector[2] = particleList.get(pair.get(0)).getCoordAtPos(2) - particleList.get(pair.get(1)).getCoordAtPos(2);
			
			
		    r = Math.sqrt(Math.pow(distanceVector[0], 2) + Math.pow(distanceVector[1], 2) + Math.pow(distanceVector[2], 2));

		    particleList.get(pair.get(0)).setForceAtPos(particleList.get(pair.get(0)).getForceAtPos(0) + getLJForce(r, distanceVector[0], eps, sig), 0);
		    particleList.get(pair.get(0)).setForceAtPos(particleList.get(pair.get(0)).getForceAtPos(1) + getLJForce(r, distanceVector[1], eps, sig), 1);
		    particleList.get(pair.get(0)).setForceAtPos(particleList.get(pair.get(0)).getForceAtPos(2) + getLJForce(r, distanceVector[2], eps, sig), 2);
		    
		    particleList.get(pair.get(1)).setForceAtPos(particleList.get(pair.get(1)).getForceAtPos(0) - getLJForce(r, distanceVector[0], eps, sig), 0);
		    particleList.get(pair.get(1)).setForceAtPos(particleList.get(pair.get(1)).getForceAtPos(1) - getLJForce(r, distanceVector[1], eps, sig), 1);
		    particleList.get(pair.get(1)).setForceAtPos(particleList.get(pair.get(1)).getForceAtPos(2) - getLJForce(r, distanceVector[2], eps, sig), 2);
		    
		}
		
		
		for (ArrayList<Integer> pair : enzPairList) {
			
			double[] pairOneCoord = particleList.get(pair.get(0)).getCoord().clone();
			double[] pairTwoCoord = particleList.get(pair.get(1)).getCoord().clone();
			double k = 0.0248907744; // 6.5 ev / angstrom squared
			double r0 = 6e-9;
			double r12 = Math.sqrt(Math.pow(pairOneCoord[0] - pairTwoCoord[0], 2) + Math.pow(pairOneCoord[1] - pairTwoCoord[1], 2) + Math.pow(pairOneCoord[2] - pairTwoCoord[2], 2));
			//System.out.println(r12);
			
			for (int i = 0; i < 3; i++) {
				particleList.get(pair.get(0)).setForceAtPos(particleList.get(pair.get(0)).getForceAtPos(i) - k * (r12 - r0) * (pairOneCoord[i] - pairTwoCoord[i]) / r12, i);
				particleList.get(pair.get(1)).setForceAtPos(particleList.get(pair.get(0)).getForceAtPos(i) + k * (r12 - r0) * (pairOneCoord[i] - pairTwoCoord[i]) / r12, i);
			}
			

		}
		
		
		for (ArrayList<Integer> pair : covPairList) {
			
			
			double[] distanceVector = new double[3];
			distanceVector[0] = particleList.get(pair.get(0)).getCoordAtPos(0) - particleList.get(pair.get(1)).getCoordAtPos(0);
			distanceVector[1] = particleList.get(pair.get(0)).getCoordAtPos(1) - particleList.get(pair.get(1)).getCoordAtPos(1);
			distanceVector[2] = particleList.get(pair.get(0)).getCoordAtPos(2) - particleList.get(pair.get(1)).getCoordAtPos(2);
			
			
			double k = 8.98755e9;// * Math.pow(1.6e-19, 2);
			double q1 =  1.60218e-19;
			double q2 = -1.60218e-19;
			double r = Math.sqrt(Math.pow(distanceVector[0], 2) + Math.pow(distanceVector[1], 2) + Math.pow(distanceVector[2], 2));


		    particleList.get(pair.get(0)).setForceAtPos(particleList.get(pair.get(0)).getForceAtPos(0) + (k * q1 * q2 / (r * r) * (distanceVector[0] / r)), 0);
		    particleList.get(pair.get(0)).setForceAtPos(particleList.get(pair.get(0)).getForceAtPos(1) + (k * q1 * q2 / (r * r) * (distanceVector[1] / r)), 1);
		    particleList.get(pair.get(0)).setForceAtPos(particleList.get(pair.get(0)).getForceAtPos(2) + (k * q1 * q2 / (r * r) * (distanceVector[2] / r)), 2);

		    particleList.get(pair.get(1)).setForceAtPos(particleList.get(pair.get(1)).getForceAtPos(0) - (k * q1 * q2 / (r * r) * (distanceVector[0] / r)), 0);
		    particleList.get(pair.get(1)).setForceAtPos(particleList.get(pair.get(1)).getForceAtPos(1) - (k * q1 * q2 / (r * r) * (distanceVector[1] / r)), 1);
		    particleList.get(pair.get(1)).setForceAtPos(particleList.get(pair.get(1)).getForceAtPos(2) - (k * q1 * q2 / (r * r) * (distanceVector[2] / r)), 2);
			
		}
		

	}
	
	public double[][] getMomentum() {
		double[][] momentum = new double[nAtoms[1]][3];
		for (int i = 0; i < nAtoms[1]; i++) {
			for (int j = 0; j < 3; j++) {
				momentum[i][j] = 0.0;
			}
		}
		return momentum;
	}
	
	public double[][] getPrime() {
		double[][] momentumPrime = new double[nAtoms[1]][3];
		for (int i = 0; i < nAtoms[1]; i++) {
			for (int j = 0; j < 3; j++) {
				momentumPrime[i][j] = 0.0;
			}
		}
		return momentumPrime;
	}
	

	
	public void setNewForce() {
		for (Particle p : particleList) {
			p.setForce(new double[] {0, 0, 0});
		}
	}
	
	
	public void setSpeed(double coord[][], int nAtoms) { //sums up the magnitudes of speed and appends to the ArrayList speed
		double speed = 0.0;
		for (int i = 0; i < nAtoms; i++) {
			for (int j = 0; j < 3; j++) {
				speed += Math.sqrt(Math.pow(coord[i][0], 2) + Math.pow(coord[i][1], 2) + Math.pow(coord[i][2], 2));
			}

		} 
	}
	
	public void setNewCoord() {
		for (Particle p : particleList) {
			//setForce();
			//TODO: Fix forces w new coord system
			//force = getForce();
			for (int j = 0; j < 3; j++) {
				//restrains particles to a box of certain dimensions, assuming box collision is completely elastic
				if(p.getCoord()[j] > 1e-6) {
					p.setCoordAtPos(1e-6, j);
					p.setMomentumAtPos(false, -p.getMomentumAtPos(false, j), j);
                	//momentum[i][j] = -momentum[i][j];
                }
                
                else if(p.getCoord()[j] < -1e-6) {
					p.setCoordAtPos(-1e-6, j);
					p.setMomentumAtPos(false, -p.getMomentumAtPos(false, j), j);
                	//momentum[i][j] = -momentum[i][j];
                }
				double gamma = baseGamma;
				double temp = baseTemp;
				//langevin dynamics algorithm -- taken from Allen and Tildesley: Computer Simulation of Liquids
				p.setMomentumAtPos(false, p.getMomentumAtPos(false, j) + 0.5 * timestep * p.getForceAtPos(j), j);
				p.setCoordAtPos(p.getCoord()[j] + 0.5 * timestep * p.getMomentumAtPos(false, j) / p.getMass(), j);
                p.setMomentumAtPos(true, Math.exp(-gamma * timestep) * p.getMomentumAtPos(false, j) + Math.sqrt(1 - Math.exp(-2 * gamma * timestep)) * Math.sqrt(p.getMass() * b * temp) * (ran.nextGaussian()), j);
				p.setCoordAtPos(p.getCoord()[j] + 0.5 * timestep * p.getMomentumAtPos(false, j) / p.getMass(), j);
				p.setCoordAtPos(p.getCoord()[j] + .5 * timestep * p.getMomentumAtPos(true, j) / p.getMass(), j);
				setForce();
				p.setMomentumAtPos(false, p.getMomentumAtPos(true, j) + .5 * timestep * p.getForceAtPos(j), j);
                
				//teleports substrates to higher conc. pool and detects substrate concentration
				if (p.getCoord()[0] < -1.2e-7 && p.getType().equals("S")) {
					p.setCoord(new double[] {p.getCoord()[0] + 2e-7, -1e-7 + 2e-7 * ran.nextDouble(), -1e-7 + 2e-7 * ran.nextDouble()});
				}
				
				if (p.getCoord()[0] > 1.2e-7 && (p.getType().equals("S"))) {
					p.setCoordAtPos(1.2e-7, 0);
					p.setMomentumAtPos(false, -p.getMomentumAtPos(false, 0), 0);
					
                }
				
				//sets up a box of dimensions 1.0e-7 by 1.0e-7 in the y and x directions
				if (p.getCoord()[j] > 1.0e-7 && (p.getType().equals("S")) && j != 0) {
					p.setCoordAtPos(1.0e-7, j);
					p.setMomentumAtPos(false, -p.getMomentumAtPos(false, j), j);
					
                }
				
				if (p.getCoord()[j] < -1.0e-7 && (p.getType().equals("S")) && j != 0) {
					p.setCoordAtPos(-1.0e-7, j);
					p.setMomentumAtPos(false, -p.getMomentumAtPos(false, j), j);
					
                }
				//detects collisions to the walls for enzyme
				if (p.getCoord()[j] > 1e-7 && (p.getType().equals("E") || p.getType().equals("A"))) {
					p.setCoordAtPos(1e-7, j);
					p.setMomentumAtPos(false, -p.getMomentumAtPos(false, j), j);
					
                }
                
                else if (p.getCoord()[j] < -1e-7 && (p.getType().equals("E") || p.getType().equals("A"))) {
					p.setCoordAtPos(-1e-7, j);
					p.setMomentumAtPos(false, -p.getMomentumAtPos(false, j), j);
                	
                }
			}

		}
		
	}
	
	

	/*
	public void setPresetCoord() {
		coord = getPresetCoord();
	}
	*/
	
	//prints in .xyz file format
	public void printCoord() throws FileNotFoundException {
		File f = new File("C:\\Users\\lochn\\OneDrive\\Documents\\output\\output" + (int)baseTemp + ".xyz");
		PrintStream x = new PrintStream(new File("C:\\Users\\lochn\\OneDrive\\Documents\\output\\enzymeKE.txt"));
		PrintStream y = new PrintStream(new File("C:\\Users\\lochn\\OneDrive\\Documents\\output\\substrateKE.txt"));
		//sets up the printing of all the x positions to files
		/*
		ArrayList<PrintStream> streamList = new ArrayList<>();
		for (int i = 0; i < nAtoms[1]; i++) {
			streamList.add(new PrintStream("C:\\Users\\lochn\\OneDrive\\Documents\\output\\xpos" + i + ".txt"));
		}
		*/
		PrintStream p = new PrintStream(f);
		p.flush();
		// sets output to file
		System.setOut(p); 
		//PrintStream console = System.out;
		setForce();
		double count = 0.0;
		//double[] newCoord = new double[3];
		//for (int i = 0; i < 3; i++) newCoord[i] = particleList.get(2).getCoordAtPos(i) + 1e-9;
		//particleList.get(0).setCoord(newCoord);
		for (double k = 0.0; k < steps; k += 1.0) {
			System.setOut(p);
			//skips 4/5 frames
			if (count >= 210000) System.exit(0);
			if //(true) {
			(k % 200 == 0) {
				System.out.println(nAtoms[0] * 2 + nAtoms[1]);
				System.out.println(count);
				count += 1.0;
				for (Particle a : particleList) {
					System.out.print(a.getType());
					for (int j = 0; j < 3; j++) {
						System.out.print(" " + (a.getCoordAtPos(j)*10e7));
						if ((j + 1) % 3 == 0) System.out.print("\n");
					}
				}
			}
			System.setOut(y);
			for (int i = 0; i < nAtoms[1]; i++) {
				System.out.println(0.5 * particleList.get(i).getMass() * Math.pow(particleList.get(i).getSpeed(), 2));
			}
			System.setOut(x); 
			for (int i = nAtoms[1]; i < nAtoms[1] + 2 * nAtoms[0]; i++) {
				System.out.println(0.5 * particleList.get(i).getMass() * Math.pow(particleList.get(i).getSpeed(), 2));
			}
			setNewCoord();
		}

	}
	
	
	
}
