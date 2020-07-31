package project1;

import java.util.*;
import java.io.*; 

import java.lang.Math;
public class LangevinAlgorithm {
	
	private double baseGamma, baseTemp, b, timestep, steps;
	private int[] nAtoms;
	private double[][] coord;
	private double[][] momentumPrime; // temp holding array for preserving the last timestep's momentum for calculation of force
	private double[][] momentum;
	private double[][] force;
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
		baseGamma = gam;
		baseTemp = t;
		b = B;
		timestep = time;
		steps = nStep;
		setMomentum();
		setPrime();
		//creates a particle list where each particle has its own coord
		particleList = new ArrayList<>();

		for (int i = 0; i < n[1]; i++) {
			particleList.add(new Particle("S", m[2], ep[2], sig[2], gam, getRandomCoord()));
		}
		for (int i = 0; i < n[0]; i++) {
			particleList.add(new Particle("E", m[0], ep[0], sig[0], gam, getRandomCoord()));
		}
		for (int i = 0; i < n[0]; i++) {
			particleList.add(new Particle("A", m[1], ep[1], sig[1], gam, getRandomCloseCoord(particleList.get(i + n[1]).getCoord())));
		}

		//0 position = enzyme #, 1 position = substrate #

	}
	
	public void setPairLists() {
		int count = 0;
		for (int i = 0; i < nAtoms[1]; i++) {
			for (int j = i + 1; j < nAtoms[1]; j++) {
				ArrayList<Integer> forcePair = new ArrayList<Integer>();
				forcePair.add(i);
				forcePair.add(j);
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
		System.out.println(Arrays.toString(enzPairList.toArray()));
	}
	
	public double[] getRandomCoord() {
		double[] coord = new double[3];
		Random rand = new Random();
		for (int i = 0; i < 3; i++) {
			coord[i] = -1e-8 + (1e-8 - (-1e-8)) * rand.nextDouble();
		}
		System.out.println(Arrays.toString(coord));
		return coord;
	}
	
	public double[] getRandomCloseCoord(double[] coord) {
		for (int i = 0; i < 3; i++) {
			coord[i] += 6e-9;
		}
		System.out.println(Arrays.toString(coord));
		return coord;
	}
	
	public void setRandomCoord() {
		for (Particle p : particleList) {
			p.setCoord(getRandomCoord());
		}
	}
	
	public double[][] getPresetCoord() {
		double[][] coord = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.1}, {0.1, 0.2, -0.3}};
		return coord;
	}
	
	public double getLJForce(double r, double vector, double epsilon, double sigma) {
		if (r != 0.0) {
			return ((48 * epsilon)/(r*r)) * ((Math.pow((sigma / r), 12)) - .5*(Math.pow((sigma / r), 6))) * vector;
		}
		else return 0.0;
	}
	
	public void setForce() {
		double[] distanceVector = new double[3];
		double r;
		
		setNewForce();
		
		for (ArrayList<Integer> pair : LJPairList) {
			double eps = particleList.get(pair.get(0)).getEpsilon();
			double sig = particleList.get(pair.get(0)).getSigma();
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
			double[] pairOneCoord = particleList.get(pair.get(0)).getCoord();
			double[] pairTwoCoord = particleList.get(pair.get(1)).getCoord();
			double k = 6.5;
			double r0 = 6e-9;
			double r12 = Math.pow(Math.pow(pairOneCoord[0] - pairTwoCoord[0], 2) + Math.pow(pairOneCoord[1] - pairTwoCoord[1], 2) + Math.pow(pairOneCoord[2] - pairTwoCoord[2], 2), 0.5);
			
			for (int i = 0; i < 3; i++) {
				particleList.get(pair.get(0)).setForceAtPos(-k * (r12 - r0) * (pairOneCoord[i] - pairTwoCoord[i]) / r12, i);
				particleList.get(pair.get(1)).setForceAtPos(k * (r12 - r0) * (pairOneCoord[i] - pairTwoCoord[i]) / r12, i);
			}

		}
		
		
		/*
		for (Particle p : particleList) {
			p.setForceAtPos(0.0, 0);
			p.setForceAtPos(0.0, 1);
			p.setForceAtPos(0.0, 2);
			
			for (Particle j : particleList) {
				
				distanceVector[0] = p.getCoordAtPos(0) - j.getCoordAtPos(0);
				distanceVector[1] = p.getCoordAtPos(1) - j.getCoordAtPos(1);
			    distanceVector[2] = p.getCoordAtPos(2) - j.getCoordAtPos(2);			    
			    
			    r = Math.sqrt(Math.pow(distanceVector[0], 2) + Math.pow(distanceVector[1], 2) + Math.pow(distanceVector[2], 2));
			    

			    
			    
			    
			    // TODO: fix LJ forces, make sure this works
			    
			    force[i][0] += getLJForce(r, distanceVector[0]);
			    force[i][1] += getLJForce(r, distanceVector[1]);
			    force[i][2] += getLJForce(r, distanceVector[2]);
			    
			    
			    
			    double forceMagnitude = Math.sqrt(Math.pow(force[i][0], 2) + Math.pow(force[i][1], 2) + Math.pow(force[i][2], 2)); // finds the magnitude of force for the atom
			    if (Math.random() * 10 >= coord[0][i]) {
			    	force[i][0] += force[i][0] * thrust / forceMagnitude;
			    	force[i][1] += force[i][1] * thrust / forceMagnitude;
			    	force[i][2] += force[i][2] * thrust / forceMagnitude;
			    { 
			    
			    
			}
		}
		//scales LJ forces appropriately
		/*
		for (int i = 0; i < nAtoms; i++) {
			for (int j = 0; j < 3; j++) {
				force[i][j] /= 1e65;
			}
		}
		*/
		//System.out.println(Arrays.deepToString(force));
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
	
	public void setMomentum() {
		momentum = getMomentum();
	}
	
	public void setNewForce() {
		for (Particle p : particleList) {
			p.setForce(new double[] {0, 0, 0});
		}
		//force = getForce();
	}
	
	public void setPrime() {
		momentumPrime = getPrime();
	}
	
	public void setSpeed(double coord[][], int nAtoms) { //sums up the magnitudes of speed and appends to the ArrayList speed
		double speed = 0.0;
		for (int i = 0; i < nAtoms; i++) {
			for (int j = 0; j < 3; j++) {
				speed += Math.sqrt(Math.pow(coord[i][0], 2) + Math.pow(coord[i][1], 2) + Math.pow(coord[i][2], 2));
			}
			//speedList.add(speed);
			//kineticEList.add(0.5*speed*speed*mass);
			//speed = 0.0;
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
			}

			//kineticEList.add(0.5*mass*Math.sqrt(Math.pow(momentum[i][0]/mass, 2) + Math.pow(momentum[i][1]/mass, 2) + Math.pow(momentum[i][2]/mass, 2)));
		}
		
		/*for (int i = 0; i < nAtoms; i++) {
			for (int j = 0; j < 3; j++) {
				setForce();
				momentum[i][j] = momentumPrime[i][j]  + .5 * timestep * force[i][j];
			}
		}*/
	}
	
	

	
	public void setPresetCoord() {
		coord = getPresetCoord();
	}
	
	//prints in .xyz file format
	public void printCoord() throws FileNotFoundException {
		File f = new File("C:\\Users\\lochn\\OneDrive\\Documents\\output\\output" + (int)baseTemp + ".xyz");
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
		setMomentum();
		setForce();
		double count = 0.0;
		for (double k = 0.0; k < steps; k += 1.0) {
			System.setOut(p);
			//skips 4/5 frames
			if (k % 2000 == 0) {
				System.out.println(nAtoms[0] * 2 + nAtoms[1]);
				System.out.println(count);
				count += 1.0;
				//System.out.println(k);
				for (Particle a : particleList) {
					//System.setOut(streamList.get(i));
					//System.out.println(p.getCoordAtPos(0));
							//coord[i][0]);
					//System.setOut(p);
					System.out.print(a.getType());
					for (int j = 0; j < 3; j++) {
						System.out.print(" " + (a.getCoordAtPos(j)*10e7));
						if ((j + 1) % 3 == 0) System.out.print("\n");
					}
				}
				//setSpeed(coord, nAtoms);
			}
			setNewCoord();
		}
		//double total = 0.0;
		/*for (int i = 0; i < kineticEList.size(); i++) {
			total += kineticEList.get(i);
		}
		double avg = total / steps;
		System.out.println(avg); */
		//System.setOut(console); // sets output back to console and prints when done
		//System.out.println("Done!");
	}
	
	
	
}
