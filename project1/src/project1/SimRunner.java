package project1;
import java.util.*;

public class SimRunner {
	public static void main(String[] args) {
		
		double mass = 39.948;
		double epsilon = 0.0661;
		double sigma = 0.3345;
		double gamma = 0.1;
		double b = 1.38064852e-23;
		double temp = 1000;
		double stepLength = 0.5;
		int n;
		int stepNum = 500;
		
		Scanner scan = new Scanner(System.in);
		System.out.print("how many atoms? :: ");
		n = scan.nextInt();
		//argon
		VerletAlgorithm sim = new VerletAlgorithm(mass, epsilon, sigma, gamma, b, temp, stepLength, n, stepNum);
		//sim.setRandomCoord();
		sim.setPresetCoord();
		sim.printCoord();
		sim.setForce();
		sim.getNewCoord();
		
		
	}
}
