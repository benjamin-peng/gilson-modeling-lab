package project1;
import java.io.*;
import java.util.*;

public class SimRunner {
	public static void main(String[] args) throws FileNotFoundException {
		
		double mass = 3.9948e-2;
		double epsilon = 0.0661;
		double sigma = 0.3345;
		double gamma = 0.1;
		double b = 1.38064852e-23;
		double temp = 300;
		double stepLength = 1e-3;
		int n;
		int stepNum = 10000;
		
		Scanner scan = new Scanner(System.in);
		System.out.print("how many atoms? :: ");
		n = scan.nextInt();
		//argon
		LangevinAlgorithm sim = new LangevinAlgorithm(mass, epsilon, sigma, gamma, b, temp, stepLength, n, stepNum);
		//System.out.println((sim.ran.nextGaussian() * .5 + .5));
		//sim.setRandomCoord();
		sim.setRandomCoord();
		sim.printCoord();
		sim.setForce();
		sim.getNewCoord();
		PrintStream console = System.out;
		System.setOut(console); 
		System.out.println(sim.speed);
	}
}
