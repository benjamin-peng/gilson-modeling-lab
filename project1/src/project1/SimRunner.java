package project1;
import java.io.*;
import java.util.*;

public class SimRunner {
	public static void main(String[] args) throws FileNotFoundException {
		
		double mass = 3.9948e-2;
		double epsilon = 0.0661;
		double sigma = 0.3345;
		double gamma = 0.1;
		double b = 6e-5;
				//1.38064852e-23;
		double temp;
		double stepLength = 1e-3;
		int n;
		int stepNum = 30000;
		
		Scanner scan = new Scanner(System.in);
		System.out.print("what temp? :: ");
		temp = scan.nextDouble();
		System.out.print("how many atoms? :: ");
		n = scan.nextInt();
		//argon
		LangevinAlgorithm sim = new LangevinAlgorithm(mass, epsilon, sigma, gamma, b, temp, stepLength, /*thrust,*/ n, stepNum);
		//sim.setRandomCoord();
		sim.setRandomCoord();
		sim.printCoord();
		sim.setForce();
		sim.getNewCoord();
		//PrintStream console = System.out;
		//System.setOut(console); 
		for (int i = 0; i < sim.speedList.size(); i++) 
			System.out.println(sim.speedList.get(i));
	}
}
