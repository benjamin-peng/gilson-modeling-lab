package project1;
import java.io.*;
import java.util.*;

public class SimRunner {
	public static void main(String[] args) throws FileNotFoundException {
		
		/*
		double mass = 3.9948e-2;
		double epsilon = 2.873e17;
				//0.0661;
		double sigma = 0.3345;
		double gamma = 3;//0000000000;
		double b = //6e-5;
				1.38064852e-23;
		double temp;
		double stepLength = 1e-3;
		int n;
		int stepNum = 1000;
		*/
		
		double[] mass = {6.644e-22, 1.661e-22, 1.661e-24}; //100,000 gram/mol
		double[] epsilon = {6.948e-21, 6.948e-21, 6.948e-21}; //epsilon = 1kcal/mol
		double[] sigma = {4e-9, 2.52e-9, 8e-10}; //meters 
		double gamma = 2e11;//5e-12; //seconds
		double b = 1.38064852e-23;
		double temp;
		double stepLength = 5e-14; //seconds
		
		double stepNum = //1000000.0;
				10000000000.0;
		
		int[] n = new int[2];
		Scanner scan = new Scanner(System.in);
		System.out.print("what temp? :: ");
		temp = scan.nextDouble();
		System.out.print("how many enzyme/active sites? :: ");
		n[0] = scan.nextInt();
		System.out.print("how many substrate particles? :: ");
		n[1] = scan.nextInt();
		scan.close();
		//argon
		LangevinAlgorithm sim = new LangevinAlgorithm(mass, epsilon, sigma, gamma, b, temp, stepLength, /*thrust,*/ n, stepNum);
		//sim.setRandomCoord();
		sim.setPairLists();
		sim.setRandomCoord();
		sim.printCoord();
		sim.setNewForce();
		sim.setNewCoord();

		//PrintStream console = System.out;
		//System.setOut(console); 
		//for (int i = 0; i < sim.speedList.size(); i++) 
			//System.out.println(sim.speedList.get(i));
	}
}
