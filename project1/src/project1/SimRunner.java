package project1;
import java.util.*;

public class SimRunner {
	public static void main(String[] args) {
		
		int n;
		int timestep = 10;
		
		Scanner scan = new Scanner(System.in);
		System.out.print("how many atoms? :: ");
		n = scan.nextInt();
		//argon
		VerletAlgorithm sim = new VerletAlgorithm(39.948, 0.0661, 0.3345, .1, 1.38064852e-23, n, timestep);
		sim.setRandomCoord();
		sim.printCoord(timestep);
		sim.getForce();
		
		
		
		
	}
}
