package pso;

/**
 * Stores the result of the fitness calculation process of a candidate,
 * including a string representation of the solution graph.
 * 
 * @author sawczualex
 *
 */
public class FitnessResult {
	public double fitness;
	public String graphString;

	/**
	 * Create a new object holding the fitness value for a given candidate and
	 * a string representation of the corresponding solution graph.
	 * 
	 * @param fitness
	 * @param graphString
	 */
	public FitnessResult(double fitness, String graphString) {
		this.fitness = fitness;
		this.graphString = graphString;
	}
}
