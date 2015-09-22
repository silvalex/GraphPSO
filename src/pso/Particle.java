package pso;


import java.util.Arrays;
import java.util.Random;

/**
 * Represents a single particle within the swarm, including
 * its dimensions and its fitness.
 *
 * @author sawczualex
 */
public class Particle {
	//public static List<Character> dimensionsLabel = new ArrayList<Character>();
	public float[] dimensions = new float[GraphPSO.numDimensions];
	public float[] velocity = new float[GraphPSO.numDimensions];
	public double fitness = 0.0; // The higher, the fitter
	public String graphString;

	// personal best values
	public double bestFitness = Double.NEGATIVE_INFINITY;
	public float[] bestDimensions = new float[GraphPSO.numDimensions];

	// global best values
	public static double globalBestFitness = Double.NEGATIVE_INFINITY;
	public static float[] globalBestDimensions = new float[GraphPSO.numDimensions];
	public static String globalGraphString;

	/**
	 * Creates a particle with null dimensions.
	 *
	 * @param serviceMap
	 */
	public Particle(Random random) {
		for (int i=0; i < dimensions.length; i++) {
			dimensions[i] = random.nextFloat();
		}
		Arrays.fill(bestDimensions, 0.0f);
		Arrays.fill(globalBestDimensions, 0.0f);
	}

	/**
	 * Resets static fields.
	 */
	public static void reset() {
		// global best values
		globalBestFitness = Double.NEGATIVE_INFINITY;
		globalBestDimensions = new float[GraphPSO.numDimensions];
		globalGraphString = null;
	}

	/**
	 * Checks if the current solution is fitter than the overall
	 * personal best. If it is, set personal best as values from
	 * current solution.
	 */
	public void updatePersonalBest() {
		if (fitness > bestFitness) {
			bestFitness = fitness;
			bestDimensions = Arrays.copyOf(dimensions, dimensions.length);
		}
	}
	@Override
	/**
	 * {@inheritDoc}
	 */
	public String toString() {
		return graphString;
	}


}
