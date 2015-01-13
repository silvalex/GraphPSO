package pso;


import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
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

	// personal best values
	public double bestFitness = Double.NEGATIVE_INFINITY;
	public float[] bestDimensions = new float[GraphPSO.numDimensions];

	// global best values
	public static double globalBestFitness = Double.NEGATIVE_INFINITY;
	public static float[] globalBestDimensions = new float[GraphPSO.numDimensions];
	public static Graph globalBestWorkflow;

	// overall global best values (for all runs)
	public static double overallGlobalBestFitness = Double.NEGATIVE_INFINITY;
	public static float[] overallGlobalBestDimensions = new float[GraphPSO.numDimensions];
	public static Graph overallGlobalBestWorkflow;

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
		globalBestWorkflow = null;

		// overall global best values (for all runs)
		overallGlobalBestFitness = Double.NEGATIVE_INFINITY;
		overallGlobalBestDimensions = new float[GraphPSO.numDimensions];
		overallGlobalBestWorkflow = null;
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
		StringBuilder builder = new StringBuilder();
		builder.append("Particle: ");
		for (int i = 0; i < GraphPSO.numDimensions; i++) {
			builder.append(String.format("%f ", dimensions[i]));
		}
		return builder.toString();
	}


}
