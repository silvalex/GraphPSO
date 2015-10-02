package pso;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;


public class SinglePassGraphPSO extends GraphPSO {

	/**
	 * Application's entry point.
	 *
	 * @param args
	 */
	public static void main(String[] args) {
		new SinglePassGraphPSO(args[0], args[1], args[2], args[3], args[4], Long.valueOf(args[5]));
	}

	/**
	 * Creates a functionally correct workflow, and runs the PSO to discover the
	 * optimal services to be used in it.
	 */
	public SinglePassGraphPSO(String logName, String histogramLogName, String taskFileName, String serviceFileName, String taxonomyFileName, long seed) {
		super(logName, histogramLogName, taskFileName, serviceFileName, taxonomyFileName, seed);
	}

	//==========================================================================================================
	//
    //                                              PSO METHODS
    //
	//==========================================================================================================

//	/**
//	 * Conducts the particle swarm optimization.
//	 */
//	@Override
//	public String runPSO() {
//		// 1. Initialize the swarm
//		initializeRandomSwarm();
//
//		int i = 0;
//		Particle p;
//		Graph workflow = null;
//		long initialization = System.currentTimeMillis() - initialisationStartTime;
//
//		while (i < MAX_NUM_ITERATIONS) {
//			long startTime = System.currentTimeMillis();
//			System.out.println("ITERATION " + i);
//
//			// Go through all particles
//			for (int j = 0; j < NUM_PARTICLES; j++) {
//				System.out.println("\tPARTICLE " + j);
//				p = swarm.get(j);
//				workflow = createNewGraph(startNode.clone(), endNode.clone(), relevant, p.dimensions);
//				// 2. Evaluate fitness of particle
//				if (workflow != null) {
//				    FitnessResult result = calculateFitness(workflow);
//				    p.fitness = result.fitness;
//				    p.graphString = result.graphString;
//				}
//				else {
//				    p.fitness = 0;
//				    p.graphString = "Partially-built graph";
//				}
//				// 3. If fitness of particle is better than Pbest, update the Pbest
//				p.updatePersonalBest();
//				// 4. If fitness of Pbest is better than Gbest, update the Gbest
//				if (p.bestFitness > Particle.globalBestFitness) {
//					Particle.globalBestFitness = p.bestFitness;
//					Particle.globalGraphString = p.graphString;
//					Particle.globalBestDimensions = Arrays.copyOf(p.bestDimensions, p.bestDimensions.length);
//				}
//				// 5. Update the velocity of particle
//				updateVelocity(p);
//				// 6. Update the position of particle
//				updatePosition(p);
//			}
//
//			fitness.add(Particle.globalBestFitness);
//			time.add((System.currentTimeMillis() - startTime) + initialization);
//			initialization = 0;
//			i++;
//		}
//
//		return Particle.globalGraphString;
//	}


	//==========================================================================================================
	//
    //                                              GRAPH METHODS
    //
	//==========================================================================================================

	@Override
	public Graph createNewGraph(Node start, Node end, Set<Node> relevant, float[] weights) {

		Graph newGraph = new Graph();

		Set<String> currentEndInputs = new HashSet<String>();
		Map<String,Edge> connections = new HashMap<String,Edge>();

		// Connect start node
		connectCandidateToGraphByInputs(start, connections, newGraph, currentEndInputs);

		Set<Node> seenNodes = new HashSet<Node>();
		List<ListItem> candidateList = new ArrayList<ListItem>();

		populateCandidateList(serviceToIndexMap, relevant, candidateList, weights);
		Collections.sort(candidateList);

		finishConstructingGraph(currentEndInputs, end, candidateList, connections, newGraph, seenNodes, relevant);

        // Keep track of nodes and edges for statistics
        for (String nodeName : newGraph.nodeMap.keySet())
            addToCountMap(nodeCount, nodeName);
        for (Edge edge : newGraph.edgeList)
            addToCountMap(edgeCount, edge.toString());
//        if (newGraph.nodeMap.containsKey( "end" ))
            return newGraph;
//        else
//            return null;
	}

	@Override
	public void finishConstructingGraph(Set<String> currentEndInputs, Node end, List<ListItem> candidateList, Map<String,Edge> connections,
	        Graph newGraph, Set<Node> seenNodes, Set<Node> relevant) {

	    // While end cannot be connected to graph
	    int index = 0;
	    boolean satisfied;
		while(!(satisfied = checkCandidateNodeSatisfied(connections, newGraph, end, end.getInputs(), null))){
			connections.clear();

            candidateLoop:
            for (; index < candidateList.size(); index++) {
            	ListItem item = candidateList.get(index);
            	Node candidate = serviceMap.get(item.serviceName).clone();
                // For all of the candidate inputs, check that there is a service already in the graph
                // that can satisfy it

                if (!checkCandidateNodeSatisfied(connections, newGraph, candidate, candidate.getInputs(), null)) {
                    connections.clear();
                	continue candidateLoop;
                }

                // Connect candidate to graph, adding its reachable services to the candidate list
                connectCandidateToGraphByInputs(candidate, connections, newGraph, currentEndInputs);
                connections.clear();

                break;
            }

			if (index < candidateList.size()) {
			    candidateList.remove(index);
			}
			else {
			    break;
			}
        }

		if (satisfied) {
		    connectCandidateToGraphByInputs(end, connections, newGraph, currentEndInputs);
		    connections.clear();
		    removeDanglingNodes(newGraph);
		}
		else {
			connections.clear();
			checkCandidateNodeSatisfied(connections, newGraph, end, end.getInputs(), null);
			newGraph.connectionsToEndNode = connections.size();
		}
	}

	@Override
	public FitnessResult calculateFitness(Graph graph) {
		if (graph.nodeMap.containsKey("end")) {
			double a = 1.0;
			double r = 1.0;
			double t = 0.0;
			double c = 0.0;

			for (Node n : graph.nodeMap.values()) {
				double[] qos = n.getQos();
				a *= qos[AVAILABILITY];
				r *= qos[RELIABILITY];
				c += qos[COST];
			}

			// Calculate longest time
			t = findLongestPath(graph);

			a = normaliseAvailability(a);
			r = normaliseReliability(r);
			t = normaliseTime(t);
			c = normaliseCost(c);

			double fitness = 1.0 + (W1 * a + W2 * r + W3 * t + W4 * c);
			String graphString = graph.toString();

			return new FitnessResult(fitness, graphString);
		}
		else {
			double a = 1.0;
			double r = 1.0;
			double c = 0.0;

			for (Node n : graph.nodeMap.values()) {
				double[] qos = n.getQos();
				a *= qos[AVAILABILITY];
				r *= qos[RELIABILITY];
				c += qos[COST];
			}

			a = normaliseAvailability(a);
			r = normaliseReliability(r);
			c = normaliseCost(c);

			double fitness = (1.0 - (1.0/(graph.connectionsToEndNode + 1))) + (W1 * a + W2 * r + W4 * c); // Time cannot be calculated with an incomplete graph
			String graphString = graph.toString();

			return new FitnessResult(fitness, graphString);
		}
	}
}