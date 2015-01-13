package pso;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.Random;
import java.util.Set;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;
import org.w3c.dom.Text;
import org.xml.sax.SAXException;


public class GraphPSO {
	// PSO settings
	private List<Particle> swarm = new ArrayList<Particle>();
	public static final int MAX_NUM_ITERATIONS = 30;
	public static final int NUM_PARTICLES = 200;
	public static final float C1 = 1;
	public static final float C2 = 1;
	public static final float W = 1;
	public static int numDimensions;
	public static long[] time = new long[MAX_NUM_ITERATIONS];
	public static double[] fitness = new double[MAX_NUM_ITERATIONS];
	public static String logName;
	public static Long initialisationStartTime;

	// Fitness function weights
	public static final double W1 = 0.25;
	public static final double W2 = 0.25;
	public static final double W3 = 0.25;
	public static final double W4 = 0.25;

	// Values for normalisation
	public final double minAvailability = 0.0;
	public double maxAvailability = -1.0;
	public final double minReliability = 0.0;
	public double maxReliability = -1.0;
	public double minTime = Double.MAX_VALUE;
	public double maxTime = -1.0;
	public double minCost = Double.MAX_VALUE;
	public double maxCost = -1.0;

	// Constants with of order of QoS attributes
	public static final int TIME = 0;
	public static final int COST = 1;
	public static final int AVAILABILITY = 2;
	public static final int RELIABILITY = 3;

	public Map<String, Node> serviceMap = new HashMap<String, Node>();
	public Set<Node> relevant;
	public Map<String, Integer> serviceToIndexMapping = new HashMap<String, Integer>();
	public Map<String, TaxonomyNode> taxonomyMap = new HashMap<String, TaxonomyNode>();
	public Set<String> taskInput;
	public Set<String> taskOutput;
	public Node startNode;
	public Node endNode;
	private Random random;
	private Graph masterGraph;


	/**
	 * Application's entry point.
	 *
	 * @param args
	 */
	public static void main(String[] args) {
		new GraphPSO(args[0], args[1], args[2], args[3], Long.valueOf(args[4]));
	}

	/**
	 * Creates a functionally correct workflow, and runs the PSO to discover the
	 * optimal services to be used in it.
	 */
	public GraphPSO(String logName, String taskFileName, String serviceFileName, String taxonomyFileName, long seed) {
		initialisationStartTime = System.currentTimeMillis();
		this.logName = logName;
		random = new Random(seed);

		parseWSCServiceFile(serviceFileName);
		parseWSCTaskFile(taskFileName);
		parseWSCTaxonomyFile(taxonomyFileName);
		findConceptsForInstances();

		double[] mockQos = new double[4];
		mockQos[TIME] = 0;
		mockQos[COST] = 0;
		mockQos[AVAILABILITY] = 1;
		mockQos[RELIABILITY] = 1;
		Set<String> startOutput = new HashSet<String>();
		startOutput.addAll(taskInput);
		startNode = new Node("start", mockQos, new HashSet<String>(), taskInput);
		endNode = new Node("end", mockQos, taskOutput ,new HashSet<String>());

		populateTaxonomyTree();
		relevant = getRelevantServices(serviceMap, taskInput, taskOutput);
		numDimensions = relevant.size() + 1;

		int i = 0;
		for (Node r : relevant) {
			serviceToIndexMapping.put(r.getName(), i++);
		}
		serviceToIndexMapping.put(endNode.getName(), i++);

		calculateNormalisationBounds(relevant);
		masterGraph = createMasterGraph(relevant);

		runPSO();
		writeLog();
	}

	private Graph createMasterGraph(Set<Node> relevant) {
		Graph masterGraph = new Graph();
		Node start = startNode.clone();
		Node end   = endNode.clone();

		Set<String> currentEndInputs = new HashSet<String>();
		Map<String,Edge> connections = new HashMap<String,Edge>();

		// Connect start node
		connectCandidateToGraphByInputs(start, connections, masterGraph, currentEndInputs);

		Set<Node> seenNodes = new HashSet<Node>();
		List<Node> candidateList = new ArrayList<Node>();


		addToCandidateList(start, seenNodes, relevant, candidateList);

		while(!candidateList.isEmpty()) {

			// Select node
			int index;

			candidateLoop:
			for (index = 0; index < candidateList.size(); index++) {
				Node candidate = candidateList.get(index).clone();
				// For all of the candidate inputs, check that there is a service already in the graph
				// that can satisfy it
				connections.clear();

				for (String input : candidate.getInputs()) {
					boolean found = false;
					 for (Node s : taxonomyMap.get(input).servicesWithOutput) {
						 if (masterGraph.nodeMap.containsKey(s.getName())) {
							 Set<String> intersect = new HashSet<String>();
							 intersect.add(input);

							 Edge mapEdge = connections.get(s.getName());
							 if (mapEdge == null) {
								 Edge e = new Edge(intersect);
								 e.setFromNode(masterGraph.nodeMap.get(s.getName()));
								 e.setToNode(candidate);
								 connections.put(e.getFromNode().getName(), e);
							 }
							 else
								 mapEdge.getIntersect().addAll(intersect);

							 found = true;
						 }
					 }
					 // If that input cannot be satisfied, move on to another candidate node to connect
					 if (!found) {
						 // Move on to another candidate
						 continue candidateLoop;
					 }
				}

				// Connect candidate to graph, adding its reachable services to the candidate list
				connectCandidateToGraphByInputs(candidate, connections, masterGraph, currentEndInputs);
				addToCandidateList(candidate, seenNodes, relevant, candidateList);

				break;
			}

			candidateList.remove(index);
		}

		// Connect end node to graph
		connections.clear();

		for (Node s : masterGraph.nodeMap.values()) {
			for (String o : s.getOutputs()) {
				Set<String> endNodeInputs = taxonomyMap.get(o).endNodeInputs;
				if (!endNodeInputs.isEmpty()) {
					Edge e = new Edge(endNodeInputs);
					e.setFromNode(s);
					e.setToNode(end);
					connections.put(e.getFromNode().getName(), e);
				}
			}
		}

		connectCandidateToGraphByInputs(end, connections, masterGraph, currentEndInputs);
		removeDanglingNodes(masterGraph);
		return masterGraph;
	}

	/**
	 * Goes through the service list and retrieves only those services which
	 * could be part of the composition task requested by the user.
	 *
	 * @param serviceMap
	 * @return relevant services
	 */
	private Set<Node> getRelevantServices(Map<String,Node> serviceMap, Set<String> inputs, Set<String> outputs) {
		// Copy service map values to retain original
		Collection<Node> services = new ArrayList<Node>(serviceMap.values());

		Set<String> cSearch = new HashSet<String>(inputs);
		Set<Node> sSet = new HashSet<Node>();
		Set<Node> sFound = discoverService(services, cSearch);
		while (!sFound.isEmpty()) {
			sSet.addAll(sFound);
			services.removeAll(sFound);
			for (Node s: sFound) {
				cSearch.addAll(s.getOutputs());
			}
			sFound.clear();
			sFound = discoverService(services, cSearch);
		}

		if (isSubsumed(outputs, cSearch)) {
			return sSet;
		}
		else {
			String message = "It is impossible to perform a composition using the services and settings provided.";
			System.out.println(message);
			System.exit(0);
			return null;
		}
	}

	/**
	 * Discovers all services from the provided collection whose
	 * input can be satisfied either (a) by the input provided in
	 * searchSet or (b) by the output of services whose input is
	 * satisfied by searchSet (or a combination of (a) and (b)).
	 *
	 * @param services
	 * @param searchSet
	 * @return set of discovered services
	 */
	private Set<Node> discoverService(Collection<Node> services, Set<String> searchSet) {
		Set<Node> found = new HashSet<Node>();
		for (Node s: services) {
			if (isSubsumed(s.getInputs(), searchSet))
				found.add(s);
		}
		return found;
	}

	/**
	 * Conducts the particle swarm optimization.
	 */
	private void runPSO() {
		// 1. Initialize the swarm
		initializeRandomSwarm();

		int i = 0;
		Particle p;
		Graph workflow;
		long initialization = System.currentTimeMillis() - initialisationStartTime;

		while (i < MAX_NUM_ITERATIONS) {
			long startTime = System.currentTimeMillis();
			System.out.println("ITERATION " + i);

			// Go through all particles
			for (int j = 0; j < NUM_PARTICLES; j++) {
				System.out.println("\tPARTICLE " + j);
				p = swarm.get(j);
				workflow = extractGraph(p);
				// 2. Evaluate fitness of particle
				p.fitness = calculateFitness(workflow);
				// 3. If fitness of particle is better than Pbest, update the
				// Pbest
				p.updatePersonalBest();
				// 4. If fitness of Pbest is better than Gbest, update the Gbest
				if (p.bestFitness > Particle.globalBestFitness) {
					Particle.globalBestFitness = p.bestFitness;
					Particle.globalBestWorkflow = workflow;
					Particle.globalBestDimensions = Arrays.copyOf(
							p.bestDimensions, p.bestDimensions.length);

					if (Particle.globalBestFitness > Particle.overallGlobalBestFitness) {
						Particle.overallGlobalBestFitness = Particle.globalBestFitness;
						Particle.overallGlobalBestWorkflow = Particle.globalBestWorkflow;
						Particle.overallGlobalBestDimensions = Arrays.copyOf(
								Particle.globalBestDimensions,
								Particle.globalBestDimensions.length);
					}
				}
				// 5. Update the velocity of particle
				updateVelocity(p);
				// 6. Update the position of particle
				updatePosition(p);
			}

			fitness[i] = Particle.globalBestFitness;
			time[i] = (System.currentTimeMillis() - startTime) + initialization;
			initialization = 0;
			i++;
		}
	}

	/**
	 * Updates the velocity vector of a particle.
	 *
	 * @param p
	 */
	private void updateVelocity(Particle p) {
		float[] vel = p.velocity;
		float[] dim = p.dimensions;
		float[] bestDim = p.bestDimensions;
		float[] globalBestDim = Particle.globalBestDimensions;

		for (int i = 0; i < vel.length; i++) {
			vel[i] = (W * vel[i])
					+ (C1 * random.nextFloat() * (bestDim[i] - dim[i]))
					+ (C2 * random.nextFloat() * (globalBestDim[i] - dim[i]));
		}
	}

	/**
	 * Initialises the swarm with random positions and velocities.
	 */
	private void initializeRandomSwarm() {
		swarm.clear();
		for (int i = 0; i < NUM_PARTICLES; i++) {
			swarm.add(new Particle(random));
		}
	}

	/**
	 * Updates the position (i.e. dimension vector) of a particle.
	 *
	 * @param p
	 */
	private void updatePosition(Particle p) {
		float newValue;
		for (int i = 0; i < numDimensions; i++) {
			// Calculate new position for that dimension
			newValue = p.dimensions[i] + p.velocity[i];
			// Ensure new position is within bounds
			if (newValue < 0.0)
				newValue = 0.0f;
			else if (newValue > 1.0)
				newValue = 1.0f;
			// Update dimension array with new value
			p.dimensions[i] = newValue;
		}
	}

	/**
	 * Calculates the fitness of a candidate particle in the swarm.
	 *
	 * @param p
	 * @return fitness
	 */
	private double calculateFitness(Graph workflow) {
		double a = 1.0;
		double r = 1.0;
		double t = 0.0;
		double c = 0.0;

		for (Node n : workflow.nodeMap.values()) {
			double[] qos = n.getQos();
			a *= qos[AVAILABILITY];
			r *= qos[RELIABILITY];
			c += qos[COST];
		}

		// Calculate longest time
		t = findLongestPath(workflow);

		a = normaliseAvailability(a);
		r = normaliseReliability(r);
		t = normaliseTime(t);
		c = normaliseCost(c);

		return W1 * a + W2 * r + W3 * t + W4 * c;
	}

	private double normaliseAvailability(double availability) {
		if (maxAvailability - minAvailability == 0.0)
			return 1.0;
		else
			return (availability - minAvailability)/(maxAvailability - minAvailability);
	}

	private double normaliseReliability(double reliability) {
		if (maxReliability - minReliability == 0.0)
			return 1.0;
		else
			return (reliability - minReliability)/(maxReliability - minReliability);
	}

	private double normaliseTime(double time) {
		if (maxTime - minTime == 0.0)
			return 1.0;
		else
			return (maxTime - time)/(maxTime - minTime);
	}

	private double normaliseCost(double cost) {
		if (maxCost - minCost == 0.0)
			return 1.0;
		else
			return (maxCost - cost)/(maxCost - minCost);
	}

	/**
	 * Uses the Bellman-Ford algorithm with negative weights to find the longest
	 * path in an acyclic directed graph.
	 *
	 * @param g
	 * @return list of edges composing longest path
	 */
	private double findLongestPath(Graph g) {
		Map<String, Double> distance = new HashMap<String, Double>();
		Map<String, Node> predecessor = new HashMap<String, Node>();

		// Step 1: initialize graph
		for (Node node : g.nodeMap.values()) {
			if (node.getName().equals("start"))
				distance.put(node.getName(), 0.0);
			else
				distance.put(node.getName(), Double.POSITIVE_INFINITY);
		}

		// Step 2: relax edges repeatedly
		for (int i = 1; i < g.nodeMap.size(); i++) {
			for (Edge e : g.edgeList) {
				if ((distance.get(e.getFromNode().getName()) - e.getToNode().getQos()[TIME]) < distance.get(e.getToNode().getName())) {
					distance.put(e.getToNode().getName(), (distance.get(e.getFromNode().getName()) - e.getToNode().getQos()[TIME]));
					predecessor.put(e.getToNode().getName(), e.getFromNode());
				}
			}
		}

		// Now retrieve total cost
		Node pre = predecessor.get("end");
		double totalTime = 0.0;

		while (pre != null) {
			totalTime += pre.getQos()[TIME];
			pre = predecessor.get(pre.getName());
		}

		return totalTime;
	}

	/**
	 * Parses the WSC Web service file with the given name, creating Web
	 * services based on this information and saving them to the service map.
	 *
	 * @param fileName
	 */
	private void parseWSCServiceFile(String fileName) {
        Set<String> inputs = new HashSet<String>();
        Set<String> outputs = new HashSet<String>();
        double[] qos = new double[4];

        try {
        	File fXmlFile = new File(fileName);
        	DocumentBuilderFactory dbFactory = DocumentBuilderFactory.newInstance();
        	DocumentBuilder dBuilder = dbFactory.newDocumentBuilder();
        	Document doc = dBuilder.parse(fXmlFile);

        	NodeList nList = doc.getElementsByTagName("service");

        	for (int i = 0; i < nList.getLength(); i++) {
        		org.w3c.dom.Node nNode = nList.item(i);
        		Element eElement = (Element) nNode;

        		String name = eElement.getAttribute("name");
				qos[TIME] = Double.valueOf(eElement.getAttribute("Res"));
				qos[COST] = Double.valueOf(eElement.getAttribute("Pri"));
				qos[AVAILABILITY] = Double.valueOf(eElement.getAttribute("Ava"));
				qos[RELIABILITY] = Double.valueOf(eElement.getAttribute("Rel"));

				// Get inputs
				org.w3c.dom.Node inputNode = eElement.getElementsByTagName("inputs").item(0);
				NodeList inputNodes = ((Element)inputNode).getElementsByTagName("instance");
				for (int j = 0; j < inputNodes.getLength(); j++) {
					org.w3c.dom.Node in = inputNodes.item(j);
					Element e = (Element) in;
					inputs.add(e.getAttribute("name"));
				}

				// Get outputs
				org.w3c.dom.Node outputNode = eElement.getElementsByTagName("outputs").item(0);
				NodeList outputNodes = ((Element)outputNode).getElementsByTagName("instance");
				for (int j = 0; j < outputNodes.getLength(); j++) {
					org.w3c.dom.Node out = outputNodes.item(j);
					Element e = (Element) out;
					outputs.add(e.getAttribute("name"));
				}

                Node ws = new Node(name, qos, inputs, outputs);
                serviceMap.put(name, ws);
                inputs = new HashSet<String>();
                outputs = new HashSet<String>();
                qos = new double[4];
        	}
        }
        catch(IOException ioe) {
            System.out.println("Service file parsing failed...");
        }
        catch (ParserConfigurationException e) {
            System.out.println("Service file parsing failed...");
		}
        catch (SAXException e) {
            System.out.println("Service file parsing failed...");
		}
    }

	/**
	 * Parses the WSC task file with the given name, extracting input and
	 * output values to be used as the composition task.
	 *
	 * @param fileName
	 */
	private void parseWSCTaskFile(String fileName) {
		try {
	    	File fXmlFile = new File(fileName);
	    	DocumentBuilderFactory dbFactory = DocumentBuilderFactory.newInstance();
	    	DocumentBuilder dBuilder = dbFactory.newDocumentBuilder();
	    	Document doc = dBuilder.parse(fXmlFile);

	    	org.w3c.dom.Node provided = doc.getElementsByTagName("provided").item(0);
	    	NodeList providedList = ((Element) provided).getElementsByTagName("instance");
	    	taskInput = new HashSet<String>();
	    	for (int i = 0; i < providedList.getLength(); i++) {
				org.w3c.dom.Node item = providedList.item(i);
				Element e = (Element) item;
				taskInput.add(e.getAttribute("name"));
	    	}

	    	org.w3c.dom.Node wanted = doc.getElementsByTagName("wanted").item(0);
	    	NodeList wantedList = ((Element) wanted).getElementsByTagName("instance");
	    	taskOutput = new HashSet<String>();
	    	for (int i = 0; i < wantedList.getLength(); i++) {
				org.w3c.dom.Node item = wantedList.item(i);
				Element e = (Element) item;
				taskOutput.add(e.getAttribute("name"));
	    	}
		}
		catch (ParserConfigurationException e) {
            System.out.println("Task file parsing failed...");
		}
		catch (SAXException e) {
            System.out.println("Task file parsing failed...");
		}
		catch (IOException e) {
            System.out.println("Task file parsing failed...");
		}
	}

	public void writeLog() {
		try {
			FileWriter writer = new FileWriter(new File(logName));
			for (int i = 0; i < MAX_NUM_ITERATIONS; i++) {
				writer.append(String.format("%d %d %f\n", i, time[i], fitness[i]));
			}
			writer.close();
		}
		catch (IOException e) {
			e.printStackTrace();
		}
	}

	/**
	 * Parses the WSC taxonomy file with the given name, building a
	 * tree-like structure.
	 *
	 * @param fileName
	 */
	private void parseWSCTaxonomyFile(String fileName) {
		try {
	    	File fXmlFile = new File(fileName);
	    	DocumentBuilderFactory dbFactory = DocumentBuilderFactory.newInstance();
	    	DocumentBuilder dBuilder = dbFactory.newDocumentBuilder();
	    	Document doc = dBuilder.parse(fXmlFile);
	    	Element taxonomy = (Element) doc.getChildNodes().item(0);

	    	processTaxonomyChildren(null, taxonomy.getChildNodes());
		}

		catch (ParserConfigurationException e) {
            System.err.println("Taxonomy file parsing failed...");
		}
		catch (SAXException e) {
            System.err.println("Taxonomy file parsing failed...");
		}
		catch (IOException e) {
            System.err.println("Taxonomy file parsing failed...");
		}
	}

	/**
	 * Recursive function for recreating taxonomy structure from file.
	 *
	 * @param parent - Nodes' parent
	 * @param nodes
	 */
	private void processTaxonomyChildren(TaxonomyNode parent, NodeList nodes) {
		if (nodes != null && nodes.getLength() != 0) {
			for (int i = 0; i < nodes.getLength(); i++) {
				org.w3c.dom.Node ch = nodes.item(i);

				if (!(ch instanceof Text)) {
					Element currNode = (Element) nodes.item(i);

					TaxonomyNode taxNode = new TaxonomyNode(currNode.getAttribute("name"), parent);
					taxonomyMap.put(taxNode.value, taxNode);
					if (parent != null)
						parent.children.add(taxNode);

					NodeList children = currNode.getChildNodes();
					processTaxonomyChildren(taxNode, children);
				}
			}
		}
	}

	private void calculateNormalisationBounds(Set<Node> services) {
		for(Node service: services) {
			double[] qos = service.getQos();

			// Availability
			double availability = qos[AVAILABILITY];
			if (availability > maxAvailability)
				maxAvailability = availability;

			// Reliability
			double reliability = qos[RELIABILITY];
			if (reliability > maxReliability)
				maxReliability = reliability;

			// Time
			double time = qos[TIME];
			if (time > maxTime)
				maxTime = time;
			if (time < minTime)
				minTime = time;

			// Cost
			double cost = qos[COST];
			if (cost > maxCost)
				maxCost = cost;
			if (cost < minCost)
				minCost = cost;
		}
		// Adjust max. cost and max. time based on the number of services in shrunk repository
		maxCost *= services.size();
		maxTime *= services.size();

	}

	/**
	 * Populates the taxonomy tree by associating services to the
	 * nodes in the tree.
	 */
	private void populateTaxonomyTree() {
		for (Node s: serviceMap.values())
			addServiceToTaxonomyTree(s);

		// Add input and output nodes
		addServiceToTaxonomyTree(startNode);
		addEndNodeToTaxonomyTree();
	}

	private void addServiceToTaxonomyTree(Node s) {
		// Populate outputs
		for (String outputVal : s.getOutputs()) {
			TaxonomyNode n = taxonomyMap.get(outputVal);
			n.servicesWithOutput.add(s);
			s.getTaxonomyOutputs().add(n);

			// Also add output to all parent nodes
			TaxonomyNode current = n.parent;
			while (current != null) {
				current.servicesWithOutput.add(s);
				s.getTaxonomyOutputs().add(current);
				current = current.parent;
			}

		}
		// Populate inputs
		for (String inputVal : s.getInputs()) {
			TaxonomyNode n = taxonomyMap.get(inputVal);
			n.servicesWithInput.add(s);

			// Also add input to all children nodes
			Queue<TaxonomyNode> queue = new LinkedList<TaxonomyNode>();
			queue.addAll(n.children);

			while(!queue.isEmpty()) {
				TaxonomyNode current = queue.poll();
				current.servicesWithInput.add(s);
				queue.addAll(current.children);
			}
		}
	}

	private void addEndNodeToTaxonomyTree() {
		for (String inputVal : endNode.getInputs()) {
			TaxonomyNode n = taxonomyMap.get(inputVal);
			n.endNodeInputs.add(inputVal);

			// Also add input to all children nodes
			Queue<TaxonomyNode> queue = new LinkedList<TaxonomyNode>();
			queue.addAll(n.children);

			while(!queue.isEmpty()) {
				TaxonomyNode current = queue.poll();
				current.endNodeInputs.add(inputVal);
				queue.addAll(current.children);
			}
		}

	}

	/**
	 * Converts input, output, and service instance values to their corresponding
	 * ontological parent.
	 */
	private void findConceptsForInstances() {
		Set<String> temp = new HashSet<String>();

		for (String s : taskInput)
			temp.add(taxonomyMap.get(s).parent.value);
		taskInput.clear();
		taskInput.addAll(temp);

		temp.clear();
		for (String s : taskOutput)
				temp.add(taxonomyMap.get(s).parent.value);
		taskOutput.clear();
		taskOutput.addAll(temp);

		for (Node s : serviceMap.values()) {
			temp.clear();
			Set<String> inputs = s.getInputs();
			for (String i : inputs)
				temp.add(taxonomyMap.get(i).parent.value);
			inputs.clear();
			inputs.addAll(temp);

			temp.clear();
			Set<String> outputs = s.getOutputs();
			for (String o : outputs)
				temp.add(taxonomyMap.get(o).parent.value);
			outputs.clear();
			outputs.addAll(temp);
		}
	}

	/**
	 * Checks whether set of inputs can be completely satisfied by the search
	 * set, making sure to check descendants of input concepts for the subsumption.
	 *
	 * @param inputs
	 * @param searchSet
	 * @return true if search set subsumed by input set, false otherwise.
	 */
	public boolean isSubsumed(Set<String> inputs, Set<String> searchSet) {
		boolean satisfied = true;
		for (String input : inputs) {
			Set<String> subsumed = taxonomyMap.get(input).getSubsumedConcepts();
			if (!isIntersection( searchSet, subsumed )) {
				satisfied = false;
				break;
			}
		}
		return satisfied;
	}

    private static boolean isIntersection( Set<String> a, Set<String> b ) {
        for ( String v1 : a ) {
            if ( b.contains( v1 ) ) {
                return true;
            }
        }
        return false;
    }

	private Graph extractGraph(Particle p) {
		Graph newGraph = new Graph();
		Node start = startNode.clone();
		Node end   = endNode.clone();

		Set<String> currentEndInputs = new HashSet<String>();
		Map<String,Edge> connections = new HashMap<String,Edge>();

		// Connect start node
		connectCandidateToGraphByInputs(start, connections, newGraph, currentEndInputs);

		Set<Node> seenNodes = new HashSet<Node>();
		List<ListItem> candidateList = new ArrayList<ListItem>();

		addToCandidateListFromMasterEdges(start, seenNodes, candidateList, p.dimensions);
		Collections.sort(candidateList);

		// While end cannot be connected to graph
		while(!currentEndInputs.containsAll(end.getInputs())) {

			// Select node
			int index;

			candidateLoop:
			for (index = 0; index < candidateList.size(); index++) {
				ListItem candidateItem = candidateList.get(index);
				Node candidate = candidateItem.node.clone();
				// For all of the candidate inputs, check that there is a service already in the graph
				// that can satisfy it
				connections.clear();

				for (String input : candidate.getInputs()) {
					boolean found = false;
					 for (Node s : taxonomyMap.get(input).servicesWithOutput) {
						 if (newGraph.nodeMap.containsKey(s.getName())) {
							 Set<String> intersect = new HashSet<String>();
							 intersect.add(input);

							 Edge mapEdge = connections.get(s.getName());
							 if (mapEdge == null) {
								 Edge e = new Edge(intersect);
								 e.setFromNode(newGraph.nodeMap.get(s.getName()));
								 e.setToNode(candidate);
								 connections.put(e.getFromNode().getName(), e);
							 }
							 else
								 mapEdge.getIntersect().addAll(intersect);

							 found = true;
							 break;
						 }
					 }
					 // If that input cannot be satisfied, move on to another candidate node to connect
					 if (!found) {
						 // Move on to another candidate
						 continue candidateLoop;
					 }
				}

				// Connect candidate to graph, adding its reachable services to the candidate list
				connectCandidateToGraphByInputs(candidate, connections, newGraph, currentEndInputs);
				addToCandidateListFromMasterEdges(candidate, seenNodes, candidateList, p.dimensions);

				break;
			}

			candidateList.remove(index);
			Collections.sort(candidateList);
		}

		// Connect end node to graph
		connections.clear();
		Iterator<Node> it = newGraph.nodeMap.values().iterator();
		Node s;

		while (!currentEndInputs.isEmpty() && it.hasNext()) {
			s = it.next();

			Set<String> intersection = new HashSet<String>();

			for (String o : s.getOutputs()) {

				Set<String> endNodeInputs = taxonomyMap.get(o).endNodeInputs;
				if (!endNodeInputs.isEmpty()) {

					for (String i : endNodeInputs) {
						if (currentEndInputs.contains(i)) {
							intersection.add(i);
							currentEndInputs.remove(i);
						}
					}

				}
			}

			if (!intersection.isEmpty()) {
				Edge e = new Edge(intersection);
				e.setFromNode(s);
				e.setToNode(end);
				connections.put(e.getFromNode().getName(), e);
			}
		}
		connectCandidateToGraphByInputs(end, connections, newGraph, currentEndInputs);
		removeDanglingNodes(newGraph);
		return newGraph;
	}

	public void removeDanglingNodes(Graph graph) {
	    List<Node> dangling = new ArrayList<Node>();
	    for (Node g : graph.nodeMap.values()) {
	        if (!g.getName().equals("end") && g.getOutgoingEdgeList().isEmpty())
	            dangling.add( g );
	    }

	    for (Node d: dangling) {
	        removeDangling(d, graph);
	    }
	}

	private void removeDangling(Node n, Graph graph) {
	    if (n.getOutgoingEdgeList().isEmpty()) {
	        graph.nodeMap.remove( n.getName() );
	        for (Edge e : n.getIncomingEdgeList()) {
	            e.getFromNode().getOutgoingEdgeList().remove( e );
	            graph.edgeList.remove( e );
	            removeDangling(e.getFromNode(), graph);
	        }
	    }
	}

	private void addToCandidateList(Node n, Set<Node> seenNode, Set<Node> relevant, List<Node> candidateList) {
		seenNode.add(n);
		List<TaxonomyNode> taxonomyOutputs;
		if (n.getName().equals("start"))
			taxonomyOutputs = startNode.getTaxonomyOutputs();
		else
			taxonomyOutputs = serviceMap.get(n.getName()).getTaxonomyOutputs();

		for (TaxonomyNode t : taxonomyOutputs) {
			// Add servicesWithInput from taxonomy node as potential candidates to be connected
			for (Node current : t.servicesWithInput) {
				if (!seenNode.contains(current) && relevant.contains(current)) {
					candidateList.add(current);
					seenNode.add(current);
				}
			}
		}
	}

	private void addToCandidateListFromMasterEdges (Node n, Set<Node> seenNode, List<ListItem> candidateList, float[] dimensions) {
		seenNode.add(n);

		Node original = masterGraph.nodeMap.get(n.getName());

		for (Edge e : original.getOutgoingEdgeList()) {
			// Add servicesWithInput from taxonomy node as potential candidates to be connected
			Node current = e.getToNode();
			if (!seenNode.contains(current)) {
				int index = serviceToIndexMapping.get(current.getName());
				candidateList.add(new ListItem(current, dimensions[index]));
				seenNode.add(current);
			}
		}
	}

	public void connectCandidateToGraphByInputs(Node candidate, Map<String,Edge> connections, Graph graph, Set<String> currentEndInputs) {
		graph.nodeMap.put(candidate.getName(), candidate);
		graph.edgeList.addAll(connections.values());
		candidate.getIncomingEdgeList().addAll(connections.values());

		for (Edge e : connections.values()) {
			Node fromNode = graph.nodeMap.get(e.getFromNode().getName());
			fromNode.getOutgoingEdgeList().add(e);
		}
		for (String o : candidate.getOutputs()) {
			currentEndInputs.addAll(taxonomyMap.get(o).endNodeInputs);
		}
		structureValidator( graph );
	}

	//==========================================================================================================================
	//                                                 Debugging Routines
	//==========================================================================================================================

    private void structureValidator( Graph graph ) {
        for ( Edge e : graph.edgeList ) {
            Node fromNode = graph.nodeMap.get(e.getFromNode().getName());

            boolean isContained = false;
            List<Edge> outgoingEdgeList = fromNode.getOutgoingEdgeList();
            for ( Edge outEdge : outgoingEdgeList ) {

                if ( e == outEdge ) {
                    isContained = true;
                    break;
                }
            }

            if ( !isContained ) {
                System.out.println( "Outgoing edge for node " + fromNode.getName() + " not detected." );
            }

            Node toNode = graph.nodeMap.get(e.getToNode().getName());

            isContained = false;
            List<Edge> incomingEdgeList = toNode.getIncomingEdgeList();
            for ( Edge inEdge : incomingEdgeList ) {
                if ( e == inEdge ) {
                    isContained = true;
                    break;
                }
            }

            if ( !isContained ) {
                System.out.println( "Incoming edge for node " + toNode.getName() + " not detected." );
            }
        }
    }

}
