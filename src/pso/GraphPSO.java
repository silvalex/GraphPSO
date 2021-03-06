package pso;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
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

/**
 * Optimises solutions for the Web service composition problem using
 * Particle Swarm Optimisation (PSO) coupled with a graph-based solution
 * representation. This class contains the application's entry point,
 * as well as methods for carrying out the optimisation process.
 * 
 * @author sawczualex
 */
public class GraphPSO {

	// PSO settings
	public List<Particle> swarm = new ArrayList<Particle>();
	public static final int MAX_NUM_ITERATIONS = 100;
	public static final int NUM_PARTICLES = 30;
	public static final float C1 = 1.49618f;
	public static final float C2 = 1.49618f;
	public static final float W = 0.7298f;
	public static final boolean dynamicNormalisation = true;
	public static int numDimensions;
	public static ArrayList<Long> initTime = new ArrayList<Long>();
	public static ArrayList<Long> time = new ArrayList<Long>();
	public static ArrayList<Double> meanFitness = new ArrayList<Double>();
	public static ArrayList<Double> bestFitnessThisGen = new ArrayList<Double>();
	public static ArrayList<Double> bestFitnessSoFar = new ArrayList<Double>();
	public static String logName;
	public static String histogramLogName;
	public static Long initialisationStartTime;

	// Fitness function weights
	public static final double W1 = 0.25;
	public static final double W2 = 0.25;
	public static final double W3 = 0.25;
	public static final double W4 = 0.25;

	public static double MINIMUM_COST = Double.MAX_VALUE;
	public static double MINIMUM_TIME = Double.MAX_VALUE;
	public static double MINIMUM_RELIABILITY = 0;
	public static double MINIMUM_AVAILABILITY = 0;
	public static double MAXIMUM_COST = Double.MIN_VALUE;
	public static double MAXIMUM_TIME = Double.MIN_VALUE;
	public static double MAXIMUM_RELIABILITY = Double.MIN_VALUE;
	public static double MAXIMUM_AVAILABILITY = Double.MIN_VALUE;

	public static List<Double> meanAvailPerGen = new ArrayList<Double>();
	public static List<Double> meanReliaPerGen = new ArrayList<Double>();
	public static List<Double> meanTimePerGen = new ArrayList<Double>();
	public static List<Double> meanCostPerGen = new ArrayList<Double>();
	public static List<Double> bestAvailThisGen = new ArrayList<Double>();
	public static List<Double> bestReliaThisGen = new ArrayList<Double>();
	public static List<Double> bestTimeThisGen = new ArrayList<Double>();
	public static List<Double> bestCostThisGen = new ArrayList<Double>();
	public static List<Double> bestAvailSoFar = new ArrayList<Double>();
	public static List<Double> bestReliaSoFar = new ArrayList<Double>();
	public static List<Double> bestTimeSoFar = new ArrayList<Double>();
	public static List<Double> bestCostSoFar = new ArrayList<Double>();

	// Constants with the order of QoS attributes
	public static final int TIME = 0;
	public static final int COST = 1;
	public static final int AVAILABILITY = 2;
	public static final int RELIABILITY = 3;

	public Map<String, Node> serviceMap = new HashMap<String, Node>();
	public Set<Node> relevant;
	public Map<String, Integer> serviceToIndexMap = new HashMap<String, Integer>();
	public Map<String, TaxonomyNode> taxonomyMap = new HashMap<String, TaxonomyNode>();
	public Set<String> taskInput;
	public Set<String> taskOutput;
	public Node startNode;
	public Node endNode;
	private Random random;

	// Statistics tracking
	Map<String, Integer> nodeCount = new HashMap<String, Integer>();
	Map<String, Integer> edgeCount = new HashMap<String, Integer>();

	/**
	 * Application's entry point.
	 *
	 * @param args
	 */
	public static void main(String[] args) {
		new GraphPSO(args[0], args[1], args[2], args[3], args[4], Long.valueOf(args[5]));
	}

	/**
	 * Creates a functionally correct workflow, and runs PSO to discover the most optimal services
	 * to be used in it.
	 * 
	 * @param lName - Name of log file to write to
	 * @param hlName - Name of histogram file to write to
	 * @param taskFileName - Name of file containing Web service composition task/request 
	 * @param serviceFileName - Name of file containing information on Web services
	 * @param taxonomyFileName - Name of file containing information on relationship between
	 * output and input types
	 * @param seed - Seed for random number generator, used for reproducing a specific run
	 */
	public GraphPSO(String lName, String hlName, String taskFileName, String serviceFileName, String taxonomyFileName, long seed) {
		initialisationStartTime = System.currentTimeMillis();

		logName = lName;
		histogramLogName = hlName;
		random = new Random(seed);

		// Parse files and figure out taxonomic relationship between inputs and outputs
		parseWSCServiceFile(serviceFileName);
		parseWSCTaskFile(taskFileName);
		parseWSCTaxonomyFile(taxonomyFileName);
		findConceptsForInstances();

		// Create mock nodes to act as the start and end points of the service workflow
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
		// relevant = getRelevantServices(serviceMap, taskInput, taskOutput); XXX
		relevant = new HashSet<Node>(serviceMap.values());
		numDimensions = relevant.size();

		mapServicesToIndices(relevant,serviceToIndexMap);
		if (!dynamicNormalisation)
			calculateNormalisationBounds(relevant);

		String finalGraph = runPSO();
		writeLogs(finalGraph);
	}

	//==========================================================================================================
	//
    //                                              PSO METHODS
    //
	//==========================================================================================================

	/**
	 * Conducts the particle swarm optimisation.
	 * 
	 * @return String representation of the service graph corresponding
	 * to the globally best solution encountered
	 */
	public String runPSO() {
		// 1. Initialize the swarm
		initializeRandomSwarm();

		int i = 0;
		Particle p;
		Graph workflow = null;
		long initialization = System.currentTimeMillis() - initialisationStartTime;

		while (i < MAX_NUM_ITERATIONS) {
			long startTime = System.currentTimeMillis();
			System.out.println("ITERATION " + i);

			double minAvailability = 2.0;
			double maxAvailability = -1.0;
			double minReliability = 2.0;
			double maxReliability = -1.0;
			double minTime = Double.MAX_VALUE;
			double maxTime = -1.0;
			double minCost = Double.MAX_VALUE;
			double maxCost = -1.0;

			// Keep track of means
			double meanAvailability = 0.0;
			double meanReliability = 0.0;
			double meanTime = 0.0;
			double meanCost = 0.0;
			double meanFit = 0.0;
			double bestFitThisGen = 0.0;
			double bestAThisGen = 0.0;
			double bestRThisGen = 0.0;
			double bestTThisGen = Double.MAX_VALUE;
			double bestCThisGen = Double.MAX_VALUE;

			// Go through all particles
			for (int j = 0; j < NUM_PARTICLES; j++) {
				System.out.println("\tPARTICLE " + j);
				p = swarm.get(j);
				workflow = createNewGraph(startNode.clone(), endNode.clone(), relevant, p.dimensions);

				// 2. Evaluate fitness of particle
				calculateOverallQoS(workflow, p);

				meanAvailability += p.availability;
				meanReliability += p.reliability;
				meanTime += p.time;
				meanCost += p.cost;

				if (dynamicNormalisation) {
					if (p.availability < minAvailability)
						minAvailability = p.availability;
					if (p.availability > maxAvailability)
						maxAvailability = p.availability;
					if (p.reliability < minReliability)
						minReliability = p.reliability;
					if (p.reliability > maxReliability)
						maxReliability = p.reliability;
					if (p.time < minTime)
						minTime = p.time;
					if (p.time > maxTime)
						maxTime = p.time;
					if (p.cost < minCost)
						minCost = p.cost;
					if (p.cost > maxCost)
						maxCost = p.cost;
				}

				if (!dynamicNormalisation) {
					double fit = calculateFitness(p);
					meanFit += fit;
					if (fit > bestFitThisGen) {
						bestFitThisGen = fit;
						bestAThisGen = p.availability;
						bestRThisGen = p.reliability;
						bestTThisGen = p.time;
						bestCThisGen = p.cost;
					}
					// 3. If fitness of particle is better than Pbest, update the Pbest
					p.updatePersonalBest();
					// 4. If fitness of Pbest is better than Gbest, update the Gbest
					if (p.bestFitness > Particle.globalBestFitness) {
						Particle.globalBestFitness = p.bestFitness;
						Particle.globalGraphString = p.graphString;
						Particle.globalBestDimensions = Arrays.copyOf(p.bestDimensions, p.bestDimensions.length);
						Particle.globalBestAvailability = p.availability;
						Particle.globalBestReliability = p.reliability;
						Particle.globalBestTime = p.time;
						Particle.globalBestCost = p.cost;
					}
					// 5. Update the velocity of particle
					updateVelocity(p);
					// 6. Update the position of particle
					updatePosition(p);
				}
			}

			// Mean QoS
			meanAvailPerGen.add(meanAvailability / NUM_PARTICLES);
			meanReliaPerGen.add(meanReliability / NUM_PARTICLES);
			meanTimePerGen.add(meanTime / NUM_PARTICLES);
			meanCostPerGen.add(meanCost / NUM_PARTICLES);

			// If normalisation is dynamic, go through particles again to finish fitness calculations
			if (dynamicNormalisation) {
				// Update the normalisation bounds with the newly found values
				MINIMUM_AVAILABILITY = minAvailability;
				MAXIMUM_AVAILABILITY = maxAvailability;
				MINIMUM_RELIABILITY = minReliability;
				MAXIMUM_RELIABILITY = maxReliability;
				MINIMUM_COST = minCost;
				MAXIMUM_COST = maxCost;
				MINIMUM_TIME = minTime;
				MAXIMUM_TIME = maxTime;

				Particle p2;
				// Finish calculating the fitness of each candidate
				for (int j = 0; j < NUM_PARTICLES; j++) {
					p2 = swarm.get(j);
					double fit = calculateFitness(p2);
					meanFit += fit;
					if (fit > bestFitThisGen) {
						bestFitThisGen = fit;
						bestAThisGen = p2.availability;
						bestRThisGen = p2.reliability;
						bestTThisGen = p2.time;
						bestCThisGen = p2.cost;
					}
					// 3. If fitness of particle is better than Pbest, update the Pbest
					p2.updatePersonalBest();
					// 4. If fitness of Pbest is better than Gbest, update the Gbest
					if (p2.bestFitness > Particle.globalBestFitness) {
						Particle.globalBestFitness = p2.bestFitness;
						Particle.globalGraphString = p2.graphString;
						Particle.globalBestDimensions = Arrays.copyOf(p2.bestDimensions, p2.bestDimensions.length);
						Particle.globalBestAvailability = p2.availability;
						Particle.globalBestReliability = p2.reliability;
						Particle.globalBestTime = p2.time;
						Particle.globalBestCost = p2.cost;
					}
					// 5. Update the velocity of particle
					updateVelocity(p2);
					// 6. Update the position of particle
					updatePosition(p2);
				}
			}

			meanFitness.add(meanFit/NUM_PARTICLES);
			bestFitnessThisGen.add(bestFitThisGen);
			bestAvailThisGen.add(bestAThisGen);
			bestReliaThisGen.add(bestRThisGen);
			bestTimeThisGen.add(bestTThisGen);
			bestCostThisGen.add(bestCThisGen);
			bestFitnessSoFar.add(Particle.globalBestFitness);
			bestAvailSoFar.add(Particle.globalBestAvailability);
			bestReliaSoFar.add(Particle.globalBestReliability);
			bestTimeSoFar.add(Particle.globalBestTime);
			bestCostSoFar.add(Particle.globalBestCost);

			initTime.add(initialization);
			time.add(System.currentTimeMillis() - startTime);
			initialization = 0;
			i++;
		}

		return Particle.globalGraphString;
	}

	/**
	 * Updates the velocity vector of a particle.
	 *
	 * @param p - Particle whose vector is to be updated
	 */
	public void updateVelocity(Particle p) {
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
	public void initializeRandomSwarm() {
		swarm.clear();
		for (int i = 0; i < NUM_PARTICLES; i++) {
			swarm.add(new Particle(random));
		}
	}

	/**
	 * Updates the position (i.e. dimension vector) of a particle.
	 *
	 * @param p - Particle whose position is to be updated
	 */
	public void updatePosition(Particle p) {
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
	 * Given a set of services deemed as relevant, this method assigns a unique integer value
	 * for each service. This information is used during the decoding process.
	 * 
	 * @param relevant - Set of relevant services
	 * @param serviceToIndexMap - Map of service name to index value
	 */
	private void mapServicesToIndices(Set<Node> relevant, Map<String,Integer> serviceToIndexMap) {
		int i = 0;
		for (Node r : relevant) {
			serviceToIndexMap.put(r.getName(), i++);
		}
	}

	//==========================================================================================================
	//
    //                                              PARSING METHODS
    //
	//==========================================================================================================

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
            e.printStackTrace();
		}
		catch (SAXException e) {
            System.out.println("Task file parsing failed...");
            e.printStackTrace();
		}
		catch (IOException e) {
            System.out.println("Task file parsing failed...");
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
	    	NodeList taxonomyRoots = doc.getChildNodes();

	    	processTaxonomyChildren(null, taxonomyRoots);
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
	 * Helper function that recursively recreates the taxonomy structure from a file,
	 * connecting children taxonomy nodes to their parents.
	 *
	 * @param parent - Taxonomy node parent
	 * @param nodes - List of children nodes of this parent 
	 */
	private void processTaxonomyChildren(TaxonomyNode parent, NodeList nodes) {
		if (nodes != null && nodes.getLength() != 0) {
			for (int i = 0; i < nodes.getLength(); i++) {
				org.w3c.dom.Node ch = nodes.item(i);

			if (!(ch instanceof Text)) {
				Element currNode = (Element) nodes.item(i);
				String value = currNode.getAttribute("name");
					TaxonomyNode taxNode = taxonomyMap.get( value );
					if (taxNode == null) {
					    taxNode = new TaxonomyNode(value);
					    taxonomyMap.put( value, taxNode );
					}
					if (parent != null) {
					    taxNode.parents.add(parent);
						parent.children.add(taxNode);
					}

					NodeList children = currNode.getChildNodes();
					processTaxonomyChildren(taxNode, children);
				}
			}
		}
	}

	//==========================================================================================================
	//
    //                                              TAXONOMY METHODS
    //
	//==========================================================================================================

	/**
	 * Populates the taxonomy tree by associating services to the
	 * nodes in the tree.
	 */
	private void populateTaxonomyTree() {
		for (Node s: serviceMap.values()) {
			addServiceToTaxonomyTree(s);
		}
	}

	/**
	 * Add references to this service node throughout the taxonomy tree,
	 * based on the node's inputs and output values.
	 * 
	 * @param s - Service node to be added to taxonomy tree
	 */
	private void addServiceToTaxonomyTree(Node s) {
		// Populate outputs
	    Set<TaxonomyNode> seenConceptsOutput = new HashSet<TaxonomyNode>();
		for (String outputVal : s.getOutputs()) {
			TaxonomyNode n = taxonomyMap.get(outputVal);
			s.getTaxonomyOutputs().add(n);

			// Also add output to all parent nodes
			Queue<TaxonomyNode> queue = new LinkedList<TaxonomyNode>();
			queue.add( n );

			while (!queue.isEmpty()) {
			    TaxonomyNode current = queue.poll();
		        seenConceptsOutput.add( current );
		        current.servicesWithOutput.add(s);
		        for (TaxonomyNode parent : current.parents) {
		            if (!seenConceptsOutput.contains( parent )) {
		                queue.add(parent);
		                seenConceptsOutput.add(parent);
		            }
		        }
			}
		}
		// Populate inputs
		Set<TaxonomyNode> seenConceptsInput = new HashSet<TaxonomyNode>();
		for (String inputVal : s.getInputs()) {
			TaxonomyNode n = taxonomyMap.get(inputVal);

			// Also add input to all children nodes
			Queue<TaxonomyNode> queue = new LinkedList<TaxonomyNode>();
			queue.add( n );

			while(!queue.isEmpty()) {
				TaxonomyNode current = queue.poll();
				seenConceptsInput.add( current );

			    Set<String> inputs = current.servicesWithInput.get(s);
			    if (inputs == null) {
			    	inputs = new HashSet<String>();
			    	inputs.add(inputVal);
			    	current.servicesWithInput.put(s, inputs);
			    }
			    else {
			    	inputs.add(inputVal);
			    }

			    for (TaxonomyNode child : current.children) {
			        if (!seenConceptsInput.contains( child )) {
			            queue.add(child);
			            seenConceptsInput.add( child );
			        }
			    }
			}
		}
		return;
	}

	/**
	 * Converts input, output, and service instance values to their corresponding
	 * ontological parent. This method should be called only once, before the
	 * optimisation starts.
	 */
	private void findConceptsForInstances() {
		Set<String> temp = new HashSet<String>();

		for (String s : taskInput)
			temp.add(taxonomyMap.get(s).parents.get(0).value);
		taskInput.clear();
		taskInput.addAll(temp);

		temp.clear();
		for (String s : taskOutput)
				temp.add(taxonomyMap.get(s).parents.get(0).value);
		taskOutput.clear();
		taskOutput.addAll(temp);

		for (Node s : serviceMap.values()) {
			temp.clear();
			Set<String> inputs = s.getInputs();
			for (String i : inputs)
				temp.add(taxonomyMap.get(i).parents.get(0).value);
			inputs.clear();
			inputs.addAll(temp);

			temp.clear();
			Set<String> outputs = s.getOutputs();
			for (String o : outputs)
				temp.add(taxonomyMap.get(o).parents.get(0).value);
			outputs.clear();
			outputs.addAll(temp);
		}
	}

	//==========================================================================================================
	//
    //                                              GRAPH METHODS
    //
	//==========================================================================================================

	/**
	 * Calculates the overall Quality of Service (QoS) attributes for a given
	 * candidate. These QoS values are stored within the particle provided as
	 * an argument.
	 * 
	 * @param graph - The candidate's corresponding service workflow
	 * @param p - The candidate particle
	 */
	public void calculateOverallQoS(Graph graph, Particle p) {

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

        p.availability = a;
        p.reliability = r;
        p.time = t;
        p.cost = c;
        p.graphString = graph.toString();
	}

	/**
	 * Calculates the fitness score of a particle based on its
	 * Quality of Service (QoS) attributes. The fitness consists
	 * of a weighted sum of normalised QoS attributes, where the
	 * weights are set as constants in the GraphPSO class.
	 * 
	 * @param p - The given particle
	 * @return fitness score
	 */
	public double calculateFitness(Particle p) {
		double a = p.availability;
        double r = p.reliability;
        double t = p.time;
        double c = p.cost;

        a = normaliseAvailability(a);
        r = normaliseReliability(r);
        t = normaliseTime(t);
        c = normaliseCost(c);

        p.fitness = W1 * a + W2 * r + W3 * t + W4 * c;
        return p.fitness;
	}

	/**
	 * Normalises availability score, using minimum and maximum availability
	 * values defined as constants in the GraphPSO class. If the minimum and
	 * maximum availability are the same, 1.0 is returned. The higher the
	 * normalised value, the better.
	 * 
	 * @param availability - The score to be normalised
	 * @return normalised availability value
	 */
	public double normaliseAvailability(double availability) {
		if (MAXIMUM_AVAILABILITY - MINIMUM_AVAILABILITY == 0.0)
			return 1.0;
		else
			return (availability - MINIMUM_AVAILABILITY)/(MAXIMUM_AVAILABILITY - MINIMUM_AVAILABILITY);
	}

	/**
	 * Normalises the reliability score, using minimum and maximum reliability
	 * values defined as constants in the GraphPSO class. If the minimum and
	 * maximum reliability are the same, 1.0 is returned. The higher the
	 * normalised value, the better.
	 * 
	 * @param reliability - The score to be normalised
	 * @return normalised reliability score
	 */
	public double normaliseReliability(double reliability) {
		if (MAXIMUM_RELIABILITY - MINIMUM_RELIABILITY == 0.0)
			return 1.0;
		else
			return (reliability - MINIMUM_RELIABILITY)/(MAXIMUM_RELIABILITY - MINIMUM_RELIABILITY);
	}

	/**
	 * Normalises the time score, using minimum and maximum time values defined
	 * as constants in the GraphPSO class. If the minimum and maximum time are
	 * the same, 1.0 is returned. As with availability and reliability, we want
	 * the highest normalised value to be the better. Thus, the raw time value
	 * (the lower, the better) is offset when calculating the normalised score
	 * (the higher, the better).
	 * 
	 * @param time - The score to be normalised
	 * @return normalised time score
	 */
	public double normaliseTime(double time) {
		if (MAXIMUM_TIME - MINIMUM_TIME == 0.0)
			return 1.0;
		else
			return (MAXIMUM_TIME - time)/(MAXIMUM_TIME - MINIMUM_TIME);
	}

	/**
	 * Normalises the cost store, using minimum and maximum cost values defined
	 * as constants in the GraphPSO class. If the minimum and maximum time are
	 * the same, 1.0 is returned. As with availability and reliability, we want
	 * the highest normalised value to be better. Thus, the raw cost value (the
	 * lower, the better) is offset when calculating the normalised score (the
	 * higher, the better).
	 * 
	 * @param cost - The score to be normalised
	 * @return normalised cost score
	 */
	public double normaliseCost(double cost) {
		if (MAXIMUM_COST - MINIMUM_COST == 0.0)
			return 1.0;
		else
			return (MAXIMUM_COST - cost)/(MAXIMUM_COST - MINIMUM_COST);
	}

	/**
	 * Uses the Bellman-Ford algorithm with negative weights to find the longest
	 * path in an directed acyclic graph.
	 *
	 * @param g - The directed acyclic graph
	 * @return list of edges composing the longest path
	 */
	public double findLongestPath(Graph g) {
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
				if ((distance.get(e.getFromNode().getName()) -
				        e.getToNode().getQos()[TIME])
				        < distance.get(e.getToNode().getName())) {
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
	 * Decodes a new functionally correct composition graph according to the
	 * particle weights provided.
	 * 
	 * @param start - Start node for the new graph
	 * @param end - End node for the new graph
	 * @param relevant - List of relevant service nodes, i.e. those considered for this composition task
	 * @param weights - Weights from a given particle
	 * @return
	 */
	public Graph createNewGraph(Node start, Node end, Set<Node> relevant, float[] weights) {

		Graph newGraph = new Graph();

		Set<String> currentEndInputs = new HashSet<String>();
		Map<String,Edge> connections = new HashMap<String,Edge>();

		// Connect start node
		connectCandidateToGraphByInputs(start, connections, newGraph, currentEndInputs);

		List<ListItem> candidateList = new ArrayList<ListItem>();

		populateCandidateList(serviceToIndexMap, relevant, candidateList, weights);
		Collections.sort(candidateList);

		finishConstructingGraph(currentEndInputs, end, candidateList, connections, newGraph, relevant);

        // Keep track of nodes and edges for statistics
        for (String nodeName : newGraph.nodeMap.keySet())
            addToCountMap(nodeCount, nodeName);
        for (Edge edge : newGraph.edgeList)
            addToCountMap(edgeCount, edge.toString());
		return newGraph;
	}

	/**
	 * Helper method to increment the count value for a string key
	 * in a given map. If string is not already a key in map, it will
	 * be added wih a count of 1.
	 *  
	 * @param map - Map to be updated
	 * @param item - String item
	 */
   public void addToCountMap(Map<String,Integer> map, String item) {
        if (map.containsKey( item )) {
            map.put( item, map.get( item ) + 1 );
        }
        else {
            map.put( item, 1 );
        }
    }


   /**
    * Adds relevant service nodes to the candidate list, which is provided as a parameter.
    * This method associates each candidate to its corresponding particle weight.
    * 
    * @param serviceToIndexMap - Map that allows a service's index to be retrieved based
    * on its name
    * @param relevant - Set of relevant services, i.e. those considered for this composition
    * task
    * @param candidateList - The candidate list to be populated
    * @param weights - The particle weights to associate with each given service
    */
	public void populateCandidateList(Map<String, Integer> serviceToIndexMap, Set<Node> relevant, List<ListItem> candidateList, float[] weights) {
		// Go through all relevant nodes
		for (Node n : relevant) {
			// Find the index for that node
			int index = serviceToIndexMap.get(n.getName());
			// Retrieve weight associated with service using the index, and create list item
			ListItem item = new ListItem(n.getName(), weights[index]);
			// Add item to list
			candidateList.add(item);
		}
	}

	/**
	 * Helper method to finish constructing the decoded graph. It iteratively
	 * adds services to the graph until the end node has all of its inputs fulfilled.
	 * 
	 * @param currentEndInputs - Set of end inputs to be fulfilled
	 * @param end - The end node
	 * @param candidateList - The list of candidate services to be added to the graph
	 * @param connections - Map to track which values should be added to outgoing edges
	 * @param newGraph - The decoded graph that is being built
	 * @param relevant - Set of relevant services, i.e. those considered for this 
	 * composition task
	 */
	public void finishConstructingGraph(Set<String> currentEndInputs, Node end, List<ListItem> candidateList, Map<String,Edge> connections,
	        Graph newGraph, Set<Node> relevant) {

	 // While end cannot be connected to graph
		while(!checkCandidateNodeSatisfied(connections, newGraph, end, end.getInputs())){
			connections.clear();

            // Select node
            int index;

            candidateLoop:
            for (index = 0; index < candidateList.size(); index++) {
            	ListItem item = candidateList.get(index);
            	Node candidate = serviceMap.get(item.serviceName).clone();
                // For all of the candidate inputs, check that there is a service already in the graph
                // that can satisfy it

                if (!checkCandidateNodeSatisfied(connections, newGraph, candidate, candidate.getInputs())) {
                    connections.clear();
                	continue candidateLoop;
                }

                // Connect candidate to graph, adding its reachable services to the candidate list
                connectCandidateToGraphByInputs(candidate, connections, newGraph, currentEndInputs);
                connections.clear();

                break;
            }

            candidateList.remove(index);
        }

        connectCandidateToGraphByInputs(end, connections, newGraph, currentEndInputs);
        connections.clear();
        removeDanglingNodes(newGraph);
	}

	/**
	 * Checks whether all of the inputs of a given candidate can be satisfied by the
	 * outputs of nodes already in the graph.
	 * 
	 * @param connections - Map to track which values should be added to outgoing edges
	 * @param newGraph - The decoded graph that is being built
	 * @param candidate - The candidate node being checked
	 * @param candInputs - Set of inputs for the candidate node
	 * @return true if all inputs are satisfied, false otherwise
	 */
	public boolean checkCandidateNodeSatisfied(Map<String, Edge> connections, Graph newGraph,
			Node candidate, Set<String> candInputs) {

		Set<String> candidateInputs = new HashSet<String>(candInputs);
		Set<String> startIntersect = new HashSet<String>();

		// Check if the start node should be considered
		Node start = newGraph.nodeMap.get("start");

		for(String output : start.getOutputs()) {
			Set<String> inputVals = taxonomyMap.get(output).servicesWithInput.get(candidate);
			if (inputVals != null) {
				candidateInputs.removeAll(inputVals);
				startIntersect.addAll(inputVals);
			}
		}

		if (!startIntersect.isEmpty()) {
			Edge startEdge = new Edge(startIntersect);
			startEdge.setFromNode(start);
			startEdge.setToNode(candidate);
			connections.put(start.getName(), startEdge);
		}
		

		for (String input : candidateInputs) {
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
					} else
						mapEdge.getIntersect().addAll(intersect);

					found = true;
					break;
				}
			}
			// If that input cannot be satisfied, move on to another candidate
			// node to connect
			if (!found) {
				// Move on to another candidate
				return false;
			}
		}
		return true;
	}

	/**
	 * Connects candidate node to graph, including the appropriate edges.
	 *  
	 * @param candidate - Node to be added to graph
	 * @param connections - Map containing edges to be added
	 * @param graph - The decoded graph being constructed
	 * @param currentEndInputs - Set of inputs for the end node that
	 * have already been satisfied
	 */
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
	}

	/**
	 * Removes the nodes from a decoded graph that do not contribute
	 * (directly or indirectly) to fulfilling the inputs of the end
	 * node, i.e. they are 'dangling'.
	 * 
	 * @param graph - The graph to be cleaned
	 */
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

	/**
	 * Helper recursive method for removing nodes that do not contribute
	 * to fulfilling the inputs of the end node.
	 * 
	 * @param n - Current node being checked for deletion
	 * @param graph - The graph to be cleaned
	 */
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

	//==========================================================================================================
	//
    //                                              AUXILIARY METHODS
    //
	//==========================================================================================================
	
	/**
	 * Estimates Quality of Service (QoS) normalisation bounds based on
	 * a provided set of services.
	 *  
	 * @param services - The set of services
	 */
	private void calculateNormalisationBounds(Set<Node> services) {
		for(Node service: services) {
			double[] qos = service.getQos();

			// Availability
			double availability = qos[AVAILABILITY];
			if (availability > MAXIMUM_AVAILABILITY)
				MAXIMUM_AVAILABILITY = availability;

			// Reliability
			double reliability = qos[RELIABILITY];
			if (reliability > MAXIMUM_RELIABILITY)
				MAXIMUM_RELIABILITY = reliability;

			// Time
			double time = qos[TIME];
			if (time > MAXIMUM_TIME)
				MAXIMUM_TIME = time;
			if (time < MINIMUM_TIME)
				MINIMUM_TIME = time;

			// Cost
			double cost = qos[COST];
			if (cost > MAXIMUM_COST)
				MAXIMUM_COST = cost;
			if (cost < MINIMUM_COST)
				MINIMUM_COST = cost;
		}
		// Adjust max. cost and max. time based on the number of services in shrunk repository
		MAXIMUM_COST *= services.size();
		MAXIMUM_TIME *= services.size();

	}

	//==========================================================================================================
	//
    //                                              LOGGING METHODS
    //
	//==========================================================================================================

	/**
	 * Writes information and statistics of the PSO run to the log and
	 * histogram files.
	 * 
	 * @param finalGraph - A string representation of the final graph
	 * identified as the best solution.
	 */
	public void writeLogs(String finalGraph) {
		try {
			FileWriter writer = new FileWriter(new File(logName));
			for (int i = 0; i < bestFitnessSoFar.size(); i++) {
				writer.append(String.format("%d %d %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",
						i, initTime.get(i), time.get(i), meanFitness.get(i), bestFitnessThisGen.get(i), bestFitnessSoFar.get(i),
						meanAvailPerGen.get(i),   meanReliaPerGen.get(i),    meanTimePerGen.get(i),    meanCostPerGen.get(i),
						bestAvailThisGen.get(i),  bestReliaThisGen.get(i),   bestTimeThisGen.get(i),   bestCostThisGen.get(i),
						bestAvailSoFar.get(i),    bestReliaSoFar.get(i),     bestTimeSoFar.get(i),     bestCostSoFar.get(i)));
			}
			writer.append(finalGraph);
			writer.close();

			FileWriter histogramWriter = new FileWriter(new File(histogramLogName));

			// Write node histogram
			List<String> keyList = new ArrayList<String>(nodeCount.keySet());
			Collections.sort( keyList );

			for (String key : keyList)
			    histogramWriter.append( key + " " );
			histogramWriter.append( "\n" );
			for (String key : keyList)
			    histogramWriter.append( String.format("%d ", nodeCount.get( key )) );
			histogramWriter.append( "\n" );

			// Write edge histogram
	        List<String> edgeList = new ArrayList<String>(edgeCount.keySet());
	        Collections.sort( edgeList );

            for (String key : edgeList)
                histogramWriter.append( key + " " );
            histogramWriter.append( "\n" );
            for (String key : edgeList)
                histogramWriter.append( String.format("%d ", edgeCount.get( key )) );
            histogramWriter.append( "\n" );
            histogramWriter.close();
		}
		catch (IOException e) {
			e.printStackTrace();
		}
	}
}
