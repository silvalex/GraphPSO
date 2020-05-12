package pso;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

/**
 * Represents a service node to be added to a composition workflow.
 * The node contains list of incoming and outgoing edges, a name,
 * sets of inputs and outputs, and mappings to taxonomy information.
 * 
 * @author sawczualex
 *
 */
public class Node implements Cloneable {
	private List<Edge> incomingEdgeList = new ArrayList<Edge>();
	private List<Edge> outgoingEdgeList = new ArrayList<Edge>();
	private List<TaxonomyNode> taxonomyOutputs = new ArrayList<TaxonomyNode>();
	private String name;
	private double[] qos;
	private Set<String> inputs;
	private Set<String> outputs;

	/**
	 * Creates a new service node with the given information.
	 * 
	 * @param name - Service name
	 * @param qos - Quality of Service (QoS) attributes for service
	 * @param inputs - Set of service inputs
	 * @param outputs - Set of service outputs
	 */
	public Node(String name, double[] qos, Set<String> inputs, Set<String> outputs) {
		this.name = name;
		this.qos = qos;
		this.inputs = inputs;
		this.outputs = outputs;
	}

	/**
	 * Returns list of incoming edges to this node.
	 * 
	 * @return incoming edge list
	 */
	public List<Edge> getIncomingEdgeList() {
		return incomingEdgeList;
	}

	/**
	 * Returns list of outgoing edges from this node.
	 * 
	 * @return outgoing edge list
	 */
	public List<Edge> getOutgoingEdgeList() {
		return outgoingEdgeList;
	}

	/**
	 * Gets array of QoS attributes from this service node.
	 * Constants with the order (i.e. corresponding index)
	 * of QoS attributes are defined in the GraphPSO class.
	 * 
	 * @return QoS attributes
	 */
	public double[] getQos() {
		return qos;
	}

	/**
	 * Gets set of inputs for this service node.
	 * 
	 * @return input set
	 */
	public Set<String> getInputs() {
		return inputs;
	}

	/**
	 * Gets set of outputs for this service node.
	 * 
	 * @return output set
	 */
	public Set<String> getOutputs() {
		return outputs;
	}

	/**
	 * Retuns this service's name.
	 * 
	 * @return service name
	 */
	public String getName() {
		return name;
	}

	/**
	 * Set this service's name.
	 * 
	 * @param name - The new name
	 */
	public void setName(String name) {
		this.name = name;
	}

	/**
	 * Shallow-clones this node, creating a new instance
	 * that shares QoS, inputs, and outputs data structures.
	 * 
	 * @return cloned node
	 */
	public Node clone() {
		return new Node(name, qos, inputs, outputs);
	}

	/**
	 * Gets the list of taxonomy nodes associated with
	 * this service.
	 * 
	 * @return taxonomy node list
	 */
	public List<TaxonomyNode> getTaxonomyOutputs() {
		return taxonomyOutputs;
	}

	@Override
	/**
	 * {@inheritDoc}
	 */
	public String toString(){
		return name;
	}

	@Override
	/**
	 * {@inheritDoc}
	 */
	public int hashCode() {
		return name.hashCode();
	}

	@Override
	/**
	 * {@inheritDoc}
	 */
	public boolean equals(Object other) {
		if (other instanceof Node) {
			Node o = (Node) other;
			return name.equals(o.name);
		}
		else
			return false;
	}
}
