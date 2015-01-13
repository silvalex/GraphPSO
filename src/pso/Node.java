package pso;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;


public class Node implements Cloneable {
	private List<Edge> incomingEdgeList = new ArrayList<Edge>();
	private List<Edge> outgoingEdgeList = new ArrayList<Edge>();
	private List<TaxonomyNode> taxonomyOutputs = new ArrayList<TaxonomyNode>();
	private String name;
	private double[] qos;
	private Set<String> inputs;
	private Set<String> outputs;

	public Node(String name, double[] qos, Set<String> inputs, Set<String> outputs) {
		this.name = name;
		this.qos = qos;
		this.inputs = inputs;
		this.outputs = outputs;
	}

	public List<Edge> getIncomingEdgeList() {
		return incomingEdgeList;
	}

	public List<Edge> getOutgoingEdgeList() {
		return outgoingEdgeList;
	}

	public double[] getQos() {
		return qos;
	}

	public Set<String> getInputs() {
		return inputs;
	}

	public Set<String> getOutputs() {
		return outputs;
	}

	public String getName() {
		return name;
	}

	public Node clone() {
		return new Node(name, qos, inputs, outputs);
	}

	public List<TaxonomyNode> getTaxonomyOutputs() {
		return taxonomyOutputs;
	}

	@Override
	public String toString(){
		return name;
	}

	@Override
	public int hashCode() {
		return name.hashCode();
	}

	@Override
	public boolean equals(Object other) {
		if (other instanceof Node) {
			Node o = (Node) other;
			return name.equals(o.name);
		}
		else
			return false;
	}
}
