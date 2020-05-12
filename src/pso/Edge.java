package pso;

import java.util.Set;

/**
 * Represents a directed edge that connects two nodes in a service workflow. This class
 * records the intersect between input and output values for the origin and destination
 * edges.
 * 
 * @author sawczualex
 */
public class Edge {
	private Node fromNode;
	private Node toNode;
	private Set<String> intersect;

	/**
	 * Creates a new edge without setting origin and destination nodes.
	 * 
	 * @param intersect - The output values from the origin node that are
	 * consumed by the inputs of the destination node.
	 */
	public Edge(Set<String> intersect) {
		this.intersect = intersect;
	}

	/**
	 * Gets the origin node for this edge.
	 * 
	 * @return origin node
	 */
	public Node getFromNode() {
		return fromNode;
	}

	/**
	 * Gets the destination node for this edge.
	 * 
	 * @return destination node
	 */
	public Node getToNode() {
		return toNode;
	}

	/**
	 * Sets the origin node for this edge.
	 * 
	 * @param fromNode - The origin node
	 */
	public void setFromNode(Node fromNode) {
		this.fromNode = fromNode;
	}

	/**
	 * Sets the destination node for this edge.
	 * 
	 * @param toNode - The destination node
	 */
	public void setToNode(Node toNode) {
		this.toNode = toNode;
	}

	/**
	 * Returns the intersect between the output values of
	 * the origin node and the input values of the destination
	 * node, i.e. the common values connecting the two nodes.
	 * 
	 * @return intersect
	 */
	public Set<String> getIntersect() {
		return intersect;
	}

	@Override
	/**
	 * {@inheritDoc}
	 */
	public String toString() {
		return String.format("%s->%s", fromNode, toNode);
	}

	@Override
	/**
	 * {@inheritDoc}
	 */
	public int hashCode() {
		return (fromNode.getName() + toNode.getName()).hashCode();
	}

	@Override
	/**
	 * {@inheritDoc}
	 */
	public boolean equals(Object other) {
		if (other instanceof Edge) {
			Edge o = (Edge) other;
			return fromNode.getName().equals(o.fromNode.getName()) && toNode.getName().equals(o.toNode.getName());
		}
		else
			return false;
	}
}
