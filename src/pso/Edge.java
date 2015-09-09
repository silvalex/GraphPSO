package pso;

import java.util.Set;

public class Edge {
	private Node fromNode;
	private Node toNode;
	private Set<String> intersect;
	private boolean consider = true;

	public Edge(Set<String> intersect) {
		this.intersect = intersect;
	}

	public Node getFromNode() {
		return fromNode;
	}

	public Node getToNode() {
		return toNode;
	}

	public void setFromNode(Node fromNode) {
		this.fromNode = fromNode;
	}

	public void setToNode(Node toNode) {
		this.toNode = toNode;
	}

	public Set<String> getIntersect() {
		return intersect;
	}

	public boolean isConsidered() {
		return consider;
	}

	public void setConsidered(boolean consider) {
		this.consider = consider;
	}

	@Override
	public String toString() {
		if (consider)
			return String.format("%s->%s", fromNode, toNode);
		else
			return String.format("%s **> %s", fromNode, toNode);
	}

	@Override
	public int hashCode() {
		return (fromNode.getName() + toNode.getName()).hashCode();
	}

	@Override
	public boolean equals(Object other) {
		if (other instanceof Edge) {
			Edge o = (Edge) other;
			return fromNode.getName().equals(o.fromNode.getName()) && toNode.getName().equals(o.toNode.getName());
		}
		else
			return false;
	}
}
