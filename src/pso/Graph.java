package pso;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


public class Graph {

	public Map<String, Node> nodeMap = new HashMap<String, Node>();
	public List<Edge> edgeList = new ArrayList<Edge>();

	@Override
	public boolean equals(Object other) {
		if (other instanceof Graph) {
			return toString().equals(other.toString());
		}
		return false;
	}

	@Override
	public int hashCode() {
		return toString().hashCode();
	}

	@Override
	public String toString() {
		StringBuilder builder = new StringBuilder();
		for(Edge e: edgeList) {
			builder.append(e);
			builder.append(" ");
		}
		return builder.toString();
	}

}
