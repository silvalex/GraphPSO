package pso;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Represents a node in the input/output taxonomy used by the datasets.
 * A taxonomy node records which services produce outputs that are compatible
 * with the given value, as well as which service require inputs that
 * are compatible with the given value. Taxonomy nodes are part of a graph/tree
 * structure, so each node record its parent(s) and its children. End
 * inputs satisfied by this node are tracked as a separate set.
 *
 * @author sawczualex
 */
public class TaxonomyNode {
	public Set<String> endNodeInputs = new HashSet<String>();
	public List<Node> servicesWithOutput = new ArrayList<Node>();
	public Map<Node,Set<String>> servicesWithInput = new HashMap<Node,Set<String>>();
	public String value;
	public List<TaxonomyNode> parents = new ArrayList<TaxonomyNode>();
	public List<TaxonomyNode> children = new ArrayList<TaxonomyNode>();

	/**
	 * Creates a new taxonomy node with an associated
	 * concept. All other data structures are initialised
	 * in an empty state.
	 * 
	 * @param value - The taxonomy concept this node represents
	 */
	public TaxonomyNode(String value) {
		this.value = value;
	}

	/**
	 * Gets all concepts in the subtree of this node.
	 *
	 * @return Set of concept strings.
	 */
	public Set<String> getSubsumedConcepts() {
		Set<String> concepts = new HashSet<String>();
        _getSubsumedConcepts( concepts );
		return concepts;
	}

	/**
	 * Helper method for recursively retrieving the concepts
	 * in a taxonomy subtree.
	 * 
	 * @param concepts - Set of concept strings retrieved so far
	 */
    private void _getSubsumedConcepts(Set<String> concepts) {
        if (!concepts.contains( value )) {
            concepts.add(value);
            for (TaxonomyNode child : children) {
                child._getSubsumedConcepts(concepts);
            }
        }
    }

    @Override
	/**
	 * {@inheritDoc}
	 */
    public boolean equals(Object other) {
        if (other instanceof TaxonomyNode) {
            return ((TaxonomyNode)other).value.equals( value );
        }
        return false;
    }

    @Override
	/**
	 * {@inheritDoc}
	 */
    public int hashCode() {
        return value.hashCode();
    }

    @Override
	/**
	 * {@inheritDoc}
	 */
    public String toString() {
        return value;
    }
}
