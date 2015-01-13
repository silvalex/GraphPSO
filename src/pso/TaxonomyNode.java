package pso;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Represents a node in the input/output taxonomy
 * used by the WSC dataset.
 *
 * @author sawczualex
 */
public class TaxonomyNode {
	public TaxonomyNode parent;
	public Set<String> endNodeInputs = new HashSet<String>();
	public List<Node> servicesWithOutput = new ArrayList<Node>();
	public List<Node> servicesWithInput = new ArrayList<Node>();
	public String value;
	public List<TaxonomyNode> children = new ArrayList<TaxonomyNode>();

	public TaxonomyNode(String value) {
		this.value = value;
	}

	public TaxonomyNode(String value, TaxonomyNode parent) {
		this.value = value;
		this.parent = parent;
	}

	/**
	 * Gets all concepts subsumed by this node (i.e. all
	 * concepts in its subtree).
	 *
	 * @return Set of concepts
	 */
	public Set<String> getSubsumedConcepts() {
		Set<String> concepts = new HashSet<String>();
        _getSubsumedConcepts( concepts );
		return concepts;
	}

    private void _getSubsumedConcepts(Set<String> concepts) {
        concepts.add(value);
        for (TaxonomyNode child : children) {
            child._getSubsumedConcepts(concepts);
        }
    }
}
