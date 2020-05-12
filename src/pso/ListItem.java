package pso;

/**
 * Class that represents an item to be added to the candidate
 * list of nodes during the decoding process. It associates
 * a given service to a given score, which allows the candidate
 * list to be sorted by importance.
 * 
 * @author sawczualex
 *
 */
public class ListItem implements Comparable<ListItem> {
	public double score;
	public String serviceName;

	/**
	 * Creates a new list item, pairing the given service
	 * name and score.
	 * 
	 * @param serviceName
	 * @param score
	 */
	public ListItem(String serviceName, double score) {
		this.serviceName = serviceName;
		this.score = score;
	}

	@Override
	/**
	 * {@inheritDoc}
	 */
	public int compareTo(ListItem o) {
		if (score > o.score)
			return -1;
		else if (score < o.score)
			return 1;
		else
			return 0;
	}

	@Override
	/**
	 * {@inheritDoc}
	 */
	public String toString() {
		return "(" + score + ", " + serviceName + ")";
	}
}
