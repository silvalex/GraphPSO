package pso;

public class ListItem implements Comparable<ListItem> {
	public double score;
	public String serviceName;

	public ListItem(String serviceName, double score) {
		this.serviceName = serviceName;
		this.score = score;
	}

	@Override
	public int compareTo(ListItem o) {
		if (score > o.score)
			return -1;
		else if (score < o.score)
			return 1;
		else
			return 0;
	}

	@Override
	public String toString() {
		return "(" + score + ", " + serviceName + ")";
	}
}
