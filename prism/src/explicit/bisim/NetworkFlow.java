/**
 * @author William Fiset, william.alexandre.fiset@gmail.com
 */
package explicit.bisim;

import prism.Evaluator;

import java.util.ArrayList;
import java.util.List;

public abstract class NetworkFlow<Value> {

    protected Evaluator<Value> evaluator;

    public class Edge {
        public int from, to;
        public Edge residual;
        public Value flow;
        public final Value capacity;
        public final Value cost;

        public Edge(int from, int to, Value capacity, Value cost) {
            this.from = from;
            this.to = to;
            this.capacity = capacity;
            this.cost = cost;
            this.flow = evaluator.zero();
        }

        public boolean isResidual() {
            return evaluator.isZero(capacity);
        }

        public Value remainingCapacity() {
            return evaluator.subtract(capacity, flow);
        }

        public void augment(Value bottleNeck) {
            flow = evaluator.add(flow, bottleNeck);
            residual.flow = evaluator.subtract(residual.flow, bottleNeck);
        }

        public String toString(int s, int t) {
            String u = (from == s) ? "s" : ((from == t) ? "t" : String.valueOf(from));
            String v = (to == s) ? "s" : ((to == t) ? "t" : String.valueOf(to));
            return String.format(
                    "Edge %s -> %s | flow = %f | capacity = %f | cost = %f | is residual: %s",
                    u, v, evaluator.toDouble(flow), evaluator.toDouble(capacity), evaluator.toDouble(cost), isResidual());
        }
    }

    // Inputs: n = number of nodes, s = source, t = sink
    protected final int n, s, t;

    protected Value maxFlow;
    protected Value minCost;

    protected List<List<Edge>> graph;

    // 'visited' and 'visitedToken' are variables used for graph sub-routines to
    // track whether a node has been visited or not. In particular, node 'i' was
    // recently visited if visited[i] == visitedToken is true. This is handy
    // because to mark all nodes as unvisited simply increment the visitedToken.
    private int visitedToken = 1;
    private int[] visited;

    // Indicates whether the network flow algorithm has run. We should not need to
    // run the solver multiple times, because it always yields the same result.
    private boolean solved;

    /**
     * Creates an instance of a flow network solver. Use the {@link #addEdge} method to add edges to
     * the graph.
     *
     * @param n - The number of nodes in the graph including source and sink nodes.
     * @param s - The index of the source node, 0 <= s < n
     * @param t - The index of the sink node, 0 <= t < n, t != s
     * @param e - An evaluator
     */
    public NetworkFlow(int n, int s, int t, Evaluator<Value> e) {
        this.n = n;
        this.s = s;
        this.t = t;
        initializeGraph();
        visited = new int[n];
        this.evaluator = e;
        this.maxFlow = evaluator.zero();
        this.minCost = evaluator.zero();
    }

    // Construct an empty graph with n nodes including the source and sink nodes.
    private void initializeGraph() {
        graph = new ArrayList<List<Edge>>(n);
        for (int i = 0; i < n; i++) graph.add(new ArrayList<Edge>());
    }

    /**
     * Adds a directed edge (and residual edge) to the flow graph.
     *
     * @param from - The index of the node the directed edge starts at.
     * @param to - The index of the node the directed edge ends at.
     * @param capacity - The capacity of the edge.
     */
    public void addEdge(int from, int to, Value capacity) {
        addEdge(from, to, capacity, evaluator.zero());
    }

    /** Cost variant of {@link #addEdge(int, int, Value)} for min-cost max-flow */
    public void addEdge(int from, int to, Value capacity, Value cost) {
        if (evaluator.gt(evaluator.zero(), capacity)) throw new IllegalArgumentException("Capacity < 0");
        if (evaluator.gt(evaluator.zero(), cost)) throw new IllegalArgumentException("Cost < 0");
        Edge e1 = new Edge(from, to, capacity, cost);
        Edge e2 = new Edge(to, from, evaluator.zero(), evaluator.subtract(evaluator.zero(), cost));
        e1.residual = e2;
        e2.residual = e1;
        graph.get(from).add(e1);
        graph.get(to).add(e2);
    }

    // Marks node 'i' as visited.
    public void visit(int i) {
        visited[i] = visitedToken;
    }

    // Returns whether node 'i' has been visited.
    public boolean visited(int i) {
        return visited[i] == visitedToken;
    }

    // Resets all nodes as unvisited. This is especially useful to do
    // between iterations of finding augmenting paths, O(1)
    public void markAllNodesAsUnvisited() {
        visitedToken++;
    }

    /**
     * Returns the graph after the solver has been executed. This allows you to inspect the {@link
     * Edge#flow} compared to the {@link Edge#capacity} in each edge. This is useful if you want to
     * figure out which edges were used during the max flow.
     */
    public List<List<Edge>> getGraph() {
        execute();
        return graph;
    }

    // Returns the maximum flow from the source to the sink.
    public Value getMaxFlow() {
        execute();
        return maxFlow;
    }

    // Returns the min cost from the source to the sink.
    // NOTE: This method only applies to min-cost max-flow algorithms.
    public Value getMinCost() {
        execute();
        return minCost;
    }

    // Wrapper method that ensures we only call solve() once
    public void execute() {
        if (solved) return;
        solved = true;
        solve();
    }

    // Method to implement which solves the network flow problem.
    protected abstract void solve();

    // Method to reset the problem between perf runs
    public void reset() {
        maxFlow = evaluator.zero();
        minCost = evaluator.zero();
        solved = false;
        for (List<Edge> list : graph) {
            for (Edge edge : list) {
                edge.flow = evaluator.zero();
            }
        }
        markAllNodesAsUnvisited();
    }
}
