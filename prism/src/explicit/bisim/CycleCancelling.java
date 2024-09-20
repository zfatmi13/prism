/**
 * Implementation of Dinic's network flow algorithm. The algorithm works by first constructing a
 * level graph using a BFS and then finding augmenting paths on the level graph using multiple DFSs.
 * Time Complexity: O(EV²)
 *
 * @author William Fiset, william.alexandre.fiset@gmail.com
 */
package explicit.bisim;

import prism.Evaluator;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

public class CycleCancelling<Value> extends MaximumFlow<Value> {

    /**
     * Creates an instance of a flow network solver. Use the {@link #addEdge} method to add edges to
     * the graph.
     *
     * @param n - The number of nodes in the graph including source and sink nodes.
     * @param s - The index of the source node, 0 <= s < n
     * @param t - The index of the sink node, 0 <= t < n, t != s
     * @param e - An evaluator
     */
    public CycleCancelling(int n, int s, int t, Evaluator<Value> e) {
        super(n, s, t, e);
    }

    @Override
    public void solve() {
        super.solve();
        if (evaluator.isOne(maxFlow) && evaluator.isOne(minCost)) {
            System.out.println("     Cycle cancelling");
            // Reduce the cost
            List<Edge> path = new LinkedList<>();
            // Run only one iteration, since we do not need the minimum cost
            getNegativeCycle(path);
            // while (getNegativeCycle(path))
            //    path.clear();
            // printGraph();
        }
    }

    private void printGraph() {
        System.out.println("Graph:");
        List<List<Edge>> g = getGraph();
        for (List<Edge> edges : g) {
            for (Edge e : edges) {
                System.out.println(e.toString(s, t));
                // System.out.println(e.residual.toString(s, t));
            }
        }
    }

    private boolean getNegativeCycle(List<Edge> path){
        List<Value> distance = new ArrayList<Value>(n);
        List<Edge> previous = new ArrayList<Edge>(n);
        for (int i = 0; i < n; i++) {
            distance.add(evaluator.infinity());
            previous.add(null);
        }
        // possibly change this because s and t are n-2 and n-1
        distance.set(s, evaluator.zero());
        distance.set(t, evaluator.one());
        for (Edge edge : graph.get(s)) {
            distance.set(edge.to, evaluator.add(distance.get(s), edge.cost));
        }
        // For each vertex, relax all the edges in the graph, O(VE)
        for (int i = 0; i < n - 1; i++) {
            for (int from = 0; from < n; from++) {
                for (Edge edge : graph.get(from)) {
                    Value d = evaluator.add(distance.get(from), edge.cost);
                    if (evaluator.gt(edge.remainingCapacity(), evaluator.zero()) && evaluator.gt(distance.get(edge.to), d)) {
                        distance.set(edge.to, d);
                        previous.set(edge.to, edge);
                    }
                }
            }
        }
        // Check for negative-weight cycles
        for (int from = 0; from < n; from++) {
            for (Edge edge : graph.get(from)) {
                if (evaluator.gt(edge.remainingCapacity(), evaluator.zero()) && evaluator.gt(distance.get(edge.to), evaluator.add(distance.get(from), edge.cost))) {
                    previous.set(edge.to, edge);
                    // Negative cycle exists
                    markAllNodesAsUnvisited();
                    visit(edge.to);
                    while (!visited(edge.from)) {
                        visit(edge.from);
                        edge = previous.get(edge.from);
                    }
                    // edge.from is a vertex in a negative cycle, find the cycle
                    path.add(edge);
                    Value bottleNeck = edge.remainingCapacity();
                    from = edge.from;
                    edge = previous.get(from);
                    // System.out.print(from + " ");
                    while (edge.from != from) {
                        // System.out.print(edge.from + " ");
                        path.add(edge);
                        bottleNeck = min(bottleNeck, edge.remainingCapacity());
                        edge = previous.get(edge.from);
                    }
                    // System.out.println();
                    // Retrace path while augmenting the flow
                    for (Edge e : path) {
                        e.augment(bottleNeck);
                        minCost = evaluator.multiply(bottleNeck, edge.cost);
                    }
                    return true;
                }
            }
        }
        return false;
    }
}
