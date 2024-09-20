/**
 * Implementation of Dinic's network flow algorithm. The algorithm works by first constructing a
 * level graph using a BFS and then finding augmenting paths on the level graph using multiple DFSs.
 * Time Complexity: O(EV²)
 *
 * @author William Fiset, william.alexandre.fiset@gmail.com
 */
package explicit.bisim;

import prism.Evaluator;

import java.util.ArrayDeque;
import java.util.Arrays;
import java.util.Deque;

public class MaximumFlow<Value> extends NetworkFlow<Value> {

    private final int[] level;

    /**
     * Creates an instance of a flow network solver. Use the {@link #addEdge} method to add edges to
     * the graph.
     *
     * @param n - The number of nodes in the graph including source and sink nodes.
     * @param s - The index of the source node, 0 <= s < n
     * @param t - The index of the sink node, 0 <= t < n, t != s
     * @param e - An evaluator
     */
    public MaximumFlow(int n, int s, int t, Evaluator<Value> e) {
        super(n, s, t, e);
        level = new int[n];
    }

    @Override
    protected void solve() {
        // next[i] indicates the next unused edge index in the adjacency list for node i. This is part
        // of the Shimon Even and Alon Itai optimization of pruning dead ends as part of the DFS phase.
        int[] next = new int[n];

        while (bfs()) {
            Arrays.fill(next, 0);
            // Find max flow by adding all augmenting path flows.
            for (Value f = dfs(s, next, evaluator.infinity()); !evaluator.isZero(f); f = dfs(s, next, evaluator.infinity())) {
                maxFlow = evaluator.add(maxFlow, f);
            }
        }
    }

    // Do a BFS from source to sink and compute the depth/level of each node
    // which is the minimum number of edges from that node to the source.
    private boolean bfs() {
        Arrays.fill(level, -1);
        level[s] = 0;
        Deque<Integer> q = new ArrayDeque<>(n);
        q.offer(s);
        while (!q.isEmpty()) {
            int node = q.poll();
            for (Edge edge : graph.get(node)) {
                Value cap = edge.remainingCapacity();
                if (evaluator.gt(cap, evaluator.zero()) && level[edge.to] == -1) {
                    level[edge.to] = level[node] + 1;
                    q.offer(edge.to);
                }
            }
        }
        return level[t] != -1;
    }

    private Value dfs(int at, int[] next, Value flow) {
        if (at == t) return flow;
        final int numEdges = graph.get(at).size();

        for (; next[at] < numEdges; next[at]++) {
            Edge edge = graph.get(at).get(next[at]);
            Value cap = edge.remainingCapacity();
            if (evaluator.gt(cap, evaluator.zero()) && level[edge.to] == level[at] + 1) {
                Value bottleNeck = dfs(edge.to, next, min(flow, cap));
                if (evaluator.gt(bottleNeck, evaluator.zero())) {
                    edge.augment(bottleNeck);
                    minCost = evaluator.multiply(bottleNeck, edge.cost);
                    return bottleNeck;
                }
            }
        }
        return evaluator.zero();
    }

    /**
     * Returns min(a, b).
     *
     * @param a an argument.
     * @param b another argument.
     * @return the smaller of the two values, a and b.
     */
    protected Value min(Value a, Value b) {
        if (evaluator.geq(a, b)) return b; else return a;
    }
}
