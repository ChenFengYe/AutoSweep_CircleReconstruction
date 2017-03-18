//https://www.quora.com/What-are-the-C-simplest-examples-of-implementation-of-breadth-first-search-and-depth-first-search

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GraphClass
{
    public class GraphAdjList
    {
        private readonly int V; //vertex num
        private readonly List<int>[] Adj;

        public int GetVNum { get { return V; } }
        public void AddEdge(int v, int w)
        {
            Adj[v].Add(w);
            Adj[w].Add(v);
        }
        public List<int> GetAdj(int v)
        {
            return Adj[v];
        }

        public GraphAdjList(int v)
        {
            V = v;
            Adj = new List<int>[V];

            for (int i = 0; i < V; i++)
            {
                Adj[i] = new List<int>();
            }
        }
    }

    public class DepthFirstSearch
    {
        public DepthFirstSearch()
        {

        }

        public void DFS(GraphAdjList G, ref bool[] visited, int s, ref int[] group, int groupIdx)
        {
            if (visited[s])
                return;
            visited[s] = true;
            group[s] = groupIdx;

            foreach (var w in G.GetAdj(s))
            {
                if (!visited[w])
                {
                    DFS(G, ref visited, w, ref group, groupIdx);

                    //visited[w] = true;
                    //group[w] = groupIdx;
                }
            }
        }
    }

    //public class DepthFirstSearch
    //{
    //    private bool[] visited;
    //    private int[] edgeTo;
    //    private int s;

    //    public DepthFirstSearch(GraphAdjList G, int s)
    //    {
    //        visited = new bool[G.GetVNum];
    //        edgeTo = new int[G.GetVNum];
    //        this.s = s;
    //    }

    //    public void DFS(GraphAdjList G, int v)
    //    {
    //        visited[s] = true;

    //        foreach (var w in G.GetAdj(v))
    //        {
    //            if (!visited[w])
    //            {
    //                DFS(G, w);
    //                edgeTo[w] = v;
    //            }
    //        }
    //    }

    //    public bool HasPathTo(int v)
    //    {
    //        return visited[v];
    //    }

    //    public IEnumerable<int> GetPathTo(int v)
    //    {
    //        if (!HasPathTo(v))
    //            return null;

    //        var stack = new Stack<int>();

    //        for (var x = v; x != s; x = edgeTo[x])
    //        {
    //            stack.Push(x);
    //        }

    //        stack.Push(s);

    //        return stack;
    //    }
    //}

    //public class BreadthFirstSearch
    //{
    //    public int[] edgeTo;
    //    public int[] distTo;
    //    public int s;

    //    public BreadthFirstSearch(GraphAdjList G, int s)
    //    {
    //        edgeTo = new int[G.GetVNum];
    //        distTo = new int[G.GetVNum];

    //        for (int i = 0; i < G.GetVNum; i++)
    //        {
    //            distTo[i] = -1;
    //            edgeTo[i] = -1;
    //        }

    //        this.s = s;

    //        BFS(G, s);
    //    }

    //    void BFS(GraphAdjList G, int s)
    //    {
    //        var queue = new Queue<int>();
    //        queue.Enqueue(s);
    //        distTo[s] = 0;

    //        while (queue.Count != 0)
    //        {
    //            int v = queue.Dequeue();

    //            foreach (var w in G.GetAdj(v))
    //            {
    //                if (distTo[w] == -1)
    //                {
    //                    queue.Enqueue(w);
    //                    distTo[w] = distTo[v] + 1;
    //                    edgeTo[w] = v;
    //                }
    //            }
    //        }
    //    }
    //}
}
