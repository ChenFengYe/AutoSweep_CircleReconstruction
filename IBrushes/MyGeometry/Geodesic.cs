using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using MyGeometry;

namespace MyGeometry
{
	/// <summary>
	/// This class is created for computing the geodesic distance on 3D mesh.
	/// Here the subdivision strategy is employed to first subdivids the mesh and
	/// then dijkstra algorithm is applied to find the distance.
	/// NoteS: 
	///			when subdivLevel = 0 the algorithm revert to dijkstra; when 
	///			when subdivLevel -> INF, the algorithm finds the exact geo-distance.
	///	CopyRight, HKUST
	///			YOUYI ZHENG @cse.ust.hk
	///			May 17, 2011
	/// </summary>
	
	public class ZyyGeodesic
	{
		/// <summary>
		/// Graph node
		/// </summary>
		public class NODE : PriorityQueueElement
		{
			public int index = -1;
			public double distance = double.MaxValue;
			public Vector3d pos;
			public NODE parent = null;
			public List<EDGE> adjEdges = new List<EDGE>();
			public NODE(int index, Vector3d p)
			{
				this.index = index;
				this.pos = p;
			}

			#region PriorityQueueElement Members
			private int pqIndex = -1;
			public int PQIndex
			{
				get
				{
					return pqIndex;
				}
				set
				{
					pqIndex = value;
				}
			}

			#endregion
			#region IComparable<PriorityQueueElement> Members

			public int CompareTo(PriorityQueueElement other)
			{
				NODE rec = other as NODE;
				if (this.distance < rec.distance) return -1;
				if (this.distance > rec.distance) return 1;
				return 0;
			}

			#endregion
		}
		/// <summary>
		/// Mesh Face, used as intermediate class
		/// </summary>
		public class FACE
		{
			public List<EDGE> edges = new List<EDGE>();
			public int index = -1;
			public FACE(int id)
			{
				this.index = id;
			}
		}
		/// <summary>
		/// Graph edge.
		/// </summary>
		public class EDGE
		{
			public List<NODE> onEdgeNodes = null;
			public int I = -1;
			public int J = -1;
			public double length = 0;
			public EDGE(int i, int j, double len)
			{
				if (i > j) { int tmp = i; i = j; j = tmp; }
				this.I = i;
				this.J = j;
				this.length = len;
			}
			public override bool Equals(object obj)
			{
				EDGE rec = obj as EDGE;
				return I == rec.I && J == rec.J;
			}
			public override int GetHashCode()
			{
				return I + J;
			}
		}


		/// <summary>
		/// Constructors & Destructors
		/// </summary>
		public ZyyGeodesic(Mesh mesh) {
			this.ConstructGraph(mesh, 0);
		}
		public ZyyGeodesic(Mesh mesh, int subdivLevel) {
			this.ConstructGraph(mesh, subdivLevel);
		}
		~ZyyGeodesic() { }


		/// <summary>
		/// Members
		/// </summary>
		private Mesh mesh = null;
		private List<NODE> nodes = null;
		private HashSet<EDGE> edges = null;


		/// <summary>
		/// This function gets geodesic distance from a single source to all vertices of the mesh.
		/// NOTES:
		///		Inserted points are ignored on return.
		/// </summary>
		/// <param name="sourceIndex"></param>
		/// <returns></returns>
		public double[] GetGeodesicDistance(int sourceIndex)
		{
			int n = this.nodes.Count;
            
			PriorityQueue queue = new PriorityQueue();

			for (int i = 0; i < n; i++)
			{
				this.nodes[i].parent = null;
				this.nodes[i].distance = double.MaxValue;
			}
            this.nodes[sourceIndex].distance = 0;

            for (int i = 0; i < n; i++)
                queue.Insert(this.nodes[i]);

            while (queue.IsEmpty() == false)
            {
                NODE minRecord = queue.DeleteMin() as NODE;

                foreach (EDGE adjEdge in minRecord.adjEdges)
                {
					double dis = adjEdge.length + minRecord.distance;
					int adj = adjEdge.I == minRecord.index ? adjEdge.J : adjEdge.I;
                    if (dis < this.nodes[adj].distance)
                    {
                        this.nodes[adj].distance = dis;
						this.nodes[adj].parent = this.nodes[minRecord.index];
						if (!queue.IsEmpty())
							queue.Update(this.nodes[adj]);
                    }
                }
            }

            double[] distance = new double[n];
            for (int i = 0; i < mesh.VertexCount; i++)
                distance[i] = this.nodes[i].distance;

            return distance;
		}

		public List<Vector3d> GetPath(int source, int target)
		{
			List<Vector3d> path = new List<Vector3d>();
			NODE node = this.nodes[target];
			while (node != null)
			{
				path.Add(node.pos);
				node = node.parent;
			}
			return path;
		}

		/// <summary>
		/// This function gets geodesic distance from multiple sources to the rest of the mesh.
		///	The points in sources are set as 0s, and the propagation goes to the entire mesh.
		/// NOTES:
		///		Inserted points are ignored on return.
		///		Once can use this function to employ segmentation by propagating from 
		///		different sources, like floodfill.
		/// </summary>
		/// <param name="sources"></param>
		/// <returns></returns>
		public double[] GetGeodesicDistance(int[] sources, double[] sourceDistances)
		{
			int n = this.nodes.Count;
			
			PriorityQueue queue = new PriorityQueue(n);

			for (int i = 0; i < n; i++)
			{
				this.nodes[i].parent = null;
				this.nodes[i].distance = double.MaxValue;
			}

			if (sourceDistances == null)
				foreach (int v in sources)
				{
					this.nodes[v].distance = 0;
				}
			else
			{
				int index = 0;
				foreach (int v in sources)
				{
					this.nodes[v].distance = sourceDistances[index++];
				}
			}
            for (int i = 0; i < n; i++)
                queue.Insert(this.nodes[i]);

            while (queue.IsEmpty() == false)
            {
                NODE minRecord = queue.DeleteMin() as NODE;

                foreach (EDGE adjEdge in minRecord.adjEdges)
                {
					double dis = adjEdge.length + minRecord.distance;
					int adj = adjEdge.I == minRecord.index ? adjEdge.J : adjEdge.I;
                    if (dis < this.nodes[adj].distance)
                    {
                        this.nodes[adj].distance = dis;
						this.nodes[adj].parent = this.nodes[minRecord.index];
						if (!queue.IsEmpty())
							queue.Update(this.nodes[adj]);
                    }
                }
            }

            double[] distance = new double[n];
            for (int i = 0; i < mesh.VertexCount; i++)
                distance[i] = this.nodes[i].distance;

            return distance;
		}


		/// <summary>
		/// This function construct the underlying graph for Dijkstra algorithm.
		/// The graph is connected with each node represent a point on mesh or mesh edges.
		/// the parameter subdiviLevel controls the subdivision resolution, i.e., how many
		/// points are inserted in each mesh edge. The larger the "subdiviLevel", the
		/// better accuracy the computed geodesic distances.
		/// NOTES:
		///		The indices for the original mesh is remained as 0 to n, the added nodes
		///		are starting with index n.
		/// </summary>
		/// <param name="mesh"></param>
		/// <param name="subdivLevel"></param>
		private void ConstructGraph(Mesh mesh, int subdivLevel)
		{
			this.mesh = mesh;

			this.nodes = new List<NODE>();

			int n = mesh.VertexCount;
			for (int i = 0; i < n; ++i)
			{
				Vector3d p = new Vector3d(mesh.VertexPos, i*3);
				this.nodes.Add(new NODE(i, p));
			}


			HashSet<EDGE> meshEdges = new HashSet<EDGE>();
			List<FACE> meshFaces = new List<FACE>();
			for (int i = 0, j = 0; i < mesh.FaceCount; ++i,j+=3)
			{
				meshFaces.Add(new FACE(i));

				int c0 = mesh.FaceIndex[j];
				int c1 = mesh.FaceIndex[j+1];
				int c2 = mesh.FaceIndex[j+2];

				Vector3d v0 = new Vector3d(mesh.VertexPos, c0*3);
				Vector3d v1 = new Vector3d(mesh.VertexPos, c1*3);
				Vector3d v2 = new Vector3d(mesh.VertexPos, c2*3);

				EDGE e0 = new EDGE(c0, c1, (v1-v0).Length());
				EDGE e1 = new EDGE(c1, c2, (v2-v1).Length());
				EDGE e2 = new EDGE(c2, c0, (v2-v0).Length());

				meshEdges.Add(e0);
				meshEdges.Add(e1);
				meshEdges.Add(e2);
			}

			/// now build the graph, first subdvids edges, then connect
			if (subdivLevel > 0)
			{ 
				foreach (EDGE e in meshEdges)
				{
					List<int> adjFaces = mesh.AdjVF[e.I].Intersect(mesh.AdjVF[e.J]).ToList();
					foreach (int f in adjFaces)
					{
						meshFaces[f].edges.Add(e);
					}
				}

				this.edges = new HashSet<EDGE>();
				
				int index = n;
				foreach (EDGE e in meshEdges)
				{
					NODE node1 = this.nodes[e.I], node2 = this.nodes[e.J];
					e.onEdgeNodes = new List<NODE>();
					for (int i = 0; i < subdivLevel; ++i)
					{
						double ratio = (double)(i + 1) / (subdivLevel + 1);
						Vector3d p = node2.pos * ratio + node1.pos*(1-ratio);
						NODE node = new NODE(index++, p);
						this.nodes.Add(node);
						e.onEdgeNodes.Add(node);
					}
				}

				foreach (EDGE e in meshEdges)
				{
					NODE start = this.nodes[e.I], end = this.nodes[e.J];
					int k = e.onEdgeNodes.Count;
					if (k > 0)
					{
						NODE n1 = e.onEdgeNodes[0], n2 = e.onEdgeNodes[k - 1];

						EDGE e1 = new EDGE(start.index, n1.index, (start.pos - n1.pos).Length());
						EDGE e2 = new EDGE(n2.index, end.index, (end.pos - n2.pos).Length());
						this.edges.Add(e1);
						this.edges.Add(e2);
					}
					for (int i = 0; i < k - 1; ++i)
					{
						NODE n1 = e.onEdgeNodes[i], n2 = e.onEdgeNodes[i + 1];
						EDGE ee = new EDGE(n1.index, n2.index, (n1.pos - n2.pos).Length());
						this.edges.Add(ee);
					}
				}
				foreach (FACE face in meshFaces)
				{
					// connect points on edges
					foreach (EDGE e1 in face.edges)
					{
						foreach (NODE node1 in e1.onEdgeNodes)
						{
							foreach (EDGE e2 in face.edges)
							{
								if (e2 == e1) continue;

								// find the node opposite to e1
								int opposite = e2.I;
								if (opposite == e1.I || opposite == e1.J)
								{
									opposite = e2.J;
								}
								this.edges.Add(
									new EDGE(opposite, node1.index,
										(node1.pos - this.nodes[opposite].pos).Length()
									)
								);

								// connect points on one edge to points on another
								foreach (NODE node2 in e2.onEdgeNodes)
								{
									EDGE ee = new EDGE(node1.index, node2.index,
										(node1.pos - node2.pos).Length()
									);
									this.edges.Add(ee);
								}
							}
						}
					}
				}
			}
			else
			{
				this.edges = meshEdges;
			}

			foreach (EDGE e in this.edges)
			{
				NODE node1 = this.nodes[e.I], node2 = this.nodes[e.J];
				node1.adjEdges.Add(e);
				node2.adjEdges.Add(e);
			}
		}
		private void Sort(int i, int j)
		{
			if (i > j)
			{
				int tmp = i;
				i = j;
				j = tmp;
			}
		}
	}
}
