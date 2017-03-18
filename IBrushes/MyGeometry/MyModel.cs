using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Threading.Tasks;

using Accord.Math;
using Accord.Math.Distances;
using Accord.MachineLearning;
using Accord.Statistics.Distributions.DensityKernels;

namespace MyGeometry
{
	public class MyDistance : IMetric<double[]>
	{
		public double Distance (double[] a, double[] b)
		{
			double squaresum = 0;
			for (int i = 0; i < a.Length; i++)
			{
				squaresum += Math.Pow(a[i] - b[i], 2);
			}
			return squaresum;
		}
	}
	public class MyModel
	{
		private List<MyVector3> vertices = new List<MyVector3>();

		public MyModel (string filename)
		{
			List<MyVector3> model = new List<MyVector3>();
			using (StreamReader sr = new StreamReader(filename, Encoding.Default))
			{
				int point_num = int.Parse(sr.ReadLine());
				string[] arr = null;
				for (int i = 0; i < point_num; i++)
				{
					arr = sr.ReadLine().Split(new char[] { ' ' });
					model.Add(new MyVector3(Convert.ToDouble(arr[0]), Convert.ToDouble(arr[1]), Convert.ToDouble(arr[2])));
				}
			}
		}
		public MyModel (List<MyVector3> vs)
		{
			vertices.AddRange(vs);
		}
		public List<MyVector3> Vertices
		{
			get { return vertices; }
			set { vertices = value; }
		}
		public void AddVertexToModel (MyVector3 vert)
		{
			vertices.Add(vert);
		}
		public void SaveModelToFile (string commonprefix)
		{
			using (StreamWriter sw = new StreamWriter(commonprefix + ".model", false))
			{
				sw.WriteLine(vertices.Count);
				foreach (MyVector3 vert in vertices)
				{
					sw.WriteLine(vert);
				}
				sw.Flush();
				sw.Close();
			}
		}
		public void SaveModelToPly (string commonprefix)
		{
			using (StreamWriter sw = new StreamWriter(commonprefix + ".ply", false))
			{
				sw.WriteLine("ply");
				sw.WriteLine("format ascii 1.0");
				sw.WriteLine("comment");
				sw.WriteLine("element vertex {0}", vertices.Count);
				sw.WriteLine("property float x");
				sw.WriteLine("property float y");
				sw.WriteLine("property float z");
				sw.WriteLine("element face 0");
				sw.WriteLine("property list uchar int vertex_indices");
				sw.WriteLine("end_header");
				foreach (MyVector3 vert in vertices)
				{
					sw.WriteLine(vert);
				}
				sw.Flush();
				sw.Close();
			}
		}
		public static double ComputeHeight(MyModel m, MyPlane p)
		{
			double height = double.MinValue;
			for (int i = 0; i < m.Vertices.Count; i++)
			{
				double dist = p.DistanceToPoint(m.Vertices[i]);
				if (dist > height)
				{
					height = dist;
				}
			}
			return height;
		}
		public void Slicing(List<MyPolygon> planes)
		{

		}
	}
}
