using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Drawing;
using System.Threading.Tasks;

using Loyc.Geometry;

using OpenTK.Graphics.OpenGL;

using SmartCanvas;

namespace MyGeometry
{
	public class MyGCuboid
	{
		public List<MyPolygon> polyList = new List<MyPolygon>();

		public MyGCuboid () { }
		public MyGCuboid (MyPolygon poly)
		{
			polyList.Add(poly);
		}
		public MyGCuboid (List<MyPolygon> list)
		{
			polyList.AddRange(list);
		}
		public List<MyPolygon> PolygonList
		{
			get { return polyList; }
			set { polyList = value; }
		}
		public void AddPolygon (MyPolygon poly)
		{
			polyList.Add(poly);
		}
		public void SliceModel (MyModel m)
		{
			List<MyVector3> vertices = m.Vertices;
			int slice_num = this.polyList.Count;
			foreach (MyVector3 vert in vertices)
			{
				double dist_min = int.MaxValue;
				int dist_min_idx = 0;
				for (int i = 0; i < slice_num; i++)
				{
					double dist = polyList[i].BelongPlane.DistanceToPoint(vert);
					if (dist < dist_min)
					{
						dist_min = dist;
						dist_min_idx = i;
					}
				}
				polyList[dist_min_idx].AddSlicePoint(vert);
			}
			RemovePolygonsInTail();
		}
		public void RemovePolygonsInTail ()
		{
			// remove last several slices with few vertices
			for (int i = polyList.Count - 1; i > 0; i--)
			{
				if (polyList[i].Neighbors.Count == 1)
				{
					if (polyList[i].SlicePoints3d.Count < 10)
					{
						polyList[i - 1].Neighbors.Remove(polyList[i]);
						polyList.RemoveAt(i);
					}
				}
				else
					break;
			}
		}
		public void Initializing ()
		{
			this.TopPlaneProcess();
			this.MidPlaneProcess();
		}
		public void Optimizing ()
		{ 
			// do optimize in 2d space
			// each gcuboid has a list of polygon(defined in MyPolygon.cs)
			// input data : polygon's 2d information, including center2d, cornerPoints2d and slicePoints2d 
			// output data : new position of cornerPoints2d
		}
		private void TopPlaneProcess ()
		{ 
			// top plane info estimating
			this.polyList.Last().ApproxCenter(); // center -get
			this.polyList.Last().BuildLocalFrame(); //frame -get
			this.polyList.Last().ProjSliceTo2d();  // slicePoints2d -get
			this.polyList.Last().FindBoundPolygon(); // cornerPoints3d, cornerPoints2d -get

			this.polyList.Last().ResetCenter();
			this.polyList.Last().ProjCenterTo2d();
		}
		private void MidPlaneProcess ()
		{
			for (int k = this.polyList.Count - 2; k >= 0; k--)
			{
				// last plane corners
				MyPolygon p_last = this.polyList[k + 1];
				List<MyVector3> corners_last = p_last.CornerPoints3d;
				// this plane info
				MyPolygon p_now = this.polyList[k];
				double offset = p_now.BelongPlane.planeEquation.Offset;
				double z = p_now.BelongPlane.planeEquation.C;
				MyVector3 p_in_plane = new MyVector3(0.0, 0.0, -offset / z);

				foreach (MyVector3 c_last in corners_last)
				{
					double bias = (c_last - p_in_plane).Dot(p_now.Normal);
					MyVector3 c_now = c_last - p_now.Normal * bias;
					p_now.AddCornerPoint(c_now);
				}

				// local optimization
				p_now.LocalOptimizing();

				p_now.ResetCenter();
				p_now.BuildLocalFrame();
				p_now.ProjCenterTo2d();
				p_now.ProjSliceTo2d();
				p_now.ProjCornerTo2d();
			}
		}
		public void Draw ()
		{
			foreach (MyPolygon poly in polyList)
			{
				poly.Draw();
			}
		}
	}
}
