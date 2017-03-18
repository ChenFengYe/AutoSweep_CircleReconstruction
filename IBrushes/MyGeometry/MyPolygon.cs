using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Drawing;
using System.Threading.Tasks;

using Loyc.Geometry;

using Emgu.CV;
using Emgu.CV.CvEnum;
using Emgu.CV.Structure;
using Emgu.CV.UI;
using Emgu.CV.Util;

using OpenTK.Graphics.OpenGL;

using SmartCanvas;

namespace MyGeometry
{
	public class MyPolygon
	{
		private MyVector3 origin;
		private MyVector3 normal;
		private MyVector3 center;
		private MyPlane belongPlane = null;
		private CoordinateFrame frame = null;
		private double area;
		private List<MyVector3> cornerPoints3d = new List<MyVector3>(); // in counterclockwise order
		private List<MyVector3> slicePoints3d = new List<MyVector3>();
		private List<MyPolygon> neighbors = new List<MyPolygon>();

		private MyVector2 center2d;
		private List<MyVector2> cornerPoints2d = new List<MyVector2>(); // in counterclockwise order
		private List<MyVector2> slicePoints2d = new List<MyVector2>(); // model points which are neareat to this polygon
		
		public MyPolygon () { }
		public MyPolygon (MyPlane p)
		{
			normal = p.Normal();
			belongPlane = p;
		}
		public MyPolygon (MyPlane p, List<MyVector3> slice_points)
		{
			normal = p.Normal();
			belongPlane = p;
			slicePoints3d.AddRange(slice_points);
		}
		public MyPolygon (MyVector3 v0, MyVector3 v1, MyPlane p)
		{
			normal = p.Normal();
			belongPlane = p;
			cornerPoints3d.Add(v0);
			cornerPoints3d.Add(v1);
		}
		public MyPolygon (MyVector3 c, MyPlane p, List<MyVector3> corners)
		{
			normal = p.Normal();
			origin = corners.First();
			belongPlane = p;
			BuildLocalFrame();
			center = c;
			double width = (corners[3] - corners[0]).Length();
			double depth = (corners[1] - corners[0]).Length();
			cornerPoints3d.AddRange(corners);
			ComputeArea();
		}
		public MyVector3 Center
		{
			get { return center; }
			set { center = value; }
		}
		public MyVector2 Center2d
		{
			get { return center2d; }
			set { center2d = value; }
		}
		public MyVector3 Origin
		{
			get { return origin; }
			set { origin = value; }
		}
		public List<MyVector3> CornerPoints3d
		{
			get { return cornerPoints3d; }
			set { cornerPoints3d = value; }
		}
		public List<MyVector2> CornerPoints2d
		{
			get { return cornerPoints2d; }
			set { cornerPoints2d = value; }
		}

		public MyPlane BelongPlane
		{
			get { return belongPlane; }
			set { belongPlane = value; }
		}
		public MyVector3 Normal
		{
			get { return normal; }
		}
		public List<MyPolygon> Neighbors
		{
			get { return neighbors; }
			set { neighbors = value; }
		}
		public List<MyVector3> SlicePoints3d
		{
			get { return slicePoints3d; }
			set { slicePoints3d = value; }
		}
		public List<MyVector2> SlicePoints2d
		{
			get { return slicePoints2d; }
			set { slicePoints2d = value; }
		}
		public void AddSlicePoint (MyVector3 sample)
		{
			this.slicePoints3d.Add(sample);
		}
		public void AddCornerPoint (MyVector3 corner)
		{
			this.cornerPoints3d.Add(corner);
		}
		public void ApproxCenter ()
		{
			// use slice points info
			MyVector3 center_tmp = new MyVector3();
			foreach (MyVector3 vert in this.slicePoints3d)
			{
				center_tmp += vert;
			}
			center_tmp /= this.slicePoints3d.Count;

			double offset = this.belongPlane.planeEquation.Offset;
			double z = this.belongPlane.planeEquation.C;
			MyVector3 p_in_plane = new MyVector3(0.0, 0.0, -offset / z);

			double bias = (center_tmp - p_in_plane).Dot(this.Normal);
			this.center = center_tmp - this.Normal * bias;
		}
		public void ResetCenter ()
		{
			// use corner info
			MyVector3 center_tmp = new MyVector3();
			foreach (MyVector3 vert in this.cornerPoints3d)
			{
				center_tmp += vert;
			}
			this.center = center_tmp / this.cornerPoints3d.Count;
		}
		public void BuildLocalFrame () 
		{
			MyVector3 X = new MyVector3(1, 0, 0);
			MyVector3 Y = new MyVector3(0, 1, 0);
			MyVector3 Z = new MyVector3(0, 0, 1);
			MyVector3 rotAxis = Z.Cross(normal).Normalize();
			double angle_cos = Z.Dot(normal);
			if (angle_cos < -1) angle_cos = -1;
			if (angle_cos > 1) angle_cos = 1;
			double rotAngle = Math.Acos(angle_cos);

			MyMatrix4d Mat = MyMatrix4d.RotationMatrix(rotAxis, rotAngle);
			MyVector3 X_new = (Mat * new MyVector4(X)).XYZ().Normalize();
			MyVector3 Y_new = (Mat * new MyVector4(Y)).XYZ().Normalize();
			frame = new CoordinateFrame(center, X_new, Y_new, normal);
		}
		public void Sort3dCornerPoints ()
		{ 
		
		}
		public void FindBoundPolygon () 
		{
			// find convexhull
			PointF[] ps = new PointF[this.slicePoints2d.Count];
			for (int i = 0; i < this.slicePoints2d.Count; i++)
			{
				PointF p = new PointF((float)this.slicePoints2d[i].x, (float)this.slicePoints2d[i].y);
				ps[i] = p;
			}
			PointF[] hull = CvInvoke.ConvexHull(ps);

			// find boundary polygon
			VectorOfPointF hull2 = new VectorOfPointF();
			hull2.Push(hull);
			VectorOfPointF poly = new VectorOfPointF();
			// when inferring # of polygon edge, the 3-rd param can be [0.0005,0.0015], than choose the best(how to define "best"??)
			CvInvoke.ApproxPolyDP(hull2, poly, 0.0003, true);
			for (int i = 0; i < poly.Size; i++)
			{
				this.cornerPoints2d.Add(new MyVector2(poly[i].X, poly[i].Y));
			}

			// unproject to 3d
			foreach (MyVector2 corner2d in this.cornerPoints2d)
			{
				MyVector3 corner3d = frame.GetPointSpaceCoord(new MyVector3(corner2d, 0.0));
				this.cornerPoints3d.Add(corner3d);
			}
		}
		public void ProjCenterTo2d ()
		{
			MyVector3 tmp = frame.GetPointLocalCoord(center);
			center2d = tmp.XY();
		}
		public void ProjSliceTo2d ()
		{
			foreach (MyVector3 vert in slicePoints3d)
			{
				MyVector3 point = frame.GetPointLocalCoord(vert);
				slicePoints2d.Add(point.XY());
			}
		}
		public void ProjCornerTo2d () 
		{ 
			foreach (MyVector3 vert in cornerPoints3d)
			{
				MyVector3 point = frame.GetPointLocalCoord(vert);
				cornerPoints2d.Add(point.XY());	
			}
		}
		public void LocalOptimizing ()
		{ 
		
		}

		public void ComputeArea ()
		{
			
		}


		public List<MyVector2> ProjectToScreen (Camera cam)
		{
			List<MyVector2> corners2d = new List<MyVector2>();
			foreach (MyVector3 vert in cornerPoints3d)
			{
				MyVector3 corner2d = cam.Project(vert.x, vert.y, vert.z);
				corners2d.Add(new MyVector2(corner2d.x, corner2d.y));
			}
			return corners2d;
		}
		public MyVector2[] ComputeLeftRightPoint (Camera cam)
		{
			MyVector2[] mp = new MyVector2[2];
			List<MyVector2> points = this.ProjectToScreen(cam);
			points.Sort((a, b) => a.x.CompareTo(b.x));
			mp[0] = points.First();
			mp[1] = points.Last();
			return mp;
		}
		public void AddNeighbor (MyPolygon r)
		{
			neighbors.Add(r);
		}
		public void Draw ()
		{
			GL.Disable(EnableCap.Lighting);
			GL.LineWidth(3.0f);
			GL.Color3(Color.Blue.R, Color.Blue.G, Color.Blue.B);
			GL.Begin(PrimitiveType.LineLoop);
			foreach (MyVector3 vert in cornerPoints3d)
			{
				GL.Vertex3(vert.ToArray());
			}
			GL.End();
			GL.PointSize(5.0f);
			GL.Begin(PrimitiveType.Points);
			GL.Vertex3(center.ToArray());
			GL.End();
			GL.Enable(EnableCap.Lighting);
		}
	}
}
