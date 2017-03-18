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
using PolygonDecompose;

namespace MyGeometry
{
	public enum CTYPE 
	{ 
		FEATURE, MIDDLE, UNCERTAIN
	};

	public class MyCompare2dX : IComparer<MyVector2>
	{
		public int Compare (MyVector2 a, MyVector2 b)
		{
			return a.x.CompareTo(b.x);
		}
	}
	public class MyCompare2dY : IComparer<MyVector2>
	{
		public int Compare (MyVector2 a, MyVector2 b)
		{
			return a.y.CompareTo(b.y);
		}
	}
	public class MyCompare3dY : IComparer<MyVector3>
	{
		public int Compare (MyVector3 a, MyVector3 b)
		{
			return a.y.CompareTo(b.y);
		}
	}
	public class MyCompareSegmentLowY : IComparer<CLineSegment>
	{
		public int Compare (CLineSegment a, CLineSegment b)
		{
			return a.StartPoint.y.CompareTo(b.StartPoint.y);
		}
	}
	public class MyCircle
	{
		private MyVector3 center;
		private MyVector2 center2d;
		private MyVector3 normal;
		private double radius;
		private double area;
		private MyPlane belongPlane = null;
		private List<MyVector3> circlePoints3d = new List<MyVector3>(); // better rename to be "circleSamples", since they are all sampled points
		private List<MyVector3> circleSlice3d = new List<MyVector3>();
		private List<MyVector2> circleSlice2d = new List<MyVector2>();
		private List<MyCircle> neighbors = new List<MyCircle>();

		public MyCircle ()
		{
		}
		public MyCircle (MyVector3 c, double r)
		{
			center = c;
			radius = r;
			ComputeArea();
		}
		public MyCircle (MyVector3 c, double r, MyVector3 n, List<MyVector3> cp)
		{
			center = c;
			radius = r;
			normal = n;
			circlePoints3d = cp;
			ComputeArea();
		}
		public MyCircle (MyVector3 c, double r, MyPlane p)
		{
			center = c;
			radius = r;
			normal = p.Normal();
			belongPlane = p;
			ComputeArea();
			SampleCirclePoints3d(c, p.Normal(), r, 100);
		}
		public MyCircle (MyVector3 c, double r, MyPlane p, List<MyVector3> cp)
		{
			center = c;
			radius = r;
			normal = p.Normal();
			belongPlane = p;
			circlePoints3d = cp;
			ComputeArea();
		}
		public MyCircle (MyCircle circle)
		{
			center = circle.Center;
			radius = circle.Radius;
			belongPlane = circle.BelongPlane;
			circlePoints3d = circle.CirclePoints;
			ComputeArea();
			neighbors = circle.Neighbors;
		}
		public MyVector3 Center   // If use same name,"public MyVector3 Center { get; set; }" is ok
		{
			get { return center; }
			set { center = value; }
		}
		public MyVector2 Center2d
		{
			get { return center2d; }
			set { center2d = value; }
		}
		public MyVector3 Normal
		{
			get { return normal; }
		}
		public double Radius
		{
			get { return radius; }
			set { radius = value; }
		}
		public double Area
		{
			get { return area; }
			set { area = value; }
		}
		public MyPlane BelongPlane
		{
			get { return belongPlane; }
			set { belongPlane = value; }
		}
		public List<MyVector3> CirclePoints
		{
			get { return circlePoints3d; }
			set { circlePoints3d = value; }
		}
		public List<MyVector3> CircleSlice3d
		{
			get { return circleSlice3d; }
			set { circleSlice3d = value; }
		}
		public List<MyVector2> CircleSlice2d
		{
			get { return circleSlice2d; }
			set { circleSlice2d = value; }
		}
		public List<MyCircle> Neighbors
		{
			get { return neighbors; }
			set { neighbors = value; }
		}

        private void SampleCirclePoints3d(MyVector3 center, MyVector3 normal, double radius, int sample)
        {
            MyPlane circleplane = new MyPlane(center, normal);
            //MyVector3 pointonplane = circleplane.ProjectPoint(center - new MyVector3(0.01, 0.01, 0.01));
            //MyVector3 u = pointonplane - center;
            MyVector3 pointonplane = circleplane.ProjectPoint(center + new MyVector3(-0.1, 0, 0));
            MyVector3 u = pointonplane - center;
            MyVector3 v = u.Cross(normal);
            u = u.Normalize();
            v = v.Normalize();
            int num = sample;

            double theta = 2.0 * Math.PI / (num - 1);

            for (int j = 0; j < num; j++)
            {
                double x = center.x + radius * (u.x * Math.Cos(j * theta) + v.x * Math.Sin(j * theta));
                double y = center.y + radius * (u.y * Math.Cos(j * theta) + v.y * Math.Sin(j * theta));
                double z = center.z + radius * (u.z * Math.Cos(j * theta) + v.z * Math.Sin(j * theta));
                this.circlePoints3d.Add(new MyVector3(x, y, z));
            }
        }

		public void AddToCircleSlice (MyVector3 sample)
		{
			this.circleSlice3d.Add(sample);
		}
		public void Project3dToBelongPlane ()
		{ 
			// prepare 2d coordinate
			MyVector3 circle_norm = this.belongPlane.Normal();
			MyVector3 X = new MyVector3(1, 0, 0);
			MyVector3 Y = new MyVector3(0, 1, 0);
			MyVector3 Z = new MyVector3(0, 0, 1);
			MyVector3 rotAxis = Z.Cross(circle_norm).Normalize();
			double angle_cos = Z.Dot(circle_norm);
			if (angle_cos > 1) angle_cos = 1;
			if (angle_cos < -1) angle_cos = -1;
			double rotAngle = Math.Acos(angle_cos);

			MyMatrix4d Mat = MyMatrix4d.RotationMatrix(rotAxis, rotAngle);
			MyVector3 X_new = (Mat * new MyVector4(X)).XYZ().Normalize();
			MyVector3 Y_new = (Mat * new MyVector4(Y)).XYZ().Normalize();
			CoordinateFrame frame = new CoordinateFrame(this.belongPlane.planeCenter, X_new, Y_new, circle_norm);

			// projection and denoise(leave out far away points)
			MyVector3 tmpc = frame.GetPointLocalCoord(this.center);
			this.center2d = new MyVector2(tmpc.x, tmpc.y);
			for (int i = 0; i < this.circleSlice3d.Count; i++)
			{
				MyVector3 vert = this.circleSlice3d[i];
				MyVector3 tmp = frame.GetPointLocalCoord(vert);
				MyVector2 tmp2 = new MyVector2(tmp.x, tmp.y);
				double dist = Math.Abs((tmp2 - this.center2d).Length() - this.radius);
				if (dist > 0.5 * this.radius && dist < 0.75 * this.radius)
				{
					this.circleSlice3d.RemoveAt(i);
					i--;
				}
				else
					this.circleSlice2d.Add(tmp2);
			}
		}
		public void ComputeArea ()
		{
			this.area = Math.PI * Math.Pow(radius, 2.0);
		}
		public List<MyVector2> ProjectToScreen (Camera cam)
		{ 
			List<MyVector2> points2d = new List<MyVector2>();
			foreach (MyVector3 vert in circlePoints3d)
			{
				MyVector3 point2d = cam.Project(vert.x, vert.y, vert.z);
				points2d.Add(new MyVector2(point2d.x, point2d.y));
			}
			return points2d;
		}
		// no use
		public MyVector2[] ComputeLeftRightPoint (Camera cam)
		{
			MyVector2[] mp = new MyVector2[2];
			List<MyVector2> points = this.ProjectToScreen(cam);
			//MyCompare2dX mc = new MyCompare2dX();
			points.Sort((a, b) => a.x.CompareTo(b.x));
			//points.Sort(mc);
			mp[0] = points.First(); mp[1] = points.Last();
			return mp;
		}
		public void AddNeighbor (MyCircle c)
		{
			neighbors.Add(c);
		}
		CTYPE CircleType ()
		{
			CTYPE tmp = CTYPE.UNCERTAIN;
			switch (neighbors.Count)
			{
				case 2:
					tmp = CTYPE.MIDDLE;
					break;
				case 1:
					tmp = CTYPE.FEATURE;
					break;
				default:
					break;
			}
			return tmp;
		}
		public void Draw ()
		{
			GL.Disable(EnableCap.Lighting);
			GL.LineWidth(3.0f);
			GL.Color3(Color.Blue.R, Color.Blue.G, Color.Blue.B);
			GL.Begin(PrimitiveType.LineLoop);
			foreach (MyVector3 vert in circlePoints3d)
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
