using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using MyGeometry;

namespace i3_ImageManipulator
{
	public enum CurveStyle {circle = 0, halfcircle, ellipse, rectangle, polyline, cubicspline, none}

	public class CurvePoint
	{
		//public double lambda = -1;
		//public double dist2skel = -1;
		//public Vector2d controllermappoint;

		public int index = -1;
		public Vector2d pos;
		public CurvePoint(Vector2d p, int id)
		{
			this.index = id;
			this.pos = p;
		}
		public CurvePoint(Vector2d p)
		{
			this.pos = p;
		}
	}

	public class Curve
	{
		public bool isClosed = true;
		public CurveStyle type = CurveStyle.none;
		public List<CurvePoint> controlPoints = null;
		public List<CurvePoint> curvePoints = null;
		public int selectedControlPointIndex = -1;
		public Vector2d selectedControlPointPos;
		public double SelectControlPoint(Vector2d point)
		{
			double mindis = double.MaxValue; int index = 0;
			CurvePoint pt = null;
			foreach (CurvePoint cp in this.controlPoints)
			{
				double d = (cp.pos - point).Length();
				if (d < mindis)
				{
					this.selectedControlPointIndex = index;
					mindis = d;
					pt = cp;
				}
				index++;
			}
			this.selectedControlPointPos = new Vector2d(pt.pos.x, pt.pos.y);

			return mindis;
		}
		public virtual void GetCurvePoints() { }
		public virtual void Update(Vector2d displacement) { }
		public Curve()
		{
			this.controlPoints = new List<CurvePoint>();
		}
		~Curve() {}
	}

	public class Circle : Curve
	{
		public double radius;
		public Vector2d center;
		public int numofPoints;
		public Circle(Vector2d c, double r, int slices)
		{
			this.type = CurveStyle.circle;
			this.center = c;
			this.radius = r;
			this.numofPoints = slices;
			this.isClosed = true;
			this.GetCurvePoints();
		}
		public void Update(Vector2d c, double r)
		{
			this.center = c;
			this.radius = r;
			double dtheta = Math.PI * 2 / this.numofPoints;
			int K = this.numofPoints / 4;
			for (int i = 0; i < this.numofPoints; ++i)
			{
				double theta = i * dtheta;
				double x = this.radius * Math.Cos(theta);
				double y = this.radius * Math.Sin(theta);
				Vector2d p = this.center + new Vector2d(x, y);
				this.curvePoints[i].pos = p;
			}
		}
		public override void GetCurvePoints()
		{
			this.curvePoints = new List<CurvePoint>();
			this.controlPoints = new List<CurvePoint>();
			double dtheta = Math.PI * 2 / this.numofPoints;
			int K = this.numofPoints / 4;
			for (int i = 0; i < this.numofPoints; ++i)
			{
				double theta = i * dtheta;
				double x = this.radius * Math.Cos(theta);
				double y = this.radius * Math.Sin(theta);
				Vector2d p = this.center + new Vector2d(x, y);
				CurvePoint cp = new CurvePoint(p, i);
				this.curvePoints.Add(cp);
				if (i % K == 0)
					this.controlPoints.Add(cp);
			}
		}
		public override void Update(Vector2d displacement)
		{
			if (this.selectedControlPointIndex == -1) return;

			Vector2d oldpos = new Vector2d(this.selectedControlPointPos.x, this.selectedControlPointPos.y);
			if (this.selectedControlPointIndex == 0 || this.selectedControlPointIndex == 2)
			{
				oldpos.x += displacement.x;
			}
			else
			{
				oldpos.y += displacement.y;
			}
			this.radius = (oldpos - this.center).Length();
			this.GetCurvePoints();
		}
	}

	public class HalfCircle : Curve
	{
		public double radius;
		public Vector2d center;
		public Vector2d xaxis = new Vector2d(1,0);
		public int numofPoints = 30;
		public HalfCircle(int slices)
		{
			this.type = CurveStyle.halfcircle;
			this.numofPoints = slices;
			this.isClosed = false;
			this.CreateDefault();
		}
		public HalfCircle(Vector2d c, Vector2d xaxis, double radius, int slices)
		{
			this.isClosed = false;
			this.type = CurveStyle.halfcircle;
			this.numofPoints = slices;
			this.Create(c, xaxis, radius);
		}
		public void CreateDefault()
		{
			this.center = new Vector2d(0,0);
			this.radius = 1.0;
			
			this.curvePoints = new List<CurvePoint>();
			this.controlPoints = new List<CurvePoint>();

			double dtheta = Math.PI / this.numofPoints;
			int K = this.numofPoints / 2;
			for (int i = 0; i < this.numofPoints; ++i)
			{
				double theta = i * dtheta;
				double x = this.radius * Math.Cos(theta);
				double y = this.radius * Math.Sin(theta);
				Vector2d p = this.center + new Vector2d(x, y);
				CurvePoint cp = new CurvePoint(p, i);
				this.curvePoints.Add(cp);
				if (i % K == 0)
					this.controlPoints.Add(cp);
			}
		}
		public void Create(Vector2d c, Vector2d xaxis, double radius)
		{
			this.center = c;
			this.xaxis = xaxis;
			this.radius = radius;

			this.curvePoints = new List<CurvePoint>();
			this.controlPoints = new List<CurvePoint>();

			double dtheta = Math.PI / this.numofPoints;
			double cos = Math.Cos(dtheta);
			double sin = Math.Sin(dtheta);
			Matrix2d R = new Matrix2d(cos,-sin,sin,cos);
			Vector2d uv = this.xaxis * radius;
			int K = this.numofPoints / 2;
			for (int i = 0; i < this.numofPoints; ++i)
			{
				uv = R * uv;
				Vector2d p = this.center + uv;
				CurvePoint cp = new CurvePoint(p, i);
				this.curvePoints.Add(cp);
				if (i % K == 0)
					this.controlPoints.Add(cp);
			}
		}
		public void Update(Vector2d c, Vector2d xaxis, double radius)
		{
			this.center = c;
			this.xaxis = xaxis;
			this.radius = radius;

			double dtheta = Math.PI / this.numofPoints;
			double cos = Math.Cos(dtheta);
			double sin = Math.Sin(dtheta);
			Matrix2d R = new Matrix2d(cos, -sin, sin, cos);
			Vector2d uv = this.xaxis * radius;
			int K = this.numofPoints / 2;
			for (int i = 0; i < this.numofPoints; ++i)
			{
				uv = R * uv;
				this.curvePoints[i].pos = this.center + uv;
			}
		}
		public void Transform(Vector2d c, Vector2d xaxis, double radius)
		{
			Vector2d translate = c - this.center;
			Vector2d norm1 = this.xaxis.Normalize();
			Vector2d norm2 = xaxis.Normalize();
			double theta1 = GetTheta(norm1);
			double theta2 = GetTheta(norm2);
			double theta = theta2 - theta1;
			double cos = Math.Cos(theta);
			double sin = Math.Sin(theta);

			double scale = radius / this.radius;

			Matrix2d M = new Matrix2d(cos, -sin, sin, cos);
			Matrix2d S = new Matrix2d(scale, 0, 0, scale);
			Matrix2d Q = M * S;
			foreach (CurvePoint cp in this.curvePoints)
			{
				cp.pos = Q * cp.pos;
			}
			foreach (CurvePoint cp in this.curvePoints)
			{
				cp.pos += translate;
			}
			this.xaxis = xaxis;
		}
		private double GetTheta(Vector2d u)
		{
			double theta = Math.Acos(u.x);
			return u.y < 0 ? Math.PI * 2 - theta : theta;
		}
		public override void GetCurvePoints()
		{
			this.CreateDefault();
		}
		public override void Update(Vector2d displacement)
		{
			if (this.selectedControlPointIndex == -1) return;

			Vector2d oldpos = new Vector2d(this.selectedControlPointPos.x, this.selectedControlPointPos.y);
			Vector2d newpos = oldpos + displacement;
			Vector2d newaxis = (newpos - this.center);
			double newradius = newaxis.Length();
			newaxis.Normalize();
			if (this.selectedControlPointIndex == 1) // the mid control point
			{
				newaxis = new Matrix2d(0, 1, -1, 0) * newaxis;
			}

			this.Update(this.center, newaxis, newradius);
		}
	}
	
	public class Ellipse : Curve
	{
		public double xRadius, yRadius;
		public Vector2d center;
		public int numofPoints;
		public Ellipse(Vector2d c, double rx, double ry, int slices)
		{
			this.type = CurveStyle.ellipse;
			this.center = c;
			this.xRadius = rx;
			this.yRadius = ry;
			this.numofPoints = slices;
			this.isClosed = true;
			this.GetCurvePoints();
		}
		public override void GetCurvePoints()
		{
			this.curvePoints = new List<CurvePoint>();
			this.controlPoints = new List<CurvePoint>();
			double dtheta = Math.PI * 2 / this.numofPoints;
			int K = this.numofPoints / 4;
			for (int i = 0; i < this.numofPoints; ++i)
			{
				double theta = i * dtheta;
				double x = this.xRadius * Math.Cos(theta);
				double y = this.yRadius * Math.Sin(theta);
				Vector2d p = this.center + new Vector2d(x, y);
				CurvePoint cp = new CurvePoint(p, i);
				this.curvePoints.Add(cp);
				if (i % K == 0)
					this.controlPoints.Add(cp);
			}
		}
		public override void Update(Vector2d displacement)
		{
			if (this.selectedControlPointIndex == -1) return;

			Vector2d oldpos = new Vector2d(this.selectedControlPointPos.x, this.selectedControlPointPos.y);
			if (this.selectedControlPointIndex == 0 || this.selectedControlPointIndex == 2)
			{
				oldpos.x += displacement.x;
				this.xRadius = (oldpos - this.center).Length();
			}
			else
			{
				oldpos.y += displacement.y;
				this.yRadius = (oldpos - this.center).Length();
			}
			this.GetCurvePoints();
		}
	}
	
	public class Rect : Curve
	{
		public double xLen, yLen;
		public Vector2d center;
		public int numofPoints;
		public Rect(Vector2d c, double rx, double ry, int slices)
		{
			this.type = CurveStyle.rectangle;
			this.center = c;
			this.xLen = rx;
			this.yLen = ry;
			this.numofPoints = slices;
			this.isClosed = true;
			this.GetCurvePoints();
		}
		public override void GetCurvePoints()
		{
			this.curvePoints = new List<CurvePoint>();
			this.controlPoints = new List<CurvePoint>();
			
			int K = this.numofPoints / 4;

			double cx = this.center.x, cy = this.center.y;
			Vector2d[] B = new Vector2d[4]
			{
				new Vector2d(cx - xLen/2, cy-yLen/2),
				new Vector2d(cx + xLen/2, cy-yLen/2),
				new Vector2d(cx + xLen/2, cy+yLen/2),
				new Vector2d(cx - xLen/2, cy+yLen/2)
			};

			int index = 0;
			for (int i = 0; i < 4; ++i)
			{
				Vector2d p = B[i], q = B[(i + 1) % 4];
				double step = 1.0 / K;
				for (int j = 0; j < K; ++j)
				{
					double r = step * j;
					Vector2d t = p * (1 - r) + q * r;
					CurvePoint cp = new CurvePoint(t,index++);
					this.curvePoints.Add(cp);
					if (j == 0)
					{
						this.controlPoints.Add(cp);
					}
				}
			}

		}
		public override void Update(Vector2d displacement)
		{
			if (this.selectedControlPointIndex == -1) return;

			Vector2d oldpos = new Vector2d(this.selectedControlPointPos.x, this.selectedControlPointPos.y);
			oldpos += displacement;

			Vector2d oppose = this.controlPoints[(this.selectedControlPointIndex + 2) % 4].pos;
			this.center = (oppose + oldpos)/2;
			this.xLen = 2*Math.Abs((oldpos - this.center).x);
			this.yLen = 2*Math.Abs((oldpos - this.center).y);

			this.GetCurvePoints();
		}
	}
	
	public class PolyLine : Curve
	{
		public int numofPoints;
		public int subdivisionLevel = 8;
		public PolyLine(List<Vector2d> inputPoints)
		{
			this.type = CurveStyle.polyline;
			this.controlPoints = new List<CurvePoint>();
			foreach (Vector2d p in inputPoints)
			{
				this.controlPoints.Add(new CurvePoint(p));
			}
			this.GetCurvePoints();
		}
		public void SetSubdivisionLevel(int level)
		{
			this.subdivisionLevel = level;
		}
		public override void GetCurvePoints()
		{
			this.curvePoints = new List<CurvePoint>();

			int n = this.controlPoints.Count;
			int k = this.isClosed ? n : n - 1;
			for (int i = 0; i < k; ++i)
			{
				Vector2d p = this.controlPoints[i].pos, q = this.controlPoints[(i + 1)%n].pos;
				double step = 1.0 / this.subdivisionLevel;
				for (int j = 0; j < this.subdivisionLevel; ++j)
				{
					double r = step * j;
					Vector2d t = p * (1 - r) + q * r;
					this.curvePoints.Add(new CurvePoint(t,i));
				}
			}

			this.numofPoints = this.curvePoints.Count;

		}
		public override void Update(Vector2d displacement)
		{
			if (this.selectedControlPointIndex == -1) return;
			Vector2d oldpos = new Vector2d(this.selectedControlPointPos.x, this.selectedControlPointPos.y);
			oldpos += displacement;
			this.controlPoints[this.selectedControlPointIndex].pos = oldpos;
			this.GetCurvePoints();
		}
	}
	
	public class CubicSpline : Curve
	{
		public int numofPoints;
		public int subdivisionLevel = 8;
		public CubicSpline(List<Vector2d> ctrlPoints)
		{
			this.type = CurveStyle.cubicspline;
			this.controlPoints = new List<CurvePoint>();
			foreach (Vector2d p in ctrlPoints)
			{
				this.controlPoints.Add(new CurvePoint(p));
			}
			this.isClosed = false;
			this.GetCurvePoints();
		}
		public CubicSpline(List<Vector2d> ctrlPoints, int subdivLevel)
		{
			this.type = CurveStyle.cubicspline;
			this.subdivisionLevel = subdivLevel;
			this.controlPoints = new List<CurvePoint>();
			foreach (Vector2d p in ctrlPoints)
			{
				this.controlPoints.Add(new CurvePoint(p));
			}
			this.isClosed = false;
			this.GetCurvePoints();
		}
		public void SetSubdivisionLevel(int level)
		{
			this.subdivisionLevel = level;
		}
		public override void GetCurvePoints()
		{
			int n = this.controlPoints.Count;
			Vector2d[] points = new Vector2d[n];
			for (int i = 0; i < n; ++i)
			{
				points[i] = this.controlPoints[i].pos;
			}

			CubicSplineFitter2d cubicSplineFitter =
				new CubicSplineFitter2d(points);
			
			Vector2d[] outPoints = 
				cubicSplineFitter.GetCurvePoints(this.subdivisionLevel);

			this.curvePoints = new List<CurvePoint>();
			int index = 0;
			foreach (Vector2d p in outPoints)
			{
				this.curvePoints.Add(new CurvePoint(p, index++));
			}

			this.numofPoints = outPoints.Length;

		}
		public override void Update(Vector2d displacement)
		{
			if (this.selectedControlPointIndex == -1) return;
			Vector2d oldpos = new Vector2d(this.selectedControlPointPos.x, this.selectedControlPointPos.y);
			oldpos += displacement;
			this.controlPoints[this.selectedControlPointIndex].pos = oldpos;
			this.GetCurvePoints();
		}
	}

}
