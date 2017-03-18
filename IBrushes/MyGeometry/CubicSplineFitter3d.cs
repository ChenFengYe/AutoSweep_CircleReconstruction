using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using MyCholmodSolver;

namespace MyGeometry
{
	public unsafe class CubicSplineFitter3d : IDisposable
	{
		private int dimension = 3; // default 2d
		private int numOfPoints = 0;
		private Vector3d[] points = null;
		private Vector3d[] B = null; // bezier control points

		// constructors and destructors
		public CubicSplineFitter3d(Vector3d[] pts)
		{
			this.points = pts;
			this.numOfPoints = pts.Length;
			this.B = new Vector3d[numOfPoints];
			for (int i = 0; i < numOfPoints; ++i)
			{
				this.B[i] = new Vector3d();
			}
			this.ComputrBezierControlPoints();
		}
		~CubicSplineFitter3d()
		{

		}


		// interfaces
		public Vector3d[] GetPoints(int numPerInterval)
		{
			List<Vector3d> curvePoints = new List<Vector3d>();
			int num = numOfPoints + numPerInterval * (numOfPoints - 1);
			double lamb = 1.0 / (numPerInterval + 1);
			for (int i = 0; i < numOfPoints - 1; ++i)
			{
				Vector3d b0 = points[i];
				Vector3d b1 = 2*B[i]/3.0 + B[i+1]/3.0;
				Vector3d b2 = B[i]/3.0 + 2*B[i+1]/3.0;;
				Vector3d b3 = points[i+1];

				for (int j = 0; j <= numPerInterval; ++j) // start point + inner points
				{
					double t = j * lamb; 
					double dt = 1.0 - t;
					double t0 = dt*dt*dt, t1 = 3*dt*dt*t, t2 = 3*dt*t*t, t3 = t*t*t;
					Vector3d pt = t0*b0 + t1*b1 + t2*b2 + t3*b3;

					curvePoints.Add(pt);
				}
			}
			Vector3d et = points[numOfPoints - 1];
			curvePoints.Add( // last point
				new Vector3d(et.x, et.y, et.z)
			);

			return curvePoints.ToArray();
		}
		public Vector3d[] GetPoints2(int numPerInterval)
		{
			List<Vector3d> curvePoints = new List<Vector3d>();
			int num = numOfPoints + numPerInterval * (numOfPoints - 1);
			double lamb = 1.0 / (numPerInterval + 1);
			for (int i = 0; i < numOfPoints - 1; ++i)
			{
				Vector3d u = (points[i] - points[i + 1]) * 2 +
				B[i] + B[i + 1];
				Vector3d v = (points[i + 1] - points[i]) * 3 -
					(2 * B[i] + B[i + 1]);
				Vector3d w = B[i];
				Vector3d r = points[i];

				Vector4d p = new Vector4d(u.x, v.x, w.x, r.x);
				Vector4d q = new Vector4d(u.y, v.y, w.y, r.y);
				Vector4d o = new Vector4d(u.z, v.z, w.z, r.z);

				for (int j = 0; j <= numPerInterval; ++j) // start point + inner points
				{
					double t = j * lamb;
					double t0 = 1, t1 = t, t2 = t * t, t3 = t * t * t;
					Vector4d g = new Vector4d(t3, t2, t1, t0);
					double x = g.Dot(p);
					double y = g.Dot(q);
					double z = g.Dot(o);

					curvePoints.Add(
						new Vector3d(x, y, z)
					);
				}
			}
			Vector3d pt = points[numOfPoints - 1];
			curvePoints.Add( // last point
				new Vector3d(pt.x, pt.y, pt.z)
			);

			return curvePoints.ToArray();
		}
		public Vector3d GetPointAt(int i, double t)
		{
			Vector3d u = (points[i] - points[i + 1]) * 2 +
				B[i] + B[i + 1];
			Vector3d v = (points[i + 1] - points[i]) * 3 -
				(2 * B[i] + B[i + 1]);
			Vector3d w = B[i];
			Vector3d r = points[i];
			double t0 = 1, t1 = t, t2 = t * t, t3 = t * t * t;
			Vector4d g = new Vector4d(t3, t2, t1, t0);
			Vector4d p = new Vector4d(u.x, v.x, w.x, r.x);
			Vector4d q = new Vector4d(u.y, v.y, w.y, r.y);
			Vector4d o = new Vector4d(u.z, v.z, w.z, r.z);

			double x = g.Dot(p);
			double y = g.Dot(q);
			double z = g.Dot(o);

			return new Vector3d(x, y, z);
		}
		public Vector3d GetDerivAt(int i, double t)
		{
			Vector3d u = (points[i] - points[i + 1]) * 2 +
				B[i] + B[i + 1];
			Vector3d v = (points[i + 1] - points[i]) * 3 -
				(2 * B[i] + B[i + 1]);
			Vector3d w = B[i];

			double t0 = 1, t1 = t, t2 = t * t;
			Vector3d g = new Vector3d(t2, t1, t0);
			Vector3d p = new Vector3d(u.x, v.x, w.x);
			Vector3d q = new Vector3d(u.y, v.y, w.y);
			Vector3d o = new Vector3d(u.z, v.z, w.z);

			double x = g.Dot(p);
			double y = g.Dot(q);
			double z = g.Dot(o);

			return new Vector3d(x, y, z);
		}

		// members
		private void ComputrBezierControlPoints()
		{
			CholmodSolver cholmodSolver = new CholmodSolver();

			int n = numOfPoints-2;
			cholmodSolver.InitializeMatrixA(n, n);
			for (int i = 0; i < n; ++i)
			{
				if (i == 0)
				{
					cholmodSolver.Add_Coef(i, i, 4);
					cholmodSolver.Add_Coef(i, i + 1, 1);
				}
				else if (i == n - 1)
				{
					cholmodSolver.Add_Coef(i, i, 4);
					cholmodSolver.Add_Coef(i, i - 1, 1);
				}
				else
				{
					cholmodSolver.Add_Coef(i, i - 1, 1);
					cholmodSolver.Add_Coef(i, i, 4);
					cholmodSolver.Add_Coef(i, i + 1, 1);
				}
			}

			// factorize
			cholmodSolver.InitializeSolver();
			cholmodSolver.SetFinalPack(0);
			cholmodSolver.Factorize();

			// rhs
			double[] b = new double[n]; double[] x = new double[n];
			for (int k = 0; k < dimension; ++k)
			{
				for (int i = 0; i < n; ++i)
				{
					double r = 6.0;
					Vector3d pt;
					if (i == 0)
					{
						pt = r * points[1] - points[0];
					}
					else if (i == n - 1)
					{
						pt = r * points[i+1] - points[i+2];
					}
					else
					{
						pt = r * points[i + 1];
					}
					switch (k)
					{
						case 0:
							b[i] = pt.x;// (points[t].x - points[s].x);
							break;
						case 1:
							b[i] = pt.y;//(points[t].y - points[s].y);
							break;
						case 2:
							b[i] = pt.z;//(points[t].z - points[s].z);
							break;
					}
				}

				// solve
				fixed (double* _b = b)
				{
					cholmodSolver.InitializeMatrixB(_b, n, 1);
				}
				fixed (double* _x = x)
				{
					cholmodSolver.Linear_Solve(_x, false);
				}

				// assign
				for (int i = 0; i < n; ++i)
				{
					int tt = i + 1;
					switch (k)
					{
						case 0:
							this.B[tt].x = x[i];
							break;
						case 1:
							this.B[tt].y = x[i];
							break;
						case 2:
							this.B[tt].z = x[i];
							break;
					}
				}
			}
			B[0] = points[0];
			B[B.Length - 1] = points[points.Length - 1];
			// release solver
			cholmodSolver.Release();
		}

		#region IDisposable Members
		public void Dispose()
		{
		}
		#endregion
	}
}
