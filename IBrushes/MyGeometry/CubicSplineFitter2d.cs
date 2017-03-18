using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using MyCholmodSolver;

namespace MyGeometry
{
	public unsafe class CubicSplineFitter2d
	{
		private int dimension = 2;		// default 2d
		private int numOfPoints = 0;
		private Vector2d[] inputPoints = null;
		private Vector2d[] B = null;	// bezier control inputPoints

		// interfaces
		public Vector2d[] GetCurvePoints(int numPerInterval)
		{
			List<Vector2d> curvePoints = new List<Vector2d>();
			int num = numOfPoints + numPerInterval * (numOfPoints - 1);
			double lamb = 1.0 / (numPerInterval + 1);
			for (int i = 0; i < numOfPoints - 1; ++i)
			{
				Vector2d b0 = inputPoints[i];
				Vector2d b1 = 2*B[i]/3.0 + B[i+1]/3.0;
				Vector2d b2 = B[i]/3.0 + 2*B[i+1]/3.0;;
				Vector2d b3 = inputPoints[i+1];

				for (int j = 0; j <= numPerInterval; ++j) // start point + inner inputPoints
				{
					double t = j * lamb; 
					double dt = 1.0 - t;
					double t0 = dt*dt*dt, t1 = 3*dt*dt*t, t2 = 3*dt*t*t, t3 = t*t*t;
					Vector2d pt = t0*b0 + t1*b1 + t2*b2 + t3*b3;

					curvePoints.Add(pt);
				}
			}
			Vector2d et = inputPoints[numOfPoints - 1];
			curvePoints.Add(
				new Vector2d(et.x, et.y)
			);

			return curvePoints.ToArray();
		}

		// constructors and destructors
		public CubicSplineFitter2d(Vector2d[] pts)
		{
			this.inputPoints = pts;
			this.numOfPoints = pts.Length;
			this.B = new Vector2d[numOfPoints];
			for (int i = 0; i < numOfPoints; ++i)
			{
				this.B[i] = new Vector2d();
			}
			this.ComputeBezierControlPoints();
		}
		~CubicSplineFitter2d()
		{

		}

		// main function
		private void ComputeBezierControlPoints()
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
					Vector2d pt;
					if (i == 0)
					{
						pt = r * inputPoints[1] - inputPoints[0];
					}
					else if (i == n - 1)
					{
						pt = r * inputPoints[i+1] - inputPoints[i+2];
					}
					else
					{
						pt = r * inputPoints[i + 1];
					}
					switch (k)
					{
						case 0:
							b[i] = pt.x;// (inputPoints[t].x - inputPoints[s].x);
							break;
						case 1:
							b[i] = pt.y;//(inputPoints[t].y - inputPoints[s].y);
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
					}
				}
			}
			B[0] = inputPoints[0];
			B[B.Length - 1] = inputPoints[inputPoints.Length - 1];
			
			cholmodSolver.Release();
		}
	}
}

