using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using MyGeometry;
using MyCholmodSolver;

namespace SmartCanvas
{
	public unsafe class CubicSpline2
	{
		private int dimension = 2;		// default 2d
		private int numOfPoints = 0;
		private MyVector2[] inputPoints = null;
		private MyVector2[] B = null;	// bezier control inputPoints

		public List<MyVector2> sampledPoints_ = null;
		public List<MyVector2> sampledPointsTangents_ = null;
		public List<MyVector2> sampledPointsNormals_ = null;
		public List<MyVector2> jointPoints_ = new List<MyVector2>();
		public void SamplePoints(int numpoints)
		{
			// handle degenerated case
			if (numpoints < 2) {
				Console.WriteLine(@" < 2 sample points");
				this.sampledPoints_ = this.GetCurvePoints(1);
				this.sampledPointsTangents_ = this.GetCurvePointsTagent(1);
				this.GetSamplePointsNormals();
				return;
			}

			// sample uniformly distributed points on the curves (#of numofPoints)
			this.sampledPoints_ = new List<MyVector2>();
			this.sampledPointsTangents_ = new List<MyVector2>();


			// get a temp rough curve length
			double tmpLength = 0;
			int K = 5;
			List<MyVector2> tmppoints = this.GetCurvePoints(5);
			for (int i = 0; i < tmppoints.Count-1; ++i)
			{
				tmpLength += (tmppoints[i + 1] - tmppoints[i]).Length();
			}
			
			// resampling
			double len = tmpLength / (numpoints-1);
			int L = 10;
			int[] numi = new int[this.numOfPoints-1];
			int count = 0;
			for (int i = 0; i < tmppoints.Count-1; i+=(K+1))
			{
				double intv_len = 0;
				for (int j = 0; j <= K; ++j)
				{
					intv_len += (tmppoints[i + j + 1] - tmppoints[i + j]).Length();
				}
				numi[count++] = (int)((intv_len / len) * L);
			}
			
			// assign & get samples
			List<MyVector2> resampledPoints = this.GetCurvePoints(numi);
			List<MyVector2> resampledPointsTagent = this.GetCurvePointsTagent(numi);
			tmpLength = 0;
			for (int i = 0; i < resampledPoints.Count - 1; ++i)
			{
				tmpLength += (resampledPoints[i + 1] - resampledPoints[i]).Length();
			}

			// get samples
			double step_len = tmpLength / (numpoints - 1);
			double curr_len = 0;
			int indx = 0;
			for (int i = 0; i < resampledPoints.Count-1; ++i)
			{
				curr_len += (resampledPoints[i + 1] - resampledPoints[i]).Length();
				if (curr_len > indx * step_len)
				{
					this.sampledPoints_.Add(resampledPoints[i]);
					this.sampledPointsTangents_.Add(resampledPointsTagent[i]);
					indx++;
				}
			}
			
			// the last point
			this.sampledPoints_.Add(resampledPoints[resampledPoints.Count - 1]);
			this.sampledPointsTangents_.Add(resampledPointsTagent[resampledPointsTagent.Count - 1]);

			this.GetSamplePointsNormals();
		}
		public void SamplePoints(double dis_thres)
		{
			// get sample points w.r.t. the distance threshold
			// a ||s_i-s_i+1|| < dis_thres is required
			this.sampledPoints_ = new List<MyVector2>();
			this.sampledPointsTangents_ = new List<MyVector2>();

			for (int i = 0; i < this.numOfPoints-1; ++i)
			{
				double len = (this.inputPoints[i+1] - this.inputPoints[i]).Length();
				int npointsPerSegment = (int)(len / dis_thres);
				bool stop = false;
				while (!stop)
				{
					stop = true;

					List<MyVector2> points = new List<MyVector2>();
					List<MyVector2> pointsTagent = new List<MyVector2>();
					for (int j = 0; j < npointsPerSegment; ++j)
					{
						double t1 = (double)j / npointsPerSegment;
						double t2 = (double)(j + 1) / npointsPerSegment;
						MyVector2 p1 = this.GetCurvePoint(i, t1);
						MyVector2 p2 = this.GetCurvePoint(i, t2);
						
						if ((p2 - p1).Length() > dis_thres)
						{
							stop = false;
							npointsPerSegment++;
							break;
						}
						
						MyVector2 pt1 = this.GetCurvePointTagent(i, t1);
						points.Add(p1);
						pointsTagent.Add(pt1);
					}

					if (stop == true)
					{
						this.sampledPoints_.AddRange(points);
						this.sampledPointsTangents_.AddRange(pointsTagent);
					}
				}
			}

			// end point
			this.sampledPoints_.Add(this.GetCurvePoint(this.numOfPoints - 2, 1));
			this.sampledPointsTangents_.Add(this.GetCurvePointTagent(this.numOfPoints - 2, 1));

			// get normals
			this.GetSamplePointsNormals();

		}
		public void GetSamplePointsNormals()
		{
			if (this.sampledPointsTangents_ == null) return;

			this.sampledPointsNormals_ = new List<MyVector2>();
			foreach (MyVector2 tagent in this.sampledPointsTangents_)
			{
				MyVector2 normal = new MyVector2(tagent.y, -tagent.x);
				this.sampledPointsNormals_.Add(normal);
			}
		}

		// interfaces
		public List<MyVector2> GetCurvePoints(int numPerInterval)
		{
			// numPerInterval -- number of points lies in-between two curve points in each interval
			List<MyVector2> curvePoints = new List<MyVector2>();
			int num = numOfPoints + numPerInterval * (numOfPoints - 1);
			double lamb = 1.0 / (numPerInterval + 1);
			for (int i = 0; i < numOfPoints - 1; ++i)
			{
				for (int j = 0; j <= numPerInterval; ++j) // start point + inner inputPoints
				{
					double t = j * lamb; 
					curvePoints.Add(this.GetCurvePoint(i,t));
				}
			}

			// last point
			MyVector2 et = inputPoints[numOfPoints - 1];
			curvePoints.Add(
				new MyVector2(et.x, et.y)
			);

			return curvePoints;
		}
		public List<MyVector2> GetCurvePointsTagent(int numPerInterval)
		{
			// numPerInterval -- number of points lies in-between two curve points in each interval
			List<MyVector2> curvePointsTagent = new List<MyVector2>();
			int num = numOfPoints + numPerInterval * (numOfPoints - 1);
			double lamb = 1.0 / (numPerInterval + 1);
			for (int i = 0; i < numOfPoints - 1; ++i)
			{
				for (int j = 0; j <= numPerInterval; ++j) // start point + inner inputPoints
				{
					double t = j * lamb;
					curvePointsTagent.Add(this.GetCurvePointTagent(i, t));
				}
			}
			
			// last point
			MyVector2 et = this.GetCurvePointTagent(numOfPoints-2, 1);
			curvePointsTagent.Add(et);
			
			return curvePointsTagent;
		}
		public List<MyVector2> GetCurvePoints(int[] numPerInterval)
		{
			List<MyVector2> curvePoints = new List<MyVector2>();
			for (int i = 0; i < numOfPoints - 1; ++i)
			{
				double lamb = 1.0 / (numPerInterval[i] + 1);
				for (int j = 0; j <= numPerInterval[i]; ++j) // start point + inner inputPoints
				{
					double t = j * lamb;
					curvePoints.Add(this.GetCurvePoint(i, t));
				}
			}

			// last point
			MyVector2 et = inputPoints[numOfPoints - 1];
			curvePoints.Add(
				new MyVector2(et.x, et.y)
			);

			return curvePoints;
		}
		public List<MyVector2> GetCurvePointsTagent(int[] numPerInterval)
		{
			List<MyVector2> curvePointsTagent = new List<MyVector2>();
			for (int i = 0; i < numOfPoints - 1; ++i)
			{
				double lamb = 1.0 / (numPerInterval[i] + 1);
				for (int j = 0; j <= numPerInterval[i]; ++j) // start point + inner inputPoints
				{
					double t = j * lamb;
					curvePointsTagent.Add(this.GetCurvePointTagent(i, t));
				}
			}
			
			// last point
			MyVector2 et = this.GetCurvePointTagent(numOfPoints - 2, 1);
			curvePointsTagent.Add(et);

			return curvePointsTagent;
		}
		private MyVector2 GetCurvePoint(int i, double t)
		{
			MyVector2 b0 = inputPoints[i];
			MyVector2 b1 = 2 * B[i] / 3.0 + B[i + 1] / 3.0;
			MyVector2 b2 = B[i] / 3.0 + 2 * B[i + 1] / 3.0;
			MyVector2 b3 = inputPoints[i + 1];
			double dt = 1.0 - t;
			double t0 = dt * dt * dt, t1 = 3 * dt * dt * t, t2 = 3 * dt * t * t, t3 = t * t * t;
			return t0 * b0 + t1 * b1 + t2 * b2 + t3 * b3;
		}
		private MyVector2 GetCurvePointTagent(int i, double t)
		{
			MyVector2 b0 = inputPoints[i];
			MyVector2 b1 = 2 * B[i] / 3.0 + B[i + 1] / 3.0;
			MyVector2 b2 = B[i] / 3.0 + 2 * B[i + 1] / 3.0; ;
			MyVector2 b3 = inputPoints[i + 1];
			double dt = 1.0 - t;
			double t0 = - dt * dt * 3, t1 = 3 * (-2 * dt * t + dt * dt), t2 = 3 * (-t * t + 2*t*dt), t3 = 3 * t * t;
			return t0 * b0 + t1 * b1 + t2 * b2 + t3 * b3;
		}

		// constructors and destructors
		public CubicSpline2(MyVector2[] pts)
		{
			if (pts.Length < 3)
			{
				Console.WriteLine("err: must have > 3 points to interpolate!");
			}

			this.inputPoints = pts;
			this.numOfPoints = pts.Length;
			this.B = new MyVector2[numOfPoints];
			for (int i = 0; i < numOfPoints; ++i)
			{
				this.B[i] = new MyVector2();
			}
			this.ComputeBezierControlPoints();
		}
		~CubicSpline2()
		{

		}

		// main function
		private void ComputeBezierControlPoints()
		{
			int n = this.numOfPoints-2;
            
            double[,] a = new double[n,n];
			for (int i = 0; i < n; ++i)
			{
				if (i == 0)
				{
					a[i, i]     = 4;
					a[i, i + 1] = 1;
				}
				else if (i == n - 1)
				{
                    a[i, i]     = 4;
					a[i, i - 1] = 1;
				}
				else
				{
                    a[i, i]     = 4;
					a[i, i - 1] = 1;
                    a[i, i + 1] = 1;
				}
			}

            int info;
            alglib.densesolverreport rep;

            double[] b = new double[n]; double[] x = new double[n];
			// rhs
			for (int k = 0; k < dimension; ++k)
			{
				for (int i = 0; i < n; ++i)
				{
					double r = 6.0;
					MyVector2 pt;
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

                alglib.rmatrixsolve(a, n, b, out info, out rep, out x);

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
			
		}
	}
}

