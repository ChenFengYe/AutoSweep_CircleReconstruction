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

            #region a theoretically fast method, but there are mystery bugs in it for now =.=
            K = this.polyList.Count;
            N = this.polyList[0].CornerPoints2d.Count;

            //// compute every triangle's initial height
            //h = new double[K * N]; // using Heron's formula to compute height
            //CLineSegment[] bottomEdges = new CLineSegment[K * N]; // bottom edges (used to calculate dis)
            //int triangleIdx = 0;
            int idx = 0;
            center2d = new MyVector2[K];
            foreach (MyPolygon mp in this.polyList)
            {
                center2d[idx] = new MyVector2();
                foreach (MyVector2 vert in mp.CornerPoints2d)
                    center2d[idx] += vert;
                center2d[idx] /= mp.CornerPoints2d.Count;

                //double edgeLengthA = (mp.CornerPoints2d[0] - center2d[idx]).Length();
                //for (int i = 0; i < mp.CornerPoints2d.Count; i++)
                //{
                //    bottomEdges[triangleIdx] = new CLineSegment(mp.CornerPoints2d[i], mp.CornerPoints2d[(i + 1) % N]);

                //    double edgeLengthB = (mp.CornerPoints2d[(i + 1) % N] - center2d).Length();
                //    double edgeLengthBottom = bottomEdges[triangleIdx].GetLineSegmentLength();
                    
                //    double p = (edgeLengthA + edgeLengthB + edgeLengthBottom) / 2;
                //    double S = Math.Sqrt(p * (p - edgeLengthA) * (p - edgeLengthB) * (p - edgeLengthBottom));
                //    h[triangleIdx] = 2 * S / edgeLengthBottom;
                    
                //    triangleIdx++;
                //    edgeLengthA = edgeLengthB;
                //}
                idx++;
            }

            //// compute initial distance between slice points and its nearest bottom edge
            //dis = new double[M_sum]; // initial distances
            //segmentIdx = new int[M_sum]; // belong to which bottom edge
            //int disIdx = 0, mpIdx = 0;
            //foreach (MyPolygon mp in this.polyList)
            //{
            //    foreach (MyVector2 mv2 in mp.SlicePoints2d)
            //    {
            //        double minDis = double.MaxValue;
            //        for (int i = mpIdx; i < mpIdx + N; i++)
            //        {
            //            // must be counterclockwise, and when point is on the left of the segment, dis > 0
            //            double tmpDis = bottomEdges[i].GetDistance(mv2) * -bottomEdges[i].GetPointLocation(mv2);
            //            if (Math.Abs(minDis) > Math.Abs(tmpDis))
            //            {
            //                minDis = tmpDis;
            //                segmentIdx[disIdx] = i;
            //            }
            //        }
            //        dis[disIdx++] = minDis;
            //    }
            //    mpIdx += N;
            //}

            //// saved rho (to compute deltaRho in each iteration)
            //rho = new double[K];
            //for (int i = 0; i < K; i++)
            //    rho[i] = 1.0;
            #endregion

            // begin to optimize
            double[] x = new double[K];
            for (int i = 0; i < K; i++)
                x[i] = 1.0;
            double epsg = 1e-4;
            double epsf = 0;
            double epsx = 0;
            int maxits = 100; //300
            alglib.minlmstate state;
            alglib.minlmreport rep;

            alglib.minlmcreatev(2, x, 0.01, out state);
            alglib.minlmsetcond(state, epsg, epsf, epsx, maxits);
            alglib.minlmoptimize(state, funcFitPolygon1, null, null);
            alglib.minlmresults(state, out x, out rep);

            idx = 0;
            foreach (MyPolygon mp in this.polyList)
            {
                MyVector3 center = mp.Center;
                for (int i = 0; i < mp.CornerPoints3d.Count; i++)
                    mp.CornerPoints3d[i] = x[idx] * (mp.CornerPoints3d[i] - center) + center;
                idx++;
            }
		}
        int K, N;
        //int M_sum;
        //double[] h;
        //double[] dis;
        //int[] segmentIdx;
        //double[] rho;

        //private void funcFitPolygon(double[] x, double[] f, object obj)
        //{
        //    double residual = 0.0;
        //    f[K] = 0.0;

        //    double[] deltaX = new double[K];
        //    for (int i = 0; i < K; i++)
        //        deltaX[i] = x[i] - rho[i];

        //    double[] deltaDis = new double[K * N];
        //    int idx = 0;
        //    for (int i = 0; i < K; i++)
        //    {
        //        for (int j = 0; j < N; j++)
        //        {
        //            deltaDis[idx] = h[idx] * deltaX[i] / rho[i];
        //            h[idx] += deltaDis[idx];
        //            idx++;
        //        }
        //    }

        //    //// update rho
        //    //for (int i = 0; i < K; i++)
        //    //    rho[i] = x[i]; // update rho

        //    // update dis
        //    idx = 0;
        //    for (int i = 0; i < K; i++)
        //    {
        //        f[i] = 0.0;
        //        for (int j = 0; j < this.polyList[i].SlicePoints2d.Count; j++)
        //        {
        //            dis[idx] += deltaDis[segmentIdx[idx]];
        //            f[i] += Math.Abs(dis[idx]);
        //            idx++;
        //        }
        //        residual += f[i];
        //    }

        //    //// compute f[1]
        //    //for (int i = 0; i < K; i++)
        //    //{
        //    //    if (i == 0)
        //    //        f[1] += Math.Pow(x[i] - x[i + 1], 2);
        //    //    else if (i == K - 1)
        //    //        f[1] += Math.Pow(x[i] - x[i - 1], 2);
        //    //    else
        //    //        f[1] += Math.Pow(2 * x[i] - x[i - 1] - x[i + 1], 2);
        //    //}
        //    //f[1] *= 1e-4;

        //    if (outputIdx == 0)
        //        Console.WriteLine(residual);
        //    outputIdx = (outputIdx + 1) % 1000;
        //}
        MyVector2[] center2d;
        int outputIdx = 0;
        private void funcFitPolygon1(double[] x, double[] f, object obj)
        {
            f[0] = 0.0;
            f[1] = 0.0;

            int idx = 0;
            foreach (MyPolygon mp in this.polyList)
            {
                MyVector2[] vertices = new MyVector2[N];
                for (int i = 0; i < mp.CornerPoints2d.Count; i++)
                    vertices[i] = center2d[idx] + (mp.CornerPoints2d[i] - center2d[idx]) * x[idx];

                CLine[] segments = new CLine[N];
                for (int i = 0; i < mp.CornerPoints2d.Count; i++)
                    segments[i] = new CLine(vertices[i], vertices[(i + 1) % N]);

                foreach(MyVector2 mv2 in mp.SlicePoints2d)
                {
                    double minDis = double.MaxValue;
                    for (int i = 0; i < N; i++)
                    {
                        double disTmp = segments[i].GetDistance(mv2);
                        if (minDis > disTmp)
                            minDis = disTmp;
                    }
                    f[0] += minDis;
                }
                idx++;
            }

            // compute f[1]
            for (int i = 0; i < K; i++)
            {
                if (i == 0)
                    f[1] += Math.Pow(x[i] - x[i + 1], 2);
                else if (i == K - 1)
                    f[1] += Math.Pow(x[i] - x[i - 1], 2);
                else
                    f[1] += Math.Pow(2 * x[i] - x[i - 1] - x[i + 1], 2);
            }
            f[1] *= 1e-1;

            if (outputIdx == 0)
                Console.WriteLine(f[0].ToString("0.000") + "\t" + f[1].ToString("0.000"));
            outputIdx = (outputIdx + 1) % 1000;
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
