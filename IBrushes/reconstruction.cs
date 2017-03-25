using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;

using MyGeometry;
using System.Drawing;
using System.Windows.Forms;

using PolygonDecompose;
using OpenTK.Graphics;
using OpenTK.Graphics.OpenGL;

using Emgu.Util;
using Emgu.CV;
using Emgu.CV.CvEnum;
using Emgu.CV.Structure;
using Emgu.CV.UI;
using Emgu.CV.Util;
using Emgu.CV.Features2D;

using Accord;
using Accord.Math;
using Accord.Math.Geometry;
using Accord.MachineLearning;
using Accord.MachineLearning.Geometry;
using Accord.Statistics.Distributions.DensityKernels;

using NumericalRecipes;


namespace SmartCanvas
{
    public partial class CanvasEngine
    {
        private Image<Gray, byte> mark;                                 // Mark of Image
        List<MyVector2> boundaryPoints_2d = new List<MyVector2>();      // The boundary of the Mark of the Image
        MyCircle topCircle = new MyCircle();                            // The Top Circle of the cylinder

        bool drawProjectedPoint = false;                                // To Draw the Project Points
        double offset = 0;                                              // Scale between 2D and 3D

        // Optima Param
        private List<MyVector2> CirclePoints_2d = new List<MyVector2>();// The Circle Sample Points projected on the screen
        List<List<double>> DistanceMap = new List<List<double>>();      // Distance Map
        double[,] DisMap;                                               // An array for distance map
        List<int> CorMap = new List<int>();                             // The correspondence of the BoundaryPoints and CirclePoints
        int Inter_Num = -1;                                             // The interaction time of optima
        int Inter_DMap = 1;                                             // While this time Update DistanceMap
        MyVector2 center_xy = new MyVector2();
        double center_z = 0;

        public void EstimatePlane()
        {
            mark = GetMarkImgae(this.Canvas);
            DisMap = new double[mark.Height, mark.Width];
            boundaryPoints_2d = GetBoundaryPoints(mark);
            //DistanceMap = BuildDistanceMap(mark.Rows, mark.Cols, boundaryPoints_2d);
            //drawProjectedPoint = true;

            // Init Optima Params
            center_z = 1.25;
            center_xy = GetMarkCenter(boundaryPoints_2d, center_z);
            //center_xy = GetMarkCenter2dPixel(boundaryPoints_2d, center_z);


            //double[] bndl = new double[] { -1, -1, -1, 0.02, -5, -5, 0 };
            //double[] x = new double[] { 0, -0.707, -0.707, 0.5, 0, 0, 1 }; //add center
            //double[] bndu = new double[] { 1, 1, 1, 5, 5, 5, 20 };

            ////double[] bndl = new double[] { -1, -1, -1, 0.02 };
            ////double[] x = new double[] { 0, -0.707, -0.707, 0.5 }; //add center
            ////double[] bndu = new double[] { 1, 1, 1, 5 };

            ////double[] scale = new double[] { 1, 1, 1, 1, 10, 10, 10 };
            //int paramsNum = 4;
            //int funNum = 3;

            double[] bndl = new double[] { -1, -1, -1, 0.02 };
            double[] x = new double[] { 0, 0, 0, 0.5 }; //add center
            double[] bndu = new double[] { 1, 1, 1, 5 };

            //double[] scale = new double[] { 1, 1, 1, 1, 10, 10, 10 };
            int paramsNum = 4;
            int funNum = 2;


            double diffstep = 0.00001;
            double epsg = 0.000000000001;
            double epsf = 0;
            double epsx = 0;
            int maxits = 0;
            alglib.minlmstate state;
            alglib.minlmreport rep;


            //set timer
            System.Diagnostics.Stopwatch stopwatch = new System.Diagnostics.Stopwatch();
            stopwatch.Start();

            // Do the optima
            alglib.minlmcreatev(funNum, x, diffstep, out state);
            alglib.minlmsetbc(state, bndl, bndu);
            alglib.minlmsetcond(state, epsg, epsf, epsx, maxits);
            //alglib.minlmsetscale(state, scale);
            alglib.minlmoptimize(state, function_project, null, null);
            alglib.minlmresults(state, out x, out rep);


            stopwatch.Stop();
            Console.WriteLine("Stop Type: {0}, Total time: {1}s", rep.terminationtype, stopwatch.ElapsedMilliseconds / 1000.0);

            // Update Circle
            //MyVector3 center = new MyVector3(x[4], x[5], x[6]);
            MyVector3 center = new MyVector3(center_xy.x, center_xy.y, center_z);


            MyVector3 normal = new MyVector3(FromRotationToNormal(x[0], x[1], x[2])); normal = normal.Normalize();
            double radius = x[3];

            topCircle = new MyCircle(center, radius, new MyPlane(center, normal));
            CirclePoints_2d = GetProjectionPoints_2D(topCircle.CirclePoints);

            ReOptimize();
            this.view.Refresh();
            topCircle.Save();
        }
        private void function_project(double[] x, double[] fi, object obj)
        {
            // Step 1: Init Params
            Inter_Num++;
            MyVector3 center = new MyVector3(center_xy.x, center_xy.y, center_z);
            //MyVector3 center = new MyVector3(x[4], x[5], x[6]);

            MyVector3 normal = new MyVector3(FromRotationToNormal(x[0], x[1], x[2])); normal = normal.Normalize();

            double radius = x[3];
            MyCircle myC = new MyCircle(center, radius, new MyPlane(center, normal), this.boundaryPoints_2d.Count());

            // Step 2: Do projection
            CirclePoints_2d = GetProjectionPoints_2D(myC.CirclePoints);

            // Step 3: Build DistanceMap
            if ((Inter_Num % Inter_DMap) == 0)
                CorMap = BuildCorrespondingMap(CirclePoints_2d, boundaryPoints_2d);

            // Step 4: Do the evaluation
            double dist_all = 0;
            double dist_max = 0;
            for (int i = 0; i < CirclePoints_2d.Count; i++)
            {
                //------------------------------------------------------------------------------------
                //// Cor Map
                //int corIndex = CorMap[i];
                //double dist1 = (CirclePoints_2d[i] - boundaryPoints_2d[corIndex]).SquareLength();

                //------------------------------------------------------------------------------------
                // Distance Map
                //int i_x = (int)Math.Round(CirclePoints_2d[i].x);
                //int i_y = (int)Math.Round(CirclePoints_2d[i].y);
                //i_x = Math.Min(i_x, mark.Cols - 1);
                //i_x = Math.Max(0, i_x);
                //i_y = Math.Min(i_y, mark.Rows - 1);
                //i_y = Math.Max(0, i_y);
                ////double dist2 = DistanceMap[i_y][i_x];

                //------------------------------------------------------------------------------------
                // Distance Interpolation from distance map
                //double dist2 = this.InterpolateDistanceWithDisMap(CirclePoints_2d[i]);

                //------------------------------------------------------------------------------------
                ////Do not use distance map, not accurate enough!
                double dist = double.MaxValue;
                for (int j = 0; j < boundaryPoints_2d.Count; j++)
                {
                    dist = Math.Min(dist, (CirclePoints_2d[i] - boundaryPoints_2d[j]).SquareLength());
                }


                //Console.WriteLine("2: {0}", Math.Abs(dist2 - dist));
                //Console.WriteLine("1: {0}", Math.Abs(dist1 - dist));

                dist_all += dist;
                //dist_max = Math.Max(dist, dist_max);
            }

            // Set Cost function
            fi[0] = dist_all;

            //////center 3d match 2d
            //double centerdis = (this.Compute2D(center) - center_xy).SquareLength();
            //fi[1] = centerdis;
            fi[1] = normal.SquareLength() - 1;

            if (Inter_Num % 50 == 0)
                System.Console.WriteLine("{0}|| N: {1},{2},{3} r: {4} cost:{5}",
                    Inter_Num, x[0], x[1], x[2], x[3], dist_all);

            // System.Console.WriteLine("{0}",Inter_Num);
        }

        public void ReOptimize()
        {
            double[] bndl = new double[] { -5, -5, 0 };
            double[] x = new double[] { topCircle.Center.x, topCircle.Center.y, topCircle.Center.z }; //add center
            double[] bndu = new double[] { 5, 5, 20 };

            int paramsNum = 3;
            int funNum = 1;

            double diffstep = 0.00001;
            double epsg = 0.000000000001;
            double epsf = 0;
            double epsx = 0;
            int maxits = 0;
            alglib.minlmstate state;
            alglib.minlmreport rep;

            //set timer
            System.Diagnostics.Stopwatch stopwatch = new System.Diagnostics.Stopwatch();
            stopwatch.Start();

            // Do the optima
            alglib.minlmcreatev(funNum, x, diffstep, out state);
            alglib.minlmsetbc(state, bndl, bndu);
            alglib.minlmsetcond(state, epsg, epsf, epsx, maxits);
            alglib.minlmoptimize(state, function_withcenter, null, null);
            alglib.minlmresults(state, out x, out rep);
            stopwatch.Stop();
            Console.WriteLine("Stop Type: {0}, Total time: {1}s", rep.terminationtype, stopwatch.ElapsedMilliseconds / 1000.0);

            // Update Circle
            MyVector3 center = new MyVector3(x[0], x[1], x[2]);
            MyVector3 normal = topCircle.Normal;
            double radius = topCircle.Radius;

            topCircle = new MyCircle(center, radius, new MyPlane(center, normal));
            CirclePoints_2d = GetProjectionPoints_2D(topCircle.CirclePoints);
            this.view.Refresh();
        }
        private void function_withcenter(double[] x, double[] fi, object obj)
        {
            Inter_Num++;
            MyVector3 center = new MyVector3(x[0], x[1], x[2]);
            MyVector3 normal = topCircle.Normal;
            double radius = topCircle.Radius;
            MyCircle myC = new MyCircle(center, radius, new MyPlane(center, normal), this.boundaryPoints_2d.Count());

            CirclePoints_2d = GetProjectionPoints_2D(myC.CirclePoints);
            double dist_all = 0;
            for (int i = 0; i < CirclePoints_2d.Count; i++)
            {
                double dist = double.MaxValue;
                for (int j = 0; j < boundaryPoints_2d.Count; j++)
                {
                    dist = Math.Min(dist, (CirclePoints_2d[i] - boundaryPoints_2d[j]).SquareLength());
                }
                dist_all += dist;
            }

            // Set Cost function
            fi[0] = dist_all;

            if (Inter_Num % 50 == 0)
                System.Console.WriteLine("{0}|| C: {1},{2},{3} cost:{4}",
                    Inter_Num, x[0], x[1], x[2], dist_all);
        }

        private MyVector3 FromRotationToNormal(double thetaX, double thetaY, double thetaZ)
        {
            MyMatrix4d RotationX = MyMatrix4d.RotationMatrix(new MyVector3(1, 0, 0), thetaX * Math.PI);
            MyMatrix4d RotationY = MyMatrix4d.RotationMatrix(new MyVector3(0, 1, 0), thetaY * Math.PI);
            MyMatrix4d RotationZ = MyMatrix4d.RotationMatrix(new MyVector3(0, 0, 1), thetaZ * Math.PI);
            MyVector3 normal_init = new MyVector3(0, -1, 0);
            MyVector3 normal_cur = (RotationZ * RotationY * RotationX * normal_init.ToMyVector4()).XYZ();
            return normal_cur.Normalize();
        }
        private double GetOtherNormalParams(double n0, double n1)
        {
            return Math.Sqrt(1 - n0 * n0 - n1 * n1);
        }

        private MyVector2 GetMarkCenter(List<MyVector2> bPoints, double CircleDistance)
        {
            offset = 1.0 / (Compute2D(new MyVector3(1, 0, CircleDistance)) - Compute2D(new MyVector3(0, 0, CircleDistance))).x;
            MyVector2 pic_c = Compute2D(new MyVector3(0, 0, CircleDistance));
            MyVector2 c = new MyVector2(0, 0);
            for (int i = 0; i < bPoints.Count; i++)
            {
                c += bPoints[i];
            }
            c = (c / bPoints.Count - pic_c) * offset;
            return c;
        }

        private MyVector2 GetMarkCenter2dPixel(List<MyVector2> bPoints, double CircleDistance)
        {
            offset = 1.0 / (Compute2D(new MyVector3(1, 0, CircleDistance)) - Compute2D(new MyVector3(0, 0, CircleDistance))).x;
            MyVector2 pic_c = Compute2D(new MyVector3(0, 0, CircleDistance));
            MyVector2 c = new MyVector2(0, 0);
            for (int i = 0; i < bPoints.Count; i++)
            {
                c += bPoints[i];
            }
            c = c / bPoints.Count;
            return c;
        }

        private List<List<double>> BuildDistanceMap(int rows, int cols, List<MyVector2> BoundaryPoints)
        {

            Image<Gray, byte> disimg = new Image<Gray, byte>(mark.Width, mark.Height, new Gray(0));

            List<List<double>> Dmap = new List<List<double>>();
            for (int i = 0; i < rows; i++)
            {
                List<double> rowDmap = new List<double>();
                for (int j = 0; j < cols; j++)
                {
                    double dist = double.MaxValue;
                    for (int k = 0; k < BoundaryPoints.Count; k++)
                    {
                        dist = Math.Min(dist, (new MyVector2(j, i) - BoundaryPoints[k]).SquareLength());
                    }
                    rowDmap.Add(dist);
                    disimg[i, j] = new Gray(dist);
                    this.DisMap[i, j] = dist;
                }
                Dmap.Add(rowDmap);
            }
            new ImageViewer(disimg, "dis map").Show();
            return Dmap;
        }

        private List<int> BuildCorrespondingMap(List<MyVector2> PrePoints, List<MyVector2> RealPoints)
        {
            List<int> DMap = new List<int>();
            for (int i = 0; i < PrePoints.Count; i++)
            {
                double dist = double.MaxValue;
                int index = 0;
                for (int j = 0; j < RealPoints.Count; j++)
                {
                    double tempdist = Math.Min(dist, Math.Pow((CirclePoints_2d[i] - boundaryPoints_2d[j]).Length(), 2));
                    if (dist > tempdist)
                    {
                        index = j;
                        dist = tempdist;
                    }
                }
                DMap.Add(index);
            }
            return DMap;
        }

        public List<MyVector2> GetBoundaryPoints(Image<Gray, byte> markImg, double cannyThreshold = 60, double cannyThresholdLinking = 100)
        {
            Image<Gray, Byte> cannyimg = markImg.Canny(cannyThreshold, cannyThresholdLinking);
            new ImageViewer(cannyimg, "boundary").Show();
            double backGround_threshold = 50;
            List<MyVector2> b_points = new List<MyVector2>();
            for (int i = 0; i < cannyimg.Height; i++)
            {
                for (int j = 0; j < cannyimg.Width; j++)
                {
                    if (cannyimg[i, j].Intensity > backGround_threshold)
                        b_points.Add(new MyVector2(j, i));
                }
            }
            return b_points;
        }

        public Image<Gray, byte> GetMarkImgae(Image<Bgr, byte> img)
        {
            Image<Gray, byte> mark = new Image<Gray, byte>(img.Width, img.Height, new Gray(0));
            mark.SetZero();

            double backGround_threshold = 10;
            for (int i = 0; i < mark.Height; i++)
            {
                for (int j = 0; j < mark.Width; j++)
                {
                    double grayValue =
                        (img[i, j].Red + img[i, j].Green + img[i, j].Blue) / 3.0;
                    if (grayValue > backGround_threshold)
                        mark[i, j] = new Gray(255);
                }
            }
            return mark;
        }

        public List<MyVector2> GetProjectionPoints_2D(List<MyVector3> points_3d)
        {
            List<MyVector2> points_2d = new List<MyVector2>();
            for (int i = 0; i < points_3d.Count; i++)
            {
                points_2d.Add(Compute2D(points_3d[i]));
            }
            return points_2d;
        }

        public double InterpolateDistanceWithDisMap(MyVector2 point_2d)
        {
            int lowerx = (int)Math.Floor(point_2d.x);
            int upperx = (int)Math.Ceiling(point_2d.x);
            int lowery = (int)Math.Floor(point_2d.y);
            int uppery = (int)Math.Ceiling(point_2d.y);

            lowerx = Math.Min(lowerx, mark.Cols - 1);
            lowerx = Math.Max(0, lowerx);
            upperx = Math.Min(upperx, mark.Cols - 1);
            upperx = Math.Max(0, upperx);

            lowery = Math.Min(lowery, mark.Rows - 1);
            lowery = Math.Max(0, lowery);
            uppery = Math.Min(uppery, mark.Rows - 1);
            uppery = Math.Max(0, uppery);


            double distance11 = DisMap[lowery, lowerx];
            double distance12 = DisMap[lowery, upperx];
            double distance21 = DisMap[uppery, lowerx];
            double distance22 = DisMap[uppery, upperx];

            if (lowerx == upperx && lowery == uppery)
                return distance11;

            //interpolation
            if (lowerx == upperx)
                return distance21 * (point_2d.y - lowery) / (uppery - lowery) + distance11 * (uppery - point_2d.y) / (uppery - lowery);
            if (lowery == uppery)
                return distance12 * (point_2d.x - lowerx) / (upperx - lowerx) + distance11 * (upperx - point_2d.x) / (upperx - lowerx);


            //bilateral interpolation
            double distance1 = (point_2d.x - lowerx) * distance12 / (upperx - lowerx) + (upperx - point_2d.x) * distance11 / (upperx - lowerx);
            double distance2 = (point_2d.x - lowerx) * distance22 / (upperx - lowerx) + (upperx - point_2d.x) * distance21 / (upperx - lowerx);
            double distance = (point_2d.y - lowery) * distance2 / (uppery - lowery) + (uppery - point_2d.y) * distance1 / (uppery - lowery);
            return distance;
        }

        public void ReadTopCircle(string filename)
        {
            topCircle = new MyCircle(filename);
        }

        //---------------------------------------------------------------------------
        //sweep----------------------------------------------------------------------
        //---------------------------------------------------------------------------
        public SweepMesh body = null;
        public MyPlane targetplane = null;
        public List<MyVector3> trajpoints_3d = new List<MyVector3>();
        public List<MyVector2> trajpoints = new List<MyVector2>();
        public List<MyVector2> trajpoints2 = new List<MyVector2>();
        public void Sweep(Image<Gray, byte> trajectoryimg_)
        {

            for (int i = 0; i < trajectoryimg_.Height; i++)
            {
                for (int j = 0; j < trajectoryimg_.Width; j++)
                {
                    if (trajectoryimg_[i, j].Intensity > 80)
                        trajpoints.Add(new MyVector2(j, i));
                }
            }

            trajpoints = CurveFitting(trajpoints);
            trajpoints = this.ResetPath(trajpoints, 10);

            trajpoints_3d = this.ProjectTrajAccording2Profile(trajpoints);

            body = new SweepMesh(this.topCircle, trajpoints_3d, null);


        }

        public List<MyVector2> CurveFitting(List<MyVector2> points)
        {
            Console.WriteLine("Fit curve");
            double[] x = new double[points.Count];
            double[] y = new double[points.Count];

            for (int i = 0; i < points.Count; i++)
            {
                x[i] = points[i].x;
                y[i] = points[i].y;
            }

            int outinfo = 0;
            alglib.spline1dfitreport rep;
            alglib.spline1dinterpolant p;
            alglib.spline1dfitpenalized(x, y, points.Count, 2.0, out outinfo, out p, out rep);

            List<MyVector2> output = new List<MyVector2>();
            for (int i = 0; i < p.innerobj.n; i++)
            {
                output.Add(new MyVector2(p.innerobj.x[i], alglib.spline1dcalc(p, p.innerobj.x[i])));
            }
            return output;

        }

        //---------------------------------------------------------------------------
        //sweep  with  Curve --------------------------------------------------
        //---------------------------------------------------------------------------
        public Image<Bgr, byte> edgeImage;  // Edge Image
        List<MyVector3> boundary3 = new List<MyVector3>();
        MyVector3 cur_p = new MyVector3();
        MyVector3 cur_dire = new MyVector3();
        MyPlane cutPlane = new MyPlane();
        MyVector3 Insection1 = new MyVector3();
        MyVector3 Insection2 = new MyVector3();
        public void CylinderSnapping()
        {
            // Get Boundary2
            List<MyVector2> boundary2 = new List<MyVector2>();
            for (int i = 0; i < edgeImage.Height; i++)
            {
                for (int j = 0; j < edgeImage.Width; j++)
                {
                    if (edgeImage[i, j].Red > 200)
                    {
                        boundary2.Add(new MyVector2(j, i));
                    }
                }
            }

            // Project  2D edge points
            MyVector3 normal = topCircle.Normal.Cross(this.camera.target).Cross(topCircle.Normal);
            MyPlane sectionPlane = new MyPlane(topCircle.Center, normal);
            boundary3 = Proj2dToPlane(sectionPlane, boundary2);

            double offset = 0.0001;
            cur_p = topCircle.Center - offset * topCircle.Normal;
            cur_dire = -1.0 * topCircle.Normal;

            // Step 1: Get IntersectionPoint

            cutPlane = new MyPlane(cur_p, sectionPlane.Normal().Cross(cur_dire).Normalize());
            int norInsec = -1;
            int notNorInsec = -1;
            RayTracein3DPlane(boundary3, new Line3(cur_p, cutPlane.Normal()), cutPlane, cur_p, out norInsec, out notNorInsec);
            Insection1 = boundary3[norInsec];
            Insection2 = boundary3[notNorInsec];
            this.view.Refresh();
        }

        //private MyVector3 GetLocalTangential(int p, List<MyVector3> boundary3)
        //{
        //    // Get Nearnest Pointss
        //    List<NearPoint> nearPs = new List<NearPoint>();
        //    for (int i = 0; i < boundary3.Count; i++)
        //    {
        //        NearPoint p_i = new NearPoint();
        //        p_i.v = boundary3[i];
        //        p_i.dist = (boundary3[p] - p_i.v).Length();
        //        p_i.index = i;
        //        nearPs.Add(p_i);
        //    }

        //    // Fit Line with k nearest points
        //    int k = 10;
        //    List<NearPoint> kNearPs = (from a in nearPs orderby a.dist ascending select a).Take(k).ToList();
        //    List<MyVector3> kNearPv = new List<MyVector3>();
        //    foreach (var kNearP in kNearPs)
        //    {
        //        kNearPv.Add(kNearP.v);
        //    }
        //    var fitline = new RansacLine3d(0.00005, 0.9);
        //    Line3 line3d = fitline.Estimate(kNearPv);

        //    // Get Local Tangential
        //    return line3d.dir;
        //}

        private struct NearPoint
        {
            public MyVector3 v;
            public double dist;
            public int index;
        }

        private void RayTracein3DPlane(List<MyVector3> points, Line3 ray, MyPlane CutPlane,MyVector3 curp, out int norInsec, out int notNorInsec)
        {
            norInsec = -1; // Normal side
            notNorInsec = -1; // Not Normal side
            double dist_left = double.MaxValue;
            double dist_right = double.MaxValue;

            for (int i = 0; i < points.Count; i++)
            {
                double dist_temp = ray.DistanceToLine(points[i]);
                if ((points[i] - curp).Dot(CutPlane.Normal()) > 0)
                {
                    // Normal side
                    if (dist_left > dist_temp)
                    {
                        dist_left = dist_temp;
                        norInsec = i;
                    }
                }
                else
                {
                    // Not Normal side
                    if (dist_right > dist_temp)
                    {
                        dist_right = dist_temp;
                        notNorInsec = i;
                    }
                }
            }
        }


        public void BuildBoundaryPoints(List<MyVector2> Boundary2)
        {
            List<int> mark = new List<int>();
            for (int i = 0; i < mark.Count; i++)
                mark.Add(0);

            for (int i = 0; i < Boundary2.Count; i++)
            {
                double dist = double.MaxValue;
                int nextIndex = -1;
                for (int j = 0; j < Boundary2.Count; j++)
                {
                    if (mark[i] == 1)
                        return;
                    double dist_temp = (Boundary2[i] - Boundary2[j]).Length();
                    if (dist > dist_temp)
                    {
                        dist = dist_temp;
                        nextIndex = j;
                    }
                }
                mark[i] = 1;
            }
        }

        List<MyVector3> ProjectTrajAccording2Profile(List<MyVector2> trajpoints)
        {
            List<MyVector3> output = new List<MyVector3>();
            MyVector3 profilenormal = this.topCircle.Normal;
            MyVector3 profilecenter = this.topCircle.Center;
            MyVector3 viewvector = this.camera.target; //suppose to be (0,0,1)

            //find a plane that is almost vertical to view vector and the normal line of profile is on the plane
            MyVector3 temp = profilenormal.Cross(viewvector);
            MyVector3 targetplanenormal = profilenormal.Cross(temp);

            targetplane = new MyPlane(profilecenter, targetplanenormal);

            output = this.Proj2dToPlane(targetplane, trajpoints);
            return output;
        }

        private List<MyVector2> ResetPath(List<MyVector2> pts, int step)
        {
            int i;
            double left;
            List<MyVector2> list = new List<MyVector2>();

            list.Add(pts[0]);
            left = step;
            for (i = 0; i < pts.Count() - 1; i++)
            {
                double tmp = (pts[i] - pts[(i + 1) % pts.Count()]).Length();
                if (tmp < left)
                {
                    left -= tmp;
                }
                else
                {
                    double offset = left;
                    while (offset < tmp)
                    {
                        MyVector2 p = pts[i] + (pts[(i + 1) % pts.Count()] - pts[i]) * (offset / tmp);
                        list.Add(p);
                        offset += step;
                    }
                    left = offset - tmp;
                }
            }
            return list;

        }

        public List<MyVector3> Proj2dToPlane(MyPlane plane, List<MyVector2> points)
        {
            List<MyVector3> output = new List<MyVector3>();
            foreach (MyVector2 p in points)
            {
                MyVector3 raynear = this.UnProject(p.x, p.y, -5);
                MyVector3 rayfar = this.UnProject(p.x, p.y, 5);
                output.Add(plane.LineIntersection(raynear, rayfar));
            }
            return output;
        }

        public class RansacLine3d
        {
            public List<MyVector3> inliers = new List<MyVector3>();
            public Line3 bestline;

            private double thres;
            private double probability;
            //* thres: if the distance between points and line are below this threshold, we regard points to be in this line;
            //* probability: the condition when we should stop current iteration. */
            public RansacLine3d(double thres_, double probability_)
            {
                this.thres = thres_;
                this.probability = probability_;
            }
            public RansacLine3d()
            {
                this.thres = 0.00005;
                this.probability = 0.9;
            }

            //public Line3 Estimate(List<MyVector3> points)
            //{
            //    int iter = 200;
            //    if (points.Count == 0) return null;
            //    Random rd = new Random();

            //    MyVector3 A = new MyVector3();
            //    MyVector3 B = new MyVector3();

            //    int maxpointinline = int.MinValue;

            //    for (int i = 0; i < iter; i++)
            //    {
            //        A = points[rd.Next(points.Count)];
            //        B = points[rd.Next(points.Count)];
            //        if (A == B) continue;  //if can't generate line

            //        Line3 testline = new Line3(A, (B - A).Normalize());
            //        List<MyVector3> tempinliers = new List<MyVector3>();
            //        int inlierscount = 0;
            //        for (int j = 0; j < points.Count; j++)
            //        {
            //            if (points[j] != A && points[j] != B)
            //            {
            //                if (testline.DistanceToLine(points[j]) < thres)
            //                {
            //                    tempinliers.Add(points[j]);
            //                    inlierscount++;
            //                }
            //            }
            //        }

            //        if (inlierscount > maxpointinline)
            //        {
            //            maxpointinline = inlierscount;
            //            this.bestline = testline.Copy();
            //            this.inliers.Clear();
            //            foreach (MyVector3 p in tempinliers)
            //                this.inliers.Add(p);
            //        }

            //        if (inlierscount >= probability * points.Count)
            //            break;
            //    }

            //    if (this.inliers.Count != 0)
            //    {
            //        double mint = double.MaxValue;
            //        double maxt = double.MinValue;
            //        foreach (MyVector3 p in this.inliers)
            //        {
            //            double t = this.bestline.ComputeT(p);
            //            if (t > maxt)
            //                maxt = t;
            //            if (t < mint)
            //                mint = t;
            //        }
            //        bestline.SetPoints(mint, maxt);
            //    }

            //    if (bestline != null && bestline.startpoint.x != double.NaN && bestline.startpoint.y != double.NaN && bestline.startpoint.z != double.NaN &&
            //    bestline.endpoint.x != double.NaN && bestline.endpoint.y != double.NaN && bestline.endpoint.z != double.NaN &&
            //        (!bestline.startpoint.IsNull() && !bestline.endpoint.IsNull()))
            //        return bestline;
            //    else
            //        return null;

            //}

        }
    }
}