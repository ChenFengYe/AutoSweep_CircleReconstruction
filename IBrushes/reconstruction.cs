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
        List<int> CorMap = new List<int>();                             // The correspondence of the BoundaryPoints and CirclePoints
        int Inter_Num = -1;                                             // The interaction time of optima
        int Inter_DMap = 1;                                             // While this time Update DistanceMap
        MyVector2 center_xy = new MyVector2();
        double center_z = 0;
        public void EstimatePlane()
        {
            mark = GetMarkImgae(this.Canvas);
            boundaryPoints_2d = GetBoundaryPoints(mark);
            DistanceMap = BuildDistanceMap(mark.Rows, mark.Cols, boundaryPoints_2d);
            drawProjectedPoint = true;

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


            MyVector3 normal = new MyVector3(FromRotationToNormal(x[0], x[1], x[2])); normal.Normalize();
            double radius = x[3];

            topCircle = new MyCircle(center, radius, new MyPlane(center, normal));
            CirclePoints_2d = GetProjectionPoints_2D(topCircle.CirclePoints);

            this.view.Refresh();
        }

        private MyVector3 FromRotationToNormal(double thetaX, double thetaY, double thetaZ)
        {
            MyMatrix4d RotationX = MyMatrix4d.RotationMatrix(new MyVector3(1, 0, 0), thetaX * Math.PI);
            MyMatrix4d RotationY = MyMatrix4d.RotationMatrix(new MyVector3(0, 1, 0), thetaY * Math.PI);
            MyMatrix4d RotationZ = MyMatrix4d.RotationMatrix(new MyVector3(0, 0, 1), thetaZ * Math.PI);
            MyVector3 normal_init = new MyVector3(0,-1,0);
            MyVector3 normal_cur = (RotationZ * RotationY * RotationX * normal_init.ToMyVector4()).XYZ();
            return normal_cur.Normalize();
        }

        private void function_project(double[] x, double[] fi, object obj)
        {
            // Step 1: Init Params
            Inter_Num++;
            MyVector3 center = new MyVector3(center_xy.x, center_xy.y, center_z);
            //MyVector3 center = new MyVector3(x[4], x[5], x[6]);

            MyVector3 normal = new MyVector3(FromRotationToNormal(x[0], x[1], x[2])); normal.Normalize();

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
                // Cor Map
                //int corIndex = CorMap[i];
                //double dist = Math.Pow((CirclePoints_2d[i] - boundaryPoints_2d[corIndex]).Length(), 2);

                //------------------------------------------------------------------------------------
                // Distance Map
                //int i_x = (int)Math.Round(CirclePoints_2d[i].x);
                //int i_y = (int)Math.Round(CirclePoints_2d[i].y);
                //i_x = Math.Min(i_x, mark.Cols - 1);
                //i_x = Math.Max(0, i_x);
                //i_y = Math.Min(i_y, mark.Rows - 1);
                //i_y = Math.Max(0, i_y);
                //double dist = DistanceMap[i_y][i_x];

                //------------------------------------------------------------------------------------
                // Distance Interpolation from distance map
                double dist = this.InterpolateDistanceWithDisMap(CirclePoints_2d[i]);



                //------------------------------------------------------------------------------------
                //double dist = double.MaxValue;
                //for (int j = 0; j < boundaryPoints_2d.Count; j++)
                //{
                //    dist = Math.Min(dist, Math.Pow((CirclePoints_2d[i] - boundaryPoints_2d[j]).Length(), 2));
                //}


                dist_all += dist;
                //dist_max = Math.Max(dist, dist_max);
            }

            // Set Cost function
            fi[0] = dist_all;

            //////center 3d match 2d
            //double centerdis = (this.Compute2D(center) - center_xy).SquareLength();
            //fi[1] = centerdis;
            fi[1] = normal.SquareLength() - 1;

            System.Console.WriteLine("{0}|| N: {1},{2},{3} r: {4} cost:{5}",
                Inter_Num, x[0], x[1], x[2], x[3], dist_all);
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
            List<List<double>> Dmap = new List<List<double>>();
            for (int i = 0; i < rows; i++)
            {
                List<double> rowDmap = new List<double>();
                for (int j = 0; j < cols; j++)
                {
                    double dist = double.MaxValue;
                    for (int k = 0; k < BoundaryPoints.Count; k++)
                    {
                        dist = Math.Min(dist, (new MyVector2(j, i) - BoundaryPoints[k]).Length());
                    }
                    rowDmap.Add(dist);
                }
                Dmap.Add(rowDmap);
            }
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

        public List<MyVector2> GetBoundaryPoints(Image<Gray, byte> markImg)
        {
            double cannyThreshold = 60;
            double cannyThresholdLinking = 100;
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
            lowery = Math.Min(lowery, mark.Rows - 1);
            lowery = Math.Max(0, lowery);

            upperx = Math.Min(upperx, mark.Cols - 1);
            upperx = Math.Max(0, upperx);
            uppery = Math.Min(uppery, mark.Rows - 1);
            uppery = Math.Max(0, uppery);


            double distance11 = DistanceMap[lowery][lowerx];
            double distance12 = DistanceMap[lowery][upperx];
            double distance21 = DistanceMap[uppery][lowerx];
            double distance22 = DistanceMap[uppery][upperx];

            if (lowerx == upperx && lowery == uppery) return distance11;

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
    }
}