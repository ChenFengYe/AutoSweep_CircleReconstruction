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
        private List<MyVector2> CirclePoints_2d = new List<MyVector2>();// The Circle Sample Points projected on the screem
        List<List<double>> DistanceMap = new List<List<double>>();      // Distance Map
        List<int> CorMap = new List<int>();                             // The correspndence of the BoundaryPoints and CirclePoints
        int Inter_Num = -1;                                             // The interaction time of optima
        int Inter_DMap = 10;                                            // While this time Update DistanceMap
        double center_z = 0;
        MyVector2 center_xy = new MyVector2();


        public void EstimatePlane()
        {
            mark = GetMarkImgae(this.rgbdBuilder.CurrFrame().Image);
            boundaryPoints_2d = GetBoundaryPoints(mark);
            DistanceMap = BuildDistanceMap(mark.Rows, mark.Cols, boundaryPoints_2d);
            drawProjectedPoint = true;

            //--------------------------- Do the Optimization----------------------------------//
            // Init Optima Params                             //                //
            center_z = 0.05;
            center_xy = GetMarkCenter(boundaryPoints_2d, center_z);
            double[] bndl = new double[] { -1, -1, -1, 0.0001 };
            double[] x = new double[] { 0, -0.5, -0.5, 0.02 };
            double[] bndu = new double[] { 1, 0.1, 0.1, 0.1 };
            int paramsNum = 4;
            double diffstep = 0.000001;
            double epsg = 0.000000000001;
            double epsf = 0;
            double epsx = 0;
            int maxits = 0;
            alglib.minlmstate state;
            alglib.minlmreport rep;

            // Do the optima
            alglib.minlmcreatev(paramsNum, x, diffstep, out state);
            alglib.minlmsetbc(state, bndl, bndu);
            alglib.minlmsetcond(state, epsg, epsf, epsx, maxits);
            alglib.minlmoptimize(state, function_project, null, null);
            alglib.minlmresults(state, out x, out rep);

            // Update Circle
            MyVector3 center = new MyVector3(center_xy.x, center_xy.y, center_z);
            MyVector3 normal = new MyVector3(x[0], x[1], x[2]); normal.Normalize();
            double radius = x[3];

            topCircle = new MyCircle(center, radius, new MyPlane(center, normal));
            CirclePoints_2d = GetProjectionPoints_2D(topCircle.CirclePoints);

            // Cout Param
            //System.Console.WriteLine("{0}", alglib.ap.format(x, 2));    // EXPECTED: [-3,+3]
        }

        private void function_project(double[] x, double[] fi, object obj)
        {
            // Step 1: Init Params
            Inter_Num++;
            MyVector3 center = new MyVector3(center_xy.x, center_xy.y, center_z);
            MyVector3 normal = new MyVector3(x[0], x[1], x[2]); normal.Normalize();
            double radius = x[3];
            MyCircle myC = new MyCircle(center, radius, new MyPlane(center, normal));

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
                int i_x = (int)Math.Round(CirclePoints_2d[i].x);
                int i_y = (int)Math.Round(CirclePoints_2d[i].y);
                i_x = Math.Min(i_x, mark.Cols - 1);
                i_x = Math.Max(0, i_x);
                i_y = Math.Min(i_y, mark.Rows - 1);
                i_y = Math.Max(0, i_y);
                double dist = DistanceMap[i_y][i_x];
                //------------------------------------------------------------------------------------
                // Caculate Every time
                //double dist = double.MaxValue;
                //for (int j = 0; j < boundaryPoints_2d.Count; j++)
                //{
                //    dist = Math.Min(dist, Math.Pow((CirclePoints_2d[i] - boundaryPoints_2d[j]).Length(), 2));
                //}
                dist_all += dist;
                dist_max = Math.Max(dist, dist_max);
            }

            // Set Cost function
            fi[0] = dist_all;
            //fi[1] = dist_max;

            System.Console.WriteLine("{0}|| N: {0},{1},{2} r: {4} cost:{5}",
                Inter_Num, x[0], x[1], x[2], x[3], dist_all);
        }

        double offsety = 0;
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
            //new ImageViewer(cannyimg, "line and corner").Show();
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
    }
}