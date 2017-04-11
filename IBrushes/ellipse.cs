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
        //identify if a ellipse

        public double ComputeBlobScore()
        {
            return ErrorScore(GetMaskPoints());
        }


        public List<MyVector2> GetMaskPoints()
        {
            //new ImageViewer(this.mark, "ori mask img").Show();

            Image<Gray, byte> blobimg = new Image<Gray, byte>(this.mark.Width, this.mark.Height);
            blobimg.SetZero();
            double backGround_threshold = 50;
            List<MyVector2> b_points = new List<MyVector2>();
            for (int i = 0; i < mark.Height; i++)
            {
                for (int j = 0; j < mark.Width; j++)
                {
                    if (mark[i, j].Intensity > backGround_threshold)
                    {
                        b_points.Add(new MyVector2(j, i));
                        blobimg[i, j] = new Gray(255);
                    }
                }
            }
            //new ImageViewer(blobimg, "blob extraction").Show();

            blobpoints = b_points;
            return b_points;

        }

        List<MyVector2> blobpoints = new List<MyVector2>();
        MyVector2 mean = new MyVector2();
        MyVector2 majorendp = new MyVector2();
        MyVector2 majorstartp = new MyVector2();
        MyVector2 minorstartp = new MyVector2();
        MyVector2 minorendp = new MyVector2();
        public void BlobDraw()
        {
            GL.PointSize(1.0f);
            GL.Begin(PrimitiveType.Points);
            GL.Color3(Color.Yellow);
            foreach (MyVector2 p in this.blobpoints)
            {
                GL.Vertex2(p.x, p.y);
            }
            GL.End();

            GL.Begin(PrimitiveType.Lines);
            GL.Color3(Color.Red);
            GL.Vertex2(majorstartp.x, majorstartp.y);
            GL.Vertex2(majorendp.x, majorendp.y);
            GL.Color3(Color.Blue);
            GL.Vertex2(minorstartp.x, minorstartp.y);
            GL.Vertex2(minorendp.x, minorendp.y);
            GL.End();
        }



        public double ErrorScore(List<MyVector2> points)
        {
            foreach (MyVector2 p in points)
                mean += p;

            mean /= points.Count;
            double W = points.Count;


            double[,] C = new double[2, 2];
            foreach (MyVector2 point in points)
            {
                C[0, 0] += (point.x - mean.x) * (point.x - mean.x);
                C[0, 1] += (point.x - mean.x) * (point.y - mean.y);

                C[1, 0] += (point.y - mean.y) * (point.x - mean.x);
                C[1, 1] += (point.y - mean.y) * (point.y - mean.y);
            }

            C[0, 0] /= W;
            C[0, 1] /= W;
            C[1, 0] /= W;
            C[1, 1] /= W;

            Matrix2d CM = new Matrix2d(C);
            SVD svd = new SVD(C);
            //svd.w - eigen value, start from 1
            //svd.u - eigen vector, start from 1

            int max = 1, min = 2;
            if (svd.w[max] < svd.w[min])
            {
                int temp = max;
                max = min;
                min = temp;
            }

            double major = 2 * Math.Sqrt(svd.w[max]);
            MyVector2 majoraxis = new MyVector2(svd.u[1, max], svd.u[2, max]);
            majoraxis = majoraxis.Normalize();
            double minor = 2 * Math.Sqrt(svd.w[min]);
            MyVector2 minoraxis = new MyVector2(svd.u[1, min], svd.u[2, min]);
            minoraxis = minoraxis.Normalize();

            majorendp = mean + majoraxis * major;
            majorstartp = mean - majoraxis * major;
            minorendp = mean + minoraxis * minor;
            minorstartp = mean - minoraxis * minor;

            double error = Math.Abs(W - 4 * Math.PI * Math.Sqrt(CM.Det()));
            error /= W;
            Console.WriteLine("Like a ellipse error: {0}", error); //10^-2 may be a good threshold

            return error;
        }
    }
}