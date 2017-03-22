using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Drawing;

using Loyc.Geometry;

using OpenTK.Graphics.OpenGL;

using Emgu.CV;
using Emgu.CV.CvEnum;

using Accord.Math;
using Accord.MachineLearning.Geometry;
using Accord.MachineLearning;
using Accord.Math.Geometry;

using PolygonDecompose;
using SmartCanvas;

namespace MyGeometry
{
    public class MyPlane
    {
        //liyuwei - 1217 - draw ori plane center while moving
        public MyVector3 oricenterformoveingcanvas = new MyVector3();
        public Plane planeEquation;
        public MyVector3 planeCenter;
        public Quad3D boundQuad = null;
        public List<MyVector3> planePoints = new List<MyVector3>();
        public List<MyVector3> planeVertices = new List<MyVector3>();

        public List<MyVector2> planeVertices2d = new List<MyVector2>();

        public Line3 normalline;
        //public MyVector3 planenormal;

        public List<Line3> axes = new List<Line3>();

        public Color planeColor;
        //Random rnd = new Random(127);
        //for candidate generation
        public int from = -1;
        public int to = -1;
        public int meanshiftId = -1;

        //tag = 0 initial plane
        //tag = 1 sweep plane
        //tag = -1 sweep curved surface
        public int tag = 0;

        public bool choosen = false;
        public double area;
        public MyPlane()
        {
            Random rnd = new Random();
            planeColor = Color.FromArgb(150, rnd.Next(100, 255), rnd.Next(100, 255), rnd.Next(100, 255));
        }
        public MyPlane(MyVector3 v1, MyVector3 v2, MyVector3 v3)
        {
            planeVertices.Add(v1);
            planeVertices.Add(v2);
            planeVertices.Add(v3);

            planePoints.Add(v1);
            planePoints.Add(v2);
            planePoints.Add(v3);

            MyVector3 planenormal = (v1 - v2).Cross(v1 - v3);
            planenormal.Normalize();
            double offset = -v1.Dot(planenormal);
            planeEquation = new Plane((float)planenormal.x, (float)planenormal.y, (float)planenormal.z, (float)offset);

            planeCenter = (v1 + v2 + v3) / 3;
            planeColor = Color.FromArgb(150, 255, 255, 0);

            // Accord.Math.Point3 _v1 = new Point3((float)v1.x, (float)v1.y, (float)v1.z);
            // Accord.Math.Point3 _v2 = new Point3((float)v2.x, (float)v2.y, (float)v2.z);
            //  Accord.Math.Point3 _v3 = new Point3((float)v3.x, (float)v3.y, (float)v3.z);
            // planeEquation = Plane.FromPoints(_v1, _v2, _v3);

            //ComputePlaneArea();
            ComputeBoundQuad();
        }
        public MyPlane(Plane _p, MyVector3 _center, List<Accord.Math.Point3> _points, List<MyVector3> _vertices)
        {
            Random rnd = new Random();
            planeColor = Color.FromArgb(150, rnd.Next(100, 255), rnd.Next(100, 255), rnd.Next(100, 255));
            planeEquation = _p;
            planeCenter = _center;
            foreach (Accord.Math.Point3 point in _points)
            {
                planePoints.Add(new MyVector3(point.X, point.Y, point.Z));
            }
            planeVertices = _vertices;
            ComputeBoundQuad();
        }
        public MyPlane(Plane _p, MyVector3 _center, List<Accord.Math.Point3> _points, int _from, int _to)
        {
            Random rnd = new Random();
            planeColor = Color.FromArgb(150, rnd.Next(100, 255), rnd.Next(100, 255), rnd.Next(100, 255));
            planeEquation = _p;
            planeCenter = _center;
            foreach (Accord.Math.Point3 point in _points)
            {
                planePoints.Add(new MyVector3(point.X, point.Y, point.Z));
            }
            from = _from;
            to = _to;
            ComputeBoundQuad();
        }
        public MyPlane(Plane _p, MyVector3 _center, int _from, int _to)
        {
            Random rnd = new Random();
            planeColor = Color.FromArgb(150, rnd.Next(100, 255), rnd.Next(100, 255), rnd.Next(100, 255));
            planeEquation = _p;
            planeCenter = _center;
            from = _from;
            to = _to;
        }
        public MyPlane(MyPlane p)
        {
            planeCenter = p.planeCenter;
            planeColor = p.planeColor;
            planeEquation = new Plane(p.planeEquation.A, p.planeEquation.B, p.planeEquation.C, p.planeEquation.Offset);
            planePoints = new List<MyVector3>(p.planePoints);
            planeVertices = new List<MyVector3>(p.planeVertices);
            choosen = p.choosen;
            ComputeBoundQuad();
            this.ProjectVerticesTo2d();
        }
        public MyPlane(Plane _p, List<MyVector3> _points)
        {
            Random rnd = new Random();
            planeColor = Color.FromArgb(150, rnd.Next(100, 255), rnd.Next(100, 255), rnd.Next(100, 255));
            planePoints = _points;
            planeEquation = _p;
            ComputeCenter();
            ComputeBoundQuad();
            this.ProjectVerticesTo2d();
        }
        public MyPlane(MyVector3 point, MyVector3 normal)
        {
            double d = -point.Dot(normal);
            planeEquation = new Plane((float)normal.x, (float)normal.y, (float)normal.z, (float)d);
            planeCenter = point;
            planeColor = Color.FromArgb(150, 255, 255, 0);

            // add three plane points
            MyVector3 x = new MyVector3(normal.y, -normal.x, 0);
            if (x.x == 0 && x.y == 0)
                x = new MyVector3(1, 0, 0);
            MyVector3 y = normal.Cross(x).Normalize();
            CoordinateFrame frame = new CoordinateFrame(point, x, y, normal);
            double s = 0.8;
            MyVector3 v1 = frame.GetPointLocalCoord(new MyVector3(s, s, 0));
            MyVector3 v2 = frame.GetPointLocalCoord(new MyVector3(-s, s, 0));
            MyVector3 v3 = frame.GetPointLocalCoord(new MyVector3(-s, -s, 0));
            MyVector3 v4 = frame.GetPointLocalCoord(new MyVector3(s, -s, 0));

            planeVertices.Add(v1);
            planeVertices.Add(v2);
            planeVertices.Add(v3);
            planeVertices.Add(v4);

            planePoints.Add(v1);
            planePoints.Add(v2);
            planePoints.Add(v3);
            planePoints.Add(v4);

            this.ComputeVertices();
            this.ComputeBoundQuad();
            this.ComputeCenter();
        }

        public MyPlane(List<MyVector3> points, bool SP = true)//sweep points
        {
            planePoints = new List<MyVector3>(points);

            double[,] A = new double[3, 3];
            double[,] B = new double[3, 1];
            foreach (MyVector3 point in points)
            {
                A[0, 0] += point.x * point.x;
                A[0, 1] += point.x * point.y;
                A[0, 2] += point.x;

                A[1, 0] += point.x * point.y;
                A[1, 1] += point.y * point.y;
                A[1, 2] += point.y;

                A[2, 0] += point.x;
                A[2, 1] += point.y;
                A[2, 2] = points.Count();

                B[0, 0] += point.x * point.z;
                B[1, 0] += point.y * point.z;
                B[2, 0] += point.z;
            }

            double[,] x = Matrix.Solve(A, B, leastSquares: true);
            MyVector3 ori_normal = new MyVector3(x[0, 0] / x[2, 0], x[1, 0] / x[2, 0], -1 / x[2, 0]);
            MyVector3 normal = ori_normal.Normalize();
            double offset = 1.0 / ori_normal.Length();
            planeEquation = new Plane((float)normal.x, (float)normal.y, (float)normal.z, (float)offset);
            Random rnd = new Random();
            planeColor = Color.FromArgb(150, rnd.Next(100, 255), rnd.Next(100, 255), rnd.Next(100, 255));
            //this.planenormal = new MyVector3(this.planeEquation.A, this.planeEquation.B, this.planeEquation.C);


            if (SP)
            {
                for (int i = 0; i < planePoints.Count; i++)
                {
                    if (this.DistanceToPoint(planePoints[i]) > 0.001)
                    {
                        tag = -1;//curved plane
                        return;
                    }
                }
                planeColor = Color.FromArgb(150, 255, 255, 0);
                tag = 1;//sweep plane
            }

            for (int i = 0; i < planePoints.Count; i++)
            {
                planePoints[i] = this.ProjectPoint(planePoints[i]);
            }
            ComputeCenter();

            choosen = true;
            if (planePoints.Count > 3)
            {
                ComputeVertices();
            }
            else
                planeVertices.AddRange(planePoints);
            //this.Scale(1.2);
            this.ProjectVerticesTo2d();
            ComputeBoundQuad();
        }
        public MyPlane(string filename)
        {
            using (StreamReader sr = new StreamReader(File.Open(filename, FileMode.Open)))
            {
                planeEquation = new Plane(float.Parse(sr.ReadLine()), float.Parse(sr.ReadLine()), float.Parse(sr.ReadLine()), float.Parse(sr.ReadLine()));
                planeCenter = new MyVector3(float.Parse(sr.ReadLine()), float.Parse(sr.ReadLine()), float.Parse(sr.ReadLine()));
                planeColor = Color.FromArgb(int.Parse(sr.ReadLine()), int.Parse(sr.ReadLine()), int.Parse(sr.ReadLine()), int.Parse(sr.ReadLine()));

                double x, y, z;
                int verticescount = int.Parse(sr.ReadLine());
                for (int i = 0; i < verticescount; i++)
                {
                    x = double.Parse(sr.ReadLine());
                    y = double.Parse(sr.ReadLine());
                    z = double.Parse(sr.ReadLine());
                    planeVertices.Add(new MyVector3(x, y, z));
                }
                this.Scale(1.1);
                ProjectVerticesTo2d();
                // ComputePlaneArea();
                //int pointcount = int.Parse(sr.ReadLine());
                //for (int i = 0; i < pointcount; i++)
                //{
                //    x = double.Parse(sr.ReadLine());
                //    y = double.Parse(sr.ReadLine());
                //    z = double.Parse(sr.ReadLine());
                //    planePoints.Add(new MyVector3(x, y, z));
                //}
                sr.Close();
            }
            //this.planenormal = (pointToPlaneClustering.RANSACNormal(this.planeVertices)).Normalize();
            //double offset = -this.planenormal.Dot(this.planeCenter);
            //this.planeEquation = new Plane((float)this.planenormal.x, (float)this.planenormal.y, (float)this.planenormal.z, (float)offset);
            this.choosen = true;
            ComputeBoundQuad();
            // foreach (MyVector3 p in planePoints)
            //   System.Console.WriteLine(planeEquation.DistanceToPoint(new Accord.Math.Point3((float)p.x, (float)p.y, (float)p.z)));
        }


        public void Scale(double ratio)
        {
            for (int i = 0; i < this.planeVertices.Count; i++)
            {
                Line3 line = new Line3(this.planeCenter, (this.planeVertices[i] - this.planeCenter).Normalize());
                double t = line.ComputeT(this.planeVertices[i]);

                //t *= ratio;
                t += (ratio / 1000);
                this.planeVertices[i] = line.GetPointwithT(t);
            }
        }


        public int SaveMyPlane(string filename)
        {
            if (File.Exists(filename + ".plane"))
                return 0;
            string dir = Environment.CurrentDirectory;
            //实例化一个文件流--->与写入文件相关联
            FileStream fs = new FileStream(filename + ".plane", FileMode.Create);
            //实例化一个StreamWriter-->与fs相关联
            StreamWriter sw = new StreamWriter(fs);
            //开始写入
            sw.WriteLine(planeEquation.A);
            sw.WriteLine(planeEquation.B);
            sw.WriteLine(planeEquation.C);
            sw.WriteLine(planeEquation.Offset);

            sw.WriteLine(planeCenter.x);
            sw.WriteLine(planeCenter.y);
            sw.WriteLine(planeCenter.z);

            sw.WriteLine(planeColor.A);
            sw.WriteLine(planeColor.R);
            sw.WriteLine(planeColor.G);
            sw.WriteLine(planeColor.B);

            sw.WriteLine(planeVertices.Count);
            foreach (MyVector3 vertex in planeVertices)
            {
                sw.WriteLine(vertex.x);
                sw.WriteLine(vertex.y);
                sw.WriteLine(vertex.z);
            }

            sw.WriteLine(planePoints.Count);
            foreach (MyVector3 point in planePoints)
            {
                sw.WriteLine(point.x);
                sw.WriteLine(point.y);
                sw.WriteLine(point.z);
            }

            //清空缓冲区
            sw.Flush();
            //关闭流
            sw.Close();
            fs.Close();
            return 1;
        }

        public MyVector3 Normal()
        {
            //return planenormal;

            return new MyVector3(this.planeEquation.A, this.planeEquation.B, this.planeEquation.C);
        }
        public void ComputeBoundQuad()
        {
            MyVector3 normal = Normal();
            MyVector3 xaxis = new MyVector3(1, 0, 0);
            MyVector3 yaxis = new MyVector3(0, 1, 0);
            MyVector3 zaxis = new MyVector3(0, 0, 1);
            MyVector3 rotAxis = zaxis.Cross(normal).Normalize();

            double cos = zaxis.Dot(normal);
            if (cos > 1) cos = 1; if (cos < -1) cos = -1;
            double rotAngle = Math.Acos(cos);
            MyMatrix4d R = MyMatrix4d.RotationMatrix(rotAxis, rotAngle);
            MyVector3 new_xaxis = (R * new MyVector4(xaxis)).XYZ().Normalize();
            MyVector3 new_yaxis = (R * new MyVector4(yaxis)).XYZ().Normalize();
            CoordinateFrame frame = new CoordinateFrame(this.planeCenter, new_xaxis, new_yaxis, normal);

            double xmin = double.MaxValue; double ymin = double.MaxValue;
            double xmax = double.MinValue; double ymax = double.MinValue;
            foreach (MyVector3 v in this.planeVertices)
            {
                MyVector3 u = frame.GetPointLocalCoord(v);
                xmin = Math.Min(u.x, xmin);
                xmax = Math.Max(u.x, xmax);
                ymin = Math.Min(u.y, ymin);
                ymax = Math.Max(u.y, ymax);
            }

            MyVector3 v1 = new MyVector3(xmin, ymin, 0);
            MyVector3 v2 = new MyVector3(xmax, ymin, 0);
            MyVector3 v3 = new MyVector3(xmax, ymax, 0);
            MyVector3 v4 = new MyVector3(xmin, ymax, 0);

            SmartCanvas.Point3[] pts3 =
                new SmartCanvas.Point3[4] {
				new SmartCanvas.Point3(frame.GetPointSpaceCoord(v1)),
				new SmartCanvas.Point3(frame.GetPointSpaceCoord(v2)),
				new SmartCanvas.Point3(frame.GetPointSpaceCoord(v3)),
				new SmartCanvas.Point3(frame.GetPointSpaceCoord(v4))
			};
            this.boundQuad = new Quad3D(pts3);

        }
        public void ComputeVertices()
        {
            ////use convex hull to find vertices
            planeVertices.Clear();
            List<Point<double>> points2 = new List<Point<double>>();
            foreach (MyVector3 p in planePoints)
            {
                Point<double> pt = new Point<double>(p.x, p.y);
                points2.Add(pt);
            }
            List<Point<double>> convex_hull = PointMath.ComputeConvexHull(points2, false).ToList();

            foreach (Point<double> point in convex_hull)
            {
                MyVector3 u = new MyVector3(point.X, point.Y, ComputeZ(point.X, point.Y));
                planeVertices.Add(u);
            }
            ProjectVerticesTo2d();
            //ComputePlaneArea();
        }
        public double ComputeZ(double x, double y)
        {
            return (x * planeEquation.A + y * planeEquation.B + planeEquation.Offset) / (-planeEquation.C);
        }
        public void ComputeCenter()
        {
            planeCenter = new MyVector3(0, 0, 0);
            if (planePoints.Count() > 0)
            {
                foreach (MyVector3 point in planePoints)
                {
                    planeCenter += point;
                }
                planeCenter = planeCenter / planePoints.Count;
            }
            else
            {
                foreach (MyVector3 point in planeVertices)
                {
                    planeCenter += point;
                }
                planeCenter = planeCenter / planeVertices.Count;
            }
        }
        public MyVector3 ProjectPoint(MyVector3 p)
        {
            double a = planeEquation.A;
            double b = planeEquation.B;
            double c = planeEquation.C;
            double t = -(p.x * a + b * p.y + c * p.z + planeEquation.Offset) / (planeEquation.Normal.Square);
            //MyVector3 r = new MyVector3(a * t + p.x, b * t + p.y, c * t + p.z);
            return new MyVector3(a * t + p.x, b * t + p.y, c * t + p.z);
        }
        public MyVector3 ProjectPointWithNewcenter(MyVector3 p, MyVector3 newcenter)
        {
            double a = planeEquation.A;
            double b = planeEquation.B;
            double c = planeEquation.C;
            double noffset = -(a * newcenter.x + b * newcenter.y + c * newcenter.z);

            double t = -(p.x * a + b * p.y + c * p.z + noffset) / (planeEquation.Normal.Square);
            MyVector3 r = new MyVector3(a * t + p.x, b * t + p.y, c * t + p.z);
            return new MyVector3(a * t + p.x, b * t + p.y, c * t + p.z);
        }
        public bool InBoundQuad(MyVector3 p)
        {
            return this.boundQuad.InQuad(p);

        }
        public bool InPlanePolygon(MyVector3 p)//the sum of inner angle is 360
        {
            //this.ProjectVerticesTo2d();
            if (this.planeVertices.Count == 0) return false;

            double theta = 0;
            MyVector3 v1 = p - planeVertices[0];
            MyVector3 v2 = p - planeVertices[planeVertices.Count - 1];
            theta = Math.Acos(v1.Dot(v2) / (v1.Length() * v2.Length()));
            for (int i = 0; i < planeVertices.Count - 1; i++)
            {
                v1 = p - planeVertices[i];
                v2 = p - planeVertices[i + 1];
                theta += Math.Acos(v1.Dot(v2) / (v1.Length() * v2.Length()));
            }

            if (Math.Abs(2 * Math.PI - theta) < 1e-10)
                return true;
            else
                return false;
        }
        public bool InPlanePolygon2d(MyVector2 p)
        {
            this.ProjectVerticesTo2d();

            if (this.planeVertices.Count == 0) return false;
            double theta = 0;
            MyVector2 v1 = p - planeVertices2d[0];
            MyVector2 v2 = p - planeVertices2d[planeVertices2d.Count - 1];
            theta = Math.Acos(v1.Dot(v2) / (v1.Length() * v2.Length()));
            for (int i = 0; i < planeVertices.Count - 1; i++)
            {
                v1 = p - planeVertices2d[i];
                v2 = p - planeVertices2d[i + 1];
                theta += Math.Acos(v1.Dot(v2) / (v1.Length() * v2.Length()));
            }

            if (Math.Abs(2 * Math.PI - theta) < 1e-10)
                return true;
            else
                return false;
        }
        public bool PointInPlane2d(MyVector2 p)
        {
            this.ProjectVerticesTo2d();
            bool c = false;
            List<MyVector2> points = this.planeVertices2d;
            int n = points.Count;
            for (int i = 0, j = n - 1; i < n; j = i++)
            {
                if (((points[i].y > p.y) != (points[j].y > p.y)) &&
                    (p.x < (points[j].x - points[i].x) * (p.y - points[i].y) / (points[j].y - points[i].y) + points[i].x))
                    c = !c;
            }
            return c;
        }
        public double DistanceToPoint(MyVector3 a)
        {
            Accord.Math.Point3 point = new Accord.Math.Point3((float)a.x, (float)a.y, (float)a.z);
            return this.planeEquation.DistanceToPoint(point);
        }
        public void ProjectVerticesTo2d()
        {
            double[] KINECT_LOCAL_PROJ = new double[]{ // a large dataset of object scans
                    525, 0.0, 319.5, 0.0,
                    0.0, 525, 239.5, 0.0,
                    0.0, 0.0, 1.0, 0.0 };
            for (int i = 0; i < planeVertices.Count; i++)
            {

                double fx = KINECT_LOCAL_PROJ[0];
                double fy = KINECT_LOCAL_PROJ[5];
                double u0 = KINECT_LOCAL_PROJ[2];
                double v0 = KINECT_LOCAL_PROJ[6];

                planeVertices2d.Add(new MyVector2(fx * planeVertices[i].x / planeVertices[i].z + u0, fy * planeVertices[i].y / planeVertices[i].z + v0));
            }
        }


        public MyVector3 LineIntersection(MyVector3 v1, MyVector3 v2)
        {
            if (Math.Abs((v1 - v2).Normalize().Dot(this.Normal())) < 1e-7)
                return new MyVector3();

            double a = (v1 - v2).Normalize().x;
            double b = (v1 - v2).Normalize().y;
            double c = (v1 - v2).Normalize().z;

            double A = this.planeEquation.A;
            double B = this.planeEquation.B;
            double C = this.planeEquation.C;
            double D = this.planeEquation.Offset;

            double t = (-D - A * v1.x - B * v1.y - C * v1.z) / (A * a + B * b + C * c);

            MyVector3 result = new MyVector3();
            result.x = a * t + v1.x;
            result.y = b * t + v1.y;
            result.z = c * t + v1.z;

            return result;
        }

        public MyVector3 RayIntersection(MyVector3 v1, MyVector3 dir)
        {
            if (Math.Abs(dir.Dot(this.Normal())) < 1e-7)
                return new MyVector3();

            double a = dir.x;
            double b = dir.y;
            double c = dir.z;

            double A = this.planeEquation.A;
            double B = this.planeEquation.B;
            double C = this.planeEquation.C;
            double D = this.planeEquation.Offset;

            double t = (-D - A * v1.x - B * v1.y - C * v1.z) / (A * a + B * b + C * c);

            MyVector3 result = new MyVector3();
            result.x = a * t + v1.x;
            result.y = b * t + v1.y;
            result.z = c * t + v1.z;

            return result;
        }

        public static void Translate(MyPlane last, out MyPlane next, float Dist)
        {
            Accord.Math.Point3 new_normal = last.planeEquation.Normal;
            float new_offset = last.planeEquation.Offset - Dist;
            Plane new_p = new Plane(new_normal, new_offset);

            MyVector3 new_center = new MyVector3();
            List<Accord.Math.Point3> new_points = new List<Accord.Math.Point3>();
            List<MyVector3> new_vertices = new List<MyVector3>();

            next = new MyPlane(new_p, new_center, new_points, new_vertices);
        }
        public float Offset()
        {
            return planeEquation.Offset;
        }

        public void SaveMyPlane2(string filename)
        {
            using (StreamWriter sw = new StreamWriter(filename + ".plane", false)) //if filename already exist, rewrite it
            {
                sw.WriteLine("{0} {1} {2} {3}", planeEquation.A, planeEquation.B, planeEquation.C, planeEquation.Offset);
                sw.WriteLine(planeCenter);
                sw.WriteLine("{0} {1} {2} {3}", planeColor.A, planeColor.R, planeColor.G, planeColor.B);

                sw.WriteLine(planeVertices.Count);
                foreach (MyVector3 vertex in planeVertices)
                {
                    sw.WriteLine(vertex);
                }

                sw.WriteLine(planePoints.Count);
                foreach (MyVector3 point in planePoints)
                {
                    sw.WriteLine(point);
                }
                sw.Flush();
                sw.Close();
            }
        }



        public void DrawMyPlane(byte opacity = 50)
        {
            // draw bounding quad
            //if (this.boundQuad != null)
            //{
            //    GL.Disable(EnableCap.Lighting);
            //    this.boundQuad.DrawBlended(this.planeColor, opacity);
            //   this.boundQuad.DrawOutLine(Color.Black);
            //    GL.Enable(EnableCap.Lighting);
            //}
            //return;


            //////draw convex hull
            GL.Disable(EnableCap.Lighting);
            //polygon
            GL.Begin(PrimitiveType.Polygon);
            GL.Color4(this.planeColor.R, this.planeColor.G, this.planeColor.B, opacity);
            for (int i = 0; i < this.planeVertices.Count(); i++)
            {
                GL.Vertex3(this.planeVertices[i].x, this.planeVertices[i].y, this.planeVertices[i].z);
            }
            GL.End();
            GL.Enable(EnableCap.Lighting);
        }



        public static bool IsVertical(MyPlane p1, MyPlane p2, double angle = 85)
        {
            MyVector3 _v1 = p1.Normal();
            MyVector3 _v2 = p2.Normal();
            if (Math.Abs(_v1.Dot(_v2) / (_v1.Length() * _v2.Length())) <= Math.Cos(angle * Math.PI / 180))
                return true;
            else
                return false;
        }
        public static bool IsParallel(MyPlane p1, MyPlane p2, double angle = 5)
        {
            MyVector3 _v1 = p1.Normal();
            MyVector3 _v2 = p2.Normal();
            if (Math.Abs(_v1.Dot(_v2) / (_v1.Length() * _v2.Length())) >= Math.Cos(angle * Math.PI / 180))
                return true;
            else
                return false;
        }

        public static bool IsParallel(MyVector3 _v1, MyVector3 _v2, double angle = 1)
        {
            if ((_v1.Dot(_v2) / (_v1.Length() * _v2.Length())) >= Math.Cos(angle * Math.PI / 180))
                return true;
            else
                return false;
        }
        public static MyVector3 RotateAroundAxis(double angle, MyVector3 v, MyVector3 a)
        {
            double s = Math.Sin(angle);
            double c = Math.Cos(angle);
            MyMatrix3d m = new MyMatrix3d();

            m[0, 0] = a.x * a.x * (1 - c) + c;
            m[0, 1] = a.x * a.y * (1 - c) + a.z * s;
            m[0, 2] = a.x * a.z * (1 - c) - a.y * s;

            m[1, 0] = a.x * a.y * (1 - c) - a.z * s;
            m[1, 1] = a.y * a.y * (1 - c) + c;
            m[1, 2] = a.y * a.z * (1 - c) + a.x * s;

            m[2, 0] = a.x * a.z * (1 - c) + a.y * s;
            m[2, 1] = a.y * a.z * (1 - c) - a.x * s;
            m[2, 2] = a.z * a.z * (1 - c) + c;

            return m * v;


        }



    }
}