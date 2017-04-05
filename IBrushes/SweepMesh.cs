using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Runtime.InteropServices;
using System.IO;
using System.Drawing;
using MyGeometry;
using SmartCanvas;

namespace SmartCanvas
{
    public class SweepMesh
    {
        List<MyVector3> profile = new List<MyVector3>();
        List<MyVector3> trajectory1 = new List<MyVector3>();
        List<MyVector3> trajectory2 = new List<MyVector3>();

        //List<List<MyVector3>> trajectory = new List<MyVector3>();
        Mesh objectmesh = null;
        public List<MyVector3> Profile
        {
            get { return this.profile; }
        }
        public List<MyVector3> Trajectory1
        {
            get { return this.trajectory1; }
        }
        public List<MyVector3> Trajectory2
        {
            get { return this.trajectory2; }
        }
        public Mesh Objectmesh
        {
            get { return this.objectmesh; }
        }
        public SweepMesh()
        {
        }
        public SweepMesh(MyCircle profile_, List<MyVector3> trajectory1_, List<MyVector3> trajectory2_)
        {
            this.profile = profile_.CirclePoints;
            this.trajectory1 = trajectory1_;
            this.trajectory2 = trajectory2_;

            //this.CreateMesh();
            this.CreateSymmetryMesh(profile_.Normal, profile_.Radius);
        }

        public SweepMesh(List<MyCircle> CircleList)
        {
            this.CreateCylinderMesh(CircleList);
        }

        private void CreateSymmetryMesh(MyVector3 profilenormal_, double radius_)
        {
            int n = this.profile.Count;   // base stroke point count
            int m = this.trajectory1.Count;   // ref stroke point count


            // sweep the points, duplicate n times and offsetting
            List<double> vertices = new List<double>();
            foreach (MyVector3 points in profile)
            {
                vertices.AddRange(points.ToArray());
            }

            // set base stroke main axis
            MyVector3 dir = profilenormal_.Normalize();
            MyVector3 basecenter3 = new MyVector3();
            for (int i = 0; i < n; i++)
            {
                basecenter3 += this.profile[i];
            }
            basecenter3 /= n;

            int closestidx = -1;
            //3d
            double dis = double.MaxValue;
            for (int i = 0; i < n; i++)
            {
                double d = (this.profile[i] - this.trajectory1.First()).Length();
                if (d < dis)
                {
                    closestidx = i;
                    dis = d;
                }
            }
            double basescale = radius_;


            MyVector3 o3 = basecenter3;
            double tx = 0, ty = 0, tz = 0;
            double sn = 1;

            List<MyVector3> onerow = new List<MyVector3>();
            onerow.AddRange(this.profile);

            for (int i = 0; i < m; i++)
            {
                Line3 normalline = new Line3(o3, dir);
                MyVector3 refp1 = normalline.ProjectToLine(this.trajectory1[i]);
                MyVector3 refp0 = normalline.ProjectToLine(onerow[closestidx]);
                MyVector3 tv = refp1 - refp0;
                tx += tv.x;
                ty += tv.y;
                tz += tv.z;
                o3 += tv;
                double d1 = (this.trajectory1[i] - refp1).Length();
                sn = d1 / basescale;

                double[] new_points = new double[n * 3];//one row
                onerow.Clear();
                for (int j = 0, k = 0; j < n; ++j, k += 3)//move base point
                {
                    new_points[k] = vertices[k] + tx;
                    new_points[k + 1] = vertices[k + 1] + ty;
                    new_points[k + 2] = vertices[k + 2] + tz;
                    MyVector3 tp = new MyVector3(new_points, k);

                    double oplength = (o3 - tp).Length();
                    oplength = oplength * sn;
                    MyVector3 opdir = (tp - o3).Normalize();
                    double t = oplength;

                    new_points[k] = opdir.x * t + o3.x;
                    new_points[k + 1] = opdir.y * t + o3.y;
                    new_points[k + 2] = opdir.z * t + o3.z;

                    onerow.Add(new MyVector3(new_points, k));
                }
                vertices.AddRange(new_points);
            }

            List<int> findices = new List<int>();

            for (int i = 0; i < m; ++i)
            {
                for (int j = 0; j < n - 1; ++j)
                {
                    int s = i * n + j, t = i * n + j + 1;
                    int p = (i + 1) * n + j, q = (i + 1) * n + j + 1;
                    // s-t-p, t-p-q
                    findices.Add(s); findices.Add(p); findices.Add(t);
                    findices.Add(t); findices.Add(p); findices.Add(q);
                }
            }
            this.objectmesh = new Mesh(vertices, findices);
        }

        private void CreateCylinderMesh(List<MyCircle> CircleList)
        {
            int n = CircleList[0].CirclePoints.Count;   // base stroke point count
            int m = CircleList.Count;   // ref stroke point count

            List<double> vertices = new List<double>();
            for (int i = 0; i < CircleList.Count; i++)
            {
                foreach (var point in CircleList[i].CirclePoints)
                {
                    vertices.AddRange(point.ToArray());
                }
            }

            List<int> findices = new List<int>();
            for (int i = 0; i < m-1; ++i)
            {
                for (int j = 0; j < n - 1; ++j)
                {
                    int s = i * n + j, t = i * n + j + 1;
                    int p = (i + 1) * n + j, q = (i + 1) * n + j + 1;
                    // s-p-t, t-p-q
                    findices.Add(s); findices.Add(t); findices.Add(p);
                    findices.Add(t); findices.Add(q); findices.Add(p);
                }
            }
            this.objectmesh = new Mesh(vertices, findices);
        }

        private void CreateMesh()
        {
            int n = this.profile.Count;   // base stroke point count
            int m = this.trajectory1.Count;   // ref stroke point count

            // sweep the points, duplicate n times and offsetting
            List<double> vertices = new List<double>();
            foreach (MyVector3 pos3 in this.profile)
            {
                vertices.AddRange(pos3.ToArray());
            }
            double dx = 0, dy = 0, dz = 0;
            for (int i = 0; i < m - 1; ++i)
            {
                MyVector3 dv = this.trajectory1[i + 1] - this.trajectory1[i];
                dx += dv.x;
                dy += dv.y;
                dz += dv.z;
                double[] new_points = new double[n * 3];//one row
                for (int j = 0, k = 0; j < n; ++j, k += 3)//move base point
                {
                    new_points[k] = vertices[k] + dx;
                    new_points[k + 1] = vertices[k + 1] + dy;
                    new_points[k + 2] = vertices[k + 2] + dz;
                }
                vertices.AddRange(new_points);
            }
            // face index
            List<int> findices = new List<int>();
            for (int i = 0; i < m - 1; ++i)
            {
                for (int j = 0; j < n - 1; ++j)
                {
                    int s = i * n + j, t = i * n + j + 1;
                    int p = (i + 1) * n + j, q = (i + 1) * n + j + 1;
                    // s-t-p, t-p-q
                    findices.Add(s); findices.Add(t); findices.Add(p);
                    findices.Add(t); findices.Add(p); findices.Add(q);
                }
            }

            this.objectmesh = new Mesh(vertices, findices);
        }


        public void Draw(Color c, byte opacity = 120) //0 - white / 1- chosen / 2- fixed chosen / 3- different color
        {
            if (this.objectmesh != null)
            {
                this.objectmesh.DrawSmoothShaded(c);
                //this.objectmesh.DrawTransparent(c);
            }
        }
    }
}
