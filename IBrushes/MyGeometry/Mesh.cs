using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;

using OpenTK.Graphics.OpenGL;
using System.Drawing;

namespace MyGeometry
{
    public unsafe class Mesh
    {
        private int vertexCount;
        private int faceCount;
        private double[] vertexPos = null;
        private double[] vertexNormal = null;
        private int[] faceIndex = null;
        private double[] faceNormal = null;
        private double[] dualVertexPos = null;
        private double[] gaussianCurvature = null;
        private byte[] flag = null;
        private bool[] isBoundary = null;
        private float[] color = null;
        private int[] singleVertexGroup = null;
        private int[][] adjVV = null;
        private int[][] adjVF = null;
        private int[][] adjFF = null;

        public double[] faceSDF = null;

        //// NEW_ADD
        //private int[] segVertexIndex = null;

        public int VertexCount
        {
            get { return vertexCount; }
        }
        public int FaceCount
        {
            get { return faceCount; }
        }
        public double[] VertexPos
        {
            get { return vertexPos; }
            set
            {
                if (value.Length < vertexCount * 3)
                    throw new Exception();
                vertexPos = value;
            }
        }
        public double[] VertexNormal { get { return vertexNormal; } set { vertexNormal = value; } }
        public int[] FaceIndex { get { return faceIndex; } }
        public double[] FaceNormal { get { return faceNormal; } set { faceNormal = value; } }
        public double[] DualVertexPos { get { return dualVertexPos; } }
        public double[] GaussianCurvature { get { return gaussianCurvature; } }
        public byte[] Flag
        {
            get { return flag; }
            set
            {
                if (value.Length < vertexCount)
                    throw new Exception();
                flag = value;
            }
        }
        public bool[] IsBoundary
        {
            get { return isBoundary; }
            set { isBoundary = value; }
        }
        public int[] SingleVertexGroup
        {
            get { return singleVertexGroup; }
        }
        public float[] Color
        {
            get { return color; }
            set { color = value; }
        }
        public int[][] AdjVV
        {
            get { return adjVV; }
            set { adjVV = value; }
        }
        public int[][] AdjVF
        {
            get { return adjVF; }
            set { adjVF = value; }
        }
        public int[][] AdjFF
        {
            get { return adjFF; }
            set { adjFF = value; }
        }
        public Mesh(StreamReader sr)
        {
            ArrayList vlist = new ArrayList();
            ArrayList flist = new ArrayList();
            char[] delimiters = { ' ', '\t' };
            string s = "";

            while (sr.Peek() > -1)
            {
                s = sr.ReadLine();
                string[] tokens = s.Split(delimiters);
                switch (tokens[0].ToLower())
                {
                    case "v":
                        for (int i = 1; i < tokens.Length; i++)
                        {
                            if (tokens[i].Equals("")) continue;
                            vlist.Add(Double.Parse(tokens[i]));
                        }
                        break;
                    case "f":
                        for (int i = 1; i < tokens.Length; i++)
                        {
                            if (tokens[i].Equals("")) continue;
                            string[] tokens2 = tokens[i].Split('/');
                            int index = Int32.Parse(tokens2[0]);
                            if (index <= 0) index = vlist.Count + index + 1;
                            flist.Add(index - 1);
                        }
                        break;
                }
            }

            this.vertexCount = vlist.Count / 3;
            this.faceCount = flist.Count / 3;
            this.vertexPos = new double[vertexCount * 3];
            this.vertexNormal = new double[vertexCount * 3];
            this.faceIndex = new int[faceCount * 3];
            this.faceNormal = new double[faceCount * 3];
            this.flag = new byte[vertexCount];
            this.isBoundary = new bool[vertexCount];
            this.color = new float[faceCount * 3];
            this.gaussianCurvature = new double[vertexCount];

            for (int i = 0; i < vlist.Count; i++) vertexPos[i] = (double)vlist[i];
            for (int i = 0; i < flist.Count; i++) faceIndex[i] = (int)flist[i];
            //			for (int i = 1; i < vlist.Count; i += 3) vertexPos[i] /= 2.0;

            ScaleToUnitBox();
            MoveToCenter();
            ComputeFaceNormal();
            ComputeVertexNormal();
            this.adjVV = BuildAdjacentMatrix().GetRowIndex();
            this.adjVF = BuildAdjacentMatrixFV().GetColumnIndex();
            this.adjFF = BuildAdjacentMatrixFF().GetRowIndex();
            FindBoundaryVertex();
            ComputeDualPosition();

            for (int i = 0; i < FaceCount; i++)
            {
                double area = ComputeFaceArea(i);
                if (double.IsNaN(area))
                    //FormMain.CurrForm.OutputText("bad tri: " + i);
                    if (AdjFF[i].Length != 3)
                    {
                        //FormMain.CurrForm.OutputText("bad FF adj: " + i + " " + AdjFF[i].Length);
                    }
            }
        }

        public Mesh(StreamReader sr, String type, bool no_use)
        {
            if (type.Equals("cgal") == false) throw new Exception();

            char[] delimiters = { ' ', '\t', '\n', '\r' };
            String[] tokens = sr.ReadToEnd().Split(delimiters);

            //FormMain.CurrForm.OutputText(tokens[0]+" "+tokens[1]+" "+tokens[2]);


            this.vertexCount = Int32.Parse(tokens[0]);
            this.faceCount = Int32.Parse(tokens[1]);
            this.vertexPos = new double[vertexCount * 3];
            this.faceIndex = new int[faceCount * 3];

            int k = 3;
            for (int i = 0, j = 3; i < vertexCount - 1; i++, j += 3)
            {
                while (tokens[k].Equals("")) k++;
                vertexPos[j + 0] = Double.Parse(tokens[k++]);
                while (tokens[k].Equals("")) k++;
                vertexPos[j + 1] = Double.Parse(tokens[k++]);
                vertexPos[j + 2] = 0;
            }


            for (int i = 0, j = 0; i < faceCount; i++, j += 3)
            {
                while (tokens[k].Equals("")) k++;
                faceIndex[j + 0] = Int32.Parse(tokens[k++]);
                while (tokens[k].Equals("")) k++;
                faceIndex[j + 1] = Int32.Parse(tokens[k++]);
                while (tokens[k].Equals("")) k++;
                faceIndex[j + 2] = Int32.Parse(tokens[k++]);
            }

            RemoveVertex(0);

            this.vertexNormal = new double[vertexCount * 3];
            this.faceNormal = new double[faceCount * 3];
            this.flag = new byte[vertexCount];
            this.isBoundary = new bool[vertexCount];

            ScaleToUnitBox();
            MoveToCenter();
            ComputeFaceNormal();
            ComputeVertexNormal();
            SparseMatrix adjMatrix = BuildAdjacentMatrix();
            SparseMatrix adjMatrixFV = BuildAdjacentMatrixFV();
            this.adjVV = adjMatrix.GetRowIndex();
            this.adjVF = adjMatrixFV.GetColumnIndex();
            FindBoundaryVertex();
            ComputeDualPosition();
        }
        public Mesh(StreamReader sr, String type)
        {
            if (type.Equals("obj"))
                LoadObjFile(sr);
            else if (type.Equals("off"))
                LoadOffFile(sr);
            FaceIndexReOrder();
        }
        public Mesh(List<double> vertices, List<int> findices)
        {
            this.vertexCount = vertices.Count / 3;
            this.faceCount = findices.Count / 3;
            this.vertexPos = vertices.ToArray();
            this.faceIndex = findices.ToArray();

            this.vertexNormal = new double[vertexCount * 3];
            this.faceNormal = new double[faceCount * 3];
            this.flag = new byte[vertexCount];
            this.isBoundary = new bool[vertexCount];

            ComputeFaceNormal();
            ComputeVertexNormal();
            SparseMatrix adjMatrix = BuildAdjacentMatrix();
            SparseMatrix adjMatrixFV = BuildAdjacentMatrixFV();
            this.adjVV = adjMatrix.GetRowIndex();
            this.adjVF = adjMatrixFV.GetColumnIndex();
            FindBoundaryVertex();
            ComputeDualPosition();
        }
        private void LoadObjFile(StreamReader sr)
        {
            ArrayList vlist = new ArrayList();
            ArrayList flist = new ArrayList();
            char[] delimiters = { ' ', '\t' };
            string s = "";

            while (sr.Peek() > -1)
            {
                s = sr.ReadLine();
                string[] tokens = s.Split(delimiters);
                switch (tokens[0].ToLower())
                {
                    case "v":
                        for (int i = 1; i < tokens.Length; i++)
                        {
                            if (tokens[i].Equals("")) continue;
                            vlist.Add(Double.Parse(tokens[i]));
                        }
                        break;
                    case "f":
                        for (int i = 1; i < tokens.Length; i++)
                        {
                            if (tokens[i].Equals("")) continue;
                            string[] tokens2 = tokens[i].Split('/');
                            int index = Int32.Parse(tokens2[0]);
                            if (index <= 0) index = vlist.Count + index + 1;
                            flist.Add(index - 1);
                        }
                        break;
                }
            }

            this.vertexCount = vlist.Count / 3;
            this.faceCount = flist.Count / 3;
            this.vertexPos = new double[vertexCount * 3];
            this.vertexNormal = new double[vertexCount * 3];
            this.faceIndex = new int[faceCount * 3];
            this.faceNormal = new double[faceCount * 3];
            this.flag = new byte[vertexCount];
            this.isBoundary = new bool[vertexCount];
            this.color = new float[faceCount * 3];
            this.gaussianCurvature = new double[vertexCount];

            for (int i = 0; i < vlist.Count; i++) vertexPos[i] = (double)vlist[i];
            for (int i = 0; i < flist.Count; i++) faceIndex[i] = (int)flist[i];
            //			for (int i = 1; i < vlist.Count; i += 3) vertexPos[i] /= 2.0;

            ScaleToUnitBox();
            MoveToCenter();
            ComputeFaceNormal();
            ComputeVertexNormal();
            this.adjVV = BuildAdjacentMatrix().GetRowIndex();
            this.adjVF = BuildAdjacentMatrixFV().GetColumnIndex();
            this.adjFF = BuildAdjacentMatrixFF().GetRowIndex();
            FindBoundaryVertex();
            ComputeDualPosition();

            for (int i = 0; i < FaceCount; i++)
            {
                double area = ComputeFaceArea(i);
                if (double.IsNaN(area))
                    //FormMain.CurrForm.OutputText("bad tri: " + i);
                    if (AdjFF[i].Length != 3)
                    {
                        //FormMain.CurrForm.OutputText("bad FF adj: " + i + " " + AdjFF[i].Length);
                    }
            }
            //this.faceArea = ComputeFaceArea();
        }
        private void LoadOffFile(StreamReader sr)
        {
            ArrayList vlist = new ArrayList();
            ArrayList flist = new ArrayList();
            char[] delimiters = { ' ', '\t' };
            string s = "";
            int vn = 0; // number of vertex
            int fn = 0; // number of face
            int en = 0; // number of edge

            while (sr.Peek() > -1)
            {
                s = sr.ReadLine();
                if (s.ToUpper().Equals("OFF"))
                {
                    s = sr.ReadLine();
                    string[] tokens = s.Split(delimiters);
                    if (tokens.Length < 3)
                        return;
                    vn = int.Parse(tokens[0]);
                    fn = int.Parse(tokens[1]);
                    en = int.Parse(tokens[2]);
                    break;
                }
            }
            for (int i = 0; i < vn; ++i)
            {
                s = sr.ReadLine();
                string[] tokens = s.Split(delimiters);
                for (int j = 0; j < tokens.Length; ++j)
                {
                    if (tokens[j] != "")
                        vlist.Add(double.Parse(tokens[j]));
                }
            }
            for (int i = 0; i < fn; ++i)
            {
                s = sr.ReadLine();
                string[] tokens = s.Split(delimiters);
                int k = int.Parse(tokens[0]);
                for (int j = 1; j <= k; ++j)
                    flist.Add(int.Parse(tokens[j]));
            }

            this.vertexCount = vlist.Count / 3;
            this.faceCount = flist.Count / 3;
            this.vertexPos = new double[vertexCount * 3];
            this.vertexNormal = new double[vertexCount * 3];
            this.faceIndex = new int[faceCount * 3];
            this.faceNormal = new double[faceCount * 3];
            this.flag = new byte[vertexCount];
            this.isBoundary = new bool[vertexCount];
            this.color = new float[faceCount * 3];
            this.gaussianCurvature = new double[vertexCount];

            for (int i = 0; i < vlist.Count; i++) vertexPos[i] = (double)vlist[i];
            for (int i = 0; i < flist.Count; i++) faceIndex[i] = (int)flist[i];
            //			for (int i = 1; i < vlist.Count; i += 3) vertexPos[i] /= 2.0;

            ScaleToUnitBox();
            MoveToCenter();
            ComputeFaceNormal();
            ComputeVertexNormal();
            this.adjVV = BuildAdjacentMatrix().GetRowIndex();
            this.adjVF = BuildAdjacentMatrixFV().GetColumnIndex();
            this.adjFF = BuildAdjacentMatrixFF().GetRowIndex();
            FindBoundaryVertex();
            ComputeDualPosition();

            for (int i = 0; i < FaceCount; i++)
            {
                double area = ComputeFaceArea(i);
                if (double.IsNaN(area))
                    //FormMain.CurrForm.OutputText("bad tri: " + i);
                    if (AdjFF[i].Length != 3)
                    {
                        //FormMain.CurrForm.OutputText("bad FF adj: " + i + " " + AdjFF[i].Length);
                    }
            }
            //this.faceArea = ComputeFaceArea();
        }
        public bool FaceIndexReOrder() // the input model must be manifold, if manifold, return true, else return false
        {
            if (this.FaceIndex == null)
                return false;

            int fn = this.faceCount;
            int vn = this.vertexCount;

            bool[] tag = new bool[fn];
            for (int i = 0; i < fn; ++i)
                tag[i] = false;

            Set<int>[] tag_e = new Set<int>[vn];
            for (int i = 0; i < vn; ++i)
                tag_e[i] = new Set<int>();

            Queue q = new Queue();
            q.Enqueue(0); // use face[0] as the first element of queue

            while (q.Count > 0)
            {
                int f = (int)q.Dequeue();
                if (tag[f])
                    continue;
                int v1 = this.faceIndex[f * 3];
                int v2 = this.faceIndex[f * 3 + 1];
                int v3 = this.faceIndex[f * 3 + 2];
                if (tag_e[v1].Contains(v2) || tag_e[v2].Contains(v3) || tag_e[v3].Contains(v1))
                {
                    tag_e[v1].Add(v3);
                    tag_e[v2].Add(v1);
                    tag_e[v3].Add(v2);

                    this.faceIndex[f * 3] = v3;
                    this.faceIndex[f * 3 + 1] = v2;
                    this.faceIndex[f * 3 + 2] = v1;
                }
                else
                {
                    tag_e[v1].Add(v2);
                    tag_e[v2].Add(v3);
                    tag_e[v3].Add(v1);

                    this.faceIndex[f * 3] = v1;
                    this.faceIndex[f * 3 + 1] = v2;
                    this.faceIndex[f * 3 + 2] = v3;
                }
                // set tag[f] to be true
                tag[f] = true;

                // add neighbour face
                foreach (int f2 in this.adjFF[f])
                {
                    if (!tag[f2])
                        q.Enqueue(f2);
                }
            }
            this.ComputeFaceNormal();
            this.ComputeVertexNormal();
            return true;
        }
        public void Write(StreamWriter sw)
        {
            for (int i = 0, j = 0; i < vertexCount; i++, j += 3)
            {
                sw.Write("v ");
                sw.Write(vertexPos[j].ToString() + " ");
                sw.Write(vertexPos[j + 1].ToString() + " ");
                sw.Write(vertexPos[j + 2].ToString());
                sw.WriteLine();
            }

            for (int i = 0, j = 0; i < faceCount; i++, j += 3)
            {
                sw.Write("f ");
                sw.Write((faceIndex[j] + 1).ToString() + " ");
                sw.Write((faceIndex[j + 1] + 1).ToString() + " ");
                sw.Write((faceIndex[j + 2] + 1).ToString());
                sw.WriteLine();
            }
        }
        public void LoadSelectedVertexPositions(StreamReader sr)
        {
            ArrayList vlist = new ArrayList();
            ArrayList flist = new ArrayList();
            char[] delimiters = { ' ', '\t' };
            string s = "";

            while (sr.Peek() > -1)
            {
                s = sr.ReadLine();
                string[] tokens = s.Split(delimiters);

                int index = Int32.Parse(tokens[0]);
                double x = Double.Parse(tokens[1]);
                double y = Double.Parse(tokens[2]);
                double z = Double.Parse(tokens[3]);
                if (index >= vertexCount) continue;
                index *= 3;
                this.vertexPos[index++] = x;
                this.vertexPos[index++] = y;
                this.vertexPos[index] = z;
            }
        }
        public void SaveSelectedVertexPositions(StreamWriter sw)
        {
            for (int i = 0, j = 0; i < vertexCount; i++, j += 3)
            {
                if (this.flag[i] == 0) continue;
                sw.Write(i.ToString() + " ");
                sw.Write(vertexPos[j].ToString() + " ");
                sw.Write(vertexPos[j + 1].ToString() + " ");
                sw.Write(vertexPos[j + 2].ToString());
                sw.WriteLine();
            }
        }
        public MyVector3 maxCoord;
        public MyVector3 minCoord;
        public MyVector3 MaxCoord()
        {
            MyVector3 maxCoord = new MyVector3(double.MinValue, double.MinValue, double.MinValue);
            for (int i = 0, j = 0; i < vertexCount; i++, j += 3)
            {
                MyVector3 v = new MyVector3(vertexPos, j);
                maxCoord = MyVector3.Max(maxCoord, v);
            }
            return maxCoord;
        }
        public MyVector3 MinCoord()
        {
            MyVector3 minCoord = new MyVector3(double.MaxValue, double.MaxValue, double.MaxValue);
            for (int i = 0, j = 0; i < vertexCount; i++, j += 3)
            {
                MyVector3 v = new MyVector3(vertexPos, j);
                minCoord = MyVector3.Min(minCoord, v);
            }
            return minCoord;
        }
        public double Volume()
        {
            double totVolume = 0;
            for (int i = 0, j = 0; i < faceCount; i++, j += 3)
            {
                int c1 = faceIndex[j] * 3;
                int c2 = faceIndex[j + 1] * 3;
                int c3 = faceIndex[j + 2] * 3;
                MyVector3 a = new MyVector3(vertexPos, c1);
                MyVector3 b = new MyVector3(vertexPos, c2);
                MyVector3 c = new MyVector3(vertexPos, c3);
                totVolume +=
                    a.x * b.y * c.z -
                    a.x * b.z * c.y -
                    a.y * b.x * c.z +
                    a.y * b.z * c.x +
                    a.z * b.x * c.y -
                    a.z * b.y * c.x;
            }
            return totVolume;
        }
        public void MoveToCenter()
        {
            this.maxCoord = MaxCoord();
            this.minCoord = MinCoord();

            MyVector3 center = (this.maxCoord + this.minCoord) / 2.0;

            for (int i = 0, j = 0; i < vertexCount; i++, j += 3)
            {
                vertexPos[j] -= center.x;
                vertexPos[j + 1] -= center.y;
                vertexPos[j + 2] -= center.z;
            }
        }
        public void ScaleToUnitBox()
        {
            this.maxCoord = MaxCoord();
            this.minCoord = MinCoord();

            MyVector3 d = this.maxCoord - this.minCoord;
            double s = (d.x > d.y) ? d.x : d.y;
            s = (s > d.z) ? s : d.z;
            if (s <= 0) return;
            for (int i = 0; i < vertexPos.Length; i++)
                vertexPos[i] /= s;
        }
        public void Transform(MyMatrix4d tran)
        {
            for (int i = 0, j = 0; i < vertexCount; i++, j += 3)
            {
                MyVector4 v = new MyVector4(vertexPos[j], vertexPos[j + 1], vertexPos[j + 2], 1.0);
                v = tran * v;
                vertexPos[j] = v.x;
                vertexPos[j + 1] = v.y;
                vertexPos[j + 2] = v.z;
            }
        }
        public void ComputeFaceNormal()
        {
            for (int i = 0, j = 0; i < faceCount; i++, j += 3)
            {
                int c1 = faceIndex[j] * 3;
                int c2 = faceIndex[j + 1] * 3;
                int c3 = faceIndex[j + 2] * 3;
                MyVector3 v1 = new MyVector3(vertexPos, c1);
                MyVector3 v2 = new MyVector3(vertexPos, c2);
                MyVector3 v3 = new MyVector3(vertexPos, c3);
                MyVector3 normal = (v2 - v1).Cross(v3 - v1).Normalize();
                faceNormal[j] = normal.x;
                faceNormal[j + 1] = normal.y;
                faceNormal[j + 2] = normal.z;
            }
        }
        public void ComputeVertexNormal()
        {
            Array.Clear(vertexNormal, 0, vertexNormal.Length);
            for (int i = 0, j = 0; i < faceCount; i++, j += 3)
            {
                int c1 = faceIndex[j] * 3;
                int c2 = faceIndex[j + 1] * 3;
                int c3 = faceIndex[j + 2] * 3;
                vertexNormal[c1] += faceNormal[j];
                vertexNormal[c2] += faceNormal[j];
                vertexNormal[c3] += faceNormal[j];
                vertexNormal[c1 + 1] += faceNormal[j + 1];
                vertexNormal[c2 + 1] += faceNormal[j + 1];
                vertexNormal[c3 + 1] += faceNormal[j + 1];
                vertexNormal[c1 + 2] += faceNormal[j + 2];
                vertexNormal[c2 + 2] += faceNormal[j + 2];
                vertexNormal[c3 + 2] += faceNormal[j + 2];
            }
            for (int i = 0, j = 0; i < vertexCount; i++, j += 3)
            {
                MyVector3 n = new MyVector3(vertexNormal, j);
                n = n.Normalize();
                vertexNormal[j] = n.x;
                vertexNormal[j + 1] = n.y;
                vertexNormal[j + 2] = n.z;
            }
        }
        public void ComputeGaussianCurvature()
        {
            int n = this.vertexCount;
            int fn = this.faceCount;
            double[] area = new double[n];

            for (int i = 0; i < n; i++)
            {
                this.gaussianCurvature[i] = 0;
                area[i] = 0;
            }

            for (int i = 0; i < fn; i++)
            {
                int b = i * 3;
                int c1 = this.faceIndex[b];
                int c2 = this.faceIndex[b + 1];
                int c3 = this.faceIndex[b + 2];
                MyVector3 v1 = new MyVector3(this.vertexPos, c1 * 3);
                MyVector3 v2 = new MyVector3(this.vertexPos, c2 * 3);
                MyVector3 v3 = new MyVector3(this.vertexPos, c3 * 3);
                double d1 = (v2 - v3).Length();
                double d2 = (v3 - v1).Length();
                double d3 = (v1 - v2).Length();
                //double a1 = 
            }
        }
        public SparseMatrix BuildAdjacentMatrix()
        {
            SparseMatrix m = new SparseMatrix(vertexCount, vertexCount, 6);

            for (int i = 0, j = 0; i < faceCount; i++, j += 3)
            {
                int c1 = faceIndex[j];
                int c2 = faceIndex[j + 1];
                int c3 = faceIndex[j + 2];
                m.AddElementIfNotExist(c1, c2, 1.0);
                m.AddElementIfNotExist(c2, c3, 1.0);
                m.AddElementIfNotExist(c3, c1, 1.0);
                m.AddElementIfNotExist(c2, c1, 1.0);
                m.AddElementIfNotExist(c3, c2, 1.0);
                m.AddElementIfNotExist(c1, c3, 1.0);
            }

            m.SortElement();
            return m;
        }
        public SparseMatrix BuildAdjacentMatrixFV()
        {
            SparseMatrix m = new SparseMatrix(faceCount, vertexCount, 6);

            for (int i = 0, j = 0; i < faceCount; i++, j += 3)
            {
                m.AddElementIfNotExist(i, faceIndex[j], 1.0);
                m.AddElementIfNotExist(i, faceIndex[j + 1], 1.0);
                m.AddElementIfNotExist(i, faceIndex[j + 2], 1.0);
            }

            m.SortElement();
            return m;
        }
        public SparseMatrix BuildAdjacentMatrixFF()
        {
            SparseMatrix m = new SparseMatrix(faceCount, faceCount, 3);

            for (int i = 0; i < faceCount; i++)
            {
                int v1 = faceIndex[i * 3];
                int v2 = faceIndex[i * 3 + 1];
                int v3 = faceIndex[i * 3 + 2];

                foreach (int j in adjVF[v1])
                    if (j != i && IsContainVertex(j, v2))
                        m.AddElementIfNotExist(i, j, 1.0);

                foreach (int j in adjVF[v2])
                    if (j != i && IsContainVertex(j, v3))
                        m.AddElementIfNotExist(i, j, 1.0);

                foreach (int j in adjVF[v3])
                    if (j != i && IsContainVertex(j, v1))
                        m.AddElementIfNotExist(i, j, 1.0);
            }

            return m;
        }
        public void FindBoundaryVertex()
        {
            for (int i = 0; i < vertexCount; i++)
            {
                int nAdjV = adjVV[i].Length;
                int nAdjF = adjVF[i].Length;
                this.isBoundary[i] = (nAdjV != nAdjF);
                if (nAdjV != nAdjF)
                {
                    //FormMain.CurrForm.OutputText("bad: " + i);
                    this.flag[i] = 1;
                }
            }
        }
        public void GroupingFlags()
        {
            for (int i = 0; i < flag.Length; i++)
                if (flag[i] != 0) flag[i] = 255;

            byte id = 0;
            Queue queue = new Queue();
            List<int> singleVertexGroupList = new List<int>();
            for (int i = 0; i < vertexCount; i++)
                if (flag[i] == 255)
                {
                    id++;
                    flag[i] = id;
                    queue.Enqueue(i);
                    bool found = false;
                    while (queue.Count > 0)
                    {
                        int curr = (int)queue.Dequeue();
                        foreach (int j in adjVV[curr])
                        {
                            if (flag[j] == 255)
                            {
                                flag[j] = id;
                                queue.Enqueue(j);
                                found = true;
                            }
                        }
                    }

                    if (!found) singleVertexGroupList.Add(i);
                }

            this.singleVertexGroup = singleVertexGroupList.ToArray();
        }
        public void FindSingleVertexGroup()
        {
            Set<int> s = new Set<int>();

            for (int i = 0; i < vertexCount; i++)
            {
                if (flag[i] != 0)
                {
                    bool found = false;
                    foreach (int j in adjVV[i])
                        if (flag[j] == flag[i])
                        {
                            found = true;
                            break;
                        }
                    if (!found)
                        s.Add(i);
                }
            }

            int[] arr = s.ToArray();
            Array.Sort(arr);
            this.singleVertexGroup = arr;
        }
        public void RemoveVertex(int index)
        {
            _RemoveVertex(index);

            this.vertexNormal = new double[vertexCount * 3];
            this.faceNormal = new double[faceCount * 3];
            this.flag = new byte[vertexCount];
            this.isBoundary = new bool[vertexCount];

            ComputeFaceNormal();
            ComputeVertexNormal();
            SparseMatrix adjMatrix = BuildAdjacentMatrix();
            SparseMatrix adjMatrixFV = BuildAdjacentMatrixFV();
            this.adjVV = adjMatrix.GetRowIndex();
            this.adjVF = adjMatrixFV.GetColumnIndex();
            FindBoundaryVertex();
            ComputeDualPosition();
        }
        public void RemoveVertex(ArrayList indice)
        {
            for (int i = 0; i < indice.Count; i++)
                _RemoveVertex(((int)indice[i]) - i);

            this.vertexNormal = new double[vertexCount * 3];
            this.faceNormal = new double[faceCount * 3];
            this.flag = new byte[vertexCount];
            this.isBoundary = new bool[vertexCount];

            ComputeFaceNormal();
            ComputeVertexNormal();
            SparseMatrix adjMatrix = BuildAdjacentMatrix();
            SparseMatrix adjMatrixFV = BuildAdjacentMatrixFV();
            this.adjVV = adjMatrix.GetRowIndex();
            this.adjVF = adjMatrixFV.GetColumnIndex();
            FindBoundaryVertex();
            ComputeDualPosition();
        }
        public bool IsContainVertex(int fIndex, int vIndex)
        {
            int b = fIndex * 3;
            int v1 = faceIndex[b];
            int v2 = faceIndex[b + 1];
            int v3 = faceIndex[b + 2];
            return (v1 == vIndex) || (v2 == vIndex) || (v3 == vIndex);
        }
        public void ComputeDualPosition()
        {
            if (dualVertexPos == null)
                dualVertexPos = new double[faceCount * 3];

            for (int i = 0; i < dualVertexPos.Length; i++)
                dualVertexPos[i] = 0.0;

            for (int i = 0, j = 0; i < vertexCount; i++, j += 3)
                foreach (int k in adjVF[i])
                {
                    int b = k * 3;
                    dualVertexPos[b] += vertexPos[j];
                    dualVertexPos[b + 1] += vertexPos[j + 1];
                    dualVertexPos[b + 2] += vertexPos[j + 2];
                }

            for (int i = 0; i < dualVertexPos.Length; i++)
                dualVertexPos[i] /= 3.0;
        }
        public MyVector3 ComputeDualPosition(int fIndex)
        {
            return new MyVector3(dualVertexPos, fIndex * 3);
        }
        public double ComputeFaceArea()
        {
            double tot = 0;
            for (int i = 0, j = 0; i < faceCount; i++, j += 3)
            {
                MyVector3 v1 = new MyVector3(VertexPos, faceIndex[j] * 3);
                MyVector3 v2 = new MyVector3(VertexPos, faceIndex[j + 1] * 3);
                MyVector3 v3 = new MyVector3(VertexPos, faceIndex[j + 2] * 3);
                tot += ((v2 - v1).Cross(v3 - v1)).Length() / 2.0;
            }
            return tot;
        }
        public double ComputeFaceArea(int fIndex)
        {
            int b = fIndex * 3;
            MyVector3 v1 = new MyVector3(VertexPos, faceIndex[b] * 3);
            MyVector3 v2 = new MyVector3(VertexPos, faceIndex[b + 1] * 3);
            MyVector3 v3 = new MyVector3(VertexPos, faceIndex[b + 2] * 3);
            return ((v2 - v1).Cross(v3 - v1)).Length() / 2.0;
        }
        public double ComputeFaceEdgeLength(int fIndex)
        {
            int b = fIndex * 3;
            MyVector3 v1 = new MyVector3(VertexPos, faceIndex[b] * 3);
            MyVector3 v2 = new MyVector3(VertexPos, faceIndex[b + 1] * 3);
            MyVector3 v3 = new MyVector3(VertexPos, faceIndex[b + 2] * 3);
            return (v1 - v2).Length() + (v2 - v3).Length() + (v3 - v1).Length();
        }
        public double AverageFaceArea()
        {
            double tot = 0;
            for (int i = 0, j = 0; i < faceCount; i++, j += 3)
            {
                MyVector3 v1 = new MyVector3(VertexPos, faceIndex[j] * 3);
                MyVector3 v2 = new MyVector3(VertexPos, faceIndex[j + 1] * 3);
                MyVector3 v3 = new MyVector3(VertexPos, faceIndex[j + 2] * 3);
                tot += ((v2 - v1).Cross(v3 - v1)).Length() / 2.0;
            }
            return (tot / faceCount);
        }
        public void SwapFlags(byte n1, byte n2)
        {
            for (int i = 0; i < vertexCount; i++)
            {
                if (flag[i] == n1) flag[i] = n2;
                else if (flag[i] == n2) flag[i] = n1;
            }
        }
        public void RandomSwapFlags()
        {
            int maxFlag = -1;
            for (int i = 0; i < vertexCount; i++)
                if (flag[i] > maxFlag)
                    maxFlag = flag[i];

            Random rand = new Random(System.DateTime.Now.Second);
            int[] map = new int[maxFlag + 1];
            for (int i = 0; i < map.Length; i++)
            {
                map[i] = i;
            }
            for (int i = 0; i < map.Length; i++)
            {
                int s = rand.Next(maxFlag + 1);
                int tmp = map[i];
                map[i] = map[s];
                map[s] = tmp;
            }

            for (int i = 0; i < vertexCount; i++)
                flag[i] = (byte)map[flag[i]];
        }
        private void _RemoveVertex(int index)
        {
            double[] vp = new double[(vertexCount - 1) * 3];
            for (int i = 0, j = 0, k = 0; i < vertexCount; i++, j += 3)
            {
                if (i == index) continue;
                vp[k++] = vertexPos[j];
                vp[k++] = vertexPos[j + 1];
                vp[k++] = vertexPos[j + 2];
            }

            ArrayList flist = new ArrayList(faceCount * 3);
            for (int i = 0, j = 0; i < faceCount; i++, j += 3)
            {
                int c1 = faceIndex[j];
                int c2 = faceIndex[j + 1];
                int c3 = faceIndex[j + 2];
                if (c1 == index || c2 == index || c3 == index) continue;
                if (c1 > index) c1--; flist.Add(c1);
                if (c2 > index) c2--; flist.Add(c2);
                if (c3 > index) c3--; flist.Add(c3);
            }

            this.vertexCount--;
            this.vertexPos = vp;
            this.faceCount = flist.Count / 3;
            this.faceIndex = new int[flist.Count];

            for (int i = 0; i < flist.Count; i++) faceIndex[i] = (int)flist[i];
        }
        private void _RemoveVertex(ArrayList indice)
        {

        }

        public void DrawFlatShaded(Color c)
        {
            GL.ShadeModel(ShadingModel.Flat);
            GL.PolygonMode(MaterialFace.FrontAndBack, PolygonMode.Fill);
            GL.Enable(EnableCap.Lighting);
            GL.Enable(EnableCap.Normalize);
            GL.Color3(c.R, c.G, c.B);
            fixed (double* vp = this.VertexPos)
            fixed (double* np = this.FaceNormal)
            fixed (int* fi = this.FaceIndex)
            {
                GL.Begin(BeginMode.Triangles);
                for (int i = 0, j = 0; i < this.FaceCount; i++, j += 3)
                {
                    GL.Normal3(np + j);
                    GL.Vertex3(vp + fi[j] * 3);
                    GL.Vertex3(vp + fi[j + 1] * 3);
                    GL.Vertex3(vp + fi[j + 2] * 3);
                }
                GL.End();
            }
            GL.Disable(EnableCap.Lighting);
        }
        public void DrawWireFrame(Color c)
        {

        }
        public void DrawSmoothShaded(Color c)
        {
            GL.Disable(EnableCap.Blend);
            GL.Enable(EnableCap.DepthTest);
            GL.Clear(ClearBufferMask.DepthBufferBit);
            GL.ShadeModel(ShadingModel.Smooth);
            GL.Enable(EnableCap.ColorMaterial);
            GL.PolygonMode(MaterialFace.FrontAndBack, PolygonMode.Fill);
            GL.Enable(EnableCap.Lighting);
            GL.Enable(EnableCap.Normalize);

            GL.Color3(c);
            GL.Begin(PrimitiveType.Triangles);
            for (int i = 0, j = 0; i < this.faceCount; i++, j += 3)
            {
                //GL.Normal3(faceNormal[j], faceNormal[j + 1], faceNormal[j + 2]);
                GL.Normal3(VertexNormal[faceIndex[j] * 3], VertexNormal[faceIndex[j] * 3 + 1], VertexNormal[faceIndex[j] * 3 + 2]);
                GL.Vertex3(VertexPos[faceIndex[j] * 3], VertexPos[faceIndex[j] * 3 + 1], VertexPos[faceIndex[j] * 3 + 2]);
                
                //GL.Normal3(faceNormal[j], faceNormal[j + 1], faceNormal[j + 2]);
                GL.Normal3(VertexNormal[faceIndex[j + 1] * 3], VertexNormal[faceIndex[j + 1] * 3 + 1], VertexNormal[faceIndex[j + 1] * 3 + 2]);
                GL.Vertex3(VertexPos[faceIndex[j + 1] * 3], VertexPos[faceIndex[j + 1] * 3 + 1], VertexPos[faceIndex[j + 1] * 3 + 2]);
                
                //GL.Normal3(faceNormal[j], faceNormal[j + 1], faceNormal[j + 2]);
                GL.Normal3(VertexNormal[faceIndex[j + 2] * 3], VertexNormal[faceIndex[j + 2] * 3 + 1], VertexNormal[faceIndex[j + 2] * 3 + 2]);
                GL.Vertex3(VertexPos[faceIndex[j + 2] * 3], VertexPos[faceIndex[j + 2] * 3 + 1], VertexPos[faceIndex[j + 2] * 3 + 2]);
            }
            GL.End();
            GL.Disable(EnableCap.Lighting);
            GL.Disable(EnableCap.Normalize);
            GL.Disable(EnableCap.DepthTest);
            GL.Enable(EnableCap.Blend);

        }
        public void DrawTransparent(Color c, byte opacity = 92)
        {
            GL.Enable(EnableCap.Blend);
            GL.BlendFunc(BlendingFactorSrc.SrcAlpha, BlendingFactorDest.OneMinusSrcAlpha);
            //	GL.PolygonMode(MaterialFace.FrontAndBack, PolygonMode.Fill);
            GL.Disable(EnableCap.Lighting);
            GL.Disable(EnableCap.Normalize);
            //GL.Color3(c.R, c.G, c.B);


            for (int i = 0, j = 0; i < this.faceCount; i++, j += 3)
            {
                GL.Color4(c.R, c.G, c.B, (byte)opacity);
                GL.Begin(PrimitiveType.Triangles);
                //GL.Normal3(VertexPos[faceIndex[j] * 3], VertexPos[faceIndex[j] * 3 + 1], VertexPos[faceIndex[j] * 3 + 2]);
                GL.Vertex3(VertexPos[faceIndex[j] * 3], VertexPos[faceIndex[j] * 3 + 1], VertexPos[faceIndex[j] * 3 + 2]);
                GL.Vertex3(VertexPos[faceIndex[j + 1] * 3], VertexPos[faceIndex[j + 1] * 3 + 1], VertexPos[faceIndex[j + 1] * 3 + 2]);
                GL.Vertex3(VertexPos[faceIndex[j + 2] * 3], VertexPos[faceIndex[j + 2] * 3 + 1], VertexPos[faceIndex[j + 2] * 3 + 2]);
                GL.End();
            }
        }
    }
}
