using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Diagnostics;

using MyGeometry;

using Emgu.CV.Structure;
using Emgu.CV;
using Emgu.CV.UI;

using Accord.Math;
using Accord.Math.Distances;
using Accord.MachineLearning;
using Accord.Statistics.Distributions.DensityKernels;

using GraphClass;

namespace SmartCanvas
{
    // store each plane's information
    public class pointPlaneClass
    {
        public pointPlaneClass(int idx, int value = 0)
        {
            this.idx = idx;
            this.value = value;
            points = new List<MyVector3>();
            pointsIdx = new List<int>();
        }

        int idx;
        int value;
        List<MyVector3> points;
        List<int> pointsIdx;

        MyVector3 normalABC;
        double d;

        public int Idx { get { return idx; } set { this.idx = value; } }
        public int Value { get { return value; } set { this.value = value; } }
        public List<MyVector3> Points { get { return points; } set { this.points = value; } }
        public List<int> PointsIdx { get { return pointsIdx; } set { this.pointsIdx = value; } }
        public MyVector3 NormalABC { get { return normalABC; } set { this.normalABC = value; } }
        public double D { get { return d; } set { this.d = value; } }
    }

    // customized distance function
    public class myDistanceClass : IMetric<double[]>
    {
        // data dimension should be 9: normalizedXYZ normals xyz
        public double Distance(double[] x, double[] y)
        {
            // xyz distance
            double dx = x[0] - y[0];
            double dy = x[1] - y[1];
            double dz = x[2] - y[2];

            // normal distance
            //MyVector3 v1 = new MyVector3(x[3], x[4], x[5]);
            //MyVector3 v2 = new MyVector3(y[3], y[4], y[5]);
            //double dn = 1 - (v1.Dot(v2) + 1) / 2;
            double d3 = x[3] - y[3];
            double d4 = x[4] - y[4];
            double d5 = x[5] - y[5];
            
            return (dx * dx + dy * dy + dz * dz) * 9 + (d3 * d3 + d4 * d4 + d5 * d5) * 1;
        }
    }

    public class pointToPlaneClustering
    {
        private int dataDimension;
        private double msSearchRadius;
        private int w, h;

        public int[] pointLabels;

        public List<pointPlaneClass> clusteringPlaneRec = null;
        public List<pointPlaneClass> extractionPlaneRec = null;
        public List<pointPlaneClass> mergedPlaneRec = null;


        public pointToPlaneClustering(int dataDimension, double radius, int dataWidth, int dataHeight)
        {
            this.dataDimension = dataDimension;
            this.msSearchRadius = radius;
            this.w = dataWidth;
            this.h = dataHeight;
        }

        public void RunProcess(double[][] inputDataMS, bool displayResult = false)
        {
            Stopwatch sw = new Stopwatch();
            sw.Start();

            MeanShift clusterMS = new MeanShift(dataDimension, new UniformKernel(), msSearchRadius);
            clusterMS.Distance = new myDistanceClass();

            MeanShiftClusterCollection clustering = clusterMS.Learn(inputDataMS);
            pointLabels = clustering.Decide(inputDataMS);

            clusteringPlaneRec = new List<pointPlaneClass>();
            for (int i = 0; i < clustering.Count; i++)
                clusteringPlaneRec.Add(new pointPlaneClass(i, 0));
            for (int i = 0; i < h * w; i++)
            {
                MyVector3 vector3T = new MyVector3(inputDataMS[i][6], inputDataMS[i][7], inputDataMS[i][8]);
                if (vector3T.x == 0 && vector3T.y == 0 && vector3T.z == 0)
                    continue;

                int idx = pointLabels[i];
                clusteringPlaneRec[idx].Points.Add(vector3T);
                clusteringPlaneRec[idx].PointsIdx.Add(i);
                clusteringPlaneRec[idx].Value++;
            }
            clusteringPlaneRec.Sort((x, y) => y.Value.CompareTo(x.Value));
            #region visualization
            if (displayResult)
            {
                int loop = 0;
                Image<Bgr, byte> image2 = new Image<Bgr, byte>(w, h);
                image2.SetZero();
                for (int i = 0; i < h; i++)
                {
                    for (int j = 0; j < w; j++)
                    {
                        if (pointLabels[loop] >= 0)
                        {
                            byte r = (byte)(Utils.ColorMall[pointLabels[loop] % 30].R);
                            byte g = (byte)(Utils.ColorMall[pointLabels[loop] % 30].G);
                            byte b = (byte)(Utils.ColorMall[pointLabels[loop] % 30].B);
                            image2[i, j] = new Bgr(b, g, r);
                        }
                        loop++;
                    }
                }
                new ImageViewer(image2, "2 - MeanShiftClustering").Show();
            }
            #endregion
            sw.Stop();
            Console.WriteLine(clusteringPlaneRec.Count + " labels\tin" + sw.ElapsedMilliseconds / 1000 + "s");
            sw.Restart();

            // extract planes from clustered data
            SceondPlaneExtraction();
            #region visualization
            if (displayResult)
            {
                int loop = 0;
                Image<Bgr, byte> image3 = new Image<Bgr, byte>(w, h);
                image3.SetZero();
                for (int i = 0; i < h; i++)
                {
                    for (int j = 0; j < w; j++)
                    {
                        if (pointLabels[loop] >= 0)
                        {
                            byte r = (byte)(Utils.ColorMall[pointLabels[loop] % 30].R);
                            byte g = (byte)(Utils.ColorMall[pointLabels[loop] % 30].G);
                            byte b = (byte)(Utils.ColorMall[pointLabels[loop] % 30].B);
                            image3[i, j] = new Bgr(b, g, r);
                        }
                        loop++;
                    }
                }
                new ImageViewer(image3, "3 - PlaneExtraction").Show();
            }
            #endregion
            sw.Stop();
            Console.WriteLine(extractionPlaneRec.Count + " labels\tin" + sw.ElapsedMilliseconds / 1000 + "s");
            sw.Restart();

            // merge planes if necessary
            MergePlanes();
            #region visualization
            if (displayResult)
            {
                int loop = 0;
                Image<Bgr, byte> image4 = new Image<Bgr, byte>(w, h);
                image4.SetZero();
                for (int i = 0; i < h; i++)
                {
                    for (int j = 0; j < w; j++)
                    {
                        if (pointLabels[loop] >= 0)
                        {
                            byte r = (byte)(Utils.ColorMall[pointLabels[loop] % 30].R);
                            byte g = (byte)(Utils.ColorMall[pointLabels[loop] % 30].G);
                            byte b = (byte)(Utils.ColorMall[pointLabels[loop] % 30].B);
                            image4[i, j] = new Bgr(b, g, r);
                        }
                        loop++;
                    }
                }
                new ImageViewer(image4, "4 - MergedPlanes").Show();
            }
            #endregion

            sw.Stop();
            Console.WriteLine(mergedPlaneRec.Count + " labels\tin" + sw.ElapsedMilliseconds / 1000 + "s");
        }

        public void SceondPlaneExtraction()
        {
            if (clusteringPlaneRec == null || clusteringPlaneRec.Count == 0)
                return;

            bool[] labelsFlag = new bool[pointLabels.Length];
            for (int i = 0; i < pointLabels.Length; i++)
                labelsFlag[i] = false;

            extractionPlaneRec = new List<pointPlaneClass>();
            for (int i = 0; i < 20; i++)
            {
                if (clusteringPlaneRec[0].Value < w * h * 0.0005) // too few points
                    break;

                List<int> validPointIdx = new List<int>();
                MyVector4 returnedValue = RANSACNormal(clusteringPlaneRec[0].Points, validPointIdx, 0.008, 1.0, 200);

                extractionPlaneRec.Add(new pointPlaneClass(i, validPointIdx.Count));
                extractionPlaneRec[i].NormalABC = new MyVector3(returnedValue.x, returnedValue.y, returnedValue.z);
                extractionPlaneRec[i].D = returnedValue.w;
                for (int j = validPointIdx.Count - 1; j >= 0; j--)
                {
                    MyVector3 mv3T = clusteringPlaneRec[0].Points[validPointIdx[j]];
                    extractionPlaneRec[i].Points.Add(mv3T);
                    clusteringPlaneRec[0].Points.RemoveAt(validPointIdx[j]);

                    extractionPlaneRec[i].PointsIdx.Add(clusteringPlaneRec[0].PointsIdx[validPointIdx[j]]);
                    labelsFlag[clusteringPlaneRec[0].PointsIdx[validPointIdx[j]]] = true;
                    pointLabels[clusteringPlaneRec[0].PointsIdx[validPointIdx[j]]] = i;
                    clusteringPlaneRec[0].PointsIdx.RemoveAt(validPointIdx[j]);

                    clusteringPlaneRec[0].Value--;
                }
                // re sort
                if (clusteringPlaneRec.Count > 1)
                {
                    for (int j = 1; j < clusteringPlaneRec.Count; j++)
                    {
                        if (clusteringPlaneRec[j - 1].Value >= clusteringPlaneRec[j].Value)
                            break;
                        pointPlaneClass lcT = clusteringPlaneRec[j - 1];
                        clusteringPlaneRec[j - 1] = clusteringPlaneRec[j];
                        clusteringPlaneRec[j] = lcT;
                    }
                }
            }

            for (int i = 0; i < pointLabels.Length; i++)
                if (labelsFlag[i] == false)
                    pointLabels[i] = -1;
        }

        public void MergePlanes()
        {
            GraphAdjList planeGraph = new GraphAdjList(extractionPlaneRec.Count);
            for (int i = 0; i < extractionPlaneRec.Count; i++)
            {
                for (int j = i + 1; j < extractionPlaneRec.Count; j++)
                {
                    double dis1 = 0;
                    for (int loop = 0; loop < extractionPlaneRec[i].Value; loop++)
                    {
                        dis1 += Math.Abs(extractionPlaneRec[i].Points[loop].x * extractionPlaneRec[j].NormalABC.x +
                            extractionPlaneRec[i].Points[loop].y * extractionPlaneRec[j].NormalABC.y +
                            extractionPlaneRec[i].Points[loop].z * extractionPlaneRec[j].NormalABC.z -
                            extractionPlaneRec[j].D) / extractionPlaneRec[j].D;
                    }
                    dis1 /= extractionPlaneRec[i].Value;

                    double dis2 = 0;
                    for (int loop = 0; loop < extractionPlaneRec[j].Value; loop++)
                    {
                        dis2 += Math.Abs(extractionPlaneRec[j].Points[loop].x * extractionPlaneRec[i].NormalABC.x +
                            extractionPlaneRec[j].Points[loop].y * extractionPlaneRec[i].NormalABC.y +
                            extractionPlaneRec[j].Points[loop].z * extractionPlaneRec[i].NormalABC.z -
                            extractionPlaneRec[i].D) / extractionPlaneRec[i].D;
                    }
                    dis2 /= extractionPlaneRec[j].Value;

                    double dis = (dis1 * extractionPlaneRec[j].Value + dis2 * extractionPlaneRec[i].Value) /
                        (extractionPlaneRec[i].Value + extractionPlaneRec[j].Value);

                    if (extractionPlaneRec[i].NormalABC.Dot(extractionPlaneRec[j].NormalABC) > Math.Cos(10 * Math.PI / 180)
                        && dis < 0.04)
                    {
                        planeGraph.AddEdge(i, j);
                    }
                }
            }

            int groupIdx = 0;
            int[] clusterInd = new int[extractionPlaneRec.Count];
            bool[] visited = new bool[extractionPlaneRec.Count];
            DepthFirstSearch dfs = new DepthFirstSearch();
            for (int i = 0; i < extractionPlaneRec.Count; i++)
            {
                if (!visited[i])
                {
                    dfs.DFS(planeGraph, ref visited, i, ref clusterInd, groupIdx);
                    groupIdx++;
                }
            }

            mergedPlaneRec = new List<pointPlaneClass>();
            for (int i = 0; i < groupIdx; i++)
                mergedPlaneRec.Add(new pointPlaneClass(i, 0));
            for (int i = 0; i < extractionPlaneRec.Count; i++)
            {
                mergedPlaneRec[clusterInd[i]].Value += extractionPlaneRec[i].Value;
                mergedPlaneRec[clusterInd[i]].Points.AddRange(extractionPlaneRec[i].Points);
                mergedPlaneRec[clusterInd[i]].PointsIdx.AddRange(extractionPlaneRec[i].PointsIdx);

                foreach (int pi in extractionPlaneRec[i].PointsIdx)
                    pointLabels[pi] = clusterInd[i];
            }
            for (int i = 0; i < mergedPlaneRec.Count; i++)
            {
                MyVector4 returnedvalue = pointToPlaneClustering.RANSACNormal(mergedPlaneRec[i].Points, null, 0.0002, 0.9, 100);
                mergedPlaneRec[i].NormalABC = new MyVector3(returnedvalue.x, returnedvalue.y, returnedvalue.z);
                mergedPlaneRec[i].D = returnedvalue.w;
            }
        }

        // save plane for visualization, may be called after RunProcess()
        public List<MyPlane> savePlaneForVisualization()
        {
            List<MyPlane> returnedPlane = new List<MyPlane>();
            foreach (pointPlaneClass ppc in mergedPlaneRec)
            {
                MyPlane plane = new MyPlane();
                plane.planeVertices = ppc.Points;
                plane.planeEquation = new Plane((float)ppc.NormalABC.x, (float)ppc.NormalABC.y, (float)ppc.NormalABC.z, (float)ppc.D);

                plane.ComputeBoundQuad();
                returnedPlane.Add(plane);
            }

            return returnedPlane;
        }
        
        public static MyVector4 RANSACNormal(List<MyVector3> input, List<int> validPointIdx = null, double thres = 0.0005, double probability = 0.99, int iterationTimes = 50)
        /* Parameter explaination:
         * thres: if the distance between points and plane are below this threshold, we regard points to be in this plane;
         * probability: the condition when we should stop current iteration, can be larger than 1 */
        {
            Random rd = new Random();
            double[,] A = new double[3, 3];
            double[] B = new double[3];
            B[0] = 1; B[1] = 1; B[2] = 1;

            MyVector3 NormalABC = new MyVector3(0, 0, 0);
            double D = 0;

            double maxValue = 0, value1 = 0, value2 = 0;
            List<int> validPointIdxT = new List<int>();
            for (int i = 0; i < iterationTimes; i++) // the maximium iteration times
            {
                int a = rd.Next(input.Count);
                int b = rd.Next(input.Count);
                int c = rd.Next(input.Count);
                A[0, 0] = input[a].x; A[0, 1] = input[a].y; A[0, 2] = input[a].z;
                A[1, 0] = input[b].x; A[1, 1] = input[b].y; A[1, 2] = input[b].z;
                A[2, 0] = input[c].x; A[2, 1] = input[c].y; A[2, 2] = input[c].z;

                if (Math.Abs(A.Determinant()) < 1e-15) // if noninvertible
                    continue;

                double[] x = A.Inverse().Dot(B);
                int num = 0;
                double averageError = 0;
                validPointIdxT.Clear();
                double xLength = Math.Sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
                for (int j = 0; j < input.Count; j++)
                {
                    double errorT = Math.Abs(input[j].x * x[0] + input[j].y * x[1] + input[j].z * x[2] - 1); // / xLength;
                    if (errorT < thres) // a threshold
                    {
                        num++;
                        averageError += errorT;
                        if (validPointIdx != null)
                            validPointIdxT.Add(j);
                    }
                }
                value1 = (double)num / input.Count;
                value2 = (1 - (averageError / input.Count / thres)) * 0;
                if (maxValue < value1 + value2)
                {
                    maxValue = value1 + value2;
                    NormalABC.x = x[0];
                    NormalABC.y = x[1];
                    NormalABC.z = x[2];
                    NormalABC = NormalABC.Normalize();
                    D = 1.0 / xLength;
                    if (validPointIdx != null)
                    {
                        validPointIdx.Clear();
                        foreach (int vpi in validPointIdxT)
                            validPointIdx.Add(vpi);
                    }
                }
                if (maxValue > probability)
                    break;
            }

            return new MyVector4(NormalABC, D);
        }
        
        public List<double> changeColor()
        {
            List<double> point_colors = new List<double>();
            for (int i = 0; i < h * w; i++)
            {
                if (pointLabels[i] == -1)
                    point_colors.AddRange(new double[3] { 0, 0, 0 });
                else
                    point_colors.AddRange(new double[3] { Utils.ColorMall[pointLabels[i] % 30].R / 255.0,
                            Utils.ColorMall[pointLabels[i] % 30].G / 255.0, Utils.ColorMall[pointLabels[i] % 30].B / 255.0});
            }

            return point_colors;
        }
    }


    // some static compute functions
    public class compute
    {
        public static void normalizeXYZ(ref double[][] inputDataMS, double dMin, double dMax, int w, int h)
        {
            dMax -= dMin;
            int loop = 0;
            for (int i = 0; i < h; i++)
            {
                for (int j = 0; j < w; j++)
                {
                    // change normalizeXYZ to 0~1
                    inputDataMS[loop][0] = (inputDataMS[loop][0] - dMin) / dMax;
                    inputDataMS[loop][1] = (inputDataMS[loop][1] - dMin) / dMax;
                    inputDataMS[loop][2] = (inputDataMS[loop][2] - dMin) / dMax;
                    loop++;
                }
            }
        }

        public static void computeVertexNormal(int areaWidth, ref List<MyVector3> points, ref double[][] inputDataMS, int w, int h, bool displayResult = false)
        {
            List<double> pointNormals = new List<double>();
            int loop = 0;
            for (int i = 0; i < h; i++)
            {
                for (int j = 0; j < w; j++)
                {
                    List<MyVector3> candidatePoints = new List<MyVector3>();
                    for (int i1 = -areaWidth; i1 <= areaWidth; i1++)
                    {
                        for (int j1 = -areaWidth; j1 <= areaWidth; j1++)
                        {
                            int iT = i + i1, jT = j + j1;
                            if (iT >= 0 && iT < h && jT >= 0 && jT < w)
                                if (!(points[iT * w + jT].x == 0 && points[iT * w + jT].y == 0 && points[iT * w + jT].z == 0))
                                    candidatePoints.Add(points[iT * w + jT]);
                        }
                    }
                    MyVector4 planeABCD;
                    if (candidatePoints.Count < areaWidth * areaWidth)
                        planeABCD = new MyVector4(0, 0, 0, 0);
                    else
                    {
                        planeABCD = pointToPlaneClustering.RANSACNormal(candidatePoints);
                    }
                    inputDataMS[loop][3] = planeABCD.x;
                    inputDataMS[loop][4] = planeABCD.y;
                    inputDataMS[loop][5] = planeABCD.z;
                    loop++;

                    pointNormals.AddRange(new double[3] { planeABCD.x, planeABCD.y, planeABCD.z });
                }
            }
            if (displayResult)
            {
                loop = 0;
                Image<Bgr, byte> image1 = new Image<Bgr, byte>(w, h);
                image1.SetZero();
                for (int i = 0; i < h; i++)
                {
                    for (int j = 0; j < w; j++)
                    {
                        double rd = pointNormals[3 * loop] > 0 ? pointNormals[3 * loop] : -pointNormals[3 * loop] * 0;
                        double gd = pointNormals[3 * loop + 1] > 0 ? pointNormals[3 * loop + 1] : -pointNormals[3 * loop + 1] * 0;
                        double bd = pointNormals[3 * loop + 2] > 0 ? pointNormals[3 * loop + 2] : -pointNormals[3 * loop + 2] * 0;

                        byte r = (byte)(rd * 255);
                        byte g = (byte)(gd * 255);
                        byte b = (byte)(bd * 255);
                        image1[i, j] = new Bgr(b, g, r);
                        loop++;
                    }
                }
                new ImageViewer(image1, "1 - PointNormal").Show();
            }
        }
    }
}
