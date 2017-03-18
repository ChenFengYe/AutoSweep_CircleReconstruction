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

using Emgu.CV;
using Emgu.CV.CvEnum;
using Emgu.CV.Structure;
using Emgu.CV.UI;
using Emgu.CV.Util;

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
        public void LoadPlane(string planefile)
        {
            this.PlaneList.Add(new MyPlane(planefile));
            this.optplane = this.PlaneList.First();
        }
        public void SaveAllPlane(string commonprefix)
        {
            for (int i = 0; i < PlaneList.Count; i++)
            {
                PlaneList[i].SaveMyPlane2(commonprefix + "_" + (i + 1));
            }
        }
        public void LoadModel(string commonprefix)
        {
            this.ModelList.Add(new MyModel(commonprefix));
        }
        public void SaveAllModel(string commonprefix)
        {
            for (int i = 0; i < ModelList.Count; i++)
            {
                ModelList[i].SaveModelToFile(commonprefix + "_" + (i + 1));
                ModelList[i].SaveModelToPly(commonprefix + "_" + (i + 1));
            }
        }
        public void LoadOutline(string commonprefix)
        {
            using (StreamReader sr = new StreamReader(commonprefix, Encoding.Default))
            {
                int outline_num = int.Parse(sr.ReadLine());
                string[] arr = null;
                for (int i = 0; i < outline_num; i++)
                {
                    List<CLineSegment> outline = new List<CLineSegment>();
                    int segment_num = int.Parse(sr.ReadLine());
                    for (int j = 0; j < segment_num; j++)
                    {
                        arr = sr.ReadLine().Split(new char[] { ' ' });
                        MyVector2 s = new MyVector2(Convert.ToDouble(arr[0]), Convert.ToDouble(arr[1]));
                        MyVector2 t = new MyVector2(Convert.ToDouble(arr[2]), Convert.ToDouble(arr[3]));
                        outline.Add(new CLineSegment(s, t));
                    }
                    this.OutlineList.Add(outline);
                }
            }
            this.ShowOutline();
        }
        public void SaveAllOutline(string commonprefix)
        {
            for (int i = 0; i < OutlineList.Count; i++)
            {
                this.SaveOutline(commonprefix + "_" + (i + 1));
            }
        }
        private void SaveOutline(string commonprefix)
        {
            using (StreamWriter sw = new StreamWriter(commonprefix + ".outline", false))
            {
                sw.WriteLine(this.OutlineList.Count);
                foreach (List<CLineSegment> list in this.OutlineList)
                {
                    sw.WriteLine(list.Count);
                    foreach (CLineSegment cls in list)
                    {
                        sw.WriteLine("{0} {1}", cls.StartPoint, cls.EndPoint);
                    }
                }
                sw.Flush();
                sw.Close();
            }
        }
    }
}
