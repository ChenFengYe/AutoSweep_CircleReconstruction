using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Drawing;
using System.Windows.Media.Imaging;

using Emgu.CV;
using Emgu.CV.Structure;
using Emgu.Util;
using Emgu.CV.UI;

using System.IO;
using MyGeometry;

using OpenTK.Graphics;
using OpenTK.Graphics.OpenGL;

using Accord.Math;
using Accord.MachineLearning.Geometry;
using Accord.MachineLearning;

namespace SmartCanvas
{
    public class RGBDBuilder
    {
        public static double[] KINECT_LOCAL_PROJ_V1 = new double[]{
            5.1930334103339817e+02, 0.0, 3.2850951551345941e+02, 0.0,
            0.0, 5.1816401430246583e+02, 2.5282555217253503e+02, 0.0,
            0.0, 0.0, 1.0, 0.0
		};

        public static double[] KINECT_LOCAL_PROJ_V2 = new double[] {
			376.6518, 0.0, 265.5583, 0.0,
            0.0, 371.4936, 206.6131, 0.0,
            0.0, 0.0, 1.0, 0.0
		};

        public static double[] KINECT_LOCAL_PROJ_V3 = new double[]{ //  a large dataset of object scans
                    525, 0.0, 319.5, 0.0,
                    0.0, 525, 239.5, 0.0,
                    0.0, 0.0, 1.0, 0.0
		};

        public static double[] KINECT_LOCAL_PROJ = null;
        public static void SetKinectLocalProjModel(int kinect_ver = 1)
        {
            if (kinect_ver == 1)
                KINECT_LOCAL_PROJ = KINECT_LOCAL_PROJ_V1;
            else if (kinect_ver == 2)
                KINECT_LOCAL_PROJ = KINECT_LOCAL_PROJ_V2;
            else if (kinect_ver == 3)
                KINECT_LOCAL_PROJ = KINECT_LOCAL_PROJ_V3;
        }

        private List<RGBDDepthFrame> allframes = null;
        private RGBDDepthFrame currFrame = null;

        public RGBDBuilder() { }
        ~RGBDBuilder() { }

        public RGBDDepthFrame CurrFrame()
        {
            return this.currFrame;
        }
        public Quad3D bbox = null; //bounding box of the scene
        static public double SceneDiagnal = 0;

		 // ////////new
		public int ReadRGBDImageDepthFileV1_new (string file)//return kinect_ver
		{
			int kinect_ver = 3;
			this.allframes = new List<RGBDDepthFrame>();

			string prefix = file.Substring(0, file.LastIndexOf('.'));

            if (prefix[prefix.Length-1]=='d')
                prefix = prefix.Remove(prefix.Length - 1);
             
			string imgdata = prefix + ".png";
			string dptdata = prefix + ".raw";

			RGBDDepthFrame frame = new RGBDDepthFrame();
			frame.Image = new Image<Bgr, byte>(imgdata);
			frame.imagename = prefix;
			float[,] depth = null;
			int w = frame.Image.Width, h = frame.Image.Height;

			if (File.Exists(prefix + "d.png"))
			{
				kinect_ver = 3;
				RGBDBuilder.SetKinectLocalProjModel(3);

				Stream imageStreamSource = new FileStream(prefix + "d.png", FileMode.Open, FileAccess.Read, FileShare.Read);
				PngBitmapDecoder decoder = new PngBitmapDecoder(imageStreamSource, BitmapCreateOptions.PreservePixelFormat, BitmapCacheOption.Default);
				BitmapSource bitmapSource = decoder.Frames[0];

				int height = bitmapSource.PixelHeight;
				int width = bitmapSource.PixelWidth;
				depth = new float[height, width];

				int stride = width * (bitmapSource.Format.BitsPerPixel / 8);
				byte[] bytes = new byte[height * stride];
				bitmapSource.CopyPixels(bytes, stride, 0);

				int counter = 0;
				for (int x = 0; x < height; x++)
				{
					for (int y = 0; y < width; y++)
					{
						byte low = bytes[counter++];
						byte high = bytes[counter++];

						ushort bit16 = (ushort)((high << 8) | low);

						depth[x, y] = bit16 / 65535.0f;
					}
				}
			}
			else if (File.Exists(dptdata))
			{
				kinect_ver = 3;
				RGBDBuilder.SetKinectLocalProjModel(1);
				Stream imageStreamSource = new FileStream(prefix + "d.png", FileMode.Open, FileAccess.Read, FileShare.Read);
				PngBitmapDecoder decoder = new PngBitmapDecoder(imageStreamSource, BitmapCreateOptions.PreservePixelFormat, BitmapCacheOption.Default);
				BitmapSource bitmapSource = decoder.Frames[0];

				int height = bitmapSource.PixelHeight;
				int width = bitmapSource.PixelWidth;
				depth = new float[height, width];

				int stride = width * (bitmapSource.Format.BitsPerPixel / 8);
				byte[] bytes = new byte[height * stride];
				bitmapSource.CopyPixels(bytes, stride, 0);

				int counter = 0;
				for (int x = 0; x < height; x++)
				{
					for (int y = 0; y < width; y++)
					{
						byte low = bytes[counter++];
						byte high = bytes[counter++];

						ushort bit16 = (ushort)((high << 8) | low);

						depth[x, y] = bit16 / 32767.0f;
					}
				}
			}

            if (!File.Exists(prefix + "d.png") && !File.Exists(dptdata))
            {

                this.currFrame = frame;
                this.allframes.Add(frame);

                return kinect_ver;
            }
			// depth & depth color
			Image<Gray, byte> test_depth_image = new Image<Gray, byte>(w, h);	// for test

			List<MyVector3> points = new List<MyVector3>();
			List<double> pointsd = new List<double>();
			List<double> point_colors = new List<double>();
			for (int i = 0; i < h; ++i)
			{
				for (int j = 0; j < w; ++j)
				{
					float dpt = depth[i, j];

					/*--------------------------------------------------------------------
					/ to enable a [i,j] map to 3d points, we disable the following code
					---------------------------------------------------------------------*/
					//if (Math.Abs(dpt) < 1e-6) continue; // no depth; 

					//double x = (j - u) * dpt / foc;
					//double y = (i - v) * dpt / foc;

					//according to kinectmodeling proj.
					double x = j * dpt;
					double y = i * dpt;

					MyMatrix3d KR = new MyMatrix3d();
					KR = MyMatrix3d.IdentityMatrix();
					KR[0, 0] = KINECT_LOCAL_PROJ[0];
					KR[0, 2] = KINECT_LOCAL_PROJ[2];
					KR[1, 1] = KINECT_LOCAL_PROJ[5];
					KR[1, 2] = KINECT_LOCAL_PROJ[6];

					MyVector3 pt = new MyVector3(KR.Inverse() * new MyVector3(x, y, dpt));

					points.Add(pt);//3d
					pointsd.AddRange(pt.ToArray());

					Bgr bgr = frame.Image[i, j];
					point_colors.AddRange(new double[3] { bgr.Red / 255.0, 
						 bgr.Green / 255.0, bgr.Blue / 255.0});

					test_depth_image[i, j] = new Gray(dpt * 255);

				}
			}

			//ImageViewer viewer = new ImageViewer(test_depth_image, "depth");
			//viewer.Show();

			frame.SetPointsV(points.ToArray());
			frame.SetPoints(pointsd.ToArray());
			frame.SetPointsColor(point_colors.ToArray());


            this.bbox = Utils.FindBoundingBox(frame.GetDepthPoints());
            RGBDBuilder.SceneDiagnal = this.bbox.DiagnalLength();

			this.currFrame = frame;
			this.allframes.Add(frame);
			return kinect_ver;
		}

        public void Draw()
        {
			if (this.currFrame != null)
			{ 
                this.currFrame.Draw();
			}

        }

        public MyVector2 ProjectPoint32Image(MyVector3 pt)//project 3d to 2d
        {
            if (pt.z == 0) return new MyVector2();

            double fx = KINECT_LOCAL_PROJ[0];
            double fy = KINECT_LOCAL_PROJ[5];
            double u0 = KINECT_LOCAL_PROJ[2];
            double v0 = KINECT_LOCAL_PROJ[6];
            return new MyVector2(fx * pt.x / pt.z + u0, fy * pt.y / pt.z + v0);
        }

        public void GetAllDepthPointspColor(out List<MyVector3> dpoints, out List<Color> dPointColor)
        {
            int pointcount = currFrame.Image.Width * currFrame.Image.Height;
            dpoints = new List<MyVector3>();
            dPointColor = new List<Color>();
            //int k=0;
            for (int i = 0; i < currFrame.Image.Rows; i++)
                for (int j = 0; j < currFrame.Image.Cols; j++)
                {
                    MyVector3 tempp = currFrame.GetPointV(j, i);
                    if (tempp.x == 0 && tempp.y == 0 && tempp.z == 0)
                        continue;
                    dpoints.Add(tempp);
                    dPointColor.Add(currFrame.GetPointColor(j, i));
                }
        }
    }

    public unsafe class RGBDDepthFrame
    {
        private Image<Bgr, byte> image = null;
        private MyVector3[] depthPointsV3 = null;
        private double[] depthPoints = null;
        private double[] depthPointsColor = null;
        private MyVector3 center;
        public string imagename;

        public RGBDDepthFrame() { }
        ~RGBDDepthFrame() { }

        public MyVector3[] GetDepthPoints()
        {
            return depthPointsV3;
        }

        public MyVector3 Center()
        {
            return this.center;
        }

        public Image<Bgr, byte> Image
        {
            get { return image; }
            set { image = value; }
        }

        public void SegImage(Image<Bgr, byte> img)
        {
            this.image = img;
        }

        public void SetPointsV(MyVector3[] depth_points_v)
        {
            this.depthPointsV3 = depth_points_v;
            this.ComputeCenter();
        }

        public void SetPoints(double[] depth_points)
        {
            this.depthPoints = depth_points;
        }

        public MyVector3 GetPointV(int i, int j)
        {

            return this.depthPointsV3[j * this.image.Width + i];
        }

        public void SetPointsColor(double[] depth_point_colors)
        {
            this.depthPointsColor = depth_point_colors;
        }
        public Color GetPointColor(int i, int j)
        {
            Color c = Color.FromArgb(255,
                (int)depthPointsColor[3 * (j * this.image.Width + i)],
                (int)depthPointsColor[3 * (j * this.image.Width + i) + 1],
                (int)depthPointsColor[3 * (j * this.image.Width + i) + 2]);
            return c;
        }

        private void ComputeCenter()
        {
            if (this.depthPointsV3 == null)
            {
                this.center = new MyVector3();
                return;
            }

            this.center = new MyVector3();
            foreach (MyVector3 v3 in this.depthPointsV3)
            {
                this.center += v3;
            }
            this.center /= this.depthPointsV3.Length;
        }

        public void Draw()
        {
            if (this.depthPointsColor == null || this.depthPoints == null) return;
            GL.Disable(EnableCap.Lighting);
            GL.PointSize(2.0f);
            GL.EnableClientState(ArrayCap.VertexArray);
            GL.EnableClientState(ArrayCap.ColorArray);
            fixed (double* vp = this.depthPoints)
            fixed (double* cp = this.depthPointsColor)
            {
                GL.ColorPointer(3, ColorPointerType.Double, 0, new IntPtr(cp));
                GL.VertexPointer(3, VertexPointerType.Double, 0, new IntPtr(vp));
                GL.DrawArrays(PrimitiveType.Points, 0, this.depthPoints.Length / 3);
            }
            GL.DisableClientState(ArrayCap.ColorArray);
            GL.DisableClientState(ArrayCap.VertexArray);
            GL.Enable(EnableCap.Lighting);
        }
    }

}
