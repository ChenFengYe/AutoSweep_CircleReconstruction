using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using System.Runtime.InteropServices;
using MyGeometry;

using Emgu.CV;
using Emgu.CV.Structure;
using Emgu.Util;
using Emgu.CV.UI;
using Emgu.CV.Util;

using System.Drawing;
using System.IO;
using System.Windows.Forms;

using MRFWrapper;
using PolygonDecompose;

using Accord.Math;
using Accord.Math.Geometry;

namespace SmartCanvas
{
	public struct LMPARAMS
	{
		public double epsg;
		public double epsf;
		public double epsx;
		public double diffstep;
		public int maxits;
		public LMPARAMS (double eg, double ef, double ex, double d, int m)
		{
			epsg = eg;
			epsf = ef;
			epsx = ex;
			diffstep = d;
			maxits = m;
		}
	}

    public unsafe partial class CanvasEngine : IDisposable
    {
        public MyVector2 mouseDownPos = new MyVector2();
        
        // bool variables
        public bool ShowCanvas = false;
		public bool ShowOutLines = false;
		public bool ShowPainting = true;
		public bool ShowOcclusion = false;
		public bool showBackground = true;
        public bool DrawTransparent = true;
        public bool showaxis = false;
        public bool showrgbdpoints = false;

        // private highlight elements
        private bool fixSelection = false;
        private bool isAdjustingRoom = false;
        private MyVector2 OldGlTrans, GlTrans;
        private double glScale = 1.0;
        public float drawingRadius = 8;
		public float SelectionRadius = 16;

        // canvas & sketch
        public Image<Bgr, byte> Canvas;
        public Image<Gray, byte> Sketch;
        public uint CanvasTextureId;
        public uint SketchTextureId;
        public List<CanvasEngine> GuestImageRecords = new List<CanvasEngine>();


		public List<List<List<int>>> candidateGroups = new List<List<List<int>>>();

        // geometry & camera
      
        public Camera camera;	// each image has a camera
        private SketchView view;
		private List<double> currStrokeSpeeds = new List<double>();

        private Trackball arcBall;
        private MyMatrix4d prevTransformation, currTransformation;
		private MyMatrix4d modelTransformation = MyMatrix4d.IdentityMatrix();
        private MyVector3 currTransCenter = new MyVector3();
        private double scaleRatio;

		// RGB-D
		private RGBDBuilder rgbdBuilder = null;

		// file name prefix(for IO)
		public string file_prefix = null;

        // grab cut
		private MyVector2 rect_s = new MyVector2();
		private MyVector2 rect_t = new MyVector2();
		public bool rect_shown = false; 
        private Rectangle rectangle;
        public Image<Gray, Byte> imageMask = null;
		private List<MyModel> ModelList = new List<MyModel>();

		// RANSAC
		private double ransac_thresh = 0.0005;
		private double ransac_probab = 0.99;
		private List<MyPlane> PlaneList = new List<MyPlane>();
		private MyPlane optplane = null;

		// ground feature curve
		double INPLANE_THRESH = 0.0002;
		double INPLANE_THRESH2 = 0.0001;
		private List<MyVector3> groundcurve = null; // for visualization and other processes
		
		// ground circle(need modification)
		private double ransac_thresh_c = 0.001;
		private double ransac_probab_c = 0.8;
		MyVector2 left2;
		MyVector2 right2;

		// ground line
		private double ransac_thresh_l = 0.00005;
		private double ransac_probab_l = 0.8;
		MyVector3 left3l = new MyVector3();
		MyVector3 right3l = new MyVector3();

		// outline
		private List<List<CLineSegment>> OutlineList = new List<List<CLineSegment>>();
		private double ratio = 0; //  = 3d/2d

		// reconstructed models
		private List<MyGCylinder> GCylinderList = new List<MyGCylinder>();
		private List<MyGCuboid> GCuboidList = new List<MyGCuboid>();

		// sweep parameters
		// (both two intervals define the density of 
		// sweeping circles. Choose the better one in a specific case.)
		private double Interval_3d = 2e-4;
		private int Interval_2d = 5; 
		private int THRESH_slicepointnumber = 10;

		// left/right outline start point
		MyVector3 left3 = new MyVector3();
		MyVector3 right3 = new MyVector3();

		 
		// render to test
		List<MyVector3> outlier = new List<MyVector3>();
		List<MyPlane> modelplanes = new List<MyPlane>();
		List<MyVector3> topplanevertices1 = new List<MyVector3>();
		List<MyVector3> topplanevertices2 = new List<MyVector3>();
		List<MyVector3> topplanevertices3 = new List<MyVector3>();

		public void SetRGBDBuilder(RGBDBuilder builder)
		{
			this.rgbdBuilder = builder;
			this.SetTransformCenter();
		}
		public void SetTransformCenter()
		{
			if (this.rgbdBuilder.CurrFrame() != null)
				this.currTransCenter = this.rgbdBuilder.CurrFrame().Center();
		}


		public string sketch2DnameAuto = "";
		public string sketch3DnameAuto = "";
		public string sketch3Dname = "";


		public void CreateKinectDefaultCamera(int kinect_ver)
		{
            if (this.camera != null)
            {
                this.camera.CreateKinectDefault(kinect_ver);
                this.GetGLMatrices();
            }
            //if (this.camera != null)
            //{
            //    this.camera.CreateKinectDefault();
            //    this.GetGLMatrices();
            //}
		}
		public void GetGLMatrices()
		{
			if (this.camera != null)
				this.camera.GetGLMatrices(this.Canvas.Width, this.Canvas.Height, Camera.zNear, Camera.zFar);
		}




		public List<MyVector3> testLines = new List<MyVector3>();
   
        private void CreateCamera()
        {
            int w = this.Canvas.Width, h = this.Canvas.Height;
            this.camera = new Camera();
            this.camera.Init(
                new MyVector2(-200, h / 2 - 200),
                new MyVector2(w + 200, h / 2 - 200),
                w, h
            );
			this.camera.SetViewport(this.Canvas, this.view);
        }
     
   
		// ctrs
		public CanvasEngine(SketchView view, Image<Bgr, byte> canvas, Image<Gray, byte> sketch, uint canvasTxtId, uint sketchTxtId)
        {
            this.view = view;
            this.Canvas = canvas;
            this.Sketch = sketch;
            this.CanvasTextureId = canvasTxtId;
            this.SketchTextureId = sketchTxtId;

            this.CreateCamera();
            this.InitTransform();

            // set stroke radius -- efficiency issue raises if the sampling rate is high
            float scaleFactor = Math.Min(
				(float)this.view.Width / this.Canvas.Width, 
				(float)this.view.Height / this.Canvas.Height
			);
            this.drawingRadius /= scaleFactor;

        }
        public void Dispose()
        {

        }
        
    }

}
