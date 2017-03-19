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
    public unsafe partial class CanvasEngine : IDisposable
    {


        public MyVector2 mouseDownPos = new MyVector2();

        // bool variables
        public bool showBackground = true;
        public bool showaxis = true;
        public bool DrawTransparent = true;
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

        // file name prefix(for IO)
        public string file_prefix = null;

        public void SetTransformCenter()
        {
            this.currTransCenter = new MyVector3(0, 0, 1);
        }

        public string sketch2DnameAuto = "";
        public string sketch3DnameAuto = "";
        public string sketch3Dname = "";


        public void CreateDefaultCamera(int kinect_ver)
        {
            if (this.camera != null)
            {
                this.camera.CreateKinectDefault();
                this.GetGLMatrices();
            }
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
