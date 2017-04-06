using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Drawing;
using System.Windows.Forms;
using System.IO;

using System.Runtime.InteropServices;
using MyGeometry;

using Emgu.CV;
using Emgu.CV.Structure;
using Emgu.Util;
using Emgu.CV.UI;

using System.Threading;
using System.Threading.Tasks;

using OpenTK;
using OpenTK.Graphics;
using OpenTK.Graphics.OpenGL;

using Microsoft.Kinect;

namespace SmartCanvas
{
    public unsafe class SketchView : OpenTK.GLControl
    {
        public interface IMeshDisplay
        {
            void Display();
            void SetData();
        }

        public enum EnumOperationMode
        {
            Viewing, NONE
        };

        // brush attributes
        public static uint canvasTextureId, drawingTextureId, sketchTextureId,
            brushTextureId, gridTextureId, circleTextureId,
            pencilTextureId, crayonTextureId, inkTextureId, waterColorTextureId, charcoalTextureId,
            loadCanvasTexutreId;
        public static bool EnableDebugMode = false;
        public static Color ShadingColor = Color.DarkGray;
        public static Color PaintColor = Color.Orchid;
        public static Color StrokeColor = Color.Black;
        public static Color HighlightColor = Color.Blue;
        public static Color PencilColor = Color.Gray;
        public static Color ConnectorColor = Color.Gold; //Color.FromArgb(0,128,255);
        public static float BrushSize;
        public static float StrokeSize = 2.0f;
        public static float StrokeDecomposeThreshold = 8.0f;
        public static float Opacity = 0.8f;
        public static double offset = 0.05f;
        public static int nCandidates = 4;
        public static float connectorSize = 8.0f;
        public static int StrokeStyle = 2;

        private uint currCanvasTextureId;
        public MyVector2 prevMousePosition;
        public MyVector2 mouseDownPosition;
        public MyVector2 currMousePosition;
        public bool isMouseDown;

        // sketching
        private const float brushPointSpace = 4;

        private List<Quad3D> candidateQuads;
        private Quad3D groundQuad;

        // member, for image rendering
        private CanvasEngine currCanvasEngine;
        public CanvasEngine CurrCanvasEngine
        {
            get { return currCanvasEngine; }
            set { currCanvasEngine = value; }
        }

        // for 3d rendering
        private EnumOperationMode currentMode;
        public EnumOperationMode CurrentMode { get { return this.currentMode; } set { this.currentMode = value; } }

        private MyVector3 lightPosition;
        private double scaleRatio;

        private MyVector3 eyePos;
        private OpenGLProjector glProjector;


        // constructor
        public SketchView()
        {
            InitializeComponent();

            this.MakeCurrent();

        }


        // public entries
        public void Clear()
        {
            // clear the screen
            if (this.currCanvasEngine != null)
                this.currCanvasEngine.Clear();
            // reset 
            BrushSize = 16;
            this.currCanvasTextureId = 0;

            this.Refresh();

        }
        public void SetMode(EnumOperationMode mode)
        {
            //if (this.currentMode == mode)
            //	this.currentMode = EnumOperationMode.Viewing;
            //else
            this.currentMode = mode;
            this.Focus();
            this.Refresh();
        }
        public bool TurnOnSketchOverlay { get; set; }


        public void Init()
        {
            BrushSize = 16;										// brush size
            this.prevMousePosition = new MyVector2();			// prev mouse pos
            this.mouseDownPosition = new MyVector2();			// mouse down pos
            this.currMousePosition = new MyVector2();			// curr mouse pos

            // initialize textures, airbrushes, and parameters
            this.LoadTextures();								// load canvas, brush textures


            // 3d entities
            //this.eyePos = new MyVector3(0, 0.1, 1);				// eye (camera) position
            this.eyePos = new MyVector3(0, 0, -0.5);				// eye (camera) position
            this.lightPosition = new MyVector3(5, 5, 5);			// light position


            // util
            Utils.InitialRandomColors(53);						// initial random colors

            // opengl
            GLBasicDraw.InitGlMaterialLights();


            // mode
            //this.currentMode = EnumOperationMode.Adjusting;
            this.currentMode = EnumOperationMode.Viewing;

            // opengl projector
            this.glProjector = new OpenGLProjector();

        }
        public void CreateCanvasEngine(Image<Bgr, byte> image)
        {

            Image<Gray, byte> gray = image.Convert<Gray, byte>();
            this.CreateTexture(image, out SketchView.canvasTextureId);

            this.currCanvasTextureId = SketchView.canvasTextureId;

            CanvasEngine engine = new CanvasEngine(this, image, gray, canvasTextureId, sketchTextureId);
            this.currCanvasEngine = engine;

            //ImageViewer iw = new ImageViewer(image, "hello");
            //iw.Show();
        }


        public void Reset()
        {
            if (this.currCanvasEngine != null)
            {
                this.currCanvasEngine.ResetView();
                this.Refresh();
            }
        }

        private void InitializeComponent()			// sketch view components
        {
            this.SuspendLayout();
            // 
            // SketchView
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 12F);
            this.Name = "SketchView";
            this.Load += new System.EventHandler(this.SketchView_Load);
            this.KeyDown += new System.Windows.Forms.KeyEventHandler(this.SketchView_KeyUp);
            this.KeyUp += new System.Windows.Forms.KeyEventHandler(this.SketchView_KeyUp);
            this.ResumeLayout(false);

        }
        private void LoadTextures()					// load textures for canvas and brush
        {
            this.CreateTexture(@"data\canvas.jpg", out SketchView.drawingTextureId);
            //	this.CreateTexture(@"..\..\data\output.png", out brushTextureId);
            this.CreateTexture(@"data\grid.png", out gridTextureId);
            this.CreateTexture(@"data\circle.png", out circleTextureId);
            this.CreateTexture(@"data\pencil.png", out SketchView.pencilTextureId);
            this.CreateTexture(@"data\crayon.png", out crayonTextureId);
            this.CreateTexture(@"data\ink.jpg", out inkTextureId);
            this.CreateTexture(@"data\watercolor.png", out waterColorTextureId);
            this.CreateTexture(@"data\charcoal.jpg", out charcoalTextureId);
        }

        private void CreateTexture(string imagefile, out uint textureid)
        {
            Bitmap image = new Bitmap(imagefile);

            // to gl texture
            Rectangle rect = new Rectangle(0, 0, image.Width, image.Height);
            //	image.RotateFlip(RotateFlipType.RotateNoneFlipY);
            System.Drawing.Imaging.BitmapData bitmapdata = image.LockBits(rect,
                System.Drawing.Imaging.ImageLockMode.ReadOnly, System.Drawing.Imaging.PixelFormat.Format32bppArgb);
            //System.Drawing.Imaging.ImageLockMode.ReadOnly, System.Drawing.Imaging.PixelFormat.Format32bppArgb);

            GL.GenTextures(1, out textureid);
            GL.BindTexture(TextureTarget.Texture2D, textureid);

            GL.TexParameter(TextureTarget.Texture2D, TextureParameterName.TextureMinFilter, (int)TextureMinFilter.Linear);
            GL.TexParameter(TextureTarget.Texture2D, TextureParameterName.TextureMagFilter, (int)TextureMagFilter.Linear);
            GL.TexParameter(TextureTarget.Texture2D, TextureParameterName.TextureWrapS, (int)TextureWrapMode.ClampToEdge);
            GL.TexParameter(TextureTarget.Texture2D, TextureParameterName.TextureWrapT, (int)TextureWrapMode.ClampToEdge);
            GL.TexParameter(TextureTarget.Texture2D, TextureParameterName.TextureWrapR, (int)TextureWrapMode.ClampToEdge);

            GL.TexImage2D(TextureTarget.Texture2D, 0, PixelInternalFormat.Rgba, image.Width, image.Height, 0, PixelFormat.Bgra,
                PixelType.UnsignedByte, bitmapdata.Scan0);
        }
        private void CreateTexture(Image<Bgr, byte> img, out uint textureid) // create gltexture
        {
            Bitmap image = img.ToBitmap();
            // to gl texture
            Rectangle rect = new Rectangle(0, 0, image.Width, image.Height);
            //	image.RotateFlip(RotateFlipType.RotateNoneFlipY);
            System.Drawing.Imaging.BitmapData bitmapdata = image.LockBits(rect,
                System.Drawing.Imaging.ImageLockMode.ReadOnly, System.Drawing.Imaging.PixelFormat.Format32bppArgb);

            GL.GenTextures(1, out textureid);
            GL.BindTexture(TextureTarget.Texture2D, textureid);

            GL.TexParameter(TextureTarget.Texture2D, TextureParameterName.TextureMinFilter, (int)TextureMinFilter.Linear);
            GL.TexParameter(TextureTarget.Texture2D, TextureParameterName.TextureMagFilter, (int)TextureMagFilter.Linear);
            GL.TexParameter(TextureTarget.Texture2D, TextureParameterName.TextureWrapS, (int)TextureWrapMode.ClampToEdge);
            GL.TexParameter(TextureTarget.Texture2D, TextureParameterName.TextureWrapT, (int)TextureWrapMode.ClampToEdge);
            GL.TexParameter(TextureTarget.Texture2D, TextureParameterName.TextureWrapR, (int)TextureWrapMode.ClampToEdge);

            GL.TexImage2D(TextureTarget.Texture2D, 0, PixelInternalFormat.Rgba, image.Width, image.Height, 0, PixelFormat.Bgra,
                PixelType.UnsignedByte, bitmapdata.Scan0);

        }
        private void CreateTexture(Image<Bgra, byte> img, out uint txtid) // create gltexture
        {
            Bitmap image = img.ToBitmap();
            // to gl texture
            Rectangle rect = new Rectangle(0, 0, image.Width, image.Height);
            //	image.RotateFlip(RotateFlipType.RotateNoneFlipY);
            System.Drawing.Imaging.BitmapData bitmapdata = image.LockBits(rect,
                System.Drawing.Imaging.ImageLockMode.ReadOnly, System.Drawing.Imaging.PixelFormat.Format32bppArgb);

            GL.GenTextures(1, out txtid);	// Create The Texture

            // Typical Texture Generation Using Data From The Bitmap
            GL.BindTexture(TextureTarget.Texture2D, txtid);
            GL.TexImage2D(TextureTarget.Texture2D, 0, PixelInternalFormat.Rgba, (int)img.Width, (int)img.Height,
                0, PixelFormat.Bgra, PixelType.UnsignedByte, bitmapdata.Scan0);

            GL.TexParameter(TextureTarget.Texture2D, TextureParameterName.TextureMinFilter, (int)TextureMinFilter.Linear);
            GL.TexParameter(TextureTarget.Texture2D, TextureParameterName.TextureMagFilter, (int)TextureMagFilter.Linear);
            GL.TexParameter(TextureTarget.Texture2D, TextureParameterName.TextureWrapS, (int)TextureWrapMode.ClampToEdge);
            GL.TexParameter(TextureTarget.Texture2D, TextureParameterName.TextureWrapT, (int)TextureWrapMode.ClampToEdge);
            GL.TexParameter(TextureTarget.Texture2D, TextureParameterName.TextureWrapR, (int)TextureWrapMode.ClampToEdge);

        }

        // override
        protected override void OnSizeChanged(EventArgs e)
        {
            base.OnSizeChanged(e);

            Size s = this.Size;
            if (s.Width == 0 || s.Height == 0) return;

            this.scaleRatio = (s.Width > s.Height) ? s.Height : s.Width;
            if (this.currCanvasEngine != null)
                this.currCanvasEngine.ObtainGLProjector();

            this.Refresh();
        }
        protected override void OnMouseDown(System.Windows.Forms.MouseEventArgs e)
        {
            base.OnMouseDown(e);

            this.currMousePosition = new MyVector2(e.X, e.Y);

            this.prevMousePosition = this.mouseDownPosition = currMousePosition;
            this.isMouseDown = true;

            if (this.currCanvasEngine != null)
            {
                this.currCanvasEngine.mouseDownPos = this.mouseDownPosition;
            }
            switch (this.currentMode)
            {
                case EnumOperationMode.Viewing:
                    {
                        if (this.currCanvasEngine != null)
                        {
                            this.currCanvasEngine.ViewingMouseDown(mouseDownPosition, e);
                            this.Refresh();
                        }
                    }
                    break;
            }
        }
        protected override void OnMouseMove(System.Windows.Forms.MouseEventArgs e)
        {
            base.OnMouseMove(e);

            this.currMousePosition = new MyVector2(e.X, e.Y);

            switch (this.currentMode)
            {
                case EnumOperationMode.Viewing:
                    {
                        if (this.currCanvasEngine != null)
                        {
                            if (this.isMouseDown)
                            {
                                this.currCanvasEngine.ViewingMouseMove(currMousePosition, e);
                                this.Refresh();
                            }
                        }
                    }
                    break;
            }

        }
        protected override void OnMouseUp(System.Windows.Forms.MouseEventArgs e)
        {
            base.OnMouseUp(e);

            this.currMousePosition = new MyVector2(e.X, e.Y);
            this.isMouseDown = false;

            switch (this.currentMode)
            {
                case EnumOperationMode.Viewing:
                    {
                        if (this.currCanvasEngine != null)
                        {
                            if (this.currMousePosition == mouseDownPosition) break;

                            this.currCanvasEngine.ViewingMouseUp();

                            this.Refresh();
                        }
                    }
                    break;
            }
        }
        protected override void OnMouseDoubleClick(MouseEventArgs e)
        {
            base.OnMouseDoubleClick(e);
        }
        protected override void OnMouseWheel(System.Windows.Forms.MouseEventArgs e)
        {
            base.OnMouseWheel(e);
            MyVector2 mousepos = new MyVector2(e.Location.X, e.Location.Y);
            switch (this.currentMode)
            {
                case EnumOperationMode.Viewing:
                    {
                        if (this.currCanvasEngine != null)
                        {
                            this.currCanvasEngine.ViewingMouseWheel(mousepos, e);
                            this.Refresh();
                        }
                    }
                    break;
            }
        }

        protected override void OnPaint(PaintEventArgs e)
        {

            if (this.DesignMode == true)
            {
                base.OnPaint(e);
                return;
            }

            this.MakeCurrent();

            GL.ClearColor(1, 1, 1, 1);
            GL.Clear(ClearBufferMask.ColorBufferBit | ClearBufferMask.DepthBufferBit);
            GL.LoadIdentity();

            // draw here
            this.DrawScene();

            this.SwapBuffers();

        }

        // draw
        private void DrawScene()
        {
            if (this.currCanvasEngine != null)
                this.currCanvasEngine.Draw();
        }

        // key events
        private void SketchView_KeyUp(object sender, KeyEventArgs e)
        {
            switch (e.KeyData)
            {

                case Keys.R:
                    {
                        if (this.currCanvasEngine != null)
                        {
                            this.Reset();
                        }
                    }
                    break;
                case Keys.P:
                    {
                        if (this.currCanvasEngine != null)
                            this.currCanvasEngine.Snap("img");
                    }
                        break;
            }
        }

        private void SketchView_Load(object sender, EventArgs e)
        {

        }
    }
}
