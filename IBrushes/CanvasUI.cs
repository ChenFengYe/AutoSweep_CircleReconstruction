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
        private OpenGLProjector glProjector = null;
        public void Clear()
        {
        }
        private MyMatrix4d meshT = MyMatrix4d.IdentityMatrix();

        // view
        public void ViewingMouseDown(MyVector2 mousedpos, System.Windows.Forms.MouseEventArgs e)
        {
            // mouse down event for viewing
            switch (e.Button)
            {
                case MouseButtons.Left:
                    this.arcBall.Click(mousedpos, Trackball.MotionType.Rotation);
                    break;
                case MouseButtons.Right:
                    this.arcBall.Click(mousedpos / scaleRatio, Trackball.MotionType.Pan);
                    break;
            }
        }
        public void ViewingMouseMove(MyVector2 mousepos, System.Windows.Forms.MouseEventArgs e)
        {
            // mouse move event for viewing
            switch (e.Button)
            {
                case System.Windows.Forms.MouseButtons.Left: this.arcBall.Drag(mousepos); break;
                case System.Windows.Forms.MouseButtons.Right: this.arcBall.Drag(mousepos / scaleRatio); break;
            }
        }
        public void ViewingMouseUp()
        {
            // mouse up event for viewing
            MyMatrix4d m = this.arcBall.GetMatrix();
            this.modelTransformation = m * this.modelTransformation;
            this.arcBall.End();
            this.ObtainGLProjector();
            this.Update2D();
        }
        public void ViewingMouseWheel(MyVector2 mousepos, System.Windows.Forms.MouseEventArgs e)
        {
            this.arcBall.Click(mousepos / scaleRatio, Trackball.MotionType.Scale);
            mousepos += mousepos * 0.6 * e.Delta;
            this.arcBall.Drag(mousepos / scaleRatio);
            MyMatrix4d ms = this.arcBall.GetMatrix();
            this.modelTransformation = ms * this.modelTransformation;
            this.arcBall.End();
            this.ObtainGLProjector();
            this.Update2D();
        }
        private void InitTransform()
        {
            this.arcBall = new Trackball(100, 100);				// trackball for viewing
            this.arcBall.SetBounds(this.Canvas.Width, this.Canvas.Height);
            this.currTransformation =							// tranfsformation matrix
                MyMatrix4d.IdentityMatrix();
            this.prevTransformation = MyMatrix4d.IdentityMatrix();
            this.scaleRatio = (this.Canvas.Width > this.Canvas.Height) ?
                this.Canvas.Height : this.Canvas.Width;
        }
        public void ResetView()
        {
            this.modelTransformation = this.currTransformation = MyMatrix4d.IdentityMatrix();
            this.GlTrans.x = this.GlTrans.y = 0;
            this.glScale = 1.0;
            this.Update2D();
        }



        // opengl utils
        /* -------------------------------------------------------------------------------
        / the following code should be intergrated into the camera class
        / however, due to some awkward bugs, don't have time yet to check -,-
        / -------------------------------------------------------------------------------*/
        public MyVector3 ProjectPointToCanvas3(MyVector2 screenpt, MyVector3 c, MyVector3 nor)
        {
            if (this.glProjector == null) this.ObtainGLProjector();

            screenpt.y = this.view.Height - screenpt.y;
            // the out parameter r measures how close the intersecting point is to the near Canvas3'
            // 1. get transformed Canvas3 normal and center (due to view point change)
            MyVector3 s = this.glProjector.UnProject(screenpt.x, screenpt.y, -1);
            MyVector3 t = this.glProjector.UnProject(screenpt.x, screenpt.y, 1);
            double r = (c - t).Dot(nor) / ((s - t).Dot(nor));
            return r * s + (1 - r) * t;
        }
        private void ProjectPointToCanvas3(MyVector2 screenpt, MyVector3 c, MyVector3 nor, out MyVector3 v, out double r)
        {
            if (this.glProjector == null) this.ObtainGLProjector();

            screenpt.y = this.view.Height - screenpt.y;
            // the out parameter r measures how close the intersecting point is to the near Canvas3'
            // 1. get transformed Canvas3 normal and center (due to view point change)
            MyVector3 s = this.glProjector.UnProject(screenpt.x, screenpt.y, -1);
            MyVector3 t = this.glProjector.UnProject(screenpt.x, screenpt.y, 1);
            r = (c - t).Dot(nor) / ((s - t).Dot(nor));
            v = r * s + (1 - r) * t;
        }
        private MyVector2 Compute2D(MyVector3 point)
        {
            if (this.glProjector == null) this.ObtainGLProjector();
            MyVector2 pt2 = this.glProjector.Project(point).ToMyVector2();
            pt2.y = this.view.Height - pt2.y;
            return pt2;
        }
        public MyVector2 WindowPoint2GLPoint(MyVector2 p)
        {
            return new MyVector2(p.x, this.view.Height - p.y);
        }
        public MyVector2 GLPoint2WindowPoint(MyVector2 p)
        {
            return new MyVector2(p.x, this.view.Height - p.y);
        }
        public MyVector3 UnProject(double winx, double winy, double winz)
        {
            if (this.glProjector == null) this.ObtainGLProjector();
            winy = this.view.Height - winy; // win to opengl
            return this.glProjector.UnProject(winx, winy, winz);
        }

        public void ObtainGLProjector()
        {
            if (!this.camera.Calibrated) return;

            // -------------------------------------------------------------------------
            // draw 3d sense
            int w = this.Canvas.Width, h = this.Canvas.Height;
            GL.Viewport(0, this.view.Height - h, w, h);
            GL.PushMatrix();
            GL.MatrixMode(MatrixMode.Projection);
            GL.LoadIdentity();
            GL.LoadMatrix(this.camera.glprojMatrix);
            GL.MatrixMode(MatrixMode.Modelview);
            GL.LoadIdentity();
            GL.LoadMatrix(this.camera.glmodelViewMatrix);

            MyMatrix4d m = this.arcBall.GetMatrix() * this.modelTransformation;
            m = MyMatrix4d.TranslationMatrix(this.currTransCenter) * m * MyMatrix4d.TranslationMatrix(
                new MyVector3() - this.currTransCenter);
            GL.PushMatrix();
            GL.MultMatrix((m.Transpose()).ToArray());
            this.glProjector = new OpenGLProjector();
            GL.PopMatrix();
            GL.PopMatrix();
        }
        private void Update2D()
        {
            this.ObtainGLProjector(); // update gl projector
        }

    }
}

