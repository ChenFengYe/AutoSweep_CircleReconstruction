using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using MyGeometry;
using System.Drawing;
using System.Windows.Forms;

using OpenTK.Graphics;
using OpenTK.Graphics.OpenGL;

namespace SmartCanvas
{
    public partial class CanvasEngine
    {
        // Rendering

        public void Draw()
        {
            this.GlClearScreen();
            this.DrawCanvas2d(SketchView.Opacity);

            this.DrawScene3d();			// draw 3d objects
            this.DrawScene2d();			// overlaid onto 3d
        }

        private void GlClearScreen()
        {
            GL.Disable(EnableCap.Dither);
            GL.Disable(EnableCap.Blend);
            GL.Disable(EnableCap.DepthTest);
            GL.Disable(EnableCap.Normalize);
            GL.Disable(EnableCap.Lighting);
            GL.Disable(EnableCap.ColorMaterial);
            GL.ClearColor(1, 1, 1, 0);
            GL.ClearDepth(1);
            GL.Clear(ClearBufferMask.ColorBufferBit | ClearBufferMask.DepthBufferBit);
        }

        private void DrawCanvas2d(float opacity)
        {
            if (this.Canvas == null || !this.showBackground)
            {
                return;
            }
            GL.Disable(EnableCap.DepthTest);
            GL.Enable(EnableCap.Blend);
            GL.BlendFunc(BlendingFactorSrc.SrcAlpha, BlendingFactorDest.OneMinusSrcAlpha);
            int w = this.Canvas.Width, h = this.Canvas.Height;
            GL.PushMatrix();
            GL.Viewport(0, this.view.Height - h, w, h);
            GL.MatrixMode(MatrixMode.Projection);
            GL.LoadIdentity();
            GL.Ortho(0, w, h, 0, -1, 1);
            GL.MatrixMode(MatrixMode.Modelview);
            GL.PushMatrix();

            // translation and scaling
            GL.Translate(this.GlTrans.x, this.GlTrans.y, 0.0);
            GL.Scale(this.glScale, this.glScale, 1);

            // draw canvas
            GLBasicDraw.DrawTexture(w, h, this.CanvasTextureId, opacity);

            GL.PopMatrix();
            GL.PopMatrix();
            GL.Disable(EnableCap.Blend);
        }
        private void DrawScene2d()
        {
            if (this.Canvas == null) return;

            int w = this.Canvas.Width, h = this.Canvas.Height;

            GL.PushMatrix();
            GL.Viewport(0, this.view.Height - h, w, h);
            GL.MatrixMode(MatrixMode.Projection);
            GL.LoadIdentity();
            GL.Ortho(0, w, h, 0, -1, 1);
            GL.MatrixMode(MatrixMode.Modelview);
            GL.PushMatrix();
            // translation and scaling
            GL.Translate(this.GlTrans.x, this.GlTrans.y, 0.0);
            GL.Scale(this.glScale, this.glScale, 1);

            GL.Enable(EnableCap.PointSmooth);
            GL.Enable(EnableCap.LineSmooth);
            GL.Enable(EnableCap.Blend);



            // draw rectangle
            if (this.rect_shown)
            {
                this.DrawRectangle(); //grab cut rectangle
            }

            //// draw ground circle in canvas
            //if (this.circleInCanvas != null)
            //	this.DrawCircleInCanvas();

            this.DrawEndpoints();

            // draw highlights
            this.DrawHighlights2D();

            this.DrawLight();

            GL.PopMatrix();
            GL.PopMatrix();
        }
        private void DrawScene3d()
        {
            if (this.camera.Calibrated == false) return;

            GLBasicDraw.InitGlMaterialLights();	// lights and materials
            int w = this.Canvas.Width, h = this.Canvas.Height;
            GL.Viewport(0, this.view.Height - h, w, h);

            // -------------------------------------------------------------------------
            // draw 3d sense
            MyMatrix4d m = this.arcBall.GetMatrix() * this.modelTransformation;
            m = MyMatrix4d.TranslationMatrix(this.currTransCenter) * m * MyMatrix4d.TranslationMatrix(
                new MyVector3() - this.currTransCenter);

            GL.PushMatrix();
            GL.MatrixMode(MatrixMode.Projection);
            GL.LoadIdentity();
            GL.LoadMatrix(this.camera.glprojMatrix);
            GL.MatrixMode(MatrixMode.Modelview);
            GL.LoadIdentity();
            GL.LoadMatrix(this.camera.glmodelViewMatrix);

            GL.PushMatrix();
            GL.MultMatrix(m.Transpose().ToArray());

            // Enable antialiased lines
            GL.Enable(EnableCap.PointSmooth);
            GL.Enable(EnableCap.LineSmooth);
            GL.Enable(EnableCap.Blend);
            GL.BlendFunc(BlendingFactorSrc.SrcAlpha, BlendingFactorDest.OneMinusSrcAlpha);




            // draw all geometric faces (depth test enabled, blend disabled)
            if (!this.DrawTransparent)
            {
                GL.Enable(EnableCap.DepthTest);
                GL.Disable(EnableCap.Blend);
                GL.DepthFunc(DepthFunction.Less);
                GL.DepthMask(true);
            }

            if (!this.DrawTransparent)
            {
                GL.DepthFunc(DepthFunction.Lequal);
                GL.DepthMask(false); // Do not remove me, else get dotted lines
            }

            GL.Enable(EnableCap.Blend);
            /* --------------------------------------------------------------------------*/
            // draw the lines on top
            // painting goes here!
            /* --------------------------------------------------------------------------*/
            // draw RGB-D
            GL.Disable(EnableCap.Blend);
            GL.Enable(EnableCap.DepthTest);
            GL.Disable(EnableCap.Blend);
            GL.DepthFunc(DepthFunction.Less);
            GL.DepthMask(true);
            if (this.rgbdBuilder != null && this.showrgbdpoints)
                this.rgbdBuilder.Draw();
            GL.Disable(EnableCap.DepthTest);
            GL.Enable(EnableCap.Blend);


            // draw plane
            if (this.PlaneList != null)
            {
                for (int i = 0; i < this.PlaneList.Count; i++)
                {
                    this.PlaneList[i].DrawMyPlane();
                }
            }

            // draw ground feature curve
            if (this.groundcurve != null)
                this.DrawGroundFeatureCurve();

            // draw ground feature line
            this.DrawGroundFeatureLine();

            this.DrawTopPlane();

            // draw model plane
            this.DrawModelPlane();

            // draw two points
            this.DrawHighLow();

            // draw gcylinder
            if (this.GCylinderList != null)
            {
                foreach (MyGCylinder mgc in this.GCylinderList)
                {
                    mgc.Draw();
                }
            }

            // draw gcuboid
            if (this.GCuboidList != null)
            {
                foreach (MyGCuboid mgc in this.GCuboidList)
                {
                    mgc.Draw();
                }
            }

            this.DrawOutliers();


            DrawAxis();

            // draw highlights
            this.DrawHighlights3D();

            GL.Disable(EnableCap.LineSmooth);
            GL.Disable(EnableCap.PointSmooth);
            GL.Disable(EnableCap.Blend);
            GL.DepthMask(true);
            GL.PopMatrix();
            GL.PopMatrix();
        }

        private void DrawAxis()
        {
            if (this.showaxis)
            {
                GL.Disable(EnableCap.Lighting);
                GL.Disable(EnableCap.Normalize);
                GL.LineWidth(2.0f);
                GL.Begin(PrimitiveType.Lines);

                GL.Color3(Color.Red);
                GL.Vertex3(0, 0, 0);
                GL.Vertex3(1, 0, 0);

                GL.Color3(Color.Green);
                GL.Vertex3(0, 0, 0);
                GL.Vertex3(0, 1, 0);

                GL.Color3(Color.Blue);
                GL.Vertex3(0, 0, 0);
                GL.Vertex3(0, 0, 1);

                GL.End();
                GL.LineWidth(1.0f);
                GL.Enable(EnableCap.Lighting);
            }

        }

        private void DrawLight()
        {
            for (int i = 0; i < GLBasicDraw.lightPositions.Count; ++i)
            {
                MyVector3 pos3 = new MyVector3(GLBasicDraw.lightPositions[i][0],
                    GLBasicDraw.lightPositions[i][1],
                    GLBasicDraw.lightPositions[i][2]);
                MyVector3 pos2 = this.camera.Project(pos3.x, pos3.y, pos3.z);
                GLBasicDraw.DrawCircle(new MyVector2(pos2.x, pos2.y), Color.Yellow, 0.2f);
            }
        }

        private void DrawRectangle()
        {
            GL.Disable(EnableCap.Lighting);
            GL.Color3(Color.Black.R, Color.Black.G, Color.Black.B);
            GL.LineWidth(3.0f);
            GL.Begin(PrimitiveType.LineLoop);
            GL.Vertex2(rect_s.ToArray());
            GL.Vertex2(rect_t.x, rect_s.y);
            GL.Vertex2(rect_t.ToArray());
            GL.Vertex2(rect_s.x, rect_t.y);
            GL.End();
            GL.Enable(EnableCap.Lighting);
        }
        //private void DrawCircleInCanvas ()
        //{
        //	GL.Disable(EnableCap.Lighting);
        //	GL.Color3(Color.Cyan.R, Color.Cyan.G, Color.Cyan.B);
        //	GL.LineWidth(3.0f);
        //	GL.Begin(PrimitiveType.LineLoop);
        //	foreach (MyVector2 pixel in circleInCanvas)
        //	{
        //		GL.Vertex2(pixel.x, pixel.y);
        //	}
        //	GL.End();
        //	GL.Enable(EnableCap.Lighting);
        //}
        private void DrawEndpoints()
        {
            GL.Disable(EnableCap.Lighting);
            GL.Color3(Color.Cyan.R, Color.Cyan.G, Color.Cyan.B);
            GL.LineWidth(3.0f);
            GL.Begin(PrimitiveType.Points);
            GL.Vertex2(left2.ToArray());
            GL.Vertex2(right2.ToArray());
            GL.End();
            GL.Enable(EnableCap.Lighting);
        }
        private void DrawGroundFeatureCurve()
        {
            GL.Disable(EnableCap.Lighting);
            GL.PointSize(2.0f);
            GL.Color3(Color.Green.R, Color.Green.G, Color.Green.B);
            GL.Begin(PrimitiveType.Points);
            foreach (MyVector3 vert in groundcurve)
            {
                GL.Vertex3(vert.ToArray());
            }
            GL.End();
            GL.Enable(EnableCap.Lighting);
        }
        private void DrawGroundFeatureLine()
        {
            GL.Disable(EnableCap.Lighting);
            GL.PointSize(2.0f);
            GL.Color3(Color.Blue.R, Color.Blue.G, Color.Blue.B);
            GL.Begin(PrimitiveType.Lines);

            GL.Vertex3(left3l.ToArray());
            GL.Vertex3(right3l.ToArray());

            GL.End();
            GL.Enable(EnableCap.Lighting);
        }
        private void DrawTopPlane()
        {
            GL.Disable(EnableCap.Lighting);
            GL.PointSize(2.0f);
            GL.Color3(Color.Red.R, Color.Red.G, Color.Red.B);
            GL.Begin(PrimitiveType.Points);
            foreach (MyVector3 vert in this.topplanevertices1)
            {
                GL.Vertex3(vert.ToArray());
            }
            GL.End();
            GL.Color3(Color.Blue.R, Color.Blue.G, Color.Blue.B);
            GL.Begin(PrimitiveType.Points);
            foreach (MyVector3 vert in this.topplanevertices2)
            {
                GL.Vertex3(vert.ToArray());
            }
            GL.End();
            GL.Color3(Color.Green.R, Color.Green.G, Color.Green.B);
            GL.Begin(PrimitiveType.Points);
            foreach (MyVector3 vert in this.topplanevertices3)
            {
                GL.Vertex3(vert.ToArray());
            }
            GL.End();
            GL.Enable(EnableCap.Lighting);
        }
        private void DrawModelPlane()
        {
            for (int i = 0; i < this.modelplanes.Count; i++)
            {
                this.modelplanes[i].DrawMyPlane();
            }
        }
        private void DrawOutliers()
        {
            GL.Disable(EnableCap.Lighting);
            GL.PointSize(2.0f);
            GL.Color3(Color.Orange.R, Color.Orange.G, Color.Orange.B);
            GL.Begin(PrimitiveType.Points);
            foreach (MyVector3 vert in outlier)
            {
                GL.Vertex3(vert.ToArray());
            }
            GL.End();
            GL.Enable(EnableCap.Lighting);
        }
        private void DrawHighLow()
        {
            GL.Disable(EnableCap.Lighting);
            GL.PointSize(6.0f);
            GL.Color3(Color.Red.R, Color.Red.G, Color.Red.B);
            GL.Begin(PrimitiveType.Points);
            GL.Vertex3(left3.ToArray());
            GL.Vertex3(right3.ToArray());
            GL.End();
            GL.Enable(EnableCap.Lighting);
        }


        private void DrawHighlights2D()
        {
            //mark.Draw();
            if (drawProjectedPoint)
            {
                GL.Disable(EnableCap.Lighting);
                GL.Color3(Color.Cyan.R, Color.Cyan.G, Color.Cyan.B);
                GL.LineWidth(3.0f);
                GL.Begin(PrimitiveType.Points);
                for (int i = 0; i < CirclePoints_2d.Count; i++)
                {
                    GL.Vertex2(CirclePoints_2d[i].ToArray());
                    GL.Vertex2(CirclePoints_2d[i].ToArray());
                }
                GL.End();
                GL.Enable(EnableCap.Lighting);
            }
        }

        private void DrawHighlights3D()
        {

            topCircle.Draw();
        }

    }// GLViewer
}
