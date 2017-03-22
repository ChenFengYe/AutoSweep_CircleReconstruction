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

            // draw highlights
            this.DrawHighlights2D();

           // this.DrawLight();

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
                GL.Vertex3(0, 0, 1);
                GL.Vertex3(0.5, 0, 1);

                GL.Color3(Color.Green);
                GL.Vertex3(0, 0, 1);
                GL.Vertex3(0, 0.5, 1);

                GL.Color3(Color.Blue);
                GL.Vertex3(0, 0, 1);
                GL.Vertex3(0, 0, 1.5);

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


        private void DrawHighlights2D()
        {
            GL.Disable(EnableCap.Lighting);
            //mark.Draw();
            if (drawProjectedPoint)
            {
                GL.Color3(Color.Cyan.R, Color.Cyan.G, Color.Cyan.B);
                GL.LineWidth(3.0f);
                GL.Begin(PrimitiveType.Points);
                for (int i = 0; i < CirclePoints_2d.Count; i++)
                {
                    GL.Vertex2(CirclePoints_2d[i].ToArray());
                    GL.Vertex2(CirclePoints_2d[i].ToArray());
                }
                GL.End();
            }

            //GL.Color3(Color.Salmon);
            //GL.Begin(PrimitiveType.Points);
            //int c = 0;

            //for (int i = 0; i < this.trajpoints.Count; i++)
            //{
            //    GL.Color3(c++, c++, c++);
            //    GL.Vertex2(trajpoints[i].ToArray());
            //}
            //GL.End();


            GL.Enable(EnableCap.Lighting);


        }

        private void DrawHighlights3D()
        {
            GL.Disable(EnableCap.Lighting);
            //topCircle.DrawCapped(Color.Salmon);


            if (this.body != null)
                body.Draw(Color.CornflowerBlue);



            GL.PointSize(1.0f);
            GL.Begin(PrimitiveType.Points);

            int c = 0;
            foreach (MyVector3 p in this.trajpoints_3d)
            {
                c = c + 10;
                if (c > 255)
                    c = 0;
                GL.Color3(Color.FromArgb(c, 0, 0));
                GL.Vertex3(p.ToArray());
            }
            GL.End();
            GL.PointSize(1.0f);



            GL.Enable(EnableCap.Lighting);
        }

    }// GLViewer
}
