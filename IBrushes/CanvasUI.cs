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
    

        //choose object by drawing rectangle 
        //chang
        public void ChoosingMouseDown(MyVector2 mousepos)
        {
            rect_s = mousepos;
        }
        public void ChoosingMouseMove(MyVector2 mousepos)
        {
            rect_t = mousepos;
        }
        public void ChoosingMouseUp(MyVector2 mousepos)
        {
            rect_t = mousepos;
            rectangle = new Rectangle((int)rect_s.x, (int)rect_s.y, (int)(rect_t.x - rect_s.x), (int)(rect_t.y - rect_s.y));
            this.Segmentation();
            this.SaveModelcut();
            this.rect_shown = false;
            rectangle = new Rectangle();

            // find outline
            this.FindOutline();

            // reconstruct either g-cuboid or g-cylinder
            this.Reconstruct();


            this.view.CurrentMode = SketchView.EnumOperationMode.Viewing;
            this.view.Focus();
            this.view.Refresh();
        }
        private void Segmentation()
        {
            imageMask = Canvas.GrabCut(rectangle, 5);
            ScalarArray one = new ScalarArray(1d);
            CvInvoke.BitwiseAnd(imageMask, one, imageMask);
            Image<Bgr, Byte> result = Canvas.Copy(imageMask);
            ImageViewer iv = new ImageViewer(result, "segmentation");
            iv.Show();
        }
        public void SaveModelcut()
        {
            List<MyVector3> modelcut = new List<MyVector3>();
            for (int i = 0; i < imageMask.Height; i++)
            {
                for (int j = 0; j < imageMask.Width; j++)
                {
                    double k = imageMask[i, j].Intensity;
                    MyVector3 vertex = this.rgbdBuilder.CurrFrame().GetPointV(j, i);
                    if (k == 1d && !vertex.IsNull()) // delete some noise
                        modelcut.Add(vertex);
                }
            }
            ModelList.Add(new MyModel(modelcut));
        }
        private void FindOutline()
        {
            VectorOfVectorOfPoint contourList = new VectorOfVectorOfPoint();
            Image<Gray, Byte> hierarchy = new Image<Gray, Byte>(imageMask.Size);
            CvInvoke.FindContours(imageMask, contourList, hierarchy, RetrType.External, ChainApproxMethod.ChainApproxSimple);
            #region Visualization
            Image<Bgr, Byte> imageContour = new Image<Bgr, Byte>(imageMask.Size);
            CvInvoke.DrawContours(imageContour, contourList, 0, new MCvScalar(255, 0, 255), 2);
            ImageViewer iv = new ImageViewer(imageContour, "Contour");
            iv.Show();
            #endregion
            // transfer to my data structure
            List<CLineSegment> outline = new List<CLineSegment>();
            int N = contourList[0].Size;
            for (int i = 0; i < N; i++)
            {
                MyVector2 s = new MyVector2(contourList[0][i].X, contourList[0][i].Y);
                MyVector2 t = new MyVector2(contourList[0][(i + 1) % N].X, contourList[0][(i + 1) % N].Y);
                if (s.y > t.y)
                    outline.Add(new CLineSegment(s, t));
                else
                    outline.Add(new CLineSegment(t, s));
            }

            // sort from bottom to up
            MyCompareSegmentLowY mcs = new MyCompareSegmentLowY();
            outline.Sort(mcs);
            outline.Reverse();

            this.OutlineList.Add(outline);
            //this.ShowOutline();
        }
        private void ShowOutline()
        {
            Image<Bgr, Byte> imageContour = new Image<Bgr, Byte>(this.Canvas.Size);
            for (int i = 0; i < this.OutlineList.Count; i++)
            {
                System.Drawing.Point[] pts = new System.Drawing.Point[this.OutlineList[i].Count];
                for (int j = 0; j < this.OutlineList[i].Count; j++)
                {
                    CLineSegment cls = this.OutlineList[i][j];
                    pts[j] = new System.Drawing.Point((int)cls.StartPoint.x, (int)cls.StartPoint.y);
                }
                imageContour.DrawPolyline(pts, true, new Bgr(Color.Magenta), 2);
            }
            ImageViewer iv = new ImageViewer(imageContour, "Outlines");
            iv.Show();
        }
        public void Reconstruct()
        {
            #region Find ground feature curve
            Dictionary<int, int> CurveModelMapping = new Dictionary<int, int>();
            // find ground curve
            if (optplane == null)
                this.FindOptimalPlanes2();
            MyVector3 opti_norm = new MyVector3(optplane.Normal());
            float offset = optplane.Offset();
            groundcurve = new List<MyVector3>();
            MyModel m = this.ModelList.First();
            for (int i = 0; i < m.Vertices.Count; i++)
            {
                MyVector3 vert = m.Vertices[i];
                double dist = Math.Abs(opti_norm.Dot(vert) + offset);
                if (dist <= INPLANE_THRESH)
                {
                    CurveModelMapping.Add(groundcurve.Count, i);
                    groundcurve.Add(vert);
                }
            }
            // project ground curve to plane found
            MyVector3 p_in_plane = new MyVector3(0.0, 0.0, -offset / opti_norm.z);
            for (int i = 0; i < groundcurve.Count; i++)
            {
                double bias = (groundcurve[i] - p_in_plane).Dot(opti_norm);
                groundcurve[i] = groundcurve[i] - opti_norm * bias;
            }
            #endregion

            #region Project to 2d-plane
            // prepare 2d coordinate
            MyVector3 plane_norm = new MyVector3(optplane.Normal());
            MyVector3 X = new MyVector3(1, 0, 0);
            MyVector3 Y = new MyVector3(0, 1, 0);
            MyVector3 Z = new MyVector3(0, 0, 1);
            MyVector3 rotAxis = Z.Cross(plane_norm).Normalize();
            double angle_cos = Z.Dot(plane_norm);
            if (angle_cos > 1) angle_cos = 1;
            if (angle_cos < -1) angle_cos = -1;
            double rotAngle = Math.Acos(angle_cos);

            MyMatrix4d Mat = MyMatrix4d.RotationMatrix(rotAxis, rotAngle);
            MyVector3 X_new = (Mat * new MyVector4(X)).XYZ().Normalize();
            MyVector3 Y_new = (Mat * new MyVector4(Y)).XYZ().Normalize();
            CoordinateFrame frame = new CoordinateFrame(optplane.planeCenter, X_new, Y_new, plane_norm);

            // prepare 2d data : project 3d points to plane with uv coordinates
            Accord.Point[] points = new Accord.Point[groundcurve.Count];
            int counter = 0;
            foreach (MyVector3 vert in groundcurve)
            {
                MyVector3 vert2d = frame.GetPointLocalCoord(vert);
                points[counter].X = (float)vert2d.x;
                points[counter].Y = (float)vert2d.y;
                counter++;
            }
            #endregion

            //#region Fit as a circle
            //RansacCircle rc = new RansacCircle(ransac_thresh_c, ransac_probab_c);
            //Circle c = rc.Estimate(points);
            //int[] c_inliers = rc.Inliers;
            //#endregion

            #region Fit as line
            RansacLine rl = new RansacLine(ransac_thresh_l, ransac_probab_l);
            Line l = rl.Estimate(points);
            int[] l_inliers = rl.Inliers;
            #endregion

            //#region Assessment
            //if (CircleDominate(ref points, ref c_inliers, ref l_inliers))
            //{
            //	System.Console.WriteLine("circle");
            //	this.ReconstructGCylinder(c, ref frame);
            //}
            //else
            //{
            System.Console.WriteLine("line");
            //this.DetectVisiblePlanes(); // may not be plane => give up
            this.PreSweeping(ref frame, ref points, ref l_inliers);
            //this.InitGCuboid();
            //this.RefineGCuboid();
            //this.ReconstructGCuboid2(l, ref frame, ref points, ref l_inliers);
            //}
            //#endregion
        }
        public void PreSweeping(ref CoordinateFrame frame, ref Accord.Point[] ps, ref int[] inliers)
        {
            #region Compute params : ratio(=3d/2d), sweep_height, sweep_interval
            // ratio
            List<MyVector3> inlierlist3d = new List<MyVector3>();
            List<MyVector2> inlierlist2d = new List<MyVector2>();
            Dictionary<MyVector2, MyVector3> map = new Dictionary<MyVector2, MyVector3>();
            foreach (int idx in inliers)
            {
                MyVector2 p = new MyVector2(ps[idx].X, ps[idx].Y);
                MyVector3 v = frame.GetPointSpaceCoord(new MyVector3(p, 0.0));
                inlierlist3d.Add(v);
                MyVector2 pix = this.rgbdBuilder.ProjectPoint32Image(v);
                inlierlist2d.Add(pix);
                map.Add(pix, v);
            }
            inlierlist2d.Sort((a, b) => a.x.CompareTo(b.x));
            MyVector3 v0; MyVector3 v1;
            bool flag1 = map.TryGetValue(inlierlist2d.Last(), out v0);
            bool flag2 = map.TryGetValue(inlierlist2d.First(), out v1);
            if (flag1 && flag2)
            {
                double dist3 = (v1 - v0).Length();
                double dist2 = (inlierlist2d.First() - inlierlist2d.Last()).Length();
                this.ratio = dist3 / dist2;
            }

            // sweep_height, sweep_interval
            int sweep_height;
            double sweep_interval_3d;
            this.FindSweepParams(out sweep_height, out sweep_interval_3d);
            #endregion

            #region Pre-Sweeping
            MyGCuboid mgc = new MyGCuboid();
            mgc.AddPolygon(new MyPolygon(optplane));
            int sweep_level = 0;
            while (sweep_level < sweep_height)
            {
                this.FindNextPlane(ref mgc, ref sweep_interval_3d);
                sweep_level++;
            }
            #endregion

            #region Slicing
            mgc.SliceModel(this.ModelList[0]);
            #endregion

            topplanevertices1.AddRange(mgc.polyList.Last().SlicePoints3d);
            topplanevertices2.AddRange(mgc.polyList[mgc.polyList.Count - 2].SlicePoints3d);
            topplanevertices3.AddRange(mgc.polyList[mgc.polyList.Count - 3].SlicePoints3d);

            this.GCuboidList.Add(mgc);
        }

        //public void ReconstructGCuboid2 (Line l, ref CoordinateFrame frame, ref Accord.Point[] ps, ref int[] inliers)
        //{
        //	MyGCuboid mgc = new MyGCuboid();
        //	#region Add ground rectangle
        //	List<MyVector3> inlierlist3d = new List<MyVector3>();
        //	List<MyVector2> inlierlist2d = new List<MyVector2>();
        //	Dictionary<MyVector2, MyVector3> map = new Dictionary<MyVector2, MyVector3>();
        //	foreach (int idx in inliers)
        //	{
        //		MyVector2 p = new MyVector2(ps[idx].X, ps[idx].Y);
        //		MyVector3 v = frame.GetPointSpaceCoord(new MyVector3(p, 0.0));
        //		inlierlist3d.Add(v); 
        //		MyVector2 pix = this.rgbdBuilder.ProjectPoint32Image(v);
        //		inlierlist2d.Add(pix);
        //		map.Add(pix, v);
        //	}
        //	inlierlist2d.Sort((a, b) => a.x.CompareTo(b.x));
        //	MyVector3 v0; MyVector3 v1;
        //	bool flag1 = map.TryGetValue(inlierlist2d.Last(), out v0);
        //	bool flag2 = map.TryGetValue(inlierlist2d.First(), out v1);
        //	if (flag1 && flag2)
        //	{
        //		mgc.AddPolygon(new MyPolygon(v0, v1, optplane));
        //	}
        //	#endregion

        //	#region Prepare sweeping
        //	int sweep_height;
        //	double sweep_interval_3d;
        //	this.FindSweepParams(ref mgc, out sweep_height, out sweep_interval_3d);
        //	#endregion

        //	int sweep_level = 0;
        //	while (sweep_level < sweep_height)
        //	{
        //		this.FindNextPolygon(ref mgc, ref sweep_interval_3d);
        //		sweep_level++;
        //	}

        //	// global optimization
        //	//this.OptimizeGCuboid(ref mgc);

        //	this.GCuboidList.Add(mgc);
        //}
        private void FindSweepParams(out int sweep_height, out double sweep_interval_3d)
        {
            double height = MyModel.ComputeHeight(this.ModelList.First(), this.optplane);
            int from_3d = (int)Math.Ceiling(height / this.Interval_3d);
            int from_2d = (int)Math.Ceiling(height / (this.Interval_2d * this.ratio));
            // Choose denser one
            if (from_3d >= from_2d)
            {
                sweep_height = from_3d;
                sweep_interval_3d = this.Interval_3d;
            }
            else
            {
                sweep_height = from_2d;
                sweep_interval_3d = this.Interval_2d * this.ratio;
            }
        }
        private void FindNextPlane(ref MyGCuboid mgc, ref double interval_3d)
        {
            MyPolygon poly_last = mgc.PolygonList.Last();

            // find next plane(simply translate last plane)
            MyPlane p_next;
            MyPlane.Translate(poly_last.BelongPlane, out p_next, (float)interval_3d);

            MyPolygon poly_next = new MyPolygon(p_next);
            poly_next.AddNeighbor(poly_last);
            poly_last.AddNeighbor(poly_next);

            // update
            mgc.AddPolygon(poly_next);
        }
        public void FindOptimalPlanes2()
        {
            // prepare RANSAC data
            List<MyVector3> dPoints = new List<MyVector3>();
            List<Color> dPointsColor = new List<Color>();
            this.rgbdBuilder.GetAllDepthPointspColor(out dPoints, out dPointsColor);
            List<int> inlier_i = new List<int>();

            // do RANSAC and complete information of the plane
            MyVector4 p = pointToPlaneClustering.RANSACNormal(dPoints, inlier_i);
            PlaneList.Clear();
            List<Accord.Math.Point3> inlier_p = new List<Accord.Math.Point3>();
            MyVector3 center = new MyVector3();
            foreach (int i in inlier_i)
            {
                center.x += dPoints[i].x;
                center.y += dPoints[i].y;
                center.z += dPoints[i].z;
                inlier_p.Add(new Accord.Math.Point3((float)dPoints[i].x, (float)dPoints[i].y, (float)dPoints[i].z));
            }
            center = center / inlier_i.Count;
            PlaneList.Add(new MyPlane(new Plane(-(float)p.x, -(float)p.y, -(float)p.z, (float)p.w), center, inlier_p, PlaneList.Count - 1, PlaneList.Count - 1));
            PlaneList.Last().ComputeVertices();
            PlaneList.Last().ComputeBoundQuad(); // necessary!
            PlaneList.Last().SaveMyPlane(file_prefix);

            // find optimal plane( design algorithm here )
            optplane = PlaneList.First();
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

