using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using System.IO;

using System.Drawing;
using MyGeometry;

using Emgu.CV;
using Emgu.CV.CvEnum;
using Emgu.CV.Structure;

using Loyc.Geometry;

using OpenTK.Graphics;
using OpenTK.Graphics.OpenGL;

namespace SmartCanvas
{
    public class Utils
    {
        // create a blur texture with gaussian alpha channel
        // the color components of the texture are less useful for airbrush
        // the texture size is [2*size, 2*size]
        public static Color4[] CreateGaussianTexture(int size)
        {
            int w = size * 2, size2 = size * size;
            Color4[] pixels = new Color4[w * w];
            for (int i = 0; i < w; ++i)
            {
                int dx = i - size;
                for (int j = 0; j < w; ++j)
                {
                    int J = i * w + j;

                    int dy = j - size;
                    double dist2 = (dx * dx + dy * dy);

                    byte gray = 0;
                    if (dist2 <= size2)	// -- not necessary actually, similar effects
                    {
                        // set gaussian values for the alphas
                        // modify the denominator to get different over-paiting effects
                        double gau_val = Math.Exp(-dist2 / (2 * size2 / 8));
                        gray = Math.Min((byte)255, (byte)(gau_val * 255));
                    }

                    pixels[J] = new Color4(gray, gray, gray, gray);

                }
            }

            return pixels;
        }

        public static double thresh = 1e-6;
        public static double pixel_thresh = 20.0;
        public static bool isNull(MyVector2 v)
        {
            if (Double.IsNaN(v.x) || Double.IsNaN(v.y))
            {
                return true;
            }
            return false;
        }

        public static bool isNull(MyVector3 v)
        {
            if (Double.IsNaN(v.x) || Double.IsNaN(v.y) || Double.IsNaN(v.z))
            {
                return true;
            }
            return false;
        }

        // find if a point in inside a polygon
        public static bool PointInPoly(MyVector2 p, MyVector2[] points)
        {
            bool c = false;
            int n = points.Length;
            for (int i = 0, j = n - 1; i < n; j = i++)
            {
                if (((points[i].y > p.y) != (points[j].y > p.y)) &&
                    (p.x < (points[j].x - points[i].x) * (p.y - points[i].y) / (points[j].y - points[i].y) + points[i].x))
                    c = !c;
            }
            return c;
        }
        public static bool PointInPoly(MyVector2 p, Point2[] points)
        {
            bool c = false;
            int n = points.Length;
            for (int i = 0, j = n - 1; i < n; j = i++)
            {
                if (((points[i].pos.y > p.y) != (points[j].pos.y > p.y)) &&
                    (p.x < (points[j].pos.x - points[i].pos.x) * (p.y - points[i].pos.y) / (points[j].pos.y - points[i].pos.y) + points[i].pos.x))
                    c = !c;
            }
            return c;
        }
        public static bool PointInPoly(MyVector2 p, List<Point2> points)
        {
            bool c = false;
            int n = points.Count;
            for (int i = 0, j = n - 1; i < n; j = i++)
            {
                if (((points[i].pos.y > p.y) != (points[j].pos.y > p.y)) &&
                    (p.x < (points[j].pos.x - points[i].pos.x) * (p.y - points[i].pos.y) / (points[j].pos.y - points[i].pos.y) + points[i].pos.x))
                    c = !c;
            }
            return c;
        }

        public static MyVector2 FindLinesegmentCircleIntersection(MyVector2 u, MyVector2 v, MyVector2 c, double radii)
        {
            if (!IsLineSegmentIntersectWithCircle(u, v, c, radii))
            {
                return u;
            }
            MyVector2 v1 = new MyVector2();
            MyVector2 v2 = new MyVector2();
            MyVector2 fp = u;
            if ((v - c).Length() > (u - c).Length())
            {
                fp = v;
            }
            // line 
            if (Math.Abs(u.x - v.x) < thresh)
            {
                double x = (u.x + v.x) / 2;
                double y = Math.Sqrt(radii * radii - Math.Pow(x - c.x, 2)) + c.y;
                double ny = -Math.Sqrt(radii * radii - Math.Pow(x - c.x, 2)) + c.y;
                v1 = new MyVector2(x, y);
                v2 = new MyVector2(x, ny);
            }
            else if (Math.Abs(u.y - v.y) < thresh)
            {
                double y = (u.x + v.y) / 2;
                double x = Math.Sqrt(radii * radii - Math.Pow(y - c.y, 2)) + c.x;
                double nx = -Math.Sqrt(radii * radii - Math.Pow(y - c.y, 2)) + c.x;
                v1 = new MyVector2(x, y);
                v2 = new MyVector2(nx, y);
            }
            else
            {
                double k = (u.y - v.y) / (u.x - v.x);
                double b = u.y - k * u.x;
                double constant = c.x * c.x + b * b + c.y * c.y - 2 * b * c.y;
                constant = radii * radii - constant;
                double coef1 = 1 + k * k;
                double coef2 = 2 * k * b - 2 * c.x - 2 * k * c.y;
                coef2 /= coef1;
                constant /= coef1;
                constant += coef2 * coef2 / 4;
                double x = Math.Sqrt(constant) - coef2 / 2;
                double nx = -Math.Sqrt(constant) - coef2 / 2;
                v1 = new MyVector2(x, k * x + b);
                v2 = new MyVector2(nx, k * nx + b);
            }
            if ((v1 - fp).Length() < (v2 - fp).Length())
            {
                return v1;
            }
            else
            {
                return v2;
            }
        }//FindLinesegmentCircleIntersection

        public static MyVector2 FindPointTolineFootPrint(MyVector2 pt, MyVector2 u, MyVector2 v)
        {
            MyVector2 uv = (v - u).Normalize();
            if (double.IsNaN(uv.x)) return pt;
            return u + (pt - u).Dot(uv) * uv;
        }

        public static bool IsLineSegmentIntersectWithCircle(MyVector2 u, MyVector2 v, MyVector2 c, double radius)
        {
            if ((u - c).Length() < radius || (v - c).Length() < radius) return true;
            MyVector2 uv = v - u;
            double r = (c - v).Dot(uv) / (u - v).Dot(uv);
            if (r < 0 || r > 1)
                return false;
            MyVector2 p = u * r + (1 - r) * v;
            return (p - c).Length() < radius;
        }
        // find Cumcentre of a triangle
        public static MyVector2 FindcirCumcentre(MyVector2 v1, MyVector2 v2, MyVector2 v3)
        {
            double a = (v2.x - v1.x) * 2, b = (v2.y - v1.y) * 2;
            double c = (v3.x - v1.x) * 2, d = (v3.y - v1.y) * 2;
            double p = v2.SquareLength() - v1.SquareLength();
            double q = v3.SquareLength() - v1.SquareLength();
            Matrix2d M = new Matrix2d(a, b, c, d);
            return M.Inverse() * new MyVector2(p, q);
        }

        public static bool IsSameOrientation(MyVector2 A, MyVector2 B, MyVector2 C)
        {
            return (C.y - A.y) * (B.x - A.x) > (B.y - A.y) * (C.x - A.x);
        }
        public static bool IsLineInterseted(MyVector2 p1, MyVector2 p2, MyVector2 q1, MyVector2 q2)
        {
            return IsSameOrientation(p1, q1, q2) != IsSameOrientation(p2, q1, q2) &&
                IsSameOrientation(p1, p2, q1) != IsSameOrientation(p1, p2, q2);
        }
        public static MyVector2 FindLinesIntersectPoint(MyVector2 p1, MyVector2 p2, MyVector2 q1, MyVector2 q2)
        {
            MyVector3 u1 = new MyVector3(p1, 1), v1 = new MyVector3(p2, 1);
            MyVector3 u2 = new MyVector3(q1, 1), v2 = new MyVector3(q2, 1);
            MyVector3 uv = u1.Cross(v1).Cross(u2.Cross(v2)).HomogenousNormalize();
            return uv.ToMyVector2();
        }

        public static bool IsLineCanvas3Intersect(MyVector3 u, MyVector3 v, MyVector3 n,
            Point3[] points, out MyVector3 ipoint)
        {
            MyVector3 dir = (u - v).Normalize();
            ipoint = new MyVector3();
            if (Math.Abs(dir.Dot(n)) < Utils.thresh)
            {
                MyVector3 mid = (u + v) / 2;
                if (PointInPoly3D(mid, points))
                    return true;
                else
                    return false;
            }

            // intersection point
            MyVector3 w = points[0].pos - v;
            MyVector3 d = u - v;
            double s = w.Dot(n) / d.Dot(n);
            ipoint = v + s * dir;
            if (PointInPoly3D(ipoint, points))
                return true;
            else
                return false;
        }

        public static bool PointInPoly3D(MyVector3 v, Point3[] points)
        {
            int N = points.Length;
            double sum = 0;
            for (int i = 0; i < N; ++i)
            {
                MyVector3 v1 = points[i].pos;
                MyVector3 v2 = points[(i + 1) % N].pos;
                v1 = (v1 - v).Normalize();
                v2 = (v2 - v).Normalize();
                double angle = Math.Acos(v1.Dot(v2));
                sum += angle;
            }
            if (Math.Abs(sum - Math.PI * 2) < 0.01)
                return true;
            return false;
        }

        // find bounding box of a set of 3d points
        public static CalibrateCuboid FindBoundingBox(List<MyVector3> points)
        {
            double minx = double.MaxValue, maxx = double.MinValue;
            double miny = double.MaxValue, maxy = double.MinValue;
            double minz = double.MaxValue, maxz = double.MinValue;
            foreach (MyVector3 p in points)
            {
                if (minx > p.x) minx = p.x;
                if (maxx < p.x) maxx = p.x;
                if (miny > p.y) miny = p.y;
                if (maxy < p.y) maxy = p.y;
                if (minz > p.z) minz = p.z;
                if (maxz < p.z) maxz = p.z;
            }
            Point3[] boxpoints = new Point3[8] {
				new Point3(new MyVector3(maxx,miny,minz)),
				new Point3(new MyVector3(maxx,maxy,minz)),
				new Point3(new MyVector3(minx,maxy,minz)),
				new Point3(new MyVector3(minx,miny,minz)),
				new Point3(new MyVector3(maxx,miny,maxz)),
				new Point3(new MyVector3(maxx,maxy,maxz)),
				new Point3(new MyVector3(minx,maxy,maxz)),
				new Point3(new MyVector3(minx,miny,maxz))
			};

            return new CalibrateCuboid(boxpoints);

        }
        public static Quad2D FindBoundingBox(List<MyVector2> points)
        {
            double minx = double.MaxValue, maxx = double.MinValue;
            double miny = double.MaxValue, maxy = double.MinValue;
            foreach (MyVector2 p in points)
            {
                if (minx > p.x) minx = p.x;
                if (maxx < p.x) maxx = p.x;
                if (miny > p.y) miny = p.y;
                if (maxy < p.y) maxy = p.y;
            }
            Point2[] boxpoints = new Point2[4] {
				new Point2(new MyVector2(maxx,miny)),
				new Point2(new MyVector2(maxx,maxy)),
				new Point2(new MyVector2(minx,maxy)),
				new Point2(new MyVector2(minx,miny)),
			};

            return new Quad2D(boxpoints.ToList());

        }
        public static Quad3D FindBoundingBox(MyVector3[] points)
        {
            double minx = double.MaxValue, maxx = double.MinValue;
            double miny = double.MaxValue, maxy = double.MinValue;
            double minz = double.MaxValue, maxz = double.MinValue;
            foreach (MyVector3 p in points)
            {
                if (minx > p.x) minx = p.x;
                if (maxx < p.x) maxx = p.x;
                if (miny > p.y) miny = p.y;
                if (maxy < p.y) maxy = p.y;
                if (minz > p.z) minz = p.z;
                if (maxz < p.z) maxz = p.z;
            }
            Point3[] boxpoints = new Point3[8] {
				new Point3(new MyVector3(maxx,miny, minz)),
				new Point3(new MyVector3(maxx,maxy, minz)),
				new Point3(new MyVector3(minx,maxy, minz)),
				new Point3(new MyVector3(minx,miny, minz)),
				new Point3(new MyVector3(maxx,miny, maxz)),
				new Point3(new MyVector3(maxx,maxy, maxz)),
				new Point3(new MyVector3(minx,maxy, maxz)),
				new Point3(new MyVector3(minx,miny, maxz))
			};

            return new Quad3D(boxpoints);
        }


        // color mall
        public static Color[] ColorMall = null;
        public static Random Random = new Random();
        public static void InitialRandomColors(int K)
        {
            Random rand = new Random();
            ColorMall = new Color[K];
            for (int i = 0; i < K; ++i)
            {
                ColorMall[i] = Color.FromArgb(
                    rand.Next(0, 255),
                    rand.Next(0, 255),
                    rand.Next(0, 255)
                );
            }

            // 8-data, color scheme
            if (true)
            {
                int off_set = 0;
                ColorMall[0 + off_set] = Color.FromArgb(228, 26, 28);
                ColorMall[1 + off_set] = Color.FromArgb(77, 175, 74);
                ColorMall[2 + off_set] = Color.FromArgb(55, 126, 184);
                //ColorMall[1 + off_set] = Color.FromArgb(55, 126, 184);
                //ColorMall[2 + off_set] = Color.FromArgb(77, 175, 74);

                ColorMall[3 + off_set] = Color.FromArgb(152, 78, 163);
                ColorMall[4 + off_set] = Color.FromArgb(255, 127, 0);
                ColorMall[5 + off_set] = Color.FromArgb(0, 128, 255);
                ColorMall[6 + off_set] = Color.FromArgb(166, 86, 40);
                ColorMall[7 + off_set] = Color.FromArgb(247, 129, 191);
            }
        }
        public static bool IsVectorCoincident(MyVector3 u, MyVector3 v)
        {
            double cos = u.Normalize().Dot(v.Normalize());
            return Math.Abs(cos) > 0.95;
        }
        public static void ComputeVanishingColor(PlanarShape3 shape)
        {
            MyVector3[] axes = new MyVector3[3] {
				new MyVector3(1,0,0),
				new MyVector3(0,1,0),
				new MyVector3(0,0,1)
			};
            for (int i = 0; i < 3; ++i)
            {
                double cos = Math.Abs(axes[i].Dot(shape.normal));
                if (cos > 0.8)
                {
                    shape.color = Utils.ColorMall[i];
                    break;
                }
            }
        }

        // others
        public static double PointToShapeDistance2D(MyVector2 point, PlanarShape3 shape)
        {
            int n = shape.points.Length;
            double mindis = double.MaxValue;
            for (int i = 0; i < n; ++i)
            {
                MyVector2 u = shape.points2[i].pos, v = shape.points2[(i + 1) % n].pos;
                double dis = PointTolineSegmentDistance2(point, u, v);
                if (dis < mindis)
                    mindis = dis;
            }
            return mindis;
        }
        public static double PointToShapeDistance(MyVector3 point, PlanarShape3 shape)
        {
            int n = shape.points.Length;
            double mindis = double.MaxValue;
            for (int i = 0; i < n; ++i)
            {
                MyVector3 u = shape.points[i].pos, v = shape.points[(i + 1) % n].pos;
                double dis = PointTolineSegmentDistance(point, u, v);
                if (dis < mindis)
                    mindis = dis;
            }
            return mindis;
        }
        public static double PointTolineSegmentDistance(MyVector3 p, MyVector3 u, MyVector3 v)
        {
            bool isProjectionIn = !((p - u).Dot(v - u) < 0 || (p - v).Dot(u - v) < 0);
            if (!isProjectionIn)
            {
                return Math.Min((p - u).Length(), (p - v).Length());
            }
            else
            {
                MyVector3 uv = (v - u).Normalize();
                MyVector3 up = p - u;
                return (up - up.Dot(uv) * uv).Length();
            }
        }
        public static double PointTolineSegmentDistance2(MyVector2 p, MyVector2 u, MyVector2 v)
        {
            bool isProjectionIn = !((p - u).Dot(v - u) < 0 || (p - v).Dot(u - v) < 0);
            if (!isProjectionIn)
            {
                return Math.Min((p - u).Length(), (p - v).Length());
            }
            else
            {
                MyVector2 uv = (v - u).Normalize();
                MyVector2 up = p - u;
                return (up - up.Dot(uv) * uv).Length();
            }
        }
        public static double PointToLineDistance(MyVector2 p, MyVector2 u, MyVector2 v)
        {
            if (Math.Abs(u.x - u.x) < Utils.thresh)
            {
                return Math.Abs(p.x - u.x);
            }
            else if (Math.Abs(u.y - v.y) < Utils.thresh)
            {
                return Math.Abs(p.y - u.y);
            }
            double k = (v.y - u.y) / (v.x - u.x);
            double b = u.y - k * u.x;
            return Math.Abs(k * p.x - p.y + b) / Math.Sqrt(1 + k * k);
        }
        public static MyMatrix3d ComputeHomographyMatrix(MyVector2[] imgpts, MyVector2[] refpts)
        {
            // compute mapping from imgpts to refpts
            if (refpts.Length != imgpts.Length || imgpts.Length < 4) return null;

            PointF[] source = new PointF[imgpts.Length];
            PointF[] target = new PointF[imgpts.Length];
            for (int i = 0; i < imgpts.Length; ++i)
            {
                source[i] = new PointF((float)imgpts[i].x, (float)imgpts[i].y);
                target[i] = new PointF((float)refpts[i].x, (float)refpts[i].y);
            }
            Matrix<double> mat = new Matrix<double>(4, 4);
            CvInvoke.FindHomography(source, target, mat, HomographyMethod.Ransac);

            MyMatrix3d H = new MyMatrix3d();
            for (int i = 0; i < 3; ++i)
            {
                for (int j = 0; j < 3; ++j)
                {
                    H[i, j] = mat[i, j];
                }
            }
            return H;
        }

        // compute the angle formed by u-o-v
        public static double ComputeTriAngle(MyVector2 u, MyVector2 o, MyVector2 v)
        {
            MyVector2 uo = (u - o).Normalize();
            MyVector2 vo = (v - o).Normalize();
            double cos = uo.Dot(vo);
            double angle = Math.Acos(cos);
            if (double.IsNaN(angle)) return 0;
            return angle;
        }
        public static double gauss(double x, double mean, double sigma)
        {
            return (Math.Exp((-(x - mean) * (x - mean)) / (2 * sigma * sigma)) / (Math.Sqrt(Math.PI * 2.0) * sigma));
        }
        public static bool IsSegmentsIntersecting(MyVector2 a, MyVector2 b, MyVector2 c, MyVector2 d, out double r1)
        {
            MyVector2 cd = c - d;
            MyVector2 ba = b - a;
            Matrix2d M = new Matrix2d(cd, ba);
            MyVector2 r = M.Inverse() * (b - d);
            r1 = r.x;
            return 0 < r.x && r.x < 1 && 0 < r.y && r.y < 1;
        }
        public static bool IsSegmentsIntersecting(MyVector2 a, MyVector2 b, MyVector2 c, MyVector2 d)
        {
            MyVector2 cd = c - d;
            MyVector2 ba = b - a;
            Matrix2d M = new Matrix2d(cd, ba);
            MyVector2 r = M.Inverse() * (b - d);
            return 0 < r.x && r.x < 1 && 0 < r.y && r.y < 1;
        }
        public static MyMatrix4d RotationMatrixU2V(MyVector3 u, MyVector3 v)
        {
            // rotate u to v
            u = u.Normalize();
            v = v.Normalize();
            MyVector3 rot_axis = u.Cross(v).Normalize();
            if (double.IsNaN(rot_axis.x)) return MyMatrix4d.IdentityMatrix();
            double cos = u.Dot(v);
            double rot_angle = Math.Acos(cos);
            return MyMatrix4d.RotationMatrix(rot_axis, rot_angle);
        }

        public static List<MyVector2> SmoothCurve(List<MyVector2> curve)
        {
            // smooth curves by checking the slope between neighboring points
            List<MyVector2> tempcurve = new List<MyVector2>();
            tempcurve.AddRange(curve);
            if (curve.Count() < 3)
            {
                return curve;
            }
            List<MyVector2> v = new List<MyVector2>();
            for (int it = 0; it < 2; it++)
            {
                v.Clear();
                List<double> shapeDs = computeShapeDescriptor(tempcurve);
                v.Add(tempcurve[0]);
                int count = 0;
                for (int i = 0; i < shapeDs.Count(); i++)
                {
                    if (shapeDs[i] > Math.PI / 4.0f && count > 0/* || asum > PI / 6.0f*/)
                    {
                        v.Add(tempcurve[i + 1]);
                        count = 0;
                    }
                    else
                    {
                        v.Add(tempcurve[i + 1] +
                            (0.5f * ((tempcurve[i] + tempcurve[i + 2]) / 2.0f - tempcurve[i + 1])));
                        count++;
                    }
                }
                v.Add(curve[curve.Count - 1]);

                tempcurve.Clear();
                tempcurve.AddRange(v);
            }
            return v;
        }
        private static List<double> computeShapeDescriptor(List<MyVector2> curve)
        {
            List<double> v = new List<double>();
            for (int i = 1; i < curve.Count - 1; i++)
            {
                MyVector2 e1 = curve[i] - curve[i - 1];
                MyVector2 e2 = curve[i + 1] - curve[i];
                e1 /= (e1.Length() + 1e-20f);
                e2 /= (e2.Length() + 1e-20f);
                double dotV = Math.Min(1.0, Math.Max(-1.0, e1.Dot(e2)));
                v.Add((double)Math.Acos(dotV));
            }
            return v;
        }
        public static MyVector2 GetIntersectPoint(MyVector2 u, MyVector2 v, MyVector2 p, MyVector2 q)
        {
            MyVector3 uu = new MyVector3(u, 1);
            MyVector3 vv = new MyVector3(v, 1);
            MyVector3 pp = new MyVector3(p, 1);
            MyVector3 qq = new MyVector3(q, 1);

            MyVector3 it = uu.Cross(vv).Cross(pp.Cross(qq));
            it.HomogenousNormalize();

            return it.ToMyVector2();
        }
    }

    interface Editable
    {
        void Transform(MyMatrix4d T);
    }

    public class CoordinateFrame
    {
        public MyVector3 old_o, old_x, old_y, old_z;
        public MyVector3 o, x, y, z;
        public double xBound, yBound;
        // xyzframe is the rotation matrix, transfer from world to local
        public MyMatrix3d XYZFrame = null, oldXYZFrame = null;
        public CoordinateFrame(MyVector3 o_, MyVector3 x_, MyVector3 y_, MyVector3 z_)
        {
            this.old_o = this.o = o_;
            this.old_x = this.x = x_;
            this.old_y = this.y = y_;
            this.old_z = this.z = z_;
            this.oldXYZFrame = new MyMatrix3d(x_, y_, z_);
            this.XYZFrame = new MyMatrix3d(x_, y_, z_);
        }
        // the bouding region
        public void SetXYBound(double xbound, double ybound)
        {
            this.xBound = xbound;
            this.yBound = ybound;
        }
        public void ShrinkToBound(ref MyVector3 p)
        {
            // if p exceeds the boundary of the region, bring it back
            MyVector3 local_p = this.GetPointLocalCoord(p);
            if (local_p.x < -xBound) local_p.x = -xBound;
            if (local_p.x > xBound) local_p.x = xBound;
            if (local_p.y < -yBound) local_p.y = -yBound;
            if (local_p.y > yBound) local_p.y = yBound;
            p = this.GetPointSpaceCoord(local_p);
        }
        public MyVector3 this[int index]
        {
            get
            {
                if (index == 0) return x;
                if (index == 1) return y;
                if (index == 2) return z;
                throw new ArgumentException();
            }
            set
            {
                if (index == 0) x = value;
                if (index == 1) y = value;
                if (index == 2) z = value;
            }
        }

        public MyVector3 GetPointLocalCoord(MyVector3 p)
        {
            MyVector3 po = (p - this.o);
            return new MyVector3(po.Dot(this.x), po.Dot(this.y), po.Dot(this.z));
        }
        public MyVector3 GetPointSpaceCoord(MyVector3 p)
        {
            return this.XYZFrame * p + this.o;
        }
        public MyVector3 PointAtXYCoord(MyVector2 coord)
        {
            return this.XYZFrame * new MyVector3(coord, 0) + this.o;
        }

        public void TransformOld(MyMatrix4d T)
        {
            this.o = (T * new MyVector4(this.old_o, 1)).XYZ();
            this.x = (T * new MyVector4(this.old_x, 0)).XYZ().Normalize();
            this.y = (T * new MyVector4(this.old_y, 0)).XYZ().Normalize();
            this.z = (T * new MyVector4(this.old_z, 0)).XYZ().Normalize();
            this.XYZFrame = new MyMatrix3d(this.x, this.y, this.z);
        }
        public void UpdateFrame()
        {
            this.old_o = this.o;
            this.old_x = this.x;
            this.old_y = this.y;
            this.old_z = this.z;
            this.oldXYZFrame = new MyMatrix3d(this.XYZFrame);
        }
    }




    public class PlanarShape3 : Editable, IComparable<PlanarShape3>
    {
        public Geometry3 HostGeometry_;
        public bool isClosed;
        public bool snapped;
        public bool clipped = false;
        public uint textureId;
        public double depth = -1;
        public Color color = Color.Gray;//Color.Cyan;//
        public MyVector3 center;
        public MyVector3 normal;
        public Point2[] points2;	// points in screen space
        public Point3[] points;		// points in 3d space
        public MyVector3[] pointsLocalCoord;		// points in local coordinate space
        public Point3[] edgePoints;
        public Point3[] skethPoints;
        public MyVector2[] texCoord;
        public PlanarShape3() { }
        public PlanarShape3(Point3[] points3)
        {
            this.points = points3;
            this.ComputeCenter();
            this.ComputeNormal();
            this.ComputePointsLocalCoord();
        }
        public CoordinateFrame GetCoordFrame()
        {
            MyVector3 u = (this.points[1].pos - this.points[0].pos).Normalize();
            MyVector3 v = this.normal.Cross(u).Normalize();
            double xlen = (this.points[1].pos - this.points[0].pos).Length();
            double ylen = (this.points[2].pos - this.points[1].pos).Length();
            CoordinateFrame frame = new CoordinateFrame(this.center, u, v, this.normal);
            frame.SetXYBound(xlen / 2, ylen / 2);
            return frame;
        }
        public void Transform(MyMatrix4d T)
        {
            foreach (Point3 pt in this.points)
            {
                pt.Transform(T);
            }
            this.ComputeCenter();
            this.ComputeNormal();
        }
        public void Transform_from_origin(MyMatrix4d T)
        {
            foreach (Point3 pt in this.points)
            {
                pt.Transform_from_origin(T);
            }
            this.ComputeCenter();
            this.ComputeNormal();
        }
        public void Transform_to_origin(MyMatrix4d T)
        {
            foreach (Point3 pt in this.points)
            {
                pt.Transform_to_origin(T);
            }
            this.ComputeCenter();
            this.ComputeNormal();
        }
        public void Scale(float scale_factor)
        {
            MyMatrix4d T1 = MyMatrix4d.TranslationMatrix(new MyVector3() - this.center);
            MyMatrix4d T2 = MyMatrix4d.TranslationMatrix(this.center);
            MyMatrix4d S = MyMatrix4d.ScalingMatrix(scale_factor, scale_factor, scale_factor);
            this.Transform(T2 * S * T1);
        }
        public void MakeGroundTouching()
        {
            double z_min = double.MaxValue;
            foreach (Point3 pt in this.points)
            {
                if (pt.pos.z < z_min)
                    z_min = pt.pos.z;
            }
            MyVector3 off = new MyVector3(0, 0, -z_min);
            this.Transform(MyMatrix4d.TranslationMatrix(off));
        }
        public void ComputeNormal()
        {
            MyVector3 nor = (this.points[2].pos - this.points[1].pos).Cross(this.points[1].pos - this.points[0].pos);
            this.normal = nor.Normalize();
        }
        public void ComputeCenter()
        {
            MyVector3 c = new MyVector3();
            foreach (Point3 pt in this.points)
            {
                c += pt.pos;
            }
            this.center = c / this.points.Length;
        }
        public double GetArea()
        {
            /* get the area of a planar 3d shape*/
            // A 2D polygon can be decomposed into triangles. For computing area, 
            // there is a very easy decomposition method for simple polygons (i.e. ones without self intersections). 
            // Let a polygon OMEGA be defined by its vertices Vi=(xi,yi) for i=0,n with Vn=V0. 
            // Also, let P be any point; and for each edge ei=ViVi+1 of OMEGA, form the triangle 
            // DELTA-i=DELTA-PViVi+1. Then, the area of OMEGA is equal to the sum of the signed areas 
            // of all the triangles DELTA-i for i=0,n-1; and we have:
            //if (this.pointsLocalCoord == null)
            this.ComputePointsLocalCoord();

            // assume a random point p is at (0,0);
            int n = this.points.Length;
            double sum = 0;
            for (int i = 0; i < n; ++i)
            {
                MyVector3 u = this.pointsLocalCoord[i], v = this.pointsLocalCoord[(i + 1) % n];
                double signed_area = u.x * v.y - u.y * v.x;
                sum += signed_area;
            }
            return Math.Abs(sum);
        }
        public PlanarShape3 Copy()
        {
            Point3[] newPoints = new Point3[this.points.Length];
            for (int i = 0; i < this.points.Length; ++i)
            {
                newPoints[i] = new Point3(this.points[i].pos);
            }

            PlanarShape3 shape = new PlanarShape3(newPoints);
            shape.points2 = new Point2[this.points2.Length];
            for (int i = 0; i < this.points.Length; ++i)
            {
                shape.points2[i] = new Point2(this.points2[i].pos);
            }
            shape.ComputeCenter();
            shape.ComputeNormal();
            return shape;
        }
        public void ComputeSketchLines()
        {
            int N = this.points.Length;
            int NK = 4;
            this.edgePoints = new Point3[N * 2];
            this.skethPoints = new Point3[NK * N * 2];
            for (int i = 0; i < N; ++i)
            {
                MyVector3 p0 = this.points[i].pos;
                MyVector3 p1 = this.points[(i + 1) % N].pos;
                MyVector3 dir = (p1 - p0).Normalize();
                double len = (p0 - p1).Length() / 10;
                MyVector3 ep0 = p0 - dir * len;
                MyVector3 ep1 = p1 + dir * len;
                this.edgePoints[i * 2] = new Point3(ep0);
                this.edgePoints[i * 2 + 1] = new Point3(ep1);
                MyVector3 vnorm = this.normal.Cross(dir).Normalize();
                MyVector3 center = (p0 + p1) / 2;
                for (int j = 0; j < NK; ++j)
                {
                    Random rand = new Random();
                    int degree = rand.Next() % 6;
                    double angle = degree / 180.0 * Math.PI;
                    MyMatrix4d mat = MyMatrix4d.RotationMatrix(this.normal, angle);
                    Point3 cp0 = new Point3(this.points[i].pos - center);
                    Point3 cp1 = new Point3(this.points[(i + 1) % N].pos - center);
                    MyVector3 rotP0 = cp0.CalTransform(mat);
                    MyVector3 rotP1 = cp1.CalTransform(mat);
                    rotP0 += center;
                    rotP1 += center;
                    dir = (rotP1 - rotP0).Normalize();
                    len = (p0 - p1).Length() / 30;
                    double rlen = rand.NextDouble();
                    rlen *= len;
                    MyVector3 sp0 = rotP0 + dir * rlen;
                    MyVector3 sp1 = rotP1 - dir * rlen;
                    rlen = rand.NextDouble() * len / 2;
                    if (rand.Next() % 2 == 0)
                    {
                        sp0 += vnorm * rlen;
                        sp1 += vnorm * rlen;
                    }
                    else
                    {
                        sp0 -= vnorm * rlen;
                        sp1 -= vnorm * rlen;
                    }
                    this.skethPoints[NK * 2 * i + 2 * j] = new Point3(sp0);
                    this.skethPoints[NK * 2 * i + 2 * j + 1] = new Point3(sp1);
                }
            }
        }//ComputeSketchLines
        public void ComputeImageTextureCoord(int imgwidth, int imgheight)
        {
            if (this.points2 == null) return;
            int N = this.points2.Length;
            this.texCoord = new MyVector2[N];
            for (int i = 0; i < N; ++i)
            {
                this.texCoord[i] = new MyVector2(this.points2[i].pos.x / imgwidth,
                    this.points2[i].pos.y / imgheight);
            }
        }
        private void ComputePointsLocalCoord()
        {
            CoordinateFrame frame = this.GetCoordFrame();
            int n = this.points.Length;
            this.pointsLocalCoord = new MyVector3[n];
            for (int i = 0; i < n; ++i)
            {
                this.pointsLocalCoord[i] = frame.GetPointLocalCoord(this.points[i].pos);
            }
        }

        public void DrawBlended(Color c, byte opacity = 30)
        {
            GL.Disable(EnableCap.DepthTest);
            GL.Enable(EnableCap.Blend);
            GL.BlendFunc(BlendingFactorSrc.SrcAlpha, BlendingFactorDest.OneMinusSrcAlpha);

            GL.PolygonMode(MaterialFace.FrontAndBack, PolygonMode.Fill);

            // draw transparent quad
            if (this.isClosed)
            {
                GL.Color4(c.R, c.G, c.B, opacity);
                GL.Begin(PrimitiveType.Polygon);
                foreach (Point3 pt3 in this.points)
                    GL.Vertex3(pt3.pos.x, pt3.pos.y, pt3.pos.z);
                GL.End();
            }

            // draw outlines
            GL.Enable(EnableCap.LineSmooth);
            GL.Color4((byte)32, (byte)32, (byte)32, (byte)120);
            GL.LineWidth(1.0f);
            GL.Begin(PrimitiveType.Lines);
            int l = this.isClosed ? this.points.Length : this.points.Length - 1;
            for (int i = 0; i < l; ++i)
            {
                Point3 p = this.points[i], q = this.points[(i + 1) % this.points.Length];
                GL.Vertex3(p.pos.x, p.pos.y, p.pos.z);
                GL.Vertex3(q.pos.x, q.pos.y, q.pos.z);
            }
            GL.End();
            //	GL.Disable(EnableCap.Blend);
        }
        public void Draw(Color c)
        {
            GL.PolygonMode(MaterialFace.Front, PolygonMode.Fill);

            // draw transparent quad
            if (this.isClosed)
            {
                GL.Color3(c.R, c.G, c.B);
                GL.Begin(PrimitiveType.Polygon);
                foreach (Point3 pt3 in this.points)
                    GL.Vertex3(pt3.pos.x, pt3.pos.y, pt3.pos.z);
                GL.End();
            }

            // draw outlines
            GL.Enable(EnableCap.LineSmooth);
            GL.Color3(128, 128, 128);
            GL.LineWidth(1.0f);
            GL.Begin(PrimitiveType.Lines);
            int l = this.isClosed ? this.points.Length : this.points.Length - 1;
            for (int i = 0; i < l; ++i)
            {
                Point3 p = this.points[i], q = this.points[(i + 1) % this.points.Length];
                GL.Vertex3(p.pos.x, p.pos.y, p.pos.z);
                GL.Vertex3(q.pos.x, q.pos.y, q.pos.z);
            }
            GL.End();
        }
        public void DrawShaded()
        {
            GL.Enable(EnableCap.Normalize);

            MyVector3 nor = this.normal;
            GL.Begin(PrimitiveType.Polygon);
            foreach (Point3 pt3 in this.points)
            {
                GL.Normal3(nor.x, nor.y, nor.z);
                GL.Vertex3(pt3.pos.x, pt3.pos.y, pt3.pos.z);
            }
            GL.End();

            GL.Disable(EnableCap.Normalize);
        }
        public void Draw()
        {
            GL.Begin(PrimitiveType.Polygon);
            foreach (Point3 pt3 in this.points)
            {
                GL.Vertex3(pt3.pos.x, pt3.pos.y, pt3.pos.z);
            }
            GL.End();
        }
        public void DrawOutLine(Color c, byte opacity = 255)
        {
            GL.Color4(c.R, c.G, c.B, opacity);
            GL.Begin(PrimitiveType.Lines);
            int N = this.points.Length;
            for (int i = 0; i < N; ++i)
            {
                Point3 p = this.points[i], q = this.points[(i + 1) % N];
                GL.Vertex3(p.pos.x, p.pos.y, p.pos.z);
                GL.Vertex3(q.pos.x, q.pos.y, q.pos.z);
            }
            GL.End();
        }
        public void DrawOutLinePlus(Color c, double ext = 0.1, byte opacity = 255)
        {
            // draw pencil sketchy style lines
            GL.Begin(PrimitiveType.Lines);
            int N = this.points.Length;
            // draw original
            GL.Color3(c.R, c.G, c.B);
            for (int i = 0; i < N; ++i)
            {
                MyVector3 p = this.points[i].pos, q = this.points[(i + 1) % N].pos;
                GL.Vertex3(p.x, p.y, p.z);
                GL.Vertex3(q.x, q.y, q.z);
            }
            // draw extended
            Color c2 = Color.FromArgb(Math.Min(c.R + 30, 255), Math.Min(c.G + 30, 255), Math.Min(c.B + 30, 255));
            GL.Color3(c2.R, c2.G, c2.B);
            for (int i = 0; i < N; ++i)
            {
                MyVector3 p = this.points[i].pos, q = this.points[(i + 1) % N].pos;
                MyVector3 pq = (q - p).Normalize();
                MyVector3 u = p - ext * pq;
                MyVector3 v = q + ext * pq;
                GL.Vertex3(u.x, u.y, u.z);
                GL.Vertex3(v.x, v.y, v.z);
            }
            GL.End();
        }
        public void DrawPencilSketchLine(Color c, byte opacity = 255)
        {
            GL.Color4(c.R, c.G, c.B, opacity);
            // Canvas3 edge
            GL.Begin(PrimitiveType.Lines);
            int N = this.edgePoints.Length / 2;
            for (int i = 0; i < N; ++i)
            {
                Point3 p = this.edgePoints[2 * i], q = this.edgePoints[2 * i + 1];
                GL.Vertex3(p.pos.x, p.pos.y, p.pos.z);
                GL.Vertex3(q.pos.x, q.pos.y, q.pos.z);
            }
            GL.End();
            // pencil sketch
            N = this.skethPoints.Length / 2;
            GL.Begin(PrimitiveType.Lines);
            for (int i = 0; i < N; ++i)
            {
                Point3 p = this.skethPoints[2 * i], q = this.skethPoints[2 * i + 1];
                GL.Vertex3(p.pos.x, p.pos.y, p.pos.z);
                GL.Vertex3(q.pos.x, q.pos.y, q.pos.z);
            }
            GL.End();
        }
        public void DrawTextured()
        {
            if (this.texCoord == null) return;

            GL.TexParameter(TextureTarget.Texture2D, TextureParameterName.TextureMagFilter, (int)TextureMagFilter.Nearest);

            // draw transparent quad
            GL.Enable(EnableCap.Texture2D);
            GL.BindTexture(TextureTarget.Texture2D, this.textureId);
            GL.Color4(255, 255, 255, 255); // white
            GL.Begin(PrimitiveType.Triangles);
            {
                for (int i = 0; i < this.texCoord.Length; ++i)
                {
                    GL.TexCoord2(this.texCoord[i].x, this.texCoord[i].y);
                    GL.Vertex3(this.points[i].pos.x, this.points[i].pos.y, this.points[i].pos.z);
                }
            }
            GL.End();
            GL.Disable(EnableCap.Texture2D);
        }
        public void DrawTextureBlended(Color c, byte opacity = 255)
        {
            if (this.texCoord == null) return;

            GL.Enable(EnableCap.Blend);
            GL.BlendFunc(BlendingFactorSrc.SrcAlpha, BlendingFactorDest.OneMinusSrcAlpha);
            GL.TexParameter(TextureTarget.Texture2D, TextureParameterName.TextureMagFilter, (int)TextureMagFilter.Nearest);

            // draw transparent quad
            GL.Enable(EnableCap.Texture2D);
            GL.BindTexture(TextureTarget.Texture2D, this.textureId);
            GL.Color4(c.R, c.G, c.B, opacity); // white
            GL.Begin(PrimitiveType.Triangles);
            {
                for (int i = 0; i < this.texCoord.Length; ++i)
                {
                    GL.TexCoord2(this.texCoord[i].x, this.texCoord[i].y);
                    GL.Vertex3(this.points[i].pos.x, this.points[i].pos.y, this.points[i].pos.z);
                }
            }
            GL.End();
            GL.Disable(EnableCap.Texture2D);

            GL.Disable(EnableCap.Blend);
        }

        public static PlanarShape3 UnitCircle3D(int sample_rate = 100)
        {
            double radius = 1.0;
            double d_theta = Math.PI * 2 / sample_rate;
            List<Point3> points3 = new List<Point3>();
            for (int i = 0; i < sample_rate; ++i)
            {
                double theta = i * d_theta;
                double x = radius * Math.Cos(theta);
                double y = radius * Math.Sin(theta);
                MyVector3 p = new MyVector3(x, y, 0);
                points3.Add(new Point3(p));
            }
            return new PlanarShape3(points3.ToArray());
        }
        public static bool IsSameShape(PlanarShape3 shape1, PlanarShape3 shape2)
        {
            if (shape1.points.Length != shape2.points.Length) return false;
            int N = shape1.points.Length;
            for (int i = 0; i < N; ++i)
            {
                if ((shape1.points[i].pos - shape2.points[i].pos).Length() > 1e-3)
                    return false;
            }
            return true;
        }
        // Implement the generic CompareTo method with the Temperature  
        // class as the Type parameter.  
        // 
        public int CompareTo(PlanarShape3 other)
        {
            // If other is not a valid object reference, this instance is greater. 
            if (other == null) return 1;

            // The temperature comparison depends on the comparison of  
            // the underlying Double values.  
            return -this.depth.CompareTo(other.depth);
        }
    }	// planar shapes in 3d
    public class Point3 : Editable
    {
        public Color color;
        public MyVector3 pos;
        public MyVector3 origin_pos;
        public Point3(double x, double y, double z)
        {
            this.pos = new MyVector3(x, y, z);
            this.origin_pos = new MyVector3(x, y, z);
        }
        public Point3(MyVector3 p)
        {
            this.pos = new MyVector3(p);
            this.origin_pos = new MyVector3(p);
        }
        public void Transform(MyMatrix4d T)
        {
            MyVector3 qt = (T * new MyVector4(this.pos, 1)).XYZ();
            this.pos = qt;
        }
        public void Transform_from_origin(MyMatrix4d T)
        {
            MyVector3 qt = (T * new MyVector4(this.origin_pos, 1)).XYZ();
            this.pos = qt;
        }
        public void Transform_to_origin(MyMatrix4d T)
        {
            MyVector3 qt = (T * new MyVector4(this.pos, 1)).XYZ();
            this.origin_pos = qt;
        }
        public MyVector3 CalTransform(MyMatrix4d T)
        {
            MyVector3 qt = (T * new MyVector4(this.pos, 1)).XYZ();
            return qt;
        }
        public object Clone()
        {
            Point3 another = new Point3(this.pos);
            return another;
        }
        public void Draw(float size = 1.0f)
        {
            GL.Color3(color.R, color.G, color.B);
            GL.PointSize(size);
            GL.Begin(PrimitiveType.Points);
            GL.Vertex3(this.pos.x, this.pos.y, this.pos.z);
            GL.End();
            GL.PointSize(1.0f);
        }
        public void Draw(Color c, float size = 1.0f)
        {
            GL.Color3(c.R, c.G, c.B);
            GL.PointSize(size);
            GL.Begin(PrimitiveType.Points);
            GL.Vertex3(this.pos.x, this.pos.y, this.pos.z);
            GL.End();
            GL.PointSize(1.0f);
        }
    }
    public class Quad3D : PlanarShape3, ICloneable
    {
        public Quad3D()
        {
            this.isClosed = true;
        }
        public Quad3D(Point3[] pts)
        {
            this.isClosed = true;
            this.points = pts;
            this.ComputeCenter();
            this.ComputeNormal();
        }
        public static Quad3D UnitQuad()		// x-z plance
        {
            double r = 1.0 / 2, h = 0;
            Point3[] points = new Point3[4] {
				new Point3(-r,h,-r), new Point3(-r, h, r), new Point3(r, h, r), new Point3(r, h,-r)
			};

            Quad3D quad = new Quad3D(points);
            return quad;
        }
        public object Clone()
        {
            Quad3D quad = new Quad3D(new Point3[4] {
				new Point3(this.points[0].pos), new Point3(this.points[1].pos),
				new Point3(this.points[2].pos), new Point3(this.points[3].pos)
				}
            );
            quad.ComputeCenter();
            quad.ComputeNormal();
            return quad;
        }
        public double DiagnalLength()
        {
            return (this.points[6].pos - this.points[0].pos).Length();
        }
        public bool InQuad(MyVector3 p)
        {
            List<MyVector3> vertices = new List<MyVector3>();
            for (int i = 0; i < this.points.Count(); i++)
                vertices.Add(points[i].pos);

            double theta = 0;
            MyVector3 v1 = p - vertices[0];
            MyVector3 v2 = p - vertices[vertices.Count - 1];
            theta = Math.Acos(v1.Dot(v2) / (v1.Length() * v2.Length()));
            for (int i = 0; i < vertices.Count - 1; i++)
            {
                v1 = p - vertices[i];
                v2 = p - vertices[i + 1];
                theta += Math.Acos(v1.Dot(v2) / (v1.Length() * v2.Length()));
            }

            if (Math.Abs(2 * Math.PI - theta) < 1e-10)
                return true;
            else
                return false;
        }
    }
    public class Polygon3D : PlanarShape3
    {
        public Polygon3D()
        {
        }
        public Polygon3D(Point3[] pts)
        {
            this.points = pts;
        }
        public static Polygon3D UnitCircle()
        {
            double r = 1.0 / 2;
            int slices = 20;
            double dtheta = Math.PI * 2 / slices;
            double z = 0;
            List<Point3> circularPoints = new List<Point3>();
            for (int i = 0; i < slices; ++i)
            {
                double theta = dtheta * i;
                double x = r * Math.Cos(theta);
                double y = r * Math.Sin(theta);
                circularPoints.Add(new Point3(new MyVector3(x, y, z)));
            }

            Polygon3D plg = new Polygon3D(circularPoints.ToArray());
            plg.isClosed = true;

            return plg;
        }
    }
    public class RotationalCircleMetaphor : PlanarShape3
    {
        private double radius;
        private int sampleRate = 100;
        private int startIndex = -1, endIndex = -1;
        public RotationalCircleMetaphor(MyVector3 center, MyVector3 axis, double radii, int sample_rate = 100)
        {
            PlanarShape3 shape = PlanarShape3.UnitCircle3D(sample_rate);
            MyMatrix4d S = MyMatrix4d.ScalingMatrix(radii, radii, radii);
            MyMatrix4d R = MyMatrix4d.RotationMatrixU2V(shape.normal, axis);
            MyMatrix4d T = MyMatrix4d.TranslationMatrix(center);
            shape.Transform(T * R * S);
            this.points = shape.points.Clone() as Point3[];
            this.ComputeCenter();
            this.ComputeNormal();
            this.radius = radii;
            this.sampleRate = sample_rate;
        }
        public RotationalCircleMetaphor(PlanarShape3 shape, int sample_rate = 100)
        {
            double r1 = (shape.points[1].pos - shape.points[0].pos).Length();
            double r2 = (shape.points[2].pos - shape.points[1].pos).Length();
            double r = Math.Max(r1, r2) / 2;
            PlanarShape3 tmp_shape = PlanarShape3.UnitCircle3D(sample_rate);
            MyMatrix4d S = MyMatrix4d.ScalingMatrix(r, r, r);
            MyMatrix4d R = MyMatrix4d.RotationMatrixU2V(tmp_shape.normal, shape.normal);
            MyMatrix4d T = MyMatrix4d.TranslationMatrix(shape.center);
            tmp_shape.Transform(T * R * S);
            this.points = tmp_shape.points.Clone() as Point3[];
            this.ComputeCenter();
            this.ComputeNormal();
            this.radius = r;
            this.sampleRate = sample_rate;
        }
        public MyMatrix4d GetRotationMatrix(MyVector3 u, MyVector3 v)
        {
            this.startIndex = this.SelectPointIndex(u);
            this.endIndex = this.SelectPointIndex(v);
            double angle = (this.startIndex - this.endIndex) * Math.PI * 2 / this.sampleRate;
            MyMatrix4d R = MyMatrix4d.RotationMatrix(this.normal, angle);
            MyMatrix4d T1 = MyMatrix4d.TranslationMatrix(new MyVector3() - this.center);
            MyMatrix4d T2 = MyMatrix4d.TranslationMatrix(this.center);
            return T2 * R * T1;
        }
        private int SelectPointIndex(MyVector3 u)
        {
            double mindis = double.MaxValue;
            int index = 0, selected_index = -1;
            foreach (Point3 pt in this.points)
            {
                double dis = (pt.pos - u).Length();
                if (dis < mindis)
                {
                    selected_index = index;
                    mindis = dis;
                }
                index++;
            }
            return selected_index;
        }
        public void DrawMetaphorPoints()
        {
            if (this.startIndex != -1)
                this.points[this.startIndex].Draw(Color.Lime, 7.0f);
            if (this.endIndex != -1)
                this.points[this.endIndex].Draw(Color.Red, 7.0f);
        }
    }
    public class Point2 : ICloneable
    {
        public int index = -1;
        public Color color = Color.Tan;
        public MyVector2 pos, oldpos, txpos; // texture position
        public Point2(MyVector2 pt)
        {
            this.txpos = this.oldpos = this.pos = pt;
        }
        public Point2()
        {
            this.txpos = this.oldpos = this.pos = new MyVector2();
        }
        public void Translate(MyVector2 off)
        {
            this.pos = this.oldpos + off;
        }
        public void UpdateOldPos()
        {
            this.oldpos = this.pos;
        }
        public object Clone()
        {
            Point2 another = new Point2(this.pos);
            another.index = this.index;
            return another;
        }
        public void Draw(float size = 1.0f)
        {
            GL.Color3(color.R, color.G, color.B);
            GL.PointSize(size);
            GL.Begin(PrimitiveType.Points);
            GL.Vertex2(this.pos.x, this.pos.y);
            GL.End();
            GL.PointSize(1.0f);
        }
        public void Draw(Color c, float size = 1.0f)
        {
            GL.Color3(c.R, c.G, c.B);
            GL.PointSize(size);
            GL.Begin(PrimitiveType.Points);
            GL.Vertex2(this.pos.x, this.pos.y);
            GL.End();
            GL.PointSize(1.0f);
        }
    }
    public class LineSegment2
    {
        public Color color;
        public Point2 u, v;
        public Point3 u3, v3; // for rotation view
        public LineSegment2(Point2 a, Point2 b)
        {
            this.u = a;
            this.v = b;
        }
        public void Translate(MyVector2 off)
        {
            this.u.Translate(off);
            this.v.Translate(off);
        }
        public void UpdatePos()
        {
            this.u.UpdateOldPos();
            this.v.UpdateOldPos();
        }
        public MyVector2 Center()
        {
            return (this.u.pos + this.v.pos) / 2;
        }
        public MyVector2 Dir()
        {
            return (this.v.pos - this.u.pos).Normalize();
        }
        public double Length()
        {
            return (this.v.pos - this.u.pos).Length();
        }
        public bool IsEndingPointsCoIncident(LineSegment2 other, double thres)
        {
            double[] dist = new double[4] {
				(this.u.pos - other.u.pos).Length(),
				(this.u.pos - other.v.pos).Length(),
				(this.v.pos - other.u.pos).Length(),
				(this.v.pos - other.v.pos).Length()
			};
            foreach (double dis in dist)
            {
                if (dis < thres)
                    return true;
            }
            return false;
        }
        public void Draw(Color c, float linewidth = 1.0f)
        {
            GL.Color3(c.R, c.G, c.B);
            GL.LineWidth(linewidth);
            GL.Begin(PrimitiveType.Lines);
            GL.Vertex2(this.u.pos.x, this.u.pos.y);
            GL.Vertex2(this.v.pos.x, this.v.pos.y);
            GL.End();
            GL.LineWidth(1.0f);
        }
        public void Draw(float linewidth = 1.0f)
        {
            GL.Color3(this.color.R, this.color.G, this.color.B);
            GL.LineWidth(linewidth);
            GL.Begin(PrimitiveType.Lines);
            GL.Vertex2(this.u.pos.x, this.u.pos.y);
            GL.Vertex2(this.v.pos.x, this.v.pos.y);
            GL.End();
            GL.LineWidth(1.0f);
        }
        public Object Clone()
        {
            Point2 uc = u.Clone() as Point2;
            Point2 vc = v.Clone() as Point2;
            LineSegment2 line = new LineSegment2(uc, vc);
            if (this.u3 != null && this.v3 != null)
            {
                line.u3 = this.u3.Clone() as Point3;
                line.v3 = this.v3.Clone() as Point3;
            }
            line.color = this.color;
            return line;
        }
        public double ComputeT(MyVector2 p)
        {
            MyVector2 dir = u.pos - v.pos;
            dir = dir.Normalize();
            double t = (p - u.pos).Dot(dir) / dir.SquareLength();
            return t;
        }
        public double DistanceToLine(MyVector2 point)
        {
            // v = y2 u = y1
            double a = v.pos.y - u.pos.y;
            double b = u.pos.x - v.pos.x;
            double c = v.pos.x * u.pos.y - u.pos.x * v.pos.y;

            return Math.Abs((a * point.x + b * point.y + c) / Math.Sqrt(a * a + b * b));
        }
        public MyVector2 ProjToLine(MyVector2 startp, MyVector2 p)
        {
            MyVector2 dir = (this.u.pos - this.v.pos).Normalize();
            MyVector2 point = startp;

            double t = (p - point).Dot(dir) / dir.SquareLength();
            return point + dir * t;
        }//liyuwei - 1217

    }
    public class Line3
    {
        public MyVector3 startpoint, dir;
        public MyVector3 endpoint;
        private double length = -1;
        public Line3() { }
        public Line3(MyVector3 point_, MyVector3 dir_)
        {
            this.startpoint = point_;
            this.dir = dir_;
        }

        public Line3(double length, MyVector3 point_, MyVector3 dir_)
        {
            this.startpoint = point_;
            this.length = length;


            this.dir = dir_;
        }

        public Line3(MyVector3 startpoint, MyVector3 endpoint, bool tag)
        {
            this.startpoint = startpoint;
            this.endpoint = endpoint;
            this.dir = (endpoint - startpoint).Normalize();
        }
        public double DistanceToLine(MyVector3 pos)
        {
            double dis = (pos - startpoint).SquareLength() - Math.Pow(dir.Dot(pos - startpoint), 2.0) / dir.SquareLength();
            dis = Math.Sqrt(dis);
            return dis;
        }
        public MyVector3 ProjectToLine(MyVector3 p)
        {
            double t = (p - startpoint).Dot(dir) / dir.SquareLength();
            return startpoint + dir * t;
        }

        public Line3 Copy()
        {
            return new Line3(this.startpoint, this.dir);
        }

        public void SetPoints(MyVector3 start, MyVector3 end)
        {
            this.startpoint = start;
            this.endpoint = end;
            this.dir = (end - start).Normalize();
        }
        public void SetPoints(double start, double end)
        {
            this.startpoint = startpoint + start * dir;
            this.endpoint = startpoint + end * dir;
        }

        public void Draw(Color color, float p = 2.0f, double linelength = 1)
        {

            linelength *= 0.01;

            if (this.length != -1)
                linelength = this.length;


            GL.Color4(color.R, color.G, color.B, (byte)150);
            GL.LineWidth(p);
            GL.Begin(PrimitiveType.Lines);
            if (!endpoint.IsNull())
            {
                GL.Vertex3(startpoint.x, startpoint.y, startpoint.z);
                GL.Vertex3(endpoint.x, endpoint.y, endpoint.z);
            }
            else
            {
                GL.Color4(color.R, color.G, color.B, (byte)200);
                GL.Vertex3(startpoint.x - dir.x * linelength, startpoint.y - dir.y * linelength, startpoint.z - dir.z * linelength);
                GL.Vertex3(startpoint.x + dir.x * linelength, startpoint.y + dir.y * linelength, startpoint.z + dir.z * linelength);
            }
            GL.End();
            GL.LineWidth(1.0f);


        }


        public double ComputeT(MyVector3 p)
        {
            //project to line
            double t = (p - startpoint).Dot(dir) / dir.SquareLength();
            return t;

            //double tx = (p.x - point.x) / dir.x;
            //double ty = (p.y - point.y) / dir.y;
            //double tz = (p.z - point.z) / dir.z;
            //return (tx + ty + tz) / 3;
        }
        public MyVector3 GetPointwithT(double t)
        {
            return startpoint + dir * t;
        }
        public double Length()
        {
            if (this.length != -1) return this.length;
            return (this.startpoint - this.endpoint).Length();
        }
        public double LineSquareLength()
        {
            return (this.startpoint - this.endpoint).SquareLength();
        }
        public void Scale(double scale)
        {
            double t = (this.endpoint.x - this.startpoint.x) / this.dir.x;
            t *= scale;
            endpoint = this.startpoint + this.dir * t;

        }

        static public bool IsParallel(Line3 a, Line3 b, double angle = 15)
        {
            if (Math.Abs(a.dir.Dot(b.dir)) >= Math.Cos(angle * Math.PI / 180))
                return true;
            else
                return false;
        }
    }

    public class Quad2D : Editable
    {
        private uint textureId_;
        public List<Point2> points_ = new List<Point2>();
        // A-B
        // C-D
        public Point2 A_, B_, C_, D_;
        public Quad2D()
        {

        }
        public Quad2D(List<Point2> points)
        {
            if (points.Count < 4) return;

            this.points_ = points;
        }
        public Quad2D(Point2 a, Point2 b, Point2 c, Point2 d)
        {
            this.points_.Add(a);
            this.points_.Add(b);
            this.points_.Add(c);
            this.points_.Add(d);
        }
        public void BindTexture(uint txtid, int w, int h)
        {
            this.textureId_ = txtid;

            // compute texture coordinates
            foreach (Point2 p2 in this.points_)
            {
                p2.txpos.x /= w;
                p2.txpos.y /= h;
            }
        }
        public bool ContainsPoint(MyVector2 pt)
        {
            return Utils.PointInPoly(pt, this.points_);
        }
        public void Transform(MyMatrix4d T)
        {
            foreach (Point2 pt in this.points_)
            {
                MyVector4 pt4 = new MyVector4(pt.pos);
                pt4.w = 1;
                pt.pos = (T * pt4).XYZ().ToMyVector2();
            }
        }
        public MyVector2 GetCenter()
        {
            return (this.points_[0].pos + this.points_[2].pos) / 2;
        }
        public void Draw(Color c)
        {
            GL.Enable(EnableCap.Blend);
            GL.BlendFunc(BlendingFactorSrc.SrcAlpha, BlendingFactorDest.OneMinusSrcAlpha);

            GL.PolygonMode(MaterialFace.FrontAndBack, PolygonMode.Fill);

            // draw transparent quad
            GL.Color4(c.R, c.G, c.B, (byte)5);
            GL.Begin(PrimitiveType.Polygon);
            foreach (Point2 pt2 in this.points_)
                GL.Vertex2(pt2.pos.x, pt2.pos.y);
            GL.End();


            // draw outlines
            this.DrawOutline(Color.Red);

            GL.Disable(EnableCap.Blend);
        }
        public void DrawWithTexture(float opacity)
        {
            GL.Enable(EnableCap.Blend);
            GL.BlendFunc(BlendingFactorSrc.SrcAlpha, BlendingFactorDest.OneMinusSrcAlpha);

            GL.Enable(EnableCap.Texture2D);
            GL.BindTexture(TextureTarget.Texture2D, this.textureId_);
            GL.Color4(1, 1, 1, opacity);
            GL.Begin(PrimitiveType.Polygon);
            foreach (Point2 pt in this.points_)
            {
                GL.TexCoord2(pt.txpos.x, pt.txpos.y);
                GL.Vertex2(pt.pos.x, pt.pos.y);
            }
            GL.End();
            GL.Disable(EnableCap.Texture2D);

            this.DrawOutline(Color.Gray);

            GL.Disable(EnableCap.Blend);
        }
        public void DrawOutline(Color c)
        {
            // draw outlines
            GL.Color3(c.R, c.G, c.B);
            GL.LineWidth(1.0f);
            GL.Begin(PrimitiveType.Lines);
            int l = this.points_.Count;
            for (int i = 0; i < l; ++i)
            {
                Point2 pt = this.points_[i], qt = this.points_[(i + 1) % l];
                GL.Vertex2(pt.pos.x, pt.pos.y);
                GL.Vertex2(qt.pos.x, qt.pos.y);
            }
            GL.End();
        }
        public static bool IsIntersecting(Quad2D quad1, Quad2D quad2)
        {
            // corner inside rect;
            foreach (Point2 p in quad1.points_)
                if (Utils.PointInPoly(p.pos, quad2.points_))
                    return true;
            foreach (Point2 p in quad2.points_)
                if (Utils.PointInPoly(p.pos, quad1.points_))
                    return true;
            // segments intersecting
            for (int i = 0; i < 4; ++i)
            {
                MyVector2 a = quad1.points_[i].pos, b = quad1.points_[(i + 1) % 4].pos;
                for (int j = 0; j < 4; ++j)
                {
                    MyVector2 c = quad2.points_[i].pos, d = quad2.points_[(i + 1) % 4].pos;
                    if (Utils.IsSegmentsIntersecting(a, b, c, d))
                        return true;
                }
            }

            return false;
        }
    }
    public class Polygon2 : Editable
    {
        private uint textureId_;
        public List<Point2> points_ = new List<Point2>();
        public Polygon2()
        {

        }
        public Polygon2(List<Point2> points)
        {
            this.points_.AddRange(points_);
        }
        public void BindTexture(uint txtid, int w, int h)
        {
            this.textureId_ = txtid;
            // compute texture coordinates
            foreach (Point2 pt2 in this.points_)
            {
                pt2.txpos.x /= w;
                pt2.txpos.y /= h;
            }
        }
        public int GetPointCount()
        {
            return this.points_.Count;
        }
        public void AddPoint(MyVector2 pt)
        {
            this.points_.Add(new Point2(pt));
        }
        public void PopPoint()
        {
            int N = this.points_.Count;
            if (N > 0)
            {
                this.points_.RemoveAt(N - 1);
            }
        }
        public bool ContainsPoint(MyVector2 pt)
        {
            return Utils.PointInPoly(pt, this.points_);
        }
        public void Transform(MyMatrix4d T)
        {
            foreach (Point2 pt in this.points_)
            {
                MyVector4 pt4 = new MyVector4(pt.pos);
                pt4.w = 1;
                pt.pos = (T * pt4).XYZ().ToMyVector2();
            }
        }
        public Quad2D GetBoundingBox()
        {
            List<MyVector2> MyVector2 = new List<MyVector2>();
            foreach (Point2 pt in this.points_)
                MyVector2.Add(pt.pos);

            return Utils.FindBoundingBox(MyVector2);
        }
        public void Draw(Color c)
        {
            GL.Enable(EnableCap.Blend);
            GL.BlendFunc(BlendingFactorSrc.SrcAlpha, BlendingFactorDest.OneMinusSrcAlpha);

            GL.PolygonMode(MaterialFace.FrontAndBack, PolygonMode.Fill);

            // draw transparent quad
            GL.Color4(c.R, c.G, c.B, (byte)5);
            GL.Begin(PrimitiveType.Polygon);
            foreach (Point2 pt2 in this.points_)
                GL.Vertex2(pt2.pos.x, pt2.pos.y);
            GL.End();


            // draw outlines
            this.DrawOutline(Color.Red);

            GL.Disable(EnableCap.Blend);
        }
        public void DrawWithTexture(float opacity)
        {
            GL.Enable(EnableCap.Blend);
            GL.BlendFunc(BlendingFactorSrc.SrcAlpha, BlendingFactorDest.OneMinusSrcAlpha);

            GL.Enable(EnableCap.Texture2D);
            GL.BindTexture(TextureTarget.Texture2D, this.textureId_);
            GL.Color4(1, 1, 1, opacity);
            GL.Begin(PrimitiveType.Polygon);
            foreach (Point2 pt in this.points_)
            {
                GL.TexCoord2(pt.txpos.x, pt.txpos.y);
                GL.Vertex2(pt.pos.x, pt.pos.y);
            }
            GL.End();
            GL.Disable(EnableCap.Texture2D);

            this.DrawOutline(Color.Gray);

            GL.Disable(EnableCap.Blend);
        }
        public void DrawOutline(Color c)
        {
            // draw outlines
            GL.Enable(EnableCap.LineSmooth);
            GL.Color3(c.R, c.G, c.B);
            GL.LineWidth(1.0f);
            GL.Begin(PrimitiveType.Lines);
            int l = this.points_.Count;
            for (int i = 0; i < l; ++i)
            {
                Point2 pt = this.points_[i], qt = this.points_[(i + 1) % l];
                GL.Vertex2(pt.pos.x, pt.pos.y);
                GL.Vertex2(qt.pos.x, qt.pos.y);
            }
            GL.End();
        }
    }

    // vanishing points camera calibrator
    public class Camera
    {
        private double xlen, ylen, zlen;
        private MyMatrix3d camera_K = null;
        private MyMatrix3d camera_R = null;
        private MyVector3 camera_t;
        private MyMatrix3d KR = null;
        private MyVector3 Kt;
        private MyVector2 focalPoint;
        private int iSInteriorStructure = 1;		// 1 => interior, -1 => outterior

        public static double zNear = 0.0001;	// default z-near value
        public static double zFar = 20;	// default z-far value
        public MyVector2 h, vx, vy, vz;		// vanishing points
        public Point2 o, x, y, z;			// first 4 corners
        public Point2 w, u, v, r;
        public List<LineSegment2> lineSegments;
        public List<LineSegment2> xycalibrators;
        public bool Calibrated;
        public MyVector3 EyePos;				// camera pos
        public MyVector3 target;				// look at pos
        public MyVector3 upvector = new MyVector3(0, 1, 0);
        public void SwitchCorner()
        {
            this.iSInteriorStructure = this.iSInteriorStructure == 1 ? -1 : 1;
            this.UpdateXYZ();
            this.UpdateWUVRPos();
        }

        private MyVector2 world_origin_imgpos;
        private double xProjScale, yProjScale, zProjScale;
        private int wndViewHeight;

        public void CreateKinectDefault(int kinect_ver = 1)
        {
            this.camera_K = MyMatrix3d.IdentityMatrix();

            if (kinect_ver == 1)
            {
                RGBDBuilder.KINECT_LOCAL_PROJ = RGBDBuilder.KINECT_LOCAL_PROJ_V1;
            }
            else if (kinect_ver == 2) // version 2
            {
                RGBDBuilder.KINECT_LOCAL_PROJ = RGBDBuilder.KINECT_LOCAL_PROJ_V2;
            }
            else if (kinect_ver == 3)
            {
                RGBDBuilder.KINECT_LOCAL_PROJ = RGBDBuilder.KINECT_LOCAL_PROJ_V3;
            }
            this.camera_K[0, 0] = RGBDBuilder.KINECT_LOCAL_PROJ[0];
            this.camera_K[1, 1] = RGBDBuilder.KINECT_LOCAL_PROJ[5];
            this.camera_K[0, 2] = RGBDBuilder.KINECT_LOCAL_PROJ[2];
            this.camera_K[1, 2] = RGBDBuilder.KINECT_LOCAL_PROJ[6];

            this.camera_R = MyMatrix3d.IdentityMatrix();
            this.camera_t = new MyVector3();

            this.KR = this.camera_K * this.camera_R;
            this.Kt = this.camera_K * this.camera_t;
            this.Calibrated = true;
        }

        // interface
        public List<Point2> adjustablePoints;

        // gl matrices
        public double[] glprojMatrix;
        public double[] glmodelViewMatrix;
        public double[] ballMatrix;
        public double[] currObjTransMatrix;
        public int[] viewport;
        private MyMatrix4d objectSpaceTransform = MyMatrix4d.IdentityMatrix();
        public void SetBallMatrix(double[] ballMat)
        {
            this.ballMatrix = ballMat;
        }
        public void SetObjectSpaceTransform(MyMatrix4d T)
        {
            this.objectSpaceTransform = T;
        }
        public MyMatrix4d GetObjectSpaceTransform()
        {
            return this.objectSpaceTransform;
        }

        public bool IsPointInViewport(MyVector2 point)
        {
            return this.viewport[0] < point.x && point.x < this.viewport[2] &&
                this.viewport[1] < point.y && point.y < this.viewport[3];
        }

        public Camera()
        {
            this.CreateCalibratingElements();
        }

        public void Init(MyVector2 vx, MyVector2 vy, int w, int h)
        {
            this.vx = vx;
            this.vy = vy;
            this.h = new MyVector2(w / 2, h / 2);
            this.vz = ComputeTriangleV3(this.h, vx, vy);
            this.o.pos = new MyVector2(w / 2, h / 2);
            this.xlen = 100;
            this.ylen = 100;
            this.zlen = 100;
            this.UpdateXYZ();
            this.UpdateWUVRPos();

            Random rand = new Random();
            MyVector2 u0 = new MyVector2(rand.Next(0, w), rand.Next(0, h));
            this.xycalibrators[0].u.pos = u0;
            this.xycalibrators[0].v.pos = u0 + (this.vx - u0).Normalize() * this.xlen;
            MyVector2 u1 = new MyVector2(rand.Next(0, w), rand.Next(0, h));
            this.xycalibrators[1].u.pos = u1;
            this.xycalibrators[1].v.pos = u1 + (this.vx - u1).Normalize() * this.xlen;
            MyVector2 v0 = new MyVector2(rand.Next(0, w), rand.Next(0, h));
            this.xycalibrators[2].u.pos = v0;
            this.xycalibrators[2].v.pos = v0 + (this.vy - v0).Normalize() * this.ylen;
            MyVector2 v1 = new MyVector2(rand.Next(0, w), rand.Next(0, h));
            this.xycalibrators[3].u.pos = v1;
            this.xycalibrators[3].v.pos = v1 + (this.vy - v1).Normalize() * this.ylen;
            this.UpdateOldPos();
        }
        public void SetViewport(Image<Bgr, Byte> canvas, SketchView view)
        {
            if (canvas != null)
            {
                this.viewport = new int[4];
                this.viewport[0] = 0;
                this.viewport[1] = view.Height - canvas.Height;
                this.viewport[2] = canvas.Width;
                this.viewport[3] = canvas.Height;
            }
        }
        public void UpdateFocalPoint(MyVector2 fp)
        {
            this.h = fp;
        }
        public void UpdateOldPos()
        {
            this.o.UpdateOldPos(); this.x.UpdateOldPos(); this.y.UpdateOldPos(); this.z.UpdateOldPos();
            this.w.UpdateOldPos(); this.u.UpdateOldPos(); this.v.UpdateOldPos(); this.r.UpdateOldPos();
            foreach (LineSegment2 segment in this.xycalibrators)
            {
                segment.u.UpdateOldPos();
                segment.v.UpdateOldPos();
            }
        }
        public void UpdateXYZLengths()
        {
            this.xlen = (this.x.pos - this.o.pos).Length();
            this.ylen = (this.y.pos - this.o.pos).Length();
            this.zlen = (this.z.pos - this.o.pos).Length();
        }
        public void UpdateVp()
        {
            this.vx = ComputeLineIntersectPoint(this.xycalibrators[0].u.pos, this.xycalibrators[0].v.pos,
                this.xycalibrators[1].u.pos, this.xycalibrators[1].v.pos);
            this.vy = ComputeLineIntersectPoint(this.xycalibrators[2].u.pos, this.xycalibrators[2].v.pos,
                this.xycalibrators[3].u.pos, this.xycalibrators[3].v.pos);
            this.vz = ComputeTriangleV3(this.h, this.vx, this.vy);
            this.UpdateXYZ();
            this.UpdateWUVRPos();
        }

        public void UpdateWUVRPos()
        {
            // compute the positions of w, u, v, r based on 
            // o, x, y, z, vx, vy, vz
            u.pos = ComputeLineIntersectPoint(vx, z.pos, vz, x.pos);
            v.pos = ComputeLineIntersectPoint(vy, z.pos, vz, y.pos);
            r.pos = ComputeLineIntersectPoint(vy, x.pos, vx, y.pos);
            w.pos = ComputeLineIntersectPoint(vy, u.pos, vx, v.pos);
        }
        public void MoveRoomPoint(Point2 p, MyVector2 off)
        {
            if (p == this.o)
            {
                this.o.Translate(off);
                this.UpdateXYZ();
            }
            else if (p == this.x)
            {
                MyVector2 ox = (this.vx - this.o.pos).Normalize();
                this.x.pos = this.x.oldpos + off.Dot(ox) * ox;
            }
            else if (p == this.y)
            {
                MyVector2 oy = (this.vy - this.o.pos).Normalize();
                this.y.pos = this.y.oldpos + off.Dot(oy) * oy;
            }
            else if (p == this.z)
            {
                MyVector2 oz = (this.vz - this.o.pos).Normalize();
                this.z.pos = this.z.oldpos + off.Dot(oz) * oz;
            }

            this.UpdateWUVRPos();

        }
        public bool Calibrate(int[] viewport, int wndHeight, double camera_height)
        {
            this.viewport = viewport;
            this.wndViewHeight = wndHeight;

            double w = this.viewport[2], h = this.viewport[3];
            double threshold = 50 * Math.Max(w, h);

            // Check if Vanishing Points lie at infinity
            bool[] infchk = new bool[3];
            infchk[0] = (Math.Abs(vx.x) > threshold || Math.Abs(vx.y) > threshold);
            infchk[1] = (Math.Abs(vy.x) > threshold || Math.Abs(vy.y) > threshold);
            infchk[2] = (Math.Abs(vz.x) > threshold || Math.Abs(vz.y) > threshold);

            int chkcount = 0;
            for (int i = 0; i < 3; ++i)
            {
                if (infchk[i]) chkcount++;
            }

            Console.WriteLine("calibrating with " + (3 - chkcount) + " vanishing points ...");

            double f = 0, u0 = 0, v0 = 0; // focal length, principal point
            //None
            if (chkcount == 0)
            {
                double Mats_11 = vy.x + vx.x;
                double Mats_12 = vx.y + vy.y;
                double Mats_13 = vy.x * vx.x + vx.y * vy.y;
                double Mats_21 = vy.x + vz.x;
                double Mats_22 = vy.y + vz.y;
                double Mats_23 = vy.x * vz.x + vy.y * vz.y;
                double Mats_31 = vz.x + vx.x;
                double Mats_32 = vz.y + vx.y;
                double Mats_33 = vz.x * vx.x + vz.y * vx.y;

                double A_11 = Mats_11 - Mats_21; double A_12 = Mats_12 - Mats_22;
                double A_21 = Mats_11 - Mats_31; double A_22 = Mats_12 - Mats_32;
                double b_1 = Mats_13 - Mats_23; double b_2 = Mats_13 - Mats_33;
                double detA = A_11 * A_22 - A_12 * A_21;
                u0 = (A_22 * b_1 - A_12 * b_2) / detA;
                v0 = (A_11 * b_2 - A_21 * b_1) / detA;

                double temp = Mats_11 * u0 + Mats_12 * v0 - Mats_13 - u0 * u0 - v0 * v0;
                if (temp < 0)
                {
                    Console.WriteLine("Calibration failed: focal length negative!!");
                    return false;
                }
                f = Math.Sqrt(temp);

                // the geometric way
                MyVector2 cp = ComputeTriangleOrthoCenter(vx, vy, vz);
                double ff = Math.Sqrt(-(vz - cp).Dot(vx - cp));
                MyMatrix3d K = new MyMatrix3d();
                K[0, 0] = K[1, 1] = ff;
                K[2, 2] = 1;
                K[0, 2] = cp.x;
                K[1, 2] = cp.y;
                this.focalPoint = cp;

                double S = (new MyVector3(vx - vz, 0).Cross(new MyVector3(vy - vz, 0))).Length();
                double sx = (new MyVector3(vy - cp, 0).Cross(new MyVector3(vz - cp, 0))).Length();
                double sy = (new MyVector3(vx - cp, 0).Cross(new MyVector3(vz - cp, 0))).Length();
                double sz = (new MyVector3(vx - cp, 0).Cross(new MyVector3(vy - cp, 0))).Length();
                double r1 = sx / S;
                double r2 = sy / S;
                double r3 = sz / S;
                MyVector3 q1 = new MyVector3((vx - cp) / f, 1) * r1;
                MyVector3 q2 = new MyVector3((vy - cp) / f, 1) * r2;
                MyVector3 q3 = new MyVector3((vz - cp) / f, 1) * r3;
                q1 = q1.Normalize();
                q2 = q2.Normalize();
                q3 = q3.Normalize();

                MyVector3 q33 = q1.Cross(q2).Normalize();
                if (q33.Dot(q3) < 0)
                    q3 = (new MyVector3() - q33).Normalize();
                else
                    q3 = q33;
                MyVector3 q23 = q3.Cross(q1).Normalize();
                if (q23.Dot(q2) < 0)
                    q2 = (new MyVector3() - q23).Normalize();
                else
                    q2 = q23;


                double e1 = q1.Dot(q2);
                double e2 = q1.Dot(q3);
                double e3 = q2.Dot(q3);
                double ee = e1 + e2 + e3;
                Console.WriteLine("rotational " + ee);

                MyMatrix3d R = new MyMatrix3d(q1, q2, q3);

                this.camera_K = K;
                this.camera_R = R;
            }
            else if (chkcount == 1)
            {
                MyVector2 v1 = vx, v2 = vy;
                if (infchk[0] == true)
                {
                    v1 = vy;
                    v2 = vz;
                }
                else if (infchk[1] == true)
                {
                    v1 = vx;
                    v2 = vz;
                }
                double r = ((w / 2 - v1.x) * (v2.x - v1.x) + (h / 2 - v1.y) * (v2.y - v1.y)) / (Math.Pow(v2.x - v1.x, 2) + Math.Pow(v2.y - v1.y, 2));
                u0 = v1.x + r * (v2.x - v1.x);
                v0 = v1.y + r * (v2.y - v1.y);
                double temp = u0 * (v1.x + v2.x) + v0 * (v2.y + v1.y) - (v1.x * v2.x + v2.y * v1.y + u0 * u0 + v0 * v0);
                if (temp < 0)
                {
                    Console.WriteLine("Calibration failed: focal length negative!!");
                    return false;
                }
                f = Math.Sqrt(temp);
            }
            if (chkcount == 1)
            {
                MyMatrix3d K = new MyMatrix3d();
                K[0, 0] = K[1, 1] = f;
                K[2, 2] = 1;
                K[0, 2] = u0;
                K[1, 2] = v0;
                MyMatrix3d Q = K.Inverse();
                MyVector3 vecx = Q * new MyVector3(vx, 1);
                vecx = vecx.Normalize();
                //	if (vx.x < u0)
                //		vecx = new MyVector3()-vecx;
                MyVector3 vecz = Q * new MyVector3(vz, 1);
                vecz = vecz.Normalize();
                MyVector3 vecy = vecz.Cross(vecx).Normalize();
                this.camera_R = new MyMatrix3d(vecx, vecy, vecz);
                this.camera_K = K;
            }

            if (chkcount == 2)
            {
                u0 = w / 2; v0 = h / 2;
                MyVector3 vecx = new MyVector3();
                MyVector3 vecy = new MyVector3();
                MyVector3 vecz = new MyVector3();
                if (infchk[0])
                {
                    vecx = new MyVector3(vx, 0).Normalize();
                }
                if (infchk[1])
                {
                    vecy = new MyVector3(vy, 0).Normalize();
                    vecy = vecy * -Math.Sign(vecy.y);
                }
                if (infchk[2])
                {
                    vecz = new MyVector3(vz, 0).Normalize();
                }

                if (infchk[0] && infchk[1])
                {
                    u0 = vz.x;
                    v0 = vz.y;
                    vecz = vecx.Cross(vecy);
                }
                else if (infchk[1] && infchk[2])
                {
                    u0 = vx.x;
                    v0 = vx.y;
                    vecx = vecy.Cross(vecz);
                }
                else
                {
                    u0 = vy.x;
                    v0 = vy.y;
                    vecy = vecz.Cross(vecx);
                }

                this.camera_R = new MyMatrix3d(new MyVector3() - vecx, vecy, vecz);
                f = 500;
                MyMatrix3d K = new MyMatrix3d();
                K[0, 0] = K[1, 1] = f;
                K[2, 2] = 1;
                K[0, 2] = u0;
                K[1, 2] = v0;
                this.camera_K = K;
            }

            Console.WriteLine("vanishing z = (" + vz.x + " " + vz.y + ")");
            Console.WriteLine("focal point = (" + camera_K[0, 2] + "," + camera_K[1, 2] + ")");

            this.ComputeCameraT(camera_height, this.o.pos);

            this.Calibrated = true;

            return true;
        }
        public bool CubeCalibrate(int[] viewport, int wndViewHeight, out CalibrateCuboid cube) // two cases
        {
            cube = null;
            this.viewport = viewport;
            this.wndViewHeight = wndViewHeight;
            double w = viewport[2], h = viewport[3];

            // find the projection matrix
            // (u,x,r,y,v,z)
            MyVector2[] imgpts = new MyVector2[6] {
	            this.u.pos,
	            this.x.pos,
	            this.r.pos,
				this.y.pos,
				this.v.pos,
				this.z.pos,
	        };
            MyVector3[] spacepoints = null;
            spacepoints = new MyVector3[6] { // case 1
				new MyVector3(1,-1,1),
				new MyVector3(1,-1,0),
				new MyVector3(1,1,0),
				new MyVector3(-1,1,0),
				new MyVector3(-1,1,1),
				new MyVector3(-1,-1,1)
			};

            if (spacepoints == null) return false;
            double[,] mat = ComputeProjectionMatrixCV(imgpts, spacepoints);

            // get initial guess from mat M$
            MyVector3 m1 = new MyVector3(mat[0, 0], mat[1, 0], mat[2, 0]); // first column
            MyVector3 m2 = new MyVector3(mat[0, 1], mat[1, 1], mat[2, 1]); // second column
            MyVector3 m3 = new MyVector3(mat[0, 2], mat[1, 2], mat[2, 2]); // third column

            // solve directly
            double u = w / 2, v = h / 2;
            double a1 = m1[0], b1 = m1[1], c1 = m1[2];
            double a2 = m2[0], b2 = m2[1], c2 = m2[2];
            double a3 = m3[0], b3 = m3[1], c3 = m3[2];
            MyVector3 b = new MyVector3(-(a1 * a2 + b1 * b2), -(a1 * a3 + b1 * b3), -(a3 * a2 + b3 * b2));
            MyMatrix3d Q = new MyMatrix3d(
                new MyVector3(c1 * c2, c1 * c3, c3 * c2),
                new MyVector3(c1 * a2 + a1 * c2, c1 * a3 + a1 * c3, c3 * a2 + a3 * c2),
                new MyVector3(c1 * b2 + b1 * c2, c1 * b3 + b1 * c3, c3 * b2 + b3 * c2)
            );
            MyVector3 output = Q.Inverse() * b;
            u = -output[1];
            v = -output[2];
            double f2 = output[0] - u * u - v * v;
            if (f2 < 0)
            {
                Console.WriteLine("focal length^2 < 0!");
                //		return false;
            }
            double f = Math.Sqrt(Math.Abs(f2));
            // output error
            double aa = a1 * a2 + b1 * b2 + c1 * c2 * (u * u + v * v + f * f) + (c1 * a2 + a1 * c2) * (-u) + (c1 * b2 + b1 * c2) * (-v);
            double bb = a1 * a3 + b1 * b3 + c1 * c3 * (u * u + v * v + f * f) + (c1 * a3 + a1 * c3) * (-u) + (c1 * b3 + b1 * c3) * (-v);
            double cc = a3 * a2 + b3 * b2 + c3 * c2 * (u * u + v * v + f * f) + (c3 * a2 + a3 * c2) * (-u) + (c3 * b2 + b3 * c2) * (-v);
            double ee = aa * aa + bb * bb + cc * cc;
            Console.WriteLine("- direct solver error: " + ee + "   - R(", false);

            // compute W
            MyMatrix3d W = MyMatrix3d.IdentityMatrix();
            W[2, 2] = u * u + v * v;
            W[0, 2] = W[2, 0] = -u;
            W[1, 2] = W[2, 1] = -v;
            W[2, 2] += f * f;
            W *= (1 / f / f);
            double lambda = Math.Sqrt(m3.Dot(W * m3));
            double l1 = Math.Sqrt(m1.Dot(W * m1)) / lambda;
            double l2 = Math.Sqrt(m2.Dot(W * m2)) / lambda;

            Matrix<double> InvK = new Matrix<double>(3, 3);
            InvK.SetZero();
            InvK[0, 0] = InvK[1, 1] = 1.0 / f;
            InvK[0, 2] = -u / f;
            InvK[1, 2] = -v / f;
            InvK[2, 2] = 1.0;
            Matrix<double> M = new Matrix<double>(mat);
            Matrix<double> F = InvK.Mul(M);
            MyVector3 r1 = new MyVector3(F[0, 0], F[1, 0], F[2, 0]).Normalize();
            MyVector3 r2 = new MyVector3(F[0, 1], F[1, 1], F[2, 1]).Normalize();
            MyVector3 r3 = new MyVector3(F[0, 2], F[1, 2], F[2, 2]).Normalize();
            MyVector3 t = new MyVector3(F[0, 3], F[1, 3], F[2, 3]) * (1.0 / lambda);

            double e1 = r1.Dot(r2);
            double e2 = r1.Dot(r3);
            double e3 = r2.Dot(r3);

            Console.WriteLine(e1.ToString() + " " + e2.ToString() + " " + e3.ToString() + ")", true);
            MyMatrix3d K = new MyMatrix3d();
            K[0, 0] = K[1, 1] = f;
            K[0, 2] = u;
            K[1, 2] = v;
            K[2, 2] = 1.0;
            MyMatrix3d R = new MyMatrix3d(r1, r2, r3);

            this.camera_K = K;
            this.camera_R = R;
            this.camera_t = t;
            this.KR = this.camera_K * this.camera_R;
            this.Kt = this.camera_K * this.camera_t;
            this.EyePos = this.GetEyePosition();

            // get the box
            MyMatrix4d S = MyMatrix4d.ScalingMatrix(l1, l2, 1);
            Point3[] pts3 = new Point3[8];
            for (int i = 0; i < 6; ++i)
            {
                MyVector3 pos = (S * new MyVector4(spacepoints[i], 0)).XYZ();
                pts3[i] = new Point3(pos);
            }

            MyVector3 x3 = (S * new MyVector4(spacepoints[1], 0)).XYZ();
            MyVector3 g3 = (S * new MyVector4(spacepoints[2], 0)).XYZ();
            MyVector3 y3 = (S * new MyVector4(spacepoints[3], 0)).XYZ();
            MyVector3 z3 = (S * new MyVector4(spacepoints[5], 0)).XYZ();
            MyVector3 o3 = x3 + y3 - g3;
            MyVector3 h3 = z3 - o3;
            MyVector3 u3 = x3 + h3;
            MyVector3 v3 = y3 + h3;
            MyVector3 w3 = g3 + h3;

            Point3[] boxpionts = new Point3[8]
			{
				new Point3(o3), new Point3(x3), new Point3(g3), new Point3(y3),
				new Point3(z3), new Point3(u3), new Point3(w3), new Point3(v3)
			};

            cube = new CalibrateCuboid(boxpionts);

            this.Calibrated = true;
            return true;
        }
        public void GetGLMatrices(out double[] glprojmatrix, out double[] glmodelviewmatrix, double w, double h, double znear, double zfar)
        {
            double[,] mat = new double[3, 4];
            for (int i = 0; i < 3; ++i)
            {
                for (int j = 0; j < 3; ++j)
                {
                    mat[i, j] = this.KR[i, j];
                }
            }
            mat[0, 3] = this.Kt[0];
            mat[1, 3] = this.Kt[1];
            mat[2, 3] = this.Kt[2];

            MyMatrix3d LHC = new MyMatrix3d();
            LHC[0, 0] = LHC[1, 1] = LHC[2, 2] = -1;
            LHC[0, 0] = 1;
            double[,] icpara = new double[3, 4];
            double[,] trans = new double[3, 4];
            if (arParamDecompMat(mat, icpara, trans) < 0)
            {
                throw new Exception();
            }
            MyMatrix3d R = new MyMatrix3d();
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    R[i, j] = trans[i, j];
                }
            }
            MyMatrix3d LHCR = LHC * R;
            MyMatrix4d modelViewMatrix = new MyMatrix4d(LHCR);
            modelViewMatrix[0, 3] = trans[0, 3];
            modelViewMatrix[1, 3] = trans[1, 3];
            modelViewMatrix[2, 3] = trans[2, 3];
            modelViewMatrix[1, 3] = modelViewMatrix[1, 3] * (-1);
            modelViewMatrix[2, 3] = modelViewMatrix[2, 3] * (-1);
            modelViewMatrix[3, 3] = 1.0;
            glmodelviewmatrix = modelViewMatrix.Transpose().ToArray();

            MyMatrix4d H_inv = new MyMatrix4d();
            H_inv[0, 0] = 2.0 / w;
            H_inv[0, 2] = -1;
            H_inv[1, 1] = -2.0 / h;
            H_inv[1, 2] = 1.0;
            H_inv[3, 2] = 1.0;
            MyMatrix3d K = new MyMatrix3d();
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    K[i, j] = icpara[i, j] / icpara[2, 2];
                }
            }
            MyMatrix3d y = K * LHC;
            MyMatrix4d y_ = new MyMatrix4d(y);
            MyMatrix4d result = H_inv * (y_);
            double C_ = -(zfar + znear) / (zfar - znear);
            double D_ = -(2 * zfar * znear) / (zfar - znear);
            result[2, 2] = C_;
            result[2, 3] = D_;
            glprojmatrix = result.Transpose().ToArray();

            this.glmodelViewMatrix = glmodelviewmatrix;
            this.glprojMatrix = glprojmatrix;

        }
        public void GetGLMatrices(double w, double h, double znear, double zfar)
        {
            double[,] mat = new double[3, 4];
            for (int i = 0; i < 3; ++i)
            {
                for (int j = 0; j < 3; ++j)
                {
                    mat[i, j] = this.KR[i, j];
                }
            }
            mat[0, 3] = this.Kt[0];
            mat[1, 3] = this.Kt[1];
            mat[2, 3] = this.Kt[2];

            MyMatrix3d LHC = new MyMatrix3d();
            LHC[0, 0] = LHC[1, 1] = LHC[2, 2] = -1;
            LHC[0, 0] = 1;
            double[,] icpara = new double[3, 4];
            double[,] trans = new double[3, 4];
            if (arParamDecompMat(mat, icpara, trans) < 0)
            {
                throw new Exception();
            }
            MyMatrix3d R = new MyMatrix3d();
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    R[i, j] = trans[i, j];
                }
            }
            MyMatrix3d LHCR = LHC * R;
            MyMatrix4d modelViewMatrix = new MyMatrix4d(LHCR);
            modelViewMatrix[0, 3] = trans[0, 3];
            modelViewMatrix[1, 3] = trans[1, 3];
            modelViewMatrix[2, 3] = trans[2, 3];
            modelViewMatrix[1, 3] = modelViewMatrix[1, 3] * (-1);
            modelViewMatrix[2, 3] = modelViewMatrix[2, 3] * (-1);
            modelViewMatrix[3, 3] = 1.0;
            this.glmodelViewMatrix = modelViewMatrix.Transpose().ToArray();

            MyMatrix4d H_inv = new MyMatrix4d();
            H_inv[0, 0] = 2.0 / w;
            H_inv[0, 2] = -1;
            H_inv[1, 1] = -2.0 / h;
            H_inv[1, 2] = 1.0;
            H_inv[3, 2] = 1.0;
            MyMatrix3d K = new MyMatrix3d();
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    K[i, j] = icpara[i, j] / icpara[2, 2];
                }
            }
            MyMatrix3d y = K * LHC;
            MyMatrix4d y_ = new MyMatrix4d(y);
            MyMatrix4d result = H_inv * (y_);
            double C_ = -(zfar + znear) / (zfar - znear);
            double D_ = -(2 * zfar * znear) / (zfar - znear);
            result[2, 2] = C_;
            result[2, 3] = D_;
            this.glprojMatrix = result.Transpose().ToArray();

            // get glu look at parameters
            double[] vmat = this.glmodelViewMatrix;
            MyVector3 forward = new MyVector3(vmat[2], vmat[6], vmat[10]);
            MyVector3 up = new MyVector3(vmat[1], vmat[5], vmat[9]);
            this.target = this.EyePos - forward.Normalize();
            this.upvector = up;

        }

        // draw calibarating element, not good oop
        public void DrawCalibratingElements()
        {
            foreach (LineSegment2 segment in this.lineSegments)
            {
                segment.Draw(2.0f);
            }

            foreach (LineSegment2 segment in this.xycalibrators)
            {
                segment.Draw(2.0f);
            }

            foreach (LineSegment2 segment in this.xycalibrators)
            {
                segment.u.Draw(10.0f);
                segment.v.Draw(10.0f);
            }

            this.o.Draw(10.0f);
            this.x.Draw(7.0f);
            this.y.Draw(7.0f);
            this.z.Draw(7.0f);

            GLBasicDraw.DrawCross2D(this.vx, 4, Color.Red);
            GLBasicDraw.DrawCross2D(this.vy, 4, Color.Green);
        }

        // get camera position and diretion
        public MyVector3 GetEyePosition()
        {
            return new MyVector3() - this.camera_R.Transpose() * this.camera_t;
        }
        public MyVector3 GetDirection()
        {
            // the extrinsic parameters tells how to convert world coordinates into camera coordinates
            // http://ksimek.github.io/2012/08/22/extrinsic/
            // the original viewing direction in world coordinate is (0,0,-1)
            // e.g. R*C+t = 0  is converting camera world coordinates to camera coordinates
            // so to convert a camera coordinate p back to world, one can simply apply
            // R^-1*(p-t), as for vectors R^-1*v
            MyVector3 camera_world_viewing = new MyVector3(0, 0, -1);
            return this.camera_R.Transpose() * camera_world_viewing;
        }

        // determine a line w.r.t. the vanishing points
        public int FindLineVanishingPoint(MyVector2 u, MyVector2 v)
        {
            MyVector2[] vanishing_pts = new MyVector2[3] { this.vx, this.vy, this.vz };
            int t = -1; double min_angle = double.MaxValue;
            for (int i = 0; i < 3; ++i)
            {
                MyVector2 o = vanishing_pts[i];
                double angle = Utils.ComputeTriAngle(u, o, v);
                if (angle < min_angle)
                {
                    min_angle = angle;
                    t = i;
                }
            }
            return t;
        }
        public int FindLineVanishingDir(LineSegment2 line)
        {
            // line (a,b) -- y = ax + b
            MyVector2 line_dir = (line.v.pos - line.u.pos).Normalize();
            MyVector2 line_ctr = (line.v.pos + line.u.pos) / 2;
            MyVector2[] vanishing_dirs = new MyVector2[3] { 
				(this.vx - line_ctr).Normalize(), 
				(this.vy - line_ctr).Normalize(), 
				(this.vz - line_ctr).Normalize() 
			};
            int t = 0; double min_angle = double.MaxValue;
            for (int i = 0; i < 3; ++i)
            {
                MyVector2 o = vanishing_dirs[i];
                double angle = Math.Acos(Math.Abs(o.Dot(line_dir)));
                if (angle < min_angle)
                {
                    min_angle = angle;
                    t = i;
                }
            }
            return t;
        }

        // interfaces
        public CalibrateCuboid ComputeRoomBox3DCoord()
        {
            MyVector3 x3 = this.Compute3DCoordinate(this.x.pos, 0);
            MyVector3 y3 = this.Compute3DCoordinate(this.y.pos, 0);
            MyVector3 r3 = this.Compute3DCoordinate(this.r.pos, 0);
            MyVector3 o3 = new MyVector3();

            MyVector2 oo = this.Compute2D(o3);
            double error = (oo - this.o.pos).Length();
            Console.WriteLine("origin error = " + error);

            MyVector2 xx = this.Compute2D(x3);
            error = (xx - this.x.pos).Length();
            Console.WriteLine("x error = " + error);

            MyVector2 yy = this.Compute2D(y3);
            error = (yy - this.y.pos).Length();
            Console.WriteLine("y error = " + error);

            MyVector2 rr = this.Compute2D(r3);
            error = (rr - this.r.pos).Length();
            Console.WriteLine("r error = " + error);

            MyVector3 z3 = this.Compute3DCoordinate(this.z.pos, this.o.pos);
            MyVector3 h3 = z3 - o3;
            MyVector3 u3 = x3 + h3;
            MyVector3 v3 = y3 + h3;
            MyVector3 w3 = r3 + h3;

            Point3[] boxpionts = new Point3[8]
			{
				new Point3(o3), new Point3(x3), new Point3(r3), new Point3(y3),
				new Point3(z3), new Point3(u3), new Point3(w3), new Point3(v3)
			};

            return new CalibrateCuboid(boxpionts);

        }
        public void Read(string file)
        {
            StreamReader sr = new StreamReader(file);
            char[] delimiters = { ' ', '\t' };
            string s = "";

            while (sr.Peek() > -1)
            {
                s = sr.ReadLine();
                string[] tokens = s.Split(delimiters);

                int j = 0;
                this.o.pos.x = double.Parse(tokens[j++]); this.o.pos.y = double.Parse(tokens[j++]);
                this.x.pos.x = double.Parse(tokens[j++]); this.x.pos.y = double.Parse(tokens[j++]);
                this.y.pos.x = double.Parse(tokens[j++]); this.y.pos.y = double.Parse(tokens[j++]);
                this.z.pos.x = double.Parse(tokens[j++]); this.z.pos.y = double.Parse(tokens[j++]);
                this.w.pos.x = double.Parse(tokens[j++]); this.w.pos.y = double.Parse(tokens[j++]);
                this.u.pos.x = double.Parse(tokens[j++]); this.u.pos.y = double.Parse(tokens[j++]);
                this.v.pos.x = double.Parse(tokens[j++]); this.v.pos.y = double.Parse(tokens[j++]);
                this.r.pos.x = double.Parse(tokens[j++]); this.r.pos.y = double.Parse(tokens[j++]);

                this.xycalibrators[0].u.pos.x = double.Parse(tokens[j++]);
                this.xycalibrators[0].u.pos.y = double.Parse(tokens[j++]);
                this.xycalibrators[0].v.pos.x = double.Parse(tokens[j++]);
                this.xycalibrators[0].v.pos.y = double.Parse(tokens[j++]);
                this.xycalibrators[1].u.pos.x = double.Parse(tokens[j++]);
                this.xycalibrators[1].u.pos.y = double.Parse(tokens[j++]);
                this.xycalibrators[1].v.pos.x = double.Parse(tokens[j++]);
                this.xycalibrators[1].v.pos.y = double.Parse(tokens[j++]);
                this.xycalibrators[2].u.pos.x = double.Parse(tokens[j++]);
                this.xycalibrators[2].u.pos.y = double.Parse(tokens[j++]);
                this.xycalibrators[2].v.pos.x = double.Parse(tokens[j++]);
                this.xycalibrators[2].v.pos.y = double.Parse(tokens[j++]);
                this.xycalibrators[3].u.pos.x = double.Parse(tokens[j++]);
                this.xycalibrators[3].u.pos.y = double.Parse(tokens[j++]);
                this.xycalibrators[3].v.pos.x = double.Parse(tokens[j++]);
                this.xycalibrators[3].v.pos.y = double.Parse(tokens[j++]);

            }
            sr.Close();

            this.UpdateXYZLengths();
            this.UpdateVp();
            this.UpdateOldPos();
            this.UpdateXYZLengths();

        }
        public void Save(string file)
        {
            StreamWriter sw = new StreamWriter(file);

            sw.Write(this.o.pos.x + " " + this.o.pos.y + " ");
            sw.Write(this.x.pos.x + " " + this.x.pos.y + " ");
            sw.Write(this.y.pos.x + " " + this.y.pos.y + " ");
            sw.Write(this.z.pos.x + " " + this.z.pos.y + " ");
            sw.Write(this.w.pos.x + " " + this.w.pos.y + " ");
            sw.Write(this.u.pos.x + " " + this.u.pos.y + " ");
            sw.Write(this.v.pos.x + " " + this.v.pos.y + " ");
            sw.Write(this.r.pos.x + " " + this.r.pos.y + " ");

            sw.Write(this.xycalibrators[0].u.pos.x + " " + this.xycalibrators[0].u.pos.y + " ");
            sw.Write(this.xycalibrators[0].v.pos.x + " " + this.xycalibrators[0].v.pos.y + " ");
            sw.Write(this.xycalibrators[1].u.pos.x + " " + this.xycalibrators[1].u.pos.y + " ");
            sw.Write(this.xycalibrators[1].v.pos.x + " " + this.xycalibrators[1].v.pos.y + " ");
            sw.Write(this.xycalibrators[2].u.pos.x + " " + this.xycalibrators[2].u.pos.y + " ");
            sw.Write(this.xycalibrators[2].v.pos.x + " " + this.xycalibrators[2].v.pos.y + " ");
            sw.Write(this.xycalibrators[3].u.pos.x + " " + this.xycalibrators[3].u.pos.y + " ");
            sw.Write(this.xycalibrators[3].v.pos.x + " " + this.xycalibrators[3].v.pos.y + " ");

            sw.Close();
        }
        public void Serialize(StreamWriter sw)
        {
            sw.WriteLine("camera");
            double[] camera_k = this.camera_K.ToArray();
            for (int i = 0; i < 9; ++i)
            {
                sw.Write(camera_k[i] + " ");
            }
            sw.WriteLine();
            double[] camera_r = this.camera_R.ToArray();
            for (int i = 0; i < 9; ++i)
            {
                sw.Write(camera_r[i] + " ");
            }
            sw.WriteLine();
            double[] camera_t = this.camera_t.ToArray();
            for (int i = 0; i < 3; ++i)
            {
                sw.Write(camera_t[i] + " ");
            }
            sw.WriteLine();
        }
        public void DeSerialize(MyMatrix3d K, MyMatrix3d R, MyVector3 t, int[] viewport, int wndHeight)
        {
            this.camera_K = K;
            this.camera_R = R;
            this.camera_t = t;
            this.KR = this.camera_K * this.camera_R;
            this.Kt = this.camera_K * this.camera_t;
            this.viewport = viewport;
            this.wndViewHeight = wndHeight;
            this.EyePos = new MyVector3() - this.camera_R.Transpose() * this.camera_t;
            this.GetGLMatrices(viewport[2], viewport[3], Camera.zNear, Camera.zFar);
            this.Calibrated = true;
        }

        private MyVector3[] Compute3DCoordinate(List<Point2> impos, double height)
        {
            MyVector2[] pt2 = new MyVector2[4] {
				new MyVector2(0,0),
				new MyVector2(1,0),
				new MyVector2(1,1),
				new MyVector2(0,1)
			};
            MyVector2[] impt2 = new MyVector2[4] {
				this.Compute2D(new MyVector3(pt2[0], height)),
				this.Compute2D(new MyVector3(pt2[1], height)),
				this.Compute2D(new MyVector3(pt2[2], height)),
				this.Compute2D(new MyVector3(pt2[3], height))
			};
            MyMatrix3d H = Utils.ComputeHomographyMatrix(impt2, pt2);
            MyVector3[] outpt3 = new MyVector3[impos.Count];
            int index = 0;
            foreach (Point2 p in impos)
            {
                MyVector3 point = H * new MyVector3(p.pos, 1);
                point.HomogenousNormalize();
                point.z = height;
                outpt3[index++] = point;
            }
            return outpt3;
        }
        private MyVector3 Compute3DCoordinate(MyVector2 impos, double height)
        {
            MyVector2[] pt2 = new MyVector2[4] {
				new MyVector2(0,0),
				new MyVector2(1,0),
				new MyVector2(1,1),
				new MyVector2(0,1)
			};
            MyVector2[] impt2 = new MyVector2[4] {
				this.Compute2D(new MyVector3(pt2[0],height)),
				this.Compute2D(new MyVector3(pt2[1],height)),
				this.Compute2D(new MyVector3(pt2[2],height)),
				this.Compute2D(new MyVector3(pt2[3],height))
			};
            MyMatrix3d H = Utils.ComputeHomographyMatrix(impt2, pt2);
            MyVector3 outpt3 = H * new MyVector3(impos, 1);
            outpt3.HomogenousNormalize();
            outpt3.z = height;
            return outpt3;
        }
        private MyVector3 Compute3DCoordinate(MyVector2 impos, MyVector2 groundproj)
        {
            MyVector3 vanishingline = new MyVector3(this.vx, 1).Cross(new MyVector3(this.vy, 1));
            MyVector3 pointT = new MyVector3(impos, 1);
            MyVector3 pointB = new MyVector3(groundproj, 1);
            MyVector3 pointQ = new MyVector3(groundproj, 1);
            MyVector3 pointV = new MyVector3(this.vz, 1);

            MyVector3 KR_1 = new MyVector3(this.KR[0, 0], this.KR[1, 0], this.KR[2, 0]);
            MyVector3 KR_2 = new MyVector3(this.KR[0, 1], this.KR[1, 1], this.KR[2, 1]);
            MyMatrix3d homographymatrix = new MyMatrix3d(KR_1, KR_2, Kt).Inverse();
            MyVector3 pointS = homographymatrix * pointQ;
            pointS.HomogenousNormalize();
            MyVector3 BT = pointB.Cross(pointT);
            MyVector3 VT = pointV.Cross(pointT);

            MyVector3 origin = new MyVector3(this.world_origin_imgpos, 1);
            int sign = BT.Dot(VT) < 0 ? -1 : 1;
            pointS.z = -sign * origin.Dot(vanishingline) * Math.Abs(BT.Length()) /
                    (this.zProjScale * pointB.Dot(vanishingline) * Math.Abs(VT.Length()));

            return pointS;
        }

        public MyVector2 Compute2D(MyVector3 p)
        {
            // in case if p is view-transformed
            p = (this.objectSpaceTransform * new MyVector4(p, 1)).XYZ();

            // the camera projection matrix
            MyVector3 q = this.KR * p + this.Kt;
            q.HomogenousNormalize();
            return q.ToMyVector2();
        }
        private void Compute2D(Geometry3 geo)
        {
            foreach (PlanarShape3 shape in geo.planarElements)
                this.Compute2D(shape);
            if (geo is CalibrateCuboid)
            {
                ((CalibrateCuboid)geo).GetVanishingPoints();
            }
        }
        private void Compute2D(PlanarShape3 face)
        {
            // if not previously assigned
            if (face.points2 == null)
            {
                face.points2 = new Point2[face.points.Length];
                for (int i = 0; i < face.points.Length; ++i)
                {
                    face.points2[i] = new Point2();
                }
            }

            int index = 0;
            foreach (Point3 pt in face.points)
            {
                MyVector2 pos = this.Project(pt.pos.x, pt.pos.y, pt.pos.z).ToMyVector2();
                face.points2[index].oldpos = face.points2[index].pos = pos;
                index++;
            }
        }
        private void ProjectPointToCanvas32(MyVector2 screenpt, MyVector3 Canvas3center, MyVector3 Canvas3normal,
            out MyVector3 p3, out double r)
        {
            // the out parameter r measures how close the intersecting point is to the near Canvas3'
            // 1. get transformed Canvas3 normal and center (due to view point change)
            MyVector3 c = this.GetObjSpaceTransformedPoint(Canvas3center);
            MyVector3 nor = this.GetObjSpaceTransformedVector(Canvas3normal);
            double[] ss = new double[3];
            double[] tt = new double[3];
            if (this.UnProject(screenpt.x, screenpt.y, -1, ss) == 0 ||
                this.UnProject(screenpt.x, screenpt.y, 1, tt) == 0)
                p3 = new MyVector3(0, 0, 0);
            MyVector3 s = new MyVector3(ss, 0);
            MyVector3 t = new MyVector3(tt, 0);
            r = (c - t).Dot(nor) / ((s - t).Dot(nor));
            p3 = r * s + (1 - r) * t;
        }
        private MyVector3 ProjectPointToCanvas32(MyVector2 screenpt, MyVector3 Canvas3center, MyVector3 Canvas3normal)
        {
            // the out parameter r measures how close the intersecting point is to the near Canvas3'
            // 1. get transformed Canvas3 normal and center (due to view point change)
            MyVector3 c = this.GetObjSpaceTransformedPoint(Canvas3center);
            MyVector3 nor = this.GetObjSpaceTransformedVector(Canvas3normal);
            double[] ss = new double[3];
            double[] tt = new double[3];
            if (this.UnProject(screenpt.x, screenpt.y, -1, ss) == 0 ||
                this.UnProject(screenpt.x, screenpt.y, 1, tt) == 0)
                return new MyVector3(0, 0, 0);
            MyVector3 s = new MyVector3(ss, 0);
            MyVector3 t = new MyVector3(tt, 0);
            double r = (c - t).Dot(nor) / ((s - t).Dot(nor));
            return r * s + (1 - r) * t;
        }
        private MyVector3 GetObjSpaceTransformedPoint(MyVector3 point)
        {
            return (this.objectSpaceTransform * new MyVector4(point, 1)).XYZ();
        }
        private MyVector3 GetObjSpaceTransformedVector(MyVector3 vector)
        {
            return (this.objectSpaceTransform * new MyVector4(vector, 0)).XYZ();
        }
        private void CreateCalibratingElements()
        {
            this.o = new Point2(); this.u = new Point2();
            this.x = new Point2(); this.v = new Point2();
            this.y = new Point2(); this.w = new Point2();
            this.z = new Point2(); this.r = new Point2();

            Color xcolor = Utils.ColorMall[0]; Color ycolor = Utils.ColorMall[1]; Color zcolor = Utils.ColorMall[2];
            LineSegment2 ox = new LineSegment2(o, x); ox.color = xcolor;
            LineSegment2 oy = new LineSegment2(o, y); oy.color = ycolor;
            LineSegment2 oz = new LineSegment2(o, z); oz.color = zcolor;
            LineSegment2 zu = new LineSegment2(z, u); zu.color = xcolor;
            LineSegment2 zv = new LineSegment2(z, v); zv.color = ycolor;
            LineSegment2 rx = new LineSegment2(r, x); rx.color = ycolor;
            LineSegment2 ry = new LineSegment2(r, y); ry.color = xcolor;
            LineSegment2 rw = new LineSegment2(r, w); rw.color = zcolor;
            LineSegment2 wu = new LineSegment2(w, u); wu.color = ycolor;
            LineSegment2 wv = new LineSegment2(w, v); wv.color = xcolor;
            LineSegment2 ux = new LineSegment2(u, x); ux.color = zcolor;
            LineSegment2 uy = new LineSegment2(v, y); uy.color = zcolor;

            this.o.color = xcolor; this.x.color = xcolor; this.y.color = ycolor; this.z.color = zcolor;

            this.lineSegments = new List<LineSegment2>();
            this.lineSegments.Add(ox);
            this.lineSegments.Add(oy);
            this.lineSegments.Add(oz);
            this.lineSegments.Add(zu);
            this.lineSegments.Add(zv);
            this.lineSegments.Add(rx);
            this.lineSegments.Add(ry);
            this.lineSegments.Add(rw);
            this.lineSegments.Add(wu);
            this.lineSegments.Add(wv);
            this.lineSegments.Add(ux);
            this.lineSegments.Add(uy);

            this.xycalibrators = new List<LineSegment2>();
            this.xycalibrators.Add(new LineSegment2(new Point2(), new Point2()));
            this.xycalibrators.Add(new LineSegment2(new Point2(), new Point2()));
            this.xycalibrators.Add(new LineSegment2(new Point2(), new Point2()));
            this.xycalibrators.Add(new LineSegment2(new Point2(), new Point2()));
            this.xycalibrators[0].color = xcolor;
            this.xycalibrators[0].u.color = this.xycalibrators[0].v.color = xcolor;
            this.xycalibrators[1].color = xcolor;
            this.xycalibrators[1].u.color = this.xycalibrators[1].v.color = xcolor;
            this.xycalibrators[2].color = ycolor;
            this.xycalibrators[2].u.color = this.xycalibrators[2].v.color = ycolor;
            this.xycalibrators[3].color = ycolor;
            this.xycalibrators[3].u.color = this.xycalibrators[3].v.color = ycolor;


            this.adjustablePoints = new List<Point2>();
            this.adjustablePoints.Add(this.o);
            this.adjustablePoints.Add(this.x);
            this.adjustablePoints.Add(this.y);
            this.adjustablePoints.Add(this.z);

            this.ballMatrix = MyMatrix4d.IdentityMatrix().ToArray();
            this.currObjTransMatrix = MyMatrix4d.IdentityMatrix().ToArray();

        }
        private void UpdateXYZ()
        {
            // compute x y z position based on position o and vx, vy, vz
            this.x.pos = this.o.pos - this.iSInteriorStructure * (this.vx - this.o.pos).Normalize() * this.xlen;
            this.y.pos = this.o.pos - this.iSInteriorStructure * (this.vy - this.o.pos).Normalize() * this.ylen;

            int sign = Math.Sign((this.vz - this.o.pos).y);
            this.z.pos = this.o.pos - sign * (this.vz - this.o.pos).Normalize() * this.zlen;
        }
        private void ComputeCameraT(double camera_height, MyVector2 world_origin_imgpos)
        {
            // C = -R'*T (a*o=KT) => C = -R'*aK^-1*o, with C.z = cameraheight;
            MyMatrix3d RT = this.camera_R.Transpose() * this.camera_R;
            double err = Math.Abs(RT.Det() - 1);
            Console.WriteLine("--- camera R accuracy = " + err);
            //	if (err > 1e-5) {
            //		throw new Exception();
            //	}

            MyVector3 o = new MyVector3(world_origin_imgpos, 1);
            MyMatrix3d K_inv = this.camera_K.Inverse();
            MyVector3 q = (this.camera_R.Transpose() * K_inv) * o;
            double alpha = -camera_height / q.z;
            this.camera_t = K_inv * o * alpha;

            this.KR = this.camera_K * this.camera_R;
            this.Kt = this.camera_K * this.camera_t;

            MyMatrix3d X = new MyMatrix3d(new MyVector3(this.vx, 1), new MyVector3(this.vy, 1), new MyVector3(this.vz, 1));
            MyMatrix3d I = X.Inverse() * this.KR;
            this.xProjScale = I[0, 0] / alpha;
            this.yProjScale = I[1, 1] / alpha;
            this.zProjScale = I[2, 2] / alpha;

            this.world_origin_imgpos = world_origin_imgpos;
            this.EyePos = new MyVector3() - this.camera_R.Transpose() * this.camera_t;
        }

        public MyVector3 Project(double objx, double objy, double objz)
        {
            double[] modelview = this.glmodelViewMatrix;
            double[] projection = this.glprojMatrix;
            double[] ballmat = this.ballMatrix;

            //Transformation vectors
            double[] tmpb = new double[4];
            //Arcball transform
            tmpb[0] = ballmat[0] * objx + ballmat[4] * objy + ballmat[8] * objz + ballmat[12];  //w is always 1
            tmpb[1] = ballmat[1] * objx + ballmat[5] * objy + ballmat[9] * objz + ballmat[13];
            tmpb[2] = ballmat[2] * objx + ballmat[6] * objy + ballmat[10] * objz + ballmat[14];
            tmpb[3] = ballmat[3] * objx + ballmat[7] * objy + ballmat[11] * objz + ballmat[15];

            double[] fTempo = new double[8];
            //Modelview transform
            fTempo[0] = modelview[0] * tmpb[0] + modelview[4] * tmpb[1] + modelview[8] * tmpb[2] + modelview[12] * tmpb[3];  //w is always 1
            fTempo[1] = modelview[1] * tmpb[0] + modelview[5] * tmpb[1] + modelview[9] * tmpb[2] + modelview[13] * tmpb[3];
            fTempo[2] = modelview[2] * tmpb[0] + modelview[6] * tmpb[1] + modelview[10] * tmpb[2] + modelview[14] * tmpb[3];
            fTempo[3] = modelview[3] * tmpb[0] + modelview[7] * tmpb[1] + modelview[11] * tmpb[2] + modelview[15] * tmpb[3];

            //Projection transform, the final row of projection matrix is always [0 0 -1 0]
            //so we optimize for that.
            fTempo[4] = projection[0] * fTempo[0] + projection[4] * fTempo[1] + projection[8] * fTempo[2] + projection[12] * fTempo[3];
            fTempo[5] = projection[1] * fTempo[0] + projection[5] * fTempo[1] + projection[9] * fTempo[2] + projection[13] * fTempo[3];
            fTempo[6] = projection[2] * fTempo[0] + projection[6] * fTempo[1] + projection[10] * fTempo[2] + projection[14] * fTempo[3];
            fTempo[7] = -fTempo[2];
            //The result normalizes between -1 and 1
            if (fTempo[7] == 0.0)        //The w value
                return new MyVector3();
            fTempo[7] = 1.0 / fTempo[7];
            //Perspective division
            fTempo[4] *= fTempo[7];
            fTempo[5] *= fTempo[7];
            fTempo[6] *= fTempo[7];
            //Window coordinates
            //Map x, y to range 0-1
            MyVector3 windowCoordinate = new MyVector3();
            windowCoordinate[0] = (fTempo[4] * 0.5 + 0.5) * viewport[2] + viewport[0];
            windowCoordinate[1] = (fTempo[5] * 0.5 + 0.5) * viewport[3] + viewport[1];
            //This is only correct when glDepthRange(0.0, 1.0)
            windowCoordinate[2] = (1.0 + fTempo[6]) * 0.5;  //Between 0 and 1

            // convert from gl 2d coords to windows coordinates
            windowCoordinate.y = this.wndViewHeight - windowCoordinate.y;

            return windowCoordinate;
        }
        public int UnProject_mat(double winx, double winy, double winz, double[] objectCoordinate)
        {
            //Transformation matrices
            MyMatrix4d A = new MyMatrix4d(this.glprojMatrix).Transpose() * (new MyMatrix4d(this.glmodelViewMatrix).Transpose()
                * new MyMatrix4d(this.ballMatrix).Transpose());
            MyMatrix4d M = A.Inverse();

            //Transformation of normalized coordinates between -1 and 1
            double[] data_in = new double[4];
            data_in[0] = (winx - (double)this.viewport[0]) / (double)this.viewport[2] * 2.0 - 1.0;
            data_in[1] = (winy - (double)this.viewport[1]) / (double)this.viewport[3] * 2.0 - 1.0;
            data_in[2] = 2.0 * winz - 1.0;
            data_in[3] = 1.0;

            //Objects coordinates
            MyVector4 data_out = M * new MyVector4(data_in, 0);
            if (data_out.w == 0.0)
                return 0;

            double w = 1.0 / data_out.w;
            objectCoordinate[0] = data_out.x * w;
            objectCoordinate[1] = data_out.y * w;
            objectCoordinate[2] = data_out.z * w;

            return 1;
        }
        public int UnProject(double winx, double winy, double winz, double[] objectCoordinate)
        {
            // convert from windows coordinate to opengl coordinate
            winy = this.wndViewHeight - winy;

            //Transformation matrices
            double[] m = new double[16];
            double[] A = new double[16];
            double[] tmpA = new double[16];
            double[] data_in = new double[4];
            double[] data_out = new double[4];
            //Calculation for inverting a matrix, compute projection x modelview
            //and store in A[16]
            MultiplyMatrices4by4OpenGL_FLOAT(tmpA, this.glmodelViewMatrix, this.ballMatrix);
            MultiplyMatrices4by4OpenGL_FLOAT(A, this.glprojMatrix, tmpA);
            //Now compute the inverse of matrix A
            if (glhInvertMatrixf2(A, m) == 0)
                return 0;
            //Transformation of normalized coordinates between -1 and 1
            data_in[0] = (winx - (double)this.viewport[0]) / (double)this.viewport[2] * 2.0 - 1.0;
            data_in[1] = (winy - (double)this.viewport[1]) / (double)this.viewport[3] * 2.0 - 1.0;
            data_in[2] = 2.0 * winz - 1.0;
            data_in[3] = 1.0;
            //Objects coordinates
            MultiplyMatrixByVector4by4OpenGL_FLOAT(data_out, m, data_in);
            if (data_out[3] == 0.0)
                return 0;
            data_out[3] = 1.0 / data_out[3];
            objectCoordinate[0] = data_out[0] * data_out[3];
            objectCoordinate[1] = data_out[1] * data_out[3];
            objectCoordinate[2] = data_out[2] * data_out[3];
            return 1;
        }
        private static void MultiplyMatrices4by4OpenGL_FLOAT(double[] result, double[] matrix1, double[] matrix2)
        {
            result[0] = matrix1[0] * matrix2[0] +
                matrix1[4] * matrix2[1] +
                matrix1[8] * matrix2[2] +
                matrix1[12] * matrix2[3];
            result[4] = matrix1[0] * matrix2[4] +
                matrix1[4] * matrix2[5] +
                matrix1[8] * matrix2[6] +
                matrix1[12] * matrix2[7];
            result[8] = matrix1[0] * matrix2[8] +
                matrix1[4] * matrix2[9] +
                matrix1[8] * matrix2[10] +
                matrix1[12] * matrix2[11];
            result[12] = matrix1[0] * matrix2[12] +
                matrix1[4] * matrix2[13] +
                matrix1[8] * matrix2[14] +
                matrix1[12] * matrix2[15];
            result[1] = matrix1[1] * matrix2[0] +
                matrix1[5] * matrix2[1] +
                matrix1[9] * matrix2[2] +
                matrix1[13] * matrix2[3];
            result[5] = matrix1[1] * matrix2[4] +
                matrix1[5] * matrix2[5] +
                matrix1[9] * matrix2[6] +
                matrix1[13] * matrix2[7];
            result[9] = matrix1[1] * matrix2[8] +
                matrix1[5] * matrix2[9] +
                matrix1[9] * matrix2[10] +
                matrix1[13] * matrix2[11];
            result[13] = matrix1[1] * matrix2[12] +
                matrix1[5] * matrix2[13] +
                matrix1[9] * matrix2[14] +
                matrix1[13] * matrix2[15];
            result[2] = matrix1[2] * matrix2[0] +
                matrix1[6] * matrix2[1] +
                matrix1[10] * matrix2[2] +
                matrix1[14] * matrix2[3];
            result[6] = matrix1[2] * matrix2[4] +
                matrix1[6] * matrix2[5] +
                matrix1[10] * matrix2[6] +
                matrix1[14] * matrix2[7];
            result[10] = matrix1[2] * matrix2[8] +
                matrix1[6] * matrix2[9] +
                matrix1[10] * matrix2[10] +
                matrix1[14] * matrix2[11];
            result[14] = matrix1[2] * matrix2[12] +
                matrix1[6] * matrix2[13] +
                matrix1[10] * matrix2[14] +
                matrix1[14] * matrix2[15];
            result[3] = matrix1[3] * matrix2[0] +
                matrix1[7] * matrix2[1] +
                matrix1[11] * matrix2[2] +
                matrix1[15] * matrix2[3];
            result[7] = matrix1[3] * matrix2[4] +
                matrix1[7] * matrix2[5] +
                matrix1[11] * matrix2[6] +
                matrix1[15] * matrix2[7];
            result[11] = matrix1[3] * matrix2[8] +
                matrix1[7] * matrix2[9] +
                matrix1[11] * matrix2[10] +
                matrix1[15] * matrix2[11];
            result[15] = matrix1[3] * matrix2[12] +
                matrix1[7] * matrix2[13] +
                matrix1[11] * matrix2[14] +
                matrix1[15] * matrix2[15];
        }
        private static void MultiplyMatrixByVector4by4OpenGL_FLOAT(double[] resultvector, double[] matrix, double[] pvector)
        {
            resultvector[0] = matrix[0] * pvector[0] + matrix[4] * pvector[1] + matrix[8] * pvector[2] + matrix[12] * pvector[3];
            resultvector[1] = matrix[1] * pvector[0] + matrix[5] * pvector[1] + matrix[9] * pvector[2] + matrix[13] * pvector[3];
            resultvector[2] = matrix[2] * pvector[0] + matrix[6] * pvector[1] + matrix[10] * pvector[2] + matrix[14] * pvector[3];
            resultvector[3] = matrix[3] * pvector[0] + matrix[7] * pvector[1] + matrix[11] * pvector[2] + matrix[15] * pvector[3];
        }
        private static void SWAP_ROWS_FLOAT(double[] a, double[] b)
        {
            double[] _tmp = a; (a) = (b); (b) = _tmp;
        }
        private static double MAT(double[] m, int r, int c)
        {
            return m[(c) * 4 + (r)];
        }

        //This code comes directly from GLU except that it is for float
        private static int glhInvertMatrixf2(double[] m, double[] out_mat)
        {
            double[][] wtmp = new double[4][];
            for (int i = 0; i < 4; ++i) wtmp[i] = new double[8];

            double m0, m1, m2, m3, s;
            double[] r0 = wtmp[0];
            double[] r1 = wtmp[1];
            double[] r2 = wtmp[2];
            double[] r3 = wtmp[3];
            r0[0] = MAT(m, 0, 0); r0[1] = MAT(m, 0, 1);
            r0[2] = MAT(m, 0, 2); r0[3] = MAT(m, 0, 3);
            r0[4] = 1.0; r0[5] = r0[6] = r0[7] = 0.0;
            r1[0] = MAT(m, 1, 0); r1[1] = MAT(m, 1, 1);
            r1[2] = MAT(m, 1, 2); r1[3] = MAT(m, 1, 3);
            r1[5] = 1.0; r1[4] = r1[6] = r1[7] = 0.0;
            r2[0] = MAT(m, 2, 0); r2[1] = MAT(m, 2, 1);
            r2[2] = MAT(m, 2, 2); r2[3] = MAT(m, 2, 3);
            r2[6] = 1.0; r2[4] = r2[5] = r2[7] = 0.0;
            r3[0] = MAT(m, 3, 0); r3[1] = MAT(m, 3, 1);
            r3[2] = MAT(m, 3, 2); r3[3] = MAT(m, 3, 3);
            r3[7] = 1.0; r3[4] = r3[5] = r3[6] = 0.0;
            /* choose pivot - or die */
            if (Math.Abs(r3[0]) > Math.Abs(r2[0]))
                SWAP_ROWS_FLOAT(r3, r2);
            if (Math.Abs(r2[0]) > Math.Abs(r1[0]))
                SWAP_ROWS_FLOAT(r2, r1);
            if (Math.Abs(r1[0]) > Math.Abs(r0[0]))
                SWAP_ROWS_FLOAT(r1, r0);
            if (0.0 == r0[0])
                return 0;
            /* eliminate first variable     */
            m1 = r1[0] / r0[0];
            m2 = r2[0] / r0[0];
            m3 = r3[0] / r0[0];
            s = r0[1];
            r1[1] -= m1 * s;
            r2[1] -= m2 * s;
            r3[1] -= m3 * s;
            s = r0[2];
            r1[2] -= m1 * s;
            r2[2] -= m2 * s;
            r3[2] -= m3 * s;
            s = r0[3];
            r1[3] -= m1 * s;
            r2[3] -= m2 * s;
            r3[3] -= m3 * s;
            s = r0[4];
            if (s != 0.0)
            {
                r1[4] -= m1 * s;
                r2[4] -= m2 * s;
                r3[4] -= m3 * s;
            }
            s = r0[5];
            if (s != 0.0)
            {
                r1[5] -= m1 * s;
                r2[5] -= m2 * s;
                r3[5] -= m3 * s;
            }
            s = r0[6];
            if (s != 0.0)
            {
                r1[6] -= m1 * s;
                r2[6] -= m2 * s;
                r3[6] -= m3 * s;
            }
            s = r0[7];
            if (s != 0.0)
            {
                r1[7] -= m1 * s;
                r2[7] -= m2 * s;
                r3[7] -= m3 * s;
            }
            /* choose pivot - or die */
            if (Math.Abs(r3[1]) > Math.Abs(r2[1]))
                SWAP_ROWS_FLOAT(r3, r2);
            if (Math.Abs(r2[1]) > Math.Abs(r1[1]))
                SWAP_ROWS_FLOAT(r2, r1);
            if (0.0 == r1[1])
                return 0;
            /* eliminate second variable */
            m2 = r2[1] / r1[1];
            m3 = r3[1] / r1[1];
            r2[2] -= m2 * r1[2];
            r3[2] -= m3 * r1[2];
            r2[3] -= m2 * r1[3];
            r3[3] -= m3 * r1[3];
            s = r1[4];
            if (0.0 != s)
            {
                r2[4] -= m2 * s;
                r3[4] -= m3 * s;
            }
            s = r1[5];
            if (0.0 != s)
            {
                r2[5] -= m2 * s;
                r3[5] -= m3 * s;
            }
            s = r1[6];
            if (0.0 != s)
            {
                r2[6] -= m2 * s;
                r3[6] -= m3 * s;
            }
            s = r1[7];
            if (0.0 != s)
            {
                r2[7] -= m2 * s;
                r3[7] -= m3 * s;
            }
            /* choose pivot - or die */
            if (Math.Abs(r3[2]) > Math.Abs(r2[2]))
                SWAP_ROWS_FLOAT(r3, r2);
            if (0.0 == r2[2])
                return 0;
            /* eliminate third variable */
            m3 = r3[2] / r2[2];
            r3[3] -= m3 * r2[3]; r3[4] -= m3 * r2[4];
            r3[5] -= m3 * r2[5]; r3[6] -= m3 * r2[6]; r3[7] -= m3 * r2[7];
            /* last check */
            if (0.0 == r3[3])
                return 0;
            s = 1.0 / r3[3];             /* now back substitute row 3 */
            r3[4] *= s;
            r3[5] *= s;
            r3[6] *= s;
            r3[7] *= s;
            m2 = r2[3];                  /* now back substitute row 2 */
            s = 1.0 / r2[2];
            r2[4] = s * (r2[4] - r3[4] * m2); r2[5] = s * (r2[5] - r3[5] * m2);
            r2[6] = s * (r2[6] - r3[6] * m2); r2[7] = s * (r2[7] - r3[7] * m2);
            m1 = r1[3];
            r1[4] -= r3[4] * m1; r1[5] -= r3[5] * m1;
            r1[6] -= r3[6] * m1; r1[7] -= r3[7] * m1;
            m0 = r0[3];
            r0[4] -= r3[4] * m0; r0[5] -= r3[5] * m0;
            r0[6] -= r3[6] * m0; r0[7] -= r3[7] * m0;
            m1 = r1[2];                  /* now back substitute row 1 */
            s = 1.0 / r1[1];
            r1[4] = s * (r1[4] - r2[4] * m1); r1[5] = s * (r1[5] - r2[5] * m1);
            r1[6] = s * (r1[6] - r2[6] * m1); r1[7] = s * (r1[7] - r2[7] * m1);
            m0 = r0[2];
            r0[4] -= r2[4] * m0; r0[5] -= r2[5] * m0;
            r0[6] -= r2[6] * m0; r0[7] -= r2[7] * m0;
            m0 = r0[1];                  /* now back substitute row 0 */
            s = 1.0 / r0[0];
            r0[4] = s * (r0[4] - r1[4] * m0); r0[5] = s * (r0[5] - r1[5] * m0);
            r0[6] = s * (r0[6] - r1[6] * m0); r0[7] = s * (r0[7] - r1[7] * m0);

            out_mat[0] = r0[4]; out_mat[4] = r0[5]; out_mat[8] = r0[6]; out_mat[12] = r0[7];
            out_mat[1] = r1[4]; out_mat[5] = r1[5]; out_mat[9] = r1[6]; out_mat[13] = r1[7];
            out_mat[2] = r2[4]; out_mat[6] = r2[5]; out_mat[10] = r2[6]; out_mat[14] = r2[7];
            out_mat[3] = r3[4]; out_mat[7] = r3[5]; out_mat[11] = r3[6]; out_mat[15] = r3[7];

            return 1;
        }
        private static int arParamDecompMat(double[,] source, double[,] cpara, double[,] trans)
        {
            int r, c;
            double[,] Cpara = new double[3, 4];
            double rem1, rem2, rem3;

            if (source[2, 3] >= 0)
            {
                for (r = 0; r < 3; r++)
                {
                    for (c = 0; c < 4; c++)
                    {
                        Cpara[r, c] = source[r, c];
                    }
                }
            }
            else
            {
                for (r = 0; r < 3; r++)
                {
                    for (c = 0; c < 4; c++)
                    {
                        Cpara[r, c] = -(source[r, c]);
                    }
                }
            }

            for (r = 0; r < 3; r++)
            {
                for (c = 0; c < 4; c++)
                {
                    cpara[r, c] = 0.0;
                }
            }
            cpara[2, 2] = norm(Cpara[2, 0], Cpara[2, 1], Cpara[2, 2]);
            trans[2, 0] = Cpara[2, 0] / cpara[2, 2];
            trans[2, 1] = Cpara[2, 1] / cpara[2, 2];
            trans[2, 2] = Cpara[2, 2] / cpara[2, 2];
            trans[2, 3] = Cpara[2, 3] / cpara[2, 2];

            cpara[1, 2] = dot(trans[2, 0], trans[2, 1], trans[2, 2],
                               Cpara[1, 0], Cpara[1, 1], Cpara[1, 2]);
            rem1 = Cpara[1, 0] - cpara[1, 2] * trans[2, 0];
            rem2 = Cpara[1, 1] - cpara[1, 2] * trans[2, 1];
            rem3 = Cpara[1, 2] - cpara[1, 2] * trans[2, 2];
            cpara[1, 1] = norm(rem1, rem2, rem3);
            trans[1, 0] = rem1 / cpara[1, 1];
            trans[1, 1] = rem2 / cpara[1, 1];
            trans[1, 2] = rem3 / cpara[1, 1];

            cpara[0, 2] = dot(trans[2, 0], trans[2, 1], trans[2, 2],
                               Cpara[0, 0], Cpara[0, 1], Cpara[0, 2]);
            cpara[0, 1] = dot(trans[1, 0], trans[1, 1], trans[1, 2],
                               Cpara[0, 0], Cpara[0, 1], Cpara[0, 2]);
            rem1 = Cpara[0, 0] - cpara[0, 1] * trans[1, 0] - cpara[0, 2] * trans[2, 0];
            rem2 = Cpara[0, 1] - cpara[0, 1] * trans[1, 1] - cpara[0, 2] * trans[2, 1];
            rem3 = Cpara[0, 2] - cpara[0, 1] * trans[1, 2] - cpara[0, 2] * trans[2, 2];
            cpara[0, 0] = norm(rem1, rem2, rem3);
            trans[0, 0] = rem1 / cpara[0, 0];
            trans[0, 1] = rem2 / cpara[0, 0];
            trans[0, 2] = rem3 / cpara[0, 0];

            trans[1, 3] = (Cpara[1, 3] - cpara[1, 2] * trans[2, 3]) / cpara[1, 1];
            trans[0, 3] = (Cpara[0, 3] - cpara[0, 1] * trans[1, 3]
                                       - cpara[0, 2] * trans[2, 3]) / cpara[0, 0];

            for (r = 0; r < 3; r++)
            {
                for (c = 0; c < 3; c++)
                {
                    cpara[r, c] /= cpara[2, 2];
                }
            }

            return 0;
        }
        private static double norm(double a, double b, double c)
        {
            return (Math.Sqrt(a * a + b * b + c * c));
        }
        private static double dot(double a1, double a2, double a3, double b1, double b2, double b3)
        {
            return (a1 * b1 + a2 * b2 + a3 * b3);
        }
        private static MyVector2 ComputeTriangleV3(MyVector2 h, MyVector2 v1, MyVector2 v2)
        {
            // this function compute the last point of a triangle, given two points
            // and the ortho center h. The algorithm uses the orthogonalty to
            // solve for the unknow X(x,y) by two linear euqations.
            // (X-H).(V2-V1) = 0 && (X-V1).(V2-H) = 0;
            MyVector2 v2H = v2 - h;
            MyVector2 v21 = v2 - v1;
            Matrix2d A = new Matrix2d(v2H, v21).Transpose();
            MyVector2 b = new MyVector2(v2H.Dot(v1), v21.Dot(h));
            MyVector2 x = A.Inverse() * b;
            return x;
        }
        private static MyVector2 ComputeTriangleOrthoCenter(MyVector2 v1, MyVector2 v2, MyVector2 v3)
        {
            // this function compute the pendicular foot of an triangle, given three points
            // v1, v2 and v3. The algorithm uses the orthogonalty to
            // solve for the unknow H(x,y) by two linear euqations.
            // (H-V1).(V3-V2) = 0 && (H-V3).(V2-V1) = 0;
            MyVector2 v21 = v2 - v1;
            MyVector2 v32 = v3 - v2;
            Matrix2d A = new Matrix2d(v21, v32).Transpose();
            MyVector2 b = new MyVector2(v21.Dot(v3), v32.Dot(v1));
            MyVector2 h = A.Inverse() * b;
            return h;
        }
        private static MyVector2 ComputeLineIntersectPoint(MyVector2 u, MyVector2 v, MyVector2 p, MyVector2 q)
        {
            // compute the intersecting point of two lines: (u,v) and (p,q);
            MyVector3 uu = new MyVector3(u, 1);
            MyVector3 vv = new MyVector3(v, 1);
            MyVector3 pp = new MyVector3(p, 1);
            MyVector3 qq = new MyVector3(q, 1);

            MyVector3 it = uu.Cross(vv).Cross(pp.Cross(qq));
            it.HomogenousNormalize();

            return it.ToMyVector2();
        }
        public static double[,] ComputeProjectionMatrixCV(MyVector2[] imgpts, MyVector3[] spacepts)
        {
            // this function computes the projection matrix given image-to-space points correspondence
            // using DLT algorithm for approximation
            // require: cv eigendecompose, num of points >= 6
            int n = imgpts.Length;
            double[,] mat = new double[n * 2, 12];
            for (int i = 0, j = 0; i < n; ++i, j += 2)
            {
                double x = imgpts[i].x, y = imgpts[i].y;
                int jj = j + 1;
                mat[j, 4] = spacepts[i].x;
                mat[j, 5] = spacepts[i].y;
                mat[j, 6] = spacepts[i].z;
                mat[j, 7] = 1.0;
                mat[j, 8] = -y * spacepts[i].x;
                mat[j, 9] = -y * spacepts[i].y;
                mat[j, 10] = -y * spacepts[i].z;
                mat[j, 11] = -y;
                mat[jj, 0] = spacepts[i].x;
                mat[jj, 1] = spacepts[i].y;
                mat[jj, 2] = spacepts[i].z;
                mat[jj, 3] = 1.0;
                mat[jj, 8] = -x * spacepts[i].x;
                mat[jj, 9] = -x * spacepts[i].y;
                mat[jj, 10] = -x * spacepts[i].z;
                mat[jj, 11] = -x;
            }
            // perform eigen decomposition
            // if n > 6, decompose ATA, else directly decompose A
            Matrix<double> eigvec = new Matrix<double>(12, 12);
            Matrix<double> eigval = new Matrix<double>(12, 1);
            Matrix<double> cvmat = new Matrix<double>(mat);
            double[,] p = new double[3, 4];

            //	if (n > 6)
            {
                cvmat = cvmat.Transpose().Mul(cvmat);
            }
            CvInvoke.Eigen((IInputArray)cvmat, (IOutputArray)eigval, (IOutputArray)eigvec);
            //CvInvoke.cvEigenVV(cvmat.Ptr, eigvec.Ptr, eigval.Ptr, 1e-10, 0, 0);
            for (int i = 0; i < 3; ++i)
            {
                for (int j = 0; j < 4; ++j)
                {
                    p[i, j] = eigvec[11, i * 4 + j];
                }
            }
            // ||p|| = 1 from cv
            return p;
        }
    }
    // reference geometry, used to map the 2d strokes into 3d
    public class Geometry3
    {
        public string categoryName_ = "";
        public RelationGraphNode hostNode_;
        public Mesh bindedMesh = null;
        public Point3[] points;
        public PlanarShape3[] textureFaces = null; // faces for texture binding
        public PlanarShape3[] planarElements;
        public Point3 center;
        public CalibrateCuboid GetBoundingBox()
        {
            List<MyVector3> points = new List<MyVector3>();
            foreach (Point3 p3 in this.points)
            {
                points.Add(p3.pos);
            }
            return Utils.FindBoundingBox(points);
        }
        public void Transform(MyMatrix4d T)
        {
            foreach (Point3 pt in this.points)
            {
                pt.Transform(T);
            }
            this.Compute3DInfo();	// update
        }
        public void Transform_from_origin(MyMatrix4d T)
        {
            foreach (Point3 pt in this.points)
            {
                pt.Transform_from_origin(T);
            }
            this.Compute3DInfo();	// update
        }
        public void Transform_to_origin(MyMatrix4d T)
        {
            foreach (Point3 pt in this.points)
            {
                pt.Transform_to_origin(T);
            }
            this.Compute3DInfo();
        }
        public void Compute3DInfo()
        {
            // compute cuboid center
            this.ComputeCenter();

            // compute quad centers & normals
            this.ComputeElementsCenter();
            this.ComputeElementsNormal();

            // link faces
            this.LinkFaces();
        }
        public bool Has2D()
        {
            return this.planarElements[0].points2 != null;
        }
        private void ComputeCenter()
        {
            // compute center
            MyVector3 c = new MyVector3();
            foreach (Point3 pt in this.points)
            {
                c += pt.pos;
            }
            c /= this.points.Length;
            this.center = new Point3(c);
        }
        private void ComputeElementsCenter()
        {
            foreach (PlanarShape3 e in this.planarElements)
            {
                e.ComputeCenter();
            }
        }
        private void ComputeElementsNormal()
        {
            foreach (PlanarShape3 e in this.planarElements)
            {
                e.ComputeNormal();
            }
        }
        private void LinkFaces()
        {
            foreach (PlanarShape3 face in this.planarElements)
            {
                face.HostGeometry_ = this;
            }
        }

        public void Draw(Color c)
        {
            // draw transparent polygons
            GL.Color3(c.R, c.G, c.B);
            foreach (PlanarShape3 quad in this.planarElements)
            {
                GL.Begin(PrimitiveType.Polygon);
                foreach (Point3 pt3 in quad.points)
                    GL.Vertex3(pt3.pos.x, pt3.pos.y, pt3.pos.z);
                GL.End();
            }
            // draw outlines
            GL.Color3(c.R, c.G, c.B);
            GL.LineWidth(1.0f);
            GL.Begin(PrimitiveType.Lines);
            foreach (PlanarShape3 quad in this.planarElements)
            {
                for (int i = 0; i < quad.points.Length; ++i)
                {
                    Point3 p = quad.points[i], q = quad.points[(i + 1) % quad.points.Length];
                    GL.Vertex3(p.pos.x, p.pos.y, p.pos.z);
                    GL.Vertex3(q.pos.x, q.pos.y, q.pos.z);
                }
            }
            GL.End();
            GL.Disable(EnableCap.Blend);
        }
        public void DrawBlended(Color c, byte opacity)
        {
            GL.Enable(EnableCap.Blend);
            GL.BlendFunc(BlendingFactorSrc.SrcAlpha, BlendingFactorDest.OneMinusSrcAlpha);

            // draw transparent polygons
            GL.Color4(c.R, c.G, c.B, (byte)5);
            foreach (PlanarShape3 quad in this.planarElements)
            {
                GL.Begin(PrimitiveType.Polygon);
                foreach (Point3 pt3 in quad.points)
                    GL.Vertex3(pt3.pos.x, pt3.pos.y, pt3.pos.z);
                GL.End();
            }
            // draw outlines
            GL.Color4(c.R, c.G, c.B, opacity);
            GL.LineWidth(1.0f);
            GL.Begin(PrimitiveType.Lines);
            foreach (PlanarShape3 quad in this.planarElements)
            {
                for (int i = 0; i < quad.points.Length; ++i)
                {
                    Point3 p = quad.points[i], q = quad.points[(i + 1) % quad.points.Length];
                    GL.Vertex3(p.pos.x, p.pos.y, p.pos.z);
                    GL.Vertex3(q.pos.x, q.pos.y, q.pos.z);
                }
            }
            GL.End();
            GL.Disable(EnableCap.Blend);
        }
    }
    // cuboid reference geometry
    public class CalibrateCuboid : Geometry3
    {
        public PlanarShape3[] snappingFaces = new PlanarShape3[2];
        public Point2 corner;
        private MyVector2[] vp; // vanish point -- in windows coordinate
        public CalibrateCuboid(Point3[] pts)
        {
            this.points = pts;
            int l = pts.Length / 2;
            int numquads = l + 2; // top bottom quads
            this.planarElements = new Quad3D[numquads];
            for (int i = 0; i < numquads; ++i) this.planarElements[i] = new Quad3D();
            for (int i = 0; i < l; ++i)
            {
                this.planarElements[i].points = new Point3[4] {
					this.points[i], this.points[(i+1)%l], 
					this.points[(i+1)%l+l], this.points[i+l]
				};
            }
            this.planarElements[l].points = new Point3[l];
            this.planarElements[l + 1].points = new Point3[l];
            for (int i = 0; i < l; ++i)
            {
                this.planarElements[l].points[i] = this.points[l - 1 - i];
                //this.planarElements[l].points[i] = this.points[i];				// bot polygon
                this.planarElements[l + 1].points[i] = this.points[i + l];		// top polygon
                // inverse normal
                //this.planarElements[l + 1].points[i] = this.points[2 * l - 1 - i];		// top polygon
            }
            this.Compute3DInfo();
            this.SetTextureFaces();
        }
        public PlanarShape3 GetBottomFace()
        {
            return this.planarElements[4];
        }
        public PlanarShape3 GetOppositeFace(PlanarShape3 face)
        {
            foreach (PlanarShape3 f in this.planarElements)
            {
                if (face == f) continue;
                if (Math.Abs(f.normal.Dot(face.normal)) > 0.8)
                    return f;
            }
            return null;
        }
        public MyVector3 GetScale()
        {
            double xlen = (this.planarElements[2].center - this.planarElements[0].center).Length();
            double ylen = (this.planarElements[3].center - this.planarElements[1].center).Length();
            double zlen = (this.planarElements[5].center - this.planarElements[4].center).Length();	// top bot
            return new MyVector3(xlen, ylen, zlen);
        }
        public PlanarShape3[] GetYPairFaces()
        {
            return new PlanarShape3[2] {
				this.planarElements[0], this.planarElements[2]
			};
        }
        public PlanarShape3[] GetXPairFaces()
        {
            return new PlanarShape3[2] {
				this.planarElements[3], this.planarElements[1]
			};
        }
        public PlanarShape3[] GetZPairFaces()
        {
            return new PlanarShape3[2] {
				this.planarElements[4], this.planarElements[5]
			};
        }
        public MyVector3[] GetXYZAxes()
        {
            PlanarShape3[] xfaces = this.GetXPairFaces();
            PlanarShape3[] yfaces = this.GetYPairFaces();
            PlanarShape3[] zfaces = this.GetZPairFaces();
            MyVector3[] dirs = new MyVector3[3] {
				(xfaces[1].center - xfaces[0].center).Normalize(),
				(yfaces[1].center - yfaces[0].center).Normalize(),
				(zfaces[1].center - zfaces[0].center).Normalize()
			};
            return dirs;
        }
        public double GetLongestOrShortestEdge(bool isLong)
        {
            PlanarShape3[] xfaces = this.GetXPairFaces();
            PlanarShape3[] yfaces = this.GetYPairFaces();
            PlanarShape3[] zfaces = this.GetZPairFaces();
            double[] lengths = new double[3] {
				(xfaces[1].center - xfaces[0].center).Length(),
				(yfaces[1].center - yfaces[0].center).Length(),
				(zfaces[1].center - zfaces[0].center).Length()
			};
            double longest = lengths[0] > lengths[1] ? lengths[0] : lengths[1];
            longest = longest > lengths[2] ? longest : lengths[2];
            double shortest = lengths[0] < lengths[1] ? lengths[0] : lengths[1];
            shortest = shortest < lengths[2] ? shortest : lengths[2];
            if (isLong)
                return longest;
            else
                return shortest;
        }

        public void GetVanishingPoints()
        {
            // compute vanishing points vx, vy, vz for the cuboid
            PlanarShape3[] xfaces = this.GetXPairFaces();
            PlanarShape3[] yfaces = this.GetYPairFaces();
            PlanarShape3[] zfaces = this.GetZPairFaces();
            MyVector3[] dirs = new MyVector3[3] {
				(xfaces[1].center - xfaces[0].center).Normalize(),
				(yfaces[1].center - yfaces[0].center).Normalize(),
				(zfaces[1].center - zfaces[0].center).Normalize()
			};
            List<MyVector2>[] xyz_lines = new List<MyVector2>[3] {
				new List<MyVector2>(),
				new List<MyVector2>(),
				new List<MyVector2>()
			};
            foreach (PlanarShape3 planar in this.planarElements)
            {
                for (int i = 0; i < 4; ++i)
                {
                    MyVector3 u3 = planar.points[i].pos, v3 = planar.points[(i + 1) % 4].pos;
                    MyVector2 u2 = planar.points2[i].pos, v2 = planar.points2[(i + 1) % 4].pos;
                    MyVector3 uv = (v3 - u3).Normalize();
                    for (int j = 0; j < 3; ++j)
                    {
                        if (Math.Abs(uv.Dot(dirs[j])) > 0.5)
                        {
                            xyz_lines[j].Add(u2);
                            xyz_lines[j].Add(v2);
                            break;
                        }
                    }
                }
            }
            this.vp = new MyVector2[3];
            for (int i = 0; i < 3; ++i)
            {
                List<MyVector2> lines = xyz_lines[i];
                int N = lines.Count / 2;	// N lines
                MyVector2 mean = new MyVector2();
                for (int j = 0, k = 0; j < N; ++j, k += 2)
                {
                    MyVector2 u = lines[k], v = lines[k + 1];
                    MyVector2 p = lines[(k + 2) % lines.Count], q = lines[(k + 3) % lines.Count];
                    MyVector2 o = Utils.GetIntersectPoint(u, v, p, q);
                    mean += o;
                }
                this.vp[i] = mean / N;
            }

            //Console.WriteLine("cuboid vp: vx(" + 
            //	this.vp[0].x.ToString("f2") + ", " + this.vp[0].y.ToString("f2") + "), vy(" +
            //	this.vp[1].x.ToString("f2") + ", " + this.vp[1].y.ToString("f2") + "), vz(" +
            //	this.vp[2].x.ToString("f2") + ", " + this.vp[2].y.ToString("f2") + ")"
            //	); 
        }
        public MyVector2[] GetVP()
        {
            this.GetVanishingPoints();
            return this.vp;
        }
        public int FindVanishingDir(MyVector2 u, MyVector2 v)
        {
            // line (a,b) -- y = ax + b
            MyVector2 uv = (v - u).Normalize();
            MyVector2 c = (v + u) / 2;
            MyVector2[] vanishing_dirs = new MyVector2[3] { 
				(this.vp[0] - c).Normalize(), 
				(this.vp[1] - c).Normalize(), 
				(this.vp[2] - c).Normalize() 
			};
            int t = 0; double min_angle = double.MaxValue;
            for (int i = 0; i < 3; ++i)
            {
                MyVector2 o = vanishing_dirs[i];
                double angle = Math.Acos(Math.Abs(o.Dot(uv)));
                if (angle < min_angle)
                {
                    min_angle = angle;
                    t = i;
                }
            }
            return t;
        }
        public double FindVanishingAngleDistance(LineSegment2 line, out int vindex)
        {
            // line (a,b) -- y = ax + b
            MyVector2 line_dir = (line.v.pos - line.u.pos).Normalize();
            MyVector2 line_ctr = (line.v.pos + line.u.pos) / 2;
            MyVector2[] vanishing_dirs = new MyVector2[3] { 
				(this.vp[0] - line_ctr).Normalize(), 
				(this.vp[1] - line_ctr).Normalize(), 
				(this.vp[2] - line_ctr).Normalize() 
			};
            int t = 0; double min_angle = double.MaxValue;
            for (int i = 0; i < 3; ++i)
            {
                MyVector2 o = vanishing_dirs[i];
                double angle = Math.Acos(Math.Abs(o.Dot(line_dir)));
                if (angle < min_angle)
                {
                    min_angle = angle;
                    t = i;
                }
            }
            vindex = t;
            return min_angle;
        }
        public double[] FindVanishingAngleDistance(LineSegment2 line)
        {
            // line (a,b) -- y = ax + b
            if (line == null) return null;
            MyVector2 line_dir = (line.v.pos - line.u.pos).Normalize();
            MyVector2 line_ctr = (line.v.pos + line.u.pos) / 2;
            MyVector2[] vanishing_dirs = new MyVector2[3] { 
				(this.vp[0] - line_ctr).Normalize(), 
				(this.vp[1] - line_ctr).Normalize(), 
				(this.vp[2] - line_ctr).Normalize() 
			};
            double[] anglesDists = new double[3];
            for (int i = 0; i < 3; ++i)
            {
                MyVector2 o = vanishing_dirs[i];
                anglesDists[i] = Math.Abs(o.Dot(line_dir));
            }
            return anglesDists;
        }
        // how well a line agrees with each vanishing point
        public double[] FindVanishingAngleConfidence(MyVector2 u, MyVector2 v)
        {
            // line (a,b) -- y = ax + b
            MyVector2 uv = (v - u).Normalize();
            MyVector2 c = (v + u) / 2;
            MyVector2[] vanishing_dirs = new MyVector2[3] { 
				(this.vp[0] - c).Normalize(), 
				(this.vp[1] - c).Normalize(), 
				(this.vp[2] - c).Normalize() 
			};
            double[] angles = new double[3];
            for (int i = 0; i < 3; ++i)
            {
                MyVector2 o = vanishing_dirs[i];
                angles[i] = Math.Acos(Math.Abs(o.Dot(uv))); // *180 / Math.PI;
            }
            return angles;
        }
        // find candidate vanishing directions
        public void FindConfidentVanishingDirs(MyVector2 u, MyVector2 v, out List<int> candidate_dirs, out List<double> confident_scores)
        {
            double[] coefs = this.FindVanishingAngleConfidence(u, v);

            Console.WriteLine("(u: " + u.x + ", " + u.y + ")，  (v: " + v.x + ", " + v.y + ")");
            for (int j = 0; j < coefs.Length; ++j)
            {
                Console.WriteLine("(" + j + " angle: " + coefs[j] + ")");
            }

            double min_coef = Math.Min(Math.Min(coefs[0], coefs[1]), coefs[2]);
            candidate_dirs = new List<int>();
            confident_scores = new List<double>();
            for (int i = 0; i < 3; ++i)
            {
                double val = Math.Abs((coefs[i] - min_coef)) / (coefs[i] + min_coef);
                if (val < 0.15)
                {
                    candidate_dirs.Add(i);
                    confident_scores.Add(coefs[i]);
                }
            }
        }
        private void SetTextureFaces()
        {
            this.textureFaces = new PlanarShape3[3];
            this.textureFaces[0] = this.planarElements[0];
            this.textureFaces[1] = this.planarElements[3];
            this.textureFaces[2] = this.planarElements[4];
        }
        public static CalibrateCuboid UnitCube()
        {
            double r = 1.0 / 2;

            Point3[] points = new Point3[8] {
				new Point3(-r,-r,-r), new Point3(r, -r,-r), new Point3(r, r,-r), new Point3(-r, r,-r),
				new Point3(-r,-r, r), new Point3(r, -r, r), new Point3(r, r, r), new Point3(-r, r, r)
			};

            CalibrateCuboid cube = new CalibrateCuboid(points);
            cube.planarElements[0] = new Quad3D(new Point3[4] { points[0], points[1], points[2], points[3] });
            cube.planarElements[1] = new Quad3D(new Point3[4] { points[4], points[5], points[6], points[7] });
            cube.planarElements[2] = new Quad3D(new Point3[4] { points[0], points[1], points[5], points[4] });
            cube.planarElements[3] = new Quad3D(new Point3[4] { points[1], points[2], points[6], points[5] });
            cube.planarElements[4] = new Quad3D(new Point3[4] { points[2], points[3], points[7], points[6] });
            cube.planarElements[5] = new Quad3D(new Point3[4] { points[3], points[0], points[4], points[7] });
            cube.Compute3DInfo();

            return cube;
        }
        public void DrawVanishingPoints()
        {
            GLBasicDraw.DrawCross2D(this.vp[0], 7, Color.Red);
            GLBasicDraw.DrawCross2D(this.vp[1], 7, Color.Green);
            GLBasicDraw.DrawCross2D(this.vp[2], 7, Color.Blue);
        }
    }
    // folding plance reference geometry
    public class FoldingPlanar : Geometry3
    {
        public FoldingPlanar(Point3[] pts)
        {
            this.points = pts;
            int l = pts.Length / 2;
            int numquads = l - 1;
            this.planarElements = new Quad3D[numquads];
            for (int i = 0; i < numquads; ++i)
            {
                this.planarElements[i] = new Quad3D();
                this.planarElements[i].points = new Point3[4] {
					this.points[i], this.points[i+1], 
					this.points[i+l+1], this.points[i+l]
				};
            }
            this.textureFaces = new PlanarShape3[numquads];
            this.Compute3DInfo();
        }
        public FoldingPlanar(Point3[] points, Quad3D[] quads)
        {
            this.points = points;
            this.planarElements = quads;
            this.textureFaces = new PlanarShape3[this.planarElements.Length];
            this.Compute3DInfo();
        }
    }
    // Cylindical ref geometry (this could include cuboid, cylinder)
    public class Prism : Geometry3
    {
        public Prism(Point3[] pts)
        {
            this.points = pts;
            int l = pts.Length / 2;
            int numpolygons = l + 2; // top bottom polygons
            this.planarElements = new Polygon3D[numpolygons];
            for (int i = 0; i < numpolygons; ++i) this.planarElements[i] = new Polygon3D();
            for (int i = 0; i < l; ++i)
            {
                this.planarElements[i].points = new Point3[4] {
					this.points[i], this.points[(i+1)%l], 
					this.points[(i+1)%l+l], this.points[i+l]
				};
            }
            this.planarElements[l].points = new Point3[l];
            this.planarElements[l + 1].points = new Point3[l];
            for (int i = 0; i < l; ++i)
            {
                this.planarElements[l].points[i] = this.points[i];				// top polygon
                this.planarElements[l + 1].points[i] = this.points[i + l];		// bot polygon
            }
            this.textureFaces = new PlanarShape3[numpolygons];
            this.Compute3DInfo();
        }
        public Prism(Point3[] points, Polygon3D[] polygons)
        {
            this.points = points;
            this.planarElements = polygons;
            this.textureFaces = new PlanarShape3[polygons.Length];
            this.Compute3DInfo();
        }
        public static Prism UnitCylinder()
        {
            double r = 1.0 / 2;
            int slices = 20;
            double dtheta = Math.PI * 2 / slices;
            double z = r / 2;
            List<Point3> cylinderPointsTop = new List<Point3>();
            List<Point3> cylinderPointsBot = new List<Point3>();
            for (int i = 0; i < slices; ++i)
            {
                double theta = dtheta * i;
                double x = r * Math.Cos(theta);
                double y = r * Math.Sin(theta);
                cylinderPointsBot.Add(new Point3(new MyVector3(x, y, -z)));
                cylinderPointsTop.Add(new Point3(new MyVector3(x, y, z)));
            }
            cylinderPointsBot.AddRange(cylinderPointsTop);

            return new Prism(cylinderPointsBot.ToArray());

        }
    }
    // relation graph
    public enum EnumGeometricRelation { Contact, Support, Attach, None };
    public class RelationGraphNode
    {
        public int index_;
        public Geometry3 geo_;
        public List<RelationGraphEdge> edges_ = new List<RelationGraphEdge>();
        public RelationGraphNode(Geometry3 gnode)
        {
            gnode.hostNode_ = this;
            this.geo_ = gnode;
        }
        public RelationGraphNode(Geometry3 gnode, int index)
        {
            gnode.hostNode_ = this;
            this.index_ = index;
            this.geo_ = gnode;
        }
        public List<RelationGraphNode> GetSupportedNodes()
        {
            List<RelationGraphNode> supportedNodes = new List<RelationGraphNode>();
            Queue<RelationGraphNode> Q = new Queue<RelationGraphNode>();
            Q.Enqueue(this);
            while (Q.Count > 0)
            {
                RelationGraphNode node = Q.Dequeue();
                foreach (RelationGraphEdge e in node.edges_)
                {
                    // if supporting, node2 is the base
                    if (e.relation_ == EnumGeometricRelation.Support && e.geoNode2_ == node)
                    {
                        supportedNodes.Add(e.geoNode1_);
                        Q.Enqueue(e.geoNode1_);
                    }
                }
            }
            return supportedNodes;
        }
    }
    public class RelationGraphEdge
    {
        public RelationGraphNode geoNode1_, geoNode2_;
        public EnumGeometricRelation relation_;
        public RelationGraphEdge(RelationGraphNode node1, RelationGraphNode node2, EnumGeometricRelation r)
        {
            // if supporting, node2 is the base
            this.geoNode1_ = node1;
            this.geoNode2_ = node2;
            this.geoNode1_.edges_.Add(this);
            this.geoNode2_.edges_.Add(this);
            this.relation_ = r;
        }
    }
    public class RelationGraph
    {
        public List<RelationGraphNode> nodes_ = new List<RelationGraphNode>();
        public List<RelationGraphEdge> edges_ = new List<RelationGraphEdge>();
        public RelationGraph()
        {

        }
        public void AddNode(RelationGraphNode node)
        {
            this.nodes_.Add(node);
        }
        public void AddEdge(RelationGraphEdge edge)
        {
            this.edges_.Add(edge);
        }
        public void RemoveNode(RelationGraphNode node)
        {
            this.nodes_.Remove(node);

            List<RelationGraphEdge> edgesInvolved = new List<RelationGraphEdge>();
            foreach (RelationGraphEdge e in this.edges_)
            {
                if (e.geoNode1_ == node || e.geoNode2_ == node)
                {
                    edgesInvolved.Add(e);
                }
            }
            foreach (RelationGraphEdge e in edgesInvolved)
            {
                this.edges_.Remove(e);
            }
        }
    }
}
