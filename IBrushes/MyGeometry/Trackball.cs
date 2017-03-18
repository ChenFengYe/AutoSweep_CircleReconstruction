using System;

namespace MyGeometry
{
    public class Trackball
    {
        public enum MotionType { None, Rotation, Pan, Scale }

        private MotionType type = MotionType.None;
        private MyVector2 stPt, edPt;
        private MyVector3 stVec;
        private MyVector3 edVec;
        private double radius = 2.0;
        private MyVector4 quat;
        private double w, h;
        private double adjustWidth;
        private double adjustHeight;
        private double angle;
        private MyMatrix4d mat = MyMatrix4d.IdentityMatrix();

        public Trackball(double w, double h)
        {
            SetBounds(w, h);
        }

        public void SetBounds(double w, double h)
        {
            double b = (w < h) ? w : h;
            this.w = w / 2.0;
            this.h = h / 2.0;
            this.adjustWidth = 1.0 / ((w - 1.0) * 0.5);
            this.adjustHeight = 1.0 / ((h - 1.0) * 0.5);
        }
        public void SetRadius(double r)
        {
            this.radius = r;
        }
        public void Click(MyVector2 pt, MotionType type)
        {
            this.stPt = pt;
            this.edPt = pt;
            this.stVec = MapToSphere(pt);
            // this.edVec = this.stVec;
            this.type = type;
        }
        public void Wheel(double scale, MotionType type)
        {
            this.type = type;
            this.adjustHeight *= scale;
            this.adjustWidth += scale;
        }
        public void Drag(MyVector2 pt)
        {
            edPt = pt;
            edVec = MapToSphere(pt);
            //angle = Math.Acos(stVec.Dot(edVec));

            double epsilon = 1.0e-5;

            MyVector3 prep = stVec.Cross(edVec);
            if (prep.Length() > epsilon)
            {
                quat = new MyVector4();

                quat.x = prep.x;
                quat.y = prep.y;
                quat.z = prep.z;
                quat.w = stVec.x + edVec.x + stVec.y + edVec.y + stVec.z + edVec.z;
            }
            else
                quat = new MyVector4();
            //if (prep.Length() > epsilon)
            //{
            //    quat = new MyVector3();

            //    quat.x = prep.x;
            //    quat.y = prep.y;
            //    quat.z = prep.z;
            //}
            //else
            //    quat = new MyVector3();

            //angle = Math.PI / 2.0 * prep.Length();
        }
        public void End()
        {
            quat = new MyVector4();
            type = MotionType.None;
        }
        public MyMatrix4d GetMatrix()
        {
            if (type == MotionType.Rotation)
                return QuatToMyMatrix4d(quat);

            if (type == MotionType.Scale)
            {
                MyMatrix4d m = MyMatrix4d.IdentityMatrix();
                m[0, 0] = m[1, 1] = m[2, 2] = 1.0 + (edPt.x - stPt.x) * adjustWidth;
                //m[0,0] = m[1,1] = m[2,2] = 1.0 +  adjustWidth;
                return m;
            }

            if (type == MotionType.Pan)
            {
                MyMatrix4d m = MyMatrix4d.IdentityMatrix();
                if (edPt == stPt)
                {
                    return m;
                }
                m[0, 3] = 0.07 * (edPt.x - stPt.x);
                m[1, 3] = 0.07 * (edPt.y - stPt.y);
                return m;
            }

            return MyMatrix4d.IdentityMatrix();
        }

        public double GetScale()
        {
            if (type == MotionType.Scale)
                return 1.0 + (edPt.x - stPt.x) * adjustWidth;
            else
                return 1.0;
        }

        private MyVector3 MapToSphere(MyVector2 pt)
        {
            #region old map
            //MyVector2 v = new MyVector2();
            ////v.x = (w - pt.x) * adjustWidth;
            //v.x = (pt.x - this.w) * adjustWidth;
            //v.y = (this.h - pt.y) * adjustHeight;

            //double lenSq = v.Dot(v);
            //MyVector3 v3 = new MyVector3(v.x, 0, v.y);
            //if (lenSq > 1.0)
            //{
            //    double norm = 1.0 / Math.Sqrt(lenSq);
            //    return new MyVector3(v.x * norm, -v.y * norm, 0);
            //}
            //else
            //{
            //    return new MyVector3(v.x, Math.Sqrt(1.0 - lenSq), v.y);
            //}  
            #endregion

            //just. 1116
            MyVector2 v = new MyVector2();
            v.x = (this.w - pt.x) * adjustWidth / radius;
            v.y = (this.h - pt.y) * adjustHeight / radius;

            double lenSq = v.Dot(v);
            //MyVector3 v3 = new MyVector3(v.x, 0, v.y);
            if (lenSq > 1.0)
            {
                double norm = 1.0 / Math.Sqrt(lenSq);
                return new MyVector3(v.x * norm, v.y * norm, 0);
            }
            else
            {
                return new MyVector3(v.x, v.y, Math.Sqrt(1.0 - lenSq));
            }
        }
        public static MyMatrix3d QuatToMyMatrix3d(MyVector4 q)
        {
            double n = q.Dot(q);
            double s = (n > 0.0) ? (2.0 / n) : 0.0f;

            double xs, ys, zs;
            double wx, wy, wz;
            double xx, xy, xz;
            double yy, yz, zz;
            xs = q.x * s; ys = q.y * s; zs = q.z * s;
            wx = q.w * xs; wy = q.w * ys; wz = q.w * zs;
            xx = q.x * xs; xy = q.x * ys; xz = q.x * zs;
            yy = q.y * ys; yz = q.y * zs; zz = q.z * zs;

            MyMatrix3d m = new MyMatrix3d();
            m[0, 0] = 1.0 - (yy + zz); m[1, 0] = xy - wz; m[2, 0] = xz + wy;
            m[0, 1] = xy + wz; m[1, 1] = 1.0 - (xx + zz); m[2, 1] = yz - wx;
            m[0, 2] = xz - wy; m[1, 2] = yz + wx; m[2, 2] = 1.0 - (xx + yy);
            return m;
        }
        private MyMatrix4d QuatToMyMatrix4d(MyVector4 q)
        {
            double n = q.Dot(q);
            double s = (n > 0.0) ? (2.0 / n) : 0.0f;

            double xs, ys, zs;
            double wx, wy, wz;
            double xx, xy, xz;
            double yy, yz, zz;
            xs = q.x * s; ys = q.y * s; zs = q.z * s;
            wx = q.w * xs; wy = q.w * ys; wz = q.w * zs;
            xx = q.x * xs; xy = q.x * ys; xz = q.x * zs;
            yy = q.y * ys; yz = q.y * zs; zz = q.z * zs;

            MyMatrix4d m = new MyMatrix4d();
            m[0, 0] = 1.0 - (yy + zz); m[1, 0] = xy + wz; m[2, 0] = xz - wy;
            m[0, 1] = xy - wz; m[1, 1] = 1.0 - (xx + zz); m[2, 1] = yz + wx;
            m[0, 2] = xz + wy; m[1, 2] = yz - wx; m[2, 2] = 1.0 - (xx + yy);
            m[3, 3] = 1.0;
            return m;
        }
        //just. 1116
        private MyMatrix4d QuatToMyMatrix4d(MyVector3 q, double a)
        {
            q.Normalize();
            double cos_theta = Math.Cos(a);
            double sin_theta = Math.Sin(a);
            double xx = q.x * q.x;
            double yy = q.y * q.y;
            double zz = q.z * q.z;
            double xy = q.x * q.y;
            double xz = q.x * q.z;
            double yz = q.y * q.z;

            MyMatrix4d m = new MyMatrix4d();
            m[0, 0] = cos_theta + xx * (1 - cos_theta);
            m[0, 1] = xy * (1 - cos_theta) - q.z * sin_theta;
            m[0, 2] = xz * (1 - cos_theta) + q.y * sin_theta;

            m[1, 0] = xy * (1 - cos_theta) + q.z * sin_theta;
            m[1, 1] = cos_theta + yy * (1 - cos_theta);
            m[1, 2] = yz * (1 - cos_theta) - q.x * sin_theta;

            m[2, 0] = xz * (1 - cos_theta) - q.y * sin_theta;
            m[2, 1] = yz * (1 - cos_theta) + q.x * sin_theta;
            m[2, 2] = cos_theta + zz * (1 - cos_theta);
            m[3, 3] = 1.0;
            //m.Inverse();
            return m;
        }
    }
}
