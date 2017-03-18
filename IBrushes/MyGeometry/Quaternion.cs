using System;
using System.Collections.Generic;
using System.Text;

namespace MyGeometry
{
	public struct Quaternion
	{
		private const double eps = 10e-6;
		public double s;
		public MyVector3 v;

		public Quaternion(double s, MyVector3 v)
		{
			this.s = s;
			this.v = v;
		}
		public Quaternion(double a, double b, double c, double d)
		{
			this.s = a;
			this.v.x = b;
			this.v.y = c;
			this.v.z = d;
		}
		public Quaternion(MyMatrix3d m)
		{
			double T = 1.0 + m.Trace();
			double S, W, X, Y, Z;

			if (T > eps)
			{
				S = Math.Sqrt(T) * 2.0;
				X = (m[5] - m[7]) / S;
				Y = (m[6] - m[2]) / S;
				Z = (m[1] - m[3]) / S;
				W = 0.25 * S;
			}
			else if (m[0] > m[4] && m[0] > m[8])
			{
				S = Math.Sqrt(1.0 + m[0] - m[4] - m[8]) * 2.0;
				X = 0.25 * S;
				Y = (m[1] + m[3]) / S;
				Z = (m[6] + m[2]) / S;
				W = (m[5] - m[7]) / S;
			}
			else if (m[4] > m[8])
			{
				S = Math.Sqrt(1.0 + m[4] - m[0] - m[8]) * 2;
				X = (m[1] + m[3]) / S;
				Y = 0.25 * S;
				Z = (m[5] + m[7]) / S;
				W = (m[6] - m[2]) / S;
			}
			else
			{
				S = Math.Sqrt(1.0 + m[8] - m[0] - m[4]) * 2;
				X = (m[6] + m[2]) / S;
				Y = (m[5] + m[7]) / S;
				Z = 0.25 * S;
				W = (m[1] - m[3]) / S;
			}

			this.s = W;
			this.v = new MyVector3(X, Y, Z);
		}
		public MyMatrix3d ToMyMatrix3d()
		{

			double xx = v.x * v.x;
			double xy = v.x * v.y;
			double xz = v.x * v.z;
			double xw = v.x * s;
			double yy = v.y * v.y;
			double yz = v.y * v.z;
			double yw = v.y * s;
			double zz = v.z * v.z;
			double zw = v.z * s;
			MyMatrix3d m = new MyMatrix3d();
			m[0] = 1 - 2 * (yy + zz);
			m[1] = 2 * (xy + zw);
			m[2] = 2 * (xz - yw);
			m[3] = 2 * (xy - zw);
			m[4] = 1 - 2 * (xx + zz);
			m[5] = 2 * (yz + xw);
			m[6] = 2 * (xz + yw);
			m[7] = 2 * (yz - xw);
			m[8] = 1 - 2 * (xx + yy);
			return m;
		}
		public double this[int index]
		{
			get
			{
				if (index == 0) return s;
				return v[index-1];
			}
			set
			{
				if (index == 0) s = value;
				v[index - 1] = value;
			}
		}

		public Quaternion Conjugate()
		{
			return new Quaternion(s, -v.x, -v.y, -v.z);
		}
		public Quaternion Inverse()
		{
			double b = s*s + v.Dot(v);
			if (b != 0) return new Quaternion(s/b, -v.x/b, -v.y/b, -v.z/b);
			throw new DivideByZeroException();
		}
		public double Norm()
		{
			return Math.Sqrt(s * s + v.Dot(v));
		}

		static public Quaternion operator +(Quaternion a, Quaternion b)
		{
			return new Quaternion(a.s+b.s, a.v+b.v);
		}
		static public Quaternion operator -(Quaternion a, Quaternion b)
		{
			return new Quaternion(a.s-b.s, a.v-b.v);
		}
		static public Quaternion operator *(Quaternion a, Quaternion b)
		{
/*
			return new Quaternion(
				a.s*b.s - a.v.Dot(b.v),
				a[0]*b[1] + a[1]*b[0] + a[2]*b[3] - a[3]*b[2],
				a[0]*b[2] - a[1]*b[3] + a[2]*b[0] + a[3]*b[1],
				a[0]*b[3] + a[1]*b[2] - a[2]*b[1] + a[3]*b[0]
				);
*/

			return new Quaternion(
				a.s*b.s - a.v.Dot(b.v),
				a.s*b.v + b.s*a.v + a.v.Cross(b.v)
				);
		}
		static public Quaternion operator *(Quaternion q, double d)
		{
			return new Quaternion(q.s*d, q.v*d);
		}
		static public Quaternion operator /(Quaternion q, double d)
		{
			return new Quaternion(q.s / d, q.v / d);
		}
		static public Quaternion slerp(Quaternion a, Quaternion b, double t)
		{
			// assume the input quaternions are normalized
			double half_angle = Math.Acos(a.s*b.s + a.v.Dot(b.v));
			double sin_theta = Math.Sin(half_angle);
			double w1 = Math.Sin((1-t) * half_angle) / sin_theta;
			double w2 = Math.Sin(t * half_angle) / sin_theta;

			return a * w1 + b * w2;
		}
	}
}
