using System;
using System.IO;
using System.Collections;
using System.Collections.Generic;
using System.Xml.Serialization;
using NumericalRecipes;

namespace MyGeometry
{
	public struct Matrix2d
	{
		private const int len = 4;
		private const int row_size = 2;
		double a, b, c, d;

		public Matrix2d(double [] arr) 
		{
			a = arr[0];
			b = arr[1];
			c = arr[2];
			d = arr[3];
		}
		public Matrix2d(double [,] arr) 
		{
			a = arr[0,0];
			b = arr[0,1];
			c = arr[1,0];
			d = arr[1,1];
		}

		// using column vectors
		public Matrix2d(MyVector2 v1, MyVector2 v2) 
		{
			a = v1.x;
			b = v2.x;
			c = v1.y;
			d = v2.y;
		}
 
		public Matrix2d(double a, double b, double c, double d) 
		{
			this.a = a;
			this.b = b;
			this.c = c;
			this.d = d;
		}

		public static Matrix2d operator + (Matrix2d m1, Matrix2d m2)
		{
			Matrix2d ret = new Matrix2d(
				m1.a + m2.a,
				m1.b + m2.b,
				m1.c + m2.c,
				m1.d + m2.d
				);
			return ret;
		}
		public static Matrix2d operator * (Matrix2d m1, Matrix2d m2) 
		{
			Matrix2d ret = new Matrix2d(
				m1.a * m2.a + m1.b * m2.c,
				m1.a * m2.b + m1.b * m2.d,
				m1.c * m2.a + m1.d * m2.c,
				m1.c * m2.b + m1.d * m2.d
				);
			return ret;
		}
		public static MyVector2 operator * (Matrix2d m, MyVector2 v) 
		{
			return new MyVector2(m.A*v.x+m.B*v.y, m.C*v.x+m.D*v.y);
		}
		public Matrix2d Inverse()
		{
			double det = (a*d - b*c);
			if (double.IsNaN(det)) throw new ArithmeticException();
			return new Matrix2d(d/det, -b/det, -c/det, a/det);
		}
		public double Det()
		{
			return (a * d - b * c);
		}
		public Matrix2d Transpose() 
		{
			return new Matrix2d(a, c, b, d);
		}
		public double Trace()
		{
			return a + d;
		}

		public double A
		{
			get { return a; }
			set { a = value; }
		}
		public double B
		{
			get { return b; }
			set { b = value; }
		}
		public double C
		{
			get { return c; }
			set { c = value; }
		}
		public double D
		{
			get { return d; }
			set { d = value; }
		}
		public override string ToString()
		{
			return 
				a.ToString("F5") + " " + b.ToString("F5") + " " +
				c.ToString("F5") + " " + d.ToString("F5");
		}
		public double[] ToArray()
		{
			double[] array = new double[4];
			array[0] = this.a; array[1] = this.b;
			array[2] = this.c; array[3] = this.d;
			return array;
		}
	}	
	
	public class MyMatrix3d
	{
		public static bool lastSVDIsFullRank = false;
		private const int len = 9;
		private const int row_size = 3;
        private double[] e = new double[len];

        public MyMatrix3d() { }
		public MyMatrix3d(double [] arr) 
		{
			for (int i=0; i<len; i++) e[i] = arr[i];
		}
		public MyMatrix3d(double [,] arr) 
		{
            for (int i=0; i<row_size; i++)
				for (int j=0; j<row_size; j++)
					this[i,j] = arr[i,j];
		}
		public MyMatrix3d(MyMatrix3d m) : this(m.e) { }
		public MyMatrix3d(MyVector3 v1, MyVector3 v2, MyVector3 v3)
		{
			for (int i=0; i<row_size; i++) 
			{
				this[i,0] = v1[i];
				this[i,1] = v2[i];
				this[i,2] = v3[i];
			}
		}

        public void Clear()
        {
            for (int i = 0; i < len; i++) e[i] = 0;
        }
        public double this[int index]
		{
			get { return e[index]; }
			set { e[index] = value; }
		}
		public double this[int row, int column] 
		{
			get { return e[row*row_size + column]; }
			set { e[row*row_size + column] = value; }
		}
		public double [] ToArray()
		{
			return (double[])e.Clone();
		}
		public double Det()
		{
			return e[0]*(e[4]*e[8]-e[5]*e[7])
				  -e[3]*(e[1]*e[8]-e[2]*e[7])
				  +e[6]*(e[1]*e[5]-e[2]*e[4]);
		}
		public double Trace()
		{
			return e[0] + e[4] + e[8];
		}
		public double SqNorm()
		{
			double sq = 0;
			for (int i=0; i<len; i++) sq += e[i] * e[i];
			return sq;
		}
		public MyMatrix3d Transpose()
		{
			MyMatrix3d m = new MyMatrix3d();
			for (int i = 0; i < row_size; i++)
				for (int j = 0; j < row_size; j++)
					m[j, i] = this[i, j];
			return m;
		}
		public MyMatrix3d Inverse()
		{
			double a = e[0];
			double b = e[1];
			double c = e[2];
			double d = e[3];
			double E = e[4];
			double f = e[5];
			double g = e[6];
			double h = e[7];
			double i = e[8];
			double det = a * (E * i - f * h) - b * (d * i - f * g) + c * (d * h - E * g);
			if (det == 0) throw new ArithmeticException();

			MyMatrix3d inv = new MyMatrix3d();
			inv[0] = (E * i - f * h) / det;
			inv[1] = (c * h - b * i) / det;
			inv[2] = (b * f - c * E) / det;
			inv[3] = (f * g - d * i) / det;
			inv[4] = (a * i - c * g) / det;
			inv[5] = (c * d - a * f) / det;
			inv[6] = (d * h - E * g) / det;
			inv[7] = (b * g - a * h) / det;
			inv[8] = (a * E - b * d) / det;
			return inv;
		}
		public MyMatrix3d InverseSVD()
		{
			SVD svd = new SVD(e, 3, 3);
			MyMatrix3d inv = new MyMatrix3d(svd.Inverse);
			lastSVDIsFullRank = svd.FullRank;
			return inv;
		}
		public MyMatrix3d InverseTranspose()
		{
			double a = e[0];
			double b = e[1];
			double c = e[2];
			double d = e[3];
			double E = e[4];
			double f = e[5];
			double g = e[6];
			double h = e[7];
			double i = e[8];
			double det = a * (E * i - f * h) - b * (d * i - f * g) + c * (d * h - E * g);
			if (det == 0) throw new ArithmeticException();

			MyMatrix3d inv = new MyMatrix3d();
			inv[0] = (E * i - f * h) / det;
			inv[3] = (c * h - b * i) / det;
			inv[6] = (b * f - c * E) / det;
			inv[1] = (f * g - d * i) / det;
			inv[4] = (a * i - c * g) / det;
			inv[7] = (c * d - a * f) / det;
			inv[2] = (d * h - E * g) / det;
			inv[5] = (b * g - a * h) / det;
			inv[8] = (a * E - b * d) / det;
			return inv;
		}
		public MyMatrix3d OrthogonalFactor(double eps)
		{
			MyMatrix3d Q = new MyMatrix3d(this);
			MyMatrix3d Q2 = new MyMatrix3d();
			double err = 0;
			do
			{
				Q2 = (Q + Q.InverseTranspose()) / 2.0;
				err = (Q2 - Q).SqNorm();
				Q = Q2;
			} while (err > eps);

			return Q2;
		}
		public MyMatrix3d OrthogonalFactorSVD()
		{
			SVD svd = new SVD(e, 3, 3);
			lastSVDIsFullRank = svd.FullRank;
			return new MyMatrix3d(svd.OrthogonalFactor);
		}
		public MyVector3 NullVector()
		{
			SVD svd = new SVD(e, 3, 3);
			for (int i = 1; i < 4; ++i)
			{
				if (Math.Abs(svd.w[i]) < 1e-7) // ==0
				{
					return new MyVector3(svd.u[1, i], svd.u[2, i], svd.u[3, i]);
				}
			}
			return new MyVector3(0, 0, 0);
		}
		public MyVector3 SmallestEigenVector()
		{
			SVD svd = new SVD(e, 3, 3);
			double min = double.MaxValue; int j = -1;
			for (int i = 1; i < 4; ++i)
			{
				if (svd.w[i] < min) // ==0
				{
					min = svd.w[i];
					j = i;
				}
			}
			return new MyVector3(svd.u[1, j], svd.u[2, j], svd.u[3, j]);
		}
        public double[] SVDSingularMat()
        {
            SVD svd = new SVD(e, 3, 3);
            return svd.w;
        }
		public MyMatrix3d OrthogonalFactorIter()
		{
			return (this + this.InverseTranspose())/2;
		}
		public static MyMatrix3d IdentityMatrix()
		{
			MyMatrix3d m = new MyMatrix3d();
			m[0] = m[4] = m[8] =  1.0;
			return m;
		}
		public MyMatrix3d LogMatrix(out MyVector3 axis, out double angle)
		{
			MyMatrix3d m;

			double cos = (this.Trace() - 1) / 2.0;
			if (cos < -1)
			{
				cos = -1;
			}
			else if (cos > 1)
			{
				cos = 1;
			}
			double theta = Math.Acos(cos);
			if (Math.Abs(theta) < 0.0001)
			{
				if (theta >= 0.0)
					theta = 0.0001;
				else
					theta = -0.0001;
			}

			m = (new MyMatrix3d(this) - this.Transpose()) * (0.5 / Math.Sin(theta) * theta);

			MyVector3 r = new MyVector3(m[2, 1], m[0, 2], m[1, 0]);

			angle = r.Length();

			axis = r / angle;

			return m;
		}
		public static MyVector3 operator * (MyMatrix3d m, MyVector3 v) 
		{
			MyVector3 ret = new MyVector3();
			ret.x = m[0]*v.x + m[1]*v.y + m[2]*v.z;
			ret.y = m[3]*v.x + m[4]*v.y + m[5]*v.z;
			ret.z = m[6]*v.x + m[7]*v.y + m[8]*v.z;
			return ret;
		}
		public static MyMatrix3d operator * (MyMatrix3d m1, MyMatrix3d m2)
		{
			MyMatrix3d ret = new MyMatrix3d();
			for (int i = 0; i < row_size; i++)
				for (int j = 0; j < row_size; j++)
				{
					ret[i, j] = 0.0;
					for (int k = 0; k < row_size; k++)
						ret[i, j] += m1[i, k] * m2[k, j];
				}
			return ret;
		}
		public static MyMatrix3d operator + (MyMatrix3d m1, MyMatrix3d m2)
		{
			MyMatrix3d ret = new MyMatrix3d();
			for (int i=0; i<len; i++) ret[i] = m1[i] + m2[i];
			return ret;
		}
		public static MyMatrix3d operator - (MyMatrix3d m1, MyMatrix3d m2)
		{
			MyMatrix3d ret = new MyMatrix3d();
			for (int i = 0; i < len; i++) ret[i] = m1[i] - m2[i];
			return ret;
		}
		public static MyMatrix3d operator * (MyMatrix3d m, double d)
		{
			MyMatrix3d ret = new MyMatrix3d();
			for (int i = 0; i < len; i++) ret[i] = m[i] * d;
			return ret;
		}
		public static MyMatrix3d operator / (MyMatrix3d m, double d)
		{
			MyMatrix3d ret = new MyMatrix3d();
			for (int i=0; i<len; i++) ret[i] = m[i] / d;
			return ret;
		}
		public override string ToString()
		{
			return
				e[0].ToString("F5") + " " + e[1].ToString("F5") + " " + e[2].ToString("F5") +
				e[3].ToString("F5") + " " + e[4].ToString("F5") + " " + e[5].ToString("F5") +
				e[6].ToString("F5") + " " + e[7].ToString("F5") + " " + e[08].ToString("F5");
		}
	}

	[XmlRootAttribute(IsNullable = false)]
	public class MyMatrix4d
	{
		private const int len = 16;
		private const int row_size = 4;
		private double[] e = new double[len];

		public MyMatrix4d()
		{
		}
		public MyMatrix4d(double[] arr)
		{
			for (int i = 0; i < len; i++) e[i] = arr[i];
		}
		public MyMatrix4d(double[,] arr)
		{
			for (int i = 0; i < row_size; i++)
				for (int j = 0; j < row_size; j++)
					this[i, j] = arr[i, j];
		}
		public MyMatrix4d(MyVector3 v1, MyVector3 v2, MyVector3 v3, MyVector3 v4)
		{
			for (int i = 0; i < 3; i++)
			{
				this[i, 0] = v1[i];
				this[i, 1] = v2[i];
				this[i, 2] = v3[i];
				this[i, 3] = v4[i];
			}
		}
		public MyMatrix4d(MyMatrix4d m) : this(m.e) { }
		public MyMatrix4d(MyMatrix3d m)
		{
			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++)
					this[i, j] = m[i, j];
		}
		public MyMatrix4d(StreamReader sr)
		{
			int c = 0;
			char[] delimiters = { ' ', '\t' };

			while (sr.Peek() > -1)
			{
				string s = sr.ReadLine();
				string[] tokens = s.Split(delimiters);
				for (int i = 0; i < tokens.Length; i++)
				{
					e[c++] = Double.Parse(tokens[i]);
					if (c >= len) return;
				}
			}
		}


		public double this[int index]
		{
			get { return e[index]; }
			set { e[index] = value; }
		}
		public double this[int row, int column]
		{
			get { return e[row * row_size + column]; }
			set { e[row * row_size + column] = value; }
		}
		public double[] Element
		{
			get { return e; }
			set
			{
				if (value.Length < len)
					throw new Exception();
				e = value;
			}
		}
		public double[] ToArray()
		{
			return (double[])e.Clone();
		}
		public MyMatrix3d ToMyMatrix3d()
		{
			MyMatrix3d ret = new MyMatrix3d();
			ret[0] = e[0];
			ret[1] = e[1];
			ret[2] = e[2];
			ret[3] = e[4];
			ret[4] = e[5];
			ret[5] = e[6];
			ret[6] = e[8];
			ret[7] = e[9];
			ret[8] = e[10];
			return ret;
		}
		public double Trace()
		{
			return e[0] + e[5] + e[10] + e[15];
		}
		public MyMatrix4d Inverse()
		{
			SVD svd = new SVD(e, row_size, row_size);
			if (svd.State == false) throw new ArithmeticException();
			return new MyMatrix4d(svd.Inverse);
		}
		public MyMatrix4d Inverse_alglib()
		{
			double[,] data = new double[row_size, row_size];
			for (int i = 0; i < row_size; ++i)
			{
				for (int j = 0; j < row_size; ++j)
					data[i, j] = this[i, j];
			}
			int info;
			alglib.matinvreport rep;
			alglib.rmatrixinverse(ref data, out info, out rep);

			if (info == 1)
				return new MyMatrix4d(data);
			else
				throw new ArithmeticException();
			
		}
		public MyMatrix4d Transpose()
		{
			MyMatrix4d m = new MyMatrix4d();
			for (int i = 0; i < row_size; i++)
				for (int j = 0; j < row_size; j++)
					m[j, i] = this[i, j];
			return m;
		}
		public static MyMatrix4d IdentityMatrix()
		{
			MyMatrix4d m = new MyMatrix4d();
			m[0] = m[5] = m[10] = m[15] = 1.0;
			return m;
		}
		public static MyMatrix4d CompressMatrix(double alpha)
		{
			MyMatrix4d m = IdentityMatrix();

			m[3, 2] = alpha;
			return m;
		}
		public static MyMatrix4d RotationMatrix(MyVector3 axis, double angle)
		{
			MyMatrix4d m = IdentityMatrix();

			double c = Math.Cos(angle);
			double s = Math.Sin(angle);

			axis.Normalize();

			double nx = axis[0];
			double ny = axis[1];
			double nz = axis[2];

			m[0, 0] = c + (1 - c) * nx * nx;
			m[0, 1] = -s * nz + (1 - c) * nx * ny;
			m[0, 2] = s * ny + (1 - c) * nx * nz;
			m[0, 3] = 0.0;

			m[1, 0] = s * nz + (1 - c) * nx * ny;
			m[1, 1] = c + (1 - c) * ny * ny;
			m[1, 2] = -s * nx + (1 - c) * ny * nz;
			m[1, 3] = 0.0;

			m[2, 0] = -s * ny + (1 - c) * nz * nx;
			m[2, 1] = s * nx + (1 - c) * nz * ny;
			m[2, 2] = c + (1 - c) * nz * nz;
			m[2, 3] = 0.0;

			m[3, 0] = 0.0;
			m[3, 1] = 0.0;
			m[3, 2] = 0.0;
			m[3, 3] = 1.0;

			return m;
		}
		public static MyMatrix4d RotationMatrix(MyVector3 axis, double cos, double sin)
		{
			MyMatrix4d m = IdentityMatrix();

			double c = cos;
			double s = sin;

			axis.Normalize();

			double nx = axis[0];
			double ny = axis[1];
			double nz = axis[2];

			m[0, 0] = c + (1 - c) * nx * nx;
			m[0, 1] = -s * nz + (1 - c) * nx * ny;
			m[0, 2] = s * ny + (1 - c) * nx * nz;
			m[0, 3] = 0.0;

			m[1, 0] = s * nz + (1 - c) * nx * ny;
			m[1, 1] = c + (1 - c) * ny * ny;
			m[1, 2] = -s * nx + (1 - c) * ny * nz;
			m[1, 3] = 0.0;

			m[2, 0] = -s * ny + (1 - c) * nz * nx;
			m[2, 1] = s * nx + (1 - c) * nz * ny;
			m[2, 2] = c + (1 - c) * nz * nz;
			m[2, 3] = 0.0;

			m[3, 0] = 0.0;
			m[3, 1] = 0.0;
			m[3, 2] = 0.0;
			m[3, 3] = 1.0;

			return m;
		}
		public static MyMatrix4d RotationMatrixU2V(MyVector3 u, MyVector3 v)
		{
			// find the rotational matrix which rotate u to v
			// one should be extremely careful here, very small viboration
			// will make a lot of difference
			// e.g., u = (0.5*e-10, 0.1*e-10, 1), v = (0,0,-1), will make
			// an fliping around axie (1, 5, 0) with angele Pi
			MyMatrix4d R = MyMatrix4d.IdentityMatrix();
			double cos = u.Normalize().Dot(v.Normalize());
			if (Math.Abs(Math.Abs(cos) - 1.0f) < 1e-5) // coincident, do nothing
				return R;
			MyVector3 axis = u.Cross(v).Normalize();
			if (!double.IsNaN(axis.x)) {
				if (cos < -1) cos = -1;
				if (cos > 1) cos = 1;
				double angle = Math.Acos(cos);
				R = MyMatrix4d.RotationMatrix(axis, angle);
			}
			return R;
		}
		public static void FindRotAxisAngle(MyMatrix4d rotMat, out MyVector3 rotAxis, out double angle)
		{
			double trace = rotMat[0, 0] + rotMat[1, 1] + rotMat[2, 2];
			angle = Math.Acos((trace - 1) / 2);

			if (rotMat[0, 1] > 0) angle = -angle; /// this may be violated...

			double sin2 = 2 * Math.Sin(angle);
			double x = (rotMat[2, 1] - rotMat[1, 2]) / sin2;
			double y = (rotMat[0, 2] - rotMat[2, 0]) / sin2;
			double z = (rotMat[1, 0] - rotMat[0, 1]) / sin2;

			rotAxis = new MyVector3(x, y, z).Normalize();
		}
		public static MyMatrix4d TranslationMatrix(MyVector3 t)
		{
			MyMatrix4d m = IdentityMatrix();

			m[0, 3] = t[0];
			m[1, 3] = t[1];
			m[2, 3] = t[2];

			return m;
		}
		public static MyVector3 GetTranslation(MyMatrix4d T)
		{
			return new MyVector3(T[0,3], T[1,3], T[2,3]);
		}
		public static MyMatrix4d ScalingMatrix(double sx, double sy, double sz)
		{
			MyMatrix4d m = IdentityMatrix();

			m[0, 0] = sx;
			m[1, 1] = sy;
			m[2, 2] = sz;
			m[3, 3] = 1.0;

			return m;
		}
		// compute the releft point of p according to the Canvas3 (Canvas3_normal, Canvas3_center)
		public static MyVector3 Reflect(MyVector3 p, MyVector3 Canvas3_normal, MyVector3 Canvas3_center)
		{
			MyVector3 u = p - Canvas3_center;
			
			// create a coord system (x,y,z), project to that sys, reflect and project back
			MyVector3 x = Canvas3_normal;
			MyVector3 y;
			if (x.x == 0 && x.y == 0)
				y = new MyVector3(0, -x.z, x.y);
			else
				y = new MyVector3(x.y, -x.x, 0);
			MyVector3 z = x.Cross(y);
			MyMatrix3d R = new MyMatrix3d(x, y, z);
			MyMatrix3d InvR = R.InverseSVD();
			MyMatrix4d U = new MyMatrix4d(R);
			MyMatrix4d V = new MyMatrix4d(InvR);
			MyMatrix4d I = MyMatrix4d.IdentityMatrix();
			I[0, 0] = -1; // reflect matrix along yz Canvas3
			MyMatrix4d T = V * I * U; // the reflection matrix
			
			// reflect
			MyVector4 r = new MyVector4(u, 0);
			MyVector3 q = Canvas3_center+(T * r).XYZ();

			return q;
		}
		public static MyVector4 operator *(MyVector4 v, MyMatrix4d m)
		{
			return m.Transpose() * v;
		}
		public static MyMatrix4d operator *(MyMatrix4d m1, MyMatrix4d m2)
		{
			MyMatrix4d ret = new MyMatrix4d();
			for (int i = 0; i < row_size; i++)
				for (int j = 0; j < row_size; j++)
				{
					ret[i, j] = 0.0;
					for (int k = 0; k < row_size; k++)
						ret[i, j] += m1[i, k] * m2[k, j];
				}
			return ret;
		}
		public static MyVector4 operator *(MyMatrix4d m, MyVector4 v)
		{
			MyVector4 ret = new MyVector4();
			ret.x = m[0] * v.x + m[1] * v.y + m[2] * v.z + m[3] * v.w;
			ret.y = m[4] * v.x + m[5] * v.y + m[6] * v.z + m[7] * v.w;
			ret.z = m[8] * v.x + m[9] * v.y + m[10] * v.z + m[11] * v.w;
			ret.w = m[12] * v.x + m[13] * v.y + m[14] * v.z + m[15] * v.w;
			return ret;
		}
		public static MyMatrix4d operator +(MyMatrix4d m1, MyMatrix4d m2)
		{
			MyMatrix4d ret = new MyMatrix4d();
			for (int i = 0; i < len; i++) ret[i] = m1[i] + m2[i];
			return ret;
		}
		public static MyMatrix4d operator -(MyMatrix4d m1, MyMatrix4d m2)
		{
			MyMatrix4d ret = new MyMatrix4d();
			for (int i = 0; i < len; i++) ret[i] = m1[i] - m2[i];
			return ret;
		}
		public static MyMatrix4d operator *(MyMatrix4d m, double d)
		{
			MyMatrix4d ret = new MyMatrix4d();
			for (int i = 0; i < len; i++) ret[i] = m[i] * d;
			return ret;
		}

		public override string ToString()
		{
			string s = "";
			foreach (double d in e)
				s += d.ToString() + "\n";
			return s;
		}
	}

	public class MatrixNd
	{
		private int m;
		private int n;
		private double[] e;

		public MatrixNd(int m, int n)
		{
			this.m = m;
			this.n = n;
			e = new double[m * n];
			for (int i = 0; i < m * n; i++)
				e[i] = 0;
		}
		public MatrixNd(SparseMatrix right) 
			: this(right.RowSize, right.ColumnSize)
		{
			int b = 0;
			foreach (List<SparseMatrix.Element> row in right.Rows)
			{
				foreach (SparseMatrix.Element element in row)
					e[b + element.j] = element.value;
				b += n;
			}
		}
		public MatrixNd(SparseMatrix right, bool transpose)
			: this(right.ColumnSize, right.RowSize)
		{
			int b = 0;
			foreach (List<SparseMatrix.Element> col in right.Columns)
			{
				foreach (SparseMatrix.Element element in col)
					e[b + element.i] = element.value;
				b += n;
			}
		}
        public MatrixNd(MyVector3[] v)
        {
            this.m = 3;
            this.n = v.Length;
            e = new double[m * n];
            for (int i = 0; i < m; ++i)
            {
                for (int j = 0; j < n; ++j)
                {
                    this[i, j] = v[j][i];
                }
            }
        }

		public MatrixNd(double[,] val)
		{
			this.m = val.GetLength(0);
			this.n = val.GetLength(1);
			this.e = new double[this.m * this.n];
			for (int i = 0; i < this.m; ++i)
			{
				for (int j = 0; j < this.n; ++j)
				{
					e[i + j * this.n] = val[i, j];
				}
			}
		}

		public double this[int row, int column]
		{
			get { return e[row * n + column]; }
			set { e[row * n + column] = value;  }
		}
		public int RowSize
		{
			get { return m; }
		}
		public int ColumnSize
		{
			get { return n; }
		}

		public void Multiply(double[] xIn, double[] xOut)
		{
			if (xIn.Length<n || xOut.Length<m) throw new ArgumentException();

			for (int i=0,b=0; i<m; i++,b+=n)
			{
				double sum = 0;
				for (int j = 0; j < n; j++)
					sum += e[b + j] * xIn[j];
				xOut[i] = sum;
			}
		}

        public MatrixNd Mul(MatrixNd right)
        {
            if (right.RowSize != this.ColumnSize) throw new ArgumentException();
            int l = right.ColumnSize;
            MatrixNd matrix = new MatrixNd(m, l);
            for (int i = 0; i < m; ++i)
            {
                for (int j = 0; j < l; ++j)
                {
                    for (int k = 0; k < n; ++k)
                    {
                        matrix[i, j] += this[i, k] * right[k, j]; 
                    }
                }
            }
            return matrix;
        }
        public MatrixNd Transpose()
        {
            MatrixNd matrix = new MatrixNd(n, m);
            for (int i = 0; i < m; ++i)
            {
                for (int j = 0; j < n; ++j)
                {
                    matrix[j, i] = this[i, j];
                }
            }
            return matrix;
        }
        public MyMatrix3d ToMyMatrix3d()
        {
            MyMatrix3d ret = new MyMatrix3d();
            ret[0] = this[0, 0];
            ret[1] = this[0, 1];
            ret[2] = this[0, 2];
            ret[3] = this[1, 0];
            ret[4] = this[1, 1];
            ret[5] = this[1, 2];
            ret[6] = this[2, 0];
            ret[7] = this[2, 1];
            ret[8] = this[2, 2];
            return ret;
        }
		public Matrix2d ToMatrix2D()
		{
			return new Matrix2d(this[0, 0], this[0, 1], this[1, 0], this[1, 1]);
		}
    }
}
