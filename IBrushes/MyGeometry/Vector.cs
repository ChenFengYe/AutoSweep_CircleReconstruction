using System;

namespace MyGeometry
{
	public struct MyVector2
	{
		public static MyVector2 MinValue = new MyVector2(double.MinValue, double.MinValue);
		public static MyVector2 MaxValue = new MyVector2(double.MaxValue, double.MaxValue);
		public static MyVector2 Zero = new MyVector2(0, 0);


		public double x, y;

		public MyVector2(double x, double y)
		{
			this.x = x;
			this.y = y;
		}

		public MyVector2(MyVector2 v)
		{
			this.x = v.x;
			this.y = v.y;
		}

		public MyVector2(double[] arr, int index)
		{
			x = arr[index];
			y = arr[index+1];
		}

		public double[] ToArray()
		{
			return new double[] { this.x, this.y };
		}

		public double this[int index] 
		{
			get 
			{
				if (index==0) return x;
				if (index==1) return y;
				throw new ArgumentException();
			}
			set 
			{
				if (index==0) x = value;
				if (index==1) y = value;
			}
		}
		public double Dot (MyVector2 v)
		{
			return x*v.x + y*v.y;
		}
		public double Length()
		{
			return Math.Sqrt(x*x + y*y);
		}
        public double SquareLength()
		{
			return x * x + y * y;
		}

		public MyVector2 Normalize()
		{
			return this / this.Length();
		}

		public override string ToString()
		{
			return x.ToString() + " " + y.ToString();
		}


		public static MyVector2 Max(MyVector2 v1, MyVector2 v2) 
		{
			return new MyVector2( 
				(v1.x > v2.x) ? v1.x : v2.x,
				(v1.y > v2.y) ? v1.y : v2.y
				);
		}
		public static MyVector2 Min(MyVector2 v1, MyVector2 v2) 
		{
			return new MyVector2(
				(v1.x < v2.x) ? v1.x : v2.x,
				(v1.y < v2.y) ? v1.y : v2.y
				);
		}
		public bool Equals(MyVector2 v)
		{
			return (this.x == v.x) && (this.y == v.y) ? true : false;
		}

		public double Cross(MyVector2 v)
		{
			return x * v.y - y * v.x;
		}

		public MyVector2 PerpendicularLeft
		{
			get { return new MyVector2(-y, x); }
		}

		public MyVector2 PerpendicularRight
		{
			get { return new MyVector2(y, -x); }
		}

		static public MyVector2 operator+ (MyVector2 v1, MyVector2 v2)
		{
			return new MyVector2(v1.x+v2.x, v1.y+v2.y);
		}
		static public MyVector2 operator- (MyVector2 v1, MyVector2 v2)
		{
			return new MyVector2(v1.x-v2.x, v1.y-v2.y);
		}
		static public MyVector2 operator* (MyVector2 v, double s)
		{
			return new MyVector2(v.x*s, v.y*s);
		}
		static public MyVector2 operator* (double s, MyVector2 v)
		{
			return new MyVector2(v.x*s, v.y*s);
		}
		static public MyVector2 operator/ (MyVector2 v, double s)
		{
			return new MyVector2(v.x/s, v.y/s);
		}

		static public bool operator == (MyVector2 v1, MyVector2 v2)
		{
			return (v1.x == v2.x) && (v1.y == v2.y);
		}
		static public bool operator !=(MyVector2 v1, MyVector2 v2)
		{
			return !(v1 == v2);
		}
		public bool EqualsPoint(MyVector2 v)
		{
			if((this-v).Length()<1e-6)
			{
				return true;
			}
			else
			{
				return false;
			}
		}
	}

	public struct MyVector3
	{
		public static MyVector3 MinValue = new MyVector3(double.MinValue, double.MinValue, double.MinValue);
		public static MyVector3 MaxValue = new MyVector3(double.MaxValue, double.MaxValue, double.MaxValue);
		public static MyVector3 Xaxis = new MyVector3(1, 0, 0);
		public static MyVector3 Yaxis = new MyVector3(0, 1, 0);
		public static MyVector3 Zaxis = new MyVector3(0, 0, 1);
		public static MyVector3[] XYZaxis = new MyVector3[3] {
			new MyVector3(1,0,0),
			new MyVector3(0,1,0),
			new MyVector3(0,0,1),
		};
		public double x, y, z;

		public MyVector3(double x, double y, double z)
		{
			this.x = x;
			this.y = y;
			this.z = z;
		}
		public MyVector3(MyVector2 v)
		{
			this.x = v.x;
			this.y = v.y;
			this.z = 0;
		}
		public MyVector3(MyVector3 v)
		{
			this.x = v.x;
			this.y = v.y;
			this.z = v.z;
		}
		public MyVector3(MyVector2 v, double z)
		{
			this.x = v.x;
			this.y = v.y;
			this.z = z;
		}

		public MyVector3(double[] arr, int index)
		{
			x = arr[index];
			y = arr[index+1];
			z = arr[index+2];
		}


		public double this[int index] 
		{
			get 
			{
				if (index==0) return x;
				if (index==1) return y;
				if (index==2) return z;
				throw new ArgumentException();
			}
			set 
			{
				if (index==0) x = value;
				if (index==1) y = value;
				if (index==2) z = value;
			}
		}
		public double Dot (MyVector3 v)
		{
			return x*v.x + y*v.y + z*v.z;
		}
		public double Length()
		{
			return Math.Sqrt(x*x + y*y + z*z);
		}
		public double SquareLength()
		{
			return (x * x + y * y + z * z);
		}
		public MyVector3 HomogenousNormalize()
		{
			if (this.z != 0)
			{
				this.x /= this.z;
				this.y /= this.z;
				this.z = 1.0;
			}
			return this;
		}
		public MyVector2 ToMyVector2()
		{
			return new MyVector2(this.x, this.y);
		}
        public MyVector4 ToMyVector4()
        {
            return new MyVector4(this.x, this.y, this.y, 1);
        }
		public double[] ToArray()
		{
			return new double[3] { x, y, z };
		}
		public bool IsNull()
		{
			return x == 0 && y == 0 && z == 0;
		}
		public MyVector3 Cross(MyVector3 v)
		{
			return new MyVector3(
				y * v.z - v.y * z,
				z * v.x - v.z * x,
				x * v.y - v.x * y
				);
		}
		public MyVector3 Normalize()
		{
			return this / this.Length();
		}
		
		public MyVector3 Rotate(MyVector3 axis, double cos, double sin)
		{
			return this * cos + (axis.Cross(this)) * sin + 
				axis * ((1.0 - cos) * (axis.Dot(this)));
		}

		public MyMatrix3d OuterCross(MyVector3 v)
		{
			MyMatrix3d m = new MyMatrix3d();
			m[0, 0] = x * v.x;
			m[0, 1] = x * v.y;
			m[0, 2] = x * v.z;
			m[1, 0] = y * v.x;
			m[1, 1] = y * v.y;
			m[1, 2] = y * v.z;
			m[2, 0] = z * v.x;
			m[2, 1] = z * v.y;
			m[2, 2] = z * v.z;
			return m;
		}
		
		public MyVector2 XY()
		{
			return new MyVector2(x, y);
		}

		public override string ToString()
		{
			return x.ToString() + " " + y.ToString() + " " + z.ToString();
		}

		public static MyVector3 Max(MyVector3 v1, MyVector3 v2) 
		{
			return new MyVector3( (v1.x > v2.x) ? v1.x : v2.x,
				(v1.y > v2.y) ? v1.y : v2.y,
				(v1.z > v2.z) ? v1.z : v2.z );
		}
		public static MyVector3 Min(MyVector3 v1, MyVector3 v2) 
		{
			return new MyVector3( (v1.x < v2.x) ? v1.x : v2.x,
				(v1.y < v2.y) ? v1.y : v2.y,
				(v1.z < v2.z) ? v1.z : v2.z );
		}

		static public MyVector3 operator+ (MyVector3 v1, MyVector3 v2)
		{
			return new MyVector3(v1.x+v2.x, v1.y+v2.y, v1.z+v2.z);
		}
		static public MyVector3 operator- (MyVector3 v1, MyVector3 v2)
		{
			return new MyVector3(v1.x-v2.x, v1.y-v2.y, v1.z-v2.z);
		}
		static public MyVector3 operator* (MyVector3 v, double s)
		{
			return new MyVector3(v.x*s, v.y*s, v.z*s);
		}
		static public MyVector3 operator* (double s, MyVector3 v)
		{
			return new MyVector3(v.x*s, v.y*s, v.z*s);
		}
		static public MyVector3 operator/ (MyVector3 v, double s)
		{
			return new MyVector3(v.x/s, v.y/s, v.z/s);
		}
	}

	public struct MyVector4
	{
		public double x, y, z, w;

		public MyVector4(double x, double y, double z, double w)
		{
			this.x = x;
			this.y = y;
			this.z = z;
			this.w = w;
		}
		public MyVector4(MyVector2 v)
		{
			this.x = v.x;
			this.y = v.y;
			this.z = 0;
			this.w = 0;
		}
		public MyVector4(MyVector2 v, double z, double w)
		{
			this.x = v.x;
			this.y = v.y;
			this.z = z;
			this.w = w;
		}
		public MyVector4(MyVector3 v)
		{
			this.x = v.x;
			this.y = v.y;
			this.z = v.z;
			this.w = 0;
		}
		public MyVector4(MyVector3 v, double w)
		{
			this.x = v.x;
			this.y = v.y;
			this.z = v.z;
			this.w = w;
		}

		public MyVector4(double[] arr, int index)
		{
			x = arr[index];
			y = arr[index + 1];
			z = arr[index + 2];
			w = arr[index + 3];
		}


		public double Dot(MyVector4 v)
		{
			return x * v.x + y * v.y + z * v.z + w * v.w;
		}
		public double Length()
		{
			return Math.Sqrt(x * x + y * y + z * z + w * w);
		}

		public MyVector4 Normalize()
		{
			return this / this.Length();
		}
		public MyMatrix4d OuterCross(MyVector4 v)
		{
			MyMatrix4d m = new MyMatrix4d();
			m[0, 0] = x * v.x;
			m[0, 1] = x * v.y;
			m[0, 2] = x * v.z;
			m[0, 3] = x * v.w;
			m[1, 0] = y * v.x;
			m[1, 1] = y * v.y;
			m[1, 2] = y * v.z;
			m[1, 3] = y * v.w;
			m[2, 0] = z * v.x;
			m[2, 1] = z * v.y;
			m[2, 2] = z * v.z;
			m[2, 3] = z * v.w;
			m[3, 0] = w * v.x;
			m[3, 1] = w * v.y;
			m[3, 2] = w * v.z;
			m[3, 3] = w * v.w;
			return m;
		}

		public MyVector3 XYZ()
		{
			return new MyVector3(x, y, z);
		}

		public override string ToString()
		{
			return x.ToString() + " " + y.ToString() + " " + z.ToString() + " " + w.ToString();
		}
		public static MyVector4 Max(MyVector4 v1, MyVector4 v2)
		{
			return new MyVector4(
				(v1.x > v2.x) ? v1.x : v2.x,
				(v1.y > v2.y) ? v1.y : v2.y,
				(v1.z > v2.z) ? v1.z : v2.z,
				(v1.w > v2.w) ? v1.w : v2.w
				);
		}
		public static MyVector4 Min(MyVector4 v1, MyVector4 v2)
		{
			return new MyVector4(
				(v1.x < v2.x) ? v1.x : v2.x,
				(v1.y < v2.y) ? v1.y : v2.y,
				(v1.z < v2.z) ? v1.z : v2.z,
				(v1.w < v2.w) ? v1.w : v2.w
				);
		}

		static public MyVector4 operator +(MyVector4 v1, MyVector4 v2)
		{
			return new MyVector4(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z, v1.w + v2.w);
		}
		static public MyVector4 operator -(MyVector4 v1, MyVector4 v2)
		{
			return new MyVector4(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z, v1.w - v2.w);
		}
		static public MyVector4 operator *(MyVector4 v, double s)
		{
			return new MyVector4(v.x * s, v.y * s, v.z * s, v.w * s);
		}
		static public MyVector4 operator *(double s, MyVector4 v)
		{
			return new MyVector4(v.x * s, v.y * s, v.z * s, v.w * s);
		}
		static public MyVector4 operator /(MyVector4 v, double s)
		{
			return new MyVector4(v.x / s, v.y / s, v.z / s, v.w / s);
		}
	}
}