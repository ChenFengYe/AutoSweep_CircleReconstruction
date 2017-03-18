using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Drawing;
using System.Threading.Tasks;

using Loyc.Geometry;

using OpenTK.Graphics.OpenGL;

using SmartCanvas;

namespace MyGeometry
{
	class MyGCylinder
	{
		private List<MyCircle> circleList = new List<MyCircle>();

		public MyGCylinder () { } 
		public MyGCylinder (MyCircle c)
		{
			circleList.Add(c);
		}
		public MyGCylinder (List<MyCircle> list)
		{
			circleList.AddRange(list);
		}

		public List<MyCircle> CircleList
		{
			get { return circleList; }
		}
		public void AddCircle(MyCircle c)
		{
			circleList.Add(c);
		}
		public void Draw ()  // need modify
		{
			foreach (MyCircle c in circleList)
			{
				c.Draw();
			}
		}
	}
}
