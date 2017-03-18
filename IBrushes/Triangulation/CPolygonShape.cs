using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using MyGeometry;

namespace PolygonDecompose
{
	public enum VertexType
	{
		ErrorPoint,
		ConvexPoint,
		ConcavePoint
	}

	public enum PolygonType
	{
		Unknown,
		Convex,
		Concave
	}

	public enum PolygonDirection
	{
		Unknown,
		Clockwise,
		Count_Clockwise
	}



	public class CPolygonShape
	{

		private MyVector2[] m_aInputVertices;
		private MyVector2[] m_aUpdatedPolygonVertices;

		private System.Collections.ArrayList m_alEars
			= new System.Collections.ArrayList();
		private MyVector2[][] m_aPolygons;

		public int NumberOfPolygons
		{
			get
			{
				return  m_aPolygons.Length;
			}
		}

		public MyVector2[] Polygons(int index)
		{
			if (index < m_aPolygons.Length)
				return m_aPolygons[index];
			else
				return null;
		}

		public CPolygonShape(MyVector2[] vertices)
		{
			int nVertices = vertices.Length;
			if (nVertices < 3)
			{
				System.Diagnostics.Trace.WriteLine("To make a polygon, "
					+ " at least 3 points are required!");
				return;
			}

			//initialize the 2D points
			m_aInputVertices = new MyVector2[nVertices];

			for (int i = 0; i < nVertices; i++)
				m_aInputVertices[i] = vertices[i];

			//make a working copy,  m_aUpdatedPolygonVertices are
			//in count clock direction from user view
			SetUpdatedPolygonVertices();
		}

		/****************************************************
		To fill m_aUpdatedPolygonVertices array with input array.
		
		m_aUpdatedPolygonVertices is a working array that will 
		be updated when an ear is cut till m_aUpdatedPolygonVertices
		makes triangle (a convex polygon).
	   ******************************************************/
		private void SetUpdatedPolygonVertices()
		{
			int nVertices = m_aInputVertices.Length;
			m_aUpdatedPolygonVertices = new MyVector2[nVertices];

			for (int i = 0; i < nVertices; i++)
				m_aUpdatedPolygonVertices[i] = m_aInputVertices[i];

			//m_aUpdatedPolygonVertices should be in count clock wise
			if (CPolygon.PointsDirection(m_aUpdatedPolygonVertices)
				== PolygonDirection.Clockwise)
				CPolygon.ReversePointsDirection(m_aUpdatedPolygonVertices);
		}

		/**********************************************************
		To check the Pt is in the Triangle or not.
		If the Pt is in the line or is a vertex, then return true.
		If the Pt is out of the Triangle, then return false.

		This method is used for triangle only.
		***********************************************************/
		private bool TriangleContainsPoint(MyVector2[] trianglePts, MyVector2 pt)
		{
			if (trianglePts.Length != 3)
				return false;

			for (int i = trianglePts.GetLowerBound(0);
				i < trianglePts.GetUpperBound(0); i++)
			{
				if (pt.EqualsPoint(trianglePts[i]))
					return true;
			}

			bool bIn = false;

			CLineSegment line0 = new CLineSegment(trianglePts[0], trianglePts[1]);
			CLineSegment line1 = new CLineSegment(trianglePts[1], trianglePts[2]);
			CLineSegment line2 = new CLineSegment(trianglePts[2], trianglePts[0]);

			if (line0.PointInLine(pt) || line1.PointInLine(pt) || line2.PointInLine(pt))
				bIn = true;
			else //point is not in the lines
			{
				double dblArea0 = CPolygon.PolygonArea(new MyVector2[] { trianglePts[0], trianglePts[1], pt });
				double dblArea1 = CPolygon.PolygonArea(new MyVector2[] { trianglePts[1], trianglePts[2], pt });
				double dblArea2 = CPolygon.PolygonArea(new MyVector2[] { trianglePts[2], trianglePts[0], pt });

				if (dblArea0 > 0)
				{
					if ((dblArea1 > 0) && (dblArea2 > 0))
						bIn = true;
				}
				else if (dblArea0 < 0)
				{
					if ((dblArea1 < 0) && (dblArea2 < 0))
						bIn = true;
				}
			}
			return bIn;
		}


		/****************************************************************
		To check whether the Vertex is an ear or not based updated Polygon vertices

		ref. www-cgrl.cs.mcgill.ca/~godfried/teaching/cg-projects/97/Ian
		/algorithm1.html

		If it is an ear, return true,
		If it is not an ear, return false;
		*****************************************************************/
		private bool IsEarOfUpdatedPolygon(MyVector2 vertex)
		{
			CPolygon polygon = new CPolygon(m_aUpdatedPolygonVertices);

			if (polygon.PolygonVertex(vertex))
			{
				bool bEar = true;
				if (polygon.PolygonVertexType(vertex) == VertexType.ConvexPoint)
				{
					MyVector2 pi = vertex;
					MyVector2 pj = polygon.PreviousPoint(vertex); //previous vertex
					MyVector2 pk = polygon.NextPoint(vertex);//next vertex

					for (int i = m_aUpdatedPolygonVertices.GetLowerBound(0);
						i < m_aUpdatedPolygonVertices.GetUpperBound(0); i++)
					{
						MyVector2 pt = m_aUpdatedPolygonVertices[i];
						if (!(pt.EqualsPoint(pi) || pt.EqualsPoint(pj) || pt.EqualsPoint(pk)))
						{
							MyVector2[] vecs = new MyVector2[3] { pj, pi, pk };
							if (TriangleContainsPoint(vecs, pt))
								bEar = false;
						}
					}
				} //ThePolygon.getVertexType(Vertex)=ConvexPt
				else  //concave point
					bEar = false; //not an ear/
				return bEar;
			}
			else //not a polygon vertex;
			{
				System.Diagnostics.Trace.WriteLine("IsEarOfUpdatedPolygon: " +
					"Not a polygon vertex");
				return false;
			}
		}

		/****************************************************
		Set up m_aPolygons:
		add ears and been cut Polygon togather
		****************************************************/
		private void SetPolygons()
		{
			int nPolygon = m_alEars.Count + 1; //ears plus updated polygon
			m_aPolygons = new MyVector2[nPolygon][];

			for (int i = 0; i < nPolygon - 1; i++) //add ears
			{
				MyVector2[] points = (MyVector2[])m_alEars[i];

				m_aPolygons[i] = new MyVector2[3]; //3 vertices each ear
				m_aPolygons[i][0] = points[0];
				m_aPolygons[i][1] = points[1];
				m_aPolygons[i][2] = points[2];
			}

			//add UpdatedPolygon:
			m_aPolygons[nPolygon - 1] = new
				MyVector2[m_aUpdatedPolygonVertices.Length];

			for (int i = 0; i < m_aUpdatedPolygonVertices.Length; i++)
			{
				m_aPolygons[nPolygon - 1][i] = m_aUpdatedPolygonVertices[i];
			}
		}

		/********************************************************
		To update m_aUpdatedPolygonVertices:
		Take out Vertex from m_aUpdatedPolygonVertices array, add 3 points
		to the m_aEars
		**********************************************************/
		private void UpdatePolygonVertices(MyVector2 vertex)
		{
			System.Collections.ArrayList alTempPts = new System.Collections.ArrayList();

			for (int i = 0; i < m_aUpdatedPolygonVertices.Length; i++)
			{
				if (vertex.EqualsPoint(
					m_aUpdatedPolygonVertices[i])) //add 3 pts to FEars
				{
					CPolygon polygon = new CPolygon(m_aUpdatedPolygonVertices);
					MyVector2 pti = vertex;
					MyVector2 ptj = polygon.PreviousPoint(vertex); //previous point
					MyVector2 ptk = polygon.NextPoint(vertex); //next point

					MyVector2[] aEar = new MyVector2[3]; //3 vertices of each ear
					aEar[0] = ptj;
					aEar[1] = pti;
					aEar[2] = ptk;

					m_alEars.Add(aEar);
				}
				else
				{
					alTempPts.Add(m_aUpdatedPolygonVertices[i]);
				} //not equal points
			}

			if (m_aUpdatedPolygonVertices.Length
				- alTempPts.Count == 1)
			{
				int nLength = m_aUpdatedPolygonVertices.Length;
				m_aUpdatedPolygonVertices = new MyVector2[nLength - 1];

				for (int i = 0; i < alTempPts.Count; i++)
					m_aUpdatedPolygonVertices[i] = (MyVector2)alTempPts[i];
			}
		}


		/*******************************************************
		To cut an ear from polygon to make ears and an updated polygon:
		*******************************************************/
		public bool CutEar()
		{
			CPolygon polygon = new CPolygon(m_aUpdatedPolygonVertices);
			bool bFinish = false;

			//if (polygon.GetPolygonType()==PolygonType.Convex) //don't have to cut ear
			//	bFinish=true;

			if (m_aUpdatedPolygonVertices.Length == 3) //triangle, don't have to cut ear
				bFinish = true;

			MyVector2 pt = new MyVector2();
			while (bFinish == false) //UpdatedPolygon
			{
				int i = 0;
				bool bNotFound = true;
				int n = m_aUpdatedPolygonVertices.Length;
				while (bNotFound
					&& (i < m_aUpdatedPolygonVertices.Length)) //loop till find an ear
				{
					pt = m_aUpdatedPolygonVertices[i];
					if (IsEarOfUpdatedPolygon(pt))
						bNotFound = false; //got one, pt is an ear
					else
						i++;
				} //bNotFount
				//An ear found:}
				//if(bNotFound)
				//{
				//	return false;
				//}
				if (pt != null)
					UpdatePolygonVertices(pt);
				if(Double.IsNaN(pt[0]) || Double.IsNaN(pt[1]))
				{
					return false;
				}

				polygon = new CPolygon(m_aUpdatedPolygonVertices);
				//if ((polygon.GetPolygonType()==PolygonType.Convex)
				//	&& (m_aUpdatedPolygonVertices.Length==3))
				if (m_aUpdatedPolygonVertices.Length == 3)
					bFinish = true;

				if (m_aUpdatedPolygonVertices.Length == n)
					return false;
			} //bFinish=false
			SetPolygons();
			return true;
		}
	}
}
