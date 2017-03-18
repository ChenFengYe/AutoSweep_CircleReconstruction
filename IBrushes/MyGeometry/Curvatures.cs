using System;
using System.Collections;
using System.Collections.Generic;
using MyGeometry;
using NumericalRecipes;


namespace MyGeometry
{
	public class Curvatures : IDisposable
	{
		private Mesh mesh = null;

		private double[] curv1 = null;		// max principle
		private double[] curv2 = null;		// min principle
		private double[] curva = null;		// curv tensor matrix {{a, b}, {b, c}}
		private double[] curvb = null;
		private double[] curvc = null;

		private double[] dual_curva = null;	// dual curv tensor matrix {{a, b}, {b, c}}
		private double[] dual_curvb = null;
		private double[] dual_curvc = null;


		private Vector3d[] pDir1 = null;
		private Vector3d[] pDir2 = null;


		private Vector3d[] localvu = null;
		private Vector3d[] localvv = null;
		
		
		private void LocateMemory()
		{
			int nv = mesh.VertexCount, nf = mesh.FaceCount;

			// Initialize curvature tensor matrix elements.
			this.curva = new double[nv];
			this.curvb = new double[nv];
			this.curvc = new double[nv];
			this.dual_curva = new double[nf];
			this.dual_curvb = new double[nf];
			this.dual_curvc = new double[nf];

			// Initialize k1, k2, e1, e2 for each vertex.
			this.pDir1 = new Vector3d[nv];
			this.pDir2 = new Vector3d[nv];
			this.curv2 = new double[nv];
			this.curv1 = new double[nv];

			// Initialize k1, k2, e1, e2 for each face.
			this.fpDir1 = new Vector3d[nf];
			this.fpDir2 = new Vector3d[nf];
			this.fcurv2 = new double[nf];
			this.fcurv1 = new double[nf];
		}
		// Compute the principle curvatures and directions, i.e. k1, k2, e1, e2.
		public void ObtainPerVertexPCurAndDirs()
		{
			// initialize areas.
			ComputePointAreas();
			int nf = mesh.FaceCount, nv = mesh.VertexCount;
			// Set up an initial coordinate system per vertex
			localvu = new Vector3d[nv];
			localvv = new Vector3d[nv];
			for (int i = 0; i < nf; i++)
			{
				int b = i * 3;
				int c0 = mesh.FaceIndex[b];
				int c1 = mesh.FaceIndex[b + 1];
				int c2 = mesh.FaceIndex[b + 2];
				Vector3d v0 = new Vector3d(mesh.VertexPos, c0 * 3);
				Vector3d v1 = new Vector3d(mesh.VertexPos, c1 * 3);
				Vector3d v2 = new Vector3d(mesh.VertexPos, c2 * 3);
				localvu[c0] = v1 - v0;
				localvu[c1] = v2 - v1;
				localvu[c2] = v0 - v2;
			}
			for (int i = 0; i < nv; i++)
			{
				Vector3d normali = new Vector3d(mesh.VertexNormal, i * 3);
				localvu[i] = (localvu[i].Cross(normali)).Normalize();
				localvv[i] = normali.Cross(localvu[i]);
			}
			// Compute curvature per-face
			for (int i = 0; i < nf; i++) {
				// Edges
				int c = i*3;
				int d0 = mesh.FaceIndex[c] * 3;
				int d1 = mesh.FaceIndex[c + 1] * 3;
				int d2 = mesh.FaceIndex[c + 2] * 3;
				Vector3d[] e = new Vector3d[3] {
					new Vector3d(mesh.VertexPos, d2) - new Vector3d(mesh.VertexPos, d1),
					new Vector3d(mesh.VertexPos, d0) - new Vector3d(mesh.VertexPos, d2),
					new Vector3d(mesh.VertexPos, d1) - new Vector3d(mesh.VertexPos, d0)
				};
				Vector3d[] normals = new Vector3d[3] { 
					new Vector3d(mesh.VertexNormal, d0),
					new Vector3d(mesh.VertexNormal, d1),
					new Vector3d(mesh.VertexNormal, d2)
				};

				// N-T-B coordinate system per face
				Vector3d t = e[0].Normalize();
				Vector3d n = e[0].Cross(e[1]);
				Vector3d b = n.Cross(t).Normalize();

				// Estimate curvature based on variation of normals
				// along edges
				double[] m = new double[3] { 0, 0, 0 };
				double[][] w = new double[3][];
				for (int j = 0; j < 3; ++j) w[j] = new double[3];
				for (int j = 0; j < 3; ++j) {
					double u = e[j].Dot(t);
					double v = e[j].Dot(b);
					w[0][0] += u*u;
					w[0][1] += u*v;
					w[2][2] += v*v;
					Vector3d dn = normals[(j+2)%3] - normals[(j+1)%3];
					double dnu = dn.Dot(t);
					double dnv = dn.Dot(b);
					m[0] += dnu*u;
					m[1] += dnu*v + dnv*u;
					m[2] += dnv*v;
				}
				w[1][1] = w[0][0] + w[2][2];
				w[1][2] = w[0][1];

				// Least squares solution
				double[] diag = new double[3];
				if (!ldltdc(w, diag, 3)) {
					System.Console.WriteLine("ldltdc failed!\n");
					continue;
				}
				ldltsl(w, diag, m, ref m, 3);

				// Push it back out to the vertices
				for (int j = 0; j < 3; j++) {
					int vj = mesh.FaceIndex[c+j];
					double c1 = 0, c12 = 0, c2 = 0;
					proj_curv(t, b, m[0], m[1], m[2],
						  localvu[vj], localvv[vj], 
						  ref c1, ref c12, ref c2);
					double wt = cornerareas[i][j] / pointareas[vj];
					curva[vj] += wt * c1;
					curvb[vj] += wt * c12;
					curvc[vj] += wt * c2;
				}
			}

			// Compute principal directions and curvatures at each vertex
			for (int i = 0; i < nv; i++) {
				Vector3d ni = new Vector3d(mesh.VertexNormal, i * 3);
				DiagonalizeCurvature(localvu[i], localvv[i],
						 curva[i], curvb[i], curvc[i],
						 ni, ref pDir1[i], ref pDir2[i],
						 ref curv1[i], ref curv2[i]);
			}
			System.Console.WriteLine("Done.\n");
		}

		// Compute the principle curvatures and directions for a specific vertex
		// this is due to the sence that the vertex position may be changed on 
		// the fly for some geometry filtering algorithms.
		public void ObtainCurvDirAt(int index, ref Vector3d e1, ref Vector3d e2, ref double k1, ref double k2)
		{
			// Initialize curvature tensor matrix elements.
			double aa = 0;
			double bb = 0;
			double cc = 0;

			Vector3d ni = new Vector3d(mesh.VertexNormal, index * 3);
			Vector3d pdir1 = localvu[index];
			Vector3d pdir2 = localvv[index];

			// Compute curvature per-face
			foreach (int f in mesh.AdjVF[index])
			{
				// Edges
				int c = f * 3;
				int d0 = mesh.FaceIndex[c] * 3;
				int d1 = mesh.FaceIndex[c + 1] * 3;
				int d2 = mesh.FaceIndex[c + 2] * 3;
				Vector3d[] e = new Vector3d[3] {
					new Vector3d(mesh.VertexPos, d2) - new Vector3d(mesh.VertexPos, d1),
					new Vector3d(mesh.VertexPos, d0) - new Vector3d(mesh.VertexPos, d2),
					new Vector3d(mesh.VertexPos, d1) - new Vector3d(mesh.VertexPos, d0)
				};
				Vector3d[] normals = new Vector3d[3] { 
					new Vector3d(mesh.VertexNormal, d0),
					new Vector3d(mesh.VertexNormal, d1),
					new Vector3d(mesh.VertexNormal, d2)
				};

				// N-T-B coordinate system per face
				Vector3d t = e[0].Normalize();
				Vector3d n = e[0].Cross(e[1]);
				Vector3d b = n.Cross(t).Normalize();

				// Estimate curvature based on variation of normals
				// along edges
				double[] m = new double[3] { 0, 0, 0 };
				double[][] w = new double[3][];
				for (int jj = 0; jj < 3; ++jj) w[jj] = new double[3];
				for (int jj = 0; jj < 3; ++jj)
				{
					double u = e[jj].Dot(t);
					double v = e[jj].Dot(b);
					w[0][0] += u * u;
					w[0][1] += u * v;
					w[2][2] += v * v;
					Vector3d dn = normals[(jj + 2) % 3] - normals[(jj + 1) % 3];
					double dnu = dn.Dot(t);
					double dnv = dn.Dot(b);
					m[0] += dnu * u;
					m[1] += dnu * v + dnv * u;
					m[2] += dnv * v;
				}
				w[1][1] = w[0][0] + w[2][2];
				w[1][2] = w[0][1];

				// Least squares solution
				double[] diag = new double[3];
				if (!ldltdc(w, diag, 3))
				{
					System.Console.WriteLine("ldltdc failed!\n");
					continue;
				}
				ldltsl(w, diag, m, ref m, 3);

				// Push it back out to the vertices
				for (int jj = 0; jj < 3; jj++)
				{
					int pj = mesh.FaceIndex[c + jj];

					if (pj != index) continue;
					
					double c1 = 0, c12 = 0, c2 = 0;
					proj_curv(t, b, m[0], m[1], m[2],
						  pdir1, pdir2,
						  ref c1, ref c12, ref c2);
					double wt = cornerareas[f][jj] / pointareas[pj];
					aa += wt * c1;
					bb += wt * c12;
					cc += wt * c2;
				}
			}

			// Compute principal directions and curvatures at each vertex
			DiagonalizeCurvature(pdir1, pdir2,
					 aa, bb, cc,
					 ni, ref e1, ref e2, ref k1, ref k2
					 );

			System.Console.WriteLine("Done.\n");
		}



		private Vector3d[] fpDir1 = null;
		private Vector3d[] fpDir2 = null;
		private double[] fcurv1 = null;
		private double[] fcurv2 = null;

		// Compute the principle curvatures and directions for a specific face.
		public void ObtainPerFacePCurAndDirs()
		{
			// Compute curvature per-face
			int nf = mesh.FaceCount, nv = mesh.VertexCount;
			for (int i = 0; i < nf; i++)
			{
				// Edges
				int c = i * 3;
				int d0 = mesh.FaceIndex[c] * 3;
				int d1 = mesh.FaceIndex[c + 1] * 3;
				int d2 = mesh.FaceIndex[c + 2] * 3;
				Vector3d[] e = new Vector3d[3] {
					new Vector3d(mesh.VertexPos, d2) - new Vector3d(mesh.VertexPos, d1),
					new Vector3d(mesh.VertexPos, d0) - new Vector3d(mesh.VertexPos, d2),
					new Vector3d(mesh.VertexPos, d1) - new Vector3d(mesh.VertexPos, d0)
				};
				Vector3d[] normals = new Vector3d[3] { 
					new Vector3d(mesh.VertexNormal, d0),
					new Vector3d(mesh.VertexNormal, d1),
					new Vector3d(mesh.VertexNormal, d2)
				};

				// N-T-B coordinate system per face
				Vector3d t = e[0].Normalize();
				Vector3d n = e[0].Cross(e[1]);
				Vector3d b = n.Cross(t).Normalize();

				// Estimate curvature based on variation of normals
				// along edges
				double[] m = new double[3] { 0, 0, 0 };
				double[][] w = new double[3][];
				for (int j = 0; j < 3; ++j) w[j] = new double[3];
				for (int j = 0; j < 3; ++j)
				{
					double u = e[j].Dot(t);
					double v = e[j].Dot(b);
					w[0][0] += u * u;
					w[0][1] += u * v;
					w[2][2] += v * v;
					Vector3d dn = normals[(j + 2) % 3] - normals[(j + 1) % 3];
					double dnu = dn.Dot(t);
					double dnv = dn.Dot(b);
					m[0] += dnu * u;
					m[1] += dnu * v + dnv * u;
					m[2] += dnv * v;
				}
				w[1][1] = w[0][0] + w[2][2];
				w[1][2] = w[0][1];

				// Least squares solution
				double[] diag = new double[3];
				if (!ldltdc(w, diag, 3))
				{
					System.Console.WriteLine("ldltdc failed!\n");
					continue;
				}
				ldltsl(w, diag, m, ref m, 3);

				// compute the principle curvatures on face i.
				Vector3d fn = new Vector3d(mesh.FaceNormal, i * 3);
				DiagonalizeCurvature(t, b,
						 m[0], m[1], m[2],
						 fn, ref fpDir1[i], ref fpDir2[i],
						 ref fcurv1[i], ref fcurv2[i]);

			}
			System.Console.WriteLine("Done.\n");
		}
		public void ObtainDualPCurAndDirs()
		{
			ComputeDualPointAreas();
			// Compute curvature per-face
			int nf = mesh.FaceCount, nv = mesh.VertexCount;
			for (int i = 0; i < nf; i++)
			{
				// compute local uv coordinates for face i
				int cc = i * 3;
				int c0 = mesh.FaceIndex[cc] * 3;
				int c1 = mesh.FaceIndex[cc + 1] * 3;
				int c2 = mesh.FaceIndex[cc + 2] * 3;
				Vector3d[] eg = new Vector3d[3] {
					new Vector3d(mesh.VertexPos, c2) - new Vector3d(mesh.VertexPos, c1),
					new Vector3d(mesh.VertexPos, c0) - new Vector3d(mesh.VertexPos, c2),
					new Vector3d(mesh.VertexPos, c1) - new Vector3d(mesh.VertexPos, c0)
				};
				// local coordinate system per face
				Vector3d uu = eg[0].Normalize();
				Vector3d nn = new Vector3d(mesh.FaceNormal, cc);
				Vector3d vv = nn.Cross(uu).Normalize();

				// Edges
				int[] c = new int[3];
				for (int j = 0; j < 3; ++j)
				{
					c[j] = -1;
				}
				for (int j = 0; j < mesh.AdjFF[i].Length; ++j)
				{
					c[j] = mesh.AdjFF[i][j]*3;
				}
				for (int j = 0; j < 3; ++j)
				{
					int d0 = cc;
					int d1 = c[j];
					int d2 = c[(j + 1) % 3];
					if (d1 < 0 || d2 < 0) continue;

					Vector3d[] e = new Vector3d[3] {
						new Vector3d(mesh.DualVertexPos, d2) - 
						new Vector3d(mesh.DualVertexPos, d1),
						new Vector3d(mesh.DualVertexPos, d0) - 
						new Vector3d(mesh.DualVertexPos, d2),
						new Vector3d(mesh.DualVertexPos, d1) - 
						new Vector3d(mesh.DualVertexPos, d0)
					};
					Vector3d[] normals = new Vector3d[3] { 
						new Vector3d(mesh.FaceNormal, d0),
						new Vector3d(mesh.FaceNormal, d1),
						new Vector3d(mesh.FaceNormal, d2)
					};

					// N-T-B coordinate system per face
					Vector3d t = e[0].Normalize();
					Vector3d n = e[0].Cross(e[1]).Normalize();
					Vector3d b = n.Cross(t).Normalize();

					// Estimate curvature based on variation of normals
					// along edges
					double[] m = new double[3] { 0, 0, 0 };
					double[][] w = new double[3][];
					for (int k = 0; k < 3; ++k) w[k] = new double[3];
					for (int k = 0; k < 3; ++k)
					{
						double u = e[k].Dot(t);
						double v = e[k].Dot(b);
						w[0][0] += u * u;
						w[0][1] += u * v;
						w[2][2] += v * v;
						Vector3d dn = normals[(k + 2) % 3] - normals[(k + 1) % 3];
						double dnu = dn.Dot(t);
						double dnv = dn.Dot(b);
						m[0] += dnu * u;
						m[1] += dnu * v + dnv * u;
						m[2] += dnv * v;
					}
					w[1][1] = w[0][0] + w[2][2];
					w[1][2] = w[0][1];

					// Least squares solution
					double[] diag = new double[3];
					if (!ldltdc(w, diag, 3))
					{
						System.Console.WriteLine("ldltdc failed!\n");
						continue;
					}
					ldltsl(w, diag, m, ref m, 3);

					// Push it back out to the dual vertices
					double g1 = 0, g12 = 0, g2 = 0;
					proj_curv(t, b, m[0], m[1], m[2],
						  uu, vv, ref g1, ref g12, ref g2);
					double wt = dual_cornerareas[i][j] / dual_pointareas[i];
					dual_curva[i] += wt * g1;
					dual_curvb[i] += wt * g12;
					dual_curvc[i] += wt * g2;
				}

				// Compute principal directions and curvatures at each dual vertex
				DiagonalizeCurvature(uu, vv,
						 dual_curva[i], dual_curvb[i], dual_curvc[i],
						 nn, ref fpDir1[i], ref fpDir2[i],
						 ref fcurv1[i], ref fcurv2[i]);
			}
			System.Console.WriteLine("Done.\n");
		}

		// From trimesh2 source code..
		// Rotate a coordinate system to be perpendicular to the given normal
		private void RotateCoordSys(Vector3d old_u, Vector3d old_v,
					  Vector3d new_norm,
					  ref Vector3d new_u, ref Vector3d new_v)
		{
			new_u = old_u;
			new_v = old_v;
			Vector3d old_norm = old_u.Cross(old_v);
			double ndot = old_norm.Dot(new_norm);
			if ((ndot <= -1.0f)) {
				new_u = new Vector3d() - new_u;
				new_v = new Vector3d() - new_v;
				return;
			}
			Vector3d perp_old = new_norm - ndot * old_norm;
			Vector3d dperp = 1.0f / (1 + ndot) * (old_norm + new_norm);
			new_u -= dperp * (new_u.Dot(perp_old));
			new_v -= dperp * (new_v.Dot(perp_old));
		}
		// Given a curvature tensor, find principal directions and curvatures
		// Makes sure that pdir1 and pdir2 are perpendicular to normal
		private void DiagonalizeCurvature(Vector3d old_u, Vector3d old_v,
					  double ku, double kuv, double kv,
					  Vector3d new_norm,
					  ref Vector3d pdir1, ref Vector3d pdir2, ref double k1, ref double k2)
		{
			Vector3d r_old_u = new Vector3d();
			Vector3d r_old_v = new Vector3d();

			RotateCoordSys(old_u, old_v, new_norm, ref r_old_u, ref r_old_v);

			double c = 1, s = 0, tt = 0;
			if (kuv != 0.0) {
				// Jacobi rotation to diagonalize
				double h = 0.5 * (kv - ku) / kuv;
				tt = (h < 0) ?
					1.0 / (h - Math.Sqrt(1.0 + h*h)) :
					1.0 / (h + Math.Sqrt(1.0 + h*h));
				c = 1.0 / Math.Sqrt(1.0 + tt*tt);
				s = tt * c;
			}

			k1 = ku - tt * kuv;
			k2 = kv + tt * kuv;

			if (Math.Abs(k1) >= Math.Abs(k2)) {
				pdir1 = c*r_old_u - s*r_old_v;
			} else {
				double tp = k1;
				k1 = k2;
				k2 = tp;
				pdir1 = s*r_old_u + c*r_old_v;
			}
			pdir2 = new_norm.Cross(pdir1);
		}
		// Reproject a curvature tensor from the basis spanned by old_u and old_v
		// (which are assumed to be unit-length and perpendicular) to the
		// new_u, new_v basis.
		private void proj_curv(Vector3d old_u, Vector3d old_v,
				   double old_ku, double old_kuv, double old_kv,
				   Vector3d new_u, Vector3d new_v,
				   ref double new_ku, ref double new_kuv, ref double new_kv)
		{
			Vector3d r_new_u = new Vector3d();
			Vector3d r_new_v = new Vector3d();

			RotateCoordSys(new_u, new_v, old_u.Cross(old_v), ref r_new_u, ref r_new_v);

			double u1 = r_new_u.Dot(old_u);
			double v1 = r_new_u.Dot(old_v);
			double u2 = r_new_v.Dot(old_u);
			double v2 = r_new_v.Dot(old_v);
			new_ku  = old_ku * u1*u1 + old_kuv * (2.0 * u1*v1) + old_kv * v1*v1;
			new_kuv = old_ku * u1*u2 + old_kuv * (u1*v2 + u2*v1) + old_kv * v1*v2;
			new_kv  = old_ku * u2*u2 + old_kuv * (2.0 * u2*v2) + old_kv * v2*v2;
		}
		// Perform LDL^T decomposition of a symmetric positive definite matrix.
		// Like Cholesky, but no square roots.  Overwrites lower triangle of matrix.
		private bool ldltdc(double[][] A, double[] rdiag, int N)
		{
			double[] v = new double[N-1];
			for (int i = 0; i < N; i++) {
				for (int k = 0; k < i; k++)
					v[k] = A[i][k] * rdiag[k];
				for (int j = i; j < N; j++) {
					double sum = A[i][j];
					for (int k = 0; k < i; k++)
						sum -= v[k] * A[j][k];
					if (i == j) {
						if ((sum <= 0))
							return false;
						rdiag[i] = 1.0 / sum;
					} else {
						A[j][i] = sum;
					}
				}
			}
			return true;
		}

		// Solve Ax=B after ldltdc
		private void ldltsl(double[][] A, double[] rdiag, double[] B, ref double[] x, int N)
		{
			int i;
			for (i = 0; i < N; i++) {
				double sum = B[i];
				for (int k = 0; k < i; k++)
					sum -= A[i][k] * x[k];
				x[i] = sum * rdiag[i];
			}
			for (i = N - 1; i >= 0; i--) {
				double sum = 0;
				for (int k = i + 1; k < N; k++)
					sum += A[k][i] * x[k];
				x[i] -= sum * rdiag[i];
			}
		}


		private Vector3d[] cornerareas = null;
		private double[] pointareas = null;
		private Vector3d[] dual_cornerareas = null;
		private double[] dual_pointareas = null;
		private void ComputePointAreas()
		{
			System.Console.WriteLine("Computing point areas... ");

			int nf = mesh.FaceCount, nv = mesh.VertexCount;
			
			cornerareas = new Vector3d[nf];
			pointareas = new double[nv];
			
			for (int i = 0; i < nf; i++) {
				cornerareas[i] = new Vector3d();
				// Edges
				int c = i*3;
				int c0 = mesh.FaceIndex[c]*3;
				int c1 = mesh.FaceIndex[c+1]*3;
				int c2 = mesh.FaceIndex[c+2]*3;
				Vector3d[] e = new Vector3d[3] {
					new Vector3d(mesh.VertexPos, c2) - new Vector3d(mesh.VertexPos, c1),
					new Vector3d(mesh.VertexPos, c0) - new Vector3d(mesh.VertexPos, c2),
					new Vector3d(mesh.VertexPos, c1) - new Vector3d(mesh.VertexPos, c0)
				};

				// Compute corner weights
				double area = 0.5 * (e[0].Cross(e[1])).Length();
				double[] l2 = new double[3] { e[0].Dot(e[0]), e[1].Dot(e[1]), e[2].Dot(e[2]) };
				double[] ew = new double[3] { 
					l2[0] * (l2[1] + l2[2] - l2[0]),
					l2[1] * (l2[2] + l2[0] - l2[1]),
					l2[2] * (l2[0] + l2[1] - l2[2]) 
				};
				if (ew[0] <= 0.0) {
					cornerareas[i][1] = -0.25 * l2[2] * area /
								(e[0].Dot(e[2]));
					cornerareas[i][2] = -0.25 * l2[1] * area /
								(e[0].Dot(e[1]));
					cornerareas[i][0] = area - cornerareas[i][1] -
								cornerareas[i][2];
				} else if (ew[1] <= 0.0) {
					cornerareas[i][2] = -0.25 * l2[0] * area /
								(e[1].Dot(e[0]));
					cornerareas[i][0] = -0.25 * l2[2] * area /
								(e[1].Dot(e[2]));
					cornerareas[i][1] = area - cornerareas[i][2] -
								cornerareas[i][0];
				} else if (ew[2] <= 0.0) {
					cornerareas[i][0] = -0.25 * l2[1] * area /
								(e[2].Dot(e[1]));
					cornerareas[i][1] = -0.25 * l2[0] * area /
								(e[2].Dot(e[0]));
					cornerareas[i][2] = area - cornerareas[i][0] -
								cornerareas[i][1];
				} else {
					double ewscale = 0.5 * area / (ew[0] + ew[1] + ew[2]);
					for (int j = 0; j < 3; j++)
						cornerareas[i][j] = ewscale * (ew[(j+1)%3] +
										   ew[(j+2)%3]);
				}
				pointareas[mesh.FaceIndex[c]] += cornerareas[i][0];
				pointareas[mesh.FaceIndex[c+1]] += cornerareas[i][1];
				pointareas[mesh.FaceIndex[c+2]] += cornerareas[i][2];
			}
			// Done!.
		}
		private void ComputeDualPointAreas() // 'Y' shape
		{
			System.Console.WriteLine("Computing dual point areas... ");

			int nf = mesh.FaceCount, nv = mesh.VertexCount;

			dual_cornerareas = new Vector3d[nf];
			dual_pointareas = new double[nf];

			for (int i = 0; i < nf; i++)
			{
				dual_cornerareas[i] = new Vector3d();
				// Edges
				int[] c = new int[3];
				for (int j = 0; j < 3; ++j)
				{
					c[j] = -1;
				}
				for (int j = 0; j < mesh.AdjFF[i].Length; ++j)
				{
					c[j] = mesh.AdjFF[i][j]*3;
				}
				for (int j = 0; j < 3; ++j)
				{
					int c0 = i*3;
					int c1 = c[j];
					int c2 = c[(j + 1) % 3];
					if (c1 < 0 || c2 < 0) continue;
					Vector3d[] e = new Vector3d[3] {
						new Vector3d(mesh.DualVertexPos, c2) - 
						new Vector3d(mesh.DualVertexPos, c1),
						new Vector3d(mesh.DualVertexPos, c0) - 
						new Vector3d(mesh.DualVertexPos, c2),
						new Vector3d(mesh.DualVertexPos, c1) - 
						new Vector3d(mesh.DualVertexPos, c0)
					};

					// Compute corner weights
					double area = 0.5 * (e[0].Cross(e[1])).Length();
					double[] l2 = new double[3] { e[0].Dot(e[0]), e[1].Dot(e[1]), e[2].Dot(e[2]) };
					double[] ew = new double[3] { 
						l2[0] * (l2[1] + l2[2] - l2[0]),
						l2[1] * (l2[2] + l2[0] - l2[1]),
						l2[2] * (l2[0] + l2[1] - l2[2]) 
					};
					if (ew[0] <= 0.0)
					{
						double a1 = -0.25 * l2[2] * area /
									(e[0].Dot(e[2]));
						double a2 = -0.25 * l2[1] * area /
									(e[0].Dot(e[1]));
						dual_cornerareas[i][j] = area - a1 - a2;
					}
					else if (ew[1] <= 0.0)
					{
						dual_cornerareas[i][j] = -0.25 * l2[2] * area /
									(e[1].Dot(e[2]));
					}
					else if (ew[2] <= 0.0)
					{
						dual_cornerareas[i][j] = -0.25 * l2[1] * area /
									(e[2].Dot(e[1]));
					}
					else
					{
						double ewscale = 0.5 * area / (ew[0] + ew[1] + ew[2]);
						dual_cornerareas[i][j] = ewscale * (ew[1] + ew[2]);
					}
					dual_pointareas[i] += dual_cornerareas[i][j];
				}
			}
			// Done!.
		}

		#region interfaces
		public double[] PMaxCurv
		{
			get { return curv1; }
		}
		public double[] PMinCurv
		{
			get { return curv2; }
		}
		public Vector3d[] PDir1
		{
			get { return pDir1; }
		}
		public Vector3d[] PDir2
		{
			get { return pDir2; }
		}
		// for facets
		public double[] PerFacePMaxCurv
		{
			get { return fcurv1; }
		}
		public double[] PerFacePMinCurv
		{
			get { return fcurv2; }
		}
		public Vector3d[] PerFacePDir1
		{
			get { return fpDir1; }
		}
		public Vector3d[] PerFacePDir2
		{
			get { return fpDir2; }
		}
		public void Process(Mesh m)
		{
			this.mesh = m;

			this.LocateMemory();
			this.ObtainPerVertexPCurAndDirs();

		}
		#endregion
		#region IDisposable Members
		public Curvatures()
		{
			
		}
		~Curvatures()
		{
			Dispose();
		}
		public void Dispose()
		{
			
		}
		#endregion
	}
}
