using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using Emgu.CV;
using Emgu.CV.Structure;
using Emgu.Util;
using Emgu.CV.UI;

using MyGeometry;
using CySoft.Collections.Geometry;			// for kd tree search


namespace IBrushes
{
	public class ManifoldPoint : PriorityQueueElement // - priority queue for geodesic distance computation
	{
		private int index;
		private double distance;						// for geodesic computation
		public double Il;								// Il = Wl*Sl; -- Wl is the reflectance component, Sl is the shading component
		public double Wl;								// the reflectance component of Il
		public double Sl;								// the shading component of Il which equals to Il/Wl;
		public double Kw;								// the weathering parameter
		public Vector3d Pos;							// (x, y, z)  <=> (l, a , b)
		public double WeatheringDegree;					// weathering degree
		public double OldWeatheringDegree;				// 
		private List<ManifoldPoint> neighbors;			// neighbor points
		public ManifoldPoint Parent;					// for dijkstra algorithm
		public ManifoldPoint(int idx, Vector3d p)
		{
			this.index = idx;
			this.distance = double.MaxValue;			// init to maximum
			this.Il = p.x;
			this.Pos = new Vector3d(0,p.y,p.z);			// not copying the L part (Lab)
			this.neighbors = new List<ManifoldPoint>();
		}

		// public entries
		public int Index { get { return this.index; } }
		public double Distance { get { return this.distance; } set { this.distance = value; } }
		public List<ManifoldPoint> Neighbors { get { return neighbors; } }
		
		public void AddNeighbor(ManifoldPoint adj)
		{
			if (this.neighbors == null) this.neighbors = new List<ManifoldPoint>();
			if (!this.neighbors.Contains(adj))
				this.neighbors.Add(adj);
		}
		public void ClearAllNeighbors()
		{
			if (this.neighbors != null) this.neighbors.Clear();
		}
		public void RemoveNeighbor(ManifoldPoint adj)
		{
			if (this.neighbors != null) 
				this.neighbors.Remove(adj);
		}
		public double[] ToArray()
		{
			return Pos.ToArray();
		}

		#region PriorityQueueElement Members
		private int pqIndex = -1;
		public int PQIndex
		{
			get
			{
				return pqIndex;
			}
			set
			{
				pqIndex = value;
			}
		}
		public int CompareTo(PriorityQueueElement other)
		{
			ManifoldPoint rec = other as ManifoldPoint;
			if (this.distance < rec.distance) return -1;
			if (this.distance > rec.distance) return 1;
			return 0;
		}
		#endregion
	}

	public class ImageWeathering
	{
		// members
		private const int kNeighbors = 8;							// num of neighbors to build the manifold
		private Image<Lab, byte> labImage;							// image to process, in L*a*b space
		private ManifoldPoint[,] manifoldPoints;					// holds manifold points
		private List<ManifoldPoint> leastWeatheredPoints;			// holds least weathered points
		private List<ManifoldPoint> mostWeatheredPoints;			// holds most weathered points
		private List<ManifoldPoint> nonWeatheredPoints;				// holds non-weathered points
		private KDTree<ManifoldPoint> kdTree;						// storing manifold points for fast query
		private double[,][] manifoldGeodesic;
		private byte[,] userMask;

		public ImageWeathering(Image<Lab, byte> image)
		{
			this.labImage = image;
		}


		public void InitUserMask(List<BrushStroke> strokes)
		{
			int w = this.labImage.Width, h = this.labImage.Height;
			this.userMask = new byte[w, h];
			foreach (BrushStroke stroke in strokes)
			{
				int size = (int)stroke.BrushPoints[0].size;
				foreach (BrushPoint point in stroke.BrushPoints)
				{
					for (int i = -size; i < size; ++i)
					{
						for (int j = -size; j < size; ++j)
						{
							int I = (int)point.x + i;
							int J = (int)point.y + j;
							if (I >= 0 && I < w && J >= 0 && J < h && 
								(new Vector2d(I,J)-new Vector2d(point.x,point.y)).Length() <= size)
							{
								this.userMask[I, J] = stroke.Type;
							}
						}
					}
				}
			}
		}
		public void Process()
		{
			Console.WriteLine("init manifold ...");
			this.InitAppearanceManifold();
			Console.WriteLine("build neighborhood ...");
			this.BuildManifoldNeighborHood();
			Console.WriteLine("convexilization ...");
			this.ConvexilizePoints();


			Console.WriteLine("improve manifold ...");
			// iteratively construct a stable appearance manifold
			int itr = 0, maxitr = 5;
			do
			{
				this.ComputeWeatheringDegreeMap();
				this.UpdateWeatheringComponent();
				this.BuildManifoldNeighborHood();
			}
			while (itr++ < maxitr);

			// refine the weathering degree map
			Console.WriteLine("compute shading components ...");
			this.ComputeShadingComponent();

			// before edit, compute all path geodesic distances
			Console.WriteLine("computing manifold all path geodesics ...");
			this.ComputeManifoldGeodesic();

			// compute weathering degrees
			Console.WriteLine("computing weathering parameter ...");
			this.ComputeWeatheringParameter();
		}
		public void Weathering(double t) // t -> time
		{
			// divid into bins by weathering degrees -- for fast query
			List<ManifoldPoint>[] w_buckets = this.PackPointsIntoBins(1000);
			double intv = 1.0 / 1000;

			// by editing t, we edit the weathering degree (0,1) is deweathering, (1,infinity) is weathering
			int w = this.labImage.Width, h = this.labImage.Height;
			double[,] new_wl = new double[w, h];
			for (int i = 0; i < w; ++i)
			{
				for (int j = 0; j < h; ++j)
				{
					ManifoldPoint pt = this.manifoldPoints[i, j];

					// compute the new weathering degree d'
					double ds = pt.WeatheringDegree;
					double dt = 1 - Math.Exp(-pt.Kw * t);
					List<ManifoldPoint> srcpoints = this.FindPointsWithDegree(w_buckets[(int)(ds / intv)], ds);
					List<ManifoldPoint> tarpoints = this.FindPointsWithDegree(w_buckets[(int)(dt / intv)], dt);

					// compute src and tar geometric center
					Vector3d srcvc = new Vector3d();
					foreach (ManifoldPoint mp in srcpoints)
					{
						srcvc += mp.Pos;
					}
					srcvc *= (1.0 / srcpoints.Count);
					Vector3d tarvc = new Vector3d();
					foreach (ManifoldPoint mp in tarpoints)
					{
						tarvc += mp.Pos;
					}
					tarvc *= (1.0 / tarpoints.Count);

					// find points closest to geometric centers
					ManifoldPoint svc = null, tvc = null;
					double mindist = double.MaxValue;
					foreach (ManifoldPoint mp in srcpoints)
					{
						double d = (mp.Pos - srcvc).Length();
						if (d < mindist)
						{
							mindist = d;
							svc = mp;
						}
					}
					mindist = double.MaxValue;
					foreach (ManifoldPoint mp in tarpoints)
					{
						double d = (mp.Pos - tarvc).Length();
						if (d < mindist)
						{
							mindist = d;
							tvc = mp;
						}
					}


					int I = svc.Index / h, J = svc.Index % h;
					double RadSv1 = double.MinValue;
					foreach (ManifoldPoint mp in srcpoints)
					{
						double d = this.manifoldGeodesic[I, J][mp.Index];
						if (d > RadSv1)
						{
							RadSv1 = d;
						}
					}
					double RbrVs1 = this.manifoldGeodesic[I, J][pt.Index] / RadSv1;		// see paper - the distance between svc and the src point

					// find the target point with degree d and assign the weathering component
					I = tvc.Index / h; J = tvc.Index % h;
					double RadTa1 = double.MinValue;
					double[] geo = this.manifoldGeodesic[I, J];
					foreach (ManifoldPoint mp in tarpoints)
					{
						double d = geo[mp.Index];
						if (d > RadTa1)
						{
							RadTa1 = d;
						}
					}

					// find the target point with degree d and assign the weathering component
					mindist = double.MaxValue;
					ManifoldPoint target = null;
					foreach (ManifoldPoint mp in tarpoints)
					{
						double Rbrv = geo[mp.Index] / RadTa1;
						double geosv = this.manifoldGeodesic[i, j][mp.Index];
						double d = geosv * (Rbrv - RbrVs1);
						if (d < mindist)
						{
							mindist = d;
							target = mp;
						}
					}

					// apply the effect (combine with shading component)
					this.labImage[i, j] = new Lab(target.Wl * pt.Sl, target.Pos.y, target.Pos.z);

				}
			}

		}
		

		
		private void ComputeWeatheringParameter()
		{
			// D'(x,y,t) = 1 - e^(-K(x,y)*t), t is in [0,infinity] and D(x,y,1) = D(x,y)
			// first compute K(x,y)
			foreach (ManifoldPoint pt in this.manifoldPoints)
			{
				pt.Kw = -Math.Log(1 - pt.WeatheringDegree);
			}
		}
		private void InitAppearanceManifold()
		{
			int w = this.labImage.Width, h = this.labImage.Height;

			this.manifoldPoints = new ManifoldPoint[w, h];
			this.leastWeatheredPoints = new List<ManifoldPoint>();
			this.mostWeatheredPoints = new List<ManifoldPoint>();
			this.nonWeatheredPoints = new List<ManifoldPoint>();

			this.kdTree = new KDTree<ManifoldPoint>(3);

			// create manifold points
			for (int i = 0; i < w; ++i)
			{
				for (int j = 0; j < h; ++j)
				{
					Lab lab = this.labImage[j, i];
					ManifoldPoint mp = new ManifoldPoint(i * h + j,
						new Vector3d(lab.X, lab.Y, lab.Z)
						);

					this.manifoldPoints[i, j] = mp;

					switch (this.userMask[i, j])
					{
						case 1:		// least weathered -- left mouse button
							this.leastWeatheredPoints.Add(mp);
							break;
						case 2:		// least weathered -- left mouse button
							this.mostWeatheredPoints.Add(mp);
							break;
						case 3:		// unused pointes (e.g., shadows, etc)
							this.nonWeatheredPoints.Add(mp);
							break;
					}

					// fill tree
					kdTree.Add(mp, mp.Pos.x, mp.Pos.y, mp.Pos.z);

				}
			}
		}
		private void BuildManifoldNeighborHood()
		{
			// build neighborhood
			foreach (ManifoldPoint mp in this.manifoldPoints) mp.ClearAllNeighbors();	// clear all previous neighbors
			foreach (ManifoldPoint mp in this.manifoldPoints)
			{
				var neighbors = kdTree.NearestNeighbors(SquareEuclideanDistance.Instance, kNeighbors, mp.ToArray());
				var dict = new HashSet<ManifoldPoint>(neighbors.Select(n => n.Value));
				dict.ToList().ForEach(p => mp.AddNeighbor(p));
			}
			
			// do filtering
			double avgdist = 0;			// compute average distance between connected points
			int ncount = 0;
			foreach (ManifoldPoint mp in this.manifoldPoints)
			{
				foreach (ManifoldPoint adjmp in mp.Neighbors)
				{
					avgdist += (adjmp.Pos - mp.Pos).Length();
				}
				ncount += mp.Neighbors.Count;
			}
			avgdist /= ncount;
			double thres = 1.5 * avgdist;
			foreach (ManifoldPoint mp in this.manifoldPoints)
			{
				List<ManifoldPoint> points2remove = new List<ManifoldPoint>();
				foreach (ManifoldPoint adjmp in mp.Neighbors)
				{
					double dis = (adjmp.Pos - mp.Pos).Length();
					if (dis > thres)
						points2remove.Add(adjmp);
				}
				points2remove.ForEach(p => mp.Neighbors.Remove(p));
			}
		}
		private void ComputeWeatheringDegreeMap()
		{
			if (this.leastWeatheredPoints == null || this.mostWeatheredPoints == null) return;

			// compute distances from src/tar points to all other points
			double[][] src_dist = new double[this.leastWeatheredPoints.Count][];
			double[][] tar_dist = new double[this.mostWeatheredPoints.Count][];
			int I = 0;
			this.leastWeatheredPoints.ForEach(p => src_dist[I++] = this.ComputeDistanceFrom(p));
			I = 0;
			this.mostWeatheredPoints.ForEach(p => tar_dist[I++] = this.ComputeDistanceFrom(p));

			// compute the weathering degree for each manifold point
			int w = this.labImage.Width, h = this.labImage.Width;
			foreach (ManifoldPoint p in this.manifoldPoints)
			{
				// find the closest points to src/tar clusters
				double src_mind = double.MaxValue;
				foreach (double[] dist in src_dist)
				{
					double d = dist[p.Index];
					if (d < src_mind)
					{
						src_mind = d;
					}
				}
				double tar_mind = double.MaxValue;
				foreach (double[] dist in tar_dist)
				{
					double d = dist[p.Index];
					if (d < tar_mind)
					{
						tar_mind = d;
					}
				}
				p.WeatheringDegree = src_mind / (src_mind + tar_mind);
			}
		}
		private void UpdateWeatheringComponent()		// update the L component (weathering reflactance) based on the weathering degrees
		{
			// divide the (0,1) into 100 bins and average
			int numbins = 100;
			double step = 1.0/numbins;
			List<ManifoldPoint>[] bins = new List<ManifoldPoint>[numbins];
			for (int i = 0; i < numbins; ++i) bins[i] = new List<ManifoldPoint>();
			for (int i = 0; i < numbins; ++i)
			{
				foreach (ManifoldPoint mp in this.manifoldPoints)
				{
					double wd = mp.WeatheringDegree;
					int j = (int)(wd / step);
					bins[j].Add(mp);
				}
			}
			// average through each bin
			foreach (List<ManifoldPoint> bin in bins)
			{
				double avgwd = 0;
				foreach (ManifoldPoint mp in bin)
				{
					avgwd += mp.Il;
				}
				avgwd /= bin.Count;
				foreach (ManifoldPoint mp in bin)
				{
					mp.Wl = mp.Pos.x = avgwd;
				}
			}
		}
		private void ConvexilizePoints()
		{
			// dilate and make convex the least- and most- weathered points
			// compute distances from src/tar points to all other points
			double[][] src_dist = new double[this.leastWeatheredPoints.Count][];
			double[][] tar_dist = new double[this.mostWeatheredPoints.Count][];
			int I = 0;
			this.leastWeatheredPoints.ForEach(p => src_dist[I++] = this.ComputeDistanceFrom(p));
			I = 0;
			this.mostWeatheredPoints.ForEach(p => tar_dist[I++] = this.ComputeDistanceFrom(p));

			// find the minmum distance between points in source and targets
			double geoV02V1 = double.MaxValue;
			foreach (ManifoldPoint mp in this.leastWeatheredPoints)
			{
				foreach (double[] dist in tar_dist)
				{
					double d = dist[mp.Index];
					if (d < geoV02V1)
					{
						geoV02V1 = d;
					}
				}
			}
			geoV02V1 *= 0.9;

			// dilate the source (least weathered points)
			foreach (ManifoldPoint mp in this.manifoldPoints)
			{
				if (this.leastWeatheredPoints.Contains(mp)) continue;
				double geovvo = double.MaxValue;
				foreach (double[] dist in tar_dist)
				{
					double d = dist[mp.Index];
					if (d < geovvo)
					{
						geovvo = d;
					}
				}
				if (geovvo > geoV02V1)
				{
					this.leastWeatheredPoints.Add(mp);
				}
			}

			// dilate the target (most weathered points)
			foreach (ManifoldPoint mp in this.manifoldPoints)
			{
				if (this.mostWeatheredPoints.Contains(mp)) continue;
				double geovvo = double.MaxValue;
				foreach (double[] dist in src_dist)
				{
					double d = dist[mp.Index];
					if (d < geovvo)
					{
						geovvo = d;
					}
				}
				if (geovvo > geoV02V1)
				{
					this.mostWeatheredPoints.Add(mp);
				}
			}


			// convexilize
			List<ManifoldPoint> res_points = new List<ManifoldPoint>();
			foreach (ManifoldPoint mp in this.leastWeatheredPoints)
			{
				this.ComputeDistanceFrom(mp);
				foreach (ManifoldPoint mp2 in this.leastWeatheredPoints)
				{
					ManifoldPoint tmp = mp2;
					while (tmp.Parent != null)
					{
						res_points.Add(tmp.Parent);
						tmp = tmp.Parent;
					}
				}
			}
			foreach (ManifoldPoint mp in res_points)
			{
				if (!this.leastWeatheredPoints.Contains(mp))
				{
					this.leastWeatheredPoints.Add(mp);
				}
			}
			
			// the most-weathered points part
			res_points.Clear();
			foreach (ManifoldPoint mp in this.mostWeatheredPoints)
			{
				this.ComputeDistanceFrom(mp);
				foreach (ManifoldPoint mp2 in this.mostWeatheredPoints)
				{
					ManifoldPoint tmp = mp2;
					while (tmp.Parent != null)
					{
						res_points.Add(tmp.Parent);
						tmp = tmp.Parent;
					}
				}
			}
			foreach (ManifoldPoint mp in res_points)
			{
				if (!this.mostWeatheredPoints.Contains(mp))
				{
					this.mostWeatheredPoints.Add(mp);
				}
			}
		}
		private double[] ComputeDistanceFrom(ManifoldPoint rec)
		{
			int w = this.labImage.Width, h = this.labImage.Height;
			int n = w*h;
			PriorityQueue queue = new PriorityQueue(n);

			rec.Distance = 0;	// initialize seed
			foreach (ManifoldPoint p in this.manifoldPoints)
				queue.Insert(p);

			while (queue.IsEmpty() == false)
			{
				ManifoldPoint u = queue.DeleteMin() as ManifoldPoint;
				foreach (ManifoldPoint v in u.Neighbors)
				{
					double dis = (u.Pos - v.Pos).Length();
					if (dis < v.Distance)
					{
						v.Parent = u;
						v.Distance = dis;
						if (!queue.IsEmpty())
							queue.Update(v);
					}
				}
			}

			double[] distance = new double[n];
			for (int i = 0; i < n; i++)
				distance[i] = this.manifoldPoints[i/h, i%h].Distance;
			return distance;
		}
		private void ComputeShadingComponent()
		{
			// find the least- and most- bright pixel in Ic channels
			// to determine the window radious r(x,y) = p*e(-q*Il(x,y)), so that
			// r(brightiest pixel) = 1, r(darkest pixel) = 20;
			double brightiest_Il = double.MinValue;
			double darkest_Il = double.MaxValue;
			foreach (ManifoldPoint mp in this.manifoldPoints)
			{
				if (mp.Il > brightiest_Il) brightiest_Il = mp.Il;
				if (mp.Il < darkest_Il) darkest_Il = mp.Il;
			}
			double q = Math.Log(20.0) / (brightiest_Il - darkest_Il);
			double p = Math.Exp(q * brightiest_Il);

			// now compute the hist_xy for each pixel
			// first fill the pixels into buckets w.r.t. ab channels
			List<ManifoldPoint>[,] ab_buckets = new List<ManifoldPoint>[255,255];
			for (int i = 0; i < 255; ++i)
			{
				for (int j = 0; j < 255; ++j)
				{
					ab_buckets[i,j] = new List<ManifoldPoint>();
				}
			}
			foreach (ManifoldPoint mp in this.manifoldPoints)
			{
				int I = (int)mp.Pos.y + 128;
				int J = (int)mp.Pos.z + 128;
				ab_buckets[I,J].Add(mp);
			}

			int w = this.labImage.Width, h = this.labImage.Height;
			double[,][] pixel_pdfs = new double[w,h][];				// for each pixel, initialize a pdf
			double[,][] new_pixel_pdfs = new double[w, h][];		// for each pixel, initialize a new pdf
			Random rand = new Random();
			int K = 500; 
			int numbins = 100;
			double[] x_val = new double[numbins];
			for (int i = 0; i < numbins; ++i)						// the x values for pdfs
			{
				x_val[i] = (double)(i + 0.5) / numbins;
			}

			
			// copy weathering degree
			foreach (ManifoldPoint mp in this.manifoldPoints)
			{
				mp.OldWeatheringDegree = mp.WeatheringDegree;
			}


			// iteration
			int itr = 0, maxitr = 5;
			do
			{
				foreach (ManifoldPoint mp in this.manifoldPoints)
				{
					int r = (int)(p * Math.Exp(-q * mp.Il));

					// get all manifold points within the radius r
					List<ManifoldPoint> ab_points = new List<ManifoldPoint>();
					for (int i = -r; i < r; ++i)
					{
						for (int j = -r; j < r; ++j)
						{
							int I = (int)mp.Pos.y + 128 + i;
							int J = (int)mp.Pos.z + 128 + j;
							if (I >= 0 && I < 255 && J >= 0 && J < 255)
							{
								ab_points.AddRange(ab_buckets[I, J]);
							}
						}
					}

					// now radom sample K points -- might be repeated
					double[] pdf = new double[numbins];
					int ncount = ab_points.Count;
					if (ncount > 0)
					{
						List<ManifoldPoint> hist_points = new List<ManifoldPoint>();
						for (int i = 0; i < K; ++i)
						{
							hist_points.Add(ab_points[rand.Next(ncount)]);
						}
						// build the histogram -- pdf
						// divid the weathering degree into 100 uniform bins
						double step = 1.0 / numbins;
						foreach (ManifoldPoint pp in hist_points)
						{
							int bin_index = (int)(pp.WeatheringDegree / step) % numbins;
							pdf[bin_index] += 1.0 / (hist_points.Count); // normalization
						}
					}
					pixel_pdfs[mp.Index / h, mp.Index % h] = pdf;
				}

				// convolve the pixel weathering degrees
				int windowsize = 2;
				for (int i = 0; i < w; ++i)
				{
					for (int j = 0; j < h; ++j)
					{
						double[] pdf1 = pixel_pdfs[i, j];
						double[] res = pdf1.Clone() as double[];
						for (int k = -windowsize; k < windowsize; ++k)
						{
							int I = i + k, J = j + k;
							if (I >= 0 && I < w && J >= 0 && J < h)
							{
								double[] pdf2 = pixel_pdfs[I, J];
								for (int l = 0; l < pdf1.Length; ++l)
								{
									res[l] *= pdf2[l];
								}
							}
						}
						this.NormalizeData(res);
						new_pixel_pdfs[i, j] = res;
					}
				}
				// assign back the weathering degree
				for (int i = 0; i < w; ++i)
				{
					for (int j = 0; j < h; ++j)
					{
						double[] pdf1 = pixel_pdfs[i, j];
						double[] pdf2 = new_pixel_pdfs[i, j];
						double mean1, var1, mean2, var2;
						this.ComputeMeanAndVariance(pdf1, x_val, out mean1, out var1);
						this.ComputeMeanAndVariance(pdf2, x_val, out mean2, out var2);
						if (var1 < var2) // confidence increase
						{
							this.manifoldPoints[i, j].WeatheringDegree = mean2;
						}
					}
				}
			}
			while (itr++ < maxitr);

			// now update the weathering (reflectance) component & compute the shading component
			List<ManifoldPoint>[] w_buckets = this.PackPointsIntoBins(1000);

			double[,] wl = new double[w,h];	// compute the new weathering component
			for (int i = 0; i < w; ++i)
			{
				for (int j = 0; j < h; ++j)
				{
					ManifoldPoint curr_point = this.manifoldPoints[i, j];
					double wd = curr_point.WeatheringDegree;
					int btk_idx = (int)(wd / 1e-3);
					List<ManifoldPoint> mp_list = new List<ManifoldPoint>();
					foreach (ManifoldPoint mp in w_buckets[btk_idx])		// get all points with weather defree = wd
					{
						if (Math.Abs(mp.OldWeatheringDegree - wd) < 1e-5)
						{
							mp_list.Add(mp);
						}
					}
					double mincrama = double.MaxValue;	// get the one point with cloest ab value
					foreach (ManifoldPoint mp in mp_list)
					{
						double dy = (curr_point.Pos.y - mp.Pos.y);
						double dz = (curr_point.Pos.z - mp.Pos.z);
						double dd = dy * dy + dz * dz;
						if (dd < mincrama)
						{
							mincrama = dd;
							wl[i, j] = mp.Wl;			// assign the weathering component
						}
					}
				}
			}
			// compute the shading component
			for (int i = 0; i < w; ++i)
			{
				for (int j = 0; j < h; ++j)
				{
					ManifoldPoint curr_point = this.manifoldPoints[i, j];
					double _wl = wl[i, j];
					if (_wl > 0)
						curr_point.Sl = curr_point.Il / wl[i, j];
					else
						curr_point.Sl = curr_point.Il / curr_point.Wl;
				}
			}
		}
		private List<ManifoldPoint>[] PackPointsIntoBins(int numbins)
		{
			double intv = 1.0 / numbins;
			List<ManifoldPoint>[] w_buckets = new List<ManifoldPoint>[numbins];
			for (int i = 0; i < numbins; ++i) w_buckets[i] = new List<ManifoldPoint>();
			foreach (ManifoldPoint mp in this.manifoldPoints)
			{
				int btk_idx = (int)(mp.WeatheringDegree / intv);
				w_buckets[btk_idx].Add(mp);
			}
			return w_buckets;
		}
		private List<ManifoldPoint> FindPointsWithDegree(List<ManifoldPoint> list, double degree)
		{
			List<ManifoldPoint> points = new List<ManifoldPoint>();
			foreach (ManifoldPoint mp in list)
			{
				if (Math.Abs(mp.WeatheringDegree - degree) < 1e-6)
				{
					points.Add(mp);
				}
			}
			return points;
		}
		private void ComputeManifoldGeodesic()
		{
			int w = this.labImage.Width, h = this.labImage.Height;
			this.manifoldGeodesic = new double[w, h][];
			for (int i = 0; i < w; ++i)
				for (int j = 0; j < h; ++j)
				{
					this.manifoldGeodesic[i,j] = this.ComputeDistanceFrom(this.manifoldPoints[i,j]);
				}
		}
		private void NormalizeData(double[] data)
		{
			double sum = 0;
			foreach (double d in data) sum += d;
			if (sum == 0) return;
			for (int i = 0; i < data.Length; ++i)
				data[i] /= sum;
		}
		// compute the mean (estimation) and variance of a pdf
		private void ComputeMeanAndVariance(double[] pdf, double[] x, out double mean, out double variance) 
		{
			mean = 0; variance = 0;
			for (int i = 0; i < pdf.Length; ++i)
			{
				mean += pdf[i] * x[i];
			}
			for (int i = 0; i < pdf.Length; ++i)
			{
				variance += pdf[i] * x[i] * x[i];
			}
			variance -= mean*mean;
		}
		
	}
}
