using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Drawing;

using Tao.OpenGl;
using MyGeometry;


namespace SmartCanvas
{
	// pen interface for various brushes
	interface Pen
	{
		void SetTargetTexture(Texture tex);
		void NewStroke(Color4 color);
		void DrawPoint(float x, float y);
		void EndStroke();
	}

	// A fix-sized buffer for storing a bitmap.
	public interface Texture
	{
		// Return the blockWidth by blockHeight region starting at x,y. 
		// The block must fit into the used mip level. The size of the
		// returned array is blockWidth*blockHeight.
		Color4[] GetPixels(int x, int y, int blockWidth, int blockHeight);
		// Like the above, but writes a block of pixel values.
		void SetPixels(int x, int y, int blockWidth, int blockHeight, Color4[] colors);
	}

	// air brush for generating gaussian-blurred strokes
	public class AirBrush : Pen
	{
		// airbrush members
		private float oldsize, size;			// brush size
		private Color4 color;					// brush color
		private Texture canvasTexture;			// the canvas texture it will write colors to
		private Color4[] blurColors;				// the blur texture used to blend colors -- here a gaussian mask

		// ctrs
		public AirBrush()
		{
		}

		// public entries
		public void SetBrushSize(float size)	// this function must be called before drawing, to get correct brush size
		{
			this.size = size;

			// precompute the gaussian textures, to avoid computation overhead
			// the creation is called when the brush size changes
			if (this.size != this.oldsize)
			{
				this.blurColors =		// this texture is used for blending
					Utils.CreateGaussianTexture((int)this.size);
				this.oldsize = this.size;
			}
		}

		// implementation of interfaces 
		public void SetTargetTexture(Texture tex)
		{
			this.canvasTexture = tex;
		}
		public void NewStroke(Color4 color)
		{
			this.color = color;
		}
		public void DrawPoint(float x, float y)
		{
			if (this.blurColors == null || this.canvasTexture == null) return;

			// two-pass blending
			// pass 1. compute the blended color1 = a * penTex color + (1-a) * pen color
			// pass 2. blend color1 with canvas color
			int blocksize = (int)size;
			Color4[] canvasColors = canvasTexture.GetPixels((int)(x) - blocksize,
				(int)(y) - blocksize, blocksize * 2, blocksize * 2);

			int numpixels = (int)size * (int)size * 4;	// num of total pixels to blend
			for (int i = 0; i < numpixels; ++i)
			{
				Color4 t_color = Utils.GaussianBlend(this.blurColors[i], this.color);	// assign gaussian alphas
				canvasColors[i] = Utils.AlphaBlend(t_color, canvasColors[i]);			// do alpha blending
			}

			// assign color to the canvas
			this.canvasTexture.SetPixels((int)(x) - blocksize,
				(int)(y) - blocksize, blocksize * 2, blocksize * 2, canvasColors);

		}
		public void EndStroke()
		{

		}
	}

	public class AirTexture : Texture
	{
		// members: texture width, height,
		private int width, height;
		private Color4[] pixels;

		public int Width { get { return width; } }
		public int Height { get { return height; } }

		public AirTexture(int w, int h, Color4[] pixels)
		{
			this.width = w;
			this.height = h;
			this.pixels = pixels;
		}

		// implementation of interfaces
		// Return the blockWidth by blockHeight region starting at x,y. 
		// The block must fit into the used mip level. The size of the
		// returned array is blockWidth*blockHeight.
		public Color4[] GetPixels(int x, int y, int blockWidth, int blockHeight)
		{
			List<Color4> colorList = new List<Color4>();
			for (int i = 0; i < blockWidth; ++i)
			{
				int I = x + i;
				for (int j = 0; j < blockHeight; ++j)
				{
					int J = y + j;
					if (I >= 0 && I < this.width && J >= 0 && J < this.height)	// check bounds
					{
						colorList.Add(this.pixels[I * this.height + J]);
					}
					else
						colorList.Add(new Color4());
				}
			}
			return colorList.ToArray();
		}

		// Like the above, but writes a block of pixel values.
		public void SetPixels(int x, int y, int blockWidth, int blockHeight, Color4[] colors)
		{
			for (int i = 0; i < blockWidth; ++i)
			{
				int I = x + i;
				for (int j = 0; j < blockHeight; ++j)
				{
					int J = y + j;
					if (I >= 0 && I < this.width && J >= 0 && J < this.height)	// check bounds
					{
						this.pixels[I * this.height + J] = colors[i * blockHeight + j];
					}
				}
			}
		}
	}

}
