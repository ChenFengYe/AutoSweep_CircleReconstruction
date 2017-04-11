using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Windows.Forms;

using MyGeometry;
using System.IO;

using Emgu.CV;
using Emgu.CV.Structure;
using Emgu.Util;
using Emgu.CV.UI;
using System.Timers;

using System.Threading;
using System.Threading.Tasks;

namespace SmartCanvas
{
    public partial class FormAirBrush : Form
    {
        public FormAirBrush()
        {
            InitializeComponent();

            this.sketchView.Init();

            this.SetModeButtons();
        }

        protected override void OnFormClosing(FormClosingEventArgs e)
        {
            base.OnFormClosing(e);
        }

        private List<ToolStripButton> modeButtons = new List<ToolStripButton>();
        private void SetModeButtons()
        {
            this.modeButtons.Add(this.toolStripButtonView);
        }

        public void SetMode(SketchView.EnumOperationMode mode)
        {
            this.sketchView.CurrentMode = mode;

            foreach (ToolStripButton button in this.modeButtons)
                button.Checked = false;
            this.sketchView.ContextMenuStrip = null;

            switch (mode)
            {
                case SketchView.EnumOperationMode.Viewing:
                    this.toolStripButtonView.Checked = true;
                    break;
            }
            this.sketchView.Focus();
            this.sketchView.Refresh();
        }


        private void toolStripButtonView_Click(object sender, EventArgs e)
        {
            this.SetMode(SketchView.EnumOperationMode.Viewing);
        }


        private void FormAirBrush_KeyDown(object sender, KeyEventArgs e)
        {
            switch (e.KeyData)
            {

                case Keys.V:
                    this.SetMode(SketchView.EnumOperationMode.Viewing);
                    break;
            }
        }




        private void showBackgroundToolStripMenuItem_Click(object sender, EventArgs e)
        {
            this.showBackgroundToolStripMenuItem.Checked = !this.showBackgroundToolStripMenuItem.Checked;
            if (this.sketchView != null && this.sketchView.CurrCanvasEngine != null)
            {
                this.sketchView.CurrCanvasEngine.showBackground = this.showBackgroundToolStripMenuItem.Checked;
            }
            this.Refresh();
        }


        private void toolStripButtonRefresh_Click(object sender, EventArgs e)
        {
            this.sketchView.Reset();
        }

        private void showAxisToolStripMenuItem_Click(object sender, EventArgs e)
        {
            this.showAxisToolStripMenuItem.Checked = !this.showAxisToolStripMenuItem.Checked;
            if (this.sketchView != null && this.sketchView.CurrCanvasEngine != null)
            {
                this.sketchView.CurrCanvasEngine.showaxis = this.showAxisToolStripMenuItem.Checked;
            }
            this.Refresh();
        }

        private void toolStripButton2_Click_1(object sender, EventArgs e)
        {
            if (this.sketchView.CurrCanvasEngine == null) return;

            this.sketchView.CurrCanvasEngine.EstimatePlane();
        }

        private void loadImageToolStripMenuItem_Click(object sender, EventArgs e)
        {
            // the canvas texture
            OpenFileDialog d = new OpenFileDialog();
            d.FileName = "";
            d.Filter = "image File (*.png)|*.png|All files (*.*)|*.*";
            d.CheckFileExists = true;
            if (d.ShowDialog(this) == DialogResult.OK)
            {
                Image<Bgr, byte> img = new Image<Bgr, byte>(d.FileName);
                this.sketchView.CreateCanvasEngine(img);
                this.sketchView.CurrCanvasEngine.mark = this.sketchView.CurrCanvasEngine.GetMarkImgae(img);
                this.sketchView.CurrCanvasEngine.CreateDefaultCamera(1);
                this.sketchView.CurrCanvasEngine.SetTransformCenter();

                this.sketchView.CurrCanvasEngine.ComputeBlobScore();


                string imgdir = d.FileName.Substring(0, d.FileName.LastIndexOf('.') + 1);
                //this.sketchView.CurrCanvasEngine.ReadTopCircle(imgdir + "circle");
                // save file prefix for late IO
                string tmp = d.FileName.Substring(d.FileName.LastIndexOf('\\') + 1);
                this.sketchView.CurrCanvasEngine.file_prefix = tmp.Substring(0, tmp.LastIndexOf('.'));

                this.sketchView.Refresh();
                this.sketchView.Focus();
            }
        }

        private void toolStripButtonSweep_Click(object sender, EventArgs e)
        {
            OpenFileDialog d = new OpenFileDialog();
            d.FileName = "";
            d.Filter = "image File (*.png)|*.png|All files (*.*)|*.*";
            d.CheckFileExists = true;
            if (d.ShowDialog(this) == DialogResult.OK)
            {
                Image<Gray, byte> img = new Image<Gray, byte>(d.FileName);
                string imgdir = d.FileName.Substring(0, d.FileName.LastIndexOf('.') + 1);
                //this.sketchView.CurrCanvasEngine.ReadTopCircle(imgdir + "circle");
                this.sketchView.CurrCanvasEngine.Sweep(img);
                this.sketchView.Refresh();
                this.sketchView.Focus();
            }
        }

        private void toolStripButton2_Click(object sender, EventArgs e)
        {
            OpenFileDialog d = new OpenFileDialog();
            d.FileName = "";
            d.Filter = "image File (*.png)|*.png|All files (*.*)|*.*";
            d.CheckFileExists = true;
            if (d.ShowDialog(this) == DialogResult.OK)
            {
                this.sketchView.CurrCanvasEngine.edgeImage = new Image<Gray, byte>(d.FileName);
                string imgdir = d.FileName.Substring(0, d.FileName.LastIndexOf('.') + 1);
                this.sketchView.CurrCanvasEngine.ReadTopCircle(imgdir + "circle");
                this.sketchView.CurrCanvasEngine.CylinderSnapping();
                this.sketchView.Refresh();
                this.sketchView.Focus();
            }
        }

        private void toolStripButton3_Click(object sender, EventArgs e)
        {
            this.sketchView.CurrCanvasEngine.IsOut_Debug = true;
        }

    }
}

