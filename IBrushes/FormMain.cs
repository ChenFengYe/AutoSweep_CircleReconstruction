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

                case SketchView.EnumOperationMode.Choosing:
                    this.sketchView.CurrCanvasEngine.rect_shown = true;
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
                case Keys.C:

                    break;

                case Keys.D:

                    break;
                case Keys.F1:

                    break;
                case Keys.F2:

                    break;
            }
        }



        private void buttonCloseTab_Click(object sender, EventArgs e)
        {
            if (this.sketchView.CurrCanvasEngine != null)
            {
                this.sketchView.RemoveImageRecord(this.sketchView.CurrCanvasEngine);
            }
            this.sketchView.Clear();
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

        private void openRGBDV1DataToolStripMenuItem_Click(object sender, EventArgs e)
        {
            // the canvas texture
            OpenFileDialog d = new OpenFileDialog();
            d.FileName = "";
            d.Filter = "RGBD image File (*.png)|*.png|All files (*.*)|*.*";
            d.CheckFileExists = true;
            d.Multiselect = true;
            if (d.ShowDialog(this) == DialogResult.OK)
            {
                RGBDBuilder rgbd_builder = new RGBDBuilder();
                int kinect_ver = rgbd_builder.ReadRGBDImageDepthFileV1_new(d.FileName);


                Image<Bgr, byte> image = rgbd_builder.CurrFrame().Image;

                if (image != null)
                {
                    this.sketchView.CreateCanvasEngine(image);
                }
                this.sketchView.CurrCanvasEngine.CreateKinectDefaultCamera(kinect_ver);
                this.sketchView.CurrCanvasEngine.SetRGBDBuilder(rgbd_builder);

                // save file prefix for late IO
                string tmp = d.FileName.Substring(d.FileName.LastIndexOf('\\') + 1);
                this.sketchView.CurrCanvasEngine.file_prefix = tmp.Substring(0, tmp.LastIndexOf('.'));

                this.sketchView.Refresh();
                this.sketchView.Focus();

            }
        }

        private void sketchView_Load(object sender, EventArgs e)
        {

        }


        private void tsbDrawAndCut_Click(object sender, EventArgs e)
        {
            this.SetMode(SketchView.EnumOperationMode.Choosing);
        }


        private void toolStripButton2_Click(object sender, EventArgs e)
        {
            this.sketchView.CurrCanvasEngine.FindOptimalPlanes2();
            this.sketchView.Refresh();
            this.sketchView.Focus();
        }

        private void tsb_Loadplane_Click(object sender, EventArgs e)
        {
            // choose plane file
            OpenFileDialog d = new OpenFileDialog();
            d.FileName = "";
            d.Filter = "Plane File (*.plane)|*.plane";
            d.CheckFileExists = true;
            d.Multiselect = true;
            if (d.ShowDialog(this) == DialogResult.OK)
            {
                string plane_prefix = d.FileName.Substring(d.FileName.LastIndexOf('\\') + 1);
                plane_prefix = plane_prefix.Substring(0, plane_prefix.LastIndexOf('_'));
                if (plane_prefix == this.sketchView.CurrCanvasEngine.file_prefix)
                {
                    // load plane, model and outline
                    this.sketchView.CurrCanvasEngine.LoadPlane(d.FileName);
                    this.sketchView.Refresh();
                    this.sketchView.Focus();
                }
            }
        }
        private void tsb_LoadModel_Click(object sender, EventArgs e)
        {
            // choose plane file
            OpenFileDialog d = new OpenFileDialog();
            d.FileName = "";
            d.Filter = "Model File (*.model)|*.model";
            d.CheckFileExists = true;
            d.Multiselect = true;
            if (d.ShowDialog(this) == DialogResult.OK)
            {
                string model_prefix = d.FileName.Substring(d.FileName.LastIndexOf('\\') + 1);
                model_prefix = model_prefix.Substring(0, model_prefix.LastIndexOf('_'));
                if (model_prefix == this.sketchView.CurrCanvasEngine.file_prefix)
                {
                    // load model
                    this.sketchView.CurrCanvasEngine.LoadModel(d.FileName);
                    this.sketchView.Refresh();
                    this.sketchView.Focus();
                }
            }
        }

        private void stb_LoadOutline_Click(object sender, EventArgs e)
        {
            // choose plane file
            OpenFileDialog d = new OpenFileDialog();
            d.FileName = "";
            d.Filter = "outline File (*.outline)|*.outline";
            d.CheckFileExists = true;
            d.Multiselect = true;
            if (d.ShowDialog(this) == DialogResult.OK)
            {
                string outline_prefix = d.FileName.Substring(d.FileName.LastIndexOf('\\') + 1);
                outline_prefix = outline_prefix.Substring(0, outline_prefix.LastIndexOf('_'));
                if (outline_prefix == this.sketchView.CurrCanvasEngine.file_prefix)
                {
                    // load outline
                    this.sketchView.CurrCanvasEngine.LoadOutline(d.FileName);
                    this.sketchView.Refresh();
                    this.sketchView.Focus();
                }
            }
        }

        private void stb_Saveplane_Click(object sender, EventArgs e)
        {
            // save plane, model and outline
            SaveFileDialog d = new SaveFileDialog();
            d.FileName = this.sketchView.CurrCanvasEngine.file_prefix;
            //d.Filter = "Plane File (*.plane)|*.plane";
            d.FilterIndex = 1;
            if (d.ShowDialog(this) == DialogResult.OK)
            {
                string common_prefix = d.FileName;
                this.sketchView.CurrCanvasEngine.SaveAllPlane(common_prefix);
                this.sketchView.CurrCanvasEngine.SaveAllModel(common_prefix);
                this.sketchView.CurrCanvasEngine.SaveAllOutline(common_prefix);
                this.sketchView.Refresh();
                this.sketchView.Focus();
            }
        }

        private void tsb_Reconstruct_Click(object sender, EventArgs e)
        {
            //this.sketchView.CurrCanvasEngine.ReconstructGCylinder();
            this.sketchView.Refresh();
            this.sketchView.Focus();
        }

        private void showRgbdpointsToolStripMenuItem_Click(object sender, EventArgs e)
        {
            this.showRgbdpointsToolStripMenuItem.Checked = !this.showRgbdpointsToolStripMenuItem.Checked;
            if (this.sketchView != null && this.sketchView.CurrCanvasEngine != null)
            {
                this.sketchView.CurrCanvasEngine.showrgbdpoints = this.showRgbdpointsToolStripMenuItem.Checked;
            }
            this.Refresh();
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
    }
}

