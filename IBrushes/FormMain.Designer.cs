namespace SmartCanvas
{
	partial class FormAirBrush
	{
		/// <summary>
		/// Required designer variable.
		/// </summary>
		private System.ComponentModel.IContainer components = null;

		/// <summary>
		/// Clean up any resources being used.
		/// </summary>
		/// <param name="disposing">true if managed resources should be disposed; otherwise, false.</param>
		protected override void Dispose(bool disposing)
		{
			if (disposing && (components != null))
			{
				components.Dispose();
			}
			base.Dispose(disposing);
		}

		#region Windows Form Designer generated code

		/// <summary>
		/// Required method for Designer support - do not modify
		/// the contents of this method with the code editor.
		/// </summary>
		private void InitializeComponent()
		{
            this.components = new System.ComponentModel.Container();
            System.ComponentModel.ComponentResourceManager resources = new System.ComponentModel.ComponentResourceManager(typeof(FormAirBrush));
            this.toolStrip1 = new System.Windows.Forms.ToolStrip();
            this.toolStripButton1 = new System.Windows.Forms.ToolStripDropDownButton();
            this.openRGBDV1DataToolStripMenuItem = new System.Windows.Forms.ToolStripMenuItem();
            this.toolStripButtonView = new System.Windows.Forms.ToolStripButton();
            this.toolStripButtonSetVIewPoint = new System.Windows.Forms.ToolStripDropDownButton();
            this.showBackgroundToolStripMenuItem = new System.Windows.Forms.ToolStripMenuItem();
            this.showRgbdpointsToolStripMenuItem = new System.Windows.Forms.ToolStripMenuItem();
            this.showAxisToolStripMenuItem = new System.Windows.Forms.ToolStripMenuItem();
            this.toolStripButtonRefresh = new System.Windows.Forms.ToolStripButton();
            this.stb_ransac = new System.Windows.Forms.ToolStripButton();
            this.tsbGrabcut = new System.Windows.Forms.ToolStripButton();
            this.stb_SavePMO = new System.Windows.Forms.ToolStripButton();
            this.toolStripSeparator1 = new System.Windows.Forms.ToolStripSeparator();
            this.tsb_Loadplane = new System.Windows.Forms.ToolStripButton();
            this.tsb_LoadModel = new System.Windows.Forms.ToolStripButton();
            this.stb_LoadOutline = new System.Windows.Forms.ToolStripButton();
            this.tsb_Reconstruct = new System.Windows.Forms.ToolStripButton();
            this.showSketchyOutlinesToolStripMenuItem = new System.Windows.Forms.ToolStripMenuItem();
            this.panelMain = new System.Windows.Forms.Panel();
            this.sketchView = new SmartCanvas.SketchView();
            this.buttonCloseTab = new System.Windows.Forms.Button();
            this.colorDialog1 = new System.Windows.Forms.ColorDialog();
            this.showOutLineToolStripMenuItem = new System.Windows.Forms.ToolStripMenuItem();
            this.contextMenuStrip = new System.Windows.Forms.ContextMenuStrip(this.components);
            this.groupToolStripMenuItem = new System.Windows.Forms.ToolStripMenuItem();
            this.ungroupToolStripMenuItem = new System.Windows.Forms.ToolStripMenuItem();
            this.sendBackSToolStripMenuItem = new System.Windows.Forms.ToolStripMenuItem();
            this.junctionToolStripMenuItem = new System.Windows.Forms.ToolStripMenuItem();
            this.rotTime = new System.Windows.Forms.Timer(this.components);
            this.toolStripButton2 = new System.Windows.Forms.ToolStripButton();
            this.toolStrip1.SuspendLayout();
            this.panelMain.SuspendLayout();
            this.contextMenuStrip.SuspendLayout();
            this.SuspendLayout();
            // 
            // toolStrip1
            // 
            this.toolStrip1.ImageScalingSize = new System.Drawing.Size(32, 32);
            this.toolStrip1.Items.AddRange(new System.Windows.Forms.ToolStripItem[] {
            this.toolStripButton1,
            this.toolStripButtonView,
            this.toolStripButtonSetVIewPoint,
            this.toolStripButtonRefresh,
            this.stb_ransac,
            this.tsbGrabcut,
            this.stb_SavePMO,
            this.toolStripSeparator1,
            this.tsb_Loadplane,
            this.tsb_LoadModel,
            this.stb_LoadOutline,
            this.tsb_Reconstruct,
            this.toolStripButton2});
            this.toolStrip1.Location = new System.Drawing.Point(0, 0);
            this.toolStrip1.Name = "toolStrip1";
            this.toolStrip1.Size = new System.Drawing.Size(748, 39);
            this.toolStrip1.TabIndex = 0;
            this.toolStrip1.Text = "toolStrip1";
            // 
            // toolStripButton1
            // 
            this.toolStripButton1.DisplayStyle = System.Windows.Forms.ToolStripItemDisplayStyle.Image;
            this.toolStripButton1.DropDownItems.AddRange(new System.Windows.Forms.ToolStripItem[] {
            this.openRGBDV1DataToolStripMenuItem});
            this.toolStripButton1.Image = ((System.Drawing.Image)(resources.GetObject("toolStripButton1.Image")));
            this.toolStripButton1.ImageTransparentColor = System.Drawing.Color.Magenta;
            this.toolStripButton1.Name = "toolStripButton1";
            this.toolStripButton1.Size = new System.Drawing.Size(45, 36);
            this.toolStripButton1.Text = "Open";
            // 
            // openRGBDV1DataToolStripMenuItem
            // 
            this.openRGBDV1DataToolStripMenuItem.Name = "openRGBDV1DataToolStripMenuItem";
            this.openRGBDV1DataToolStripMenuItem.Size = new System.Drawing.Size(194, 22);
            this.openRGBDV1DataToolStripMenuItem.Text = "Open RGBD v1 Data";
            this.openRGBDV1DataToolStripMenuItem.Click += new System.EventHandler(this.openRGBDV1DataToolStripMenuItem_Click);
            // 
            // toolStripButtonView
            // 
            this.toolStripButtonView.DisplayStyle = System.Windows.Forms.ToolStripItemDisplayStyle.Image;
            this.toolStripButtonView.Image = ((System.Drawing.Image)(resources.GetObject("toolStripButtonView.Image")));
            this.toolStripButtonView.ImageTransparentColor = System.Drawing.Color.Magenta;
            this.toolStripButtonView.Name = "toolStripButtonView";
            this.toolStripButtonView.Size = new System.Drawing.Size(36, 36);
            this.toolStripButtonView.Text = "toolStripButton2";
            this.toolStripButtonView.ToolTipText = "View mode";
            this.toolStripButtonView.Click += new System.EventHandler(this.toolStripButtonView_Click);
            // 
            // toolStripButtonSetVIewPoint
            // 
            this.toolStripButtonSetVIewPoint.DisplayStyle = System.Windows.Forms.ToolStripItemDisplayStyle.Image;
            this.toolStripButtonSetVIewPoint.DropDownItems.AddRange(new System.Windows.Forms.ToolStripItem[] {
            this.showBackgroundToolStripMenuItem,
            this.showRgbdpointsToolStripMenuItem,
            this.showAxisToolStripMenuItem});
            this.toolStripButtonSetVIewPoint.Image = ((System.Drawing.Image)(resources.GetObject("toolStripButtonSetVIewPoint.Image")));
            this.toolStripButtonSetVIewPoint.ImageTransparentColor = System.Drawing.Color.Magenta;
            this.toolStripButtonSetVIewPoint.Name = "toolStripButtonSetVIewPoint";
            this.toolStripButtonSetVIewPoint.Size = new System.Drawing.Size(45, 36);
            this.toolStripButtonSetVIewPoint.Text = "toolStripButton6";
            this.toolStripButtonSetVIewPoint.ToolTipText = "switch original/sketch image";
            // 
            // showBackgroundToolStripMenuItem
            // 
            this.showBackgroundToolStripMenuItem.Name = "showBackgroundToolStripMenuItem";
            this.showBackgroundToolStripMenuItem.Size = new System.Drawing.Size(181, 22);
            this.showBackgroundToolStripMenuItem.Text = "show background";
            this.showBackgroundToolStripMenuItem.Click += new System.EventHandler(this.showBackgroundToolStripMenuItem_Click);
            // 
            // showRgbdpointsToolStripMenuItem
            // 
            this.showRgbdpointsToolStripMenuItem.Name = "showRgbdpointsToolStripMenuItem";
            this.showRgbdpointsToolStripMenuItem.Size = new System.Drawing.Size(181, 22);
            this.showRgbdpointsToolStripMenuItem.Text = "show rgbdpoints";
            this.showRgbdpointsToolStripMenuItem.Click += new System.EventHandler(this.showRgbdpointsToolStripMenuItem_Click);
            // 
            // showAxisToolStripMenuItem
            // 
            this.showAxisToolStripMenuItem.Name = "showAxisToolStripMenuItem";
            this.showAxisToolStripMenuItem.Size = new System.Drawing.Size(181, 22);
            this.showAxisToolStripMenuItem.Text = "show axis";
            this.showAxisToolStripMenuItem.Click += new System.EventHandler(this.showAxisToolStripMenuItem_Click);
            // 
            // toolStripButtonRefresh
            // 
            this.toolStripButtonRefresh.DisplayStyle = System.Windows.Forms.ToolStripItemDisplayStyle.Image;
            this.toolStripButtonRefresh.Image = ((System.Drawing.Image)(resources.GetObject("toolStripButtonRefresh.Image")));
            this.toolStripButtonRefresh.ImageTransparentColor = System.Drawing.Color.Magenta;
            this.toolStripButtonRefresh.Name = "toolStripButtonRefresh";
            this.toolStripButtonRefresh.Size = new System.Drawing.Size(36, 36);
            this.toolStripButtonRefresh.Text = "Refresh";
            this.toolStripButtonRefresh.ToolTipText = "Reset camera View";
            this.toolStripButtonRefresh.Click += new System.EventHandler(this.toolStripButtonRefresh_Click);
            // 
            // stb_ransac
            // 
            this.stb_ransac.DisplayStyle = System.Windows.Forms.ToolStripItemDisplayStyle.Text;
            this.stb_ransac.Image = ((System.Drawing.Image)(resources.GetObject("stb_ransac.Image")));
            this.stb_ransac.ImageTransparentColor = System.Drawing.Color.Magenta;
            this.stb_ransac.Name = "stb_ransac";
            this.stb_ransac.Size = new System.Drawing.Size(53, 36);
            this.stb_ransac.Text = "Ransac";
            this.stb_ransac.ToolTipText = "Find the optimal plane by ransac.";
            this.stb_ransac.Click += new System.EventHandler(this.toolStripButton2_Click);
            // 
            // tsbGrabcut
            // 
            this.tsbGrabcut.DisplayStyle = System.Windows.Forms.ToolStripItemDisplayStyle.Text;
            this.tsbGrabcut.Image = ((System.Drawing.Image)(resources.GetObject("tsbGrabcut.Image")));
            this.tsbGrabcut.ImageTransparentColor = System.Drawing.Color.Magenta;
            this.tsbGrabcut.Name = "tsbGrabcut";
            this.tsbGrabcut.Size = new System.Drawing.Size(58, 36);
            this.tsbGrabcut.Text = "Grabcut";
            this.tsbGrabcut.ToolTipText = "Draw rectangle and do others automatically.";
            this.tsbGrabcut.Click += new System.EventHandler(this.tsbDrawAndCut_Click);
            // 
            // stb_SavePMO
            // 
            this.stb_SavePMO.DisplayStyle = System.Windows.Forms.ToolStripItemDisplayStyle.Text;
            this.stb_SavePMO.Image = ((System.Drawing.Image)(resources.GetObject("stb_SavePMO.Image")));
            this.stb_SavePMO.ImageTransparentColor = System.Drawing.Color.Magenta;
            this.stb_SavePMO.Name = "stb_SavePMO";
            this.stb_SavePMO.Size = new System.Drawing.Size(68, 36);
            this.stb_SavePMO.Text = "SavePMO";
            this.stb_SavePMO.ToolTipText = "Save Plane, model and outline.";
            this.stb_SavePMO.Click += new System.EventHandler(this.stb_Saveplane_Click);
            // 
            // toolStripSeparator1
            // 
            this.toolStripSeparator1.Name = "toolStripSeparator1";
            this.toolStripSeparator1.Size = new System.Drawing.Size(6, 39);
            // 
            // tsb_Loadplane
            // 
            this.tsb_Loadplane.DisplayStyle = System.Windows.Forms.ToolStripItemDisplayStyle.Text;
            this.tsb_Loadplane.Image = ((System.Drawing.Image)(resources.GetObject("tsb_Loadplane.Image")));
            this.tsb_Loadplane.ImageTransparentColor = System.Drawing.Color.Magenta;
            this.tsb_Loadplane.Name = "tsb_Loadplane";
            this.tsb_Loadplane.Size = new System.Drawing.Size(72, 36);
            this.tsb_Loadplane.Text = "LoadPlane";
            this.tsb_Loadplane.ToolTipText = "Choose a plane to load.";
            this.tsb_Loadplane.Click += new System.EventHandler(this.tsb_Loadplane_Click);
            // 
            // tsb_LoadModel
            // 
            this.tsb_LoadModel.DisplayStyle = System.Windows.Forms.ToolStripItemDisplayStyle.Text;
            this.tsb_LoadModel.Image = ((System.Drawing.Image)(resources.GetObject("tsb_LoadModel.Image")));
            this.tsb_LoadModel.ImageTransparentColor = System.Drawing.Color.Magenta;
            this.tsb_LoadModel.Name = "tsb_LoadModel";
            this.tsb_LoadModel.Size = new System.Drawing.Size(79, 36);
            this.tsb_LoadModel.Text = "LoadModel";
            this.tsb_LoadModel.ToolTipText = "Choose a model to load.";
            this.tsb_LoadModel.Click += new System.EventHandler(this.tsb_LoadModel_Click);
            // 
            // stb_LoadOutline
            // 
            this.stb_LoadOutline.DisplayStyle = System.Windows.Forms.ToolStripItemDisplayStyle.Text;
            this.stb_LoadOutline.Image = ((System.Drawing.Image)(resources.GetObject("stb_LoadOutline.Image")));
            this.stb_LoadOutline.ImageTransparentColor = System.Drawing.Color.Magenta;
            this.stb_LoadOutline.Name = "stb_LoadOutline";
            this.stb_LoadOutline.Size = new System.Drawing.Size(82, 36);
            this.stb_LoadOutline.Text = "LoadOutline";
            this.stb_LoadOutline.ToolTipText = "Choose an outline to load.";
            this.stb_LoadOutline.Click += new System.EventHandler(this.stb_LoadOutline_Click);
            // 
            // tsb_Reconstruct
            // 
            this.tsb_Reconstruct.DisplayStyle = System.Windows.Forms.ToolStripItemDisplayStyle.Text;
            this.tsb_Reconstruct.Image = ((System.Drawing.Image)(resources.GetObject("tsb_Reconstruct.Image")));
            this.tsb_Reconstruct.ImageTransparentColor = System.Drawing.Color.Magenta;
            this.tsb_Reconstruct.Name = "tsb_Reconstruct";
            this.tsb_Reconstruct.Size = new System.Drawing.Size(80, 36);
            this.tsb_Reconstruct.Text = "Reconstruct";
            this.tsb_Reconstruct.ToolTipText = "Press to reconstruct gcylinder.";
            this.tsb_Reconstruct.Click += new System.EventHandler(this.tsb_Reconstruct_Click);
            // 
            // showSketchyOutlinesToolStripMenuItem
            // 
            this.showSketchyOutlinesToolStripMenuItem.Name = "showSketchyOutlinesToolStripMenuItem";
            this.showSketchyOutlinesToolStripMenuItem.Size = new System.Drawing.Size(32, 19);
            // 
            // panelMain
            // 
            this.panelMain.Anchor = ((System.Windows.Forms.AnchorStyles)((((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Bottom) 
            | System.Windows.Forms.AnchorStyles.Left) 
            | System.Windows.Forms.AnchorStyles.Right)));
            this.panelMain.Controls.Add(this.sketchView);
            this.panelMain.Controls.Add(this.buttonCloseTab);
            this.panelMain.Location = new System.Drawing.Point(0, 39);
            this.panelMain.Name = "panelMain";
            this.panelMain.Size = new System.Drawing.Size(748, 578);
            this.panelMain.TabIndex = 1;
            // 
            // sketchView
            // 
            this.sketchView.Anchor = ((System.Windows.Forms.AnchorStyles)((((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Bottom) 
            | System.Windows.Forms.AnchorStyles.Left) 
            | System.Windows.Forms.AnchorStyles.Right)));
            this.sketchView.BackColor = System.Drawing.Color.WhiteSmoke;
            this.sketchView.CurrCanvasEngine = null;
            this.sketchView.CurrentMode = SmartCanvas.SketchView.EnumOperationMode.Viewing;
            this.sketchView.Location = new System.Drawing.Point(0, 3);
            this.sketchView.Name = "sketchView";
            this.sketchView.Size = new System.Drawing.Size(745, 572);
            this.sketchView.TabIndex = 3;
            this.sketchView.TurnOnSketchOverlay = false;
            this.sketchView.VSync = false;
            this.sketchView.Load += new System.EventHandler(this.sketchView_Load);
            // 
            // buttonCloseTab
            // 
            this.buttonCloseTab.Anchor = ((System.Windows.Forms.AnchorStyles)((System.Windows.Forms.AnchorStyles.Top | System.Windows.Forms.AnchorStyles.Right)));
            this.buttonCloseTab.Image = ((System.Drawing.Image)(resources.GetObject("buttonCloseTab.Image")));
            this.buttonCloseTab.Location = new System.Drawing.Point(433, 41);
            this.buttonCloseTab.Margin = new System.Windows.Forms.Padding(3, 2, 3, 2);
            this.buttonCloseTab.Name = "buttonCloseTab";
            this.buttonCloseTab.Size = new System.Drawing.Size(36, 22);
            this.buttonCloseTab.TabIndex = 18;
            this.buttonCloseTab.UseVisualStyleBackColor = true;
            this.buttonCloseTab.Click += new System.EventHandler(this.buttonCloseTab_Click);
            // 
            // showOutLineToolStripMenuItem
            // 
            this.showOutLineToolStripMenuItem.Name = "showOutLineToolStripMenuItem";
            this.showOutLineToolStripMenuItem.Size = new System.Drawing.Size(32, 19);
            // 
            // contextMenuStrip
            // 
            this.contextMenuStrip.Items.AddRange(new System.Windows.Forms.ToolStripItem[] {
            this.groupToolStripMenuItem,
            this.ungroupToolStripMenuItem,
            this.sendBackSToolStripMenuItem,
            this.junctionToolStripMenuItem});
            this.contextMenuStrip.Name = "contextMenuStrip";
            this.contextMenuStrip.Size = new System.Drawing.Size(176, 92);
            // 
            // groupToolStripMenuItem
            // 
            this.groupToolStripMenuItem.Image = ((System.Drawing.Image)(resources.GetObject("groupToolStripMenuItem.Image")));
            this.groupToolStripMenuItem.Name = "groupToolStripMenuItem";
            this.groupToolStripMenuItem.Size = new System.Drawing.Size(175, 22);
            this.groupToolStripMenuItem.Text = "Group strokes";
            // 
            // ungroupToolStripMenuItem
            // 
            this.ungroupToolStripMenuItem.Name = "ungroupToolStripMenuItem";
            this.ungroupToolStripMenuItem.Size = new System.Drawing.Size(175, 22);
            this.ungroupToolStripMenuItem.Text = "Ungroup strokes";
            // 
            // sendBackSToolStripMenuItem
            // 
            this.sendBackSToolStripMenuItem.Name = "sendBackSToolStripMenuItem";
            this.sendBackSToolStripMenuItem.Size = new System.Drawing.Size(175, 22);
            // 
            // junctionToolStripMenuItem
            // 
            this.junctionToolStripMenuItem.Name = "junctionToolStripMenuItem";
            this.junctionToolStripMenuItem.Size = new System.Drawing.Size(175, 22);
            // 
            // toolStripButton2
            // 
            this.toolStripButton2.DisplayStyle = System.Windows.Forms.ToolStripItemDisplayStyle.Text;
            this.toolStripButton2.Image = ((System.Drawing.Image)(resources.GetObject("toolStripButton2.Image")));
            this.toolStripButton2.ImageTransparentColor = System.Drawing.Color.Magenta;
            this.toolStripButton2.Name = "toolStripButton2";
            this.toolStripButton2.Size = new System.Drawing.Size(32, 36);
            this.toolStripButton2.Text = "Rcs";
            this.toolStripButton2.Click += new System.EventHandler(this.toolStripButton2_Click_1);
            // 
            // FormAirBrush
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 12F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(748, 617);
            this.Controls.Add(this.panelMain);
            this.Controls.Add(this.toolStrip1);
            this.Icon = ((System.Drawing.Icon)(resources.GetObject("$this.Icon")));
            this.Name = "FormAirBrush";
            this.Text = " ";
            this.KeyDown += new System.Windows.Forms.KeyEventHandler(this.FormAirBrush_KeyDown);
            this.toolStrip1.ResumeLayout(false);
            this.toolStrip1.PerformLayout();
            this.panelMain.ResumeLayout(false);
            this.contextMenuStrip.ResumeLayout(false);
            this.ResumeLayout(false);
            this.PerformLayout();

		}

		#endregion

		private System.Windows.Forms.ToolStrip toolStrip1;
		private System.Windows.Forms.Panel panelMain;
		private System.Windows.Forms.ColorDialog colorDialog1;
        private System.Windows.Forms.ToolStripDropDownButton toolStripButton1;
        private System.Windows.Forms.ToolStripButton toolStripButtonView;
		private SketchView sketchView;
		private System.Windows.Forms.ContextMenuStrip contextMenuStrip;
		private System.Windows.Forms.ToolStripMenuItem sendBackSToolStripMenuItem;
        private System.Windows.Forms.ToolStripDropDownButton toolStripButtonSetVIewPoint;
        private System.Windows.Forms.ToolStripMenuItem groupToolStripMenuItem;
        private System.Windows.Forms.ToolStripMenuItem ungroupToolStripMenuItem;
		private System.Windows.Forms.ToolStripMenuItem showSketchyOutlinesToolStripMenuItem;
        private System.Windows.Forms.ToolStripMenuItem showOutLineToolStripMenuItem;
		private System.Windows.Forms.ToolStripMenuItem showBackgroundToolStripMenuItem;
        private System.Windows.Forms.ToolStripMenuItem junctionToolStripMenuItem;
        private System.Windows.Forms.Button buttonCloseTab;
        private System.Windows.Forms.Timer rotTime;
        private System.Windows.Forms.ToolStripButton toolStripButtonRefresh;
        private System.Windows.Forms.ToolStripMenuItem openRGBDV1DataToolStripMenuItem;
		private System.Windows.Forms.ToolStripButton tsbGrabcut;
		private System.Windows.Forms.ToolStripButton stb_ransac;
		private System.Windows.Forms.ToolStripButton stb_SavePMO;
		private System.Windows.Forms.ToolStripButton tsb_Loadplane;
		private System.Windows.Forms.ToolStripSeparator toolStripSeparator1;
		private System.Windows.Forms.ToolStripButton tsb_LoadModel;
		private System.Windows.Forms.ToolStripButton stb_LoadOutline;
		private System.Windows.Forms.ToolStripButton tsb_Reconstruct;
        private System.Windows.Forms.ToolStripMenuItem showRgbdpointsToolStripMenuItem;
        private System.Windows.Forms.ToolStripMenuItem showAxisToolStripMenuItem;
        private System.Windows.Forms.ToolStripButton toolStripButton2;
	}
}

