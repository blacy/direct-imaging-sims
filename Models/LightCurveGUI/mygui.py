#!/usr/bin/env pythonw

import matplotlib
matplotlib.use('WXAgg')
import wx
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams
import wx.lib.agw.floatspin as FS
rcParams['font.family'] = 'serif'
rc('mathtext', default='rm')
rc('axes', labelsize=20)
rc('figure.subplot', bottom=0.12)
rc('figure.subplot', top=0.93)
rc('legend', frameon=False)
rc('legend', numpoints=1)

def title_font(text, color, size, bold=True):
    text.SetForegroundColour(color)
    font = text.GetFont()
    font.SetPointSize(size)
    if bold: font.SetWeight(wx.FONTWEIGHT_BOLD)
    text.SetFont(font)
 
class p1(wx.Panel):
    def __init__(self,parent):
        wx.Panel.__init__(self, parent)
        self.f, self.a = plt.subplots(1,1)
        self.f.subplots_adjust(right=0.92, left=0.14)
        self.a.minorticks_on()
        self.a.set_xlabel(r'$(t - t_p)/P$')
        self.a.set_ylabel(r'$\left(F_p/F_{\!\star}\right) \times\ 10^{-9}$')
        self.canvas = FigureCanvas(self,-1, self.f)
     
    def plot(self, omega, ecc, incl, argperi, Omega, Rp, semi, models, delta):
        from LightCurve import make_lc
        self.a.cla()
        make_lc(self.a, omega, ecc, incl, argperi, Omega, Rp, semi, models, delta)
        self.canvas.draw()    

    def save(self):
        """ save figure image to file"""
        import os
        file_choices = "PDF (*.pdf)|*.pdf|" \
                       "PNG (*.png)|*.png|" \
                       "PS (*.ps)|*.ps|" \
                       "EPS (*.eps)|*.eps|" \
                       "BMP (*.bmp)|*.bmp"
        thisdir  = os.getcwd()
        dlg = wx.FileDialog(self, message='Save Plot Figure as...',
                            defaultDir = thisdir, defaultFile='plot.pdf',
                            wildcard=file_choices, style=wx.SAVE)
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            self.canvas.print_figure(path,dpi=300)
            if (path.find(thisdir) ==  0):
                path = path[len(thisdir)+1:]
 
class MyFrame(wx.Frame):
    def __init__(self,parent,title):
        wx.Frame.__init__(self, parent, title=title, size=(940,520), style=wx.MINIMIZE_BOX|wx.SYSTEM_MENU|
                  wx.CAPTION|wx.CLOSE_BOX|wx.CLIP_CHILDREN)
        self.sp = wx.SplitterWindow(self)
        self.p1 = p1(self.sp)
        self.p2 = wx.Panel(self.sp,style=wx.SUNKEN_BORDER)
        self.sp.SplitVertically(self.p1, self.p2, 625)
         
        # orbital params *****
        leftnum = 220
        base = 20.5
        dnum = 30
        orbparams=wx.StaticText(self.p2, -1, label='orbital parameters',pos=(75, base))
        title_font(orbparams, 'blue', 15)
        wx.StaticText(self.p2, -1, label='semi-major axis (AU)', pos=(10, base+dnum))
        self.semi_str = wx.TextCtrl(self.p2, -1, "5.0", size=(55,20), pos=(leftnum,base+dnum-1.5))
        wx.StaticText(self.p2, -1, label='eccentricity', pos=(10, base+2*dnum))
        self.ecc_str = wx.TextCtrl(self.p2, -1, "0.3", size=(55, 20), pos=(leftnum,base+2*dnum-1.5))
        wx.StaticText(self.p2, -1, label='inclination (deg)', pos=(10, base+3*dnum))
        self.incl_str = wx.TextCtrl(self.p2, -1, "90.0", size=(55, 20), pos=(leftnum,base+3*dnum-1.5))
        wx.StaticText(self.p2, -1, label='arg. of pericenter (deg)', size=(55,20),pos=(10, base+4*dnum))
        self.argperi_str = wx.TextCtrl(self.p2, -1, '0.0', size=(55,20),pos=(leftnum,base+4*dnum-1.5))
        #wx.StaticText(self.p2, -1, label='long. of ascending node (deg)',size=(55,20), pos=(10, 161.5))
        #self.Omega_str = wx.TextCtrl(self.p2, -1, '90.0', size=(55,20),pos=(leftnum,160))
        # *****

        # planet params *****
        basenum = 185
        dextra = 1.5
        planparams = wx.StaticText(self.p2, -1, label='planet parameters', pos=(72, basenum))
        title_font(planparams, 'blue', 15)
        wx.StaticText(self.p2, -1, label='single scattering albedo', pos=(10, basenum+dnum+dextra))
        self.omega_str = wx.TextCtrl(self.p2,-1, "1.0", size=(55,20), pos=(leftnum, basenum+dnum))
        wx.StaticText(self.p2, -1, label='radius (Jupiter radii)', pos=(10, basenum+2*dnum+dextra))
        self.Rp_str = wx.TextCtrl(self.p2, -1, "1.0", size=(55, 20), pos=(leftnum, basenum+2*dnum))
        # *****

        # check boxex *****
        baseleft = 26
        bump = 10
        self.lambox = wx.CheckBox(self.p2, -1, 'Lambert', pos=(baseleft, basenum+3.5*dnum+bump))
        self.lambox.SetValue(True)
        self.isobox = wx.CheckBox(self.p2, -1, 'isotropic', pos=(baseleft, basenum+4.5*dnum+bump))
        self.isobox.SetValue(True)
        self.raybox = wx.CheckBox(self.p2, -1, 'scalar Rayleigh', pos=(baseleft+100, basenum+3.5*dnum+bump))
        self.raybox.SetValue(True)
        self.vecbox = wx.CheckBox(self.p2, -1, 'vector Rayleigh', pos=(baseleft+100, basenum+4.5*dnum+bump))
        self.vecbox.SetValue(True)
        self.anibox = wx.CheckBox(self.p2, -1, 'anisotropic:   '+u'\N{GREEK SMALL LETTER DELTA} =', pos=(baseleft, basenum+5.5*dnum+bump))
        self.anibox.SetValue(True)
        self.anispin = FS.FloatSpin(self.p2, -1, pos=(baseleft+138, basenum+5.5*dnum+bump), min_val=0.1, max_val=1.0,
                                    increment=0.1, value=0.1, agwStyle=FS.FS_LEFT)
        self.anispin.SetFormat("%f")
        self.anispin.SetDigits(1)
        # *****

        # plot default params
        self.plot(wx.EVT_BUTTON)

        # buttons *****
        self.plotbut = wx.Button(self.p2, -1, "plot", size=(110,20), pos=(22, 425))
        self.plotbut.Bind(wx.EVT_BUTTON, self.plot)
        self.savebut = wx.Button(self.p2, -1, "save figure", size=(110,20), pos=(152, 425))
        self.savebut.Bind(wx.EVT_BUTTON, self.savefig)
        # *****

    def savefig(self, event):
        self.p1.save() 
        
    def plot(self,event):
        models = [self.lambox.GetValue(), self.isobox.GetValue(),\
                  self.raybox.GetValue(), self.vecbox.GetValue(), self.anibox.GetValue()]
        omega = float(self.omega_str.GetValue())
        ecc = float(self.ecc_str.GetValue())
        incl = float(self.incl_str.GetValue())
        argperi = float(self.argperi_str.GetValue())
        #Omega = float(self.Omega_str.GetValue())
        Omega = 90.0
        Rp = float(self.Rp_str.GetValue())
        semi = float(self.semi_str.GetValue())
        delta = self.anispin.GetValue()
        self.p1.plot(omega, ecc, incl, argperi, Omega, Rp, semi, models, delta)       

if __name__ == '__main__':
    app = wx.App(redirect=False)
    frame = MyFrame(None,"Light Curves")
    frame.Show()
    app.MainLoop()
