import sys

from java.awt import *
from java.awt.image import *
from java.io import *
from java.nio import *
from java.lang import *
from javax.swing import *
from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.ogl.Gl import *
from edu.mines.jtk.sgl import *
from edu.mines.jtk.util import *
from edu.mines.jtk.util.ArrayMath import *

from hv import *
from unc import *
from util import *

#############################################################################
# parameters
# fake image
# real image
n1,n2= 120,260
d1,d2= 0.004, 0.025
f1,f2= 110.0*d1,10.0*d2
s1 = Sampling(n1,d1,f1)
s2 = Sampling(n2,d2,f2)
pngDir = "../../../png/"
seismicDir = "../../data/unc/"
#############################################################################

def main(args):
  f = readImage("f3d148")
  u1  = zerofloat(n1,n2)
  u2  = zerofloat(n1,n2)
  u1m = zerofloat(n1,n2)
  u2m = zerofloat(n1,n2)
  us = fillfloat(0.0,n1,n2)
  ut1 = fillfloat(1.0,n1,n2)
  ut2 = fillfloat(1.0,n1,n2)
  unc(f,ut1,ut2,us)
  f = gain(f)
  unweightedLof(f,u1,u2)
  weightedLof(f,ut1,ut2,u1m,u2m)
  correctNormal(u1m,u2m,ut1)
  plot(u1,"u1",cbi=0.04,cbar="u1",cmap=jet,cmin=min(u1m),cmax=max(u1m),interp=False)
  plot(u1m,"u1m",cbi=0.04,cbar="u1",cmap=jet,cmin=min(u1m),cmax=max(u1m),interp=False)
  plot(u2,"u2",cbi=0.1,cbar="u2",cmap=jet,cmin=min(u2m),cmax=max(u2m),interp=False)
  plot(u2m,"u2m",cbi=0.1,cbar="u2",cmap=jet,cmin=min(u2m),cmax=max(u2m),interp=False)
  p =mul(div(u2m,u1m),-1.0) 
  plot(p,"slopem",cbi=0.2,cbar="Slope",cmap=jet,cmin=min(p),cmax=max(p),interp=False)
  
  #plotVectors2(f,u1,u2,u1m,u2m,"normalVectors")
  #flatten(f,us,ut2,u1m,u2m)
def unc(f,ut1,ut2,us):
  u  = zerofloat(n1,n2)
  uc1 = zerofloat(n1,n2)
  uc2 = zerofloat(n1,n2)
  ua1 = zerofloat(n1,n2)
  ua2 = zerofloat(n1,n2)
  ud = Unconformity()
  lof = LocalOrientFilter(1.0,2.0)
  et = lof.applyForTensors(f)
  #lsf = LaplacianSmoother(2.0)
  #fs = zerofloat(n1,n2)
  #lsf.apply(LaplacianSmoother.Direction2.V,et,f,fs)
  ud.setForLofu(0.9,n2*4,0.0)
  ud.uncLikelihood(et,f,u)
  ud.applyForNormals(et,f,uc1,uc2,ua1,ua2)
  #ref = RecursiveExponentialFilter(2.0)
  #ref.apply(u,us)
  copy(u,us)
  sub(us,min(us),us)
  div(us,max(us),us)
  pow(us,1.5,us)
  #plot2(s1,s2,f,us,png="Vector difference")
  f = gain(f)
  ud.thin(1,us,ut1,ut2) 
  plot2(s1,s2,f,g=sub(1.0,ut1),png="thinnedUnconformityAttribute1")
  plot2(s1,s2,f,g=sub(1.0,ut2),png="thinnedUnconformityAttribute2")
  plot2(s1,s2,f,g=sub(1.0,us),png="UnconformityAttribute")
  #plotVectors2(f,uc1,uc2,ua1,ua2,"uc(causal,red) & ua(anticausal,blue): a=0.9")
  #plot2(s1,s2,f)
def unweightedLof(f,u1,u2):
  el = zerofloat(n1,n2)
  w1 = zerofloat(n1,n2)
  lof = LocalOrientFilter(12.0,6.0)
  lof.applyForNormalLinear(f,u1,u2,el)
  el = pow(el,8.0)
  p2 = div(u2,u1)
  p2 = mul(p2,-1.0)
  fl2 = Flattener2()
  mp = fl2.getMappingsFromSlopes(s1,s2,p2,el,el,w1)
  g = mp.flatten(f)
  rgt = mp.u1
  plot(g,"flattened",cbi=1.0,cbar="Amplitude",cmap=gray,cmin=-2.0,cmax=2.0,interp=False)
  plot(rgt,"RGT",cbi=0.1,cbar="RGT",cmap=jet,interp=True)
  plot(p2,"p2",cbi=0.1,cbar="Slopes",cmap=jet,cmin=-0.48,cmax=0.18,interp=False)
  #plotVectors1(f,u1,u2,"normalVectors")

def weightedLof(f,ut1,ut2,u1,u2):
  el = zerofloat(n1,n2)
  sg = pow(ut1,2)
  sgp= pow(ut2,8)
  wlof = WeightedLocalOrientFilter(6.0,2.0)
  wlof.setScaleFactor(sg,sgp)
  wlof.setGradientSmoothing(2,2)
  wlof.applyForNormalLinear(f,u1,u2,el)
  correctNormal(u1,u2,ut1)
  et = wlof.applyForTensors(f)
  '''
  lof = LocalOrientFilter(4.0,2.0)
  lof.applyForNormalLinear(f,u1,u2,el)
  '''
  #plot(pow(ut2,4),"scale factors for LSFilter",cmap=jet)
  #plot(u1,"u1m",cbi=0.1,cbar="u1",cmap=jet,cmin=0.8,cmax=1.0,interp=False)
  #plot(u2,"u2m",cbi=0.2,cbar="u2",cmap=jet,cmin=-0.4,cmax=0.6,interp=False)
  #plotVectors1(f,u1,u2,"normalVectorsm")
  #flatten(f,ut2,u1,u2)

def flatten(f,us,ut,u1,u2):
  p = div(u2,u1)
  p = mul(p,-1.0)
  el = zerofloat(n1,n2)
  lof = LocalOrientFilter(2.0,8.0)
  lof.applyForNormalLinear(f,u1,u2,el)
  fl2 = Flattener2()
  fl2.setSmoothings(75,75)
  el = pow(el,9.0)
  w1 = zerofloat(n1,n2)
  #add(w1,0.01,w1)
  #mul(w1,el,w1)
  wp1 = pow(ut,80)
  mul(w1,wp1,w1)
  mul(el,wp1,wp1)
  plot(wp1,"w1",cmap=jet)
  mp = fl2.getMappingsFromSlopes(s1,s2,p,wp1,wp1,w1)
  ft = mp.flatten(f)
  rgt = mp.u1
  plot(f,"input",cbi=1.0,cbar="Amplitude",cmap=gray,cmin=-2.0,cmax=2.0,interp=False)
  plot(ft,"flattenedm",cbi=1.0,cbar="Amplitude",cmap=gray,cmin=-2.0,cmax=2.0,interp=False)
  plot(rgt,"RGTm",cbi=0.1,cbar="RGT",cmap=jet,interp=True)
  plot(p,"p2m",cbi=0.1,cbar="Slopes",cmap=jet,cmin=-0.48,cmax=0.18,interp=False)

def plotVectors1(f,u1,u2,title):
  np = int(n1*n2)
  x1 = zerofloat(2,np)
  x2 = zerofloat(2,np)
  i = 0
  dth = 8
  for i2 in range(1,n2-2,4):
    for i1 in range(20,n1-20,3):
      x1[i][0] = i1*d1+f1
      x1[i][1] = (i1+dth*u1[i2][i1])*d1+f1
      x2[i][0] = i2*d2+f2
      x2[i][1] = (i2+dth*u2[i2][i1])*d2+f2
      i = i+1
  x1 = copy(2,i,0,0,x1)
  x2 = copy(2,i,0,0,x2)
  orientation = PlotPanel.Orientation.X1DOWN_X2RIGHT;
  panel = PlotPanel(1,1,orientation);
  pxv = panel.addPixels(0,0,s1,s2,f);
  pxv.setColorModel(ColorMap.GRAY)
  #panel.setTitle(title)
  panel.setHLabel("Time (s)")
  panel.setVLabel("Crossline (km)")
  ptv = panel.addPoints(x1,x2)
  ptv.setLineColor(Color.RED)
  ptv.setLineWidth(2.0)
  panel.addColorBar();
  panel.setColorBarWidthMinimum(120)
  moc = panel.getMosaic();
  frame = PlotFrame(panel);
  frame.setDefaultCloseOperation(PlotFrame.EXIT_ON_CLOSE);
  frame.setFontSizeForSlide(1.0,0.9,16.0/9.0)
  frame.setVisible(True);
  frame.setSize(n2*4+120,500)

def plotVectors2(f,us1,us2,ua1,ua2,title):
  np = int(n1*n2)
  xs1 = zerofloat(2,np)
  xa1 = zerofloat(2,np)
  xs2 = zerofloat(2,np)
  xa2 = zerofloat(2,np)
  i = 0
  dth = 6
  for i2 in range(3,n2-2,5):
    for i1 in range(5,n1-6,4):
  #for i2 in range(3,n2-2,5):
    #for i1 in range(5,n1-6,10):
      xs1[i][0] = i1*d1+f1
      xa1[i][0] = i1*d1+f1
      xs1[i][1] = (i1+dth*us1[i2][i1])*d1+f1
      xa1[i][1] = (i1+dth*ua1[i2][i1])*d1+f1
      xs2[i][0] = i2*d2+f2
      xa2[i][0] = i2*d2+f2
      xs2[i][1] = (i2+dth*us2[i2][i1])*d2+f2
      xa2[i][1] = (i2+dth*ua2[i2][i1])*d2+f2
      i = i+1
  xs1 = copy(2,i,0,0,xs1)
  xa1 = copy(2,i,0,0,xa1)
  xs2 = copy(2,i,0,0,xs2)
  xa2 = copy(2,i,0,0,xa2)
  orientation = PlotPanel.Orientation.X1DOWN_X2RIGHT;
  panel = PlotPanel(1,1,orientation);
  pxv = panel.addPixels(0,0,s1,s2,f);
  pxv.setColorModel(ColorMap.GRAY)
  #panel.setTitle(title)
  panel.setHLabel("Crossline (km)")
  panel.setVLabel("Time (s)")
  panel.setHInterval(1.0)
  panel.setVInterval(0.1)
  ptvs = panel.addPoints(xs1,xs2)
  ptva = panel.addPoints(xa1,xa2)
  ptvs.setLineColor(Color.RED)
  ptva.setLineColor(Color.BLUE)
  ptva.setLineWidth(2.0)
  ptvs.setLineWidth(2.0)
  cb = panel.addColorBar();
  cb.setLabel("Amplitude")
  panel.setColorBarWidthMinimum(120)
  moc = panel.getMosaic();
  frame = PlotFrame(panel);
  frame.setDefaultCloseOperation(PlotFrame.EXIT_ON_CLOSE);
  frame.setVisible(True);
  frame.setSize(n2*4+120,500)
  frame.setFontSizeForSlide(1.0,0.9,16.0/9.0)
  frame.paintToPng(720,3.3,pngDir+title+".png")

def gain(x):
  g = mul(x,x) 
  ref = RecursiveExponentialFilter(10.0)
  ref.apply1(g,g)
  y = like(x)
  div(x,sqrt(g),y)
  return y

def like(x):
  n2 = len(x)
  n1 = len(x[0])
  return zerofloat(n1,n2)
def correctNormal(u1,u2,ut):
  ut = pow(ut,20)
  for i2 in range(n2):
    for i1 in range(n1):
      if ((ut[i2][i1]<0.2)):
        u1[i2][i1  ] = u1[i2][i1-2]
        u2[i2][i1  ] = u2[i2][i1-2]
        u1[i2][i1-1] = u1[i2][i1-2]
        u2[i2][i1-1] = u2[i2][i1-2]
        u1[i2][i1+1] = u1[i2][i1+2]
        u2[i2][i1+1] = u2[i2][i1+2]

#############################################################################
# plotting

def plot2(s1,s2,f,g=None,gmin=None,gmax=None,
          label=None,title=None,png=None):
  n1 = len(f[0])
  n2 = len(f)
  panel = panel2()
  if label:
    panel.addColorBar(label)
  else:
    cb = panel.addColorBar()
    cb.setInterval(0.2)
    cb.setLabel("Unconformity likelihood")
    #cb.setLabel("Vector difference")
  if title:
    panel.setTitle(title)
  panel.setColorBarWidthMinimum(120)
  panel.setVInterval(0.1)
  panel.setHInterval(1.0)
  panel.setHLabel("Crossline (km)")
  panel.setVLabel("Time (s)")
  pv = panel.addPixels(s1,s2,f)
  #pv = panel.addPixels(f)
  pv.setInterpolation(PixelsView.Interpolation.LINEAR)
  pv.setColorModel(ColorMap.GRAY)
  pv.setClips(-2.0,2.0)
  if g:
    alpha = 0.0
    pv = panel.addPixels(s1,s2,g)
    #pv = panel.addPixels(g)
    pv.setInterpolation(PixelsView.Interpolation.NEAREST)
    if gmin==None: gmin = min(g)
    if gmax==None: gmax = max(g)
    pv.setClips(gmin,gmax)
    pv.setColorModel(ColorMap.getJet(alpha))
    if gmin==0:
      updateColorModel(pv,1.0)
  frame2(panel,png)

def updateColorModel(pv,alpha):
    n = 256
    r = zerobyte(n)
    g = zerobyte(n)
    b = zerobyte(n)
    a = zerobyte(n)
    icm = pv.getColorModel()
    icm.getReds(r)
    icm.getGreens(g)
    icm.getBlues(b)
    for i in range(n):
      ai = int(255.0*alpha*i/n)
      if ai>127:
        ai -= 256
      a[i] = ai
    icm = IndexColorModel(8,n,r,g,b,a)
    pv.setColorModel(icm)

def panel2():
  #panel = PlotPanel(1,1,
  #  PlotPanel.Orientation.X1DOWN_X2RIGHT,
  #  PlotPanel.AxesPlacement.NONE)
  panel = PlotPanel(1,1,
    PlotPanel.Orientation.X1DOWN_X2RIGHT)
  return panel

def frame2(panel,png=None):
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setSize(n2*4+120,500);
  frame.setFontSizeForSlide(1.0,0.9,16.0/9.0)
  frame.setVisible(True)
  if png:
    frame.paintToPng(720,3.3,pngDir+png+".png")
  return frame




##################################################################
# plots
jet = ColorMap.JET
gray = ColorMap.GRAY
def plot(x,title,cbi=0.2,cbar=None,cmap=jet,px1=None,px2=None,cmin=0,cmax=0,interp=False):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  pv = sp.addPixels(s1,s2,x)
  if (interp):
    pv.setInterpolation(PixelsView.Interpolation.NEAREST);
  pv.setColorModel(cmap);
  if cmin<cmax:
    pv.setClips(cmin,cmax)
  if px1:
    ptv = sp.addPoints(px1,px2);
    ptv.setLineStyle(PointsView.Line.NONE);
    ptv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE);
    #ptv.setMarkColor(Color.red)
    ptv.setMarkSize(3.0);
  sp.setSize(n2*4+120,500);
  sp.setHInterval(1.0)
  sp.setVInterval(0.1)
  #sp.addTitle(title)
  cb=sp.addColorBar(cbar)
  cb.setWidthMinimum(120)
  cb.setInterval(cbi)
  cb.setLabel(cbar)
  sp.setHLabel("Crossline (km)")
  sp.setVLabel("Time (s)")
  sp.setFontSizeForSlide(1.0,0.9,16.0/9.0)
  sp.paintToPng(720,3.3,pngDir+title+".png")

def addImageToWorld(s1,s2,s3,world,image,cmap=gray,cmin=0.0,cmax=0.0):
  ipg = ImagePanelGroup(s1,s2,s3,image)
  ipg.setColorModel(cmap)
  if cmin<cmax:
    ipg.setClips(cmin,cmax)
  world.addChild(ipg)
  return ipg

def addImage2ToWorld(world,image1,image2):
  ipg = ImagePanelGroup2(s1,s2,s3,image1,image2)
  ipg.setColorModel1(ColorMap.getGray())
  #ipg.setColorModel2(ColorMap.getJet(0.3))
  ipg.setColorModel2(ColorMap.getHue(0.0,20.0,0.3))
  world.addChild(ipg)
  return ipg


#############################################################################
# read image
def readImage(name):
  fileName = seismicDir+name+".dat"
  n1,n2= s1.count,s2.count
  image = zerofloat(n1,n2)
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()
  return image
def writeImage(name,image):
  fileName = seismicDir+name+".dat"
  aos = ArrayOutputStream(fileName)
  aos.writeFloats(image)
  aos.close()
  return image

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
