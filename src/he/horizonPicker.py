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

from he import *

#############################################################################
# parameters
# fake image
# 3D real image
n1,n2,n3= 251,357,161
d1,d2,d3= 0.004,0.0250,0.0250
f1,f2,f3= 0.5,0.0,0.0
s1 = Sampling(n1,d1,f1)
s2 = Sampling(n2,d2,f2)
s3 = Sampling(n3,d3,f3)
pngDir = "./png/"
seismicDir = "../../data/"
#############################################################################

def main(args):
  applyFor2D()
  applyFor3D()

def applyFor3D():
  g = readImage("tpst")
  lof =  LocalSlopeFinder(4.0,2.0,10.0)
  p2 = zerofloat(n1,n2,n3) 
  p3 = zerofloat(n1,n2,n3) 
  ep = zerofloat(n1,n2,n3) 
  lof.findSlopes(g,p2,p3,ep) 
  p = zerofloat(n1,n2,n3) 
  t = zerofloat(n1,n2,n3) 
  h = zerofloat(n1,n2,n3) 
  hp = HorizonPicker()
  hp.setForPeaksAndTroughs(2)
  hp.setForHorizons(1,4,4,2,30)
  hp.peakFinder(g,ep,p)
  hp.troughFinder(g,ep,t)
  hp.applyForHorizons(g,t,h)
  world = World()
  addImageToWorld(s1,s2,s3,world,g)
  addImageToWorld(s1,s2,s3,world,p,cmap=jet)
  addImageToWorld(s1,s2,s3,world,t,cmap=jet)
  addImageToWorld(s1,s2,s3,world,h,cmap=jet)
  makeFrame(world)
def applyFor2D():
  g = readImage("tpst")
  f = zerofloat(n1,n2) 
  SimpleFloat3(g).get12(n1,n2,0,0,120,f)
  u1 = zerofloat(n1,n2)
  u2 = zerofloat(n1,n2)
  ep = zerofloat(n1,n2)
  lof = LocalOrientFilter(4.0,2.0)
  lof.applyForNormalLinear(f,u1,u2,ep)
  #f = gain(f)
  t = like(f)
  p = like(f)
  h = like(f)
  hp = HorizonPicker()
  hp.setForPeaksAndTroughs(2)
  hp.setForHorizons(4,4,8,30)
  hp.peakFinder(f,ep,p)
  hp.troughFinder(f,ep,t)
  hp.applyForHorizons(f,p,h)
  plot(f,"seismic",cbi=0.5,cbar="Amplitude",cmap=gray)
  plot(p,"peaks",cbi=0.5,cbar="Amplitude",cmap=jet,interp=True)
  plot(t,"troughs",cbi=0.5,cbar="Amplitude",cmap=jet,interp=True)
  plot(h,"peaks",cbi=10,cbar="Horizons",cmap=jet,interp=True)
  
def gain(x):
  g = mul(x,x) 
  ref = RecursiveExponentialFilter(5.0)
  ref.apply1(g,g)
  y = like(x)
  div(x,sqrt(g),y)
  return y

def like(x):
  n2 = len(x)
  n1 = len(x[0])
  return zerofloat(n1,n2)
#############################################################################
# plots
jet = ColorMap.JET
gray = ColorMap.GRAY
def plot(x,title,cbi=0.5,cbar=None,cmap=jet,px1=None,px2=None,cmin=0,cmax=0,interp=False):
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
  sp.setSize(n2*3+125,n1*3);
  sp.setHInterval(1.0)
  sp.setVInterval(0.1)
  #sp.addTitle(title)
  cb=sp.addColorBar(cbar)
  cb.setWidthMinimum(125)
  cb.setInterval(cbi)
  cb.setLabel(cbar)
  sp.setHLabel("Crossline (km)")
  sp.setVLabel("Time (s)")
  sp.setFontSize(36)
  #sp.paintToPng(720,3.3,pngDir+title+".png")

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

def makeFrame(world):
  n1,n2,n3 = s1.count,s2.count,s3.count
  d1,d2,d3 = s1.delta,s2.delta,s3.delta
  f1,f2,f3 = s1.first,s2.first,s3.first
  l1,l2,l3 = s1.last,s2.last,s3.last
  frame = SimpleFrame(world)
  view = frame.getOrbitView()
  zscale = 0.75*max(n2*d2,n3*d3)/(n1*d1)
  view.setAxesScale(1.0,1.0,zscale)
  view.setScale(1.1)
  #view.setAzimuth(azimuth)
  #view.setElevation(elevation)
  view.setWorldSphere(BoundingSphere(BoundingBox(f3,f2,f1,l3,l2,l1)))
  #frame.viewCanvas.setBackground(frame.getBackground())
  frame.setSize(1000,900)
  frame.setVisible(True)
  return frame

#############################################################################
# read image
def readImage(name):
  fileName = seismicDir+name+".dat"
  image = zerofloat(n1,n2,n3)
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
