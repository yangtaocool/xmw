import sys

from java.awt import *
from java.io import *
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
from util import *

#pngDir = None
pngDir = "../../png/"
seismicDir = "../../data/f3d/"
s1 = Sampling(155,1.0,0)
s2 = Sampling(951,1.0,0)
s3 = Sampling(550,1.0,0)
n1,n2,n3 = s1.count,s2.count,s3.count
d1,d2,d3 = s1.delta,s2.delta,s3.delta
#k1,k2,k3 = 88,60,160; azimuth=285; elevation=11 # for 3D view of all horizons
k1,k2,k3 = 154,950,540; azimuth=240; elevation=25 # for 3D view of strips
fmin,fmax = -5.5,5.5
k1f,k2f,k3f = 65,406,114
k1f,k2f,k3f = 48,406,114
k1f,k2f,k3f = 48,406,0
gmin,gmax,gint,glab = -2.0,2.0,0.5,"Amplitude"
background = Color.WHITE
def main(args):
  slopes()
  #oneControlPoint()
  multipleControlPoints()
def slopes():
  f = readImage("f3dm")
  p2 = copy(f)
  p3 = copy(f)
  ep = copy(f)
  #good to use sigma2>sigma1 for extracting sequence boundaries
  sigma1,sigma2=2.0,1.0
  lsf = LocalSlopeFinder(sigma1,sigma2,20) 
  lsf.findSlopes(f,p2,p3,ep);
  writeImage("f3dmp2",p2)
  writeImage("f3dmp3",p3)
  writeImage("f3dmep",ep)
  for g in [p2,p3,ep]:
    world = World()
    addImage2ToWorld(world,f,g)
    makeFrame(world)

def oneControlPoint():
  k11 = [100.0]
  k12 = [335.0]
  k13 = [433.0]
  f  = readImage("f3dm")
  p2 = readImage("f3dmp2")
  p3 = readImage("f3dmp3")
  ep = readImage("f3dmep")
  ep = pow(ep,6.0) 
  horizonExtraction(f,p2,p3,ep,k11,k12,k13,"surf1")
  pg = setPointGroup(k11,k12,k13,9.0)
  displayHorizon(pg,"surf1")
  
def multipleControlPoints():
  k11 = [100, 43, 35, 91, 39, 38, 82, 76, 47, 76, 86, 57, 39, 37, 35,106, 58,101, 39,  6]
  k12 = [335,706,832,624,945,920,620,620,650,640,635,519,875,821,950,370,556,365,768,940]
  k13 = [433,200,495,  0,353,  9, 95,165,286,120, 22,547, 26,150,168,280,500,380,200,530]
  f = readImage("f3dm")
  p2 = readImage("f3dmp2")
  p3 = readImage("f3dmp3")
  ep = readImage("f3dmep")
  wp = pow(ep,20.0) 
  horizonExtraction(f,p2,p3,wp,k11,k12,k13,"surfm")
  pg = setPointGroup(k11,k12,k13,12.0)
  displayHorizon(pg,"surfm")

def horizonExtraction(f,p2,p3,wp,k11,k12,k13,filename):
  lmt = n1-1
  se = SurfaceExtractorC()
  #k11=se.refineConstraints(k11,k12,k13,f)
  se.setCG(0.01,200)
  se.setWeights(0.0)
  se.setSmoothings(4.0,4.0)
  se.setCG(0.01,100)
  surf = se.surfaceInitialization(n2,n3,lmt,k11,k12,k13)
  se.surfaceUpdateFromSlopes(wp,p2,p3,k11,k12,k13,surf)
  writeImage(filename,surf) 

def display(filename):
  f = readImage(filename)
  world = World()
  ipg = addImageToWorld(world,f,cmap=gray)
  ipg.setSlices(k1,k2,k3)
  makeFrame(world)

def setPointGroup(k1,k2,k3,size):
  np  = len(k1)
  xyz = zerofloat(np*3)
  rgb = zerofloat(np*3)
  ki = 0
  for i in range(np):
    xyz[ki  ] = k3[i]
    xyz[ki+1] = k2[i]
    xyz[ki+2] = k1[i]
    rgb[ki  ]  = 0#1/225 
    rgb[ki+1]  = 1#225/225 
    rgb[ki+2]  = 0#1/225 
    ki = ki+3
  pg = PointGroup(size,xyz,rgb);
  states = StateSet();
  cs = ColorState();
  cs.setColor(Color.GREEN);
  states.add(cs);
  lms = LightModelState();
  lms.setTwoSide(True);
  states.add(lms);
  ms = MaterialState();
  ms.setColorMaterial(GL_AMBIENT_AND_DIFFUSE);
  ms.setShininess(100.0);
  states.add(ms);
  pg.setStates(states);
  return pg;

def rgbFromHeight(h,r,g,b):
  n1 = len(h[0])
  n2 = len(h)
  ht = zerofloat(n1*n2)
  mp = ColorMap(-max(h),-min(h),ColorMap.JET)
  i = 0
  for i1 in range(n1):
    for i2 in range(n2):
      ht[i] = -h[i2][i1]
      i=i+1
  htRGB = mp.getRgbFloats(ht)
  i = 0
  for i1 in range(n1):
    for i2 in range(n2):
      r[i2][i1] = htRGB[i  ] 
      g[i2][i1] = htRGB[i+1] 
      b[i2][i1] = htRGB[i+2] 
      i = i+3

def rgbFromAmplitude(f,h,r,g,b):
  amp = zerofloat(n2*n3)
  si = SincInterpolator()
  #si.setUniform(n1,1,0,n2,1,0,n3,1,0,f)
  i = 0
  for i3 in range(n3):
    for i2 in range(n2):
      amp[i] = si.interpolate(n1,1,0,n2,1,0,n3,1,0,f,h[i3][i2],i2,i3)
      i = i+1
  aMin,aMax = -2.5,2.5
  mp = ColorMap(aMin,aMax,ColorMap.RED_WHITE_BLUE)
  ampRGB = mp.getRgbFloats(amp)
  i = 0
  for i3 in range(n3):
    for i2 in range(n2):
      r[i3][i2] = ampRGB[i  ] 
      g[i3][i2] = ampRGB[i+1] 
      b[i3][i2] = ampRGB[i+2] 
      i = i+3

def displayHorizon(pg,filename):
  f = readImage("f3dm")
  f = gain(f)
  h = readImage2(filename) 
  h1 = readImage2("hc") 
  r = zerofloat(n2,n3)
  g = zerofloat(n2,n3)
  b = zerofloat(n2,n3)
  #rgbFromHeight(h,r,g,b)
  rgbFromAmplitude(f,h,r,g,b)
  world = World()
  print min(f)
  print max(f)
  ipg = addImageToWorld(world,f,cmap=rwb,cmin=-2.5,cmax=2.5)
  ipg.setSlices(k1,k2,k3)
  tg  = TriangleGroup(True,s3,s2,add(h,0.0),r,g,b)
  world.addChild(pg)
  world.addChild(tg)
  makeFrame(world,png=filename)
def gain(x):
  g = mul(x,x) 
  ref = RecursiveExponentialFilter(20.0)
  ref.apply1(g,g)
  y = zerofloat(n1,n2,n3)
  div(x,sqrt(g),y)
  return y

#############################################################################
# read/write files
def readImage2(name):
  fileName = seismicDir+name+".dat"
  image = zerofloat(n2,n3)
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()
  return image
  
def readImage(name):
  fileName = seismicDir+name+".dat"
  n1,n2,n3 = s1.count,s2.count,s3.count
  #n1,n2,n3 = 120,480,430
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

def readSlice3(name):
  fileName = seismicDir+name+".dat"
  n1,n2 = s1.count,s2.count
  image = zerofloat(n1,n2)
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()
  return image

def slice12(k3,f):
  n1,n2,n3 = len(f[0][0]),len(f[0]),len(f)
  s = zerofloat(n1,n2)
  SimpleFloat3(f).get12(n1,n2,0,0,k3,s)
  return s

def slice13(k2,f):
  n1,n2,n3 = len(f[0][0]),len(f[0]),len(f)
  s = zerofloat(n1,n3)
  SimpleFloat3(f).get13(n1,n3,0,k2,0,s)
  return s

def slice23(k1,f):
  n1,n2,n3 = len(f[0][0]),len(f[0]),len(f)
  s = zerofloat(n2,n3)
  SimpleFloat3(f).get23(n2,n3,k1,0,0,s)
  return s


#############################################################################
# graphics

gray = ColorMap.GRAY
jet = ColorMap.JET
bwr = ColorMap.BLUE_WHITE_RED
rwb = ColorMap.RED_WHITE_BLUE

def addImageToWorld(world,image,cmap=gray,cmin=0,cmax=0):
  ipg = ImagePanelGroup(s1,s2,s3,image)
  ipg.setColorModel(cmap)
  if cmin<cmax:
    ipg.setClips(cmin,cmax)
  world.addChild(ipg)
  return ipg

def addImage2ToWorld(world,image1,image2):
  ipg = ImagePanelGroup2(s1,s2,s3,image1,image2)
  ipg.setColorModel1(ColorMap.getGray())
  ipg.setColorModel2(ColorMap.getJet())
  #ipg.setColorModel2(ColorMap.getHue(0.0,20.0,0.3))
  world.addChild(ipg)
  return ipg

def makeFrame(world,png=None):
  n1,n2,n3 = s1.count,s2.count,s3.count
  d1,d2,d3 = s1.delta,s2.delta,s3.delta
  f1,f2,f3 = s1.first,s2.first,s3.first
  l1,l2,l3 = s1.last,s2.last,s3.last
  #lightPosition=[-0.18,-0.4,0.8,0.0] # good for horizons 1 and 2
  lightPosition=[0.2,0.0,0.30,0.0] # good for horizons 1 and 2
  #lightPosition=[0.,0.,1.0,0.0] #default position
  frame = SimpleFrame(world)
  view = frame.getOrbitView()
  view.setLightPosition(lightPosition)
  zscale = 0.5*max(n2*d2,n3*d3)/(n1*d1)
  view.setAxesScale(1.0,1.0,zscale)
  view.setScale(1.1)
  view.setAzimuth(azimuth)
  view.setElevation(elevation)
  view.setWorldSphere(BoundingSphere(BoundingBox(f3,f2,f1,l3,l2,l1)))
  #frame.viewCanvas.setBackground(frame.getBackground())
  frame.setSize(1000,900)
  frame.setVisible(True)
  if png:
    frame.paintToFile(png+".png")
  return frame
"""
  frame = SimpleFrame(world)
  view = frame.getOrbitView()
  zscale = 0.75*max(n2*d2,n3*d3)/(n1*d1)
  view.setAxesScale(1.0,1.0,zscale)
  view.setScale(1.3)
  #view.setAzimuth(75.0)
  #view.setAzimuth(-75.0)
  view.setAzimuth(-65.0)
"""
 
def display2(k3,s,png=None):
  pp = PlotPanel(1,1,
    PlotPanel.Orientation.X1DOWN_X2RIGHT,
    PlotPanel.AxesPlacement.LEFT_BOTTOM)
  pp.setHInterval(10.0)
  pp.setVInterval(10.0)
  pp.setHLabel("Crossline (km)")
  pp.setVLabel("Inline (km)")
  pv = pp.addPixels(s1,s2,slice12(k3,s))
  pv.setClips(fmin,fmax)
  pf = PlotFrame(pp)
  pf.setFontSizeForPrint(6.0,480)
  pf.setSize(926,510)
  pf.setVisible(True)
  if png and pngDir:
    pf.paintToPng(300,6,pngDir+png+".png")
 
def display3(s,c=None,clabel="",cmin=0,cmax=0,png=None):
  pp = PlotPanelPixels3(
    PlotPanelPixels3.Orientation.X1DOWN_X2RIGHT,
    PlotPanelPixels3.AxesPlacement.LEFT_BOTTOM,
    s1,s2,s3,s)
  pp.setSlices(k1,k2,k3)
  pp.setLabel1("Time (s)")
  pp.setLabel2("Crossline (km)")
  pp.setLabel3("Inline (km)")
  pp.setClips(fmin,fmax)
  if c:
    cb = pp.addColorBar(clabel)
    #cb.setInterval(1.0)
    pp.setColorBarWidthMinimum(140)
    pp.setLineColor(Color.BLACK)
  else:
    pp.setLineColor(Color.YELLOW)
    #cb = pp.addColorBar("Amplitude")
    #cb.setInterval(5.0)
  pp.setInterval1(0.5)
  pp.setInterval2(2.0)
  pp.setInterval3(2.0)
  pp.mosaic.setHeightElastic(0,100)
  pp.mosaic.setHeightElastic(1,200)
  if c:
    pv12 = PixelsView(s1,s2,slice12(k3,c))
    pv12.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT)
    pv13 = PixelsView(s1,s3,slice13(k2,c))
    pv13.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT)
    pv23 = PixelsView(s2,s3,slice23(k1,c))
    pv23.setOrientation(PixelsView.Orientation.X1RIGHT_X2UP)
    for pv in [pv12,pv13,pv23]:
      pv.setColorModel(ColorMap.getJet(0.5))
      if cmin!=cmax:
        pv.setClips(cmin,cmax)
    pp.pixelsView12.tile.addTiledView(pv12)
    pp.pixelsView13.tile.addTiledView(pv13)
    pp.pixelsView23.tile.addTiledView(pv23)
  pf = PlotFrame(pp)
  pf.setFontSizeForSlide(1.0,1.0)
  if c:
    pf.setSize(1036,814)
  else:
    pf.setSize(859,814)
  pf.setVisible(True)
  if png and pngDir:
    pf.paintToPng(300,6,pngDir+png+".png")

def plot3f(g,a=None,amin=None,amax=None,amap=None,alab=None,aint=None,
           png=1):
  pp = PlotPanelPixels3(
    PlotPanelPixels3.Orientation.X1DOWN_X3RIGHT,
    PlotPanelPixels3.AxesPlacement.LEFT_BOTTOM,
    s1,s2,s3,g)
  pp.setSlices(k1f,k2f,k3f)
  pp.setLabel1("Time (s)")
  pp.setLabel2("Inline (km)")
  pp.setLabel3("Crossline (km)")
  pp.mosaic.setHeightElastic(0,180)
  pp.mosaic.setHeightElastic(1, 70)
  pp.setClips(gmin,gmax)
  pp.setColorModel(gray)
  if a:
    pp.setLineColor(Color.WHITE)
    cb = pp.addColorBar(alab)
    if aint:
      cb.setInterval(aint)
  else:
    pp.setLineColor(Color.YELLOW)
    cb = pp.addColorBar("Amplitude")
    cb.setInterval(0.5)
  pp.setInterval1(0.1)
  pp.setInterval2(2.0)
  pp.setInterval3(2.0)
  if a:
    pv12 = PixelsView(s1,s2,slice12(k3,a))
    pv12.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT)
    pv12.setInterpolation(PixelsView.Interpolation.NEAREST)
    pv13 = PixelsView(s1,s3,slice13(k2,a))
    pv13.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT)
    pv13.setInterpolation(PixelsView.Interpolation.NEAREST)
    pv23 = PixelsView(s2,s3,slice23(k1f,a))
    pv23.setOrientation(PixelsView.Orientation.X1RIGHT_X2UP)
    pv23.setInterpolation(PixelsView.Interpolation.NEAREST)
    for pv in [pv12,pv13,pv23]:
      pv.setColorModel(ColorMap.getJet(0.5))
      if amin!=amax:
        pv.setClips(amin,amax)
    pp.pixelsView12.tile.addTiledView(pv12)
    pp.pixelsView13.tile.addTiledView(pv13)
    pp.pixelsView23.tile.addTiledView(pv23)
  pf = PlotFrame(pp)
  pf.setBackground(background)
  pp.setColorBarWidthMinimum(170)
  #pf.setFontSizeForSlide(1.0,0.8)
  pf.setSize(1200,700)
  pf.setVisible(True)
  if png and pngDir:
    png = pngDir+"f"+str(k1f)
    pf.paintToPng(360,7.0,png+".png")

def plotFrame(s1,s2,f,h,i3t):
  orient = PlotPanel.Orientation.X1DOWN_X2RIGHT
  panel  = PlotPanel(2,1,orient)
  pxv    = panel.addPixels(0,0,s1,s2,f)
  pxv.setColorModel(ColorMap.GRAY)
  pxv    = panel.addPixels(1,0,s1,s2,f)
  pxv.setColorModel(ColorMap.GRAY)
  ptv1 = panel.addPoints(0,0,h[0],h[2])
  ptv2 = panel.addPoints(0,0,h[1],h[2])
  ptv1.setStyle("b-")
  ptv2.setStyle("g-")
  panel.setTitle("section "+i3t)
  panel.addColorBar()
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(PlotFrame.EXIT_ON_CLOSE);
  frame.setVisible(True)
#############################################################################
# Run the function main on the Swing thread
import sys
class _RunMain(Runnable):
  def __init__(self,main):
    self.main = main
  def run(self):
    self.main(sys.argv)
def run(main):
  SwingUtilities.invokeLater(_RunMain(main)) 
run(main)
