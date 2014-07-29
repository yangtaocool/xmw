import sys
from math import *
from java.awt import *
from java.lang import *
from java.util import *
from javax.swing import *
from java.util.Random import *

from edu.mines.jtk.la import *
from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.interp import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.ogl.Gl import *
from edu.mines.jtk.util import *
from edu.mines.jtk.util.ArrayMath import *

from sf import *

#############################################################################
# parameters
n = 101
d = 0.05
f = 0.00
s = Sampling(n,1.0,1.0)
pngDir = "./png/"
#############################################################################

def main(args):
  x = zerofloat(n)
  y1 = zerofloat(n)
  y2 = zerofloat(n)
  y3 = zerofloat(n)
  x[50]=1.0
  ef = ExponentialFilter(6.0) 
  ef.applyTwoSide(x,y1)
  ef.applyOneSideCausal(x,y2)
  ef.applyOneSideAnticausal(x,y3)
  plot1D(s,x, "x","y",color=Color.black,title="input")
  plot1D(s,y1,"x","y",color=Color.red,title="towSide")
  plot1D(s,y2,"x","y",color=Color.blue,title="oneSideC")
  plot1D(s,y3,"x","y",color=Color.green,title="oneSideA")
##################################################################
# plots
jet = ColorMap.JET
gray = ColorMap.GRAY
def plot1D(s,x1,vlabel,hlabel,x2=None,color=Color.black,title=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  pv = sp.addPoints(s,x1)
  if x2:
   pv1 = sp.addPoints(s,x1)
   pv2 = sp.addPoints(s,x2)
   pv1.setLineColor(Color.blue)
   pv2.setLineColor(Color.red)
   pv1.setLineWidth(2.0)
   pv2.setLineWidth(2.0)
  else:
    pv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
    pv.setMarkColor(color)
    #pv.setLineColor(Color.magenta)
    pv.setLineColor(color)
    pv.setLineWidth(2.0)
    pv.setMarkSize(8.0)
  sp.setSize(400,1000)
  sp.setVLabel(vlabel)
  sp.setHLabel(hlabel)
  if title:
    sp.paintToPng(720,3.3,pngDir+title+".png")
#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
