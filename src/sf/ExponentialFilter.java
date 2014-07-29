package sf;

import java.io.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import edu.mines.jtk.interp.*;
import static edu.mines.jtk.util.ArrayMath.*;

import vec.*;
import util.*;

/**
 * Two-side and one-side recursive exponential filters 
 * @author Xinming Wu and Dave Hale
 * @version 2014.03.12
 */

public class ExponentialFilter {
  public ExponentialFilter(
    double sigma)
  {
    Check.argument(sigma>=0.0,"sigma is non-negative");
    _sigma = (float)sigma;
    _a = aFromSigma(sigma);
  }
 
  public void applyTwoSide(float[] x, float[] y) {
   int n = x.length;
   float b = 1.0f-_a;
   float sx = 1.0f, sy = _a;
   float yi = 0.0f;
   y[0] = yi = sy*yi+sx*x[0];
   for (int i=1; i<n-1; ++i)
     y[i] = yi = _a*yi+b*x[i];
   sx /= 1.0f+_a; sy /= 1.0f+_a;
   y[n-1] = yi = sy*yi+sx*x[n-1];
   for (int i=n-2; i>=0; --i)
    y[i] = yi = _a*yi+b*y[i]; 
  }

  public void applyCausal(float[] x, float[] y) {
   int n = x.length;
   float b = 1.0f-_a;
   float yi = y[0] = x[0];
   for (int i=1; i<n; ++i)
     y[i] = yi = _a*yi+b*x[i]; 
  }
  
  public void applyAnticausal(float[] x, float[] y) {
   int n = x.length;
   float b = 1.0f-_a;
   float yi = y[n-1] = x[n-1];
   for (int i=n-2; i>=0; --i)
     y[i] = yi = _a*yi+b*x[i]; 
  }

  private static float aFromSigma(double sigma) {
    if (sigma<=0.0f)
      return 0.0f;
    double ss = sigma*sigma;
    return (float)((1.0+ss-sqrt(1.0+2.0*ss))/ss);
  }

  private float _sigma,_a; 
}
