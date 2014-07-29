package he;

import java.util.ArrayList;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.Parallel;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Connect peaks or troughs in a seismic image to form horizons.
 * These horizons are usually not complete, just patches, but they 
 * exactly align with amplitude peaks or troughs and can be used 
 * as constraints for other methods.
 * @author Xinming Wu
 * @version 2014.03.11
 */
public class HorizonPicker {
  public void setForPeaksAndTroughs(int dd) {
    _dd = dd;
  } 
  public void setForHorizons(int d1, int d2, int dm, int np) {
    _d1 = d1;
    _d2 = d2;
    _dm = dm;
    _np = np;
  } 
  public void setForHorizons(int d1, int d2, int d3, int dm, int np) {
    _d1 = d1;
    _d2 = d2;
    _d3 = d3;
    _dm = dm;
    _np = np;
  } 

  public void peakFinder(float[] f, float[] ep, float[] p) {
    int n = f.length;
    for (int i=_dd; i<n-_dd; ++i) {
      float fi = f[i];
      float epi = ep[i];
      float sum1 = 0.0f;
      float sum2 = 0.0f;
      for (int k=1; k<=_dd; ++k) {
        float dfmi = fi - f[i-k];
        float dfpi = fi - f[i+k];
        sum1 += dfmi+dfpi;
        sum2 += abs(dfmi)+abs(dfpi);
      }
      if (sum1==sum2 && epi>0.9f) {
        p[i] = f[i];
      } 
    } 
  }
  
  public void peakFinder(final float[][] f, final float[][] ep, final float[][] p) {
    final int n2 = f.length;
    Parallel.loop(n2,new Parallel.LoopInt() {
      public void compute(int i2) {
        peakFinder(f[i2],ep[i2],p[i2]);
      }
    });
  }

  public void peakFinder(final float[][][] f, final float[][][] ep, final float[][][] p) {
    final int n3 = f.length;
    Parallel.loop(n3,new Parallel.LoopInt() {
      public void compute(int i3) {
        peakFinder(f[i3],ep[i3],p[i3]);
      }
    });
  }

  public void troughFinder(float[] f, float[] ep, float[] p) {
    int n = f.length;
    for (int i=_dd; i<n-_dd; ++i) {
      float fi = f[i];
      float epi = ep[i];
      float sum1 = 0.0f;
      float sum2 = 0.0f;
      for (int k=1; k<=_dd; ++k) {
        float dfmi = -fi + f[i-k];
        float dfpi = -fi + f[i+k];
        sum1 += dfmi+dfpi;
        sum2 += abs(dfmi)+abs(dfpi);
      }
      if (sum1==sum2 && epi>0.9f) {
        p[i] = f[i];
      } 
    } 
  }

  public void troughFinder(final float[][] f, final float[][] ep, final float[][] p) {
    final int n2 = f.length;
    Parallel.loop(n2,new Parallel.LoopInt() {
      public void compute(int i2) {
        troughFinder(f[i2],ep[i2],p[i2]);
      }
    });
  }

  public void troughFinder(final float[][][] f, final float[][][] ep, final float[][][] p) {
    final int n3 = f.length;
    Parallel.loop(n3,new Parallel.LoopInt() {
      public void compute(int i3) {
        troughFinder(f[i3],ep[i3],p[i3]);
      }
    });
  }

  public void applyForHorizons(float[][] f, float[][] pt, float[][] h) {
    int n2 = f.length;
    int n1 = f[0].length;
    int hi = 0;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float pti = pt[i2][i1];
        if(pti !=0.0f) {
          float[] dp = new float[n2];
          int np = floodFill(i1,i2,f,pt,dp);
          if(np>_np) {
            hi++;
            horizonMarker(hi,h,dp);
          }
        } 
      }
    }
  }

  public void applyForHorizons(
    float[][][] f, float[][][] pt, float[][][] h) 
  {
    int n3 = f.length;
    int n2 = f[0].length;
    int n1 = f[0][0].length;
    int hi = 0;
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float pti = pt[i3][i2][i1];
          if(pti !=0.0f) {
            float[][] dp = new float[n3][n2];
            int np = floodFill(i1,i2,i3,f,pt,dp);
            if(np>_np) {
              hi++;
              horizonMarker(hi,h,dp);
            }
          } 
        }
      }
    }
  }

  private void horizonMarker(int hi, float[][][] h, float[][] dp) {
    int n3 = h.length;
    int n2 = h[0].length;
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        int i1 = (int)dp[i3][i2];
        if (i1 != 0) h[i3][i2][i1] = (float) hi;
      }
    }
  }
  private void horizonMarker(int hi, float[][] h, float[] dp) {
    int n2 = h.length;
    for (int i2=0; i2<n2; ++i2) {
        int i1 = (int)dp[i2];
        if (i1 != 0) h[i2][i1] = (float) hi;
    }
  }

  // pi1,pi2: indexes of the start point  
  // f[n2][n1]:  seismic image
  // pt[n2][n1]: peaks or troughs
  // dp[n2]:  depth of the horizon
  // return:  number of points on the horizons
  private int floodFill(int pi1, int pi2, float[][] f, float[][] pt, float[] dp) {
    int count = 0;
    ArrayList<int[]> pile = new ArrayList<int[]>();
    int[] id = {pi1,pi2};
    pile.add(id);
    int n1 = f[0].length;
    float fi = f[pi2][pi1];
    if (pi1==0 || pi1==n1-1) 
      dp[pi2] = (float)pi1; 
    else {
      float fm = f[pi2][pi1-1];
      float fp = f[pi2][pi1+1];
      dp[pi2] = parabolicPeak(pi1,fm,fi,fp);
    }
    count ++;
    pt[pi2][pi1] = 0.0f;
    while (pile.size()>0) {
      id = pile.remove(0);
      int i1 = id[0];
      int i2 = id[1];
      count = extend(i1,i2,count,pile,f,pt,dp); 
    } 
    return count;
  }
  private int floodFill(
    int pi1, int pi2, int pi3, float[][][] f, 
    float[][][] pt, float[][] dp) 
  {
    int count = 0;
    ArrayList<int[]> pile = new ArrayList<int[]>();
    int[] id = {pi1,pi2,pi3};
    pile.add(id);
    int n1 = f[0][0].length;
    float fi = f[pi3][pi2][pi1];
    if (pi1==0 || pi1==n1-1) 
      dp[pi3][pi2] = (float)pi1; 
    else {
      float fm = f[pi3][pi2][pi1-1];
      float fp = f[pi3][pi2][pi1+1];
      dp[pi3][pi2] = parabolicPeak(pi1,fm,fi,fp);
    }
    count ++;
    pt[pi3][pi2][pi1] = 0.0f;
    while (pile.size()>0) {
      id = pile.remove(0);
      int i1 = id[0];
      int i2 = id[1];
      int i3 = id[2];
      count = extend(i1,i2,i3,count,pile,f,pt,dp); 
    } 
    return count;
  }

  
  private int extend(
    int i1t, int i2t, int count, ArrayList<int[]> pile,
    float[][] f, float[][] pt, float[] dp)
  {
    int n2 = f.length;
    int n1 = f[0].length;
    int i1b = i1t-_d1; if(i1b<0) i1b=0;
    int i2b = i2t-_d2; if(i2b<0) i2b=0;
    int i1e = i1t+_d1; if(i1e>n1-1) i1e=n1-1;
    int i2e = i2t+_d2; if(i2e>n2-1) i2e=n2-1;
    ArrayList<int[]> id = new ArrayList<int[]>(); 
    int dMin = 10000;
    int[] ind = new int[2];
    for (int i2=i2b; i2<=i2e; ++i2) { 
      for (int i1=i1b; i1<=i1e; ++i1) { 
        float pti = pt[i2][i1];
        if(pti!= 0.0f) {
          int d1i = i1-i1t;
          int d2i = i2-i2t;
          int ddi = d1i*d1i+d2i*d2i;
          if(ddi<=1) {
            int[] indt = {i1,i2};
            id.add(indt);
          }
          if(ddi<dMin) {
            dMin = ddi;
            ind[0] = i1;
            ind[1] = i2;
          }
        }
      }
    }
    if(dMin<_dm) id.add(ind);
    int np = id.size();
    count += np;
    for (int i=0; i<np; ++i) {
      ind = id.get(i);
      int i1 = ind[0];
      int i2 = ind[1];
      pile.add(ind);
      pt[i2][i1] = 0.0f;
      float fi = f[i2][i1];
      if (i1==0 || i1==n1-1) 
        dp[i2] = (float)i1; 
      else {
        float fm = f[i2][i1-1];
        float fp = f[i2][i1+1];
        dp[i2] = parabolicPeak(i1,fm,fi,fp);
      }
    }
    return count;
  }  

  private int extend(
    int i1t, int i2t, int i3t, int count, ArrayList<int[]> pile,
    float[][][] f, float[][][] pt, float[][] dp)
  {
    int n3 = f.length;
    int n2 = f[0].length;
    int n1 = f[0][0].length;
    int i1b = i1t-_d1; if(i1b<0) i1b=0;
    int i2b = i2t-_d2; if(i2b<0) i2b=0;
    int i3b = i3t-_d3; if(i3b<0) i3b=0;
    int i1e = i1t+_d1; if(i1e>n1-1) i1e=n1-1;
    int i2e = i2t+_d2; if(i2e>n2-1) i2e=n2-1;
    int i3e = i3t+_d3; if(i3e>n3-1) i3e=n3-1;
    ArrayList<int[]> id = new ArrayList<int[]>(); 
    int dMin = 10000;
    int[] ind = new int[3];
    for (int i3=i3b; i3<=i3e; ++i3) { 
      for (int i2=i2b; i2<=i2e; ++i2) { 
        for (int i1=i1b; i1<=i1e; ++i1) { 
          float fi = f[i3][i2][i1];
          float pti = pt[i3][i2][i1];
          if(pti!= 0.0f) {
            int d1i = i1-i1t;
            int d2i = i2-i2t;
            int d3i = i3-i3t;
            int ddi = d1i*d1i+d2i*d2i+d3i*d3i;
            float ri = sqrt((float)ddi);
            if(ddi<=1) {
              int[] indt = {i1,i2,i3};
              id.add(indt);
            }
            if(ddi<dMin) {
              dMin = ddi;
              ind[0] = i1;
              ind[1] = i2;
              ind[2] = i3;
            }
          }
        }
      }
    }
    if(dMin<_dm) id.add(ind);
    int np = id.size();
    count += np;
    for (int i=0; i<np; ++i) {
      ind = id.get(i);
      int i1 = ind[0];
      int i2 = ind[1];
      int i3 = ind[2];
      pile.add(ind);
      pt[i3][i2][i1] = 0.0f;
      float fi = f[i3][i2][i1];
      if (i1==0 || i1==n1-1) 
        dp[i3][i2] = (float)i1; 
      else {
        float fm = f[i3][i2][i1-1];
        float fp = f[i3][i2][i1+1];
        dp[i3][i2] = parabolicPeak(i1,fm,fi,fp);
      }
    }
    return count;
  }  
 
  
  // use 3 points to fit a parabolic curve and find its peak
  private float parabolicPeak(int i1, float um, float ui, float up) {
    float z = (float)i1;
    float a = um-up;
    float b = 2.0f*(um+up)-4.0f*ui;
    return (z+a/b);
  }
  private RecursiveGaussianFilter _rgf = new RecursiveGaussianFilter(2.0f);
  private int _dd = 1;
  private int _d1 = 1;
  private int _d2 = 1;
  private int _d3 = 1;
  private int _dm = 4;
  private int _np = 20;
}
