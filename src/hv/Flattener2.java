/****************************************************************************
Copyright (c) 2013, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package hv;

import edu.mines.jtk.dsp.*;
import static edu.mines.jtk.util.ArrayMath.*;


/**
 * Flattens and unflattens locally linear features in a 2D image.
 * <p>
 * We assume that a 2D image f(x1,x2) has coherent features that are locally
 * linear with finite (non-vertical) slope. This class computes and applies
 * coordinate mappings for flattening and unflattening such features.
 * In this version, the mappings change only the vertical coordinate.
 * <p> 
 * Let x1 denote the vertical coordinate for the image f(x1,x2). This class
 * computes a mapping function x1(u1,u2) such that features in a transformed
 * image g(u1,u2) = f(x1(u1,u2),u2) are approximately horizontal. This process
 * is often called "flattening", and the transformed image g(u1,u2) is "the
 * flattened image." Likewise, for any constant u1, the curve defined by the
 * function x1(u1,x2) is called a "horizon".
 * <p>
 * If the coordinates u1 and u2 are sampled finely enough, then the mapping
 * function x1(u1,u2) is invertible. Specifically, there exists an inverse
 * mapping u1(x1,x2) such that x1 = x1(u1(x1,x2),x2) for all x2. Currently,
 * samplings of u1 and u2 are set to be identical to those for x1 and x2. If
 * used directly, this sampling of the mapping x1(u1,u2) may cause aliasing of
 * the flattened image g(u1,u2), so that f(x1,x2) cannot be recovered by
 * unflattening.
 * <p>
 * Without additional constraints, the description above of the mapping
 * function x1(u1,u2) is ambiguous. For example, one could add any constant to
 * this mapping function, and features in the transformed image g(u1,u2) would
 * still be horizontal. The mapping function x1(u1,u2) is here constrained so
 * that x1(u1,u2) = u1 for the first and last sampled values of u1. In other
 * words, features at the top or bottom of an image are not shifted by
 * flattening or unflattening.
 * @author Dave Hale, Colorado School of Mines
 * @version 2013.01.29
 */
public class Flattener2 {

  /** Coordinate mappings u1(x1,x2) and x1(u1,u2). */
  public static class Mappings {
    
    /** Sampling for the 1st dimension (the vertical coordinate). */
    public Sampling s1;
    
    /** Sampling for the 2nd dimension. */
    public Sampling s2;

    /** Array of sampled u1(x1,x2). */
    public float[][] u1;
    
    /** Array of sampled x1(u1,u2). */
    public float[][] x1;

    /**
     * Uses these mappings to flatten the specified image.
     * @param f the image to flatten.
     * @return the flattened image.
     */
    public float[][] flatten(float[][] f) {
      return apply(x1,f);
    }

    /**
     * Uses these mappings to unflatten the specified image.
     * @param f the image to unflatten.
     * @return the unflattened image.
     */
    public float[][] unflatten(float[][] f) {
      return apply(u1,f);
    }

    /**
     * Gets the flattening shifts s(u1,u2) = u1 - x1(u1,u2).
     * @return the flattening shifts.
     */
    public float[][] getShiftsS() {
      int n1 = s1.getCount();
      int n2 = s2.getCount();
      float[][] s = new float[n2][n1];
      float d1 = (float)s1.getDelta();
      float f1 = (float)s1.getFirst();
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float u1i = f1+i1*d1;
          s[i2][i1] = u1i-x1[i2][i1];
        }
      }
      return s;
    }

    /**
     * Gets the unflattening shifts r(x1,x2) = u1(x1,x2) - x1.
     * @return the unflattening shifts.
     */
    public float[][] getShiftsR() {
      int n1 = s1.getCount();
      int n2 = s2.getCount();
      float[][] r = new float[n2][n1];
      float d1 = (float)s1.getDelta();
      float f1 = (float)s1.getFirst();
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float x1i = f1+i1*d1;
          r[i2][i1] = u1[i2][i1]-x1i;
        }
      }
      return r;
    }

    private Mappings(Sampling s1, Sampling s2, float[][] u1, float[][] x1) {
      this.s1 = s1;
      this.s2 = s2;
      this.u1 = u1;
      this.x1 = x1;
    }

    private float[][] apply(float[][] ux, float[][] f) {
      int n1 = s1.getCount();
      int n2 = s2.getCount();
      double d1 = s1.getDelta();
      double f1 = s1.getFirst();
      SincInterpolator si = new SincInterpolator();
      float[][] g = new float[n2][n1];
      for (int i2=0; i2<n2; ++i2)
        si.interpolate(n1,d1,f1,f[i2],n1,ux[i2],g[i2]);
      return g;
    }
  }

  /**
   * Sets the relative weight of the PDE dr(x1,x2)/dx1 = 0.
   * Increasing this weight will cause shifts r(x1,x2) and s(u1,u2) to vary
   * more slowly with vertical coordinates x1 and u1, respectively. A weight
   * of 1.0 will cause this equation to get as much weight as other PDEs that
   * cause contours of constant u1 = u1(x1,x2) to be aligned with coherent
   * linear image features.
   * @param w1 the weight.
   */
  public void setWeight1(double w1) {
    _weight1 = (float)w1;
  }

  /**
   * Sets half-widths for smoothings in 1st and 2nd dimensions.
   * These smoothings serve as a preconditioner; they accelerate convergence
   * of the iterative solver used to compute mappings.
   * @param sigma1 half-width for smoothing in 1st dimension, in samples.
   * @param sigma2 half-width for smoothing in 2nd dimension, in samples.
   */
  public void setSmoothings(double sigma1, double sigma2) {
    _sigma1 = (float)sigma1;
    _sigma2 = (float)sigma2;
  }

  /**
   * Sets parameters that control the number of solver iterations.
   * @param small stop iterations when error norm is reduced by this fraction.
   * @param niter maximum number of solver iterations.
   */
  public void setIterations(double small, int niter) {
    _small = (float)small;
    _niter = niter;
  }

  /**
   * Gets mappings computed from specified slopes and linearities.
   * @param s1 sampling of 1st dimension.
   * @param s2 sampling of 2nd dimension.
   * @param p2 array of slopes of image features.
   * @param el array of linearities of image features.
   */
  public Mappings getMappingsFromSlopes(
    Sampling s1, Sampling s2, float[][] p2, float[][] el1, float[][] el2,float[][] w1) 
  {
    // Sampling parameters.
    int n1 = s1.getCount();
    int n2 = s2.getCount();
    float d1 = (float)s1.getDelta();
    float d2 = (float)s2.getDelta();
    float f1 = (float)s1.getFirst();

    // If necessary, convert units for slopes to samples per sample.
    //if (d1!=d2)
      //p2 = mul(d2/d1,p2);

    // Compute shifts r(x1,x2), in samples.
    float[][] b = new float[n2][n1]; // right-hand side
    float[][] r = new float[n2][n1]; // shifts, in samples
    VecArrayFloat2 vb = new VecArrayFloat2(b);
    VecArrayFloat2 vr = new VecArrayFloat2(r);
    Smoother2 smoother2 = new Smoother2(n1,n2,_sigma1,_sigma2,el1,el2);
    A2 a2 = new A2(smoother2,w1,el2,p2);
    CgSolver cs = new CgSolver(_small,_niter);
    makeRhs(el1,p2,b);
    smoother2.applyTranspose(b);
    cs.solve(a2,vb,vr);
    smoother2.apply(r);
    cleanShifts(r);

    // Compute u1(x1,x2).
    float[][] u1 = r;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        float x1i = f1+i1*d1;
        u1[i2][i1] = x1i+r[i2][i1]*d1;
      }
    }

    // Compute x1(u1,u2).
    float[][] x1 = b;
    InverseInterpolator ii = new InverseInterpolator(s1,s1);
    for (int i2=0; i2<n2; ++i2) 
      ii.invert(u1[i2],x1[i2]);

    return new Mappings(s1,s2,u1,x1);
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private float _weight1 = 0.010f; // weight for dr/d1 = 0 equation
  private float _sigma1 = 24.0f; // precon smoothing extent for 1st dim
  private float _sigma2 = 24.0f; // precon smoothing extent for 2nd dim
  private float _small = 0.01f; // stop CG iterations if residuals small
  private int _niter = 1000; // maximum number of CG iterations

  // Conjugate-gradient operators.
  private static class A2 implements CgSolver.A {
    A2(Smoother2 s2, float[][] w1, float[][] wp, float[][] p2) {
      _s2 = s2;
      _w1 = w1;
      _wp = wp;
      _p2 = p2;
    }
    public void apply(Vec vx, Vec vy) {
      VecArrayFloat2 v2x = (VecArrayFloat2)vx;
      VecArrayFloat2 v2y = (VecArrayFloat2)vy;
      float[][] x = v2x.getArray();
      float[][] y = v2y.getArray();
      float[][] z = copy(x);
      _s2.apply(z);
      zero(y);
      applyLhs(_w1,_wp,_p2,z,y);
      _s2.applyTranspose(y);
    }
    private Smoother2 _s2;
    private float[][] _w1;
    private float[][] _wp;
    private float[][] _p2;
  }

  // Smoother used as a preconditioner. After smoothing, enforces zero-shift
  // boundary conditions at top and bottom.
  private static class Smoother2 {
    public Smoother2(
      int n1, int n2, float sigma1, float sigma2, float[][] el1, float[][] el2) 
    {
      _sigma1 = sigma1;
      _sigma2 = sigma2;
      _el1 = el1;
      _el2 = el2;
    }
    public void apply(float[][] x) {
      //removeAverage(_el1,x);
      smooth2(_sigma2,_el2,x);
      smooth1(_sigma1,_el1,x);
      zero1(x);
    }
    public void applyTranspose(float[][] x) {
      zero1(x);
      smooth1(_sigma1,_el1,x);
      smooth2(_sigma2,_el2,x);
      //removeAverage(_el1,x);
    }
    private float _sigma1,_sigma2;
    private float[][] _el1;
    private float[][] _el2;
    private void zero1(float[][] x) {
      int n1 = x[0].length;
      int n2 = x.length;
      for (int i2=0; i2<n2; ++i2) {
        x[i2][   0] = 0.0f;
        x[i2][n1-1] = 0.0f;
      }
    }
  }
  /*
  // Smoothing for dimension 1.
  private static void smooth1(float sigma, float[][] x) {
    if (sigma<=0.0f)
      return;
    RecursiveExponentialFilter.Edges edges =
      RecursiveExponentialFilter.Edges.OUTPUT_ZERO_VALUE;
    RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sigma);
    ref.setEdges(edges);
    ref.apply1(x,x);
  }
  */
  // we do not want 
  private static void removeAverage(float[][] w, float[][] x) {
    int n2 = x.length;
    int n1 = x[0].length;
    for (int i1=0; i1<n1; ++i1) {
      float sumw = 0.0f;
      float sumx = 0.0f;
      for (int i2=0; i2<n2; ++i2) {
        float wi = w[i2][i1];
        sumw += wi;
        sumx += x[i2][i1];
      }
      float avg = sumx/(float)n2;
      if (sumw<(float)n2) avg = 0.0f;
      for (int i2=0; i2<n2; ++i2) { 
        x[i2][i1] -= avg;
      }
    }
  }
  // Smoothing for dimension 1.
  private static void smooth1(float sigma, float[][] s, float[][] x) {
    if (sigma<1.0f)
      return;
    float c = 0.5f*sigma*sigma;
    int n2 = x.length;
    int n1 = x[0].length;
    float[] xt = zerofloat(n2);
    float[] yt = zerofloat(n2);
    float[] st = fillfloat(1.0f,n2);
    LocalSmoothingFilter lsf = new LocalSmoothingFilter();
    lsf.setPreconditioner(true);
    for (int i2=0; i2<n2; ++i2) {
      if (s!=null) {
        for (int i1=0; i1<n1; ++i1)
          st[i1] = s[i2][i1];
      }
      for (int i1=0; i1<n1; ++i1)
        xt[i1] = x[i2][i1];
      lsf.apply(c,st,xt,yt);
      for (int i1=0; i1<n1; ++i1)
        x[i2][i1] = yt[i1];
    }
  }

  // Smoothing for dimension 2.
  private static void smooth2(float sigma, float[][] s, float[][] x) {
    if (sigma<1.0f)
      return;
    float c = 0.5f*sigma*sigma;
    int n1 = x[0].length;
    int n2 = x.length;
    float[] st = fillfloat(1.0f,n2);
    float[] xt = zerofloat(n2);
    float[] yt = zerofloat(n2);
    LocalSmoothingFilter lsf = new LocalSmoothingFilter();
    lsf.setPreconditioner(true);
    for (int i1=0; i1<n1; ++i1) {
      if (s!=null) {
        for (int i2=0; i2<n2; ++i2)
          st[i2] = s[i2][i1];
      }
      for (int i2=0; i2<n2; ++i2)
        xt[i2] = x[i2][i1];
      //lsf.apply(c,st,xt,yt);
      lsf.apply(c,xt,yt);
      for (int i2=0; i2<n2; ++i2)
        x[i2][i1] = yt[i2];
    }
  }

  private static void makeRhs(float[][] wp, float[][] p2, float[][] y) {
    int n1 = y[0].length;
    int n2 = y.length;
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        float wpi = (wp!=null)?wp[i2][i1]:1.000f;
        float p2i = p2[i2][i1];
        float b12 = wpi*p2i;
        float b22 = wpi;
        float x2 = -wpi*p2i;
        float y1 = b12*x2;
        float y2 = b22*x2;
        float ya = 0.5f*(y1+y2);
        float yb = 0.5f*(y1-y2);
        y[i2  ][i1  ] += ya;
        y[i2  ][i1-1] -= yb;
        y[i2-1][i1  ] += yb;
        y[i2-1][i1-1] -= ya;
      }
    }
  }
  private static void applyLhs(
    float[][] w1, float[][] wp, float[][] p2, float[][] x, float[][] y) 
  {
    int n1 = x[0].length;
    int n2 = x.length;
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        float wpi = (wp!=null)?wp[i2][i1]:1.000f;
	      float w1i = w1[i2][i1];
        float p2i = p2[i2][i1];
        float wps = wpi*wpi;
        float w1s = w1i*w1i;
        //float w1s = 0.0f;
        float p2s = p2i*p2i;
        float d11 = wps*p2s+w1s;
        float d12 = wps*p2i;
        float d22 = wps;
        float xa = 0.0f;
        float xb = 0.0f;
        xa += x[i2  ][i1  ];
        xb -= x[i2  ][i1-1];
        xb += x[i2-1][i1  ];
        xa -= x[i2-1][i1-1];
        float x1 = 0.5f*(xa+xb);
        float x2 = 0.5f*(xa-xb);
        float y1 = d11*x1+d12*x2;
        float y2 = d12*x1+d22*x2;
        float ya = 0.5f*(y1+y2);
        float yb = 0.5f*(y1-y2);
        y[i2  ][i1  ] += ya;
        y[i2  ][i1-1] -= yb;
        y[i2-1][i1  ] += yb;
        y[i2-1][i1-1] -= ya;
      }
    }
  }

  // Post-processing of computed shifts to ensure monotonic u1.
  private static void cleanShifts(float[][] r) {
    int n1 = r[0].length;
    int n2 = r.length;
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        if (r[i2][i1]<=r[i2][i1-1]-0.99f)
          r[i2][i1] = r[i2][i1-1]-0.99f;
      }
    }
  }

  private static void printStats(String s, int i1, float[][] a) {
    int n2 = a.length;
    float amin = a[0][i1];
    float amax = a[0][i1];
    for (int i2=1; i2<n2; ++i2) {
      if (a[i2][i1]<amin) amin = a[i2][i1];
      if (a[i2][i1]>amax) amax = a[i2][i1];
    }
    System.out.println(s+": i1="+i1+" min="+amin+" max="+amax);
  }
}
